import multiprocessing
import os
import time

from graph import add_edges_to_G, add_edges_to_D, get_HNet_G, get_HNet_D
from parsers.kgml_parser import parse_kgml, add_ec_numbers


class NetworkBuilder(object):
    """Handles all required information for trail finding (CoMetGeNe)."""
    def __init__(self, G_init, kgml, args, exclusions, ec_numbers, dev_out):
        """Initializes networks for CoMetGeNe trail finding and stores
        additional required information.

        If network initialization takes too long, the process is aborted and the
        pathway is blacklisted for the gap parameters with which CoMetGeNe was
        executed.

        :param G_init: undirected graph representing gene neighborhood for a
            given species
        :param kgml: metabolic pathway in KGML format
        :param args: command-line arguments for CoMetGeNe.py
        :param exclusions: Exclusions object representing blacklisted pathways
        :param ec_numbers: dict associating a list of EC numbers (values) to R
            numbers (keys)
        :param dev_out: output device for results (file or stdout)
        """
        self.delta_G = args.delta_G
        self.delta_D = args.delta_D
        self.blacklisted = False
        self.start_pw = time.time()
        self.pathway = os.path.join(os.path.abspath(args.DIR), kgml)
        self.reactions, self.compounds = parse_kgml(self.pathway)

        # Initialize the required data structures. The initialization is started
        # as a process that is terminated if it takes longer than a set timeout.
        # Terminating the process also results in blacklisting the pathway for
        # which initialization took too long (such pathways are added to the
        # exclusions list).
        init_queue = multiprocessing.Queue()
        init_process = multiprocessing.Process(
            target=self.initialize_networks,
            args=(init_queue, G_init, ec_numbers,)
        )
        init_process.start()
        init_process.join(timeout=args.timeout)
        if init_process.is_alive():
            init_process.terminate()
            init_process.join()
            aborted = kgml + ": Aborted (initialization takes longer than "
            aborted += str(args.timeout) + " seconds)\n"
            dev_out.write(aborted + '\n')
            exclusions.blacklist(
                self, 
                os.path.join(os.path.abspath(args.DIR), kgml))
            self.blacklisted = True
        else:
            elapsed = kgml + ": Initialization took %.2f seconds" % (
                time.time() - self.start_pw)
            dev_out.write(elapsed + '\n')
            self.reactions = init_queue.get()
            self.G_reduced = init_queue.get()
            self.G = init_queue.get()
            self.D = init_queue.get()

    def initialize_networks(self, queue, G_init, ec_numbers):
        """Adds EC number information to reactions present in the pathway to be
        analyzed and creates the input graphs for the HNet algorithm.

        The queue argument is modified.

        :param queue: multiprocessing queue storing reaction information, gene
            neighborhood graph (with additional edges allowing to skip genes, if
            required by the gap parameter), and the input graphs for the HNet
            algorithm
        :param G_init: undirected graph representing gene neighborhood for a
            given species
        :param ec_numbers: dict associating a list of EC numbers (values) to R
            numbers (keys)
        """
        # Add EC number information for the reactions in the metabolic pathway.
        add_ec_numbers(self.reactions, ec_numbers)

        # Add edges to the initial undirected graph G_init.
        G_reduced = add_edges_to_G(G_init, self.reactions, self.delta_G)

        # Build the initial directed graph D_init.
        D_init = get_HNet_D(self.reactions, self.compounds)

        # Build the undirected graph G on the same vertex set as D_init.
        G = get_HNet_G(self.reactions, D_init, G_reduced)

        # Build the final directed graph D by adding edges to D_init as
        # specified by delta_D.
        D = add_edges_to_D(D_init, self.delta_D)

        # Add the required data structures to the multiprocessing queue. They
        # must be retrieved in the same order as they are added to the queue.
        queue.put(self.reactions)
        queue.put(G_reduced)
        queue.put(G)
        queue.put(D)

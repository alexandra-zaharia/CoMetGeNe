import os

from ..definitions import EXCLUSIONS_FILENAME
from parsers.kgml_parser import get_pathway_id
import datetime


class Exclusions(object):
    """Handles blacklisted pathways."""
    def __init__(self, org):
        """Initializes exclusions for the specified species.

        :param org: species for which exclusions will be read
        """
        self.exclusions_file = EXCLUSIONS_FILENAME
        self.org = org
        self.exclusions = self.get_excluded_pathways()

    def get_excluded_pathways(self):
        """Returns excluded pathways for a given species.

        :return: dict of excluded pathways for a given species storing
            combinations of the gap parameters (delta_D for skipped reactions
            and delta_G for skipped genes) for which the pathways have been
            blacklisted
        """
        exclusions = dict()
        if os.path.isfile(self.exclusions_file):
            with open(self.exclusions_file, 'r') as f_in:
                for line in f_in:
                    if line.startswith(self.org):
                        # Parse the line corresponding to the query organism.
                        fields = line.split()
                        assert len(fields) == 4
                        for field in fields[1:4]:
                            assert field.isdigit()
                        delta_G = int(fields[1])
                        delta_D = int(fields[2])
                        pathway = fields[3]

                        if pathway not in exclusions:
                            exclusions[pathway] = list()
                        exclusions[pathway].append((delta_G, delta_D))
        return exclusions

    def can_analyze(self, kgml, delta_G, delta_D):
        """Determines whether the pathway designated by 'kgml' is analyzable,
        according to the gap parameters in the exclusions list.

        :param kgml: filename for input metabolic pathway
        :param delta_G: integer designating the value of the gap parameter
            delta_G (maximum number of skipped genes) with which CoMetGeNe is
            executed
        :param delta_D: integer designating the value of the gap parameter
            delta_D (maximum number of skipped reactions) with which CoMetGeNe
            is executed
        :return: True if kgml has not already been blacklisted for the given
            species and smaller or equal values of delta_G and delta_D, False
            otherwise
        """
        kegg_id = get_pathway_id(kgml)

        if int(kegg_id) >= 1100:  # exclude meta-pathways
            return False
        
        if kegg_id not in self.exclusions:
            return True
        
        for (dG, dD) in self.exclusions[kegg_id]:
            if delta_G >= dG and delta_D >= dD:
                return False
            
        return True

    def blacklist(self, network_instance, kgml):
        """Blacklists the specified pathway (adds it to the list of pathways
        excluded from analysis).

        The excluded pathway is also displayed on stdout, along with the gap
        parameters and the timestamp when it was added to the list of excluded
        pathways.

        :param network_instance: NetworkBuilder object
        :param kgml: filename for input metabolic pathway
        """
        pathway = get_pathway_id(kgml)
        with open(self.exclusions_file, 'a') as blacklist:
            entry = "%s %d %d %s" % (
                self.org,
                network_instance.delta_G,
                network_instance.delta_D,
                pathway)
            blacklist.write(entry + '\n')
            print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), entry

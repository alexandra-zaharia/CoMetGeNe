"""
Provides basic error handling for CoMetGeNe.

Version: 1.0 (May 2018)
License: MIT
Author: Alexandra Zaharia (contact@alexandra-zaharia.org)
"""

error = {
    # Errors related to kegg_import.py
    # --------------------------------
    'import_org_code':   1,  # invalid organism code or name
    'import_not_dir':    2,  # specified path is not a directory
    'import_not_r':      3,  # specified output directory is not readable
    'import_not_w':      4,  # specified output directory is not writable
    'import_not_x':      5,  # specified output directory is not executable
    'import_not_found':  6,  # specified organism not present in KEGG database
    'import_mkdir':      7,  # unknown error when creating the output directory
    'import_ec':         8,  # error retrieving reaction/EC numbers associations

    # Errors related to HNet.py
    # --------------------------------
    'hnet_node':       101,  # node in path not in graph D
    'hnet_edge':       102   # edge in path not in graph D
}


class CoMetGeNeError(Exception):
    """Creates instances of errors specific to CoMetGeNe."""
    def __init__(self, value, filename, arg2, arg3=None, arg4=None):
        """Creates a CoMetGeNeError object.

        :param value: string representing the type of error
        :param filename: file where the error is created and raised
        :param arg2: argument for error
        :param arg3: argument for error (optional)
        :param arg4: argument for error (optional)
        """
        self.value = value
        self.text = filename + ': '
        if error[value] <= 100:
            self.text += self.import_error_string(arg2)
        else:
            self.text += self.hnet_error_string(arg2, arg3, arg4)
        self.text += "Aborting."

    def import_error_string(self, arg2):
        """This CoMetGeNeError instance is related to downloads and permissions
        on directories storing downloads in kegg_import.py.

        :param arg2: string representing a KEGG organism code or a directory
        :return: error string
        """
        text = ''
        if self.value == 'import_org_code':
            text += "'%s' is not a valid KEGG organism code " % arg2
            text += "(three or four letters expected). "
        elif self.value == 'import_not_dir':
            text += "'%s' is not a directory. " % arg2
        elif self.value == 'import_not_r':
            text += "'%s' is not readable. " % arg2
        elif self.value == 'import_not_w':
            text += "'%s' is not writable. " % arg2
        elif self.value == 'import_not_x':
            text += "'%s' is not executable. " % arg2
        elif self.value == 'import_not_found':
            text += "Cannot retrieve metabolic pathways for '%s': " % arg2
            text += "organism not found in the KEGG database. "
        elif self.value == 'import_mkdir':
            text += "Could not create directory '%s'. " % arg2
        elif self.value == 'import_ec':
            text += "Could not retrieve associations between reactions and " + \
                    "EC numbers."
        return text

    def hnet_error_string(self, arg2, arg3=None, arg4=None):
        """This CoMetGeNeError instance is related to errors raised when
        performing trail finding (HNet algorithm) in HNet.py.

        :param arg2: string representing a vertex in a graph
        :param arg3: string representing a vertex in a graph (optional)
        :param arg4: list of strings representing a path in a graph (optional)
        :return: error string
        """
        text = ''
        if self.value == 'hnet_node':
            text += "Graph D does not contain node %s. " % arg2
        elif self.value == 'hnet_edge':
            assert arg3 is not None and arg4 is not None
            text += "There is no edge from %s to %s " % (arg2, arg3)
            text += "in path %s. " % arg4
        return text

    def __str__(self):
        """Returns a textual representation of this CoMetGeNeError instance.

        :return: textual representation of the error
        """
        return self.text

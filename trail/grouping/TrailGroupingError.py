import sys

from ..definitions import PICKLE_RN_FILENAME, PICKLE_RS_FILENAME, \
    PICKLE_GEN_FILENAME, PICKLE_GENOME_FILENAME


class TrailGroupingError(Exception):
    """Creates instances of errors specific to trail grouping."""
    @staticmethod
    def display_error():
        """Displays an error message on stderr instructing the user to rename
        pickle files that were automatically created on another data set."""
        sys.stderr.write(
            "\n\n###########\n#  ERROR\n###########\n#\n"
            "# It would appear that you are trying to perform trail grouping "
            "on\n"
            "# CoMetGeNe results different than those used to generate the\n"
            "# pickle files in the pickle/ directory.\n#\n"
            "# Try the following:\n"
            "# \t* rename " + PICKLE_RN_FILENAME + " to something else\n"
            "# \t* rename " + PICKLE_RS_FILENAME + " to something else\n"
            "# \t* rename " + PICKLE_GEN_FILENAME + " to something else\n"
            "# \t* rename " + PICKLE_GENOME_FILENAME + " to something else and"
            "\n# \t  use the correct file instead (created by CoMetGeNe.py)\n\n"
        )

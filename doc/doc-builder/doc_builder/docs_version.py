"""
A class defining characteristics of a documentation version
"""


class DocsVersion:
    """
    A class defining characteristics of a documentation version
    """

    # pylint: disable=too-few-public-methods,too-many-arguments
    def __init__(
        self,
        *,
        short_name,
        display_name,
        ref,
        landing_version=False,
    ):
        # The name of this version in file/URL paths
        self.short_name = short_name

        # What gets shown in the dropdown menu
        self.display_name = display_name

        # Whether this version should be the one on the landing page (i.e., default version)
        self.landing_version = landing_version

        # Branch, tag, or commit SHA
        self.ref = ref

    def subdir(self):
        """
        Get the subdirectory under --publish-dir where this version's HTML will be moved
        """
        if self.landing_version:
            return ""
        return self.short_name

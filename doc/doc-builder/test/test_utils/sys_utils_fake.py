"""Fake implementations of some system utilities"""

def make_fake_isdir(dirs_exist):
    """Return a fake replacement for os.path.isdir, which returns true
    if the given directory is in the list dirs_exist, false
    otherwise.
    
    dirs_exist: list of strings (paths)
    """

    dirs_exist_no_trailing_slashes = [thedir.rstrip("/") for thedir in dirs_exist]

    def isdir(thedir):
        thedir_no_trailing_slashes = thedir.rstrip("/")
        return thedir_no_trailing_slashes in dirs_exist_no_trailing_slashes

    return isdir

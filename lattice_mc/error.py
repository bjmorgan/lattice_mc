class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class BlockedLatticeError(Error):
    """Exception raised if there are no possible moves.

    Attributes:
        message: explanation of the error
    """

    def __init__(self, message):
        self.message = message

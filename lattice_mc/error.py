"""Exceptions raised by lattice_mc."""

from __future__ import annotations


class Error(Exception):
    """Base class for exceptions in this module."""

    pass


class BlockedLatticeError(Error):
    """Exception raised if there are no possible moves.

    Attributes:
        message: Explanation of the error.
    """

    def __init__(self, message: str) -> None:
        """
        Initialise a BlockedLatticeError.

        Args:
            message: Explanation of the error.
        """
        super().__init__(message)
        self.message = message

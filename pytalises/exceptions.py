"""Public exceptions for pyTALISES."""


class PyTalisesError(Exception):
    """Base exception for all pyTALISES errors."""


class BackendError(PyTalisesError):
    """Backend operation failed."""


class GridMismatchError(PyTalisesError):
    """Wavefunction and potential have incompatible grids."""


class ValidationError(PyTalisesError):
    """Invalid user input."""


class SerializationError(PyTalisesError):
    """Failed to save or load wavefunction."""

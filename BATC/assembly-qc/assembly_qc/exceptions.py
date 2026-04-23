"""Custom exception hierarchy for assembly-qc.

All public functions raise subclasses of AssemblyQCError so callers can
catch a single base class when fine-grained distinction is unnecessary.
"""


class AssemblyQCError(Exception):
    """Base exception for all assembly-qc errors."""


class FastaParseError(AssemblyQCError):
    """Raised when a FASTA file cannot be parsed or is structurally invalid."""


class FilterError(AssemblyQCError):
    """Raised when filter criteria are inconsistent or cannot be applied."""


class ComparisonError(AssemblyQCError):
    """Raised when assembly comparison inputs are invalid."""


class PlotError(AssemblyQCError):
    """Raised when a plot cannot be generated or saved."""

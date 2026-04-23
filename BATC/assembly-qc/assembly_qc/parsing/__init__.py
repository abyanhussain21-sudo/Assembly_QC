"""FASTA parsing utilities."""

from assembly_qc.parsing.fasta_parser import (
    ContigRecord,
    parse_fasta,
    parse_fasta_streaming,
    compute_gc,
    validate_fasta_path,
)

__all__ = [
    "ContigRecord",
    "parse_fasta",
    "parse_fasta_streaming",
    "compute_gc",
    "validate_fasta_path",
]

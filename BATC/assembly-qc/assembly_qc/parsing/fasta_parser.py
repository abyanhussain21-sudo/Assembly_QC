"""FASTA parsing layer wrapping Biopython SeqIO.

This module is the single entry point for all FASTA reading in assembly-qc.
All other modules receive List[ContigRecord] — they never call Bio.SeqIO directly.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Generator, List

from Bio import SeqIO

from assembly_qc.exceptions import FastaParseError

# Recognised FASTA file extensions (case-insensitive)
_FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".fas", ".fsa"}

# Bases counted as GC; everything else (including IUPAC ambiguity codes) is
# either AT or excluded from the denominator when it is a gap/N.
_GC_BASES = frozenset("GCgc")
_N_BASES = frozenset("Nn")


@dataclass(frozen=True)
class ContigRecord:
    """Immutable record for a single contig or scaffold.

    Attributes:
        seq_id:     Sequence identifier (text before first space in header).
        length:     Sequence length in base pairs.
        gc_content: GC fraction in [0.0, 1.0]; excludes N/gap characters from
                    denominator so assembly gaps don't dilute the estimate.
        n_count:    Number of ambiguous N bases in the sequence.
        sequence:   Raw uppercase nucleotide string.
    """

    seq_id: str
    length: int
    gc_content: float
    n_count: int
    sequence: str


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def compute_gc(sequence: str) -> float:
    """Return GC fraction (0.0–1.0) for *sequence*, excluding N/gap characters.

    Pure calculation exposed for direct testing and reuse.

    Args:
        sequence: Nucleotide string; case-insensitive.

    Returns:
        GC fraction.  Returns 0.0 if the sequence consists entirely of
        ambiguous characters (avoids division-by-zero).

    Examples:
        >>> compute_gc("GGCCAATT")
        0.5
        >>> compute_gc("NNNN")
        0.0
    """
    sequence = sequence.upper()
    # Denominator excludes Ns so assembly gaps don't deflate GC estimates
    effective_bases = sum(1 for b in sequence if b not in _N_BASES)
    if effective_bases == 0:
        return 0.0
    gc_count = sum(1 for b in sequence if b in _GC_BASES)
    return gc_count / effective_bases


def validate_fasta_path(fasta_path: Path) -> None:
    """Check that *fasta_path* exists, is a file, and has a FASTA extension.

    Args:
        fasta_path: Path to validate.

    Raises:
        FastaParseError: If the path does not exist, is not a regular file, or
            has an unrecognised extension.
    """
    if not fasta_path.exists():
        raise FastaParseError(f"File not found: {fasta_path}")
    if not fasta_path.is_file():
        raise FastaParseError(f"Path is not a regular file: {fasta_path}")
    if fasta_path.suffix.lower() not in _FASTA_EXTENSIONS:
        raise FastaParseError(
            f"Unrecognised extension '{fasta_path.suffix}' for {fasta_path}. "
            f"Expected one of: {', '.join(sorted(_FASTA_EXTENSIONS))}"
        )


def parse_fasta(fasta_path: Path) -> List[ContigRecord]:
    """Parse a FASTA file and return all sequences as a list of ContigRecord.

    Loads the entire file into memory.  For assemblies >1 GB consider
    :func:`parse_fasta_streaming` instead.

    Args:
        fasta_path: Path to the input FASTA file.

    Returns:
        List of ContigRecord, one per sequence.  Empty sequences are skipped
        with a warning (they carry no biological information).

    Raises:
        FastaParseError: If the file fails path validation, is truly empty
            (zero bytes), or contains no valid nucleotide records.
    """
    validate_fasta_path(fasta_path)

    # Zero-byte file check — Bio.SeqIO silently returns an empty iterator
    if fasta_path.stat().st_size == 0:
        raise FastaParseError(f"File is empty (0 bytes): {fasta_path}")

    records: List[ContigRecord] = []
    try:
        # Use an explicit file handle so it is closed deterministically after parsing,
        # avoiding ResourceWarning from Biopython's generator-based SeqIO.parse().
        with open(fasta_path, "r") as fh:
            for bio_record in SeqIO.parse(fh, "fasta"):
                seq_str = str(bio_record.seq).upper()
                if len(seq_str) == 0:
                    # Skip zero-length records; they don't contribute to assembly stats
                    continue
                n_count = sum(1 for b in seq_str if b in _N_BASES)
                records.append(
                    ContigRecord(
                        seq_id=bio_record.id,
                        length=len(seq_str),
                        gc_content=compute_gc(seq_str),
                        n_count=n_count,
                        sequence=seq_str,
                    )
                )
    except FastaParseError:
        raise
    except Exception as exc:
        # Re-raise Biopython errors as our own so callers don't need to import Bio
        raise FastaParseError(f"Failed to parse {fasta_path}: {exc}") from exc

    if not records:
        raise FastaParseError(
            f"No valid sequences found in {fasta_path}. "
            "Is this a valid FASTA file?"
        )

    return records


def parse_fasta_streaming(fasta_path: Path) -> Generator[ContigRecord, None, None]:
    """Streaming FASTA parser for large assemblies.

    Yields one ContigRecord at a time without loading all sequences into RAM.
    Useful for assemblies that would exceed available memory as a list.

    Args:
        fasta_path: Path to the input FASTA file.

    Yields:
        ContigRecord for each non-empty sequence.

    Raises:
        FastaParseError: If the file fails path validation or is empty.
    """
    validate_fasta_path(fasta_path)

    if fasta_path.stat().st_size == 0:
        raise FastaParseError(f"File is empty (0 bytes): {fasta_path}")

    try:
        with open(fasta_path, "r") as fh:
            for bio_record in SeqIO.parse(fh, "fasta"):
                seq_str = str(bio_record.seq).upper()
                if len(seq_str) == 0:
                    continue
                n_count = sum(1 for b in seq_str if b in _N_BASES)
                yield ContigRecord(
                    seq_id=bio_record.id,
                    length=len(seq_str),
                    gc_content=compute_gc(seq_str),
                    n_count=n_count,
                    sequence=seq_str,
                )
    except FastaParseError:
        raise
    except Exception as exc:
        raise FastaParseError(f"Failed to parse {fasta_path}: {exc}") from exc

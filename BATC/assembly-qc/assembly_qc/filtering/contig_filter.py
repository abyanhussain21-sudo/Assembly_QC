"""Contig filtering by length, GC content, N-fraction, and sequence ID patterns.

Design rationale:
  FilterCriteria is a frozen dataclass so it can be safely passed between
  functions without mutation.  apply_filters() validates criteria consistency
  before processing, giving actionable error messages before any I/O occurs.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from assembly_qc.exceptions import FilterError
from assembly_qc.parsing.fasta_parser import ContigRecord


@dataclass(frozen=True)
class FilterCriteria:
    """All filtering parameters in one immutable structure.

    None values mean "not applied".  Numeric thresholds are inclusive.

    Attributes:
        min_length:      Minimum contig length in bp (inclusive).
        max_length:      Maximum contig length in bp (inclusive).
        min_gc:          Minimum GC fraction [0.0, 1.0] (inclusive).
        max_gc:          Maximum GC fraction [0.0, 1.0] (inclusive).
        max_n_fraction:  Discard contigs where N/total > this value [0.0, 1.0].
        name_include:    Regex pattern; retain only contigs whose seq_id matches.
        name_exclude:    Regex pattern; discard contigs whose seq_id matches.
    """

    min_length: Optional[int] = None
    max_length: Optional[int] = None
    min_gc: Optional[float] = None
    max_gc: Optional[float] = None
    max_n_fraction: Optional[float] = None
    name_include: Optional[str] = None
    name_exclude: Optional[str] = None


@dataclass
class FilterResult:
    """Outcome of applying FilterCriteria to a list of contigs.

    Attributes:
        retained:  Contigs that passed every active criterion.
        discarded: Contigs that failed at least one criterion.
        criteria:  The FilterCriteria applied.
    """

    retained: List[ContigRecord]
    discarded: List[ContigRecord]
    criteria: FilterCriteria

    @property
    def retention_rate(self) -> float:
        """Fraction of input contigs that passed the filter (0.0–1.0)."""
        total = len(self.retained) + len(self.discarded)
        if total == 0:
            return 0.0
        return len(self.retained) / total


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------


def _validate_criteria(criteria: FilterCriteria) -> None:
    """Check that all criteria are internally consistent.

    Args:
        criteria: FilterCriteria to validate.

    Raises:
        FilterError: On any logical inconsistency.
    """
    if criteria.min_length is not None and criteria.min_length < 0:
        raise FilterError(f"min_length must be non-negative, got {criteria.min_length}")

    if criteria.max_length is not None and criteria.max_length < 0:
        raise FilterError(f"max_length must be non-negative, got {criteria.max_length}")

    if (
        criteria.min_length is not None
        and criteria.max_length is not None
        and criteria.min_length > criteria.max_length
    ):
        raise FilterError(
            f"min_length ({criteria.min_length}) > max_length ({criteria.max_length}). "
            "No contig can satisfy both constraints."
        )

    if criteria.min_gc is not None and not 0.0 <= criteria.min_gc <= 1.0:
        raise FilterError(f"min_gc must be in [0, 1], got {criteria.min_gc}")

    if criteria.max_gc is not None and not 0.0 <= criteria.max_gc <= 1.0:
        raise FilterError(f"max_gc must be in [0, 1], got {criteria.max_gc}")

    if (
        criteria.min_gc is not None
        and criteria.max_gc is not None
        and criteria.min_gc > criteria.max_gc
    ):
        raise FilterError(
            f"min_gc ({criteria.min_gc}) > max_gc ({criteria.max_gc})."
        )

    if criteria.max_n_fraction is not None and not 0.0 <= criteria.max_n_fraction <= 1.0:
        raise FilterError(
            f"max_n_fraction must be in [0, 1], got {criteria.max_n_fraction}"
        )

    # Validate regex patterns compile without error
    for attr in ("name_include", "name_exclude"):
        pattern = getattr(criteria, attr)
        if pattern is not None:
            try:
                re.compile(pattern)
            except re.error as exc:
                raise FilterError(
                    f"Invalid regex for {attr} '{pattern}': {exc}"
                ) from exc


# ---------------------------------------------------------------------------
# Core filter logic
# ---------------------------------------------------------------------------


def _passes_criteria(contig: ContigRecord, criteria: FilterCriteria) -> bool:
    """Return True if *contig* satisfies every active criterion.

    Each criterion is checked only when its value is not None, allowing callers
    to apply partial filters without any default behaviour.
    """
    if criteria.min_length is not None and contig.length < criteria.min_length:
        return False

    if criteria.max_length is not None and contig.length > criteria.max_length:
        return False

    if criteria.min_gc is not None and contig.gc_content < criteria.min_gc:
        return False

    if criteria.max_gc is not None and contig.gc_content > criteria.max_gc:
        return False

    if criteria.max_n_fraction is not None:
        n_fraction = contig.n_count / contig.length if contig.length > 0 else 0.0
        if n_fraction > criteria.max_n_fraction:
            return False

    if criteria.name_include is not None:
        if not re.search(criteria.name_include, contig.seq_id):
            return False

    if criteria.name_exclude is not None:
        if re.search(criteria.name_exclude, contig.seq_id):
            return False

    return True


def apply_filters(
    contigs: List[ContigRecord],
    criteria: FilterCriteria,
) -> FilterResult:
    """Apply all non-None criteria to *contigs* and return a FilterResult.

    Args:
        contigs:  Input contig list to filter.
        criteria: FilterCriteria specifying thresholds.

    Returns:
        FilterResult with retained/discarded lists and the applied criteria.

    Raises:
        FilterError: If criteria are logically inconsistent (validated first,
            before any contig is processed).
    """
    _validate_criteria(criteria)

    retained: List[ContigRecord] = []
    discarded: List[ContigRecord] = []

    for contig in contigs:
        if _passes_criteria(contig, criteria):
            retained.append(contig)
        else:
            discarded.append(contig)

    return FilterResult(retained=retained, discarded=discarded, criteria=criteria)


# ---------------------------------------------------------------------------
# FASTA output
# ---------------------------------------------------------------------------


def write_filtered_fasta(
    result: FilterResult,
    output_path: Path,
    original_fasta: Path,
) -> int:
    """Write retained contigs from *result* to a new FASTA file.

    Sequences are re-read from *original_fasta* via Bio.SeqIO to preserve the
    original record metadata (description, per-base quality if present).

    Args:
        result:         FilterResult from apply_filters().
        output_path:    Destination FASTA path.
        original_fasta: Source FASTA used to produce the contigs.

    Returns:
        Number of sequences written.

    Raises:
        FilterError: If output_path is not writable or the source cannot be read.
    """
    # Build a set of retained IDs for O(1) membership testing
    retained_ids = {c.seq_id for c in result.retained}

    try:
        # Filter records from source FASTA preserving original SeqRecord objects
        kept_records = [
            record
            for record in SeqIO.parse(str(original_fasta), "fasta")
            if record.id in retained_ids
        ]
        count = SeqIO.write(kept_records, str(output_path), "fasta")
    except OSError as exc:
        raise FilterError(f"Cannot write to {output_path}: {exc}") from exc
    except Exception as exc:
        raise FilterError(f"Error writing filtered FASTA: {exc}") from exc

    return count


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def filter_summary_dataframe(result: FilterResult) -> pd.DataFrame:
    """Build a per-contig filter report as a pandas DataFrame.

    Columns: seq_id, length, gc_content (%), n_count, passed_filter.

    Useful for TSV/CSV reports written with --report.

    Args:
        result: FilterResult from apply_filters().

    Returns:
        DataFrame with one row per input contig, sorted by length descending.
    """
    rows = []
    for contig in result.retained:
        rows.append(
            {
                "seq_id": contig.seq_id,
                "length": contig.length,
                "gc_content_pct": round(contig.gc_content * 100, 4),
                "n_count": contig.n_count,
                "passed_filter": True,
            }
        )
    for contig in result.discarded:
        rows.append(
            {
                "seq_id": contig.seq_id,
                "length": contig.length,
                "gc_content_pct": round(contig.gc_content * 100, 4),
                "n_count": contig.n_count,
                "passed_filter": False,
            }
        )

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("length", ascending=False).reset_index(drop=True)
    return df

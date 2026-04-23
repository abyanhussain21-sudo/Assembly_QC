"""Assembly statistics calculation.

Core metrics:
  - Total assembly size, contig count, min/max/mean/median length
  - N50, N90, L50, L90  (Nx = length at which x% of assembly is covered;
                          Lx = number of contigs needed to reach that point)
  - Assembly-wide GC content (length-weighted average across all contigs)
  - Total N-base count

All heavy computation uses pandas; results are stored in a frozen dataclass so
downstream code can rely on attribute access rather than dict key lookups.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import pandas as pd

from assembly_qc.exceptions import AssemblyQCError
from assembly_qc.parsing.fasta_parser import ContigRecord


@dataclass
class AssemblyStatistics:
    """Full statistics for a single assembly.  All length values are in bp."""

    fasta_path: Path
    total_contigs: int
    total_length: int
    min_length: int
    max_length: int
    mean_length: float
    median_length: float
    n50: int
    n90: int
    l50: int   # number of contigs needed to cover >= 50% of assembly length
    l90: int   # number of contigs needed to cover >= 90% of assembly length
    gc_content: float   # assembly-wide GC fraction (length-weighted)
    n_count: int        # total N bases across all contigs
    # raw length list excluded from repr to keep output readable
    contig_lengths: List[int] = field(repr=False)


# ---------------------------------------------------------------------------
# Nx / Lx calculation
# ---------------------------------------------------------------------------


def compute_nx(sorted_lengths: List[int], total_length: int, x: float) -> int:
    """Compute Nx length for a genome assembly.

    The Nx statistic is the length of the shortest contig in the set of
    longest contigs whose combined length accounts for at least x% of the
    total assembly size.  N50 (x=0.5) is the most common form.

    Args:
        sorted_lengths: Contig lengths sorted in **descending** order.
        total_length:   Sum of all contig lengths (the assembly span).
        x:              Coverage threshold as a fraction, e.g. 0.5 for N50.

    Returns:
        The Nx length in bp.  Returns 0 if the list is empty.

    Example (N50 for [50, 40, 30, 20, 10], total=150):
        target = 75; after 50: cumsum=50 < 75; after 40: cumsum=90 >= 75 → N50=40
    """
    if not sorted_lengths or total_length == 0:
        return 0

    target = total_length * x
    cumulative = 0
    for length in sorted_lengths:
        cumulative += length
        if cumulative >= target:
            return length
    # Should not normally reach here if sorted_lengths sums to total_length
    return sorted_lengths[-1]


def compute_lx(sorted_lengths: List[int], total_length: int, x: float) -> int:
    """Compute Lx: the count of contigs required to reach Nx coverage.

    Args:
        sorted_lengths: Contig lengths sorted in **descending** order.
        total_length:   Sum of all contig lengths.
        x:              Coverage threshold as a fraction, e.g. 0.5 for L50.

    Returns:
        Number of contigs.  Returns 0 if the list is empty.
    """
    if not sorted_lengths or total_length == 0:
        return 0

    target = total_length * x
    cumulative = 0
    for i, length in enumerate(sorted_lengths, start=1):
        cumulative += length
        if cumulative >= target:
            return i
    return len(sorted_lengths)


# ---------------------------------------------------------------------------
# Main statistics computation
# ---------------------------------------------------------------------------


def compute_statistics(contigs: List[ContigRecord], fasta_path: Path = Path(".")) -> AssemblyStatistics:
    """Compute all assembly statistics from a list of ContigRecord.

    Args:
        contigs:    List of ContigRecord returned by parse_fasta().
        fasta_path: Source FASTA path stored for reference in the result.

    Returns:
        Populated AssemblyStatistics instance.

    Raises:
        AssemblyQCError: If *contigs* is empty.
    """
    if not contigs:
        raise AssemblyQCError(
            "Cannot compute statistics: contig list is empty. "
            "Check that the FASTA file contains at least one sequence."
        )

    # Build a pandas Series of lengths — makes percentile/sum/mean trivial
    lengths_series = pd.Series([c.length for c in contigs], dtype="int64")

    total_length = int(lengths_series.sum())
    total_n = sum(c.n_count for c in contigs)

    # Length-weighted GC: sum(gc_i * len_i) / total_length
    # This gives the true assembly GC rather than the average of per-contig GC,
    # which would be biased by short contigs.
    weighted_gc = sum(c.gc_content * c.length for c in contigs) / total_length

    # Sort lengths descending once; reused for Nx and Lx
    sorted_lengths = sorted((c.length for c in contigs), reverse=True)

    return AssemblyStatistics(
        fasta_path=fasta_path,
        total_contigs=len(contigs),
        total_length=total_length,
        min_length=int(lengths_series.min()),
        max_length=int(lengths_series.max()),
        mean_length=float(lengths_series.mean()),
        median_length=float(lengths_series.median()),
        n50=compute_nx(sorted_lengths, total_length, 0.50),
        n90=compute_nx(sorted_lengths, total_length, 0.90),
        l50=compute_lx(sorted_lengths, total_length, 0.50),
        l90=compute_lx(sorted_lengths, total_length, 0.90),
        gc_content=weighted_gc,
        n_count=total_n,
        contig_lengths=sorted_lengths,
    )


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------


def statistics_to_dataframe(stats: AssemblyStatistics) -> pd.DataFrame:
    """Convert an AssemblyStatistics object to a single-row pandas DataFrame.

    Each metric becomes a column.  The resulting shape is (1, 14).
    Suitable for --format tsv|csv|json output and for the compare subcommand.

    Args:
        stats: Computed assembly statistics.

    Returns:
        Single-row DataFrame with snake_case column names.
    """
    row = {
        "assembly":      str(stats.fasta_path),
        "total_contigs": stats.total_contigs,
        "total_length":  stats.total_length,
        "min_length":    stats.min_length,
        "max_length":    stats.max_length,
        "mean_length":   round(stats.mean_length, 1),
        "median_length": round(stats.median_length, 1),
        "n50":           stats.n50,
        "n90":           stats.n90,
        "l50":           stats.l50,
        "l90":           stats.l90,
        "gc_content":    round(stats.gc_content * 100, 4),  # store as percentage
        "n_count":       stats.n_count,
    }
    return pd.DataFrame([row])


def statistics_to_string(stats: AssemblyStatistics) -> str:
    """Render AssemblyStatistics as a human-readable two-column table.

    Metric names and values are aligned with fixed-width columns so the output
    looks clean in terminals and log files.

    Args:
        stats: Computed assembly statistics.

    Returns:
        Multi-line string ready for print() or file.write().
    """
    # Contig length distribution buckets (rough summary, not a histogram)
    lengths = stats.contig_lengths
    buckets = _length_distribution_summary(lengths)

    lines = [
        f"Assembly statistics: {stats.fasta_path}",
        "=" * 55,
        f"  {'Total contigs':<28} {stats.total_contigs:>12,}",
        f"  {'Total assembly size (bp)':<28} {stats.total_length:>12,}",
        f"  {'Largest contig (bp)':<28} {stats.max_length:>12,}",
        f"  {'Smallest contig (bp)':<28} {stats.min_length:>12,}",
        f"  {'Mean contig length (bp)':<28} {stats.mean_length:>12,.1f}",
        f"  {'Median contig length (bp)':<28} {stats.median_length:>12,.1f}",
        "-" * 55,
        f"  {'N50 (bp)':<28} {stats.n50:>12,}",
        f"  {'N90 (bp)':<28} {stats.n90:>12,}",
        f"  {'L50 (contigs)':<28} {stats.l50:>12,}",
        f"  {'L90 (contigs)':<28} {stats.l90:>12,}",
        "-" * 55,
        f"  {'GC content (%)':<28} {stats.gc_content * 100:>11.2f}%",
        f"  {'Total N bases':<28} {stats.n_count:>12,}",
        "-" * 55,
        "  Contig length distribution:",
    ]
    for label, count in buckets.items():
        lines.append(f"    {label:<26} {count:>8,} contigs")

    return "\n".join(lines)


def _length_distribution_summary(lengths: List[int]) -> dict:
    """Build an ordered dict of length-bucket → contig count for display.

    Buckets are chosen to mirror common quality thresholds used in genome
    assembly reports (e.g., >=500 bp is a common minimum for annotation).

    Args:
        lengths: Contig lengths (any order).

    Returns:
        Ordered dict mapping bucket label to count.
    """
    thresholds = [
        (">= 10 Mbp",   10_000_000),
        (">= 1 Mbp",     1_000_000),
        (">= 100 kbp",     100_000),
        (">= 10 kbp",       10_000),
        (">= 1 kbp",         1_000),
        (">= 500 bp",          500),
        ("< 500 bp",             0),
    ]
    counts: dict = {}
    for label, threshold in thresholds:
        if threshold == 0:
            # Final bucket: everything below 500 bp
            counts[label] = sum(1 for l in lengths if l < 500)
        else:
            counts[label] = sum(1 for l in lengths if l >= threshold)
    return counts

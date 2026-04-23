"""Side-by-side comparison of two or more genome assemblies.

Typical use-case: benchmarking haplotype-resolved assemblies, comparing
assembler versions, or tracking quality improvement across polishing rounds.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from assembly_qc.exceptions import ComparisonError
from assembly_qc.statistics.assembly_stats import AssemblyStatistics, statistics_to_dataframe

# Metrics that are exported to the comparison table (subset of all stats fields).
# Order determines column order in the output table.
_COMPARISON_METRICS = [
    "total_contigs",
    "total_length",
    "min_length",
    "max_length",
    "mean_length",
    "median_length",
    "n50",
    "n90",
    "l50",
    "l90",
    "gc_content",
    "n_count",
]


@dataclass
class ComparisonResult:
    """Side-by-side statistics for two or more assemblies.

    Attributes:
        assemblies: List of AssemblyStatistics, one per assembly.
        labels:     Maps fasta_path → human label used in output tables.
    """

    assemblies: List[AssemblyStatistics]
    labels: Dict[Path, str]

    def to_dataframe(self) -> pd.DataFrame:
        """Wide-format DataFrame: one row per assembly, one column per metric.

        The DataFrame index is the human label string.
        The 'assembly' column holds the file path.

        Returns:
            DataFrame of shape (n_assemblies, len(_COMPARISON_METRICS) + 1).
        """
        rows = []
        idx = []
        for stats in self.assemblies:
            df_row = statistics_to_dataframe(stats)
            row_dict = df_row.iloc[0].to_dict()
            rows.append(row_dict)
            idx.append(self.labels[stats.fasta_path])

        df = pd.DataFrame(rows, index=idx)
        # Reorder columns: assembly path first, then metrics in canonical order
        ordered_cols = ["assembly"] + [
            m for m in _COMPARISON_METRICS if m in df.columns
        ]
        return df[ordered_cols]

    def delta_series(self, metric: str) -> pd.Series:
        """Absolute and percentage difference for *metric* between two assemblies.

        Only valid when exactly two assemblies are present.

        Args:
            metric: Column name from the comparison table (e.g. "n50").

        Returns:
            Series with index ["absolute_diff", "pct_diff"].  Positive values
            mean the second assembly is larger.

        Raises:
            ComparisonError: If not exactly two assemblies, or metric unknown.
        """
        if len(self.assemblies) != 2:
            raise ComparisonError(
                f"delta_series() requires exactly 2 assemblies, got {len(self.assemblies)}"
            )

        df = self.to_dataframe()
        if metric not in df.columns:
            raise ComparisonError(
                f"Unknown metric '{metric}'. "
                f"Available: {', '.join(df.columns.tolist())}"
            )

        val_a = float(df.iloc[0][metric])
        val_b = float(df.iloc[1][metric])
        abs_diff = val_b - val_a
        # Avoid division by zero for metrics that are zero in assembly A
        pct_diff = (abs_diff / val_a * 100) if val_a != 0 else float("nan")

        return pd.Series(
            {"absolute_diff": abs_diff, "pct_diff": round(pct_diff, 4)},
            name=metric,
        )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compare_assemblies(
    stats_list: List[AssemblyStatistics],
    labels: Optional[List[str]] = None,
) -> ComparisonResult:
    """Wrap a list of AssemblyStatistics in a ComparisonResult.

    Args:
        stats_list: One AssemblyStatistics per assembly (minimum two).
        labels:     Optional human-readable labels, same length as stats_list.
                    Defaults to the stem of each fasta_path.

    Returns:
        ComparisonResult ready for to_dataframe() or delta_series().

    Raises:
        ComparisonError: If fewer than two assemblies are supplied, or if
            labels length mismatches stats_list length.
    """
    if len(stats_list) < 2:
        raise ComparisonError(
            f"At least two assemblies are required for comparison, got {len(stats_list)}."
        )

    if labels is not None and len(labels) != len(stats_list):
        raise ComparisonError(
            f"Number of labels ({len(labels)}) must match "
            f"number of assemblies ({len(stats_list)})."
        )

    label_map: Dict[Path, str] = {}
    for i, stats in enumerate(stats_list):
        if labels is not None:
            label_map[stats.fasta_path] = labels[i]
        else:
            # Use file stem as default label; append index on collision
            stem = stats.fasta_path.stem
            if stem in label_map.values():
                stem = f"{stem}_{i}"
            label_map[stats.fasta_path] = stem

    return ComparisonResult(assemblies=stats_list, labels=label_map)


def rank_assemblies(
    result: ComparisonResult,
    primary_metric: str = "n50",
    secondary_metric: str = "total_length",
    ascending: bool = False,
) -> pd.DataFrame:
    """Return the comparison DataFrame sorted by assembly quality metrics.

    Args:
        result:           ComparisonResult to rank.
        primary_metric:   Metric to sort by first (default: n50).
        secondary_metric: Tiebreaker metric (default: total_length).
        ascending:        False = highest value is rank 1 (better for N50).

    Returns:
        Sorted DataFrame with a "rank" column prepended.

    Raises:
        ComparisonError: If either metric is not present in the table.
    """
    df = result.to_dataframe()

    for metric in (primary_metric, secondary_metric):
        if metric not in df.columns:
            raise ComparisonError(
                f"Metric '{metric}' not found in comparison table. "
                f"Available: {', '.join(df.columns.tolist())}"
            )

    df_sorted = df.sort_values(
        by=[primary_metric, secondary_metric],
        ascending=ascending,
    ).copy()

    # Insert rank as first column (1-based, human-readable)
    df_sorted.insert(0, "rank", range(1, len(df_sorted) + 1))
    return df_sorted


def comparison_to_string(result: ComparisonResult) -> str:
    """Render a ComparisonResult as a human-readable aligned table.

    Args:
        result: ComparisonResult from compare_assemblies().

    Returns:
        Multi-line string suitable for terminal output.
    """
    df = result.to_dataframe()

    # Metric display names and format strings
    metric_formats = {
        "total_contigs": ("Total contigs",         "{:>12,}"),
        "total_length":  ("Total length (bp)",     "{:>12,}"),
        "max_length":    ("Largest contig (bp)",   "{:>12,}"),
        "min_length":    ("Smallest contig (bp)",  "{:>12,}"),
        "mean_length":   ("Mean length (bp)",      "{:>12,.1f}"),
        "median_length": ("Median length (bp)",    "{:>12,.1f}"),
        "n50":           ("N50 (bp)",              "{:>12,}"),
        "n90":           ("N90 (bp)",              "{:>12,}"),
        "l50":           ("L50 (contigs)",         "{:>12,}"),
        "l90":           ("L90 (contigs)",         "{:>12,}"),
        "gc_content":    ("GC content (%)",        "{:>11.2f}%"),
        "n_count":       ("Total N bases",         "{:>12,}"),
    }

    labels = df.index.tolist()
    col_width = max(len(lbl) for lbl in labels) + 4
    metric_col_width = 28

    header = f"  {'Metric':<{metric_col_width}}" + "".join(
        f" {lbl:>{col_width}}" for lbl in labels
    )
    separator = "=" * (metric_col_width + 4 + col_width * len(labels) + 2)

    lines = [separator, header, separator]

    for col, (display_name, fmt) in metric_formats.items():
        if col not in df.columns:
            continue
        values = df[col].tolist()
        row = f"  {display_name:<{metric_col_width}}"
        for val in values:
            try:
                formatted = fmt.format(val)
            except (ValueError, TypeError):
                formatted = f"{str(val):>{col_width}}"
            row += f" {formatted:>{col_width}}"
        lines.append(row)

    lines.append(separator)
    return "\n".join(lines)

"""Assembly visualisation using matplotlib.

All functions set the non-interactive Agg backend before importing pyplot,
ensuring the module can be used in headless environments (HPC clusters, CI).
The backend must be set before `import matplotlib.pyplot` — do not reorder.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Literal, Optional

import matplotlib
matplotlib.use("Agg")  # non-interactive backend; must come before pyplot import
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from assembly_qc.exceptions import PlotError
from assembly_qc.statistics.assembly_stats import AssemblyStatistics

PlotFormat = Literal["png", "svg", "pdf"]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _resolve_format(output_path: Path, fmt: Optional[str]) -> PlotFormat:
    """Infer plot format from output_path extension, falling back to fmt or 'png'.

    Args:
        output_path: Destination path (extension used for inference).
        fmt:         Explicit format override (may be None).

    Returns:
        Lowercase format string: 'png', 'svg', or 'pdf'.
    """
    if fmt is not None:
        return fmt.lower()  # type: ignore[return-value]
    suffix = output_path.suffix.lstrip(".").lower()
    if suffix in ("png", "svg", "pdf"):
        return suffix  # type: ignore[return-value]
    return "png"


def _save_figure(
    fig: plt.Figure,
    output_path: Path,
    fmt: PlotFormat,
    dpi: int,
) -> Path:
    """Save *fig* to *output_path* and close it to prevent memory leaks.

    Args:
        fig:         Matplotlib Figure to save.
        output_path: Destination file path.
        fmt:         Format string ('png', 'svg', 'pdf').
        dpi:         Raster resolution (ignored for SVG/PDF).

    Returns:
        Resolved output_path.

    Raises:
        PlotError: If the file cannot be written.
    """
    try:
        fig.savefig(str(output_path), format=fmt, dpi=dpi, bbox_inches="tight")
    except OSError as exc:
        raise PlotError(f"Cannot write plot to {output_path}: {exc}") from exc
    finally:
        # Always close to release memory, even if save failed
        plt.close(fig)
    return output_path.resolve()


def _format_bp(x: float, _pos: int) -> str:
    """Matplotlib ticker formatter: convert bp to human-readable units.

    Args:
        x:    Tick value in base pairs.
        _pos: Tick position (required by matplotlib; unused here).

    Returns:
        String like "1.5 Mbp", "500 kbp", or "200 bp".
    """
    if x >= 1_000_000:
        return f"{x / 1_000_000:.1f} Mbp"
    if x >= 1_000:
        return f"{x / 1_000:.0f} kbp"
    return f"{int(x)} bp"


# ---------------------------------------------------------------------------
# Public plot functions
# ---------------------------------------------------------------------------


def plot_length_distribution(
    stats: AssemblyStatistics,
    output_path: Path,
    bins: int = 50,
    log_scale: bool = False,
    fmt: Optional[str] = None,
    dpi: int = 150,
) -> Path:
    """Histogram of contig length distribution with N50 and N90 annotations.

    The x-axis shows contig length in bp (human-readable units); the y-axis
    shows the number of contigs in each bin.  Vertical dashed lines mark N50
    (blue) and N90 (orange) — the two most commonly reported assembly metrics.

    Args:
        stats:       AssemblyStatistics for the assembly to plot.
        output_path: Destination PNG/SVG/PDF path.
        bins:        Number of histogram bins (default 50).
        log_scale:   Apply log10 scale to the x-axis.  Useful when a few very
                     large scaffolds dominate the distribution.
        fmt:         Output format; inferred from output_path extension if None.
        dpi:         Raster resolution in dots per inch.

    Returns:
        Resolved path of the saved image.

    Raises:
        PlotError: If output_path is not writable or stats has no contigs.
    """
    if not stats.contig_lengths:
        raise PlotError("Cannot plot length distribution: no contigs in assembly.")

    resolved_fmt = _resolve_format(output_path, fmt)
    lengths = stats.contig_lengths  # already sorted descending; order doesn't matter for histogram

    fig, ax = plt.subplots(figsize=(10, 5))

    # Use log-spaced bins when log_scale is True to avoid all data in first bin
    if log_scale and min(lengths) > 0:
        bin_edges = np.logspace(
            np.log10(min(lengths)), np.log10(max(lengths)), bins + 1
        )
    else:
        bin_edges = bins

    ax.hist(lengths, bins=bin_edges, color="#4C72B0", edgecolor="white", linewidth=0.4)

    # Mark N50 and N90 with vertical lines
    ax.axvline(stats.n50, color="#DD8452", linestyle="--", linewidth=1.5,
               label=f"N50 = {stats.n50:,} bp")
    ax.axvline(stats.n90, color="#55A868", linestyle="--", linewidth=1.5,
               label=f"N90 = {stats.n90:,} bp")

    if log_scale:
        ax.set_xscale("log")

    # Apply human-readable bp tick labels on the x-axis
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(_format_bp))
    plt.xticks(rotation=30, ha="right")

    ax.set_xlabel("Contig length", fontsize=12)
    ax.set_ylabel("Number of contigs", fontsize=12)
    ax.set_title(
        f"Contig length distribution — {stats.fasta_path.name}\n"
        f"{stats.total_contigs:,} contigs · {stats.total_length:,} bp total",
        fontsize=11,
    )
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.35, linestyle=":")

    fig.tight_layout()
    return _save_figure(fig, output_path, resolved_fmt, dpi)


def plot_nx_curve(
    stats_list: List[AssemblyStatistics],
    labels: List[str],
    output_path: Path,
    fmt: Optional[str] = None,
    dpi: int = 150,
) -> Path:
    """Nx curve for one or more assemblies plotted on the same axes.

    The Nx curve shows the contig length at which x% of the assembly is
    covered (x on the x-axis, Nx length on the y-axis).  A sharper drop-off
    at high x values indicates a more contiguous assembly.

    Multiple assemblies are shown as separate coloured lines, making it easy
    to compare contiguity profiles side-by-side.

    Args:
        stats_list:  List of AssemblyStatistics (one per assembly).
        labels:      Human-readable labels for each assembly (legend entries).
        output_path: Destination image path.
        fmt:         Output format override.
        dpi:         Raster resolution.

    Returns:
        Resolved path of the saved image.

    Raises:
        PlotError: If stats_list is empty, or labels length mismatches.
    """
    if not stats_list:
        raise PlotError("Cannot plot Nx curve: stats_list is empty.")
    if len(stats_list) != len(labels):
        raise PlotError(
            f"len(stats_list)={len(stats_list)} != len(labels)={len(labels)}"
        )

    resolved_fmt = _resolve_format(output_path, fmt)

    # Compute Nx for x in [1, 100] for each assembly
    x_values = list(range(1, 101))

    fig, ax = plt.subplots(figsize=(10, 5))

    for stats, label in zip(stats_list, labels):
        from assembly_qc.statistics.assembly_stats import compute_nx
        sorted_lengths = stats.contig_lengths  # already sorted descending
        nx_values = [
            compute_nx(sorted_lengths, stats.total_length, x / 100)
            for x in x_values
        ]
        ax.plot(x_values, nx_values, linewidth=2, label=label)

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(_format_bp))
    ax.set_xlabel("x (%)", fontsize=12)
    ax.set_ylabel("Nx length", fontsize=12)
    ax.set_title("Nx curve — assembly contiguity comparison", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.35, linestyle=":")

    fig.tight_layout()
    return _save_figure(fig, output_path, resolved_fmt, dpi)

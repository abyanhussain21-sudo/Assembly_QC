"""Command-line interface for assembly-qc.

Architecture:
  build_parser()  — constructs the full argparse tree; called from main() and
                    from tests so the parser can be exercised without I/O.
  main()          — parses argv, validates inputs, dispatches to domain modules.
                    No business logic lives here.

Subcommands:
  stats    — compute and display assembly statistics
  filter   — filter sequences by length / GC / N-fraction / ID pattern
  compare  — side-by-side comparison of two or more assemblies
  plot     — generate length distribution histogram or Nx curve
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path
from typing import List, Optional

from assembly_qc import __version__
from assembly_qc.exceptions import AssemblyQCError


# ---------------------------------------------------------------------------
# Parser construction
# ---------------------------------------------------------------------------


def build_parser():
    """Construct and return the top-level ArgumentParser with all subcommands.

    Separated from main() so tests can call build_parser() without triggering
    sys.exit() when testing argument parsing edge-cases.
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="assembly-qc",
        description=textwrap.dedent(
            """\
            assembly-qc: genome assembly quality analysis toolkit.

            Analyse FASTA assemblies for contiguity, GC content, and length
            distributions.  Four subcommands cover the most common QC tasks
            used in genome assembly pipelines.
            """
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version", action="version", version=f"assembly-qc {__version__}"
    )

    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    subparsers.required = True

    _add_stats_parser(subparsers)
    _add_filter_parser(subparsers)
    _add_compare_parser(subparsers)
    _add_plot_parser(subparsers)

    return parser


def _add_stats_parser(subparsers) -> None:
    """Register the 'stats' subcommand."""
    import argparse

    p = subparsers.add_parser(
        "stats",
        help="Compute assembly statistics (N50, L50, GC content, etc.)",
        description=(
            "Parse a FASTA assembly and report total size, contig count, "
            "N50, N90, L50, L90, GC content, and a length distribution summary."
        ),
    )
    p.add_argument(
        "fasta",
        type=Path,
        help="Path to the input assembly FASTA file (.fa/.fasta/.fna/.fas).",
    )
    p.add_argument(
        "--output", "-o",
        type=Path,
        default=None,
        metavar="FILE",
        help="Write statistics to FILE instead of stdout.",
    )
    p.add_argument(
        "--format", "-f",
        choices=["human", "tsv", "csv", "json"],
        default="human",
        metavar="FORMAT",
        help="Output format: human (default), tsv, csv, or json.",
    )
    p.add_argument(
        "--min-contig-length",
        type=int,
        default=0,
        metavar="BP",
        help=(
            "Exclude contigs shorter than BP before computing statistics. "
            "Useful to match NCBI's >=200 bp reporting threshold. Default: 0."
        ),
    )


def _add_filter_parser(subparsers) -> None:
    """Register the 'filter' subcommand."""
    p = subparsers.add_parser(
        "filter",
        help="Filter sequences by length, GC content, N fraction, or ID pattern.",
        description=(
            "Apply one or more filters to a FASTA assembly and write passing "
            "sequences to a new FASTA.  A summary of retained/removed counts "
            "is printed to stderr."
        ),
    )
    p.add_argument(
        "fasta",
        type=Path,
        help="Path to the input assembly FASTA file.",
    )
    p.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        metavar="FILE",
        help="Destination FASTA file for retained sequences (required).",
    )
    p.add_argument(
        "--min-length",
        type=int,
        default=None,
        metavar="BP",
        help="Minimum contig length in bp (inclusive). Shorter contigs are removed.",
    )
    p.add_argument(
        "--max-length",
        type=int,
        default=None,
        metavar="BP",
        help="Maximum contig length in bp (inclusive). Longer contigs are removed.",
    )
    p.add_argument(
        "--min-gc",
        type=float,
        default=None,
        metavar="FRAC",
        help="Minimum GC fraction [0–1]. Contigs below this GC are removed.",
    )
    p.add_argument(
        "--max-gc",
        type=float,
        default=None,
        metavar="FRAC",
        help="Maximum GC fraction [0–1]. Contigs above this GC are removed.",
    )
    p.add_argument(
        "--max-n-fraction",
        type=float,
        default=None,
        metavar="FRAC",
        help=(
            "Maximum allowed N fraction per contig [0–1]. "
            "Contigs with more Ns than this fraction are removed."
        ),
    )
    p.add_argument(
        "--include-pattern",
        type=str,
        default=None,
        metavar="REGEX",
        help="Keep only contigs whose sequence ID matches this regex.",
    )
    p.add_argument(
        "--exclude-pattern",
        type=str,
        default=None,
        metavar="REGEX",
        help="Remove contigs whose sequence ID matches this regex.",
    )
    p.add_argument(
        "--report",
        type=Path,
        default=None,
        metavar="FILE",
        help="Write per-contig filter report (TSV) to FILE.",
    )
    p.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress filter summary output to stderr.",
    )


def _add_compare_parser(subparsers) -> None:
    """Register the 'compare' subcommand."""
    p = subparsers.add_parser(
        "compare",
        help="Compare two or more assemblies side-by-side.",
        description=(
            "Compute statistics for each FASTA and display them in a side-by-side "
            "table.  Useful for benchmarking haplotype-resolved assemblies or "
            "evaluating assembly improvements across polishing rounds."
        ),
    )
    p.add_argument(
        "fastas",
        type=Path,
        nargs="+",
        metavar="FASTA",
        help="Two or more assembly FASTA files to compare.",
    )
    p.add_argument(
        "--labels",
        type=str,
        nargs="+",
        metavar="LABEL",
        default=None,
        help=(
            "Human-readable labels for each assembly (same order as FASTA args). "
            "Defaults to the file stem of each path."
        ),
    )
    p.add_argument(
        "--output", "-o",
        type=Path,
        default=None,
        metavar="FILE",
        help="Write comparison table to FILE instead of stdout.",
    )
    p.add_argument(
        "--format", "-f",
        choices=["human", "tsv", "csv", "json"],
        default="human",
        metavar="FORMAT",
        help="Output format: human (default), tsv, csv, or json.",
    )
    p.add_argument(
        "--rank-by",
        type=str,
        default="n50",
        metavar="METRIC",
        help="Metric to sort assemblies by in the output table. Default: n50.",
    )
    p.add_argument(
        "--ascending",
        action="store_true",
        help="Sort in ascending order (default: descending, higher is better).",
    )
    p.add_argument(
        "--min-contig-length",
        type=int,
        default=0,
        metavar="BP",
        help="Exclude contigs shorter than BP before computing stats for each assembly.",
    )


def _add_plot_parser(subparsers) -> None:
    """Register the 'plot' subcommand."""
    p = subparsers.add_parser(
        "plot",
        help="Generate assembly visualisation plots (histogram or Nx curve).",
        description=(
            "Plot contig length distribution or Nx contiguity curve for one or "
            "more assemblies.  Output is saved as PNG, SVG, or PDF."
        ),
    )
    p.add_argument(
        "fastas",
        type=Path,
        nargs="+",
        metavar="FASTA",
        help="One or more assembly FASTA files to plot.",
    )
    p.add_argument(
        "--type", "-t",
        choices=["length-dist", "nx-curve"],
        default="length-dist",
        metavar="TYPE",
        help=(
            "Plot type: 'length-dist' (histogram, single assembly) or "
            "'nx-curve' (contiguity curve, one or more assemblies). "
            "Default: length-dist."
        ),
    )
    p.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        metavar="FILE",
        help="Output image path.  Extension infers format (.png/.svg/.pdf).",
    )
    p.add_argument(
        "--format", "-f",
        choices=["png", "svg", "pdf"],
        default=None,
        metavar="FORMAT",
        help="Override output format (default: inferred from --output extension).",
    )
    p.add_argument(
        "--bins",
        type=int,
        default=50,
        metavar="N",
        help="Number of histogram bins for length-dist. Default: 50.",
    )
    p.add_argument(
        "--log-scale",
        action="store_true",
        help="Use log scale on the x-axis (useful for highly fragmented assemblies).",
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=150,
        metavar="DPI",
        help="Output resolution in dots per inch (PNG only). Default: 150.",
    )
    p.add_argument(
        "--labels",
        type=str,
        nargs="+",
        metavar="LABEL",
        default=None,
        help="Legend labels for multi-assembly nx-curve plots.",
    )
    p.add_argument(
        "--min-contig-length",
        type=int,
        default=0,
        metavar="BP",
        help="Exclude contigs shorter than BP before plotting.",
    )


# ---------------------------------------------------------------------------
# Input validation helpers
# ---------------------------------------------------------------------------


def _validate_fasta_inputs(paths: List[Path], parser) -> None:
    """Check every path in *paths* for existence; emit parser.error on failure."""
    for p in paths:
        if not p.exists():
            parser.error(f"File not found: {p}")
        if not p.is_file():
            parser.error(f"Not a regular file: {p}")


def _validate_positive_int(value: int, name: str, parser) -> None:
    """Emit parser.error if *value* is negative."""
    if value < 0:
        parser.error(f"{name} must be non-negative, got {value}")


def _validate_gc_fraction(value: Optional[float], name: str, parser) -> None:
    """Emit parser.error if *value* is outside [0, 1]."""
    if value is not None and not 0.0 <= value <= 1.0:
        parser.error(f"{name} must be in [0, 1], got {value}")


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------


def _run_stats(args, parser) -> int:
    """Execute the 'stats' subcommand."""
    from assembly_qc.parsing import parse_fasta
    from assembly_qc.statistics import (
        compute_statistics,
        statistics_to_string,
        statistics_to_dataframe,
    )

    _validate_fasta_inputs([args.fasta], parser)
    _validate_positive_int(args.min_contig_length, "--min-contig-length", parser)

    contigs = parse_fasta(args.fasta)

    # Apply minimum contig length filter if requested
    if args.min_contig_length > 0:
        before = len(contigs)
        contigs = [c for c in contigs if c.length >= args.min_contig_length]
        print(
            f"[filter] Excluded {before - len(contigs)} contigs shorter than "
            f"{args.min_contig_length:,} bp; {len(contigs):,} remain.",
            file=sys.stderr,
        )

    stats = compute_statistics(contigs, fasta_path=args.fasta)

    if args.format == "human":
        output_text = statistics_to_string(stats)
    else:
        df = statistics_to_dataframe(stats)
        if args.format == "tsv":
            output_text = df.to_csv(sep="\t", index=False)
        elif args.format == "csv":
            output_text = df.to_csv(index=False)
        elif args.format == "json":
            output_text = df.to_json(orient="records", indent=2)
        else:
            output_text = statistics_to_string(stats)  # fallback

    if args.output:
        args.output.write_text(output_text)
        print(f"Statistics written to {args.output}", file=sys.stderr)
    else:
        print(output_text)

    return 0


def _run_filter(args, parser) -> int:
    """Execute the 'filter' subcommand."""
    from assembly_qc.parsing import parse_fasta
    from assembly_qc.filtering import (
        FilterCriteria,
        apply_filters,
        write_filtered_fasta,
        filter_summary_dataframe,
    )

    _validate_fasta_inputs([args.fasta], parser)

    # Validate numeric arguments before touching any file I/O
    if args.min_length is not None:
        _validate_positive_int(args.min_length, "--min-length", parser)
    if args.max_length is not None:
        _validate_positive_int(args.max_length, "--max-length", parser)
    if args.min_length is not None and args.max_length is not None:
        if args.min_length > args.max_length:
            parser.error(
                f"--min-length ({args.min_length}) > --max-length ({args.max_length})"
            )
    _validate_gc_fraction(args.min_gc, "--min-gc", parser)
    _validate_gc_fraction(args.max_gc, "--max-gc", parser)
    if args.min_gc is not None and args.max_gc is not None and args.min_gc > args.max_gc:
        parser.error(f"--min-gc ({args.min_gc}) > --max-gc ({args.max_gc})")
    _validate_gc_fraction(args.max_n_fraction, "--max-n-fraction", parser)

    contigs = parse_fasta(args.fasta)

    criteria = FilterCriteria(
        min_length=args.min_length,
        max_length=args.max_length,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        max_n_fraction=args.max_n_fraction,
        name_include=args.include_pattern,
        name_exclude=args.exclude_pattern,
    )

    result = apply_filters(contigs, criteria)
    count = write_filtered_fasta(result, args.output, args.fasta)

    if not args.quiet:
        total = len(result.retained) + len(result.discarded)
        print(
            f"Filter summary:\n"
            f"  Input sequences:    {total:>8,}\n"
            f"  Retained:           {len(result.retained):>8,} "
            f"({result.retention_rate * 100:.1f}%)\n"
            f"  Removed:            {len(result.discarded):>8,}\n"
            f"  Output written to:  {args.output}",
            file=sys.stderr,
        )

    if args.report:
        df = filter_summary_dataframe(result)
        df.to_csv(str(args.report), sep="\t", index=False)
        print(f"Per-contig report written to {args.report}", file=sys.stderr)

    return 0


def _run_compare(args, parser) -> int:
    """Execute the 'compare' subcommand."""
    from assembly_qc.parsing import parse_fasta
    from assembly_qc.statistics import compute_statistics
    from assembly_qc.comparison import compare_assemblies, rank_assemblies
    from assembly_qc.comparison.assembly_compare import comparison_to_string

    if len(args.fastas) < 2:
        parser.error("compare requires at least two FASTA files.")

    _validate_fasta_inputs(args.fastas, parser)

    if args.labels is not None and len(args.labels) != len(args.fastas):
        parser.error(
            f"Number of --labels ({len(args.labels)}) must match "
            f"number of FASTA files ({len(args.fastas)})."
        )

    _validate_positive_int(args.min_contig_length, "--min-contig-length", parser)

    stats_list = []
    for fasta_path in args.fastas:
        contigs = parse_fasta(fasta_path)
        if args.min_contig_length > 0:
            contigs = [c for c in contigs if c.length >= args.min_contig_length]
        stats_list.append(compute_statistics(contigs, fasta_path=fasta_path))

    result = compare_assemblies(stats_list, labels=args.labels)

    if args.format == "human":
        try:
            ranked = rank_assemblies(
                result,
                primary_metric=args.rank_by,
                ascending=args.ascending,
            )
            # Add rank to a string representation by injecting into the comparison string
            output_text = comparison_to_string(result)
        except AssemblyQCError as exc:
            # rank_by metric might not exist; fall back to unranked output
            print(f"Warning: {exc}", file=sys.stderr)
            output_text = comparison_to_string(result)
    else:
        df = result.to_dataframe()
        if args.format == "tsv":
            output_text = df.to_csv(sep="\t")
        elif args.format == "csv":
            output_text = df.to_csv()
        elif args.format == "json":
            output_text = df.to_json(orient="index", indent=2)
        else:
            output_text = comparison_to_string(result)

    if args.output:
        args.output.write_text(output_text)
        print(f"Comparison written to {args.output}", file=sys.stderr)
    else:
        print(output_text)

    return 0


def _run_plot(args, parser) -> int:
    """Execute the 'plot' subcommand."""
    from assembly_qc.parsing import parse_fasta
    from assembly_qc.statistics import compute_statistics
    from assembly_qc.visualisation import plot_length_distribution, plot_nx_curve

    _validate_fasta_inputs(args.fastas, parser)
    _validate_positive_int(args.min_contig_length, "--min-contig-length", parser)
    _validate_positive_int(args.bins, "--bins", parser)
    _validate_positive_int(args.dpi, "--dpi", parser)

    if args.type == "length-dist" and len(args.fastas) > 1:
        parser.error("length-dist supports exactly one FASTA file.")

    if args.labels is not None and len(args.labels) != len(args.fastas):
        parser.error(
            f"Number of --labels ({len(args.labels)}) must match "
            f"number of FASTA files ({len(args.fastas)})."
        )

    stats_list = []
    for fasta_path in args.fastas:
        contigs = parse_fasta(fasta_path)
        if args.min_contig_length > 0:
            contigs = [c for c in contigs if c.length >= args.min_contig_length]
        stats_list.append(compute_statistics(contigs, fasta_path=fasta_path))

    labels = args.labels or [s.fasta_path.stem for s in stats_list]

    if args.type == "length-dist":
        out = plot_length_distribution(
            stats_list[0],
            args.output,
            bins=args.bins,
            log_scale=args.log_scale,
            fmt=args.format,
            dpi=args.dpi,
        )
    else:  # nx-curve
        out = plot_nx_curve(
            stats_list,
            labels=labels,
            output_path=args.output,
            fmt=args.format,
            dpi=args.dpi,
        )

    print(f"Plot saved to {out}", file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main(argv: Optional[List[str]] = None) -> int:
    """Parse arguments and dispatch to the appropriate subcommand handler.

    Args:
        argv: Argument list (default: sys.argv[1:]).

    Returns:
        Integer exit code (0 = success, 1 = error).
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        dispatch = {
            "stats":   _run_stats,
            "filter":  _run_filter,
            "compare": _run_compare,
            "plot":    _run_plot,
        }
        handler = dispatch[args.command]
        return handler(args, parser)

    except AssemblyQCError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130

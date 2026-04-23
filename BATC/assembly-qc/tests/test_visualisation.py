"""Tests for assembly_qc.visualisation (assembly_plot module).

Tests:
  — test_plot_length_distribution_creates_file  — PNG file is created and non-empty
  — test_plot_nx_curve_creates_file             — Nx curve PNG is created and non-empty
  — test_plot_empty_stats_raises                — PlotError for assembly with no contigs
  — test_plot_nx_curve_empty_list_raises        — PlotError for empty stats_list
"""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch, MagicMock

from assembly_qc.exceptions import PlotError
from assembly_qc.parsing.fasta_parser import ContigRecord
from assembly_qc.statistics import compute_statistics
from assembly_qc.visualisation import plot_length_distribution, plot_nx_curve


def _make_contig(seq_id: str, length: int) -> ContigRecord:
    return ContigRecord(
        seq_id=seq_id,
        length=length,
        gc_content=0.5,
        n_count=0,
        sequence="A" * length,
    )


def _stats(lengths: list[int], name: str = "test"):
    contigs = [_make_contig(f"ctg{i}", l) for i, l in enumerate(lengths)]
    return compute_statistics(contigs, fasta_path=Path(f"{name}.fasta"))


class TestAssemblyPlot(unittest.TestCase):
    """Unit tests for visualisation functions."""

    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmpdir.name)

    def tearDown(self):
        self._tmpdir.cleanup()

    # ------------------------------------------------------------------
    # Test: plot_length_distribution creates a non-empty file
    # ------------------------------------------------------------------

    def test_plot_length_distribution_creates_file(self):
        """plot_length_distribution() saves a non-empty PNG file."""
        stats = _stats([500, 1_000, 2_000, 5_000, 10_000])
        out = self.tmpdir / "length_dist.png"

        result_path = plot_length_distribution(stats, out, bins=10, dpi=72)

        self.assertTrue(out.exists(), "Output PNG file was not created.")
        self.assertGreater(out.stat().st_size, 0, "Output PNG file is empty.")
        self.assertEqual(result_path, out.resolve())

    # ------------------------------------------------------------------
    # Test: plot_nx_curve creates a non-empty file
    # ------------------------------------------------------------------

    def test_plot_nx_curve_creates_file(self):
        """plot_nx_curve() saves a non-empty PNG when given two assemblies."""
        stats_a = _stats([100, 500, 1_000], name="asm_a")
        stats_b = _stats([1_000, 5_000, 10_000], name="asm_b")
        out = self.tmpdir / "nx_curve.png"

        result_path = plot_nx_curve(
            [stats_a, stats_b],
            labels=["Assembly A", "Assembly B"],
            output_path=out,
            dpi=72,
        )

        self.assertTrue(out.exists(), "Nx curve PNG was not created.")
        self.assertGreater(out.stat().st_size, 0, "Nx curve PNG is empty.")
        self.assertEqual(result_path, out.resolve())

    # ------------------------------------------------------------------
    # Test: PlotError when stats has no contigs
    # ------------------------------------------------------------------

    def test_plot_empty_stats_raises(self):
        """plot_length_distribution raises PlotError for an assembly with no contigs."""
        # Manually construct an AssemblyStatistics with an empty contig list
        from assembly_qc.statistics.assembly_stats import AssemblyStatistics
        empty_stats = AssemblyStatistics(
            fasta_path=Path("empty.fasta"),
            total_contigs=0,
            total_length=0,
            min_length=0,
            max_length=0,
            mean_length=0.0,
            median_length=0.0,
            n50=0, n90=0, l50=0, l90=0,
            gc_content=0.0,
            n_count=0,
            contig_lengths=[],
        )
        out = self.tmpdir / "should_not_exist.png"

        with self.assertRaises(PlotError):
            plot_length_distribution(empty_stats, out)

    # ------------------------------------------------------------------
    # Test: PlotError for empty stats_list in nx_curve
    # ------------------------------------------------------------------

    def test_plot_nx_curve_empty_list_raises(self):
        """plot_nx_curve raises PlotError when given an empty stats_list."""
        out = self.tmpdir / "nx_empty.png"
        with self.assertRaises(PlotError):
            plot_nx_curve([], labels=[], output_path=out)


if __name__ == "__main__":
    unittest.main()

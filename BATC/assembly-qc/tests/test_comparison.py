"""Tests for assembly_qc.comparison (assembly_compare module).

Tests:
  — test_compare_two_assemblies     — ComparisonResult.to_dataframe() shape and labels
  — test_rank_assemblies_by_n50    — first row in ranked table has higher N50
  — test_compare_single_raises     — ComparisonError when only one assembly supplied
  — test_delta_series              — delta_series() computes absolute and pct diff
  — test_labels_mismatch_raises    — ComparisonError if labels count != fastas count
"""

import unittest
from pathlib import Path

from assembly_qc.exceptions import ComparisonError
from assembly_qc.parsing.fasta_parser import ContigRecord
from assembly_qc.statistics import compute_statistics
from assembly_qc.comparison import compare_assemblies, rank_assemblies


def _make_contig(seq_id: str, length: int) -> ContigRecord:
    """Build a minimal ContigRecord."""
    return ContigRecord(
        seq_id=seq_id,
        length=length,
        gc_content=0.5,
        n_count=0,
        sequence="A" * length,
    )


def _stats(lengths: list[int], name: str):
    """Build AssemblyStatistics from a list of lengths."""
    contigs = [_make_contig(f"ctg{i}", l) for i, l in enumerate(lengths)]
    return compute_statistics(contigs, fasta_path=Path(f"{name}.fasta"))


class TestAssemblyComparison(unittest.TestCase):
    """Unit tests for compare_assemblies() and rank_assemblies()."""

    def setUp(self):
        """Create two assemblies with clearly different N50 values."""
        # Assembly A: small contigs, low N50
        self.stats_a = _stats([100, 200, 300, 400, 500], name="assembly_a")
        # Assembly B: large contigs, high N50
        self.stats_b = _stats([1_000, 5_000, 10_000, 50_000, 100_000], name="assembly_b")

    # ------------------------------------------------------------------
    # Test: compare_two_assemblies produces correct DataFrame
    # ------------------------------------------------------------------

    def test_compare_two_assemblies(self):
        """ComparisonResult.to_dataframe() has shape (2, N) with both labels."""
        result = compare_assemblies(
            [self.stats_a, self.stats_b],
            labels=["asm_A", "asm_B"],
        )
        df = result.to_dataframe()

        # Two rows, one per assembly
        self.assertEqual(df.shape[0], 2)
        # At least 12 metric columns (+ 'assembly' path column)
        self.assertGreaterEqual(df.shape[1], 12)

        # Both labels appear in the index
        self.assertIn("asm_A", df.index)
        self.assertIn("asm_B", df.index)

        # Key metrics are present as columns
        for col in ("n50", "n90", "l50", "total_contigs", "gc_content"):
            self.assertIn(col, df.columns, f"Missing column: {col}")

    # ------------------------------------------------------------------
    # Test: rank_assemblies sorts by N50 descending
    # ------------------------------------------------------------------

    def test_rank_assemblies_by_n50(self):
        """rank_assemblies() places the assembly with higher N50 first."""
        result = compare_assemblies([self.stats_a, self.stats_b])
        ranked = rank_assemblies(result, primary_metric="n50", ascending=False)

        # First row should have the higher N50 (assembly_b)
        first_row_n50 = ranked.iloc[0]["n50"]
        second_row_n50 = ranked.iloc[1]["n50"]
        self.assertGreater(first_row_n50, second_row_n50)

        # Rank column exists and starts at 1
        self.assertIn("rank", ranked.columns)
        self.assertEqual(ranked.iloc[0]["rank"], 1)

    # ------------------------------------------------------------------
    # Test: single assembly raises ComparisonError
    # ------------------------------------------------------------------

    def test_compare_single_raises(self):
        """compare_assemblies() raises ComparisonError with fewer than 2 assemblies."""
        with self.assertRaises(ComparisonError):
            compare_assemblies([self.stats_a])

    # ------------------------------------------------------------------
    # Test: delta_series computes differences correctly
    # ------------------------------------------------------------------

    def test_delta_series(self):
        """delta_series() returns correct absolute and percentage differences."""
        result = compare_assemblies([self.stats_a, self.stats_b])

        # Both assemblies have GC = 0.5, so gc_content% in the table is the same
        # Use total_contigs: both have 5 contigs → absolute_diff = 0
        delta = result.delta_series("total_contigs")
        self.assertEqual(delta["absolute_diff"], 0)
        self.assertAlmostEqual(delta["pct_diff"], 0.0, places=4)

        # total_length differs; verify sign and magnitude
        delta_len = result.delta_series("total_length")
        expected_diff = self.stats_b.total_length - self.stats_a.total_length
        self.assertEqual(delta_len["absolute_diff"], expected_diff)

    # ------------------------------------------------------------------
    # Test: labels count mismatch raises ComparisonError
    # ------------------------------------------------------------------

    def test_labels_mismatch_raises(self):
        """ComparisonError is raised when len(labels) != len(stats_list)."""
        with self.assertRaises(ComparisonError):
            compare_assemblies(
                [self.stats_a, self.stats_b],
                labels=["only_one_label"],   # should be 2
            )


if __name__ == "__main__":
    unittest.main()

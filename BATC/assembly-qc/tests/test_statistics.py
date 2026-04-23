"""Tests for assembly_qc.statistics (assembly_stats module).

Tests:
  4. test_n50_hand_calculated       — N50/L50 verified by hand on a small dataset
  5. test_statistics_dataframe      — to_dataframe shape and column content
  6. test_empty_contigs_raises      — compute_statistics([]) raises AssemblyQCError
  7. test_single_sequence_assembly  — single contig: N50 == L50 == contig length
"""

import unittest
from pathlib import Path

from assembly_qc.exceptions import AssemblyQCError
from assembly_qc.parsing.fasta_parser import ContigRecord
from assembly_qc.statistics import (
    AssemblyStatistics,
    compute_statistics,
    compute_nx,
    compute_lx,
    statistics_to_dataframe,
    statistics_to_string,
)


def _make_contig(seq_id: str, sequence: str) -> ContigRecord:
    """Build a ContigRecord from a raw sequence string."""
    seq = sequence.upper()
    gc = sum(1 for b in seq if b in "GC") / max(len(seq), 1)
    n_count = sum(1 for b in seq if b == "N")
    return ContigRecord(seq_id=seq_id, length=len(seq), gc_content=gc,
                        n_count=n_count, sequence=seq)


def _make_contigs_by_length(lengths: list[int]) -> list[ContigRecord]:
    """Create synthetic ContigRecords with specified lengths (all-A sequence)."""
    return [
        _make_contig(f"ctg{i}", "A" * l)
        for i, l in enumerate(lengths)
    ]


class TestAssemblyStatistics(unittest.TestCase):
    """Unit tests for N50, N90, L50, L90, and statistics formatting."""

    # ------------------------------------------------------------------
    # Test 4 — N50/L50 hand calculation
    # ------------------------------------------------------------------

    def test_n50_hand_calculated(self):
        """N50 and L50 match values derived by hand from a 5-contig assembly.

        Lengths: [10, 20, 30, 40, 50]  total = 150 bp
        Sorted descending: [50, 40, 30, 20, 10]
        50% of 150 = 75 bp
          cumsum after 50 = 50  (< 75)
          cumsum after 40 = 90  (>= 75) → N50 = 40, L50 = 2
        90% of 150 = 135 bp
          cumsum after 50+40+30 = 120 (< 135)
          cumsum after +20 = 140       (>= 135) → N90 = 20, L90 = 4
        """
        contigs = _make_contigs_by_length([10, 20, 30, 40, 50])
        stats = compute_statistics(contigs, fasta_path=Path("test.fasta"))

        self.assertEqual(stats.total_length, 150)
        self.assertEqual(stats.total_contigs, 5)
        self.assertEqual(stats.n50, 40)
        self.assertEqual(stats.l50, 2)
        self.assertEqual(stats.n90, 20)
        self.assertEqual(stats.l90, 4)
        self.assertEqual(stats.max_length, 50)
        self.assertEqual(stats.min_length, 10)

    # ------------------------------------------------------------------
    # Test 5 — statistics_to_dataframe shape and content
    # ------------------------------------------------------------------

    def test_statistics_dataframe_shape(self):
        """statistics_to_dataframe returns a (1, N) DataFrame with key columns."""
        contigs = _make_contigs_by_length([100, 200, 300, 400, 500])
        stats = compute_statistics(contigs, fasta_path=Path("asm.fasta"))
        df = statistics_to_dataframe(stats)

        # Must be a single-row DataFrame
        self.assertEqual(df.shape[0], 1)
        # Must have at least 12 metric columns
        self.assertGreaterEqual(df.shape[1], 12)

        # Key columns must exist
        for col in ("n50", "n90", "l50", "l90", "total_contigs", "gc_content"):
            self.assertIn(col, df.columns, f"Missing column: {col}")

        # The n50 value in the DataFrame must match AssemblyStatistics
        self.assertEqual(int(df.iloc[0]["n50"]), stats.n50)

    # ------------------------------------------------------------------
    # Test 6 — empty contig list raises AssemblyQCError
    # ------------------------------------------------------------------

    def test_empty_contigs_raises(self):
        """compute_statistics([]) raises AssemblyQCError."""
        with self.assertRaises(AssemblyQCError):
            compute_statistics([], fasta_path=Path("empty.fasta"))

    # ------------------------------------------------------------------
    # Test 7 — single-sequence assembly
    # ------------------------------------------------------------------

    def test_single_sequence_assembly(self):
        """For a single contig, N50 == the contig length and L50 == 1."""
        contigs = _make_contigs_by_length([1_000_000])
        stats = compute_statistics(contigs, fasta_path=Path("single.fasta"))

        self.assertEqual(stats.total_contigs, 1)
        self.assertEqual(stats.total_length, 1_000_000)
        self.assertEqual(stats.n50, 1_000_000)
        self.assertEqual(stats.l50, 1)
        self.assertEqual(stats.n90, 1_000_000)
        self.assertEqual(stats.l90, 1)

    # ------------------------------------------------------------------
    # Utility: compute_nx / compute_lx edge cases
    # ------------------------------------------------------------------

    def test_compute_nx_empty_returns_zero(self):
        """compute_nx on an empty list returns 0 without raising."""
        self.assertEqual(compute_nx([], 0, 0.5), 0)

    def test_compute_lx_empty_returns_zero(self):
        """compute_lx on an empty list returns 0 without raising."""
        self.assertEqual(compute_lx([], 0, 0.5), 0)

    def test_statistics_to_string_contains_key_labels(self):
        """statistics_to_string output contains expected metric labels."""
        contigs = _make_contigs_by_length([500, 1000, 2000])
        stats = compute_statistics(contigs, fasta_path=Path("asm.fasta"))
        text = statistics_to_string(stats)

        for label in ("N50", "N90", "L50", "GC content", "Total contigs"):
            self.assertIn(label, text, f"Missing label in output: {label}")


if __name__ == "__main__":
    unittest.main()

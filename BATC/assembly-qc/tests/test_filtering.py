"""Tests for assembly_qc.filtering (contig_filter module).

Tests:
  8.  test_min_length_filter           — contigs below threshold are removed
  9.  test_max_length_filter           — contigs above threshold are removed
  10. test_inconsistent_criteria_raises — min_length > max_length raises FilterError
  11. test_filter_retention_rate       — retention_rate property is correct
  12. test_filter_boundary_inclusive   — boundary values pass (inclusive thresholds)
"""

import tempfile
import unittest
from pathlib import Path

from assembly_qc.exceptions import FilterError
from assembly_qc.filtering import FilterCriteria, FilterResult, apply_filters
from assembly_qc.parsing.fasta_parser import ContigRecord


def _make_contig(seq_id: str, length: int, gc: float = 0.5, n_count: int = 0) -> ContigRecord:
    """Build a synthetic ContigRecord without needing real sequence data."""
    # Construct a plausible sequence of the right length and approximate GC
    gc_bases = int(length * gc)
    at_bases = length - gc_bases
    sequence = "G" * gc_bases + "A" * at_bases
    return ContigRecord(
        seq_id=seq_id,
        length=length,
        gc_content=gc,
        n_count=n_count,
        sequence=sequence,
    )


class TestContigFilter(unittest.TestCase):
    """Unit tests for apply_filters() and FilterCriteria validation."""

    def _contigs(self):
        """Return a standard set of 5 contigs with known lengths."""
        return [
            _make_contig("ctg1", 100),
            _make_contig("ctg2", 500),
            _make_contig("ctg3", 1_000),
            _make_contig("ctg4", 5_000),
            _make_contig("ctg5", 10_000),
        ]

    # ------------------------------------------------------------------
    # Test 8 — min_length filter
    # ------------------------------------------------------------------

    def test_min_length_filter(self):
        """Contigs below min_length threshold are discarded; others retained."""
        contigs = self._contigs()
        criteria = FilterCriteria(min_length=500)
        result = apply_filters(contigs, criteria)

        retained_lengths = sorted(c.length for c in result.retained)
        discarded_lengths = sorted(c.length for c in result.discarded)

        # Contigs 500, 1000, 5000, 10000 should be retained (>= 500)
        self.assertEqual(retained_lengths, [500, 1_000, 5_000, 10_000])
        # Contig 100 should be discarded
        self.assertEqual(discarded_lengths, [100])

        # Retention rate: 4 of 5
        self.assertAlmostEqual(result.retention_rate, 0.8, places=6)

    # ------------------------------------------------------------------
    # Test 9 — max_length filter
    # ------------------------------------------------------------------

    def test_max_length_filter(self):
        """Contigs above max_length threshold are discarded."""
        contigs = self._contigs()
        criteria = FilterCriteria(max_length=1_000)
        result = apply_filters(contigs, criteria)

        retained_lengths = sorted(c.length for c in result.retained)
        self.assertEqual(retained_lengths, [100, 500, 1_000])

    # ------------------------------------------------------------------
    # Test 10 — inconsistent criteria raises FilterError
    # ------------------------------------------------------------------

    def test_inconsistent_criteria_raises(self):
        """FilterError is raised when min_length > max_length."""
        contigs = self._contigs()
        bad_criteria = FilterCriteria(min_length=1_000, max_length=100)

        with self.assertRaises(FilterError):
            apply_filters(contigs, bad_criteria)

    # ------------------------------------------------------------------
    # Test 11 — retention_rate property
    # ------------------------------------------------------------------

    def test_filter_retention_rate(self):
        """retention_rate property equals len(retained) / total."""
        contigs = self._contigs()
        # Retain only contigs with length >= 1000 → 3 of 5 → rate = 0.6
        result = apply_filters(contigs, FilterCriteria(min_length=1_000))
        self.assertAlmostEqual(result.retention_rate, 0.6, places=6)

    # ------------------------------------------------------------------
    # Test 12 — boundary values are inclusive
    # ------------------------------------------------------------------

    def test_filter_boundary_inclusive(self):
        """Contigs exactly at min_length or max_length are retained."""
        contigs = [
            _make_contig("boundary_min", 500),
            _make_contig("boundary_max", 1_000),
            _make_contig("too_short",    499),
            _make_contig("too_long",    1_001),
        ]
        criteria = FilterCriteria(min_length=500, max_length=1_000)
        result = apply_filters(contigs, criteria)

        retained_ids = {c.seq_id for c in result.retained}
        discarded_ids = {c.seq_id for c in result.discarded}

        self.assertIn("boundary_min", retained_ids)
        self.assertIn("boundary_max", retained_ids)
        self.assertIn("too_short", discarded_ids)
        self.assertIn("too_long", discarded_ids)

    # ------------------------------------------------------------------
    # Additional edge cases
    # ------------------------------------------------------------------

    def test_no_criteria_retains_all(self):
        """FilterCriteria with all None retains every contig."""
        contigs = self._contigs()
        result = apply_filters(contigs, FilterCriteria())
        self.assertEqual(len(result.retained), len(contigs))
        self.assertEqual(len(result.discarded), 0)
        self.assertAlmostEqual(result.retention_rate, 1.0, places=6)

    def test_filter_empty_contig_list(self):
        """apply_filters on an empty list returns an empty FilterResult."""
        result = apply_filters([], FilterCriteria(min_length=500))
        self.assertEqual(len(result.retained), 0)
        self.assertEqual(len(result.discarded), 0)
        self.assertEqual(result.retention_rate, 0.0)

    def test_gc_filter(self):
        """GC fraction filter correctly removes contigs outside bounds."""
        contigs = [
            _make_contig("low_gc",  1000, gc=0.20),
            _make_contig("mid_gc",  1000, gc=0.50),
            _make_contig("high_gc", 1000, gc=0.80),
        ]
        criteria = FilterCriteria(min_gc=0.30, max_gc=0.70)
        result = apply_filters(contigs, criteria)

        retained_ids = {c.seq_id for c in result.retained}
        self.assertIn("mid_gc", retained_ids)
        self.assertNotIn("low_gc", retained_ids)
        self.assertNotIn("high_gc", retained_ids)

    def test_negative_min_length_raises(self):
        """FilterError is raised for a negative min_length value."""
        with self.assertRaises(FilterError):
            apply_filters([], FilterCriteria(min_length=-1))


if __name__ == "__main__":
    unittest.main()

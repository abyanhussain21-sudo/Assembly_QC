"""Tests for assembly_qc.parsing (fasta_parser module).

Tests:
  1. test_parse_valid_fasta         — valid multi-sequence FASTA round-trips correctly
  2. test_parse_empty_fasta_raises  — empty file raises FastaParseError
  3. test_compute_gc_known_sequence — GC calculation on hand-verified sequences
"""

import tempfile
import unittest
from pathlib import Path

from assembly_qc.exceptions import FastaParseError
from assembly_qc.parsing import compute_gc, parse_fasta, ContigRecord


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    """Write a minimal FASTA file from a list of (id, sequence) tuples."""
    with open(path, "w") as fh:
        for seq_id, seq in records:
            fh.write(f">{seq_id}\n{seq}\n")


# ---------------------------------------------------------------------------
# Test class
# ---------------------------------------------------------------------------

class TestFastaParser(unittest.TestCase):
    """Unit tests for FASTA parsing and GC computation."""

    def setUp(self):
        """Create a temporary directory that is cleaned up after each test."""
        self._tmpdir = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmpdir.name)

    def tearDown(self):
        self._tmpdir.cleanup()

    # ------------------------------------------------------------------
    # Test 1 — valid FASTA round-trips correctly
    # ------------------------------------------------------------------

    def test_parse_valid_fasta(self):
        """parse_fasta returns correct seq_id, length, and gc_content."""
        records_input = [
            ("contig1", "ATCGATCGATCG"),   # 12 bp, GC = 6/12 = 0.5
            ("contig2", "AAAA"),            # 4 bp, GC = 0.0
            ("contig3", "GCGCGCGC"),        # 8 bp, GC = 1.0
        ]
        fasta_path = self.tmpdir / "test.fasta"
        _write_fasta(fasta_path, records_input)

        result = parse_fasta(fasta_path)

        self.assertEqual(len(result), 3)

        # Verify first record
        r0 = result[0]
        self.assertEqual(r0.seq_id, "contig1")
        self.assertEqual(r0.length, 12)
        self.assertAlmostEqual(r0.gc_content, 0.5, places=4)

        # Verify all seq_ids are present
        ids = {r.seq_id for r in result}
        self.assertEqual(ids, {"contig1", "contig2", "contig3"})

    # ------------------------------------------------------------------
    # Test 2 — empty file raises FastaParseError
    # ------------------------------------------------------------------

    def test_parse_empty_fasta_raises(self):
        """parse_fasta raises FastaParseError for a zero-byte file."""
        empty_path = self.tmpdir / "empty.fasta"
        empty_path.write_text("")  # zero bytes

        with self.assertRaises(FastaParseError):
            parse_fasta(empty_path)

    # ------------------------------------------------------------------
    # Test 3 — GC calculation on known sequences
    # ------------------------------------------------------------------

    def test_compute_gc_known_sequence(self):
        """compute_gc returns correct fractions for manually verified inputs."""
        # 50% GC: GGCCAATT → 4 GC out of 8 effective bases
        self.assertAlmostEqual(compute_gc("GGCCAATT"), 0.5, places=6)

        # 100% GC
        self.assertAlmostEqual(compute_gc("GCGCGCGC"), 1.0, places=6)

        # 0% GC
        self.assertAlmostEqual(compute_gc("AAAATTTT"), 0.0, places=6)

        # All N: denominator is 0, should return 0.0 without division error
        self.assertEqual(compute_gc("NNNN"), 0.0)

        # Mixed with Ns: only non-N bases count in denominator
        # GCNN → 2 GC out of 2 effective bases = 1.0
        self.assertAlmostEqual(compute_gc("GCNN"), 1.0, places=6)

    # ------------------------------------------------------------------
    # Test 4 (bonus) — missing file raises FastaParseError
    # ------------------------------------------------------------------

    def test_parse_nonexistent_file_raises(self):
        """parse_fasta raises FastaParseError for a path that does not exist."""
        missing = self.tmpdir / "does_not_exist.fasta"
        with self.assertRaises(FastaParseError):
            parse_fasta(missing)


if __name__ == "__main__":
    unittest.main()

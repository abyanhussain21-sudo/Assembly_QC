# assembly-qc

A Python command-line toolkit for analysing genome assembly quality from FASTA files.  
Built for use in production genome assembly pipelines alongside tools like Merqury, QUAST, and BUSCO.

---

## What it does

`assembly-qc` provides four subcommands:

| Subcommand | Description |
|---|---|
| `stats`   | Compute N50, N90, L50, L90, GC content, contig count, length distribution |
| `filter`  | Filter contigs by length, GC fraction, N content, or sequence ID pattern |
| `compare` | Side-by-side statistics table for two or more assemblies |
| `plot`    | Contig length distribution histogram or Nx contiguity curve (PNG/SVG/PDF) |

---

## Installation

**Requirements:** Python ≥ 3.10

```bash
# Clone the repository
git clone https://github.com/example/assembly-qc.git
cd assembly-qc

# Install dependencies
pip install -r requirements.txt

# Install the package (adds the assembly-qc command to PATH)
pip install -e .
```

Or install dependencies directly:

```bash
pip install biopython>=1.83 pandas>=2.0 matplotlib>=3.8 numpy>=1.26
```

---

## Usage

### `stats` — compute assembly statistics

```bash
# Default human-readable output
assembly-qc stats genome.fasta

# Only count contigs >= 500 bp (common for NCBI submission thresholds)
assembly-qc stats genome.fasta --min-contig-length 500

# Export as TSV for downstream analysis
assembly-qc stats genome.fasta --format tsv --output stats.tsv

# Export as JSON
assembly-qc stats genome.fasta --format json --output stats.json
```

**Example output:**

```
Assembly statistics: genome.fasta
=======================================================
  Total contigs                          1,247
  Total assembly size (bp)         623,481,902
  Largest contig (bp)               45,231,018
  Smallest contig (bp)                     201
  Mean contig length (bp)              500,386.2
  Median contig length (bp)            128,442.0
-------------------------------------------------------
  N50 (bp)                          12,450,391
  N90 (bp)                           1,823,044
  L50 (contigs)                             18
  L90 (contigs)                             89
-------------------------------------------------------
  GC content (%)                         38.71%
  Total N bases                       3,421,009
-------------------------------------------------------
  Contig length distribution:
    >= 10 Mbp                               12 contigs
    >= 1 Mbp                                47 contigs
    >= 100 kbp                             203 contigs
    >= 10 kbp                              589 contigs
    >= 1 kbp                               891 contigs
    >= 500 bp                            1,102 contigs
    < 500 bp                               145 contigs
```

---

### `filter` — filter sequences

```bash
# Keep only contigs >= 500 bp
assembly-qc filter genome.fasta --output filtered.fasta --min-length 500

# Keep contigs between 1 kbp and 10 Mbp
assembly-qc filter genome.fasta --output filtered.fasta \
    --min-length 1000 --max-length 10000000

# Remove highly repetitive contigs (GC > 70%) and gap-heavy contigs (N > 5%)
assembly-qc filter genome.fasta --output filtered.fasta \
    --max-gc 0.70 --max-n-fraction 0.05

# Keep only chromosome-scale scaffolds by sequence ID pattern
assembly-qc filter genome.fasta --output chr_only.fasta \
    --include-pattern "^chr"

# Write a per-contig filter report
assembly-qc filter genome.fasta --output filtered.fasta \
    --min-length 500 --report filter_report.tsv
```

**Example stderr summary:**

```
Filter summary:
  Input sequences:        1,247
  Retained:               1,102 (88.4%)
  Removed:                  145
  Output written to:  filtered.fasta
```

---

### `compare` — compare two or more assemblies

```bash
# Compare two assemblies with default labels (file stem)
assembly-qc compare haplotype1.fasta haplotype2.fasta

# Provide human-readable labels
assembly-qc compare hap1.fasta hap2.fasta --labels "Haplotype 1" "Haplotype 2"

# Compare three polishing rounds, ranked by N50
assembly-qc compare raw.fasta polished_r1.fasta polished_r2.fasta \
    --labels "Raw" "Round 1" "Round 2" --rank-by n50

# Export comparison as CSV
assembly-qc compare asm_v1.fasta asm_v2.fasta \
    --format csv --output comparison.csv
```

**Example output:**

```
============================================================
  Metric                           Haplotype 1   Haplotype 2
============================================================
  Total contigs                          1,247         1,193
  Total length (bp)               623,481,902   621,903,451
  Largest contig (bp)              45,231,018    46,892,003
  N50 (bp)                         12,450,391    14,103,882
  N90 (bp)                          1,823,044     2,011,231
  L50 (contigs)                            18            16
  L90 (contigs)                            89            81
  GC content (%)                        38.71%        38.69%
  Total N bases                      3,421,009     2,891,044
============================================================
```

---

### `plot` — generate plots

```bash
# Contig length distribution histogram (saved as PNG)
assembly-qc plot genome.fasta --type length-dist --output length_dist.png

# Log-scale histogram for highly fragmented assemblies
assembly-qc plot genome.fasta --type length-dist --output length_dist.png \
    --log-scale --bins 100

# Nx contiguity curve comparing two assemblies
assembly-qc plot hap1.fasta hap2.fasta \
    --type nx-curve --output nx_curve.png \
    --labels "Haplotype 1" "Haplotype 2"

# High-resolution SVG for publication
assembly-qc plot genome.fasta --type length-dist \
    --output figure1.svg --format svg --dpi 300

# Filter short contigs before plotting
assembly-qc plot genome.fasta --type length-dist --output plot.png \
    --min-contig-length 500
```

---

## Running tests

```bash
# Run all 12 tests
python -m unittest discover tests/

# With verbose output
python -m unittest discover tests/ -v

# Run a specific test file
python -m unittest tests.test_statistics -v
```

---

## Project structure

```
assembly-qc/
├── assembly_qc/
│   ├── cli.py                      # argparse subcommands and dispatch
│   ├── exceptions.py               # custom exception hierarchy
│   ├── parsing/
│   │   └── fasta_parser.py         # Biopython SeqIO wrapper, GC calculation
│   ├── statistics/
│   │   └── assembly_stats.py       # N50/Nx, L50/Lx, GC, pandas formatting
│   ├── filtering/
│   │   └── contig_filter.py        # length/GC/N/regex filtering
│   ├── comparison/
│   │   └── assembly_compare.py     # side-by-side table, delta_series
│   └── visualisation/
│       └── assembly_plot.py        # matplotlib histogram and Nx curve
├── tests/
│   ├── test_parsing.py             # 4 tests: parsing correctness and edge cases
│   ├── test_statistics.py          # 6 tests: N50/L50 arithmetic, DataFrame, single-seq
│   ├── test_filtering.py           # 8 tests: boundaries, GC, N-fraction, errors
│   ├── test_comparison.py          # 5 tests: table shape, ranking, delta, errors
│   └── test_visualisation.py       # 4 tests: file creation, edge-case errors
├── requirements.txt
└── README.md
```

---

## Key design decisions

- **Biopython SeqIO** is used exclusively for FASTA I/O; all downstream modules receive `List[ContigRecord]` and never import Bio directly.
- **pandas** is used for stats formatting and tabular output; all numeric calculations are done with plain Python for transparency.
- **N/GC denominator**: N bases are excluded from the GC denominator, so assembly scaffolding gaps do not artificially lower GC estimates.
- **Length-weighted GC**: Assembly-wide GC is computed as `sum(gc_i * len_i) / total_len`, not the mean of per-contig values, to give the true nucleotide composition.
- **Non-interactive matplotlib backend** (`Agg`) is set before any pyplot import, ensuring plots work on HPC nodes without a display.

---


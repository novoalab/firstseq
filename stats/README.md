# stats

This folder contains base-resolution statistics output from sequencing data analysis.

Each file, such as `FIRST_Invitro_DMS.STATS`, contains per-position metrics useful for quality control, mutation rate analysis, and base composition.

---

## File Format
Tab-separated values (`.STATS` files)

### Columns:
| Column     | Description                                             |
|------------|---------------------------------------------------------|
| `chr`      | Chromosome name                                         |
| `pos`      | Genomic position (1-based)                              |
| `ref_nuc`  | Reference nucleotide                                    |
| `coverage` | Number of reads covering this position                  |
| `meanQual` | Average base quality score                              |
| `medianQual` | Median base quality score                            |
| `rtstop`   | Fraction of reverse transcription stops                 |
| `ins`      | Fraction of insertions                                  |
| `del`      | Fraction of deletions                                   |
| `A` `T` `C` `G` | Count of each nucleotide observed at this position |

---

## Usage
These files are used for:
- Base modification detection
- RT stop profiling
- Mutation rate analysis
- Nucleotide bias assessment
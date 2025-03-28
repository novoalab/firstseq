# readends

This directory contains `.readends.tsv` files detailing read end positions captured from sequencing datasets.  
Useful for identifying termination sites or analyzing coverage profiles.

### File Naming Convention:
`{EXPERIMENT}_{ENZYME}_{CONDITION}.readends.tsv`

Examples:
- `dcDNA_Maxima.readends.tsv`: dcDNA with Maxima
- `FIRST_Mn_SS2.readends.tsv`: FIRST synthesis with Mn and SS2 condition

### Columns (expected):
- Chromosome
- Position
- Strand (optional)
- Read count at end position

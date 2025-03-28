# chrcounts

This folder contains chromosome-level read count data files generated from sequencing experiments.  
Each file is a `.chrcount.tsv` file with tab-separated values.

### File Naming Convention:
`{EXPERIMENT}_{ENZYME}_{CONDITION}.chrcount.tsv`

Examples:
- `dcDNA_Maxima.chrcount.tsv`: dcDNA experiment using Maxima enzyme
- `FIRST_Mg_PS3.chrcount.tsv`: FIRST-strand synthesis using Mg enzyme, PS3 condition

### Data Content:
Each file contains:
- Chromosome ID
- Read count
- Possibly strand-specific or positional data

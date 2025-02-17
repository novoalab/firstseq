# FIRST-seq Data Analysis Repository
![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)  
![R](https://img.shields.io/badge/Made%20with-R-blue)

## 🧬 Project Overview  
This repository contains the code, data, and manuscript for our study on **FIRST-seq**. Our approach leverages **customised nanopore direct cDNA sequencing** to investigate reverse transcriptase fingerprints and determine RNA secondary structures in long read sequencing reads. 

## 🛠 Installation  
### 1️⃣ **Clone this repository** 

```bash
git clone https://github.com/YOUR-USERNAME/RNA-Nanopore-Analysis.git
cd RNA-Nanopore-Analysis
```

### 2️⃣ **Set up dependencies** 

Ensure R and required packages are installed. We use renv for reproducibility.

```R
install.packages("renv")
renv::restore()  # Install all required packages
```

Alternatively, install key dependencies manually:

```R
install.packages(c("tidyverse", "ggplot2", "data.table", "Biostrings", "Seurat"))
```


📄 Reference Paper
This repo is created for the research paper below:

Oguzhan Begik, Gregor Diensthuber,Ivana Borovska,John S Mattick,Danny Incarnato,Eva Maria Novoa. "Long-read transcriptome-wide RNA structure maps using DMS-FIRST-seq" Under Review



📜 License
This project is licensed under the MIT License - see the LICENSE file for details.


📧 Contact
For questions, reach out to:
📩 oguzhanbegik@gmail.com
🔬 Novoa Lab / Centre for Genomic Regulation




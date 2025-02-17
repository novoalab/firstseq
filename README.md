# FIRST-seq Data Analysis Repository
![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)  
![R](https://img.shields.io/badge/Made%20with-R-blue)

## 🧬 Project Overview  
This repository contains the code, data, and manuscript for our study on **FIRST-seq**. Our approach leverages **customised nanopore direct cDNA sequencing** to investigate reverse transcriptase fingerprints and determine RNA secondary structures in long read sequencing reads. 

## 📂 Repository Structure  
📦 RNA-Nanopore-Analysis/
│
├── 📂 analysis/           🧬 R scripts for data analysis  
│   ├── 01_data_preprocessing.R  
│   ├── 02_quality_control.R  
│   ├── 03_structure_mapping.R  
│   ├── 04_visualization.R  
│   ├── 05_statistical_analysis.R  
│
├── 📂 data/               📊 Raw and processed sequencing data  
│   ├── 📂 raw/            🔬 Unprocessed sequencing data  
│   ├── 📂 processed/      🏗️ Cleaned & structured data  
│   ├── metadata.csv       📄 Sample metadata  
│
├── 📂 results/            📈 Figures, tables, and outputs  
│   ├── 📂 plots/          🎨 Data visualizations  
│   ├── 📂 tables/         📑 Summarized results  
│   ├── 📂 supplementary/  📎 Additional supporting materials  
│
├── 📂 scripts/            🛠️ Utility functions and helpers  
│   ├── functions.R  
│   ├── plotting.R  
│
├── 📂 notebooks/          📓 Reproducible RMarkdown analyses  
│   ├── exploratory_analysis.Rmd  
│   ├── final_paper_analysis.Rmd  
│
├── 📂 paper/              📜 Manuscript, references, and figures  
│   ├── manuscript.tex     ✍️ LaTeX/Markdown for the paper  
│   ├── references.bib     📚 Bibliography  
│   ├── 📂 figures/        🖼️ Manuscript-ready figures  
│
├── 📂 environment/        🔧 Dependency management files  
│   ├── renv.lock         📌 R environment dependencies  
│   ├── renv/             📂 R package library (if using `renv`)  
│
├── 📜 .gitignore          🚫 Ignore unnecessary files  
├── 📝 LICENSE             📄 License file  
├── 📖 README.md           📘 Project documentation  
├── 📂 RNA-Nanopore-Analysis.Rproj  🎯 RStudio project file  


## 🛠 Installation  
### 1️⃣ **Clone this repository** 

```bash
git clone https://github.com/YOUR-USERNAME/RNA-Nanopore-Analysis.git
cd RNA-Nanopore-Analysis
```





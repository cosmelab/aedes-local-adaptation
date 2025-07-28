# Tool Versions Report - Aedes Local Adaptation Analysis Container

This document provides a comprehensive inventory of all bioinformatics tools, software packages, and computational environments used in the Aedes aegypti local adaptation analysis pipeline. This information is essential for reproducibility and can be referenced in publications.

**Container Information:**
- Generated on: 2025-07-28 15:47:33 UTC
- Container: Aedes Local Adaptation Analysis
- Base OS: Ubuntu (via micromamba)
- Architecture: AMD64/x86_64

---

## Core Programming Languages and Environments

| Software | Version | Purpose |
|----------|---------|---------|
| Python | 3.11.7 | Primary scripting and data analysis |
| R | 4.3.3 | Statistical computing and population genetics |
| Bash | 5.2.15 | Shell scripting and pipeline automation |
| Java | 21.0.2-internal | Runtime for bioinformatics tools |
| Ruby | 3.2.2 | Utility scripts and gems |

## Primary Genomics and Bioinformatics Tools

| Tool | Version | Description |
|------|---------|-------------|
| samtools | 1.22 | SAM/BAM file manipulation and analysis |
| bcftools | Library issue | VCF/BCF file manipulation and analysis |
| vcftools | (0.1.17) | VCF file analysis and filtering |
| bedtools | v2.31.1 | Genomic interval manipulation |
| plink | 1.9.0 | Population genetics analysis |
| plink2 | 2.0.0 | Next-generation PLINK |
| bwa | 0.7.19-r1273 | Short read alignment |
| angsd | 0.940 | NGS data analysis |
| tabix | 1.22 | Genomic file indexing |
| iqtree | 3.0.1 | Phylogenetic inference |

## Population Genetics and Local Adaptation Tools

| Tool | Version | Description |
|------|---------|-------------|
| ADMIXTURE | 1.3.0 | Population structure analysis |
| BayeScan | 2.1 | FST outlier detection |
| GEMMA | 0.98.42021012920122021 | Genome-wide association studies |
| sambada | 0.8.0 | Landscape genomics analysis |
| BA3-SNPS | 3.0.4 | Bayesian migration rate estimation |

## Python Packages

| Package | Version | Description |
|---------|---------|-------------|
| numpy | 1.26.4 | Numerical computing foundation |
| pandas | 2.3.1 | Data manipulation and analysis |
| scipy | 1.16.0 | Scientific computing library |
| matplotlib | 3.10.3 | Plotting and visualization |
| seaborn | 0.13.2 | Statistical data visualization |
| scikit-learn | 1.7.1 | Machine learning library |
| scikit-allel | 1.3.13 | Population genetics analysis |
| cyvcf2 | 0.31.0 | Fast VCF file processing |
| pysam | 0.22.1 | SAM/BAM file interface |
| biopython | 1.85 | Biological computation |
| geopandas | 1.1.1 | Geospatial data analysis |
| rasterio | 1.3.10 | Raster data I/O |
| plotly | 6.2.0 | Interactive visualization |
| networkx | 2.8.8 | Network analysis |

## R Packages

### Core R Packages
| Package | Version | Description |
|---------|---------|-------------|
| tidyverse | 2.0.0 | Data science meta-package |
| data.table | 1.17.8 | High-performance data manipulation |
| ggplot2 | 3.5.2 | Grammar of graphics plotting |
| dplyr | 1.1.4 | Data manipulation verbs |
| devtools | 2.4.5 | Package development tools |
| remotes | 2.5.0 | Package installation from remote repositories |

### Population Genetics R Packages
| Package | Version | Description |
|---------|---------|-------------|
| adegenet | 2.1.11 | Exploratory analysis of genetic data |
| vcfR | 1.15.0 | VCF file manipulation in R |
| SNPRelate | 1.36.0 | Parallel computing for SNP data |
| LEA | 3.19.8 | Landscape and ecological association studies |
| lfmm | 1.0 | Latent factor mixed models |
| OutFLANK | 0.2 | FST outlier detection |
| pcadapt | 4.4.0 | Principal component analysis for population genetics |
| lostruct | 0.0.0.9000 | Local PCA for population structure |
| admixtools | 2.0.10 | ADMIXTOOLS wrapper for R |
| admixr | 0.10.0 | ADMIXTOOLS analysis in R |

### Geospatial and Visualization R Packages
| Package | Version | Description |
|---------|---------|-------------|
| sf | 1.0.16 | Simple features for spatial data |
| raster | 3.6.32 | Raster data analysis |
| terra | 1.8.60 | Spatial data analysis |
| leaflet | 2.2.2 | Interactive web maps |
| tmap | 3.3.4 | Thematic maps |
| ggstatsplot | 0.13.1 | Statistical plots with ggplot2 |
| plotly | 4.11.0 | Interactive web graphics |

## System Information

| Component | Details |
|-----------|---------|
| Operating System | Debian GNU/Linux 12 (bookworm) |
| Kernel | 4.18.0-477.27.1.el8_8.x86_64 |
| Architecture | x86_64 |
| Container Engine | micromamba 1.5.0 |
| Shell | /bin/bash |

---

## Usage Notes

This tool inventory was automatically generated and includes all software versions used in the Aedes aegypti local adaptation analysis pipeline. For reproducibility:

1. **Container Usage**: All analyses should be performed within the containerized environment to ensure identical software versions.

2. **Version Dependencies**: Some tools may have interdependent version requirements. The container environment resolves all compatibility issues.

3. **Citation Guidelines**: When publishing results, reference this tool inventory in supplementary materials and cite individual tools as appropriate.

4. **Updates**: This report reflects the software versions at the time of generation. Container updates may include newer versions.

## Reproducibility Statement

All computational analyses were performed using the software versions documented above within a containerized environment. The complete computational environment can be reproduced using the Dockerfile and associated configuration files in this repository.

For questions about specific tool versions or installation procedures, please refer to the project documentation or contact the development team through the project repository.

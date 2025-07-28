# 🦟 Aedes Local Adaptation Analysis - Tools & Versions

<div align="center">

![Container](https://img.shields.io/badge/Environment-Containerized-orange?style=for-the-badge&logo=docker)
![Package Manager](https://img.shields.io/badge/Package%20Manager-micromamba-green?style=for-the-badge&logo=anaconda)
![Architecture](https://img.shields.io/badge/Architecture-AMD64%20HPC-blue?style=for-the-badge&logo=linux)

**📦 Comprehensive bioinformatics environment for local adaptation analysis**

</div>

---

## 🐳 **Environment Setup**

This project uses a **containerized environment** with Docker/Singularity for reproducible analysis. All packages are pre-installed and optimized using **micromamba** for fast, reliable package management.

<div align="center">

### 🎯 **Container Approach Benefits**

<table>
<tr>
<td align="center" width="25%">

**🔄 Reproducible**
<br>
Exact same environment across different systems

</td>
<td align="center" width="25%">

**⚡ Optimized**
<br>
Pre-compiled packages with dependency resolution

</td>
<td align="center" width="25%">

**🖥️ HPC-ready**
<br>
Compatible with Singularity on HPC systems

</td>
<td align="center" width="25%">

**🚀 Fast startup**
<br>
No need to install packages individually

</td>
</tr>
</table>

</div>

### 💻 **Usage Examples**

<details>
<summary><b>🐳 Container Usage Commands</b> (Click to expand)</summary>

```zsh
# Docker
docker run -it ghcr.io/cosmelab/aedes-local-adaptation:latest

# Singularity (HPC) - Interactive shell
singularity exec aedes-local-adaptation.sif zsh
```

### 💾 **Cache Management for HPC**

Before pulling containers on HPC systems, set up proper cache management:

```zsh
# Check quota first
check_quota home

# Use temporary cache to avoid quota issues
mkdir -p ./singularity_temp_cache
export SINGULARITY_CACHEDIR=$PWD/singularity_temp_cache

# Pull container
singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest

# Clean up
rm -rf ./singularity_temp_cache
unset SINGULARITY_CACHEDIR
```

</details>

---

## 🛠️ **Core System Tools**

<table>
<tr>
<th width="20%">🔧 Tool</th>
<th width="20%">📦 Version</th>
<th width="60%">🎯 Purpose</th>
</tr>
<tr>
<td><b>🐚 zsh</b></td>
<td><code>Latest</code></td>
<td>Advanced shell environment with Oh-My-Zsh</td>
</tr>
<tr>
<td><b>🐍 Python</b></td>
<td><code>3.11.7</code></td>
<td>Primary scripting and data analysis</td>
</tr>
<tr>
<td><b>📊 R</b></td>
<td><code>4.3.3</code></td>
<td>Statistical computing and population genetics</td>
</tr>
<tr>
<td><b>⚡ Bash</b></td>
<td><code>5.2.15</code></td>
<td>Shell scripting and pipeline automation</td>
</tr>
<tr>
<td><b>☕ Java</b></td>
<td><code>21.0.2-internal</code></td>
<td>Runtime for bioinformatics tools</td>
</tr>
<tr>
<td><b>💎 Ruby</b></td>
<td><code>3.2.2</code></td>
<td>Utility scripts and gems</td>
</tr>
<tr>
<td><b>📦 micromamba</b></td>
<td><code>1.5.0</code></td>
<td>Fast package management (optimized for containers)</td>
</tr>
</table>

---

## 🧬 **Primary Genomics and Bioinformatics Tools**

<table>
<tr>
<th width="25%">🔧 Tool</th>
<th width="20%">📦 Version</th>
<th width="55%">🎯 Purpose</th>
</tr>
<tr>
<td><b>🔬 samtools</b></td>
<td><code>1.22</code></td>
<td>SAM/BAM file manipulation and analysis</td>
</tr>
<tr>
<td><b>🧮 bcftools</b></td>
<td><code>Library issue</code></td>
<td>VCF/BCF file manipulation and analysis</td>
</tr>
<tr>
<td><b>📊 vcftools</b></td>
<td><code>(0.1.17)</code></td>
<td>VCF file analysis and filtering</td>
</tr>
<tr>
<td><b>🏗️ bedtools</b></td>
<td><code>v2.31.1</code></td>
<td>Genomic interval manipulation</td>
</tr>
<tr>
<td><b>🧮 PLINK</b></td>
<td><code>1.9.0</code></td>
<td>Population genetics analysis</td>
</tr>
<tr>
<td><b>🧮 PLINK2</b></td>
<td><code>2.0.0</code></td>
<td>Next-generation PLINK</td>
</tr>
<tr>
<td><b>🎯 bwa</b></td>
<td><code>0.7.19-r1273</code></td>
<td>Short read alignment</td>
</tr>
<tr>
<td><b>📈 ANGSD</b></td>
<td><code>0.940</code></td>
<td>NGS data analysis</td>
</tr>
<tr>
<td><b>📇 tabix</b></td>
<td><code>1.22</code></td>
<td>Genomic file indexing</td>
</tr>
<tr>
<td><b>🌳 IQ-TREE</b></td>
<td><code>3.0.1</code></td>
<td>Phylogenetic inference and tree construction</td>
</tr>
</table>

---

## 🧬 **Population Genetics and Local Adaptation Tools**

<table>
<tr>
<th width="25%">🔧 Tool</th>
<th width="20%">📦 Version</th>
<th width="55%">🎯 Purpose</th>
</tr>
<tr>
<td><b>🌍 ADMIXTURE</b></td>
<td><code>1.3.0</code></td>
<td>Population structure analysis and ancestry analysis</td>
</tr>
<tr>
<td><b>🔍 BayeScan</b></td>
<td><code>2.1</code></td>
<td>FST outlier detection</td>
</tr>
<tr>
<td><b>📊 GEMMA</b></td>
<td><code>0.98.42021012920122021</code></td>
<td>Genome-wide association studies</td>
</tr>
<tr>
<td><b>🎯 SamBada</b></td>
<td><code>0.8.0</code></td>
<td>Local adaptation analysis and environmental association</td>
</tr>
<tr>
<td><b>🔄 BA3-SNPS</b></td>
<td><code>3.0.4</code></td>
<td>Migration analysis and gene flow estimation</td>
</tr>
</table>

---

## 🐍 **Python Packages**

<details>
<summary><b>📊 Core Scientific Computing</b> (Click to expand)</summary>

<table>
<tr>
<th width="30%">Package</th>
<th width="15%">Version</th>
<th width="55%">Purpose</th>
</tr>
<tr>
<td><b>numpy</b></td>
<td><code>1.26.4</code></td>
<td>Numerical computing foundation</td>
</tr>
<tr>
<td><b>pandas</b></td>
<td><code>2.3.1</code></td>
<td>Data manipulation and analysis</td>
</tr>
<tr>
<td><b>scipy</b></td>
<td><code>1.16.0</code></td>
<td>Scientific computing library</td>
</tr>
<tr>
<td><b>scikit-learn</b></td>
<td><code>1.7.1</code></td>
<td>Machine learning library</td>
</tr>
<tr>
<td><b>scikit-allel</b></td>
<td><code>1.3.13</code></td>
<td>Population genetics analysis</td>
</tr>
<tr>
<td><b>networkx</b></td>
<td><code>2.8.8</code></td>
<td>Network analysis</td>
</tr>
</table>

</details>

<details>
<summary><b>📈 Visualization</b> (Click to expand)</summary>

<table>
<tr>
<th width="30%">Package</th>
<th width="15%">Version</th>
<th width="55%">Purpose</th>
</tr>
<tr>
<td><b>matplotlib</b></td>
<td><code>3.10.3</code></td>
<td>Plotting and visualization</td>
</tr>
<tr>
<td><b>seaborn</b></td>
<td><code>0.13.2</code></td>
<td>Statistical data visualization</td>
</tr>
<tr>
<td><b>plotly</b></td>
<td><code>6.2.0</code></td>
<td>Interactive visualization</td>
</tr>
</table>

</details>

<details>
<summary><b>🧬 Bioinformatics Tools</b> (Click to expand)</summary>

<table>
<tr>
<th width="30%">Package</th>
<th width="15%">Version</th>
<th width="55%">Purpose</th>
</tr>
<tr>
<td><b>pysam</b></td>
<td><code>0.22.1</code></td>
<td>SAM/BAM file processing</td>
</tr>
<tr>
<td><b>biopython</b></td>
<td><code>1.85</code></td>
<td>Biological computation and utilities</td>
</tr>
<tr>
<td><b>cyvcf2</b></td>
<td><code>0.31.0</code></td>
<td>Fast VCF parsing</td>
</tr>
</table>

</details>

<details>
<summary><b>🌍 Geospatial Analysis</b> (Click to expand)</summary>

<table>
<tr>
<th width="30%">Package</th>
<th width="15%">Version</th>
<th width="55%">Purpose</th>
</tr>
<tr>
<td><b>geopandas</b></td>
<td><code>1.1.1</code></td>
<td>Geospatial data analysis</td>
</tr>
<tr>
<td><b>rasterio</b></td>
<td><code>1.3.10</code></td>
<td>Raster data I/O</td>
</tr>
</table>

</details>

---

## 📊 **R Packages**

<details>
<summary><b>🔧 Core Data Manipulation</b> (Click to expand)</summary>

<table>
<tr>
<th width="30%">Package</th>
<th width="15%">Version</th>
<th width="55%">Purpose</th>
</tr>
<tr>
<td><b>tidyverse</b></td>
<td><code>2.0.0</code></td>
<td>Data science meta-package toolkit</td>
</tr>
<tr>
<td><b>data.table</b></td>
<td><code>1.17.8</code></td>
<td>High-performance data manipulation</td>
</tr>
<tr>
<td><b>ggplot2</b></td>
<td><code>3.5.2</code></td>
<td>Grammar of graphics plotting</td>
</tr>
<tr>
<td><b>dplyr</b></td>
<td><code>1.1.4</code></td>
<td>Data manipulation verbs</td>
</tr>
<tr>
<td><b>devtools</b></td>
<td><code>2.4.5</code></td>
<td>Package development tools</td>
</tr>
<tr>
<td><b>remotes</b></td>
<td><code>2.5.0</code></td>
<td>Package installation from remote repositories</td>
</tr>
</table>

</details>

<details>
<summary><b>🧬 Population Genetics & Bioinformatics</b> (Click to expand)</summary>

<table>
<tr>
<th width="30%">Package</th>
<th width="15%">Version</th>
<th width="55%">Purpose</th>
</tr>
<tr>
<td><b>adegenet</b></td>
<td><code>2.1.11</code></td>
<td>Exploratory analysis of genetic data</td>
</tr>
<tr>
<td><b>vcfR</b></td>
<td><code>1.15.0</code></td>
<td>VCF file manipulation in R</td>
</tr>
<tr>
<td><b>SNPRelate</b></td>
<td><code>1.36.0</code></td>
<td>Parallel computing for SNP data</td>
</tr>
<tr>
<td><b>LEA</b></td>
<td><code>3.19.8</code></td>
<td>Landscape and ecological association studies</td>
</tr>
<tr>
<td><b>lfmm</b></td>
<td><code>1.0</code></td>
<td>Latent factor mixed models</td>
</tr>
<tr>
<td><b>OutFLANK</b></td>
<td><code>0.2</code></td>
<td>FST outlier detection</td>
</tr>
<tr>
<td><b>pcadapt</b></td>
<td><code>4.4.0</code></td>
<td>Principal component analysis for population genetics</td>
</tr>
<tr>
<td><b>lostruct</b></td>
<td><code>0.0.0.9000</code></td>
<td>Local PCA for population structure</td>
</tr>
<tr>
<td><b>admixtools</b></td>
<td><code>2.0.10</code></td>
<td>ADMIXTOOLS wrapper for R</td>
</tr>
<tr>
<td><b>admixr</b></td>
<td><code>0.10.0</code></td>
<td>ADMIXTOOLS analysis in R</td>
</tr>
</table>

</details>

<details>
<summary><b>🌍 Geospatial & Visualization</b> (Click to expand)</summary>

<table>
<tr>
<th width="30%">Package</th>
<th width="15%">Version</th>
<th width="55%">Purpose</th>
</tr>
<tr>
<td><b>sf</b></td>
<td><code>1.0.16</code></td>
<td>Simple features for spatial data</td>
</tr>
<tr>
<td><b>raster</b></td>
<td><code>3.6.32</code></td>
<td>Raster data analysis</td>
</tr>
<tr>
<td><b>terra</b></td>
<td><code>1.8.60</code></td>
<td>Spatial data analysis</td>
</tr>
<tr>
<td><b>leaflet</b></td>
<td><code>2.2.2</code></td>
<td>Interactive web maps</td>
</tr>
<tr>
<td><b>tmap</b></td>
<td><code>3.3.4</code></td>
<td>Thematic maps</td>
</tr>
<tr>
<td><b>ggstatsplot</b></td>
<td><code>0.13.1</code></td>
<td>Statistical plots with ggplot2</td>
</tr>
<tr>
<td><b>plotly</b></td>
<td><code>4.11.0</code></td>
<td>Interactive web graphics</td>
</tr>
</table>

</details>

---

## 🔧 **System Dependencies**

<details>
<summary><b>🖥️ Required System Libraries</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🌍 Geospatial Libraries**
- `GDAL` - Geospatial data abstraction
- `PROJ` - Coordinate reference systems
- `GEOS` - Geometry engine
- `SQLite3` - Database support

</td>
<td width="50%">

**🛠️ Development Tools**
- `GCC7` - C++ compiler (SamBada)
- `GSL` - GNU Scientific Library
- `libcurl4-openssl-dev` - HTTP library
- `libssl-dev` - SSL/TLS library
- `libxml2-dev` - XML parsing

</td>
</tr>
</table>

</details>

---

## 📦 **Installation Methods**

<div align="center">

### 🎯 **Package Installation Strategy**

</div>

<details>
<summary><b>🐍 Micromamba Packages</b> (Click to expand)</summary>

Most packages are installed via **micromamba** from conda-forge and bioconda channels for optimal performance and dependency resolution.

**Advantages:**
- ⚡ Faster than conda
- 🔧 Better dependency resolution
- 🐳 Container-optimized
- 🔄 Reproducible environments

</details>

<details>
<summary><b>🐙 GitHub Packages</b> (Click to expand)</summary>

Specialized packages installed via devtools/remotes:

<table>
<tr>
<th width="40%">📦 Package</th>
<th width="60%">🔗 GitHub Repository</th>
</tr>
<tr>
<td><b>R.SamBada</b></td>
<td><code>SolangeD/R.SamBada</code></td>
</tr>
<tr>
<td><b>gdalUtils</b></td>
<td><code>gearslaboratory/gdalUtils</code></td>
</tr>
<tr>
<td><b>lostruct</b></td>
<td><code>petrelharp/local_pca/lostruct</code></td>
</tr>
<tr>
<td><b>TESS3_encho_sen</b></td>
<td><code>eriqande/TESS3_encho_sen</code></td>
</tr>
<tr>
<td><b>admixtools</b></td>
<td><code>uqrmaie1/admixtools</code></td>
</tr>
<tr>
<td><b>lfmm</b></td>
<td><code>bcm-uga/lfmm</code></td>
</tr>
<tr>
<td><b>LEA</b></td>
<td><code>bcm-uga/LEA</code></td>
</tr>
</table>

</details>

---

## 📋 **System Information**

| Component | Details |
|-----------|---------|
| **Operating System** | Debian GNU/Linux 12 (bookworm) |
| **Kernel** | 4.18.0-477.27.1.el8_8.x86_64 |
| **Architecture** | x86_64 |
| **Container Engine** | micromamba 1.5.0 |
| **Shell** | /bin/bash |

---

## 📝 **Usage Notes**

This tool inventory was automatically generated and includes all software versions used in the Aedes aegypti local adaptation analysis pipeline. For reproducibility:

1. **Container Usage**: All analyses should be performed within the containerized environment to ensure identical software versions.

2. **Version Dependencies**: Some tools may have interdependent version requirements. The container environment resolves all compatibility issues.

3. **Citation Guidelines**: When publishing results, reference this tool inventory in supplementary materials and cite individual tools as appropriate.

4. **Updates**: This report reflects the software versions at the time of generation. Container updates may include newer versions.

## 🔄 **Reproducibility Statement**

All computational analyses were performed using the software versions documented above within a containerized environment. The complete computational environment can be reproduced using the Dockerfile and associated configuration files in this repository.

For questions about specific tool versions or installation procedures, please refer to the project documentation or contact the development team through the project repository.

---

<div align="center">

**🐳 Container Size:** ~3GB | **🏗️ Architecture:** AMD64 (HPC Optimized) | **📦 Total Packages:** 100+

*All packages are pre-installed and tested for compatibility*

</div>
# 🦟 Aedes Local Adaptation Analysis - Package Requirements

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

```bash
# Docker
docker run -it ghcr.io/cosmelab/aedes-local-adaptation:latest

# Singularity (HPC)
singularity exec aedes-local-adaptation.sif bash
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
<td><b>🐍 python</b></td>
<td><code>3.11.7</code></td>
<td>Python programming environment</td>
</tr>
<tr>
<td><b>📊 R</b></td>
<td><code>4.3.2+</code></td>
<td>Statistical computing (rocker/tidyverse base)</td>
</tr>
<tr>
<td><b>📦 micromamba</b></td>
<td><code>Latest</code></td>
<td>Fast package management (optimized for containers)</td>
</tr>
</table>

---

## 🧬 **Bioinformatics Tools**

<details>
<summary><b>🔬 Core Analysis Tools</b> (Click to expand)</summary>

<table>
<tr>
<th width="25%">🔧 Tool</th>
<th width="75%">🎯 Purpose</th>
</tr>
<tr>
<td><b>🧮 PLINK2</b></td>
<td>Genotype data manipulation and population genetics analysis</td>
</tr>
<tr>
<td><b>📊 bcftools</b></td>
<td>VCF file manipulation and variant calling</td>
</tr>
<tr>
<td><b>🏗️ fastStructure</b></td>
<td>Population structure analysis</td>
</tr>
<tr>
<td><b>🌍 Admixture</b></td>
<td>Ancestry analysis and admixture proportions</td>
</tr>
<tr>
<td><b>🌳 IQ-TREE</b></td>
<td>Phylogenetic analysis and tree construction</td>
</tr>
<tr>
<td><b>🔄 BA3-SNPs</b></td>
<td>Migration analysis and gene flow estimation</td>
</tr>
<tr>
<td><b>🎯 SamBada</b></td>
<td>Local adaptation analysis and environmental association</td>
</tr>
</table>

</details>

---

## 📊 **R Packages**

<details>
<summary><b>🔧 Core Data Manipulation</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**📋 Data Processing**
- `tidyverse` - Data science toolkit
- `data.table` - Fast data manipulation
- `here` - File path management
- `dplyr` - Data manipulation grammar
- `readr` - Data import/export
- `purrr` - Functional programming
- `stringr` - String manipulation

</td>
<td width="50%">

**🔍 Additional Tools**
- `forcats` - Factor handling
- `lubridate` - Date/time manipulation
- `janitor` - Data cleaning
- `broom` - Statistical model tidying
- `knitr` - Dynamic report generation
- `rmarkdown` - R Markdown documents

</td>
</tr>
</table>

</details>

<details>
<summary><b>📈 Visualization Packages</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🎨 Core Visualization**
- `ggplot2` - Grammar of graphics
- `scales` - Scale functions for ggplot2
- `ggrepel` - Non-overlapping text labels
- `ggtext` - Rich text rendering
- `ggvenn` - Venn diagrams
- `ggstatsplot` - Statistical plots
- `ggforce` - Extended ggplot2 functionality
- `ggpattern` - Pattern fills

</td>
<td width="50%">

**🎯 Specialized Plots**
- `scatterpie` - Scatter pie charts
- `RColorBrewer` - Color palettes
- `extrafont` - Font management
- `flextable` - Table formatting
- `officer` - Office document creation
- `Cairo` - Graphics device
- `ggspatial` - Spatial visualization

</td>
</tr>
</table>

</details>

<details>
<summary><b>🧬 Bioinformatics & Population Genetics</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🧬 Core Bioinformatics**
- `vcfR` - VCF file handling
- `adegenet` - Population genetics
- `dartR` - Genomic analysis
- `Biostrings` - DNA/RNA sequence analysis
- `SNPRelate` - SNP analysis
- `biomaRt` - Biomart database access

</td>
<td width="50%">

**🌍 Population Genetics**
- `OutFLANK` - Outlier detection
- `pcadapt` - PCA-based adaptation detection
- `R.SamBada` - Environmental association
- `LEA` - Landscape genomics
- `qvalue` - Multiple testing correction
- `pegas` - Population genetics analysis
- `ape` - Phylogenetic analysis
- `genetics` - Genetic data analysis

</td>
</tr>
</table>

</details>

<details>
<summary><b>🌍 Spatial Analysis</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🗺️ Spatial Data**
- `raster` - Raster data analysis
- `sf` - Simple features for R
- `rgdal` - Geospatial data abstraction
- `geosphere` - Spherical geometry
- `fields` - Spatial statistics

</td>
<td width="50%">

**🌐 Mapping & Visualization**
- `rnaturalearth` - Natural Earth data
- `rnaturalearthdata` - Natural Earth datasets
- `ggspatial` - Spatial ggplot2 extensions
- `contextily` - Contextual map tiles

</td>
</tr>
</table>

</details>

---

## 🐍 **Python Packages**

<details>
<summary><b>📊 Core Python Stack</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🔢 Scientific Computing**
- `numpy` - Numerical computing
- `scipy` - Scientific computing
- `pandas` - Data manipulation
- `scikit-learn` - Machine learning
- `scikit-allel` - Population genetics
- `networkx` - Network analysis

</td>
<td width="50%">

**📈 Visualization**
- `matplotlib` - Plotting library
- `seaborn` - Statistical visualization
- `plotly` - Interactive plots
- `bokeh` - Web-based visualization
- `altair` - Declarative visualization
- `holoviews` - Data visualization

</td>
</tr>
</table>

</details>

<details>
<summary><b>🧬 Bioinformatics Python</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🧬 Genomics Tools**
- `pysam` - SAM/BAM file processing
- `biopython` - Bioinformatics utilities
- `cyvcf2` - Fast VCF parsing
- `pyfaidx` - FASTA file indexing

</td>
<td width="50%">

**🌍 Geospatial Analysis**
- `geopandas` - Geospatial data analysis
- `fiona` - Vector data I/O
- `rasterio` - Raster data I/O
- `pyproj` - Cartographic projections
- `shapely` - Geometric operations
- `folium` - Interactive maps
- `earthpy` - Earth science data

</td>
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
<summary><b>📊 R Packages (CRAN)</b> (Click to expand)</summary>

Additional R packages installed via `install.packages()`:

```r
# Core packages
install.packages(c("data.table", "tidyverse", "qqman", "qqplotr",
                   "reticulate", "broom", "readxl", "writexl"))

# Bioinformatics packages
install.packages(c("pegas", "ape", "phangorn", "vcfR", "genetics"))

# Visualization packages
install.packages(c("ggrepel", "ggtext", "ggvenn", "ggstatsplot",
                   "RColorBrewer", "Cairo"))
```

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

<details>
<summary><b>🧬 Bioconductor Packages</b> (Click to expand)</summary>

Specialized bioinformatics packages:

```r
# Install BiocManager
install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("SNPRelate", "biomaRt", "LEA",
                       "AnnotationDbi", "VariantAnnotation",
                       "Biostrings"))
```

</details>

<details>
<summary><b>🔧 External Tools</b> (Click to expand)</summary>

<table>
<tr>
<th width="30%">🛠️ Tool</th>
<th width="70%">📥 Installation Method</th>
</tr>
<tr>
<td><b>SamBada</b></td>
<td>Via R.SamBada package (automated)</td>
</tr>
<tr>
<td><b>AdmixTools</b></td>
<td>Compiled from source</td>
</tr>
<tr>
<td><b>BayeScan</b></td>
<td>Binary download</td>
</tr>
<tr>
<td><b>GEMMA</b></td>
<td>Binary download</td>
</tr>
<tr>
<td><b>BayesAss3-SNPs</b></td>
<td>Binary download</td>
</tr>
<tr>
<td><b>BA3-SNPS-autotune</b></td>
<td>Python script</td>
</tr>
</table>

</details>

---

<div align="center">

**🐳 Container Size:** ~3GB | **🏗️ Architecture:** AMD64 (HPC Optimized) | **📦 Total Packages:** 100+

*All packages are pre-installed and tested for compatibility*

</div>

# 🦟 Aedes Local Adaptation Analysis

<div align="center">

<!-- Primary Badges -->
![Species](https://img.shields.io/badge/Species-Aedes%20aegypti-brightgreen?style=for-the-badge&logo=dna)
![Analysis](https://img.shields.io/badge/Analysis-Local%20Adaptation-blue?style=for-the-badge&logo=chart-line)
![Container](https://img.shields.io/badge/Container-Ready-orange?style=for-the-badge&logo=docker)

<!-- Secondary Badges -->
![Python](https://img.shields.io/badge/Python-3.11.7-blue?style=flat-square&logo=python)
![R](https://img.shields.io/badge/R-4.3.2-blue?style=flat-square&logo=r)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)
![Platform](https://img.shields.io/badge/Platform-Linux%20HPC-lightgrey?style=flat-square)

**🧬 Local Adaptation Analysis Pipeline for Aedes aegypti Populations**

*Population genomics approach to identify local adaptation signatures using containerized bioinformatics tools*

<!-- Container Registry Badges -->
[![Docker Hub](https://img.shields.io/badge/Docker%20Hub-cosmelab%2Faedes--local--adaptation-2496ED?style=flat&logo=docker&logoColor=white)](https://hub.docker.com/r/cosmelab/aedes-local-adaptation)
[![GHCR](https://img.shields.io/badge/GHCR-ghcr.io%2Fcosmelab%2Faedes--local--adaptation-181717?style=flat&logo=github&logoColor=white)](https://github.com/cosmelab/aedes-local-adaptation/pkgs/container/aedes-local-adaptation)

<!-- Quick Stats -->
**📊 Analysis Capabilities:** 40+ Tools | **🐳 Container Size:** ~3GB | **🏗️ Architecture:** AMD64 (HPC Optimized)

</div>

---

## 📋 **Table of Contents**

- [🎯 Project Overview](#-project-overview)
- [🚀 Quick Start](#-quick-start)
- [🏗️ Project Structure](#️-project-structure)
- [🐳 Container Usage](#-container-usage)
- [🛠️ Features & Tools](#️-features--tools)
- [🔗 Remote Development Setup](#-remote-development-setup)
- [📚 Documentation](#-documentation)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)

---

## 🎯 **Project Overview**

This repository contains the analysis pipeline for studying local adaptation in *Aedes aegypti* populations across different geographic regions using population genomics approaches.

### Research Context

- **Species**: *Aedes aegypti* (Yellow fever mosquito)
- **Research Goal**: Identify genetic signatures of local adaptation across geographic populations
- **Main Question**: What genetic variants and population structure patterns are associated with local adaptation?
- **Analysis Focus**: Population structure, selection detection, and environmental association
- **Methodology**: Population genomics with environmental correlation analysis

<div align="center">

### 📊 **Project Statistics**

<table>
<tr>
<td align="center" width="16.66%">

**🧬 Analysis Type**
<br>
`Local Adaptation`

</td>
<td align="center" width="16.66%">

**📈 Geographic Scope**
<br>
`Multiple Populations`

</td>
<td align="center" width="16.66%">

**🔬 Tools Available**
<br>
`40+`

</td>
<td align="center" width="16.66%">

**🐳 Container Size**
<br>
`~3GB`

</td>
<td align="center" width="16.66%">

**🏗️ Architecture**
<br>
`AMD64 (HPC Optimized)`

</td>
<td align="center" width="16.66%">

**📚 Documentation**
<br>
`Complete`

</td>
</tr>
</table>

</div>

## 🏗️ **Project Structure**

```
aedes-local-adaptation/
├── configs/              # Configuration files
│   └── analysis_config/  # Analysis-specific configs
├── scripts/
│   ├── analysis/         # Analysis pipelines
│   └── visualization/    # Plotting and visualization
├── containers/           # Custom analysis containers
├── Dockerfile           # Docker environment setup
├── setup.sh             # Environment setup script
├── start_jupyter.sh     # Jupyter lab launcher
└── check_packages.sh    # Package verification
```

## 👥 **For New Collaborators**

<div align="center">

### 🚀 **Quick Start (3 Steps)**

**New to the project?** Get up and running in minutes:

</div>

<details>
<summary><b>📋 Step-by-Step Setup Guide</b> (Click to expand)</summary>

<table>
<tr>
<td align="center" width="60">

### 1️⃣

</td>
<td>

**🔧 Clone Repository & Setup**
```zsh
git clone https://github.com/cosmelab/aedes-local-adaptation.git
cd aedes-local-adaptation
./setup.sh  # Creates directories + shows container instructions
```
*Creates project structure and shows HPC-specific commands*

</td>
</tr>
<tr>
<td align="center" width="60">

### 2️⃣

</td>
<td>

**🐳 Download Container** (Choose either registry)
```zsh
# 🐙 GitHub Container Registry:
singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest

# 🐳 Docker Hub:
singularity pull aedes-local-adaptation.sif docker://cosmelab/aedes-local-adaptation:latest
```
*Both registries work identically - no authentication needed*

</td>
</tr>
<tr>
<td align="center" width="60">

### 3️⃣

</td>
<td>

**🧪 Test Everything Works**
```zsh
singularity shell --cleanenv --bind $PWD:/proj aedes-local-adaptation.sif
# Inside container, test all tools:
cd /proj && zsh scripts/test_container_tools.sh
```
*Comprehensive tool testing*

</td>
</tr>
</table>

</details>

<div align="center">

**🎉 That's it!** You now have:

| Status | Component |
|--------|-----------|
| ✅ | Project directory structure |
| ✅ | All analysis tools in container |
| ✅ | Scripts and configurations |
| ✅ | Ready to start analysis |

</div>

### 📚 **Next Steps for Collaborators**

<table>
<tr>
<td width="33%">

**🖥️ HPC Users**
- See [README_HPC.md](README_HPC.md)
- SLURM job templates
- Resource recommendations

</td>
<td width="33%">

**📊 Analysis Guide**
- Check `scripts/` directory
- Analysis workflows
- Example pipelines

</td>
<td width="33%">

**📖 Documentation**
- [package_requirements.md](package_requirements.md)
- Tool information
- Environment details

</td>
</tr>
</table>

---

## 🐳 **Container Usage**

This project provides a comprehensive bioinformatics environment for local adaptation analysis through Docker and Singularity containers.

### Available Containers

The analysis environment is available from two registries (choose either):

#### GitHub Container Registry (GHCR)
```bash
# Pull from GHCR (GitHub integration)
docker pull ghcr.io/cosmelab/aedes-local-adaptation:latest
```

#### Docker Hub
```bash
# Pull from Docker Hub (widely supported)
docker pull cosmelab/aedes-local-adaptation:latest
```

**Both registries work great!** Choose based on your preference or institutional requirements.

### Singularity on HPC

For HPC systems without Docker, use Singularity with either registry:

```bash
# Load Singularity module (if required)
module load singularity-ce/3.9.3

# Choose either registry:
# From GHCR:
singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest
# OR from Docker Hub:
singularity pull aedes-local-adaptation.sif docker://cosmelab/aedes-local-adaptation:latest
```

**Both registries work identically** - no authentication needed for public repositories.

### 💾 **Cache Management (Important!)**

Singularity uses a cache directory that can quickly fill up your home directory quota. Here's how to manage it:

<details>
<summary><b>🚨 Quota Issues & Solutions</b> (Click to expand)</summary>

#### Check Your Quota
```zsh
# Check home directory quota
check_quota home

# Check cache size
du -sh ~/.singularity/cache 2>/dev/null || echo "No cache found"
```

#### Use Temporary Cache (Recommended)
```zsh
# Create temporary cache in current directory (usually has more space)
mkdir -p ./singularity_temp_cache
export SINGULARITY_CACHEDIR=$PWD/singularity_temp_cache

# Pull container with temporary cache
singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest

# Clean up temporary cache after successful pull
rm -rf ./singularity_temp_cache
unset SINGULARITY_CACHEDIR
```

#### One-Line Solution for Quota Issues
```zsh
# Remove old cache, pull with temp cache, clean up
rm -rf ~/.singularity/cache && mkdir -p ./singularity_temp_cache && export SINGULARITY_CACHEDIR=$PWD/singularity_temp_cache && singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest && rm -rf ./singularity_temp_cache && unset SINGULARITY_CACHEDIR
```

#### Clean Up Existing Cache
```zsh
# Remove all cached images (frees up space)
rm -rf ~/.singularity/cache

# Or clean specific cache
singularity cache clean
```

</details>

### Running the Container

#### Interactive Shell
```zsh
# Start interactive shell with project directory mounted
singularity shell \
    --cleanenv \
    --bind /path/to/your/project:/proj \
    aedes-local-adaptation.sif

# Once inside the container, navigate to the project directory
cd /proj
```

#### Execute Commands
```bash
# Run a single command
singularity exec \
    --cleanenv \
    --bind /path/to/your/project:/proj \
    aedes-local-adaptation.sif your-command

# Run Python script
singularity exec \
    --cleanenv \
    --bind /path/to/your/project:/proj \
    aedes-local-adaptation.sif python /proj/scripts/your_script.py
```

### SLURM Batch Script Examples

#### Example 1: Local Adaptation Analysis
```bash
#!/bin/bash
#SBATCH --job-name=local-adaptation
#SBATCH --output=local_adaptation_%j.out
#SBATCH --error=local_adaptation_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=epyc

# Load Singularity module (if required)
module load singularity-ce/3.9.3

# Set paths
PROJECT_DIR="/path/to/your/project"
CONTAINER="aedes-local-adaptation.sif"

# Run local adaptation analysis
singularity exec \
    --cleanenv \
    --bind ${PROJECT_DIR}:/proj \
    ${CONTAINER} \
    bash -c "cd /proj && python scripts/analysis/local_adaptation_pipeline.py"
```

#### Example 2: Population Structure Analysis
```bash
#!/bin/bash
#SBATCH --job-name=pop-structure
#SBATCH --output=pop_structure_%j.out
#SBATCH --error=pop_structure_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=epyc

module load singularity-ce/3.9.3

PROJECT_DIR="/path/to/your/project"
CONTAINER="aedes-local-adaptation.sif"

# Run population structure analysis
singularity exec \
    --cleanenv \
    --bind ${PROJECT_DIR}:/proj \
    ${CONTAINER} \
    bash -c "cd /proj && python scripts/analysis/population_structure.py"
```

## 🛠️ **Features & Tools**

<div align="center">

### 🔬 **40+ Bioinformatics Tools Included**

*Comprehensive analysis environment for population genomics and local adaptation studies*

</div>

<details>
<summary><b>🧬 Core Genomics Tools</b> (Click to expand)</summary>

<table>
<tr>
<th width="20%">🔧 Tool</th>
<th width="15%">📦 Version</th>
<th width="65%">🎯 Purpose</th>
</tr>
<tr>
<td><b>🧬 samtools</b></td>
<td><code>1.22</code></td>
<td>SAM/BAM file manipulation and processing</td>
</tr>
<tr>
<td><b>📊 bcftools</b></td>
<td><code>1.22</code></td>
<td>VCF/BCF file manipulation and variant calling</td>
</tr>
<tr>
<td><b>🔍 vcftools</b></td>
<td><code>0.1.17</code></td>
<td>VCF file processing and population genetics statistics</td>
</tr>
<tr>
<td><b>🧮 angsd</b></td>
<td><code>0.940</code></td>
<td>Genotype likelihood analysis for low-coverage data</td>
</tr>
<tr>
<td><b>📈 plink</b></td>
<td><code>v1.9.0-b.8</code></td>
<td>Population genetics analysis and association studies</td>
</tr>
<tr>
<td><b>🚀 plink2</b></td>
<td><code>v2.0.0-a.6.9LM</code></td>
<td>Modern PLINK for large datasets and genomic analysis</td>
</tr>
<tr>
<td><b>🌍 bayenv2</b></td>
<td><code>2.0</code></td>
<td>Environmental association analysis and local adaptation</td>
</tr>
</table>

</details>

<details>
<summary><b>🐍 Python Environment</b> (Click to expand)</summary>

**Python 3.11.7** with comprehensive bioinformatics packages:

<table>
<tr>
<td width="50%">

**📊 Data Analysis**
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `scipy` - Scientific computing
- `sklearn` - Machine learning

</td>
<td width="50%">

**📈 Visualization**
- `matplotlib` - Plotting
- `seaborn` - Statistical visualization
- `plotly` - Interactive plots
- `bokeh` - Web-based visualization

</td>
</tr>
<tr>
<td>

**🧬 Bioinformatics**
- `pysam` - SAM/BAM processing
- `biopython` - Sequence analysis
- `pybedtools` - Genomic intervals
- `pyvcf` - VCF file handling

</td>
<td>

**🔢 Statistics**
- `statsmodels` - Statistical modeling
- `pingouin` - Statistical analysis
- `lifelines` - Survival analysis
- `patsy` - Statistical formulas

</td>
</tr>
</table>

</details>
<details>
<summary><b>📊 R Environment</b> (Click to expand)</summary>

**R 4.4.2** with specialized population genetics packages:

<table>
<tr>
<td width="50%">

**🧬 Population Genetics**
- `pegas` - Population genetics analysis
- `ape` - Phylogenetic analysis
- `vcfR` - VCF file handling in R
- `genetics` - Genetic data analysis
- `OutFLANK` - Outlier detection
- `pcadapt` - PCA-based outlier detection

</td>
<td width="50%">

**📈 Visualization & Stats**
- `ggplot2` - Data visualization
- `qqman` - Manhattan plots
- `qqplotr` - Q-Q plots
- `vegan` - Community ecology
- `ade4` - Multivariate analysis
- `MASS` - Statistical functions

</td>
</tr>
<tr>
<td>

**🔧 Data Manipulation**
- `data.table` - Fast data processing
- `tidyverse` - Data science toolkit
- `readxl` / `writexl` - Excel files
- `broom` - Statistical model tidying

</td>
<td>

**🌍 Environmental Analysis**
- `R.SamBada` - Environmental association
- `seqinr` - Sequence analysis
- `reticulate` - Python integration
- `qqconf` - Statistical testing

</td>
</tr>
</table>

</details>

<details>
<summary><b>⚙️ Development Environment</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🔄 Workflow Management**
- `snakemake` (9.6.2) - Pipeline orchestration
- `jupyter` (4.4.4) - Interactive notebooks
- `jupyterlab` (4.4.4) - Advanced notebook interface

</td>
<td width="50%">

**🖥️ Shell Environment**
- `zsh` (5.9) - Advanced shell
- `starship` (1.23.0) - Cross-shell prompt
- `lsd` (1.1.5) - Modern ls replacement
- `colorls` (1.5.0) - Colored ls with git status

</td>
</tr>
</table>

</details>

## 🔬 **Technical Framework**

<div align="center">

### 🔄 **Analysis Workflow**

</div>

<details>
<summary><b>📊 Primary Analysis Pipeline</b> (Click to expand)</summary>

<table>
<tr>
<td align="center" width="33%">

### 🧬 **Population Structure**

**🎯 Purpose**
<br>
Identify genetic clusters and population structure

**🔧 Tools**
<br>
PCA, ADMIXTURE, STRUCTURE

**📈 Output**
<br>
Population assignments and genetic distances

</td>
<td align="center" width="33%">

### 🌍 **Local Adaptation**

**🎯 Purpose**
<br>
Identify signatures of local adaptation

**🔧 Tools**
<br>
FST outlier analysis, environmental association

**📈 Methods**
<br>
Bayenv2, LFMM, RDA analysis

</td>
<td align="center" width="33%">

### 🔍 **Selection Detection**

**🎯 Purpose**
<br>
Detect various types of selection

**🔧 Tools**
<br>
Tajima's D, iHS, XP-EHH

**📈 Methods**
<br>
ANGSD, selscan, BEAGLE, IQ-TREE

</td>
</tr>
</table>

</details>

<details>
<summary><b>⚙️ Detailed Analysis Methods</b> (Click to expand)</summary>

#### 🧬 **Population Structure Analysis**
- **PCA**: Principal component analysis for population structure
- **ADMIXTURE**: Ancestry inference and admixture proportions
- **STRUCTURE**: Bayesian clustering of individuals

#### 🌍 **Local Adaptation Detection**
- **FST Outlier Analysis**: Identify loci under selection
- **Environmental Association**: Correlate genetic variation with environmental variables
- **Bayenv2**: Bayesian environmental association testing

#### 🔍 **Selection Detection**
- **Tajima's D**: Neutrality tests (ANGSD)
- **iHS**: Integrated haplotype score (selscan)
- **XP-EHH**: Cross-population extended haplotype homozygosity (selscan)
- **BEAGLE**: Haplotype phasing
- **IQ-TREE**: Phylogenetic inference

</details>

## 📊 **Analytical Goals**

<div align="center">

### 🎯 **Research Objectives**

</div>

<details>
<summary><b>🧬 Population Structure Analysis</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🎯 Objective**
<br>
Identify genetic clusters and population structure

**🔧 Methods**
<br>
PCA, ADMIXTURE, STRUCTURE analysis

</td>
<td width="50%">

**📈 Expected Output**
<br>
Population assignments and genetic distances

**🔍 Applications**
<br>
Understanding population history and migration patterns

</td>
</tr>
</table>

</details>

<details>
<summary><b>🌍 Local Adaptation Detection</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🎯 FST Outlier Analysis**
<br>
Identify loci under selection

**🎯 Environmental Association**
<br>
Correlate genetic variation with environmental variables

</td>
<td width="50%">

**🔧 Bayenv2 Analysis**
<br>
Bayesian environmental association testing

**📈 Statistical Methods**
<br>
RDA, LFMM, Bayenv2

</td>
</tr>
</table>

</details>

<details>
<summary><b>🔍 Selection Detection</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🧮 Tajima's D**
<br>
Detect balancing and directional selection

**📊 iHS**
<br>
Identify recent positive selection within populations

</td>
<td width="50%">

**🔬 XP-EHH**
<br>
Detect selection differences between populations

**🌍 Environmental Variables**
<br>
Temperature, precipitation, elevation, urbanization

</td>
</tr>
</table>

</details>

## 🔗 **Remote Development Setup**

### SSH Configuration for HPC Access

To access your HPC system remotely via Cursor, VS Code, or other SSH clients:

#### 1. SSH Config Setup
Add this to your `~/.ssh/config` file:

```bash
# HPC SSH Configuration
Host ucr-hpc
  HostName cluster.hpcc.ucr.edu
  User your-username
  IdentityFile ~/.ssh/id_rsa
  AddKeysToAgent yes
  UseKeychain yes
  ControlMaster auto
  ControlPath ~/.ssh/cm-%r@%h:%p
  ControlPersist 600
```

#### 2. SSH Key Setup
```bash
# Generate SSH key (if you don't have one)
ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa

# Copy public key to HPC
ssh-copy-id -i ~/.ssh/id_rsa.pub your-username@cluster.hpcc.ucr.edu

# Test connection
ssh ucr-hpc
```

### Remote Development with Cursor/VS Code

#### 1. Connect to HPC
1. Open Cursor/VS Code
2. Press `Ctrl+Shift+P` (or `Cmd+Shift+P` on Mac)
3. Type "Remote-SSH: Connect to Host"
4. Select your configured host (e.g., `ucr-hpc`)
5. Choose platform (Linux)
6. Enter password if prompted

#### 2. Open Project Directory
1. Once connected, click "Open Folder"
2. Navigate to your project directory
3. Open the project

#### 3. Use Integrated Terminal
- Open terminal in Cursor/VS Code
- Run Singularity commands directly
- Access all HPC resources

### SSH Tunneling for Jupyter

```bash
# Create SSH tunnel for Jupyter
ssh -L 8888:localhost:8888 ucr-hpc

# On HPC, start Jupyter
singularity exec --cleanenv --bind /path/to/project:/proj aedes-local-adaptation.sif \
    bash -c "cd /proj && jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root"

# Open in browser: http://localhost:8888
```

## 📚 **Documentation**

- **HPC_Setup_Guide.md**: Instructions for HPC deployment
- **package_requirements.md**: Detailed package requirements
- **USER_RULES.md**: AI assistant interaction guidelines

## 🤝 **Contributing**

This project is designed for publication purposes. For questions or issues, please open a GitHub issue or contact the maintainers.

## 📄 **License**

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 **Acknowledgments**

- **cosmelab team** for project development
- **Bioconductor** for R package ecosystem
- **Conda-forge** for Python package management
- **Singularity** for HPC containerization

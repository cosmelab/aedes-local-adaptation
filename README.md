# 🦟 Aedes Local Adaptation Analysis

<div align="center">

![Aedes aegypti](https://img.shields.io/badge/Species-Aedes%20aegypti-brightgreen)
![Analysis Type](https://img.shields.io/badge/Analysis-Local%20Adaptation-blue)
![Container](https://img.shields.io/badge/Container-Docker%20%7C%20Singularity-orange)
![License](https://img.shields.io/badge/License-MIT-green)

**Local Adaptation Analysis Pipeline for Aedes aegypti Populations**

*Population genomics approach to identify local adaptation signatures*

[![Docker Hub](https://img.shields.io/badge/Docker%20Hub-cosmelab%2Faedes--local--adaptation-blue?logo=docker)](https://hub.docker.com/r/cosmelab/aedes-local-adaptation)
[![GitHub Container Registry](https://img.shields.io/badge/GHCR-ghcr.io%2Fcosmelab%2Faedes--local--adaptation-purple?logo=github)](https://github.com/cosmelab/aedes-local-adaptation/pkgs/container/aedes-local-adaptation)

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

| Metric | Value |
|--------|-------|
| 🧬 **Analysis Type** | Local Adaptation |
| 📈 **Geographic Scope** | Multiple Populations |
| 🔬 **Tools Available** | 40+ |
| 🐳 **Container Size** | ~2.3GB |
| 🏗️ **Architectures** | AMD64, ARM64 |
| 📚 **Documentation** | Complete |

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

## 🐳 **Container Usage**

This project provides a comprehensive bioinformatics environment for local adaptation analysis through Docker and Singularity containers.

### Available Containers

The analysis environment is available from two sources:

#### Docker Hub
```bash
# Pull from Docker Hub (public)
docker pull cosmelab/aedes-local-adaptation:latest
```

#### GitHub Container Registry (GHCR)
```bash
# Pull from GHCR (requires authentication for private repos)
docker pull ghcr.io/cosmelab/aedes-local-adaptation:latest
```

### Singularity on HPC

For HPC systems without Docker, use Singularity:

#### From Docker Hub
```bash
# Load Singularity module (if required)
module load singularity-ce/3.9.3

# Pull and build Singularity container
singularity pull aedes-local-adaptation.sif docker://cosmelab/aedes-local-adaptation:latest
```

#### From GHCR (requires authentication)
```bash
# Authenticate with GHCR (requires GitHub token)
echo $GITHUB_TOKEN | singularity remote login --username cosmelab --password-stdin oras://ghcr.io

# Pull and build Singularity container
singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest
```

### Running the Container

#### Interactive Shell
```bash
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

</div>

### 🧬 **Core Genomics Tools**

| Tool | Version | Purpose |
|------|---------|---------|
| 🧬 **samtools** | 1.22 | SAM/BAM file manipulation |
| 📊 **bcftools** | 1.22 | VCF/BCF file manipulation |
| 🔍 **vcftools** | 0.1.17 | VCF file processing |
| 🧮 **angsd** | 0.940 | Genotype likelihood analysis |
| 📈 **plink** | v1.9.0-b.8 | Population genetics analysis |
| 🚀 **plink2** | v2.0.0-a.6.9LM | Modern PLINK for large datasets |
| 🌍 **bayenv2** | 2.0 | Environmental association analysis |

### Python Environment
- **Python** (3.11.7) with comprehensive bioinformatics packages:
  - `pandas`, `numpy`, `scipy` - Data analysis
  - `matplotlib`, `seaborn`, `plotly` - Visualization
  - `sklearn` - Machine learning
  - `pysam` - SAM/BAM file processing
  - `Bio` (biopython) - Bioinformatics utilities
  - `allel` (scikit-allel) - Population genetics
  - `cyvcf2` - Fast VCF parsing
  - `pyfaidx` - FASTA file indexing

### R Environment
- **R** (4.4.2) with population genetics packages:
  - `ade4`, `MASS`, `ggplot2`, `vegan`, `seqinr`, `qqconf`
  - `data.table`, `tidyverse` - Data manipulation
  - `pegas`, `ape`, `vcfR`, `genetics` - Population genetics
  - `qqman`, `qqplotr` - Manhattan plots
  - `reticulate`, `broom`, `readxl`, `writexl` - Utilities

### Workflow Management
- **snakemake** (9.6.2) - Workflow management
- **jupyter** (4.4.4) - Interactive analysis
- **jupyterlab** (4.4.4) - JupyterLab interface

### Shell Environment
- **zsh** (5.9) with Oh-My-Zsh and Powerlevel10k
- **starship** (1.23.0) - Cross-shell prompt
- **lsd** (1.1.5) - Modern ls replacement
- **colorls** (1.5.0) - Colored ls with git status

## 🔬 **Technical Framework**

### Primary Analysis Pipeline

#### Population Structure Analysis
- **Purpose**: Identify population structure and genetic clusters
- **Tools**: PCA, ADMIXTURE, STRUCTURE
- **Output**: Population assignments and genetic distances

#### Local Adaptation Detection
- **Purpose**: Identify signatures of local adaptation
- **Tools**: FST outlier analysis, environmental association
- **Methods**: Bayenv2, LFMM, RDA analysis

#### Selection Detection
- **Tajima's D**: Neutrality tests (ANGSD)
- **iHS**: Integrated haplotype score (selscan)
- **XP-EHH**: Cross-population extended haplotype homozygosity (selscan)
- **BEAGLE**: Haplotype phasing
- **IQ-TREE**: Phylogenetic inference

## 📊 **Analytical Goals**

### Population Structure
- **Objective**: Identify genetic clusters and population structure
- **Method**: PCA, ADMIXTURE, STRUCTURE analysis
- **Output**: Population assignments and genetic distances

### Local Adaptation Detection
- **FST Outlier Analysis**: Identify loci under selection
- **Environmental Association**: Correlate genetic variation with environmental variables
- **Bayenv2 Analysis**: Bayesian environmental association testing

### Selection Detection
- **Tajima's D**: Detect balancing and directional selection
- **iHS**: Identify recent positive selection within populations
- **XP-EHH**: Detect selection differences between populations

### Environmental Correlation
- **Climate Variables**: Temperature, precipitation, elevation
- **Landscape Features**: Urbanization, land use patterns
- **Statistical Methods**: RDA, LFMM, Bayenv2

## 🚀 **Quick Start**

<div align="center">

### ⚡ **Get Started in 3 Steps**

</div>

### 1️⃣ **Clone & Setup**

```bash
# Clone repository
git clone https://github.com/cosmelab/aedes-local-adaptation.git
cd aedes-local-adaptation

# Pull Singularity container
singularity pull aedes-local-adaptation.sif docker://cosmelab/aedes-local-adaptation:latest
```

### 2️⃣ **Test Environment**

```bash
# Start interactive shell
singularity shell --cleanenv --bind /path/to/your/project:/proj aedes-local-adaptation.sif

# Test all tools
cd /proj
./check_packages.sh
```

### 3️⃣ **Start Analysis**

```bash
# Your local adaptation analysis pipeline is ready!
# Check the documentation below for detailed workflows
```

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

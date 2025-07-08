# HPC Setup Guide for Aedes Local Adaptation Analysis

## Overview

This guide provides step-by-step instructions for setting up and using the Aedes Local Adaptation analysis environment on High-Performance Computing (HPC) systems. The analysis pipeline includes population structure analysis, local adaptation detection, and genomic analysis tools for Aedes aegypti mosquitoes.

## Prerequisites

- Access to UCR HPCC or Yale HPC systems
- SSH access configured
- Basic familiarity with command line operations
- Docker installed locally (for building containers)

### Docker Installation and Setup

#### Install Docker Desktop

1. **macOS**: Download and install Docker Desktop from <https://www.docker.com/products/docker-desktop>
2. **Linux**: Follow instructions at <https://docs.docker.com/engine/install/>
3. **Windows**: Download Docker Desktop from <https://www.docker.com/products/docker-desktop>

#### Enable Docker Buildx

```bash
# Enable Buildx (one time setup)
docker buildx create --use
```

#### Create Output Directory

```bash
# Create directory for Docker images
mkdir -p ~/Dropbox/docker
```

## System-Specific Instructions

### UCR HPCC Setup

#### 1. Container Preparation (Local Machine)

```bash
# Clone the repository locally
git clone <repository-url>
cd AedesLocalAdaptation

# Build AMD64 image and load it locally
docker buildx build \
    --platform linux/amd64 \
    --load \
    --progress=auto \
    -t aedes-local-adaptation:amd64 \
    ~/Projects/AedesLocalAdaptation

# Save the AMD64 image tarball into Dropbox
docker save aedes-local-adaptation:amd64 \
    | gzip > ~/Dropbox/docker/aedes-local-adaptation-amd64.tar.gz
```

#### 2. Transfer to UCR HPCC

```bash
# Transfer the compressed tarball to UCR HPCC
scp ~/Dropbox/docker/aedes-local-adaptation-amd64.tar.gz \
    lcosme@cluster.hpcc.ucr.edu:/bigdata/cosmelab/lcosme/docker/

# Transfer your project files
scp -r ~/Projects/AedesLocalAdaptation/* \
    lcosme@cluster.hpcc.ucr.edu:/bigdata/cosmelab/lcosme/docker/
```

#### 3. Connect to UCR HPCC

```bash
# SSH to UCR HPCC
ssh username@cluster.hpcc.ucr.edu

# Request compute resources for interactive work
srun -p epyc -t 5:00:00 --pty -c 4 --mem=4g bash -l

# Load VS Code module
module load vscode

# Start VS Code tunnel
code tunnel
```

#### 4. Connect via VS Code/Cursor

1. Open VS Code or Cursor on your local machine
2. Use `Cmd+Shift+P` (Mac) or `Ctrl+Shift+P` (Windows/Linux)
3. Select "Connect to Tunnel"
4. Enter the tunnel URL provided by the `code tunnel` command
5. You will now have full access to the remote file system

#### 5. Run Analysis Container

```bash
# In the VS Code/Cursor terminal
singularity shell \
    --cleanenv \
    --bind /bigdata/cosmelab/lcosme/docker:/proj \
    /bigdata/cosmelab/lcosme/docker/AedesLocalAdaptation.sif

# Your container is now running with all analysis tools available
```

#### 3. Convert to Singularity Format (on UCR HPCC)

```bash
# SSH to UCR HPCC
ssh lcosme@cluster.hpcc.ucr.edu

# Navigate to your directory
cd /bigdata/cosmelab/lcosme/docker

# Load Singularity module
module load singularity-ce/3.9.3

# Convert Docker image to Singularity format
singularity build AedesLocalAdaptation.sif docker-archive://aedes-local-adaptation-amd64.tar.gz

# If you run out of space in your home directory
singularity cache clean -f -T all

# Use your scratch space for cache & temp
export SINGULARITY_CACHEDIR=/tmp/${USER}/.singularity_cache
export SINGULARITY_TMPDIR=/tmp/${USER}/.singularity_tmp
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

# Build with scratch space (if needed)
singularity build --quiet --tmpdir "$SINGULARITY_TMPDIR" \
    AedesLocalAdaptation.sif \
    docker-archive://aedes-local-adaptation-amd64.tar.gz

# Clean up the tar file to save space
rm aedes-local-adaptation-amd64.tar.gz
```

### Yale HPC Setup

#### 1. Container Preparation (Local Machine)

```bash
# Follow the same steps as UCR HPCC for building
docker buildx build \
    --platform linux/amd64 \
    --load \
    --progress=auto \
    -t aedes-local-adaptation:amd64 \
    ~/Projects/AedesLocalAdaptation

# Save the AMD64 image tarball into Dropbox
docker save aedes-local-adaptation:amd64 \
    | gzip > ~/Dropbox/docker/aedes-local-adaptation-amd64.tar.gz
```

#### 2. Transfer to Yale HPC

```bash
# Transfer to Yale HPC
scp ~/Dropbox/docker/aedes-local-adaptation-amd64.tar.gz \
    lvc26@mccleary.ycrc.yale.edu:/ycga-gpfs/project/caccone/lvc26/docker/

# Transfer project files
scp -r ~/Projects/AedesLocalAdaptation/* \
    lvc26@mccleary.ycrc.yale.edu:/ycga-gpfs/project/caccone/lvc26/docker/
```

#### 3. Convert to Apptainer Format (on Yale HPC)

```bash
# SSH to Yale HPC
ssh lvc26@mccleary.ycrc.yale.edu

# Navigate to your directory
cd /ycga-gpfs/project/caccone/lvc26/docker

# Convert Docker image to Apptainer format
apptainer build AedesLocalAdaptation.sif docker-archive://aedes-local-adaptation-amd64.tar.gz

# Clean up the tar file to save space
rm aedes-local-adaptation-amd64.tar.gz
```

#### 4. Access via Web Portal

1. Navigate to Yale's Open OnDemand portal
2. Log in with your Yale credentials
3. Launch VS Code through the web interface
4. You will have full file system access in the browser

#### 5. Run Analysis Container

```bash
# Request compute resources (must be on compute nodes, not login nodes)
srun --partition=general --cpus-per-task=8 --mem=16G --pty bash -l

# Run your container
apptainer shell \
    --cleanenv \
    --bind /ycga-gpfs/project/caccone/lvc26/docker:/proj \
    AedesLocalAdaptation.sif
```

## Available Analysis Tools

The container includes the following tools and packages:

### R Packages

- **Population Genetics**: `adegenet`, `dartR`, `OutFLANK`, `pcadapt`
- **Local Adaptation**: `R.SamBada`, `LEA`, `TESS3`, `lfmm`
- **Visualization**: `ggplot2`, `ggrepel`, `ggstatsplot`, `RColorBrewer`
- **Geospatial**: `sf`, `raster`, `rgdal`, `rnaturalearth`
- **Admixture Analysis**: `admixtools`, `admixr`

### Python Packages

- **Bioinformatics**: `scikit-allel`
- **Statistics**: `limix`
- **Visualization**: `matplotlib`, `seaborn`, `bokeh`
- **Geospatial**: `geopandas`, `folium`, `rasterio`

### Command Line Tools

- **Population Structure**: `fastStructure`, `Admixture`
- **Phylogenetics**: `IQ-TREE`
- **Selection Analysis**: `BayeScan`, `GEMMA`
- **Migration**: `BA3-SNPS`, `BA3-SNPS-autotune`
- **Variant Analysis**: `PLINK`, `bcftools`

## Analysis Workflow

### 1. Data Organization

```bash
# Navigate to your project directory
cd /proj

# Create analysis directories
mkdir -p data/{raw,processed,metadata}
mkdir -p output/{populations,local_adaptation,structure}
mkdir -p scripts/{analysis,visualization}
```

### 2. Quality Control

```bash
# Run quality control scripts
Rscript scripts/analysis/01_quality_control.R
```

### 3. Population Structure Analysis

```bash
# Run fastStructure
faststructure -K 3 -i data/processed/genotypes.txt -o output/structure/

# Run Admixture
admixture data/processed/genotypes.bed 3
```

### 4. Local Adaptation Analysis

```bash
# Run OutFLANK
Rscript scripts/analysis/outflank_analysis.R

# Run pcadapt
Rscript scripts/analysis/pcadapt_analysis.R

# Run R.SamBada
Rscript scripts/analysis/sambada_analysis.R
```

### 5. Visualization

```bash
# Generate plots
Rscript scripts/visualization/create_plots.R
```

## Job Submission (Non-Interactive)

For long-running analyses, submit jobs to the queue:

### UCR HPCC Job Script Example

```bash
#!/bin/bash
#SBATCH --job-name=aedes_analysis
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=aedes_analysis_%j.out
#SBATCH --error=aedes_analysis_%j.err

# Load modules
module load singularity-ce/3.9.3

# Run analysis
singularity exec \
    --cleanenv \
    --bind /bigdata/cosmelab/lcosme/docker:/proj \
    AedesLocalAdaptation.sif \
    Rscript scripts/analysis/main_analysis.R
```

### Yale HPC Job Script Example

```bash
#!/bin/bash
#SBATCH --job-name=aedes_analysis
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=aedes_analysis_%j.out
#SBATCH --error=aedes_analysis_%j.err

# Run analysis
apptainer exec \
    --cleanenv \
    --bind /ycga-gpfs/project/caccone/lvc26/docker:/proj \
    AedesLocalAdaptation.sif \
    Rscript scripts/analysis/main_analysis.R
```

Submit jobs with:

```bash
sbatch job_script.sh
```

## Troubleshooting

### Common Issues

1. **Container won't start**: Ensure you're on compute nodes, not login nodes
2. **Permission denied**: Check file permissions and ownership
3. **Memory issues**: Request more memory in your job submission
4. **Connection timeout**: Use VS Code tunnels instead of Remote-SSH

### Getting Help

- **UCR HPCC**: Check documentation at <https://hpcc.ucr.edu/manuals/hpc_cluster>
- **Yale HPC**: Visit <https://docs.ycrc.yale.edu>
- **Container issues**: Check Singularity/Apptainer documentation

## Best Practices

1. **Data Management**: Use scratch space for temporary files
2. **Resource Requests**: Request appropriate CPU and memory for your analyses
3. **File Organization**: Keep your project structure organized
4. **Version Control**: Use Git for tracking changes to scripts
5. **Backup**: Regularly backup important results and scripts

## Support

For technical support:

- **UCR HPCC**: Submit tickets through the HPCC portal
- **Yale HPC**: Contact YCRC support team
- **Analysis Questions**: Consult with your lab's bioinformatics team

---

*This guide is maintained by the Aedes Local Adaptation Analysis Team. Last updated: [Current Date]*

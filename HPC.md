# 🖥️ HPC Setup Guide - Aedes Local Adaptation Analysis

<div align="center">

![HPC](https://img.shields.io/badge/HPC-Singularity%20%7C%20SLURM-red?style=for-the-badge&logo=linux)
![Architecture](https://img.shields.io/badge/Architecture-AMD64%20HPC-blue?style=for-the-badge&logo=amd)
![Container](https://img.shields.io/badge/Container-3.7GB-orange?style=for-the-badge&logo=docker)
![Package Manager](https://img.shields.io/badge/Package%20Manager-micromamba-green?style=for-the-badge&logo=anaconda)

**🚀 Complete guide for running Aedes Local Adaptation analysis on HPC systems**

*Optimized for SLURM-based HPC clusters with Singularity/Apptainer*

</div>

---

## 🚀 **Quick Setup**

<details>
<summary><b>📥 1. Download Container</b> (Click to expand)</summary>

**Choose either container registry** (both work identically):

<table>
<tr>
<td width="50%">

**🐙 GitHub Container Registry (GHCR)**
```bash
# Load Singularity module (if required)
module load singularity-ce/3.9.3

# Download from GHCR (recommended)
singularity pull aedes-local-adaptation.sif \
    docker://ghcr.io/cosmelab/aedes-local-adaptation:latest
```

</td>
<td width="50%">

**🐳 Docker Hub**
```bash
# Load Singularity module (if required)
module load singularity-ce/3.9.3

# Download from Docker Hub
singularity pull aedes-local-adaptation.sif \
    docker://cosmelab/aedes-local-adaptation:latest
```

</td>
</tr>
</table>

### 💾 **Cache Management First! (Critical for HPC)**

**⚠️ IMPORTANT:** Before pulling containers, set up proper cache management to avoid quota issues:

```bash
# Check your quota status
check_quota home
check_quota bigdata

# Set up temporary cache (choose one method):

# Method 1: Use current directory (recommended)
mkdir -p ./singularity_temp_cache
export SINGULARITY_CACHEDIR=$PWD/singularity_temp_cache

# Method 2: Use scratch space (if available)
export SINGULARITY_CACHEDIR=/scratch/$USER/.singularity_cache
mkdir -p "$SINGULARITY_CACHEDIR"

# Pull container
singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest

# Clean up afterward
rm -rf ./singularity_temp_cache
unset SINGULARITY_CACHEDIR
```

**Emergency quota cleanup:**
```bash
# If you're over quota and can't pull containers:
rm -rf ~/.singularity/cache
```

</details>

<details>
<summary><b>🧪 2. Quick Test</b> (Click to expand)</summary>

### Basic Tool Testing
```bash
# Test core tools
singularity exec aedes-local-adaptation.sif python3 --version
singularity exec aedes-local-adaptation.sif R --version
singularity exec aedes-local-adaptation.sif micromamba --version
singularity exec aedes-local-adaptation.sif plink2 --version
```

### Comprehensive Tool Test
```bash
# Test all 50+ tools comprehensively
singularity shell --cleanenv --bind $PWD:/proj aedes-local-adaptation.sif
cd /proj && bash scripts/test_all_tools.sh
```

**🎯 The test script validates all tools including:**
- Python packages (numpy, pandas, scikit-allel, geopandas)
- R packages (adegenet, pcadapt, OutFLANK, lostruct)
- Bioinformatics tools (PLINK, bcftools, ADMIXTURE)
- Local adaptation tools (BayeScan, GEMMA, AdmixTools)

</details>

<details>
<summary><b>🖥️ 3. Interactive Session</b> (Click to expand)</summary>

### Start Interactive Session
```bash
# Request interactive resources
srun --partition=general --cpus-per-task=4 --mem=8G --time=2:00:00 --pty bash

# Enter container
singularity shell --cleanenv --bind $PWD:/proj aedes-local-adaptation.sif

# Activate micromamba environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate base
```

### Verify Environment
```bash
# Check package manager
micromamba list | head -10

# Test key tools
python3 -c "import numpy, pandas, allel; print('Python packages OK')"
Rscript -e "library(adegenet); library(pcadapt); cat('R packages OK\n')"
```

</details>

---

## 📋 **SLURM Job Templates**

<details>
<summary><b>🌍 Local Adaptation Analysis</b> (Click to expand)</summary>

```bash
#!/bin/bash
#SBATCH --job-name=local-adapt
#SBATCH --output=logs/local_adapt_%j.out
#SBATCH --error=logs/local_adapt_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=epyc

# Load required modules
module load singularity-ce/3.9.3

# Set paths
PROJECT_DIR=$PWD
CONTAINER="aedes-local-adaptation.sif"

# Create logs directory
mkdir -p logs

# Run OutFLANK analysis
singularity exec \
    --cleanenv \
    --bind ${PROJECT_DIR}:/proj \
    ${CONTAINER} \
    Rscript /proj/scripts/analysis/outflank_analysis.R

# Run pcadapt analysis
singularity exec \
    --cleanenv \
    --bind ${PROJECT_DIR}:/proj \
    ${CONTAINER} \
    Rscript /proj/scripts/analysis/pcadapt_analysis.R

# Run lostruct (local PCA) analysis
singularity exec \
    --cleanenv \
    --bind ${PROJECT_DIR}:/proj \
    ${CONTAINER} \
    Rscript /proj/scripts/analysis/lostruct_analysis.R

echo "Local adaptation analysis completed: $(date)"
```

**📊 Resource Requirements:** 16 CPUs, 64GB RAM, 24 hours

</details>

<details>
<summary><b>🧬 Population Structure Analysis</b> (Click to expand)</summary>

```bash
#!/bin/bash
#SBATCH --job-name=pop-struct
#SBATCH --output=logs/pop_struct_%j.out
#SBATCH --error=logs/pop_struct_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=epyc

module load singularity-ce/3.9.3

PROJECT_DIR=$PWD
CONTAINER="aedes-local-adaptation.sif"

# Run ADMIXTURE analysis
singularity exec \
    --cleanenv \
    --bind ${PROJECT_DIR}:/proj \
    ${CONTAINER} \
    bash -c "cd /proj && admixture --cv data/processed/filtered.bed 2"

# Run PCA analysis
singularity exec \
    --cleanenv \
    --bind ${PROJECT_DIR}:/proj \
    ${CONTAINER} \
    python3 /proj/scripts/analysis/pca_analysis.py

echo "Population structure analysis completed: $(date)"
```

**📊 Resource Requirements:** 8 CPUs, 32GB RAM, 12 hours

</details>

<details>
<summary><b>🌳 Phylogenetic Analysis</b> (Click to expand)</summary>

```bash
#!/bin/bash
#SBATCH --job-name=phylo
#SBATCH --output=logs/phylo_%j.out
#SBATCH --error=logs/phylo_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --partition=epyc

module load singularity-ce/3.9.3

PROJECT_DIR=$PWD
CONTAINER="aedes-local-adaptation.sif"

# Run IQ-TREE analysis
singularity exec \
    --cleanenv \
    --bind ${PROJECT_DIR}:/proj \
    ${CONTAINER} \
    iqtree -s /proj/data/processed/alignment.phy -m MFP -bb 1000 -nt AUTO

echo "Phylogenetic analysis completed: $(date)"
```

**📊 Resource Requirements:** 20 CPUs, 128GB RAM, 48 hours

</details>

<details>
<summary><b>📊 Interactive Analysis with Jupyter Lab</b> (Click to expand)</summary>

```bash
#!/bin/bash
#SBATCH --job-name=jupyter
#SBATCH --output=logs/jupyter_%j.out
#SBATCH --error=logs/jupyter_%j.err
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=general

module load singularity-ce/3.9.3

PROJECT_DIR=$PWD
CONTAINER="aedes-local-adaptation.sif"

# Start Jupyter Lab in container
singularity exec \
    --cleanenv \
    --bind ${PROJECT_DIR}:/proj \
    ${CONTAINER} \
    bash /proj/start_jupyter.sh

# Access at: http://your-hpc-node:8888
```

**📊 Resource Requirements:** 4 CPUs, 16GB RAM, 8 hours

</details>

---

## ⚡ **Analysis Workflow**

<details>
<summary><b>📁 1. Data Organization</b> (Click to expand)</summary>

```bash
# Navigate to project directory
cd /proj

# Create analysis directories (done by setup.sh)
mkdir -p data/{raw,processed,metadata}
mkdir -p results/{populations,local_adaptation,structure}
mkdir -p scripts/{analysis,visualization}
mkdir -p logs

# Verify directory structure
./setup.sh --validate
```

</details>

<details>
<summary><b>🔍 2. Quality Control</b> (Click to expand)</summary>

```bash
# Filter VCF files
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    bcftools view -q 0.05 -Q 0.95 data/raw/variants.vcf.gz > data/processed/filtered.vcf

# Convert to PLINK format
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    plink2 --vcf data/processed/filtered.vcf --make-bed --out data/processed/filtered
```

</details>

<details>
<summary><b>🧬 3. Population Structure Analysis</b> (Click to expand)</summary>

```bash
# Run ADMIXTURE
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    admixture --cv data/processed/filtered.bed 3

# Run PCA with R
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    Rscript scripts/analysis/pca_analysis.R
```

</details>

<details>
<summary><b>🌍 4. Local Adaptation Analysis</b> (Click to expand)</summary>

```bash
# Run OutFLANK
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    Rscript scripts/analysis/outflank_analysis.R

# Run pcadapt
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    Rscript scripts/analysis/pcadapt_analysis.R

# Run lostruct (local PCA)
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    Rscript scripts/analysis/lostruct_analysis.R
```

</details>

<details>
<summary><b>📊 5. Visualization</b> (Click to expand)</summary>

```bash
# Generate comprehensive plots
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    Rscript scripts/visualization/create_plots.R

# Interactive analysis with Jupyter
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    bash start_jupyter.sh
```

</details>

---

## 📊 **Resource Recommendations**

<div align="center">

### 🎯 **Job Resources by Analysis Type**

</div>

<table>
<tr>
<th width="25%">🔬 Analysis Type</th>
<th width="15%">💻 CPUs</th>
<th width="15%">🧠 Memory</th>
<th width="15%">⏱️ Time</th>
<th width="15%">🏗️ Partition</th>
<th width="15%">📊 Priority</th>
</tr>
<tr>
<td><b>📁 VCF Processing</b></td>
<td><code>4-8</code></td>
<td><code>16-32GB</code></td>
<td><code>2-6h</code></td>
<td><code>general</code></td>
<td>🟢 Low</td>
</tr>
<tr>
<td><b>🧬 Population Structure</b></td>
<td><code>8-16</code></td>
<td><code>32-64GB</code></td>
<td><code>6-12h</code></td>
<td><code>epyc</code></td>
<td>🟡 Medium</td>
</tr>
<tr>
<td><b>🌍 Local Adaptation</b></td>
<td><code>16-32</code></td>
<td><code>64-128GB</code></td>
<td><code>12-24h</code></td>
<td><code>epyc</code></td>
<td>🔴 High</td>
</tr>
<tr>
<td><b>🌳 Phylogenetics</b></td>
<td><code>20-40</code></td>
<td><code>128-256GB</code></td>
<td><code>24-48h</code></td>
<td><code>epyc</code></td>
<td>🔴 High</td>
</tr>
<tr>
<td><b>🌍 Environmental Association</b></td>
<td><code>8-16</code></td>
<td><code>32-64GB</code></td>
<td><code>6-12h</code></td>
<td><code>general</code></td>
<td>🟡 Medium</td>
</tr>
</table>

<div align="center">

### 💾 **Storage Requirements**

<table>
<tr>
<td align="center" width="20%">

**🐳 Container**
<br>
`~3.7GB`

</td>
<td align="center" width="20%">

**📁 Input VCF**
<br>
`1-10GB`

</td>
<td align="center" width="20%">

**🔄 Intermediate**
<br>
`2-5x input`

</td>
<td align="center" width="20%">

**📊 Results**
<br>
`100MB-1GB`

</td>
<td align="center" width="20%">

**💽 Total Recommended**
<br>
`20-50GB`

</td>
</tr>
</table>

</div>

---

## 🔧 **Common Commands**

<details>
<summary><b>📁 Data Processing</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🔍 VCF Filtering**
```bash
# Filter by MAF and missing data
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    bcftools view -q 0.05 -Q 0.95 \
    data/raw/variants.vcf.gz \
    > data/processed/filtered.vcf
```

</td>
<td width="50%">

**📊 Format Conversion**
```bash
# Convert VCF to PLINK
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    plink2 --vcf data/processed/filtered.vcf \
    --make-bed --out data/processed/filtered
```

</td>
</tr>
</table>

</details>

<details>
<summary><b>🧬 Population Analysis</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🔄 ADMIXTURE**
```bash
# Run ADMIXTURE with cross-validation
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    admixture --cv data/processed/filtered.bed 3
```

</td>
<td width="50%">

**📈 PCA Analysis**
```bash
# Principal component analysis
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    Rscript scripts/analysis/pca_analysis.R
```

</td>
</tr>
</table>

</details>

<details>
<summary><b>🌍 Local Adaptation</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🎯 OutFLANK**
```bash
# Outlier detection
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    Rscript scripts/analysis/outflank_analysis.R
```

</td>
<td width="50%">

**📊 pcadapt**
```bash
# PC-based adaptation detection
singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    Rscript scripts/analysis/pcadapt_analysis.R
```

</td>
</tr>
</table>

</details>

---

## 🔧 **Performance Optimization**

<details>
<summary><b>⚙️ HPC-Specific Settings</b> (Click to expand)</summary>

### Thread Control
```bash
# Set optimal thread counts (in your job script)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Use fast local storage for temporary files
export TMPDIR=/tmp/$USER
mkdir -p $TMPDIR
```

### Container Optimization
```bash
# Use --cleanenv for reproducible environments
singularity exec --cleanenv --bind $PWD:/proj aedes-local-adaptation.sif command

# For GPU nodes (if available)
singularity exec --nv --bind $PWD:/proj aedes-local-adaptation.sif command
```

</details>

---

## 🛠️ **Troubleshooting**

<details>
<summary><b>⚠️ Common Issues</b> (Click to expand)</summary>

<table>
<tr>
<th width="30%">🚨 Issue</th>
<th width="70%">🔧 Solution</th>
</tr>
<tr>
<td><b>Container won't start</b></td>
<td>Ensure you're on compute nodes, not login nodes. Use <code>srun</code> first.</td>
</tr>
<tr>
<td><b>Permission denied</b></td>
<td>Check file permissions and ownership. Use <code>--cleanenv</code> flag.</td>
</tr>
<tr>
<td><b>Memory issues</b></td>
<td>Request more memory in your job submission. Check resource limits.</td>
</tr>
<tr>
<td><b>micromamba not found</b></td>
<td>Use <code>eval "$(micromamba shell hook --shell bash)"</code> in container.</td>
</tr>
<tr>
<td><b>Tool not found</b></td>
<td>Run <code>bash scripts/test_all_tools.sh</code> to verify installation.</td>
</tr>
<tr>
<td><b>Quota exceeded</b></td>
<td>Clean up <code>~/.singularity/cache</code> and use temporary cache directories.</td>
</tr>
</table>

</details>

<details>
<summary><b>🔍 Debugging Commands</b> (Click to expand)</summary>

```bash
# Check container environment
singularity exec aedes-local-adaptation.sif env | grep -E "(PATH|CONDA|MAMBA)"

# Test micromamba
singularity exec aedes-local-adaptation.sif micromamba list

# Check tool availability
singularity exec aedes-local-adaptation.sif bash scripts/test_all_tools.sh

# Verify R packages
singularity exec aedes-local-adaptation.sif Rscript -e "installed.packages()[,1]"

# Check Python packages
singularity exec aedes-local-adaptation.sif python3 -c "import pkg_resources; [print(d) for d in pkg_resources.working_set]"
```

</details>

---

## 🔗 **Additional Resources**

<details>
<summary><b>🐳 Container Registries</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🐙 GitHub Container Registry**
<br>
[ghcr.io/cosmelab/aedes-local-adaptation](https://ghcr.io/cosmelab/aedes-local-adaptation)
<br>
*Integrated with GitHub, fast pulls, automatic builds*

</td>
<td width="50%">

**🐳 Docker Hub**
<br>
[cosmelab/aedes-local-adaptation](https://hub.docker.com/r/cosmelab/aedes-local-adaptation)
<br>
*Widely supported, reliable, good for institutional use*

</td>
</tr>
</table>

</details>

<details>
<summary><b>📚 Documentation</b> (Click to expand)</summary>

<table>
<tr>
<td width="33%">

**📖 Project Documentation**
- [Main README](README.md)
- [Tool Versions](TOOLS.md)
- [Test Scripts](scripts/test_all_tools.sh)

</td>
<td width="33%">

**🐳 Container Resources**
- [GitHub Container Registry](https://ghcr.io/cosmelab/aedes-local-adaptation)
- [Docker Hub](https://hub.docker.com/r/cosmelab/aedes-local-adaptation)
- [Tool Testing Script](scripts/test_all_tools.sh)

</td>
<td width="33%">

**🖥️ HPC Resources**
- [UCR HPCC Manual](https://hpcc.ucr.edu/manuals/hpc_cluster)
- [Yale HPC Docs](https://docs.ycrc.yale.edu)
- [Singularity Documentation](https://sylabs.io/docs/)

</td>
</tr>
</table>

</details>

---

<div align="center">

**🦟 Happy analyzing on HPC! 🖥️**

<table>
<tr>
<td align="center" width="25%">

**🐳 Container Size**
<br>
`~3.7GB`

</td>
<td align="center" width="25%">

**🏗️ Architecture**
<br>
`AMD64 (HPC Optimized)`

</td>
<td align="center" width="25%">

**🔧 Tools Available**
<br>
`50+`

</td>
<td align="center" width="25%">

**📊 Testing Status**
<br>
`Fully Validated`

</td>
</tr>
</table>

*For comprehensive tool testing, run: `bash scripts/test_all_tools.sh`*

</div>
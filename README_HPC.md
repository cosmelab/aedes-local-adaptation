# 🖥️ HPC Quick Start Guide - Aedes Local Adaptation Analysis

<div align="center">

![HPC](https://img.shields.io/badge/HPC-Singularity%20%7C%20SLURM-red?style=for-the-badge&logo=linux)
![Architecture](https://img.shields.io/badge/Architecture-AMD64%20HPC-blue?style=for-the-badge&logo=amd)
![Container](https://img.shields.io/badge/Container-3GB-orange?style=for-the-badge&logo=docker)

**🚀 Quick reference for running Aedes Local Adaptation analysis on HPC systems**

*Optimized for SLURM-based HPC clusters with Singularity/Apptainer*

</div>

---

## 🚀 **Quick Setup**

<details>
<summary><b>📥 1. Download Container</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🐙 GitHub Container Registry (GHCR)**
```bash
# Load Singularity module (system-dependent)
module load singularity-ce/3.9.3

# Download from GHCR
singularity pull aedes-local-adaptation.sif \
    docker://ghcr.io/cosmelab/aedes-local-adaptation:latest
```

</td>
<td width="50%">

**🐳 Docker Hub**
```bash
# Load Singularity module (system-dependent)
module load singularity-ce/3.9.3

# Download from Docker Hub
singularity pull aedes-local-adaptation.sif \
    docker://cosmelab/aedes-local-adaptation:latest
```

</td>
</tr>
</table>

**💡 Both registries work identically - choose based on your preference**

</details>

<details>
<summary><b>🧪 2. Quick Test</b> (Click to expand)</summary>

### Basic Tool Testing
```bash
# Test container and check tools
singularity exec aedes-local-adaptation.sif python --version
singularity exec aedes-local-adaptation.sif R --version
singularity exec aedes-local-adaptation.sif plink2 --version
```

### Comprehensive Tool Test
```bash
# Comprehensive tool testing
singularity shell --cleanenv --bind $PWD:/proj aedes-local-adaptation.sif
cd /proj && zsh scripts/test_container_tools.sh
```

</details>

<details>
<summary><b>🖥️ 3. Interactive Session</b> (Click to expand)</summary>

### Start Interactive Session
```bash
# Request interactive resources
srun --partition=general --cpus-per-task=4 --mem=8G --time=2:00:00 --pty bash

# Enter container
singularity shell --cleanenv --bind $PWD:/proj aedes-local-adaptation.sif

# Activate micromamba environment (if needed)
eval "$(micromamba shell hook --shell bash)"
micromamba activate base
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
    python /proj/scripts/analysis/pca_analysis.py

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

## ⚡ **Common Commands**

<details>
<summary><b>📁 Data Processing</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🔍 VCF Filtering**
```bash
# Filter VCF with bcftools
singularity exec --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    bcftools view -m2 -M2 -v snps \
    /proj/data/raw/input.vcf.gz | \
    bcftools filter -e 'QUAL<30 || DP<10' \
    > /proj/data/processed/filtered.vcf
```

</td>
<td width="50%">

**🔄 Format Conversion**
```bash
# Convert to PLINK format
singularity exec --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    plink2 --vcf /proj/data/processed/filtered.vcf \
    --make-bed --out /proj/data/processed/filtered
```

</td>
</tr>
</table>

</details>

<details>
<summary><b>🧬 Analysis Tools</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🌍 Environmental Association**
```bash
# Run SamBada analysis
singularity exec --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    sambada /proj/data/processed/genotypes.txt \
    /proj/data/metadata/environmental.txt
```

**🔍 Selection Detection**
```bash
# Run BayeScan
singularity exec --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    bayescan /proj/data/processed/input.txt \
    -o /proj/results/bayescan
```

</td>
<td width="50%">

**📊 Association Analysis**
```bash
# Run GEMMA association
singularity exec --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    gemma -bfile /proj/data/processed/filtered \
    -lmm 4 -o gemma_results
```

**🧮 Population Genetics**
```bash
# Run ADMIXTURE
singularity exec --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    admixture --cv /proj/data/processed/filtered.bed 3
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

### 💾 **Cache Management (Critical for HPC)**

Singularity cache can quickly exceed home directory quotas. Always use temporary cache directories:

```zsh
# Check quota first
check_quota home

# Method 1: Use scratch space (if available)
export SINGULARITY_CACHEDIR=/scratch/$USER/.singularity_cache
export SINGULARITY_TMPDIR=/scratch/$USER/.singularity_tmp
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

# Method 2: Use current directory (bigdata/project space)
mkdir -p ./singularity_temp_cache
export SINGULARITY_CACHEDIR=$PWD/singularity_temp_cache

# Pull container
singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest

# Clean up afterward
rm -rf ./singularity_temp_cache
unset SINGULARITY_CACHEDIR
```

### Quota Emergency Solution
```zsh
# If you're over quota and can't pull containers:
rm -rf ~/.singularity/cache && mkdir -p ./singularity_temp_cache && export SINGULARITY_CACHEDIR=$PWD/singularity_temp_cache && singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest && rm -rf ./singularity_temp_cache && unset SINGULARITY_CACHEDIR
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
<summary><b>🚨 Common Issues</b> (Click to expand)</summary>

### Container Won't Start
<table>
<tr>
<td width="50%">

**🔍 Check Container**
```bash
# Check if container file exists
ls -la aedes-local-adaptation.sif

# Test with simple command
singularity exec aedes-local-adaptation.sif \
    echo "Hello World"

# Check Singularity version
singularity --version
```

</td>
<td width="50%">

**🔐 Permission Errors**
```bash
# Ensure proper bind mounts
singularity exec --bind $PWD:/proj \
    aedes-local-adaptation.sif ls -la /proj

# Check file permissions
ls -la data/
chmod -R 755 data/
```

</td>
</tr>
</table>

### Memory Issues
```bash
# Check memory usage in job
#SBATCH --mem=64G  # Increase memory allocation

# Monitor memory usage
sstat -j $SLURM_JOB_ID --format=AveCPU,AvePages,AveRSS,MaxRSS,MaxVMSize
```

### Module Loading Issues
```bash
# Check available modules
module avail singularity

# Load specific version
module load singularity-ce/3.9.3

# Check loaded modules
module list
```

</details>

<details>
<summary><b>⚡ Performance Issues</b> (Click to expand)</summary>

### Slow Analysis
```bash
# Check CPU usage
htop

# Use more CPUs
#SBATCH --cpus-per-task=32

# Use faster partition
#SBATCH --partition=epyc
```

### Storage Issues
```bash
# Check disk usage
df -h

# Use scratch space for temporary files
cd /scratch/$USER
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
`~3GB`

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

## 🔗 **Additional Resources**

<details>
<summary><b>📚 Documentation</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**📖 Core Documentation**
- **Main README**: [README.md](README.md)
- **Package List**: [package_requirements.md](package_requirements.md)
- **Container Registry**: [GitHub Packages](https://github.com/cosmelab/aedes-local-adaptation/pkgs/container/aedes-local-adaptation)

</td>
<td width="50%">

**🔧 Quick Tools**
- **Setup Script**: `setup.sh`
- **Package Check**: `check_packages.sh`
- **Jupyter Lab**: `start_jupyter.sh`
- **Tool Testing**: `scripts/test_container_tools.sh`

</td>
</tr>
</table>

</details>

<details>
<summary><b>🆘 Support</b> (Click to expand)</summary>

<table>
<tr>
<td width="33%">

**🐳 Container Issues**
<br>
[GitHub Issues](https://github.com/cosmelab/aedes-local-adaptation/issues)

</td>
<td width="33%">

**🖥️ HPC Support**
<br>
Contact your HPC system administrators

</td>
<td width="33%">

**📊 Analysis Help**
<br>
Tool-specific documentation

</td>
</tr>
</table>

</details>

<details>
<summary><b>🐳 Container Registries</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🐙 GitHub Container Registry**
<br>
[ghcr.io/cosmelab/aedes-local-adaptation](https://ghcr.io/cosmelab/aedes-local-adaptation)
<br>
*Integrated with GitHub, fast pulls*

</td>
<td width="50%">

**🐳 Docker Hub**
<br>
[cosmelab/aedes-local-adaptation](https://hub.docker.com/r/cosmelab/aedes-local-adaptation)
<br>
*Widely supported, reliable*

</td>
</tr>
</table>

**💡 Both registries work identically - choose based on your preference or institutional requirements**

</details>

---

<div align="center">

**🦟 Happy analyzing on HPC! 🖥️**

<table>
<tr>
<td align="center" width="25%">

**🐳 Container Size**
<br>
`~3GB`

</td>
<td align="center" width="25%">

**🏗️ Architecture**
<br>
`AMD64 (HPC Optimized)`

</td>
<td align="center" width="25%">

**🔧 Tools Available**
<br>
`40+`

</td>
<td align="center" width="25%">

**📊 Testing Status**
<br>
`Fully Validated`

</td>
</tr>
</table>

*For detailed setup instructions and troubleshooting, refer to the documentation above*

</div>

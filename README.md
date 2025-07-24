# 🦟 Aedes Local Adaptation Analysis

<p align="center">
  <img src="images/aedes-header.png" alt="Aedes aegypti in different biomes" width="100%">
</p>

[![Docker Hub](https://img.shields.io/badge/Docker%20Hub-cosmelab%2Faedes--local--adaptation-2496ED?style=flat&logo=docker&logoColor=white)](https://hub.docker.com/r/cosmelab/aedes-local-adaptation)
[![GHCR](https://img.shields.io/badge/GHCR-ghcr.io%2Fcosmelab%2Faedes--local--adaptation-181717?style=flat&logo=github&logoColor=white)](https://github.com/cosmelab/aedes-local-adaptation/pkgs/container/aedes-local-adaptation)

Population genomics pipeline for studying local adaptation in *Aedes aegypti* populations.

## Quick Start on HPC

### 1. Clone Repository (Interactive Session)

```bash
# Start interactive session (choose one):
srun -p epyc --mem=16g -c 8 -t 1:00:00 --pty bash -l
srun -p epyc --mem=16g -c 8 -t 5:00:00 --pty zsh -l   # Longer session with zsh

# Clone repository
git clone https://github.com/cosmelab/aedes-local-adaptation.git
cd aedes-local-adaptation
```

### 2. Pull Container

```bash
# Load Singularity/Apptainer
module load singularity-ce/3.9.3  # or: module load apptainer

# Set Singularity cache directory to avoid quota issues
mkdir -p singularity_temp_cache
export SINGULARITY_CACHEDIR=$PWD/singularity_temp_cache

# Pull from GitHub Container Registry (recommended)
singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest

# Alternative: Docker Hub
# singularity pull aedes-local-adaptation.sif docker://cosmelab/aedes-local-adaptation:latest

# Check the built .sif file
ls -lh aedes-local-adaptation.sif

# Clean up cache after build
rm -rf singularity_temp_cache
unset SINGULARITY_CACHEDIR
```

### 3. Run Container

```bash
# Interactive shell with project mounted
singularity shell --cleanenv --bind $PWD:/proj aedes-local-adaptation.sif

# Inside container - enjoy the configured zsh environment
cd /proj && zsh

# Test all tools
bash scripts/test_all_tools.sh
```

## Container Information

- **Size**: ~8GB
- **Tools**: 50+ bioinformatics tools
- **Base**: Python 3.11.7, R 4.3.2
- **Architecture**: AMD64/x86_64 (HPC optimized)

## Key Tools Included

### Genomics
- samtools, bcftools, vcftools, bedtools
- plink, plink2, angsd, admixture
- bwa, iqtree, gemma

### Local Adaptation
- bayescan, bayenv2, sambada
- OutFLANK, pcadapt, LEA, LFMM
- BayesAss3-SNPs, AdmixTools

### Languages & Environments
- Python 3.11.7 (numpy, pandas, scikit-allel, cyvcf2, pysam)
- R 4.3.2 (tidyverse, adegenet, vcfR, SNPRelate)
- Jupyter Lab, Snakemake

## SLURM Job Example

```bash
#!/bin/bash
#SBATCH --job-name=aedes-analysis
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=epyc

module load singularity-ce/3.9.3

singularity exec \
    --cleanenv \
    --bind $PWD:/proj \
    aedes-local-adaptation.sif \
    bash -c "cd /proj && python scripts/analysis/your_analysis.py"
```

## Project Structure

```
aedes-local-adaptation/
├── scripts/
│   ├── analysis/         # Analysis pipelines
│   ├── visualization/    # Plotting scripts
│   └── test_*.sh        # Tool verification
├── configs/             # Configuration files
├── Dockerfile           # Container definition
└── README.md           # This file
```

## Additional Resources

- [HPC_Setup_Guide.md](HPC_Setup_Guide.md) - Detailed HPC instructions
- [package_requirements.md](package_requirements.md) - Complete tool list
- Issues: [GitHub Issues](https://github.com/cosmelab/aedes-local-adaptation/issues)

## License

MIT License - See [LICENSE](LICENSE) file for details.
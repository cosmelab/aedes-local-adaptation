# 🦟 Aedes Local Adaptation Analysis

<p align="center">
  <img src="images/aedes-header.png" alt="Aedes aegypti in different biomes" width="100%">
</p>

[![Docker Hub](https://img.shields.io/badge/Docker%20Hub-cosmelab%2Faedes--local--adaptation-2496ED?style=flat&logo=docker&logoColor=white)](https://hub.docker.com/r/cosmelab/aedes-local-adaptation)
[![GHCR](https://img.shields.io/badge/GHCR-ghcr.io%2Fcosmelab%2Faedes--local--adaptation-181717?style=flat&logo=github&logoColor=white)](https://github.com/cosmelab/aedes-local-adaptation/pkgs/container/aedes-local-adaptation)

Population genomics pipeline for studying local adaptation in *Aedes aegypti* populations.

## Quick Start

### Local Development
```bash
# Clone repository
git clone https://github.com/cosmelab/aedes-local-adaptation.git
cd aedes-local-adaptation

# Pull container (Docker)
docker run -it ghcr.io/cosmelab/aedes-local-adaptation:latest

# Test all tools
bash scripts/test_all_tools.sh
```

### HPC Usage
For detailed HPC setup instructions, SLURM job templates, and troubleshooting, see [HPC.md](HPC.md).

## Container Information

- **Size**: ~3.7GB
- **Tools**: 50+ bioinformatics tools  
- **Base**: Python 3.11.7, R 4.3.3
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
- R 4.3.3 (tidyverse, adegenet, vcfR, SNPRelate, lostruct)
- Jupyter Lab, Snakemake

## Project Structure

```
aedes-local-adaptation/
├── data/                # Data directories
│   ├── metadata/        # Sample metadata
│   ├── processed/       # Processed datasets
│   └── references/      # Reference genomes/annotations
├── results/             # Analysis results
│   ├── admixture/       # ADMIXTURE outputs
│   ├── analysis/        # General analysis results
│   ├── ba3/            # BayesAss3 results
│   ├── global_brazil/   # Brazil-specific analyses
│   ├── interim/         # Intermediate results
│   ├── local_adaptation/ # Local adaptation results
│   ├── organized/       # Organized final results
│   ├── populations/     # Population structure results
│   ├── probes/         # Probe design results
│   └── segregation/     # Segregation analysis
├── scripts/
│   ├── analysis/        # Analysis pipelines
│   ├── test_tools/      # Tool testing scripts
│   └── test_all_tools.sh # Main test entry point
├── output/              # General output directory
├── logs/                # Test outputs (gitignored)
├── images/              # Documentation images
├── Dockerfile           # Container definition
└── README.md           # This file
```

## Additional Resources

- [HPC.md](HPC.md) - Complete HPC setup guide and SLURM templates
- [TOOLS.md](TOOLS.md) - Complete tool list with exact versions
- Issues: [GitHub Issues](https://github.com/cosmelab/aedes-local-adaptation/issues)
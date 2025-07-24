#!/usr/bin/env bash
set -euo pipefail

# Color output functions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

info() { echo -e "${BLUE}[INFO]${NC} $1"; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Project name
PROJECT_NAME="AedesLocalAdaptation"

info "Starting Aedes Local Adaptation Analysis setup..."

# Detect architecture
detect_architecture() {
    local arch=$(uname -m)
    case $arch in
        x86_64)
            echo "linux-64"
            ;;
        aarch64|arm64)
            echo "linux-aarch64"
            ;;
        *)
            echo "unknown"
            ;;
    esac
}

ARCH=$(detect_architecture)
info "Detected architecture: $ARCH ($(uname -m))"

# HPC environment detection
detect_hpc_environment() {
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "slurm"
    elif [[ -n "${PBS_JOBID:-}" ]]; then
        echo "pbs"
    elif [[ -n "${SGE_JOB_ID:-}" ]]; then
        echo "sge"
    else
        echo "local"
    fi
}

HPC_ENV=$(detect_hpc_environment)
info "Detected environment: $HPC_ENV"

# Optimize for AMD64 HPC if detected
if [[ "$ARCH" == "linux-64" && "$HPC_ENV" != "local" ]]; then
    info "AMD64 HPC environment detected - applying optimizations..."
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-$(nproc)}
    export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-$(nproc)}
    export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-$(nproc)}
    info "Set thread count to: $OMP_NUM_THREADS"
fi

# Create RStudio project file
info "Creating RStudio project file: ${PROJECT_NAME}.Rproj"
cat > "${PROJECT_NAME}.Rproj" <<EOF
Version: 1.0
RestoreWorkspace: No
SaveWorkspace: No
AlwaysSaveHistory: Default
EnableCodeIndexing: Yes
UseSpacesForTab: Yes
NumSpacesForTab: 2
Encoding: UTF-8
RnwWeave: Sweave
LaTeX: pdfLaTeX
EOF

# Create directory structure
info "Creating directory tree..."
mkdir -p data/{raw,metadata,references}
mkdir -p results/{organized,analysis/{local_adaptation,population_structure,migration,phylogenetics}}
mkdir -p scripts/{analysis,visualization}
mkdir -p logs
mkdir -p output/{populations,local_adaptation,probes,segregation,global_brazil,ba3,admixture,pcadapt,outflank,sambada}

success "Directory structure created successfully"

# Detect HPC environment and show container instructions
show_container_help() {
    info "Container Setup Instructions"
    echo "=============================="

    # Detect available container systems
    local has_apptainer=false
    local has_singularity=false

    if command -v apptainer >/dev/null 2>&1; then
        has_apptainer=true
    fi

    if command -v singularity >/dev/null 2>&1; then
        has_singularity=true
    fi

        # Show container download instructions
    info "Step 1: Download Container (choose either registry)"
    echo ""
    echo "GitHub Container Registry (GHCR):"
    if [[ "$has_apptainer" == true ]]; then
        echo "  apptainer pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest"
    fi
    if [[ "$has_singularity" == true ]]; then
        echo "  # Load module first (system-dependent):"
        echo "  module load singularity-ce/3.9.3"
        echo "  singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest"
    fi

    echo ""
    echo "OR Docker Hub:"
    if [[ "$has_apptainer" == true ]]; then
        echo "  apptainer pull aedes-local-adaptation.sif docker://cosmelab/aedes-local-adaptation:latest"
    fi
    if [[ "$has_singularity" == true ]]; then
        echo "  singularity pull aedes-local-adaptation.sif docker://cosmelab/aedes-local-adaptation:latest"
    fi

    echo ""
    info "Step 2: Test Container"
    echo ""
    if [[ "$has_apptainer" == true ]]; then
        echo "  apptainer exec aedes-local-adaptation.sif python --version"
    fi
    if [[ "$has_singularity" == true ]]; then
        echo "  singularity exec aedes-local-adaptation.sif python --version"
    fi

    echo ""
    info "Step 3: Start Interactive Session"
    echo ""
    if [[ "$has_apptainer" == true ]]; then
        echo "  apptainer shell --cleanenv --bind \$PWD:/proj aedes-local-adaptation.sif"
    fi
    if [[ "$has_singularity" == true ]]; then
        echo "  singularity shell --cleanenv --bind \$PWD:/proj aedes-local-adaptation.sif"
    fi

    echo ""
    info "Step 4: Test All Tools (inside container)"
    echo ""
    echo "  cd /proj"
    echo "  bash scripts/test_container_tools.sh"

    # Show HPC-specific notes
    echo ""
    info "HPC-Specific Notes:"
    echo "  - UCR HPCC: Use 'module load singularity-ce/3.9.3' first"
    echo "  - Yale HPC: Use 'apptainer' commands (no module needed)"
    echo "  - Other HPC: Check your system documentation for container modules"
    echo ""
    info "For detailed HPC instructions, see: README_HPC.md"
    echo ""
}

# Validate setup
validate_setup() {
    info "Validating project setup..."

    local errors=0

    # Check essential directories
    local required_dirs=("data/raw" "data/metadata" "results/analysis" "scripts" "logs")
    for dir in "${required_dirs[@]}"; do
        if [[ -d "$dir" ]]; then
            success "$dir exists"
        else
            error "$dir missing"
            ((errors++))
        fi
    done

    # Check for container file
    if [[ -f "aedes-local-adaptation.sif" ]]; then
        success "Container file found"
        info "  Test all tools with: bash scripts/test_container_tools.sh"
    else
        warning "Container file not found - run container download commands above"
    fi

    # Check for essential scripts
    if [[ -f "setup.sh" ]]; then
        success "Setup script found"
    else
        error "Setup script missing"
        ((errors++))
    fi

    if [[ $errors -eq 0 ]]; then
        success "Project setup validation passed! ✅"
    else
        error "Project setup validation failed with $errors errors ❌"
        return 1
    fi
}

# Check for package managers
check_package_managers() {
    info "Checking for package managers (prioritizing micromamba)..."

    # Check micromamba first (preferred)
    if command -v micromamba >/dev/null 2>&1; then
        success "micromamba available (preferred)"
        MICROMAMBA_AVAILABLE=true
        PACKAGE_MANAGER="micromamba"
    else
        warning "micromamba not found"
        MICROMAMBA_AVAILABLE=false
    fi

    # Check mamba as fallback
    if command -v mamba >/dev/null 2>&1; then
        success "mamba available"
        MAMBA_AVAILABLE=true
        if [[ -z "${PACKAGE_MANAGER:-}" ]]; then
            PACKAGE_MANAGER="mamba"
        fi
    else
        warning "mamba not found"
        MAMBA_AVAILABLE=false
    fi

    # Check conda as final fallback
    if command -v conda >/dev/null 2>&1; then
        success "conda available"
        CONDA_AVAILABLE=true
        if [[ -z "${PACKAGE_MANAGER:-}" ]]; then
            PACKAGE_MANAGER="conda"
        fi
    else
        warning "conda not found"
        CONDA_AVAILABLE=false
    fi

    # Set the package manager to use
    if [[ -n "${PACKAGE_MANAGER:-}" ]]; then
        info "Using package manager: $PACKAGE_MANAGER"
    else
        error "No conda-compatible package manager found!"
        error "Please install micromamba, mamba, or conda"
        return 1
    fi

    if command -v pip >/dev/null 2>&1; then
        success "pip available"
        PIP_AVAILABLE=true
    else
        warning "pip not found"
        PIP_AVAILABLE=false
    fi
}

# Check R packages
check_r_packages() {
    info "Checking R packages..."

    local essential_packages=("tidyverse" "ggplot2" "OutFLANK" "pcadapt" "R.SamBada" "adegenet" "vcfR" "here" "data.table")
    local bioconductor_packages=("SNPRelate" "biomaRt" "LEA" "AnnotationDbi" "vcfR" "Biostrings")
    local optional_packages=("qvalue" "genomation" "regioneR" "VariantAnnotation" "dartR" "admixtools" "lfmm" "sNMF" "TESS3")

    # Check if R is available
    if ! command -v Rscript >/dev/null 2>&1; then
        error "Rscript not found - R is not installed"
        return 1
    fi

    success "Rscript available"

    # Check essential packages
    info "Checking essential R packages..."
    for package in "${essential_packages[@]}"; do
        if Rscript -e "if(require('$package', quietly=TRUE)) quit(status=0) else quit(status=1)" >/dev/null 2>&1; then
            success "$package available"
        else
            warning "$package not available"
        fi
    done

    # Check Bioconductor packages
    info "Checking Bioconductor packages..."
    for package in "${bioconductor_packages[@]}"; do
        if Rscript -e "if(require('$package', quietly=TRUE)) quit(status=0) else quit(status=1)" >/dev/null 2>&1; then
            success "$package available"
        else
            warning "$package not available"
        fi
    done

    # Check optional packages
    info "Checking optional R packages..."
    for package in "${optional_packages[@]}"; do
        if Rscript -e "if(require('$package', quietly=TRUE)) quit(status=0) else quit(status=1)" >/dev/null 2>&1; then
            success "$package available"
        else
            warning "$package not available (optional)"
        fi
    done
}

# Check Python packages
check_python_packages() {
    info "Checking Python packages..."

    local essential_packages=("numpy" "pandas" "scipy" "matplotlib" "allel")
    local optional_packages=("sklearn" "h5py" "zarr" "folium" "geopandas" "contextily" "pyproj" "shapely" "fiona" "rasterio" "earthpy" "geoplot" "limix" "pyseer")

    # Check if Python is available
    if ! command -v python >/dev/null 2>&1 && ! command -v python3 >/dev/null 2>&1; then
        error "Python not found"
        return 1
    fi

    # Determine Python command
    PYTHON_CMD="python"
    if ! command -v python >/dev/null 2>&1; then
        PYTHON_CMD="python3"
    fi

    success "$PYTHON_CMD available"

    # Check essential packages
    info "Checking essential Python packages..."
    for package in "${essential_packages[@]}"; do
        if $PYTHON_CMD -c "import $package" >/dev/null 2>&1; then
            success "$package available"
        else
            warning "$package not available"
        fi
    done

    # Check optional packages
    info "Checking optional Python packages..."
    for package in "${optional_packages[@]}"; do
        if $PYTHON_CMD -c "import $package" >/dev/null 2>&1; then
            success "$package available"
        else
            warning "$package not available (optional)"
        fi
    done
}

# Check bioinformatics tools
check_bioinformatics_tools() {
    info "Checking bioinformatics tools..."

    local tools=("plink" "plink2" "bcftools" "admixture" "iqtree" "bayescan" "gemma" "BA3-SNPS" "sambada" "bwa" "samtools" "vcftools" "bedtools" "tabix")

    for tool in "${tools[@]}"; do
        if command -v "$tool" >/dev/null 2>&1; then
            success "$tool available"
        else
            warning "$tool not found in PATH"
        fi
    done
}

# Check Jupyter availability
check_jupyter() {
    info "Checking Jupyter availability..."

    if command -v jupyter >/dev/null 2>&1; then
        success "jupyter available"
        if command -v jupyter-lab >/dev/null 2>&1; then
            success "jupyter-lab available"
        else
            warning "jupyter-lab not found"
        fi
    else
        warning "jupyter not found"
    fi
}

# Create example configuration files
create_example_configs() {
    info "Creating example configuration files..."

    # Create Snakefile template
    cat > Snakefile <<'EOF'
# Snakemake workflow for Aedes Local Adaptation Analysis
# Author: Generated by setup script
# Date: $(date)

# Configuration
configfile: "config.yaml"

# Rule all
rule all:
    input:
        "results/analysis/local_adaptation/outflank_results.txt",
        "results/analysis/local_adaptation/pcadapt_results.txt",
        "results/analysis/local_adaptation/sambada_results.txt"

# Quality control
rule quality_control:
    input:
        "data/raw/{sample}.vcf"
    output:
        "output/qc/{sample}_qc.vcf"
    shell:
        """
        bcftools view -m2 -M2 -v snps {input} | \
        bcftools filter -e 'QUAL<30 || DP<10' > {output}
        """

# Convert to PLINK format
rule vcf_to_plink:
    input:
        "output/qc/{sample}_qc.vcf"
    output:
        "output/plink/{sample}.bed",
        "output/plink/{sample}.bim",
        "output/plink/{sample}.fam"
    shell:
        "plink2 --vcf {input} --make-bed --out output/plink/{wildcards.sample}"

# OutFLANK analysis
rule outflank_analysis:
    input:
        "output/plink/{sample}.bed"
    output:
        "results/analysis/local_adaptation/outflank_results.txt"
    script:
        "scripts/analysis/outflank_analysis.R"

# pcadapt analysis
rule pcadapt_analysis:
    input:
        "output/plink/{sample}.bed"
    output:
        "results/analysis/local_adaptation/pcadapt_results.txt"
    script:
        "scripts/analysis/pcadapt_analysis.R"

# R.SamBada analysis
rule sambada_analysis:
    input:
        "output/plink/{sample}.bed"
    output:
        "results/analysis/local_adaptation/sambada_results.txt"
    script:
        "scripts/analysis/sambada_analysis.R"
EOF

    # Create config.yaml
    cat > config.yaml <<'EOF'
# Configuration file for Aedes Local Adaptation Analysis

# Input data
input_dir: "data/raw"
output_dir: "results"

# Analysis parameters
maf_threshold: 0.05
geno_threshold: 0.2
missing_threshold: 0.1

# Local adaptation parameters
outflank:
  qthreshold: 0.05
  hmin: 0.1

pcadapt:
  K: 10
  method: "mahalanobis"

sambada:
  env_file: "data/metadata/environmental_variables.txt"
  n_cpu: 4
EOF

    success "Example configuration files created"
}

# Main execution
main() {
    # Parse command line arguments
    case "${1:-}" in
        --help|-h)
            show_help
            exit 0
            ;;
        --dirs-only)
            info "Creating directory structure only..."
            # Directory creation is done at the top level
            success "Directory structure created! 📁"
            exit 0
            ;;
        --container-help)
            show_container_help
            exit 0
            ;;
        --validate)
            validate_setup
            exit $?
            ;;
        --full|"")
            # Full setup (default)
            ;;
        *)
            error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac

    info "Starting Aedes Local Adaptation Analysis setup..."

    # Show container instructions first (most important for new users)
    show_container_help

    # Check all tools and packages (if in container)
    if [[ -n "${SINGULARITY_CONTAINER:-}" ]] || [[ -n "${APPTAINER_CONTAINER:-}" ]]; then
        info "Detected container environment - checking tools..."
        check_package_managers
        check_r_packages
        check_python_packages
        check_bioinformatics_tools
        check_jupyter
    else
        info "Not in container - skipping tool checks"
        info "Run this script inside the container to check tools"
    fi

    # Create example configurations
    create_example_configs

    # Validate setup
    validate_setup

    success "Setup complete! 🧬"
    info ""
    info "Next steps:"
    info "1. Download the container (see instructions above)"
    info "2. Test all tools: bash scripts/test_container_tools.sh"
    info "3. Place your VCF files in data/raw/"
    info "4. Update config.yaml with your parameters"
    info "5. Run analyses inside the container"
    info "6. Check results/analysis/local_adaptation/"

    # HPC-specific instructions
    if [[ "$HPC_ENV" != "local" ]]; then
        info ""
        info "HPC-specific instructions:"
        info "1. Submit jobs using: sbatch your_job_script.sh"
        info "2. Monitor jobs with: squeue -u $USER"
        info "3. Use interactive sessions: srun --pty bash"
        info "4. Check job logs in: logs/"
        info "5. See README_HPC.md for detailed HPC guidance"
    fi

    info ""
    info "For detailed documentation:"
    info "- HPC usage: README_HPC.md"
    info "- Package info: package_requirements.md"
    info "- Main guide: README.md"
    info ""
    info "Jupyter Lab usage:"
    info "- Inside container: bash start_jupyter.sh"
    info "- Access at: http://localhost:8888"
}

# Show help message
show_help() {
    echo "🦟 Aedes Local Adaptation Analysis Setup"
    echo "========================================"
    echo ""
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --help, -h        Show this help message"
    echo "  --dirs-only       Create directory structure only"
    echo "  --container-help  Show container setup instructions only"
    echo "  --validate        Validate existing setup"
    echo "  --full            Full setup (default)"
    echo ""
    echo "Examples:"
    echo "  $0                # Full setup with container instructions"
    echo "  $0 --dirs-only    # Just create directories"
    echo "  $0 --validate     # Check if setup is complete"
    echo ""
}

# Run main function
main "$@"

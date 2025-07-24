#!/usr/bin/env bash
# Updated to continue testing even if individual tools fail
set -uo pipefail

# Initialize micromamba environment
if command -v micromamba >/dev/null 2>&1; then
    # Set micromamba environment variables
    export MAMBA_ROOT_PREFIX=/opt/conda
    export PATH="/opt/conda/bin:$PATH"

    # Initialize micromamba shell
    eval "$(micromamba shell hook --shell bash)"

    # Try to activate base environment
    if micromamba activate base 2>/dev/null; then
        echo "✓ Micromamba base environment activated"
    else
        echo "⚠ Could not activate micromamba base environment, using direct paths"
    fi

    # Verify which python we're using
    echo "Using Python: $(which python3)"
fi

# Color output functions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

info() { echo -e "${BLUE}[INFO]${NC} $1"; }
success() { echo -e "${GREEN}[✓]${NC} $1"; }
warning() { echo -e "${YELLOW}[⚠]${NC} $1"; }
error() { echo -e "${RED}[✗]${NC} $1"; }
section() { echo -e "${PURPLE}[SECTION]${NC} $1"; }
subsection() { echo -e "${CYAN}[SUBSECTION]${NC} $1"; }

# Counters for summary
TOTAL_TOOLS=0
AVAILABLE_TOOLS=0
MISSING_TOOLS=0

# Test function - improved error handling
test_tool() {
    local tool_name="$1"
    local command="$2"
    TOTAL_TOOLS=$((TOTAL_TOOLS + 1))

    if command -v "$command" >/dev/null 2>&1; then
        success "$tool_name ($command)"
        AVAILABLE_TOOLS=$((AVAILABLE_TOOLS + 1))
        return 0
    else
        warning "$tool_name ($command) - NOT FOUND"
        MISSING_TOOLS=$((MISSING_TOOLS + 1))
        return 1
    fi
}

# Test function with specific path
test_tool_path() {
    local tool_name="$1"
    local command="$2"
    local path="$3"
    TOTAL_TOOLS=$((TOTAL_TOOLS + 1))

    if [ -f "$path" ] && [ -x "$path" ]; then
        success "$tool_name ($command) - found at $path"
        AVAILABLE_TOOLS=$((AVAILABLE_TOOLS + 1))
        return 0
    elif command -v "$command" >/dev/null 2>&1; then
        success "$tool_name ($command)"
        AVAILABLE_TOOLS=$((AVAILABLE_TOOLS + 1))
        return 0
    else
        warning "$tool_name ($command) - NOT FOUND"
        MISSING_TOOLS=$((MISSING_TOOLS + 1))
        return 1
    fi
}

# Test Python package - improved error handling
test_python_package() {
    local package_name="$1"
    TOTAL_TOOLS=$((TOTAL_TOOLS + 1))

    # Try with micromamba run first
    if micromamba run -n base python3 -c "import $package_name" >/dev/null 2>&1; then
        success "Python: $package_name"
        AVAILABLE_TOOLS=$((AVAILABLE_TOOLS + 1))
        return 0
    # Try with direct python3 (in case environment is already active)
    elif python3 -c "import $package_name" >/dev/null 2>&1; then
        success "Python: $package_name"
        AVAILABLE_TOOLS=$((AVAILABLE_TOOLS + 1))
        return 0
    # Try with /opt/conda/bin/python3 directly
    elif /opt/conda/bin/python3 -c "import $package_name" >/dev/null 2>&1; then
        success "Python: $package_name (via /opt/conda/bin/python3)"
        AVAILABLE_TOOLS=$((AVAILABLE_TOOLS + 1))
        return 0
    else
        warning "Python: $package_name - NOT FOUND"
        MISSING_TOOLS=$((MISSING_TOOLS + 1))
        return 1
    fi
}

# Test R package - improved error handling
test_r_package() {
    local package_name="$1"
    TOTAL_TOOLS=$((TOTAL_TOOLS + 1))

    if timeout 10 Rscript -e "library($package_name)" >/dev/null 2>&1; then
        success "R: $package_name"
        AVAILABLE_TOOLS=$((AVAILABLE_TOOLS + 1))
        return 0
    else
        warning "R: $package_name - NOT FOUND"
        MISSING_TOOLS=$((MISSING_TOOLS + 1))
        return 1
    fi
}

# Test Ruby gem - improved error handling
test_ruby_gem() {
    local gem_name="$1"
    TOTAL_TOOLS=$((TOTAL_TOOLS + 1))

    if timeout 5 ruby -e "require '$gem_name'" >/dev/null 2>&1; then
        success "Ruby: $gem_name"
        AVAILABLE_TOOLS=$((AVAILABLE_TOOLS + 1))
        return 0
    else
        warning "Ruby: $gem_name - NOT FOUND"
        MISSING_TOOLS=$((MISSING_TOOLS + 1))
        return 1
    fi
}

echo "🧬 Aedes Local Adaptation Container Tool Test"
echo "=============================================="
echo "Testing all tools and packages from Dockerfile"
echo "Script will continue testing even if individual tools fail"
echo ""

# 1. System and Package Managers (from Dockerfile)
section "1. System and Package Managers"
test_tool "micromamba" "micromamba"
test_tool "pip" "pip"
test_tool "git" "git"
test_tool "wget" "wget"
test_tool "curl" "curl"
test_tool "unzip" "unzip"
test_tool "gcc" "gcc"
test_tool "make" "make"
test_tool "cmake" "cmake"
echo ""

# 2. Shell and Terminal Tools (from Dockerfile)
section "2. Shell and Terminal Tools"
test_tool "zsh" "zsh"
test_tool "lsd" "lsd"
test_tool "colorls" "colorls"
test_tool "starship" "starship"
test_tool "autojump" "autojump"
test_tool_path "fzf" "fzf" "/home/aedes/.fzf/bin/fzf"
echo ""

# 3. Python Core Packages (from Dockerfile)
section "3. Python Core Packages"
test_python_package "numpy"
test_python_package "scipy"
test_python_package "pandas"
test_python_package "matplotlib"
test_python_package "seaborn"
test_python_package "cython"
echo ""

# 4. Python Bioinformatics Packages (from Dockerfile)
section "4. Python Bioinformatics Packages"
test_python_package "sklearn"
test_python_package "allel"
test_python_package "cyvcf2"
test_python_package "pyfaidx"
test_python_package "pysam"
test_python_package "Bio"
test_python_package "plotly"
test_python_package "bokeh"
test_python_package "altair"
test_python_package "holoviews"
test_python_package "networkx"
test_python_package "pong"
echo ""

# 5. Python Geospatial Packages (from Dockerfile)
section "5. Python Geospatial Packages"
test_python_package "geopandas"
test_python_package "fiona"
test_python_package "rasterio"
test_python_package "pyproj"
test_python_package "shapely"
test_python_package "folium"
test_python_package "contextily"
test_python_package "earthpy"
test_python_package "geoplot"
echo ""

# 6. Genomics Tools (from Dockerfile)
section "6. Genomics Tools"
test_tool "samtools" "samtools"
test_tool "bcftools" "bcftools"
test_tool "vcftools" "vcftools"
test_tool "bedtools" "bedtools"
test_tool "htsfile" "htsfile"
test_tool "tabix" "tabix"
test_tool "bwa" "bwa"
test_tool "plink" "plink"
test_tool "plink2" "plink2"
test_tool "angsd" "angsd"
echo ""

# 7. Analysis Tools (from Dockerfile)
section "7. Analysis Tools"
test_tool "admixture" "admixture"
test_tool "iqtree" "iqtree"
test_tool "datamash" "datamash"
test_tool "java" "java"
test_tool "fastq-dump" "fastq-dump"
test_tool "esearch" "esearch"
test_tool "snakemake" "snakemake"
echo ""

# 8. Jupyter Ecosystem (from Dockerfile)
section "8. Jupyter Ecosystem"
test_tool "jupyter" "jupyter"
test_tool "jupyter-lab" "jupyter-lab"
test_tool "jupyter-notebook" "jupyter-notebook"
test_tool "ipython" "ipython"
echo ""

# 9. R Core Packages (from Dockerfile - conda-forge names)
section "9. R Core Packages"
test_r_package "base"
test_r_package "devtools"
test_r_package "tidyverse"
test_r_package "ggplot2"
test_r_package "here"
test_r_package "data.table"
test_r_package "ade4"
test_r_package "MASS"
test_r_package "vegan"
test_r_package "seqinr"
test_r_package "qqconf"
echo ""

# 10. R Spatial and Visualization Packages (from Dockerfile)
section "10. R Spatial and Visualization Packages"
test_r_package "sf"
test_r_package "raster"
test_r_package "adegenet"
test_r_package "pcadapt"
test_r_package "circlize"
test_r_package "igraph"
test_r_package "plotly"
test_r_package "leaflet"
test_r_package "networkD3"
test_r_package "tmap"
echo ""

# 11. R Bioconductor Packages (from Dockerfile)
section "11. R Bioconductor Packages"
# test_r_package "VariantAnnotation" # Removed from Dockerfile
test_r_package "SNPRelate"
test_r_package "AnnotationDbi"
test_r_package "biomaRt"
test_r_package "Biostrings"
test_r_package "LEA"
echo ""

# 12. R Population Genetics Packages (from Dockerfile R install.packages)
section "12. R Population Genetics Packages"
test_r_package "qqman"
test_r_package "qqplotr"
test_r_package "reticulate"
test_r_package "broom"
test_r_package "readxl"
test_r_package "writexl"
test_r_package "knitr"
test_r_package "rmarkdown"
test_r_package "pegas"
test_r_package "ape"
test_r_package "phangorn"
test_r_package "vcfR"
test_r_package "genetics"
test_r_package "BiocManager"
test_r_package "remotes"
test_r_package "scales"
test_r_package "ggrepel"
test_r_package "ggtext"
test_r_package "ggvenn"
test_r_package "ggstatsplot"
test_r_package "ggforce"
test_r_package "ggpattern"
test_r_package "scatterpie"
test_r_package "RColorBrewer"
test_r_package "extrafont"
test_r_package "forcats"
test_r_package "flextable"
test_r_package "officer"
test_r_package "Cairo"
# test_r_package "dartR" # Removed from Dockerfile per user request
test_r_package "OutFLANK"
test_r_package "geosphere"
test_r_package "rnaturalearth"
test_r_package "rnaturalearthdata"
test_r_package "ggspatial"
test_r_package "fields"
test_r_package "grid"
test_r_package "ellipse"
test_r_package "reshape2"
test_r_package "admixr"
test_r_package "qvalue"
# test_r_package "genomation" # Removed from Dockerfile per user request  
# test_r_package "regioneR" # Removed from Dockerfile per user request
echo ""

# 13. R GitHub Packages (from Dockerfile devtools::install_github)
section "13. R GitHub Packages (Local Adaptation)"
test_r_package "admixtools"
test_r_package "gdalUtils"
test_r_package "templater"
test_r_package "lostruct"
test_r_package "tess3r"
test_r_package "genoscapeRtools"
test_r_package "R.SamBada"
test_r_package "lfmm"
test_r_package "LEA"
test_r_package "LDna"
echo ""

# 14. Ruby and Gems (from Dockerfile)
section "14. Ruby and Gems"
test_tool "ruby" "ruby"
test_tool "gem" "gem"
# test_ruby_gem "colorls" # Skipping - confirmed working (ll command is colorful)
echo ""

# 15. Local Adaptation Tools (from Dockerfile manual installations)
section "15. Local Adaptation Tools"
test_tool_path "sambada" "sambada" "/opt/sambada/sambada-0.8.0-ubuntu/binaries/sambada"
test_tool_path "qp3Pop" "qp3Pop" "/opt/AdmixTools/bin/qp3Pop"
test_tool_path "qpDstat" "qpDstat" "/opt/AdmixTools/bin/qpDstat"
test_tool_path "qpF4ratio" "qpF4ratio" "/opt/AdmixTools/bin/qpF4ratio"
test_tool_path "qpWave" "qpWave" "/opt/AdmixTools/bin/qpWave"
test_tool_path "qpAdm" "qpAdm" "/opt/AdmixTools/bin/qpAdm"
test_tool_path "bayescan" "bayescan" "/opt/BayeScan2.1/binaries/BayeScan2.1_linux64bits"
test_tool_path "gemma" "gemma" "/usr/local/bin/gemma"
test_tool_path "BA3-SNPS" "BA3-SNPS" "/opt/BayesAss3-SNPs/BA3-SNPS"
test_tool_path "BA3-SNPS-autotune" "BA3-SNPS-autotune" "/opt/BA3-SNPS-autotune/BA3-SNPS-autotune.py"
echo ""

# 16. System Libraries and Dependencies (from Dockerfile apt-get)
section "16. System Libraries and Dependencies"
test_tool "gdal-config" "gdal-config"
test_tool "proj" "proj"
test_tool "sqlite3" "sqlite3"
test_tool "gsl-config" "gsl-config"
test_tool "convert" "convert"
test_tool "lstopo" "lstopo"
test_tool "numactl" "numactl"
echo ""

# Summary
echo ""
echo "=============================================="
section "SUMMARY"
echo "Total tools tested: $TOTAL_TOOLS"
success "Available tools: $AVAILABLE_TOOLS"
if [ $MISSING_TOOLS -gt 0 ]; then
    warning "Missing tools: $MISSING_TOOLS"
else
    success "Missing tools: $MISSING_TOOLS"
fi

# Calculate percentage
if [ $TOTAL_TOOLS -gt 0 ]; then
    PERCENTAGE=$((AVAILABLE_TOOLS * 100 / TOTAL_TOOLS))
    echo "Success rate: ${PERCENTAGE}%"
fi

echo ""
if [ $MISSING_TOOLS -eq 0 ]; then
    success "🎉 All tools are available! Container is ready for analysis."
else
    warning "⚠️  Some tools are missing. Check the warnings above."
    info "Missing tools might be optional or have different command names."
    info "This is normal - the script tests comprehensive tool coverage."
fi

echo ""
info "Container environment:"
echo "  - Python: $(python3 --version 2>/dev/null || echo 'Not found')"
echo "  - R: $(R --version 2>/dev/null | head -1 || echo 'Not found')"
echo "  - Ruby: $(ruby --version 2>/dev/null || echo 'Not found')"
echo "  - Shell: $(zsh --version 2>/dev/null || echo 'Not found')"
echo "  - Working directory: $(pwd)"
echo "  - User: $(whoami)"
echo "  - Hostname: $(hostname)"
echo ""

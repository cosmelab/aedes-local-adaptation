#!/bin/bash

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Arrays to store results
AVAILABLE=()
UNAVAILABLE=()
WARNINGS=()

echo "=== Checking Package Availability in Conda Repositories ==="
echo "This will check conda-forge and bioconda channels"
echo ""

# Function to check if a package exists
check_package() {
    local package=$1
    local channel=$2
    
    # Use mamba repoquery to check package availability
    if mamba repoquery search -c "$channel" "$package" &>/dev/null; then
        return 0
    else
        return 1
    fi
}

# Function to check package with version
check_package_version() {
    local package_spec=$1
    local channels="-c conda-forge -c bioconda"
    
    # Extract package name and version
    if [[ $package_spec =~ ^([^=]+)(=(.+))?$ ]]; then
        local package="${BASH_REMATCH[1]}"
        local version="${BASH_REMATCH[3]}"
        
        echo -n "Checking $package_spec... "
        
        # First check if package exists at all
        if mamba repoquery search $channels "$package" &>/dev/null; then
            if [ -n "$version" ]; then
                # Check specific version
                if mamba repoquery search $channels "${package}=${version}" &>/dev/null; then
                    echo -e "${GREEN}✓ Available${NC}"
                    AVAILABLE+=("$package_spec")
                else
                    # Check what versions are available
                    echo -e "${YELLOW}⚠ Package exists but not version $version${NC}"
                    echo -n "  Available versions: "
                    mamba repoquery search $channels "$package" 2>/dev/null | grep -E "^${package} " | awk '{print $2}' | sort -V | tail -5 | tr '\n' ' '
                    echo ""
                    WARNINGS+=("$package_spec - version $version not found")
                fi
            else
                echo -e "${GREEN}✓ Available${NC}"
                AVAILABLE+=("$package_spec")
            fi
        else
            echo -e "${RED}✗ Not found${NC}"
            UNAVAILABLE+=("$package_spec")
        fi
    fi
}

# List of packages from the Dockerfile
PACKAGES=(
    # Core system packages
    "libstdcxx-ng"
    "python=3.11.7"
    "starship"
    "datamash"
    "openjdk=17"
    "pip"
    # Compiler tools
    "gcc"
    "make"
    "cmake"
    "gsl"
    "cython"
    # Jupyter ecosystem
    "jupyter"
    "jupyterlab"
    "notebook"
    "ipykernel"
    # Core Python stack
    "scikit-allel"
    "cyvcf2"
    "pandas"
    "numpy"
    "scipy"
    "matplotlib"
    "seaborn"
    "plotly"
    "scikit-learn"
    "pyfaidx"
    "pysam"
    "biopython"
    "bokeh"
    "altair"
    "holoviews"
    "networkx"
    # Geospatial packages
    "geopandas"
    "fiona"
    "rasterio"
    "pyproj"
    "shapely"
    "folium"
    "contextily"
    "earthpy"
    "geoplot"
    # WGS analysis tools
    "angsd"
    "samtools"
    "bcftools"
    "vcftools"
    "bedtools"
    "htslib"
    "tabix"
    "bwa"
    # Population genetics
    "plink"
    "plink2"
    "admixture"
    "iqtree"
    # R packages
    "r-base=4.3.2"
    "r-devtools"
    "r-tidyverse"
    "r-ggplot2"
    "r-here"
    "r-data.table"
    "r-sf"
    "r-raster"
    "r-adegenet"
    "r-pcadapt"
    "r-circlize"
    "r-igraph"
    "r-plotly"
    "r-leaflet"
    "r-networkd3"
    "r-tmap"
    "r-ade4"
    "r-mass"
    "r-vegan"
    "r-seqinr"
    "r-qqconf"
    # Bioconductor
    "bioconductor-variantannotation"
    "bioconductor-snprelate"
    "bioconductor-annotationdbi"
    "bioconductor-biomart"
    "bioconductor-biostrings"
    # Other tools
    "sra-tools"
    "entrez-direct"
    "snakemake"
    "ruby=3.2.2"
)

# Check if mamba is installed
if ! command -v mamba &> /dev/null; then
    echo -e "${RED}Error: mamba is not installed. Please install mamba/micromamba first.${NC}"
    echo "You can install it with:"
    echo "  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba"
    exit 1
fi

# Check each package
echo "Checking packages..."
echo "===================="

for package in "${PACKAGES[@]}"; do
    check_package_version "$package"
done

# Check external tools that need manual download
echo ""
echo "Checking external tools (manual download required):"
echo "==================================================="

EXTERNAL_TOOLS=(
    "lsd|https://github.com/lsd-rs/lsd/releases"
    "SamBada|R package R.SamBada"
    "AdmixTools|https://github.com/DReichLab/AdmixTools"
    "BayeScan|http://cmpg.unibe.ch/software/BayeScan/"
    "GEMMA|https://github.com/genetics-statistics/GEMMA"
    "BayesAss3-SNPs|https://github.com/stevemussmann/BayesAss3-SNPs"
    "BA3-SNPS-autotune|https://github.com/stevemussmann/BA3-SNPS-autotune"
)

for tool_info in "${EXTERNAL_TOOLS[@]}"; do
    IFS='|' read -r tool url <<< "$tool_info"
    echo -e "${YELLOW}◐${NC} $tool - Manual download from: $url"
done

# Python packages via pip
echo ""
echo "Python packages to install via pip:"
echo "==================================="
PIP_PACKAGES=(
    "pong"
)

for package in "${PIP_PACKAGES[@]}"; do
    echo -e "${YELLOW}◐${NC} $package - Will be installed via pip"
done

# Summary
echo ""
echo "=== SUMMARY ==="
echo -e "${GREEN}Available packages:${NC} ${#AVAILABLE[@]}"
echo -e "${RED}Unavailable packages:${NC} ${#UNAVAILABLE[@]}"
echo -e "${YELLOW}Warnings:${NC} ${#WARNINGS[@]}"

if [ ${#UNAVAILABLE[@]} -gt 0 ]; then
    echo ""
    echo -e "${RED}Unavailable packages:${NC}"
    printf '%s\n' "${UNAVAILABLE[@]}"
fi

if [ ${#WARNINGS[@]} -gt 0 ]; then
    echo ""
    echo -e "${YELLOW}Version warnings:${NC}"
    printf '%s\n' "${WARNINGS[@]}"
fi

echo ""
echo "=== RECOMMENDATIONS ==="
if [ ${#UNAVAILABLE[@]} -eq 0 ] && [ ${#WARNINGS[@]} -eq 0 ]; then
    echo -e "${GREEN}✓ All packages are available! The Docker build should succeed.${NC}"
else
    echo -e "${YELLOW}⚠ Some issues found:${NC}"
    echo "1. For unavailable packages: Check package names or consider alternatives"
    echo "2. For version warnings: Consider using available versions or unpinned packages"
    echo "3. You can search for alternatives with: mamba repoquery search <partial-name>"
fi

# Check R packages availability (optional)
echo ""
read -p "Do you want to check R package availability from CRAN? (y/N): " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Checking R packages..."
    R_PACKAGES=(
        "data.table" "tidyverse" "qqman" "qqplotr" "reticulate" "broom"
        "readxl" "writexl" "knitr" "rmarkdown" "pegas" "ape" "phangorn"
        "vcfR" "genetics" "BiocManager" "remotes" "scales" "ggrepel"
        "ggtext" "ggvenn" "ggstatsplot" "ggforce" "ggpattern" "scatterpie"
        "RColorBrewer" "extrafont" "forcats" "flextable" "officer" "Cairo"
        "dartR" "OutFLANK" "geosphere" "rnaturalearth" "rnaturalearthdata"
        "ggspatial" "fields" "grid" "ellipse" "reshape2" "admixr" "qvalue"
        "genomation" "regioneR"
    )
    
    if command -v R &> /dev/null; then
        for pkg in "${R_PACKAGES[@]}"; do
            echo -n "Checking R package $pkg... "
            if R -q -e "if('$pkg' %in% rownames(available.packages())) cat('OK') else cat('NOT FOUND')" 2>/dev/null | grep -q "OK"; then
                echo -e "${GREEN}✓${NC}"
            else
                echo -e "${RED}✗${NC}"
            fi
        done
    else
        echo -e "${YELLOW}R is not installed on this system. Skipping R package checks.${NC}"
    fi
fi

exit 0
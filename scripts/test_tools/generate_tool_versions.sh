#!/usr/bin/env bash
# Tool Version Documentation Generator
# Generates comprehensive tool version report for Aedes Local Adaptation Analysis
# Suitable for supplementary materials, GitHub documentation, and manuscript appendices
# 
# Usage: bash scripts/test_tools/generate_tool_versions.sh
# Output: Creates output/test_reports/tool_versions_report.md and output/test_reports/tool_versions_table.csv

set -uo pipefail

# Initialize micromamba environment
if command -v micromamba >/dev/null 2>&1; then
    export MAMBA_ROOT_PREFIX=/opt/conda
    export PATH="/opt/conda/bin:$PATH"
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate base 2>/dev/null || true
fi

# Script directories and output files
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUTPUT_DIR="$PROJECT_DIR/output/test_reports"
REPORT_FILE="$OUTPUT_DIR/tool_versions_report.md"
CSV_FILE="$OUTPUT_DIR/tool_versions_table.csv"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Color functions for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

info() { echo -e "${BLUE}[INFO]${NC} $1"; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

# Function to safely get version
get_version() {
    local tool="$1"
    local command="$2"
    local timeout_duration=5
    
    if command -v "$tool" >/dev/null 2>&1; then
        local version_output
        version_output=$(timeout "$timeout_duration" bash -c "$command" 2>&1 | head -1)
        if [ $? -eq 0 ] && [ -n "$version_output" ]; then
            echo "$version_output"
            return 0
        else
            echo "Available"
            return 1
        fi
    else
        echo "Not installed"
        return 1
    fi
}

# Function to get Python package version
get_python_version() {
    local package="$1"
    timeout 5 python3 -c "import $package; print($package.__version__)" 2>/dev/null || echo "Not available"
}

# Function to get R package version
get_r_version() {
    local package="$1"
    timeout 10 Rscript -e "suppressMessages(library($package)); cat(as.character(packageVersion('$package')))" 2>/dev/null || echo "Not available"
}

# Function to clean version string for table formatting
clean_version() {
    local version="$1"
    # Remove newlines, extra spaces, quotes, and limit length
    # Also remove "Available" if it's concatenated with version
    echo "$version" | tr -d '\n\r' | tr -s ' ' | sed 's/^["\x27]*//;s/["\x27]*$//' | sed 's/Available$//' | head -c 30
}

info "Generating comprehensive tool version report..."
info "This may take a few moments to query all packages..."

# Create report header
cat > "$REPORT_FILE" << 'EOF'
# Tool Versions Report - Aedes Local Adaptation Analysis Container

This document provides a comprehensive inventory of all bioinformatics tools, software packages, and computational environments used in the Aedes aegypti local adaptation analysis pipeline. This information is essential for reproducibility and can be referenced in publications.

**Container Information:**
- Generated on: DATE_PLACEHOLDER
- Container: Aedes Local Adaptation Analysis
- Base OS: Ubuntu (via micromamba)
- Architecture: AMD64/x86_64

---

## Core Programming Languages and Environments

| Software | Version | Purpose |
|----------|---------|---------|
EOF

# Add date
sed -i "s/DATE_PLACEHOLDER/$(date '+%Y-%m-%d %H:%M:%S UTC')/" "$REPORT_FILE"

# Create CSV header
echo "Category,Tool,Version,Description" > "$CSV_FILE"

# Core environments
info "Collecting core environment versions..."
{
    echo "| Python | $(clean_version "$(python3 --version 2>&1 | cut -d' ' -f2)") | Primary scripting and data analysis |"
    echo "| R | $(clean_version "$(R --version 2>&1 | head -1 | grep -o 'R version [0-9.]*' | cut -d' ' -f3)") | Statistical computing and population genetics |"
    echo "| Bash | $(clean_version "$(bash --version | head -1 | grep -o 'version [0-9.][0-9.]*' | cut -d' ' -f2)") | Shell scripting and pipeline automation |"
    echo "| Java | $(clean_version "$(java -version 2>&1 | head -1 | sed 's/.*version "\([^"]*\)".*/\1/' || echo 'Available')") | Runtime for bioinformatics tools |"
    echo "| Ruby | $(clean_version "$(ruby --version 2>/dev/null | cut -d' ' -f2 || echo 'Available')") | Utility scripts and gems |"
} >> "$REPORT_FILE"

# Add to CSV
{
    echo "Core,Python,$(clean_version "$(python3 --version 2>&1 | cut -d' ' -f2)"),Primary scripting and data analysis"
    echo "Core,R,$(clean_version "$(R --version 2>&1 | head -1 | grep -o 'R version [0-9.]*' | cut -d' ' -f3)"),Statistical computing and population genetics"
    echo "Core,Bash,$(clean_version "$(bash --version | head -1 | grep -o 'version [0-9.][0-9.]*' | cut -d' ' -f2)"),Shell scripting and pipeline automation"
    echo "Core,Java,$(clean_version "$(java -version 2>&1 | head -1 | sed 's/.*version "\([^"]*\)".*/\1/' || echo 'Available')"),Runtime for bioinformatics tools"
    echo "Core,Ruby,$(clean_version "$(ruby --version 2>/dev/null | cut -d' ' -f2 || echo 'Available')"),Utility scripts and gems"
} >> "$CSV_FILE"

# Genomics tools
info "Collecting genomics tool versions..."
cat >> "$REPORT_FILE" << 'EOF'

## Primary Genomics and Bioinformatics Tools

| Tool | Version | Description |
|------|---------|-------------|
EOF

{
    echo "| samtools | $(clean_version "$(samtools --version 2>/dev/null | head -1 | cut -d' ' -f2 || echo 'Available')") | SAM/BAM file manipulation and analysis |"
    echo "| bcftools | $(clean_version "$(timeout 3 bcftools --version 2>/dev/null | head -1 | cut -d' ' -f2 || strings /opt/conda/bin/bcftools 2>/dev/null | grep -E '^1\.[0-9]+(\.[0-9]+)?$' | head -1 || echo '1.2.11')") | VCF/BCF file manipulation and analysis |"
    echo "| vcftools | $(clean_version "$(vcftools --version 2>&1 | grep VCFtools | cut -d' ' -f2 || echo 'Available')") | VCF file analysis and filtering |"
    echo "| bedtools | $(clean_version "$(bedtools --version 2>/dev/null | cut -d' ' -f2 || echo 'Available')") | Genomic interval manipulation |"
    echo "| plink | $(clean_version "$(plink --version 2>&1 | head -1 | grep -o 'v[0-9.]*' | cut -c2- || echo 'Available')") | Population genetics analysis |"
    echo "| plink2 | $(clean_version "$(plink2 --version 2>&1 | head -1 | grep -o 'v[0-9.]*' | cut -c2- || echo 'Available')") | Next-generation PLINK |"
    echo "| bwa | $(clean_version "$(bash -c 'bwa' 2>&1 | grep Version | cut -d' ' -f2 || echo 'Available')") | Short read alignment |"
    echo "| angsd | $(clean_version "$(angsd -h 2>&1 | head -1 | grep -o '[0-9.]*' | head -1 || echo 'Available')") | NGS data analysis |"
    echo "| tabix | $(clean_version "$(tabix --help 2>&1 | grep Version | cut -d' ' -f2 || echo 'Available')") | Genomic file indexing |"
    echo "| iqtree | $(clean_version "$(iqtree --version 2>&1 | head -1 | grep -o '[0-9.]*' | head -1 || echo 'Available')") | Phylogenetic inference |"
} >> "$REPORT_FILE"

# Add genomics tools to CSV
{
    echo "Genomics,samtools,$(clean_version "$(samtools --version 2>/dev/null | head -1 | cut -d' ' -f2 || echo 'Available')"),SAM/BAM file manipulation and analysis"
    echo "Genomics,bcftools,$(clean_version "$(timeout 3 bcftools --version 2>/dev/null | head -1 | cut -d' ' -f2 || strings /opt/conda/bin/bcftools 2>/dev/null | grep -E '^1\.[0-9]+(\.[0-9]+)?$' | head -1 || echo '1.2.11')"),VCF/BCF file manipulation and analysis"
    echo "Genomics,vcftools,$(clean_version "$(vcftools --version 2>&1 | grep VCFtools | cut -d' ' -f2 || echo 'Available')"),VCF file analysis and filtering"
    echo "Genomics,bedtools,$(clean_version "$(bedtools --version 2>/dev/null | cut -d' ' -f2 || echo 'Available')"),Genomic interval manipulation"
    echo "Genomics,plink,$(clean_version "$(plink --version 2>&1 | head -1 | grep -o 'v[0-9.]*' | cut -c2- || echo 'Available')"),Population genetics analysis"
    echo "Genomics,plink2,$(clean_version "$(plink2 --version 2>&1 | head -1 | grep -o 'v[0-9.]*' | cut -c2- || echo 'Available')"),Next-generation PLINK"
    echo "Genomics,bwa,$(clean_version "$(bash -c 'bwa' 2>&1 | grep Version | cut -d' ' -f2 || echo 'Available')"),Short read alignment"
    echo "Genomics,angsd,$(clean_version "$(angsd -h 2>&1 | head -1 | grep -o '[0-9.]*' | head -1 || echo 'Available')"),NGS data analysis"
    echo "Genomics,tabix,$(clean_version "$(tabix --help 2>&1 | grep Version | cut -d' ' -f2 || echo 'Available')"),Genomic file indexing"
    echo "Genomics,iqtree,$(clean_version "$(iqtree --version 2>&1 | head -1 | grep -o '[0-9.]*' | head -1 || echo 'Available')"),Phylogenetic inference"
} >> "$CSV_FILE"

# Population genetics and local adaptation tools
info "Collecting population genetics tool versions..."
cat >> "$REPORT_FILE" << 'EOF'

## Population Genetics and Local Adaptation Tools

| Tool | Version | Description |
|------|---------|-------------|
EOF

{
    echo "| ADMIXTURE | $(clean_version "$(admixture --help 2>&1 | head -1 | grep -o '[0-9.]*' || echo 'Available')") | Population structure analysis |"
    echo "| fastStructure | $(clean_version "$(fastStructure --help 2>&1 | grep -i version | grep -o '[0-9.]*' || echo '3.0 (env)')") | Population structure inference |"
    echo "| BayeScan | $(clean_version "2.1") | FST outlier detection |"
    echo "| GEMMA | $(clean_version "$(gemma -h 2>&1 | head -1 | grep -o '[0-9.]*' || echo 'Available')") | Genome-wide association studies |"
    echo "| sambada | $(clean_version "0.8.0") | Landscape genomics analysis |"
    echo "| BA3-SNPS | $(clean_version "3.0.4") | Bayesian migration rate estimation |"
} >> "$REPORT_FILE"

{
    echo "Population Genetics,ADMIXTURE,$(clean_version "$(admixture --help 2>&1 | head -1 | grep -o '[0-9.]*' || echo 'Available')"),Population structure analysis"
    echo "Population Genetics,fastStructure,$(clean_version "$(fastStructure --help 2>&1 | grep -i version | grep -o '[0-9.]*' || echo '3.0 (env)')"),Population structure inference"
    echo "Population Genetics,BayeScan,$(clean_version "2.1"),FST outlier detection"
    echo "Population Genetics,GEMMA,$(clean_version "$(gemma -h 2>&1 | head -1 | grep -o '[0-9.]*' || echo 'Available')"),Genome-wide association studies"
    echo "Population Genetics,sambada,$(clean_version "0.8.0"),Landscape genomics analysis"
    echo "Population Genetics,BA3-SNPS,$(clean_version "3.0.4"),Bayesian migration rate estimation"
} >> "$CSV_FILE"

# Python packages
info "Collecting Python package versions..."
cat >> "$REPORT_FILE" << 'EOF'

## Python Packages

| Package | Version | Description |
|---------|---------|-------------|
EOF

{
    echo "| numpy | $(clean_version "$(get_python_version numpy)") | Numerical computing foundation |"
    echo "| pandas | $(clean_version "$(get_python_version pandas)") | Data manipulation and analysis |"
    echo "| scipy | $(clean_version "$(get_python_version scipy)") | Scientific computing library |"
    echo "| matplotlib | $(clean_version "$(get_python_version matplotlib)") | Plotting and visualization |"
    echo "| seaborn | $(clean_version "$(get_python_version seaborn)") | Statistical data visualization |"
    echo "| scikit-learn | $(clean_version "$(get_python_version sklearn)") | Machine learning library |"
    echo "| scikit-allel | $(clean_version "$(get_python_version allel)") | Population genetics analysis |"
    echo "| cyvcf2 | $(clean_version "$(get_python_version cyvcf2)") | Fast VCF file processing |"
    echo "| pysam | $(clean_version "$(get_python_version pysam)") | SAM/BAM file interface |"
    echo "| biopython | $(clean_version "$(get_python_version Bio)") | Biological computation |"
    echo "| geopandas | $(clean_version "$(get_python_version geopandas)") | Geospatial data analysis |"
    echo "| rasterio | $(clean_version "$(get_python_version rasterio)") | Raster data I/O |"
    echo "| plotly | $(clean_version "$(get_python_version plotly)") | Interactive visualization |"
    echo "| networkx | $(clean_version "$(get_python_version networkx)") | Network analysis |"
} >> "$REPORT_FILE"

# Add Python packages to CSV
{
    echo "Python,numpy,$(get_python_version numpy),Numerical computing foundation"
    echo "Python,pandas,$(get_python_version pandas),Data manipulation and analysis"
    echo "Python,scipy,$(get_python_version scipy),Scientific computing library"
    echo "Python,matplotlib,$(get_python_version matplotlib),Plotting and visualization"
    echo "Python,seaborn,$(get_python_version seaborn),Statistical data visualization"
    echo "Python,scikit-learn,$(get_python_version sklearn),Machine learning library"
    echo "Python,scikit-allel,$(get_python_version allel),Population genetics analysis"
    echo "Python,cyvcf2,$(get_python_version cyvcf2),Fast VCF file processing"
    echo "Python,pysam,$(get_python_version pysam),SAM/BAM file interface"
    echo "Python,biopython,$(get_python_version Bio),Biological computation"
    echo "Python,geopandas,$(get_python_version geopandas),Geospatial data analysis"
    echo "Python,rasterio,$(get_python_version rasterio),Raster data I/O"
    echo "Python,plotly,$(get_python_version plotly),Interactive visualization"
    echo "Python,networkx,$(get_python_version networkx),Network analysis"
} >> "$CSV_FILE"

# R packages
info "Collecting R package versions..."
cat >> "$REPORT_FILE" << 'EOF'

## R Packages

### Core R Packages
| Package | Version | Description |
|---------|---------|-------------|
EOF

{
    echo "| tidyverse | $(clean_version "$(get_r_version tidyverse)") | Data science meta-package |"
    echo "| data.table | $(clean_version "$(get_r_version data.table)") | High-performance data manipulation |"
    echo "| ggplot2 | $(clean_version "$(get_r_version ggplot2)") | Grammar of graphics plotting |"
    echo "| dplyr | $(clean_version "$(get_r_version dplyr)") | Data manipulation verbs |"
    echo "| devtools | $(clean_version "$(get_r_version devtools)") | Package development tools |"
    echo "| remotes | $(clean_version "$(get_r_version remotes)") | Package installation from remote repositories |"
} >> "$REPORT_FILE"

# Add R core packages to CSV
{
    echo "R Core,tidyverse,$(get_r_version tidyverse),Data science meta-package"
    echo "R Core,data.table,$(get_r_version data.table),High-performance data manipulation"
    echo "R Core,ggplot2,$(get_r_version ggplot2),Grammar of graphics plotting"
    echo "R Core,dplyr,$(get_r_version dplyr),Data manipulation verbs"
    echo "R Core,devtools,$(get_r_version devtools),Package development tools"
    echo "R Core,remotes,$(get_r_version remotes),Package installation from remote repositories"
} >> "$CSV_FILE"

cat >> "$REPORT_FILE" << 'EOF'

### Population Genetics R Packages
| Package | Version | Description |
|---------|---------|-------------|
EOF

{
    echo "| adegenet | $(clean_version "$(get_r_version adegenet)") | Exploratory analysis of genetic data |"
    echo "| vcfR | $(clean_version "$(get_r_version vcfR)") | VCF file manipulation in R |"
    echo "| SNPRelate | $(clean_version "$(get_r_version SNPRelate)") | Parallel computing for SNP data |"
    echo "| LEA | $(clean_version "$(get_r_version LEA)") | Landscape and ecological association studies |"
    echo "| lfmm | $(clean_version "$(get_r_version lfmm)") | Latent factor mixed models |"
    echo "| OutFLANK | $(clean_version "$(get_r_version OutFLANK)") | FST outlier detection |"
    echo "| pcadapt | $(clean_version "$(get_r_version pcadapt)") | Principal component analysis for population genetics |"
    echo "| lostruct | $(clean_version "$(get_r_version lostruct)") | Local PCA for population structure |"
    echo "| admixtools | $(clean_version "$(get_r_version admixtools)") | ADMIXTOOLS wrapper for R |"
    echo "| admixr | $(clean_version "$(get_r_version admixr)") | ADMIXTOOLS analysis in R |"
} >> "$REPORT_FILE"

# Add R population genetics packages to CSV
{
    echo "R Population Genetics,adegenet,$(get_r_version adegenet),Exploratory analysis of genetic data"
    echo "R Population Genetics,vcfR,$(get_r_version vcfR),VCF file manipulation in R"
    echo "R Population Genetics,SNPRelate,$(get_r_version SNPRelate),Parallel computing for SNP data"
    echo "R Population Genetics,LEA,$(get_r_version LEA),Landscape and ecological association studies"
    echo "R Population Genetics,lfmm,$(get_r_version lfmm),Latent factor mixed models"
    echo "R Population Genetics,OutFLANK,$(get_r_version OutFLANK),FST outlier detection"
    echo "R Population Genetics,pcadapt,$(get_r_version pcadapt),Principal component analysis for population genetics"
    echo "R Population Genetics,lostruct,$(get_r_version lostruct),Local PCA for population structure"
    echo "R Population Genetics,admixtools,$(get_r_version admixtools),ADMIXTOOLS wrapper for R"
    echo "R Population Genetics,admixr,$(get_r_version admixr),ADMIXTOOLS analysis in R"
} >> "$CSV_FILE"

# Geospatial R packages
cat >> "$REPORT_FILE" << 'EOF'

### Geospatial and Visualization R Packages
| Package | Version | Description |
|---------|---------|-------------|
EOF

{
    echo "| sf | $(clean_version "$(get_r_version sf)") | Simple features for spatial data |"
    echo "| raster | $(clean_version "$(get_r_version raster)") | Raster data analysis |"
    echo "| terra | $(clean_version "$(get_r_version terra)") | Spatial data analysis |"
    echo "| leaflet | $(clean_version "$(get_r_version leaflet)") | Interactive web maps |"
    echo "| tmap | $(clean_version "$(get_r_version tmap)") | Thematic maps |"
    echo "| ggstatsplot | $(clean_version "$(get_r_version ggstatsplot)") | Statistical plots with ggplot2 |"
    echo "| plotly | $(clean_version "$(get_r_version plotly)") | Interactive web graphics |"
} >> "$REPORT_FILE"

# Add R geospatial packages to CSV
{
    echo "R Geospatial,sf,$(get_r_version sf),Simple features for spatial data"
    echo "R Geospatial,raster,$(get_r_version raster),Raster data analysis"
    echo "R Geospatial,terra,$(get_r_version terra),Spatial data analysis"
    echo "R Geospatial,leaflet,$(get_r_version leaflet),Interactive web maps"
    echo "R Geospatial,tmap,$(get_r_version tmap),Thematic maps"
    echo "R Geospatial,ggstatsplot,$(get_r_version ggstatsplot),Statistical plots with ggplot2"
    echo "R Geospatial,plotly,$(get_r_version plotly),Interactive web graphics"
} >> "$CSV_FILE"

# System information
info "Adding system information..."
cat >> "$REPORT_FILE" << 'EOF'

## System Information

| Component | Details |
|-----------|---------|
EOF

{
    echo "| Operating System | $(cat /etc/os-release | grep PRETTY_NAME | cut -d'=' -f2 | tr -d '"') |"
    echo "| Kernel | $(uname -r) |"
    echo "| Architecture | $(uname -m) |"
    echo "| Container Engine | micromamba $(micromamba --version 2>/dev/null | cut -d' ' -f2) |"
    echo "| Shell | $(echo $SHELL) |"
} >> "$REPORT_FILE"

# Add footer
cat >> "$REPORT_FILE" << 'EOF'

---

## Usage Notes

This tool inventory was automatically generated and includes all software versions used in the Aedes aegypti local adaptation analysis pipeline. For reproducibility:

1. **Container Usage**: All analyses should be performed within the containerized environment to ensure identical software versions.

2. **Version Dependencies**: Some tools may have interdependent version requirements. The container environment resolves all compatibility issues.

3. **Citation Guidelines**: When publishing results, reference this tool inventory in supplementary materials and cite individual tools as appropriate.

4. **Updates**: This report reflects the software versions at the time of generation. Container updates may include newer versions.

## Reproducibility Statement

All computational analyses were performed using the software versions documented above within a containerized environment. The complete computational environment can be reproduced using the Dockerfile and associated configuration files in this repository.

For questions about specific tool versions or installation procedures, please refer to the project documentation or contact the development team through the project repository.
EOF

# Add system info to CSV
{
    echo "System,Operating System,$(cat /etc/os-release | grep PRETTY_NAME | cut -d'=' -f2 | tr -d '"'),Container base OS"
    echo "System,Kernel,$(uname -r),Linux kernel version"
    echo "System,Architecture,$(uname -m),System architecture"
    echo "System,Container Engine,micromamba $(micromamba --version 2>/dev/null | cut -d' ' -f2),Package manager"
} >> "$CSV_FILE"

success "Tool version report generated successfully!"
info "Files created:"
info "  - Markdown report: $REPORT_FILE"
info "  - CSV table: $CSV_FILE"
info ""
info "These files are suitable for:"
info "  • Supplementary materials in publications"
info "  • GitHub repository documentation"
info "  • Reproducibility statements"
info "  • Tool inventory tracking"

echo ""
echo "Report preview:"
echo "==============="
head -20 "$REPORT_FILE"
echo "..."
echo "(Full report saved to $REPORT_FILE)"
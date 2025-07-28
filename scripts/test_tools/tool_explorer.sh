#!/usr/bin/env bash
# Interactive Tool Explorer for Aedes Local Adaptation Container
# Provides user-friendly tool discovery and usage examples

set -uo pipefail

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# Functions for colored output
title() { echo -e "\n${BOLD}${PURPLE}=== $1 ===${NC}\n"; }
section() { echo -e "\n${CYAN}▶ $1${NC}"; }
info() { echo -e "${BLUE}ℹ${NC} $1"; }
success() { echo -e "${GREEN}✓${NC} $1"; }
example() { echo -e "${YELLOW}Example:${NC} $1"; }
cmd() { echo -e "  ${GREEN}\$${NC} $1"; }

# Tool categories with descriptions
declare -A CATEGORIES=(
    ["sequence"]="Sequence alignment and manipulation"
    ["variant"]="Variant calling and analysis"
    ["population"]="Population genetics analysis"
    ["adaptation"]="Local adaptation analysis"
    ["phylo"]="Phylogenetic analysis"
    ["stats"]="Statistical analysis"
    ["spatial"]="Spatial and geographic analysis"
    ["genomics"]="General genomics utilities"
    ["python"]="Python packages for genomics"
    ["r"]="R packages for analysis"
)

# Tool definitions with examples
declare -A TOOLS=(
    # Sequence tools
    ["bwa"]="sequence:BWA aligner for mapping reads to reference"
    ["samtools"]="sequence:SAM/BAM file manipulation"
    ["bcftools"]="variant:VCF/BCF file manipulation and variant calling"
    ["vcftools"]="variant:VCF file analysis and manipulation"
    ["bedtools"]="genomics:BED file manipulation and genomic arithmetic"
    
    # Population genetics
    ["plink"]="population:GWAS and population genetics analysis"
    ["plink2"]="population:Next-generation plink with improved algorithms"
    ["admixture"]="population:Ancestry estimation from SNP data"
    ["angsd"]="population:Analysis of next generation sequencing data"
    
    # Local adaptation
    ["bayescan"]="adaptation:Bayesian genome scan for selection"
    ["pcadapt"]="adaptation:PCA-based genome scan (R package)"
    ["OutFLANK"]="adaptation:F_ST outlier approach (R package)"
    ["lfmm"]="adaptation:Latent factor mixed models (R package)"
    ["LEA"]="adaptation:Landscape and ecological association (R package)"
    ["sambada"]="adaptation:Spatial analysis of molecular variance"
    
    # Phylogenetics
    ["iqtree"]="phylo:Maximum likelihood phylogenetic inference"
    
    # Python genomics
    ["allel"]="python:scikit-allel for population genetics"
    ["cyvcf2"]="python:Fast VCF parsing"
    ["pysam"]="python:Python interface to samtools"
    
    # R packages
    ["adegenet"]="r:Genetic data analysis"
    ["pcadapt"]="r:PCA-based selection scan"
    ["vegan"]="r:Community ecology analysis"
    ["sf"]="r:Spatial data analysis"
)

# Usage examples for key tools
declare -A EXAMPLES=(
    ["bwa"]="bwa index reference.fa && bwa mem reference.fa reads.fq > aligned.sam"
    ["samtools"]="samtools view -b aligned.sam > aligned.bam && samtools sort aligned.bam > sorted.bam"
    ["bcftools"]="bcftools call -mv -Oz -o variants.vcf.gz sorted.bam"
    ["vcftools"]="vcftools --gzvcf input.vcf.gz --maf 0.05 --recode --out filtered"
    ["plink"]="plink --vcf input.vcf --make-bed --out dataset"
    ["admixture"]="admixture dataset.bed 3 -j4"
    ["bayescan"]="bayescan -snp input.bayescan -threads 8"
    ["iqtree"]="iqtree -s alignment.phy -m MFP -bb 1000 -nt AUTO"
)

# Main menu
show_menu() {
    clear
    echo -e "${BOLD}${PURPLE}🧬 Aedes Local Adaptation Tool Explorer${NC}"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo
    echo "Select an option:"
    echo
    echo "  1) Browse tools by category"
    echo "  2) Search for a specific tool"
    echo "  3) Show tool usage examples"
    echo "  4) Test tool availability"
    echo "  5) Show common workflows"
    echo "  6) View Python package imports"
    echo "  7) View R package library"
    echo "  8) Quick reference guide"
    echo "  9) Exit"
    echo
    read -p "Enter your choice [1-9]: " choice
}

# Browse by category
browse_category() {
    clear
    title "Tool Categories"
    
    local i=1
    declare -a cat_array
    for cat in "${!CATEGORIES[@]}"; do
        echo "  $i) ${cat^} - ${CATEGORIES[$cat]}"
        cat_array[$i]=$cat
        ((i++))
    done
    echo "  0) Back to main menu"
    echo
    read -p "Select category [0-$((i-1))]: " cat_choice
    
    if [[ "$cat_choice" == "0" ]]; then
        return
    fi
    
    if [[ -n "${cat_array[$cat_choice]:-}" ]]; then
        show_category_tools "${cat_array[$cat_choice]}"
    fi
}

# Show tools in category
show_category_tools() {
    local category=$1
    clear
    title "${category^} Tools"
    info "${CATEGORIES[$category]}"
    echo
    
    for tool in "${!TOOLS[@]}"; do
        if [[ "${TOOLS[$tool]}" =~ ^${category}: ]]; then
            local desc="${TOOLS[$tool]#*:}"
            echo -e "${GREEN}$tool${NC} - $desc"
            if [[ -n "${EXAMPLES[$tool]:-}" ]]; then
                example "${EXAMPLES[$tool]}"
            fi
            echo
        fi
    done
    
    read -p "Press Enter to continue..."
}

# Search for tool
search_tool() {
    clear
    title "Tool Search"
    read -p "Enter tool name (partial match ok): " search_term
    echo
    
    local found=0
    for tool in "${!TOOLS[@]}"; do
        if [[ "$tool" =~ $search_term ]]; then
            local desc="${TOOLS[$tool]#*:}"
            echo -e "${GREEN}$tool${NC} - $desc"
            if command -v "$tool" &> /dev/null; then
                success "Available at: $(which $tool)"
            else
                info "Tool defined but not in PATH"
            fi
            if [[ -n "${EXAMPLES[$tool]:-}" ]]; then
                example "${EXAMPLES[$tool]}"
            fi
            echo
            found=1
        fi
    done
    
    if [[ $found -eq 0 ]]; then
        info "No tools found matching '$search_term'"
    fi
    
    read -p "Press Enter to continue..."
}

# Show examples
show_examples() {
    clear
    title "Tool Usage Examples"
    
    section "Basic Workflow Examples"
    echo
    
    info "1. Read alignment workflow:"
    cmd "bwa index reference.fa"
    cmd "bwa mem -t 8 reference.fa reads_R1.fq reads_R2.fq > aligned.sam"
    cmd "samtools view -b aligned.sam > aligned.bam"
    cmd "samtools sort -@ 8 aligned.bam > sorted.bam"
    cmd "samtools index sorted.bam"
    echo
    
    info "2. Variant calling workflow:"
    cmd "bcftools mpileup -Ou -f reference.fa sorted.bam | bcftools call -mv -Oz -o variants.vcf.gz"
    cmd "bcftools index variants.vcf.gz"
    cmd "vcftools --gzvcf variants.vcf.gz --maf 0.05 --max-missing 0.8 --recode --out filtered"
    echo
    
    info "3. Population structure analysis:"
    cmd "plink --vcf filtered.recode.vcf --make-bed --out dataset"
    cmd "plink --bfile dataset --pca --out pca_results"
    cmd "admixture dataset.bed 3 -j8 --cv"
    echo
    
    info "4. Selection scan workflow:"
    cmd "Rscript -e \"library(pcadapt); K <- read.pcadapt('dataset.bed'); x <- pcadapt(K, K=2)\""
    cmd "bayescan -snp dataset.bayescan -threads 8 -out_pilot -out_freq"
    echo
    
    read -p "Press Enter to continue..."
}

# Test tool availability
test_availability() {
    clear
    title "Testing Tool Availability"
    
    local available=0
    local total=0
    
    for tool in "${!TOOLS[@]}"; do
        ((total++))
        if command -v "$tool" &> /dev/null 2>&1; then
            success "$tool is available"
            ((available++))
        else
            # Check if it's an R package
            if [[ "${TOOLS[$tool]}" =~ ^r: ]]; then
                if R --slave -e "if(require('$tool', quietly=TRUE)) quit(status=0) else quit(status=1)" &> /dev/null; then
                    success "$tool (R package) is available"
                    ((available++))
                else
                    echo -e "${RED}✗${NC} $tool (R package) not found"
                fi
            # Check if it's a Python package
            elif [[ "${TOOLS[$tool]}" =~ ^python: ]]; then
                if python3 -c "import $tool" &> /dev/null 2>&1; then
                    success "$tool (Python package) is available"
                    ((available++))
                else
                    echo -e "${RED}✗${NC} $tool (Python package) not found"
                fi
            else
                echo -e "${RED}✗${NC} $tool not found in PATH"
            fi
        fi
    done
    
    echo
    info "Available: $available/$total tools"
    
    read -p "Press Enter to continue..."
}

# Show common workflows
show_workflows() {
    clear
    title "Common Analysis Workflows"
    
    section "1. Local Adaptation Analysis Pipeline"
    echo "   a) Prepare VCF file with quality filtering"
    echo "   b) Convert to required formats (PLINK, BayeScan, etc.)"
    echo "   c) Run environmental association tests"
    echo "   d) Perform selection scans"
    echo "   e) Integrate results"
    echo
    
    section "2. Population Structure Analysis"
    echo "   a) Filter variants for LD and MAF"
    echo "   b) Run PCA analysis"
    echo "   c) Estimate admixture proportions"
    echo "   d) Calculate FST between populations"
    echo "   e) Create visualization"
    echo
    
    section "3. Phylogenetic Analysis"
    echo "   a) Extract gene sequences"
    echo "   b) Perform multiple sequence alignment"
    echo "   c) Build phylogenetic tree"
    echo "   d) Test for selection on branches"
    echo
    
    read -p "Press Enter to continue..."
}

# Show Python packages
show_python() {
    clear
    title "Python Genomics Packages"
    
    python3 << 'EOF'
import sys
packages = [
    ('numpy', 'Numerical computing'),
    ('pandas', 'Data manipulation'),
    ('matplotlib', 'Plotting'),
    ('seaborn', 'Statistical visualization'),
    ('allel', 'Population genetics (scikit-allel)'),
    ('cyvcf2', 'Fast VCF parsing'),
    ('pysam', 'SAM/BAM manipulation'),
    ('Bio', 'Biopython - biological computation'),
    ('pyfaidx', 'FASTA indexing and retrieval')
]

for pkg, desc in packages:
    try:
        mod = __import__(pkg)
        version = getattr(mod, '__version__', 'unknown')
        print(f"✓ {pkg:<15} {version:<10} - {desc}")
    except ImportError:
        print(f"✗ {pkg:<15} {'N/A':<10} - {desc}")
EOF
    
    echo
    read -p "Press Enter to continue..."
}

# Show R packages
show_r_packages() {
    clear
    title "R Analysis Packages"
    
    R --slave << 'EOF'
packages <- list(
    c("ggplot2", "Data visualization"),
    c("tidyverse", "Data science toolkit"),
    c("adegenet", "Genetic data analysis"),
    c("pcadapt", "PCA-based selection scan"),
    c("OutFLANK", "FST outlier detection"),
    c("LEA", "Landscape genomics"),
    c("vegan", "Community ecology"),
    c("sf", "Spatial data analysis"),
    c("vcfR", "VCF file manipulation"),
    c("pegas", "Population genetics")
)

for (pkg_info in packages) {
    pkg <- pkg_info[1]
    desc <- pkg_info[2]
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
        version <- packageVersion(pkg)
        cat(sprintf("✓ %-15s %-10s - %s\n", pkg, version, desc))
    } else {
        cat(sprintf("✗ %-15s %-10s - %s\n", pkg, "N/A", desc))
    }
}
EOF
    
    echo
    read -p "Press Enter to continue..."
}

# Quick reference
quick_reference() {
    clear
    title "Quick Reference Guide"
    
    section "Essential Commands"
    cmd "samtools view -b -q 20 input.bam > filtered.bam  # Filter by quality"
    cmd "bcftools view -m2 -M2 -v snps input.vcf > snps.vcf  # Extract biallelic SNPs"
    cmd "vcftools --vcf input.vcf --maf 0.05 --max-missing 0.9 --recode  # Filter VCF"
    cmd "plink --vcf input.vcf --make-bed --out dataset  # Convert VCF to PLINK"
    echo
    
    section "File Conversions"
    cmd "plink --bfile dataset --recode vcf --out dataset  # PLINK to VCF"
    cmd "plink --vcf input.vcf --make-bed --out dataset  # VCF to PLINK"
    cmd "vcftools --vcf input.vcf --plink --out dataset  # VCF to PLINK (alt)"
    echo
    
    section "Quality Control"
    cmd "vcftools --vcf input.vcf --missing-indv  # Check missing data per sample"
    cmd "vcftools --vcf input.vcf --missing-site  # Check missing data per site"
    cmd "plink --bfile dataset --missing  # Generate missingness report"
    echo
    
    read -p "Press Enter to continue..."
}

# Main loop
while true; do
    show_menu
    case $choice in
        1) browse_category ;;
        2) search_tool ;;
        3) show_examples ;;
        4) test_availability ;;
        5) show_workflows ;;
        6) show_python ;;
        7) show_r_packages ;;
        8) quick_reference ;;
        9) echo "Goodbye!"; exit 0 ;;
        *) echo "Invalid choice. Press Enter to continue..."; read ;;
    esac
done
#!/usr/bin/env bash
# Quick Tool Help System - Get instant help for any tool
# Usage: ./tool_help.sh [tool_name]

set -uo pipefail

# Color definitions
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# Tool help database
declare -A TOOL_HELP=(
    ["samtools"]="SAM/BAM manipulation
Common commands:
  view    - Convert between SAM/BAM
  sort    - Sort alignment file
  index   - Index sorted BAM
  flagstat - Get alignment statistics
  depth   - Compute read depth
Example: samtools sort -@ 8 input.bam -o sorted.bam"

    ["bcftools"]="VCF/BCF manipulation
Common commands:
  view    - View/subset VCF
  call    - Call variants
  filter  - Filter variants
  merge   - Merge VCF files
  stats   - Get VCF statistics
Example: bcftools view -m2 -M2 -v snps input.vcf > biallelic_snps.vcf"

    ["vcftools"]="VCF analysis toolkit
Common options:
  --maf        - Minor allele frequency filter
  --max-missing - Filter by completeness
  --thin       - Thin sites by distance
  --hardy      - Hardy-Weinberg test
  --weir-fst   - Calculate FST
Example: vcftools --vcf input.vcf --maf 0.05 --recode --out filtered"

    ["plink"]="GWAS and population genetics
Common options:
  --make-bed   - Create binary fileset
  --pca        - Principal component analysis
  --genome     - IBD/IBS calculation
  --assoc      - Basic association test
  --distance   - Generate distance matrix
Example: plink --vcf input.vcf --make-bed --maf 0.05 --out dataset"

    ["bwa"]="Burrows-Wheeler Aligner
Commands:
  index   - Index reference genome
  mem     - Align with BWA-MEM algorithm
  aln     - Align with BWA-ALN (legacy)
Example: bwa mem -t 8 ref.fa reads_R1.fq reads_R2.fq > aligned.sam"

    ["admixture"]="Ancestry estimation
Usage: admixture [input.bed] [K]
Options:
  -j       - Number of threads
  --cv     - Cross-validation
  -s       - Random seed
Example: admixture dataset.bed 3 -j8 --cv"

    ["iqtree"]="Phylogenetic inference
Common options:
  -s       - Input alignment
  -m       - Substitution model
  -bb      - Bootstrap replicates
  -nt      - Number of threads
Example: iqtree -s alignment.phy -m GTR+G -bb 1000 -nt AUTO"

    ["bayescan"]="Bayesian selection scan
Options:
  -threads    - Number of threads
  -n          - Number of iterations
  -thin       - Thinning interval
  -nbp        - Number of pilot runs
Example: bayescan -snp input.bayescan -threads 8"
)

# Function to show tool help
show_help() {
    local tool="$1"
    
    if [[ -n "${TOOL_HELP[$tool]:-}" ]]; then
        echo -e "${BOLD}${GREEN}$tool${NC}"
        echo -e "${TOOL_HELP[$tool]}"
        echo
        
        # Check if tool is available
        if command -v "$tool" &> /dev/null; then
            echo -e "${CYAN}Location:${NC} $(which $tool)"
            echo -e "${CYAN}Version:${NC}"
            case "$tool" in
                samtools|bcftools) $tool --version | head -1 ;;
                vcftools) $tool --version ;;
                plink) $tool --version | head -1 ;;
                bwa) $tool 2>&1 | grep Version || echo "Version info not available" ;;
                *) echo "Version info not available" ;;
            esac
        else
            echo -e "${YELLOW}Note:${NC} Tool not found in PATH"
        fi
    else
        # Try to get help from the tool itself
        if command -v "$tool" &> /dev/null; then
            echo -e "${BOLD}${GREEN}$tool${NC}"
            echo -e "${CYAN}Location:${NC} $(which $tool)"
            echo
            echo "Getting help from tool..."
            echo
            
            case "$tool" in
                *.py) python3 "$tool" --help 2>&1 | head -20 ;;
                *.R|*.r) Rscript "$tool" --help 2>&1 | head -20 ;;
                *) $tool --help 2>&1 | head -20 || $tool -h 2>&1 | head -20 || echo "No help available" ;;
            esac
        else
            echo "Tool '$tool' not found or no help available"
            echo
            echo "Try searching with: ./tool_explorer.sh"
        fi
    fi
}

# Main
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <tool_name>"
    echo
    echo "Get quick help for bioinformatics tools"
    echo
    echo "Examples:"
    echo "  $0 samtools"
    echo "  $0 vcftools"
    echo "  $0 plink"
    echo
    echo "For interactive exploration, use: ./tool_explorer.sh"
else
    show_help "$1"
fi
#!/usr/bin/env bash

# Start Jupyter Lab for Aedes Local Adaptation Analysis
# Ensure micromamba environment is properly activated

# Initialize micromamba environment if available
if command -v micromamba >/dev/null 2>&1; then
    export MAMBA_ROOT_PREFIX=/opt/conda
    export PATH="/opt/conda/bin:$PATH"
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate base 2>/dev/null || true
fi

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

info "Starting Jupyter Lab for Aedes Local Adaptation Analysis..."

# Check if Jupyter is available
if ! command -v jupyter >/dev/null 2>&1; then
    error "Jupyter not found. This should be pre-installed in the container."
    echo "  If running outside container, install with:"
    echo "  micromamba install -c conda-forge jupyter jupyterlab"
    echo "  Or use the containerized environment."
    exit 1
fi

# Set Jupyter configuration
export JUPYTER_CONFIG_DIR="/proj/.jupyter"
export JUPYTER_DATA_DIR="/proj/.jupyter"
export JUPYTER_RUNTIME_DIR="/proj/.jupyter/runtime"

# Create Jupyter directories
mkdir -p "${JUPYTER_CONFIG_DIR}" "${JUPYTER_DATA_DIR}" "${JUPYTER_RUNTIME_DIR}"

# Generate Jupyter config if it doesn't exist
if [[ ! -f "${JUPYTER_CONFIG_DIR}/jupyter_lab_config.py" ]]; then
    info "Generating Jupyter Lab configuration..."
    jupyter lab --generate-config

    # Configure Jupyter Lab settings
    cat >> "${JUPYTER_CONFIG_DIR}/jupyter_lab_config.py" <<EOF

# Aedes Local Adaptation Analysis Jupyter Configuration
c.ServerApp.ip = '0.0.0.0'
c.ServerApp.port = 8888
c.ServerApp.open_browser = False
c.ServerApp.allow_root = True
c.ServerApp.token = ''
c.ServerApp.password = ''
c.ServerApp.notebook_dir = '/proj'
c.ServerApp.allow_origin = '*'
c.ServerApp.allow_remote_access = True

# Disable security warnings for development
c.ServerApp.disable_check_xsrf = True

# Increase max buffer size for large outputs
c.ServerApp.max_buffer_size = 1073741824  # 1GB

# Enable extensions
c.LabApp.check_for_updates_class = 'jupyterlab.NeverCheckForUpdate'
EOF
    success "Jupyter Lab configuration created"
fi

# Set environment variables for R kernel
export R_LIBS_USER="/opt/conda/lib/R/library"

# Start Jupyter Lab
info "Starting Jupyter Lab..."
info "Access URL: http://localhost:8888"
info "Working directory: /proj"
info ""
info "Available kernels:"
echo "  - Python 3 (ipykernel) - Population genomics analysis"
echo "  - R (IRkernel) - Local adaptation detection"
echo ""
info "Pre-installed packages:"
echo "  Python: pandas, numpy, matplotlib, seaborn, scikit-allel, cyvcf2, pysam"
echo "  R: tidyverse, OutFLANK, pcadapt, adegenet, vcfR, R.SamBada"
echo ""
info "Analysis workflows available:"
echo "  - Population structure analysis (PCA, ADMIXTURE)"
echo "  - Local adaptation detection (OutFLANK, pcadapt)"
echo "  - Environmental association (R.SamBada)"
echo "  - VCF processing and filtering"
echo "  - Visualization and plotting"
echo ""
info "Example notebooks:"
echo "  - Check scripts/ directory for analysis examples"
echo "  - Data should be placed in data/raw/"
echo "  - Results will be saved to results/analysis/"
echo ""
warning "Press Ctrl+C to stop Jupyter Lab"
echo ""

# Change to project directory
cd /proj

# Start Jupyter Lab
jupyter lab \
    --config="${JUPYTER_CONFIG_DIR}/jupyter_lab_config.py" \
    --no-browser \
    --allow-root \
    --port=8888 \
    --ip=0.0.0.0

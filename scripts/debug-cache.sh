#!/bin/bash
# Cache debugging helper script

echo "=== Docker Cache Debug Info ==="
echo "Date: $(date)"
echo "Cache bust value: ${CACHE_BUST:-not set}"
echo ""

echo "=== Checking installed tools ==="
for tool in bcftools vcftools plink admixtools bayescan baypass ba3 pcadapt; do
    if command -v $tool &> /dev/null; then
        echo "✓ $tool found at: $(which $tool)"
    else
        echo "✗ $tool NOT FOUND"
    fi
done

echo ""
echo "=== Python packages ==="
python -c "import numpy; print(f'numpy: {numpy.__version__}')" 2>/dev/null || echo "✗ numpy not found"
python -c "import pandas; print(f'pandas: {pandas.__version__}')" 2>/dev/null || echo "✗ pandas not found"
python -c "import allel; print(f'scikit-allel: {allel.__version__}')" 2>/dev/null || echo "✗ scikit-allel not found"

echo ""
echo "=== Conda environments ==="
micromamba env list

echo ""
echo "=== FastStructure3 check ==="
if [ -d "/opt/fastStructure3" ]; then
    echo "✓ FastStructure3 directory exists"
    if [ -f "/usr/local/bin/fastStructure" ]; then
        echo "✓ FastStructure wrapper script exists"
    else
        echo "✗ FastStructure wrapper script missing"
    fi
else
    echo "✗ FastStructure3 directory missing"
fi

echo ""
echo "=== If you see missing tools that should be installed ==="
echo "This is likely a cache issue. Solutions:"
echo "1. Trigger a manual workflow with 'Force rebuild without cache' checked"
echo "2. Add 'RUN echo \"bust cache \$(date)\"' before the problematic layer"
echo "3. Change CACHE_DATE build arg in the workflow"
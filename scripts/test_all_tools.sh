#!/usr/bin/env bash
# Master Test Script - Aedes Local Adaptation Container
# Runs all tool tests and provides comprehensive summary
set -uo pipefail

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

# Script directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
LOGS_DIR="$PROJECT_DIR/logs"

echo "🧬 Aedes Local Adaptation Container - Master Test Suite"
echo "======================================================"
echo "Running comprehensive tool testing for bioinformatics container"
echo "Date: $(date)"
echo "Container: $(hostname)"
echo ""

# Initialize counters
TOTAL_MAIN_TOOLS=0
AVAILABLE_MAIN_TOOLS=0
TOTAL_GDAL_TOOLS=0
AVAILABLE_GDAL_TOOLS=0
SEGFAULT_GDAL_TOOLS=0

# Test 1: Main container tools
section "TEST 1: Main Bioinformatics Tools"
echo "Running primary tool test suite..."
echo ""

if [ -f "$SCRIPT_DIR/test_tools/test_container_tools.sh" ]; then
    # Create logs directory if it doesn't exist
    mkdir -p "$LOGS_DIR"

    # Run main test and capture output
    if bash "$SCRIPT_DIR/test_tools/test_container_tools.sh" > "$LOGS_DIR/main_test.log" 2>&1; then
        success "Main tools test completed"
    else
        warning "Main tools test completed with warnings"
    fi

    # Parse results from main test
    if [ -f "$LOGS_DIR/main_test.log" ]; then
        # Remove ANSI color codes before parsing
        TOTAL_MAIN_TOOLS=$(grep "Total tools tested:" "$LOGS_DIR/main_test.log" | sed 's/\x1b\[[0-9;]*m//g' | awk '{print $4}' 2>/dev/null || echo "0")
        AVAILABLE_MAIN_TOOLS=$(grep "Available tools:" "$LOGS_DIR/main_test.log" | sed 's/\x1b\[[0-9;]*m//g' | awk '{print $4}' 2>/dev/null || echo "0")

        # Ensure variables are numeric
        TOTAL_MAIN_TOOLS=${TOTAL_MAIN_TOOLS//[^0-9]/}
        AVAILABLE_MAIN_TOOLS=${AVAILABLE_MAIN_TOOLS//[^0-9]/}
        TOTAL_MAIN_TOOLS=${TOTAL_MAIN_TOOLS:-0}
        AVAILABLE_MAIN_TOOLS=${AVAILABLE_MAIN_TOOLS:-0}

        # Show summary from main test
        echo "Main tools summary:"
        grep -E "(Total tools tested|Available tools|Missing tools|Success rate)" "$LOGS_DIR/main_test.log" || echo "Could not parse main test results"
    fi
else
    error "Main test script not found: $SCRIPT_DIR/test_tools/test_container_tools.sh"
fi

echo ""
echo "────────────────────────────────────────────────────────"
echo ""

# Test 2: GDAL and geospatial tools
section "TEST 2: GDAL and Geospatial Tools"
echo "Running GDAL/geospatial tool test suite..."
echo ""

if [ -f "$SCRIPT_DIR/test_tools/test_gdal_tools.sh" ]; then
        # Run GDAL test and capture output
    if bash "$SCRIPT_DIR/test_tools/test_gdal_tools.sh" > "$LOGS_DIR/gdal_test.log" 2>&1; then
        success "GDAL tools test completed"
    else
        warning "GDAL tools test completed with warnings"
    fi

    # Parse results from GDAL test
    if [ -f "$LOGS_DIR/gdal_test.log" ]; then
        # Remove ANSI color codes before parsing
        TOTAL_GDAL_TOOLS=$(grep "Total tools tested:" "$LOGS_DIR/gdal_test.log" | sed 's/\x1b\[[0-9;]*m//g' | awk '{print $4}' 2>/dev/null || echo "0")
        AVAILABLE_GDAL_TOOLS=$(grep "Working tools:" "$LOGS_DIR/gdal_test.log" | sed 's/\x1b\[[0-9;]*m//g' | awk '{print $4}' 2>/dev/null || echo "0")

        # Ensure variables are numeric
        TOTAL_GDAL_TOOLS=${TOTAL_GDAL_TOOLS//[^0-9]/}
        AVAILABLE_GDAL_TOOLS=${AVAILABLE_GDAL_TOOLS//[^0-9]/}
        TOTAL_GDAL_TOOLS=${TOTAL_GDAL_TOOLS:-0}
        AVAILABLE_GDAL_TOOLS=${AVAILABLE_GDAL_TOOLS:-0}

        # Handle segfaulting tools line which might not exist
        if grep -q "Segfaulting tools:" "$LOGS_DIR/gdal_test.log"; then
            SEGFAULT_GDAL_TOOLS=$(grep "Segfaulting tools:" "$LOGS_DIR/gdal_test.log" | awk '{print $3}' 2>/dev/null || echo "0")
            SEGFAULT_GDAL_TOOLS=${SEGFAULT_GDAL_TOOLS//[^0-9]/}
            SEGFAULT_GDAL_TOOLS=${SEGFAULT_GDAL_TOOLS:-0}
        else
            SEGFAULT_GDAL_TOOLS=0
        fi

        # Show summary from GDAL test
        echo "GDAL tools summary:"
        grep -E "(Total tools tested|Working tools|Segfaulting tools|Tools present)" "$LOGS_DIR/gdal_test.log" || echo "Could not parse GDAL test results"
    fi
else
    warning "GDAL test script not found: $SCRIPT_DIR/test_tools/test_gdal_tools.sh"
    warning "Skipping GDAL tests (this is optional)"
fi

echo ""
echo "══════════════════════════════════════════════════════════"
echo ""

# Overall Summary
section "OVERALL CONTAINER SUMMARY"
echo "Master test suite completed at $(date)"
echo ""

# Set defaults if variables are not set
TOTAL_MAIN_TOOLS=${TOTAL_MAIN_TOOLS:-0}
AVAILABLE_MAIN_TOOLS=${AVAILABLE_MAIN_TOOLS:-0}
TOTAL_GDAL_TOOLS=${TOTAL_GDAL_TOOLS:-0}
AVAILABLE_GDAL_TOOLS=${AVAILABLE_GDAL_TOOLS:-0}
SEGFAULT_GDAL_TOOLS=${SEGFAULT_GDAL_TOOLS:-0}

# Calculate totals
TOTAL_ALL_TOOLS=$((TOTAL_MAIN_TOOLS + TOTAL_GDAL_TOOLS))
AVAILABLE_ALL_TOOLS=$((AVAILABLE_MAIN_TOOLS + AVAILABLE_GDAL_TOOLS))

echo "📊 Test Results:"
echo "  Main bioinformatics tools: $AVAILABLE_MAIN_TOOLS/$TOTAL_MAIN_TOOLS"
echo "  GDAL/geospatial tools: $AVAILABLE_GDAL_TOOLS/$TOTAL_GDAL_TOOLS"
if [ $SEGFAULT_GDAL_TOOLS -gt 0 ]; then
    echo "  GDAL segfaulting tools: $SEGFAULT_GDAL_TOOLS (installed but unstable)"
fi
echo "  ─────────────────────────────────────"
echo "  Total working tools: $AVAILABLE_ALL_TOOLS/$TOTAL_ALL_TOOLS"

# Calculate overall percentage
if [ $TOTAL_ALL_TOOLS -gt 0 ]; then
    OVERALL_PERCENTAGE=$((AVAILABLE_ALL_TOOLS * 100 / TOTAL_ALL_TOOLS))
    echo "  Overall success rate: ${OVERALL_PERCENTAGE}%"
fi

echo ""

# Container readiness assessment
MAIN_PERCENTAGE=0
if [ $TOTAL_MAIN_TOOLS -gt 0 ]; then
    MAIN_PERCENTAGE=$((AVAILABLE_MAIN_TOOLS * 100 / TOTAL_MAIN_TOOLS))

    if [ $MAIN_PERCENTAGE -ge 90 ]; then
        success "🎉 Container is READY for production use!"
        success "   Main bioinformatics tools: ${MAIN_PERCENTAGE}% functional"
    elif [ $MAIN_PERCENTAGE -ge 75 ]; then
        warning "⚠️  Container is MOSTLY READY with minor issues"
        warning "   Main bioinformatics tools: ${MAIN_PERCENTAGE}% functional"
    else
        error "❌ Container has SIGNIFICANT ISSUES"
        error "   Main bioinformatics tools: only ${MAIN_PERCENTAGE}% functional"
    fi
else
    error "❌ Unable to assess container readiness - test failure"
fi

echo ""
info "📋 Container Information:"
echo "  - Architecture: $(uname -m)"
echo "  - Kernel: $(uname -r)"
echo "  - Container OS: $(cat /etc/os-release | grep PRETTY_NAME | cut -d'=' -f2 | tr -d '"' || echo 'Unknown')"
echo "  - Python: $(python3 --version 2>/dev/null || echo 'Not found')"
echo "  - R: $(R --version 2>/dev/null | head -1 || echo 'Not found')"
echo "  - Working directory: $(pwd)"
echo "  - User: $(whoami)"

echo ""
info "📁 Log files generated:"
if [ -f "$LOGS_DIR/main_test.log" ]; then
    echo "  - Main tools: $LOGS_DIR/main_test.log"
fi
if [ -f "$LOGS_DIR/gdal_test.log" ]; then
    echo "  - GDAL tools: $LOGS_DIR/gdal_test.log"
fi

echo ""
info "🚀 Next steps:"
echo "  1. If success rate >90%: Container ready for HPC deployment"
echo "  2. If success rate 75-90%: Review warnings, container usable"
echo "  3. If success rate <75%: Check container build logs"
echo "  4. For GDAL issues: Expected due to library conflicts"

# Save summary to file
{
    echo "# Container Test Summary - $(date)"
    echo "Main tools: $AVAILABLE_MAIN_TOOLS/$TOTAL_MAIN_TOOLS"
    echo "GDAL tools: $AVAILABLE_GDAL_TOOLS/$TOTAL_GDAL_TOOLS"
    echo "Total: $AVAILABLE_ALL_TOOLS/$TOTAL_ALL_TOOLS"
    if [ $TOTAL_ALL_TOOLS -gt 0 ]; then
        echo "Success rate: $((AVAILABLE_ALL_TOOLS * 100 / TOTAL_ALL_TOOLS))%"
    fi
    echo "Status: $(if [ $MAIN_PERCENTAGE -ge 90 ]; then echo "READY"; elif [ $MAIN_PERCENTAGE -ge 75 ]; then echo "MOSTLY_READY"; else echo "ISSUES"; fi)"
} > "$LOGS_DIR/container_test_summary.txt"

echo ""
success "📄 Summary saved to: $LOGS_DIR/container_test_summary.txt"
echo ""

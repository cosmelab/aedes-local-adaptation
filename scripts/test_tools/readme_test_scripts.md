# Container Testing Scripts

This directory (`scripts/test_tools/`) contains comprehensive testing scripts for validating the Aedes Local Adaptation bioinformatics container.

## Overview

The container includes 100+ bioinformatics tools and packages. These scripts systematically test tool availability and functionality to ensure the container is ready for HPC deployment.

## Test Scripts

### 1. `test_container_tools.sh`
**Purpose**: Tests the main bioinformatics tools and packages

**Coverage**:
- System and package managers (micromamba, pip, git, etc.)
- Shell and terminal tools (zsh, lsd, colorls, starship)
- Python core and bioinformatics packages
- Genomics tools (samtools, bcftools, vcftools, plink, etc.)
- Analysis tools (admixture, iqtree, snakemake)
- Jupyter ecosystem
- R packages (300+ packages)
- Local adaptation tools (sambada, AdmixTools, BayeScan, etc.)

**Usage**:
```bash
# From project root:
bash scripts/test_all_tools.sh

# Or run individual test:
bash scripts/test_tools/test_container_tools.sh
```

### 2. `test_gdal_tools.sh`
**Purpose**: Tests GDAL and geospatial tools separately due to known segfault issues

**Coverage**:
- GDAL core tools (gdalinfo, gdal_translate, gdalwarp, ogr2ogr)
- PROJ tools (proj, cs2cs, projinfo)
- Python geospatial packages (rasterio, fiona, geopandas, pyproj, shapely)
- R spatial packages (sf, raster)

**Special Features**:
- Timeout protection to handle segfaults gracefully
- Distinguishes between missing tools and segfaulting tools
- Provides context about expected GDAL/PROJ version conflicts

**Usage**:
```bash
# From project root:
bash scripts/test_all_tools.sh

# Or run individual test:
bash scripts/test_tools/test_gdal_tools.sh
```

### 3. `test_all_tools.sh` (Main Script)
**Purpose**: Master wrapper script located in `scripts/` that runs all tests and provides comprehensive summary

**Location**: `scripts/test_all_tools.sh` (main directory)

**Features**:
- Runs both main and GDAL test suites from the test_tools directory
- Aggregates results from all tests
- Provides overall container readiness assessment
- Generates summary files in the main logs directory
- Color-coded status reporting

**Usage**:
```bash
# This is the main entry point - run from project root:
bash scripts/test_all_tools.sh
```

## Running Tests

### On HPC with Singularity

```bash
# Pull the latest container
singularity pull aedes-local-adaptation.sif docker://ghcr.io/cosmelab/aedes-local-adaptation:latest

# Run all tests
singularity exec aedes-local-adaptation.sif bash /proj/scripts/test_all_tools.sh

# Or run interactively
singularity shell --cleanenv --bind $PWD:/proj aedes-local-adaptation.sif
cd /proj
bash scripts/test_all_tools.sh
```

### With Docker

```bash
# Run all tests
docker run --rm -v $PWD:/proj ghcr.io/cosmelab/aedes-local-adaptation:latest \
    bash /proj/scripts/test_all_tools.sh

# Or run interactively
docker run -it --rm -v $PWD:/proj ghcr.io/cosmelab/aedes-local-adaptation:latest /bin/zsh
cd /proj
bash scripts/test_all_tools.sh
```

## Understanding Results

### Success Indicators
- **>90% tools available**: Container is ready for production use
- **75-90% tools available**: Container is usable with minor limitations
- **<75% tools available**: Container has significant issues

### Expected Issues
- **GDAL/PROJ segfaults**: Known issue due to library version conflicts. Tools are installed but may be unstable for certain operations.
- **Missing optional packages**: Some specialized packages may not be available but don't affect core functionality.

### Output Files
- `logs/main_test.log`: Detailed results from main tool testing
- `logs/gdal_test.log`: Detailed results from GDAL tool testing
- `logs/container_test_summary.txt`: Concise summary of all tests

Note: All log files are saved to the main `logs/` directory in the project root, not in the scripts directory.

## Troubleshooting

### Common Issues

1. **"Micromamba not found"**
   - Ensure you're running inside the container
   - Check that the container image is correct

2. **Many tools showing as "NOT FOUND"**
   - Verify container was built successfully
   - Check GitHub Actions logs for build errors
   - Ensure you're using the latest container version

3. **GDAL tools segfaulting**
   - This is expected due to library conflicts
   - Core bioinformatics functionality is not affected
   - Use alternative tools when possible

### Getting Help

If you encounter issues:
1. Check the container build logs in GitHub Actions
2. Review the Dockerfile for expected tool installations
3. Run individual test scripts to isolate problems
4. Check `logs/*.log` files for detailed error messages

## Container Validation Workflow

1. **Build container** (GitHub Actions)
2. **Pull to HPC** (Singularity)
3. **Run test suite** (`test_all_tools.sh`)
4. **Review results** (check success rate)
5. **Deploy if ready** (>90% success rate)

## Notes

- Tests are non-destructive and read-only
- No data is modified during testing
- Tests can be run multiple times safely
- Results may vary slightly between runs due to network timeouts
- GDAL segfaults do not affect other tools
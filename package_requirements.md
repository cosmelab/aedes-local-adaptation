# Aedes Local Adaptation Analysis - Package Requirements

## Core System Tools
- **zsh** - Shell environment
- **python** - Python 3.x
- **R** - R 4.3.2+ (rocker/tidyverse base)
- **conda/mamba** - Package management

## Bioinformatics Tools
- **PLINK2** - Genotype data manipulation
- **bcftools** - VCF file manipulation
- **fastStructure** - Population structure analysis
- **Admixture** - Ancestry analysis
- **IQ-TREE** - Phylogenetic analysis
- **BA3-SNPs** - Migration analysis
- **SamBada** - Local adaptation analysis

## R Packages

### Core Data Manipulation
- tidyverse
- data.table
- here
- dplyr
- readr
- purrr
- stringr

### Visualization
- ggplot2
- scales
- ggrepel
- ggtext
- ggvenn
- ggstatsplot
- ggforce
- ggpattern
- scatterpie
- RColorBrewer
- extrafont
- forcats
- flextable
- officer
- Cairo

### Bioinformatics
- vcfR
- adegenet
- dartR
- Biostrings
- SNPRelate
- biomaRt

### Population Genetics
- OutFLANK
- pcadapt
- R.SamBada
- LEA
- qvalue

### Spatial Analysis
- raster
- sf
- rgdal
- geosphere
- rnaturalearth
- rnaturalearthdata
- ggspatial
- fields

### Development
- devtools
- remotes
- reticulate

### Other
- colorout
- grid
- ellipse
- reshape2
- metafor

## Python Packages
- numpy
- scipy
- pandas
- matplotlib
- seaborn
- bokeh
- altair
- holoviews
- cython
- gsl (GNU Scientific Library)

## System Dependencies
- **GDAL** - Required for R.SamBada spatial analysis
- **PROJ** - Coordinate reference system library
- **SQLite3** - Database support
- **GCC7** - Required for SamBada compilation
- **GSL** - GNU Scientific Library
- **Development libraries**: libcurl4-openssl-dev, libssl-dev, libxml2-dev, libfontconfig1-dev, libharfbuzz-dev, libfribidi-dev, libfreetype6-dev, libpng-dev, libtiff5-dev, libjpeg-dev

## GitHub Packages (via devtools/remotes)
- R.SamBada (SolangeD/R.SamBada)
- gdalUtils (gearslaboratory/gdalUtils)
- templater (petrelharp/templater)
- lostruct (petrelharp/local_pca/lostruct)
- TESS3_encho_sen (eriqande/TESS3_encho_sen)
- genoscapeRtools (eriqande/genoscapeRtools)

## Bioconductor Packages
- SNPRelate
- biomaRt
- LEA
- AnnotationDbi
- DESeq2 (if needed for RNA-seq)
- tximport (if needed for RNA-seq)
- sva (if needed for RNA-seq)
- edgeR (if needed for RNA-seq)
- limma (if needed for RNA-seq) 
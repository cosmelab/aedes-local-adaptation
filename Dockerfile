# HPC-Optimized Container for Aedes Local Adaptation Analysis
# Architecture: AMD64/x86_64 only
# Target: HPC clusters and compute nodes

FROM mambaorg/micromamba:1.5.0

# Document architecture requirement
LABEL maintainer="cosmelab@domain.com" \
      architecture="amd64" \
      description="HPC-optimized bioinformatics environment for local adaptation analysis"

ENV MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:/usr/bin:$PATH \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# HPC Performance: Prevent tools from auto-detecting all cores
ENV OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1

SHELL ["bash", "-lc"]
USER root

# Update micromamba and all packages to latest versions
RUN micromamba update --all -y && micromamba clean --all --yes

# Install system dependencies and create user in a single layer
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    ca-certificates \
    curl \
    git \
    unzip \
    zsh \
    libcairo2-dev \
    libbz2-dev \
    liblzma-dev \
    wget \
    lsb-release \
    libz-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libtinfo5 \
    gdal-bin \
    libgdal-dev \
    libproj-dev \
    proj-data \
    sqlite3 \
    libsqlite3-dev \
    libgsl-dev \
    libgslcblas0 \
    libopenblas-dev \
    autojump \
    imagemagick \
    libmagick++-dev \
    libmagickwand-dev \
    hwloc \
    numactl \
    libv8-dev \
    libpoppler-cpp-dev \
    jags \
    cargo \
    rustc \
    && apt-get clean && rm -rf /var/lib/apt/lists/* \
    && groupadd -r aedes && useradd -r -g aedes aedes \
    && mkdir -p /home/aedes && chown aedes:aedes /home/aedes

# Install lsd with error checking and cleanup
RUN wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
    https://github.com/lsd-rs/lsd/releases/download/v1.1.5/lsd-v1.1.5-x86_64-unknown-linux-gnu.tar.gz \
    && tar -xzf lsd-v1.1.5-x86_64-unknown-linux-gnu.tar.gz \
    && mv lsd-v1.1.5-x86_64-unknown-linux-gnu/lsd /usr/local/bin/ \
    && rm -rf lsd-v1.1.5-* \
    && chmod +x /usr/local/bin/lsd

# Install core Python and system packages (combined for better caching)
RUN micromamba install --channel-priority strict -c conda-forge -c bioconda \
    libstdcxx-ng \
    python=3.11.7 \
    pip \
    gcc \
    make \
    cmake \
    gsl \
    cython \
    numpy \
    scipy \
    pandas \
    matplotlib \
    seaborn \
    scikit-allel \
    cyvcf2 \
    pyfaidx \
    pysam \
    biopython \
    scikit-learn \
    plotly \
    bokeh \
    altair \
    holoviews \
    networkx \
    geopandas \
    fiona \
    rasterio \
    pyproj \
    shapely \
    folium \
    contextily \
    earthpy \
    geoplot \
    -y && micromamba clean --all --yes && \
    # Verify critical packages are installed
    python3 -c "import pysam; import allel; import cyvcf2; print('Critical Python packages verified')"

# Install genomics tools
RUN micromamba install --channel-priority strict -c conda-forge -c bioconda \
    samtools \
    bcftools \
    vcftools \
    bedtools \
    htslib \
    tabix \
    bwa \
    plink \
    plink2 \
    angsd \
    admixture \
    iqtree \
    starship \
    datamash \
    openjdk=17 \
    sra-tools \
    entrez-direct \
    snakemake \
    jupyter \
    jupyterlab \
    notebook \
    ipykernel \
    -y && micromamba clean --all --yes

# Install R and R packages (combined for efficiency)
RUN micromamba install --channel-priority strict -c conda-forge -c bioconda \
    r-base=4.3.2 \
    r-devtools \
    r-tidyverse \
    r-ggplot2 \
    r-here \
    r-data.table \
    r-ade4 \
    r-mass \
    r-vegan \
    r-seqinr \
    r-qqconf \
    r-sf \
    r-raster \
    r-adegenet \
    r-pcadapt \
    r-circlize \
    r-igraph \
    r-plotly \
    r-leaflet \
    r-networkd3 \
    r-tmap \
    r-v8 \
    bioconductor-snprelate \
    bioconductor-annotationdbi \
    bioconductor-biomart \
    bioconductor-biostrings \
    ruby=3.2.2 \
    -y && micromamba clean --all --yes

# Install R packages for population genetics and local adaptation analysis
# Note: Many packages are already installed via conda-forge above
RUN R -e "install.packages(c('qqman', 'qqplotr', 'reticulate', 'broom', \
    'readxl', 'writexl', 'knitr', 'rmarkdown', 'pegas', 'ape', 'phangorn', 'vcfR', 'genetics', \
    'BiocManager', 'remotes', 'scales', 'ggrepel', 'ggtext', 'ggvenn', 'ggstatsplot', 'ggforce', \
    'ggpattern', 'scatterpie', 'RColorBrewer', 'extrafont', 'forcats', 'flextable', 'officer', \
    'Cairo', 'geosphere', 'rnaturalearth', 'rnaturalearthdata', 'ggspatial', \
    'fields', 'ellipse', 'reshape2'), \
    repos='https://cloud.r-project.org/', dependencies=TRUE, Ncpus=4)"

# Install OutFLANK dependencies first, then OutFLANK from GitHub for better reliability
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(c('qvalue'), update = FALSE, ask = FALSE)" && \
    R -e "remotes::install_github('whitlock/OutFLANK')"

# Install Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(c('LEA'), update = FALSE, ask = FALSE)"

# Install GitHub R packages with error handling
RUN R -e "options(Ncpus = 4); \
    devtools::install_github('uqrmaie1/admixtools'); \
    devtools::install_github('gearslaboratory/gdalUtils', force = TRUE); \
    devtools::install_github('petrelharp/templater'); \
    devtools::install_github('petrelharp/local_pca/lostruct'); \
    remotes::install_github('eriqande/TESS3_encho_sen'); \
    remotes::install_github('eriqande/genoscapeRtools'); \
    devtools::install_github('SolangeD/R.SamBada', build_vignettes = FALSE); \
    devtools::install_github('bcm-uga/lfmm'); \
    devtools::install_github('bcm-uga/TESS3_encho_sen'); \
    devtools::install_github('bcm-uga/LEA'); \
    devtools::install_github('petrikemppainen/LDna', ref = 'v.2.15')"

# Install Python packages not in conda-forge
RUN pip3 install --no-cache-dir pong

# Install colorls globally as root
RUN /opt/conda/bin/gem install colorls --no-document -n /usr/local/bin && \
    chmod +x /usr/local/bin/colorls

# Install additional system dependencies for compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    autoconf \
    automake \
    libtool \
    build-essential \
    gfortran \
    libgsl-dev \
    liblapack-dev \
    libblas-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Download and install local adaptation tools with proper error handling
# Split into separate steps for better debugging
RUN R -e "library(R.SamBada); downloadSambada('/opt/sambada')" && \
    find /opt/sambada -name "sambada*" -type f -executable -exec chmod +x {} \; && \
    find /opt/sambada -name "sambada*" -type f -executable -exec ln -sf {} /usr/local/bin/ \;

# Install AdmixTools with OpenBLAS linking fix and function declaration patch
RUN cd /opt && \
    git clone --depth=1 https://github.com/DReichLab/AdmixTools.git && \
    cd AdmixTools && cd src && \
    export LDFLAGS="-L/opt/conda/lib" && \
    export CFLAGS="-I/opt/conda/include" && \
    sed -i 's/-lopenblas/-lblas/g' Makefile && \
    sed -i 's/void setgtime ();/void setgtime(double *time);/' qpgsubs.c && \
    make && make install && cd .. && \
    ln -sf /opt/AdmixTools/bin/* /usr/local/bin/ && \
    rm -rf /opt/AdmixTools/.git

# Install other binary tools  
RUN cd /opt && \
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip && \
    unzip BayeScan2.1.zip && \
    chmod +x BayeScan2.1/binaries/BayeScan2.1_linux64bits && \
    ln -sf /opt/BayeScan2.1/binaries/BayeScan2.1_linux64bits /usr/local/bin/bayescan && \
    rm -f BayeScan2.1.zip && \
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.4/gemma-0.98.4-linux-static-AMD64.gz && \
    gunzip gemma-0.98.4-linux-static-AMD64.gz && \
    chmod +x gemma-0.98.4-linux-static-AMD64 && \
    mv gemma-0.98.4-linux-static-AMD64 /usr/local/bin/gemma && \
    git clone --depth=1 https://github.com/stevemussmann/BayesAss3-SNPs.git && \
    cd BayesAss3-SNPs && \
    mv BA3-SNPS-Ubuntu64 BA3-SNPS && \
    chmod +x BA3-SNPS && \
    ln -sf /opt/BayesAss3-SNPs/BA3-SNPS /usr/local/bin/BA3-SNPS && \
    rm -rf /opt/BayesAss3-SNPs/.git && \
    cd /opt && \
    git clone --depth=1 https://github.com/stevemussmann/BA3-SNPS-autotune.git && \
    cd BA3-SNPS-autotune && \
    chmod +x BA3-SNPS-autotune.py && \
    ln -sf /opt/BA3-SNPS-autotune/BA3-SNPS-autotune.py /usr/local/bin/BA3-SNPS-autotune && \
    rm -rf /opt/BA3-SNPS-autotune/.git

# Create HPC module-compatible activation script
RUN echo '#!/bin/zsh' > /opt/conda/bin/activate-env.sh && \
    echo 'eval "$(micromamba shell hook --shell zsh)"' >> /opt/conda/bin/activate-env.sh && \
    echo 'micromamba activate base' >> /opt/conda/bin/activate-env.sh && \
    chmod +x /opt/conda/bin/activate-env.sh

# Switch to aedes user
USER aedes

# Install and configure Oh-My-Zsh, Powerlevel10k, and all plugins
RUN git clone --depth=1 https://github.com/romkatv/powerlevel10k.git /tmp/powerlevel10k && \
    cp /tmp/powerlevel10k/config/p10k.zsh /home/aedes/.p10k.zsh && \
    git clone --depth=1 https://github.com/ohmyzsh/ohmyzsh.git /home/aedes/.oh-my-zsh && \
    git clone --depth=1 https://github.com/dracula/zsh.git /home/aedes/.oh-my-zsh/themes/dracula && \
    git clone --depth=1 https://github.com/zsh-users/zsh-completions.git /home/aedes/.oh-my-zsh/custom/plugins/zsh-completions && \
    git clone --depth=1 https://github.com/zsh-users/zsh-autosuggestions.git /home/aedes/.oh-my-zsh/custom/plugins/zsh-autosuggestions && \
    git clone --depth=1 https://github.com/zsh-users/zsh-syntax-highlighting.git /home/aedes/.oh-my-zsh/custom/plugins/zsh-syntax-highlighting && \
    rm -rf /home/aedes/.oh-my-zsh/.git && \
    rm -rf /home/aedes/.oh-my-zsh/themes/dracula/.git && \
    rm -rf /home/aedes/.oh-my-zsh/custom/plugins/zsh-completions/.git && \
    rm -rf /home/aedes/.oh-my-zsh/custom/plugins/zsh-autosuggestions/.git && \
    rm -rf /home/aedes/.oh-my-zsh/custom/plugins/zsh-syntax-highlighting/.git && \
    cp /home/aedes/.oh-my-zsh/templates/zshrc.zsh-template /home/aedes/.zshrc && \
    sed -i 's/ZSH_THEME=".*"/ZSH_THEME="dracula\/dracula"/' /home/aedes/.zshrc && \
    sed -i 's/plugins=(git)/plugins=(git autojump zsh-completions zsh-autosuggestions zsh-syntax-highlighting)/' /home/aedes/.zshrc && \
    echo 'export DISABLE_AUTO_UPDATE="true"' >> /home/aedes/.zshrc && \
    echo 'export DISABLE_UPDATE_PROMPT="true"' >> /home/aedes/.zshrc && \
    echo 'POWERLEVEL9K_DISABLE_CONFIGURATION_WIZARD=true' >> /home/aedes/.zshrc && \
    sed -i 's/typeset -g POWERLEVEL9K_TIME_BACKGROUND=.*/typeset -g POWERLEVEL9K_TIME_BACKGROUND=magenta/' /home/aedes/.p10k.zsh || echo 'typeset -g POWERLEVEL9K_TIME_BACKGROUND=magenta' >> /home/aedes/.p10k.zsh && \
    # Configure shell aliases
    echo 'if command -v lsd > /dev/null; then' >> /home/aedes/.zshrc && \
    echo '  export LS_COLORS="di=1;34:ln=1;36:so=1;35:pi=1;94:ex=1;31:bd=1;95:cd=1;96:ur=0;32:uw=0;33:ux=0;31:ue=0;32:gr=0;32:gw=0;33:gx=0;31:tr=0;90:tw=0;93:tx=0;92"' >> /home/aedes/.zshrc && \
    echo '  alias ls="lsd --color=always --header"' >> /home/aedes/.zshrc && \
    echo '  alias ll="colorls --long --almost-all --sort-dirs --git-status"' >> /home/aedes/.zshrc && \
    echo '  alias la="lsd -la --color=always --header"' >> /home/aedes/.zshrc && \
    echo '  alias lt="lsd --tree --color=always --header"' >> /home/aedes/.zshrc && \
    echo 'fi' >> /home/aedes/.zshrc && \
    # Color configuration
    echo '# Color environment' >> /home/aedes/.zshrc && \
    echo 'export LSCOLORS="Gxfxcxdxbxegedabagacad"' >> /home/aedes/.zshrc && \
    echo 'export TERM="xterm-256color"' >> /home/aedes/.zshrc && \
    echo 'export COLORTERM="truecolor"' >> /home/aedes/.zshrc && \
    echo 'export EZA_COLORS="ur=0:uw=0:ux=0:ue=0:gr=0:gw=0:gx=0:tr=0:tw=0:tx=0:su=0:sf=0:xa=0"' >> /home/aedes/.zshrc && \
    echo 'export COLORFGBG="15;0"' >> /home/aedes/.zshrc && \
    # Initialize micromamba (not conda)
    echo '# Initialize micromamba' >> /home/aedes/.zshrc && \
    echo 'eval "$(micromamba shell hook --shell zsh)"' >> /home/aedes/.zshrc && \
    echo 'micromamba activate base' >> /home/aedes/.zshrc && \
    # Local adaptation tools
    echo '# Local adaptation tools' >> /home/aedes/.zshrc && \
    echo 'export PATH="/opt/BayesAss3-SNPs:/opt/BA3-SNPS-autotune:/opt/sambada/binaries:$PATH"' >> /home/aedes/.zshrc && \
    # Install fzf
    mkdir -p ~/.fzf && \
    git clone --depth=1 https://github.com/junegunn/fzf.git ~/.fzf && \
    ~/.fzf/install --all && \
    rm -rf ~/.fzf/.git && \
    echo '[ -f ~/.fzf.zsh ] && source ~/.fzf.zsh' >> /home/aedes/.zshrc && \
    echo 'export FZF_BASE=~/.fzf' >> /home/aedes/.zshrc && \
    rm -rf /tmp/powerlevel10k

# Create colorls configuration
RUN mkdir -p ~/.config/colorls && \
    echo "unrecognized_file: white" > ~/.config/colorls/dark_colors.yaml && \
    echo "recognized_file: white" >> ~/.config/colorls/dark_colors.yaml && \
    echo "executable_file: red" >> ~/.config/colorls/dark_colors.yaml && \
    echo "dir: blue" >> ~/.config/colorls/dark_colors.yaml && \
    echo "user: magenta" >> ~/.config/colorls/dark_colors.yaml && \
    echo "group: cyan" >> ~/.config/colorls/dark_colors.yaml && \
    echo "date: yellow" >> ~/.config/colorls/dark_colors.yaml && \
    echo "time: darkgreen" >> ~/.config/colorls/dark_colors.yaml && \
    echo "file_size: palegreen" >> ~/.config/colorls/dark_colors.yaml && \
    echo "read: darkgreen" >> ~/.config/colorls/dark_colors.yaml && \
    echo "write: yellow" >> ~/.config/colorls/dark_colors.yaml && \
    echo "exec: red" >> ~/.config/colorls/dark_colors.yaml && \
    echo "no_access: gray" >> ~/.config/colorls/dark_colors.yaml && \
    echo "image: magenta" >> ~/.config/colorls/dark_colors.yaml && \
    echo "video: blue" >> ~/.config/colorls/dark_colors.yaml && \
    echo "music: cyan" >> ~/.config/colorls/dark_colors.yaml && \
    echo "log: yellow" >> ~/.config/colorls/dark_colors.yaml

# Set working directory and create ALL project directories
WORKDIR /proj

# Create all project subdirectories for local adaptation analysis
RUN mkdir -p /proj/{data/{raw,processed,references,metadata},results/{interim,organized,analysis,populations,local_adaptation,probes,segregation,global_brazil,ba3,admixture},scripts/{download,analysis,visualization},configs,containers,logs,output}

# Expose Jupyter port
EXPOSE 8888

# Default command
CMD ["/bin/zsh"]

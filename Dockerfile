FROM mambaorg/micromamba:1.5.0

ENV MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:/usr/bin:$PATH \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

SHELL ["bash", "-lc"]
USER root

# Update micromamba and all packages to latest versions
RUN micromamba update --all -y

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
    && apt-get clean && rm -rf /var/lib/apt/lists/* \
    && groupadd -r aedes && useradd -r -g aedes aedes \
    && mkdir -p /home/aedes && chown aedes:aedes /home/aedes

# Install lsd manually (not available in conda-forge)
RUN wget https://github.com/lsd-rs/lsd/releases/download/v1.1.5/lsd-v1.1.5-x86_64-unknown-linux-gnu.tar.gz && \
    tar -xzf lsd-v1.1.5-x86_64-unknown-linux-gnu.tar.gz && \
    mv lsd-v1.1.5-x86_64-unknown-linux-gnu/lsd /usr/local/bin/ && \
    rm -rf lsd-v1.1.5-*

# Install core system and Python packages first
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
    -y && micromamba clean --all --yes

# Install bioinformatics Python packages
RUN micromamba install --channel-priority strict -c conda-forge -c bioconda \
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
    -y && micromamba clean --all --yes

# Install geospatial packages
RUN micromamba install --channel-priority strict -c conda-forge -c bioconda \
    geopandas \
    fiona \
    rasterio \
    pyproj \
    shapely \
    folium \
    contextily \
    earthpy \
    geoplot \
    -y && micromamba clean --all --yes

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
    -y && micromamba clean --all --yes

# Install additional analysis tools
RUN micromamba install --channel-priority strict -c conda-forge -c bioconda \
    admixture \
    iqtree \
    starship \
    datamash \
    openjdk=17 \
    sra-tools \
    entrez-direct \
    snakemake \
    -y && micromamba clean --all --yes

# Install Jupyter ecosystem
RUN micromamba install --channel-priority strict -c conda-forge -c bioconda \
    jupyter \
    jupyterlab \
    notebook \
    ipykernel \
    -y && micromamba clean --all --yes

# Install R and R packages
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
    -y && micromamba clean --all --yes

# Install R spatial and visualization packages
RUN micromamba install --channel-priority strict -c conda-forge -c bioconda \
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
    -y && micromamba clean --all --yes

# Install Bioconductor packages
RUN micromamba install --channel-priority strict -c conda-forge -c bioconda \
    bioconductor-variantannotation \
    bioconductor-snprelate \
    bioconductor-annotationdbi \
    bioconductor-biomart \
    bioconductor-biostrings \
    -y && micromamba clean --all --yes

# Install Ruby separately (with its bundled gem)
RUN micromamba install -c conda-forge ruby=3.2.2 -y && micromamba clean --all --yes

# Install R packages for population genetics and local adaptation analysis (with parallel processing)
RUN R -e "install.packages(c('data.table', 'tidyverse', 'qqman', 'qqplotr', 'reticulate', 'broom', \
    'readxl', 'writexl', 'knitr', 'rmarkdown', 'pegas', 'ape', 'phangorn', 'vcfR', 'genetics', \
    'BiocManager', 'remotes', 'scales', 'ggrepel', 'ggtext', 'ggvenn', 'ggstatsplot', 'ggforce', \
    'ggpattern', 'scatterpie', 'RColorBrewer', 'extrafont', 'forcats', 'flextable', 'officer', \
    'Cairo', 'dartR', 'OutFLANK', 'geosphere', 'rnaturalearth', 'rnaturalearthdata', 'ggspatial', \
    'fields', 'grid', 'ellipse', 'reshape2', 'admixr', 'qvalue', 'genomation', 'regioneR'), \
    repos='https://cloud.r-project.org/', dependencies=TRUE, Ncpus=4)"

# Install Bioconductor packages not available via conda
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
    BiocManager::install(c('LEA'), update = FALSE, ask = FALSE)"

# Install GitHub R packages for local adaptation
RUN R -e "options(Ncpus = 4); \
    tryCatch(devtools::install_github('uqrmaie1/admixtools'), error = function(e) print(paste('Failed:', e))); \
    tryCatch(devtools::install_github('gearslaboratory/gdalUtils', force = TRUE), error = function(e) print(paste('Failed:', e))); \
    tryCatch(devtools::install_github('petrelharp/templater'), error = function(e) print(paste('Failed:', e))); \
    tryCatch(devtools::install_github('petrelharp/local_pca/lostruct'), error = function(e) print(paste('Failed:', e))); \
    tryCatch(remotes::install_github('eriqande/TESS3_encho_sen'), error = function(e) print(paste('Failed:', e))); \
    tryCatch(remotes::install_github('eriqande/genoscapeRtools'), error = function(e) print(paste('Failed:', e))); \
    tryCatch(devtools::install_github('SolangeD/R.SamBada', build_vignettes = FALSE), error = function(e) print(paste('Failed:', e))); \
    tryCatch(devtools::install_github('bcm-uga/lfmm'), error = function(e) print(paste('Failed:', e))); \
    tryCatch(devtools::install_github('bcm-uga/TESS3_encho_sen'), error = function(e) print(paste('Failed:', e))); \
    tryCatch(devtools::install_github('bcm-uga/LEA'), error = function(e) print(paste('Failed:', e))); \
    tryCatch(devtools::install_github('petrikemppainen/LDna', ref = 'v.2.15'), error = function(e) print(paste('Failed:', e)))"

# Install Python packages not in conda-forge
RUN pip3 install --no-cache-dir pong

# Install colorls globally as root
RUN /opt/conda/bin/gem install colorls --no-document -n /usr/local/bin && \
    chmod +x /usr/local/bin/colorls

# Download and install local adaptation tools (with error handling)
RUN R -e "library(R.SamBada); downloadSambada('/opt/sambada')" && \
    find /opt/sambada -name "sambada*" -type f -executable -exec chmod +x {} \; && \
    find /opt/sambada -name "sambada*" -type f -executable -exec ln -sf {} /usr/local/bin/ \; && \
    cd /opt && \
    git clone https://github.com/DReichLab/AdmixTools.git && \
    cd AdmixTools && cd src && make && make install && cd .. && \
    ln -sf /opt/AdmixTools/bin/* /usr/local/bin/ && \
    cd /opt && \
    wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip && \
    unzip BayeScan2.1.zip && \
    chmod +x BayeScan2.1/binaries/BayeScan2.1_linux64bits && \
    ln -sf /opt/BayeScan2.1/binaries/BayeScan2.1_linux64bits /usr/local/bin/bayescan && \
    rm BayeScan2.1.zip && \
    wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.4/gemma-0.98.4-linux-static-AMD64.gz && \
    gunzip gemma-0.98.4-linux-static-AMD64.gz && \
    chmod +x gemma-0.98.4-linux-static-AMD64 && \
    mv gemma-0.98.4-linux-static-AMD64 /usr/local/bin/gemma && \
    git clone https://github.com/stevemussmann/BayesAss3-SNPs.git && \
    cd BayesAss3-SNPs && \
    mv BA3-SNPS-Ubuntu64 BA3-SNPS && \
    chmod +x BA3-SNPS && \
    ln -sf /opt/BayesAss3-SNPs/BA3-SNPS /usr/local/bin/BA3-SNPS && \
    cd /opt && \
    git clone https://github.com/stevemussmann/BA3-SNPS-autotune.git && \
    cd BA3-SNPS-autotune && \
    chmod +x BA3-SNPS-autotune.py && \
    ln -sf /opt/BA3-SNPS-autotune/BA3-SNPS-autotune.py /usr/local/bin/BA3-SNPS-autotune || true

# Switch to aedes user
USER aedes

# Install and configure Oh-My-Zsh, Powerlevel10k, and all plugins in a single layer
RUN git clone --depth=1 https://github.com/romkatv/powerlevel10k.git /tmp/powerlevel10k \
  && cp /tmp/powerlevel10k/config/p10k.zsh /home/aedes/.p10k.zsh \
  && git clone https://github.com/ohmyzsh/ohmyzsh.git /home/aedes/.oh-my-zsh \
  && git clone https://github.com/dracula/zsh.git /home/aedes/.oh-my-zsh/themes/dracula \
  && git clone https://github.com/zsh-users/zsh-completions.git /home/aedes/.oh-my-zsh/custom/plugins/zsh-completions \
  && git clone https://github.com/zsh-users/zsh-autosuggestions.git /home/aedes/.oh-my-zsh/custom/plugins/zsh-autosuggestions \
  && git clone https://github.com/zsh-users/zsh-syntax-highlighting.git /home/aedes/.oh-my-zsh/custom/plugins/zsh-syntax-highlighting \
  && cp /home/aedes/.oh-my-zsh/templates/zshrc.zsh-template /home/aedes/.zshrc \
  && sed -i 's/ZSH_THEME=".*"/ZSH_THEME="dracula\/dracula"/' /home/aedes/.zshrc \
  && sed -i 's/plugins=(git)/plugins=(git autojump zsh-completions zsh-autosuggestions zsh-syntax-highlighting)/' /home/aedes/.zshrc \
  && echo 'export DISABLE_AUTO_UPDATE="true"' >> /home/aedes/.zshrc \
  && echo 'export DISABLE_UPDATE_PROMPT="true"' >> /home/aedes/.zshrc \
  && echo 'POWERLEVEL9K_DISABLE_CONFIGURATION_WIZARD=true' >> /home/aedes/.zshrc \
  && sed -i 's/typeset -g POWERLEVEL9K_TIME_BACKGROUND=.*/typeset -g POWERLEVEL9K_TIME_BACKGROUND=magenta/' /home/aedes/.p10k.zsh || echo 'typeset -g POWERLEVEL9K_TIME_BACKGROUND=magenta' >> /home/aedes/.p10k.zsh \
  && echo 'if command -v lsd > /dev/null; then' >> /home/aedes/.zshrc \
  && echo '  sed -i "/alias ls=/d" /home/aedes/.zshrc' >> /home/aedes/.zshrc \
  && echo '  sed -i "/LS_COLORS=/d" /home/aedes/.zshrc' >> /home/aedes/.zshrc \
  && echo '  export LS_COLORS="di=1;34:ln=1;36:so=1;35:pi=1;94:ex=1;31:bd=1;95:cd=1;96:ur=0;32:uw=0;33:ux=0;31:ue=0;32:gr=0;32:gw=0;33:gx=0;31:tr=0;90:tw=0;93:tx=0;92"' >> /home/aedes/.zshrc \
  && echo '  alias ls="lsd --color=always --header"' >> /home/aedes/.zshrc \
  && echo '  alias ll="colorls --long --almost-all --sort-dirs --git-status"' >> /home/aedes/.zshrc \
  && echo '  alias la="lsd -la --color=always --header"' >> /home/aedes/.zshrc \
  && echo '  alias lt="lsd --tree --color=always --header"' >> /home/aedes/.zshrc \
  && echo 'fi' >> /home/aedes/.zshrc \
  && echo '# Match local colorls environment' >> /home/aedes/.zshrc \
  && echo 'export LSCOLORS="Gxfxcxdxbxegedabagacad"' >> /home/aedes/.zshrc \
  && echo 'export TERM="xterm-256color"' >> /home/aedes/.zshrc \
  && echo 'export COLORTERM="truecolor"' >> /home/aedes/.zshrc \
  && echo 'export EZA_COLORS="ur=0:uw=0:ux=0:ue=0:gr=0:gw=0:gx=0:tr=0:tw=0:tx=0:su=0:sf=0:xa=0"' >> /home/aedes/.zshrc \
  && echo 'export COLORFGBG="15;0"' >> /home/aedes/.zshrc \
  && echo '# Initialize conda' >> /home/aedes/.zshrc \
  && echo 'export PATH="/opt/conda/bin:${PATH}"' >> /home/aedes/.zshrc \
  && echo 'source /opt/conda/etc/profile.d/conda.sh' >> /home/aedes/.zshrc \
  && echo 'conda activate base' >> /home/aedes/.zshrc \
  && echo '# Local adaptation tools' >> /home/aedes/.zshrc \
  && echo 'export PATH="/opt/BayesAss3-SNPs:/opt/BA3-SNPS-autotune:/opt/sambada/binaries:$PATH"' >> /home/aedes/.zshrc \
  && mkdir -p ~/.fzf \
  && git clone --depth 1 https://github.com/junegunn/fzf.git ~/.fzf \
  && ~/.fzf/install --all \
  && echo '[ -f ~/.fzf.zsh ] && source ~/.fzf.zsh' >> /home/aedes/.zshrc \
  && echo 'export FZF_BASE=~/.fzf' >> /home/aedes/.zshrc

# Create colorls configuration directory and file for aedes user
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

# Clean up only cache files
RUN rm -rf /var/tmp/* 2>/dev/null || true

# Default command
CMD ["/bin/zsh"]

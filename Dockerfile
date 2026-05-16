# =============================================================================
# Dockerfile — CITE-Seq Analysis Pipeline
# Base: Rocker r-ver (Debian Bookworm, R 4.5.0)
# Reproduces: RNA + ADT multimodal integration, Symphony reference mapping,
#             T-cell progenitor identification and developmental annotation.
#
# Build:  docker build -t citeseq-analysis:1.0 .
# Run:    docker run --rm -v $(pwd)/data:/app/data \
#                         -v $(pwd)/results:/app/results \
#                         citeseq-analysis:1.0
# =============================================================================

FROM rocker/r-ver:4.5.0

LABEL maintainer="Elwin Paulose"
LABEL description="CITE-Seq multimodal analysis pipeline"
LABEL version="1.0"
LABEL org.opencontainers.image.source="https://github.com/elwin-paul/CITE-Seq-Analysis"

# ── 1. System libraries ───────────────────────────────────────────────────────
# Installed once; cached by Docker layer unless this block changes.
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    # Compression & archive
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    # Network & SSL
    libcurl4-openssl-dev \
    libssl-dev \
    libssh2-1-dev \
    # XML / HTML parsing
    libxml2-dev \
    libxslt1-dev \
    # HDF5 (SeuratDisk)
    libhdf5-dev \
    hdf5-tools \
    # Linear algebra (PCA / WNN)
    libblas-dev \
    liblapack-dev \
    libopenblas-dev \
    # Font & graphics rendering (ggplot2 / ggsave PNG)
    libfontconfig1-dev \
    libfreetype6-dev \
    libcairo2-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libxt-dev \
    # Git (for remotes::install_github)
    git \
    # Pandoc (R Markdown reports)
    pandoc \
    # Utilities
    wget \
    curl \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ── 2. Set CRAN mirror and R options ─────────────────────────────────────────
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/bookworm/latest"))' \
    >> /usr/local/lib/R/etc/Rprofile.site

# ── 3. Install renv for reproducible package restoration ─────────────────────
RUN Rscript -e "install.packages('renv', version = '1.0.7')"

# ── 4. Copy renv lockfile and restore all packages ───────────────────────────
# Only renv.lock is copied first so Docker cache is invalidated only when
# the lockfile changes, not when scripts change.
WORKDIR /app
COPY renv.lock renv.lock

RUN Rscript -e " \
    renv::consent(provided = TRUE); \
    renv::restore(lockfile = 'renv.lock', prompt = FALSE) \
    "

# ── 5. Install GitHub-only packages not on CRAN / Bioconductor ───────────────
# Pinned to exact commit SHAs that match renv.lock entries.
RUN Rscript -e " \
    remotes::install_github('mojaveazure/seurat-disk@4c37539c76f75e9fe3e03b0d1ca58ab84fef1da2', \
                             upgrade = 'never', quiet = TRUE); \
    remotes::install_github('immunogenomics/symphony@v0.1.1', \
                             upgrade = 'never', quiet = TRUE); \
    remotes::install_github('immunogenomics/BoneMarrowMap@main', \
                             upgrade = 'never', quiet = TRUE) \
    "

# ── 6. Copy pipeline scripts ──────────────────────────────────────────────────
COPY scripts/ /app/scripts/
COPY Main_Pipeline_Target.R /app/Main_Pipeline_Target.R

# ── 7. Set up expected directory structure ────────────────────────────────────
RUN mkdir -p /app/data \
             /app/results \
             /app/figures \
             /app/outputs \
             /app/docs \
             /app/workflow

# ── 8. Environment variables ──────────────────────────────────────────────────
# Allow large R objects without hitting memory limits
ENV R_MAX_VSIZE=32Gb
# Parallelism: adjust to match your host CPU count
ENV OMP_NUM_THREADS=4
# Suppress R startup messages in logs
ENV R_QUIET=1

# ── 9. Default entrypoint ─────────────────────────────────────────────────────
# Mount your data directory at /app/data before running.
# Results will be written to /app/results and /app/figures.
#
# Override entrypoint to run a single module:
#   docker run ... citeseq-analysis:1.0 Rscript scripts/Module_02_qc_filtered.R
ENTRYPOINT ["Rscript"]
CMD ["Main_Pipeline_Target.R"]

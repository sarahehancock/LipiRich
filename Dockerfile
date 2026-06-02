# ─────────────────────────────────────────────────────────────────────────────
# LipiDash — Dockerfile
# Base: rocker/shiny-verse (R 4.5.2 + Shiny Server + tidyverse pre-installed)
# Bioconductor: 3.22
# ─────────────────────────────────────────────────────────────────────────────

FROM rocker/shiny-verse:4.5.2

# ── System dependencies ───────────────────────────────────────────────────────
# libxml2-dev    → xml2 / curl
# libssl-dev     → openssl / httr
# libcurl4-*     → curl / httr / readr
# libharfbuzz-*
# libfribidi-dev → textshaping / ragg
# libfreetype6-* → ragg / systemfonts
# libpng-dev     → png graphics
# libtiff5-dev   → tiff support in ragg
# libjpeg-dev    → jpeg in ragg
# libfontconfig  → system font discovery
# libgit2-dev    → gert (optional but common dep)
# libsodium-dev  → sodium (Shiny auth helpers)
# libgmp-dev     → arbitrary precision (polynom)
# pandoc         → rmarkdown/knitr
# libglpk-dev    → igraph (optional, pulled by visNetwork)
# libabsl-dev    → s2 / sf dependency chain
# libudunits2-dev→ units / sf dependency chain
# cmake          → nloptr / lme4 / car dependency chain
# cargo / rustc  → some Bioc packages with Rust extensions

RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    libgit2-dev \
    libsodium-dev \
    libgmp-dev \
    libglpk-dev \
    libabsl-dev \
    libudunits2-dev \
    pandoc \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# ── R package installation ────────────────────────────────────────────────────
# ARG CACHE_BUST forces Docker to re-run the package install step when changed.
# Update the value here whenever you want a clean reinstall (e.g. new packages.R)
ARG CACHE_BUST=2026-05-17b
COPY packages.R /tmp/packages.R
RUN Rscript /tmp/packages.R

# ── Shiny Server configuration ────────────────────────────────────────────────
# Serves the app at / so the Cloudflare tunnel URL http://lipidash:3838
# goes straight to the app without needing a subpath.
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# ── App files ─────────────────────────────────────────────────────────────────
# Remove the default Shiny Server sample apps and deploy LipiDash
RUN rm -rf /srv/shiny-server/*
COPY app.R /srv/shiny-server/lipidash/app.R

# ── Persistent data directory ─────────────────────────────────────────────────
# Mount a host volume here in Portainer to persist any cache files across
# container restarts and image rebuilds.
RUN mkdir -p /srv/shiny-data && chmod 777 /srv/shiny-data

# ── Runtime settings ──────────────────────────────────────────────────────────
USER shiny

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]

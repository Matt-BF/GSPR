FROM rocker/tidyverse:latest

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libglpk-dev libxt-dev cmake build-essential libcurl4-openssl-dev libnode-dev libnlopt-dev

RUN R -q -e 'install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))'

# Clear R package cache to prevent possible corrupted downloads
RUN R -q -e 'pak::cache_clean()'

RUN R -q -e 'pak::pkg_install(c("brunnermunzel", "ggforce", "patchwork", "gt", "khroma", "Biostrings", "viridis", "gtools", "dunn.test", "fs", "ggokabeito", "ggtree", "ggtreeExtra", "ggnewscale", "Peptides", "lmerTest", "gluedown", "TreeDist", "Quartet", "cowplot", "ape", "phytools", "ggblend", "ggseqlogo", "bio3d", "colorspace", "plotwidgets", "PCAtools", "cluster", "dbscan", "sf", "rnaturalearth", "rnaturalearthdata", "fullRankMatrix", "caret", "iNEXT", "styler"))'
RUN R -q -e 'install.packages("geosphere")'

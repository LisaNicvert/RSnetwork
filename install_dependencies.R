# Header #############################################################
#
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
#
# Date: 2024-10-03
#
# Script Description: this scripts installs all the necessary dependencies to
# run this project.

# Dependencies for the custom package
packages <- c("ade4",
              "adegraphics",
              "dplyr",
              "ggplot2",
              "ggrepel",
              "grid",
              "lmtest",
              "mvtnorm",
              "stats",
              "stringr",
              "tibble",
              "tidyr",
              "tidyselect")
base::lapply(packages, require)

# Install custom package
devtools::install_local(upgrade = "never")

# Other required libraries (for the analyses)
packages <- c("here", "patchwork", "RColorBrewer", "viridisLite",
              "tidytext", "effectsize", "DT",
              "igraph", "tidygraph", "ggraph")
base::lapply(packages, require)

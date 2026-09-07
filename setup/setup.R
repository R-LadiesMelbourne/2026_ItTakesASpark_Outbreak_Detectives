
# ============================================================
# Outbreak Detectives - workshop setup
# ============================================================

needed_packages <- c(
  "tidyverse",
  "scales",
  "ape",
  "phangorn",
  "igraph",
  "ggraph"
)

missing_packages <- needed_packages[
  !vapply(
    needed_packages,
    requireNamespace,
    FUN.VALUE = logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste(
      "These packages need to be installed before the workshop:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(tidyverse)

source("setup/workshop_functions.R")

## setwd("~/research/SoftWare/forestCausalSearch//")
library(targets)
# ---------------------------------------------------------------------
# packages
# ---------------------------------------------------------------------
thepackages <- c("targets",
                 "future.batchtools",
                 "tarchetypes",
                 "future.callr",
                 "future",
                 "grf",
                 "ranger",
                 ## "nleqslv"
                 "data.table",
                 "scales",
                 "riskRegression",
                 "prodlim",
                 "survival",
                 "foreach",
                 "parallel",
                 "grid",
                 "ggplot2",
                 "gridExtra",
                 "ggplotify",
                 "cowplot")
library(targets)
targets::tar_option_set(packages = thepackages)
# ---------------------------------------------------------------------
# R functions
# ---------------------------------------------------------------------
for(f in list.files("R",".R$",full.names=TRUE)){source(f)}
for(f in list.files("functions",".R$",full.names=TRUE)){source(f)}
# ---------------------------------------------------------------------
# Simulation settings
# ---------------------------------------------------------------------
source("./setting/simulation_targets.R")

# ---------------------------------------------------------------------
# The target flow
# ---------------------------------------------------------------------

## MCCORES <- 5
## MCCORES <- 50
## MCCORES are set in setting/simulation-targets.R
list(varying_target,
     fixed_target,
     fixed,
     truth_varying,
     truth,
     estimates,
     ate,
     results,
     ranking,
     plotframe)





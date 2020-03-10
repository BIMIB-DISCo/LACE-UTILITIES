# set the working directory
setwd("/data_directory/BimiB/share/LACE/MELANOMA/real_data")

# load required scripts
library("timescape")
source("R/utils.R")
source("R/visualization.R")

# inference
load("results/inference.RData")
draw.B.ts(inference$B,inference$C)

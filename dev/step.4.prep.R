
list.of.packages <- c("pbapply","EpiModel", "Rcpp","ergm","btergm","texreg","rstudioapi","data.table", "haven", "networkDynamic", "intergraph", "ape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)

## load library
library("EpiModel")
library("sna")
library('magrittr')
require(data.table)
source("dev/helper-functions.R")

## re-load all required files
## initial status vector
source_lines("~/Dropbox/GitHub/NetSim/EpiModel.R", 39)
source_lines("~/Dropbox/GitHub/NetSim/EpiModel.R", 208:211)

## load network simulations
source_lines("~/Dropbox/GitHub/NetSim/EpiModel.R", 274:275)

## define custom modules
source_lines("~/Dropbox/GitHub/NetSim/EpiModel.R", 284:526)

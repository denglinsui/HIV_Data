#' This file sets the enviroment

rm(list = ls())
# Output Table
library(knitr)

# Plot
library(ggplot2)
library(ggpubr)

# Tidy Dataset
library(dplyr)
library(reshape)

# To discover gene
library(knockoff)

setwd("~/project/HIV_Data")
# Extract codes of methods
source("method/multilayer_knockoff_filter.R")
source("method/aux_mkf.R")
source("method/aux_knockoffs.R")
source("method/p_filter.R")

# Extract and clean data
source("codes/DataExtract.R")
source("codes/RunProcedure.R")
source("codes/ResultExtract.R")
source("codes/PlotData.R")


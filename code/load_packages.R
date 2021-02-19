# ============================================= #
# script: load_packages.R
# Project: POOLS Projects
# Author(s): R.N. Padgett, L. Shero, S. Jiang,
#              & T. Kettler
# ============================================= #
# Data Created: 2020-10-14
# Date Modified: 2020-10-14
# By: R. Noah Padgett                   
# ============================================= #
# Stems from work in EDP 6367                   
# ============================================= #
# Purpose:
# This R script is for loading all necessary
#   R packages
#
# No output - just loading packages into the 
#   environment
# ============================================= #
# Set up directory and libraries
rm(list=ls())
# list of packages
packages <- c("tidyverse", "readr", "patchwork",
              "tidyr","data.table", "dplyr","ggplot2",
              "lavaan", "semTools", "lavaanPlot",
              "MIIVsem", "lslx",
              "simsem", "naniar", "ggcorrplot",
              "mvtnorm", "psychometric", "psych",
              "nFactors", "coda",
              "readxl",
              "kableExtra", "xtable")   
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])] 
if(length(new.packages)) install.packages(new.packages) 
# Load packages
lapply(packages, library, character.only = TRUE)

w.d <- getwd()


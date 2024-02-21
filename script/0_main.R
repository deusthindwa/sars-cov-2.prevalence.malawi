#written by Deus
#13/02/2024
#pneumococcal carriage and serotype dynamics by adult HIV status a mature PCV program

#====================================================================

#load packages
pacman::p_load(char = c("lubridate", "tidyverse", "dplyr", "here", "rio", "scales", "boot", "magrittr",  "mvtnorm", "zoo", "patchwork", "ggplotify", "sf",
                        "PropCIs", "reshape2","purrr", "msm", "minqa", "ggridges", "timetk", "ggbreak", "ggpubr", "gridExtra", "doParallel", "igraph", "rgdal"))

#====================================================================

#set seed for entire session globally to ensure reproducibility using a task call
addTaskCallback(function(...) {set.seed(12345);TRUE})

#turn off the global task call for set seed if needed
#removeTaskCallback(1)

#getting and putting datasets in right order for analysis
source(here("script", "1_data_wrangling.R"))

#characterizing spn carriage at baseline and follow up
source(here("script", "2_carriage_char.R"))

#run SIS Markov model
source(here("script", "3_markov_modelfit.R"))

#compute carriage acquisition and clearance rates at vaccine-serotype group level
source(here("script", "4_carriage_acq_dur.R"))

#compute hazard ratios comparing different covariate levels
source(here("script", "5_carriage_hazard_ratios.R"))

#compute carriage acquisition and clearance rates at serotype level
source(here("script", "6_serotype_dynamics.R"))

#compute the number of carriage episodes 
source(here("script", "7_carriage_miscelleneous.R"))

#characterize and compute acquisition and clearance for serotypes in each pcv 
source(here("script", "8_pcv_serotype_dynamics.R"))

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

#characterizing spn carriage at baseline
source(here("script", "2_carriage_char.R"))

#run SIS Markov model in a Markov modelling framework
source(here("script", "3_markov_modelfit.R"))

#use SIS Markov model to estimate acquisition and clearance of carriage 
source(here("script", "4_carriage_acq_dur.R"))

#use SIS Markov model to estimate hazard ratios comparing covariate levels
source(here("script", "5_carriage_hazard_ratios.R"))

#use SIS Markov model to estimate the number of carriage episodes 
source(here("script", "6_serotype_dynamics.R"))

#rerun multiple chains SIS Markov model for convergence and model fit check
#source(here("script", "7_model_convergence.R"))

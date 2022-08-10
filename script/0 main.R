#Written by Deus Thindwa
#Estimating SARS-CoV-2 prevalence in Malawi
#Generalized additive model.
#01/8/2022 - 31/12/2022

#=======================================================================

#load the require packages
if (!require(pacman)){install.packages("pacman")}
pacman::p_load(char = c("tidyverse", "dplyr", "lubridate", "patchwork", "rio", "boot", "devtools", "Metrics", "PropCIs", "forecast", "here"))

#=======================================================================

#load PHIM dataset
phimage 
phim <- import(here("data", "phimserosurvey.csv")) %>%
  
#=======================================================================

#data wrangling
  #select(sid:comorbidities) %>%
  mutate(sid = sid,
         
         agegp = as_factor(if_else(agegroup == "<12 yrs", "0-12y",
                                   if_else(agegroup == "13 to 17 yrs", "13-17y",
                                           if_else(agegroup == "18 to 50 yrs", "18-50y",
                                                   if_else(agegroup == "50 yrs", "50+y", NA_character_))))),
         
         sexgp = as_factor(sexgrp),
         
         loc = as_factor(location),
         
         site = as_factor(if_else(site == "BT City", "Blantyre city",
                                  if_else(agegroup == "BT Rural", "Blantyre rural",
                                          if_else(agegroup == "LL City", "Lilongwe city",
                                                  if_else(agegroup == "LL Rural", "Lilongwe rural",
                                                          if_else(agegroup == "ZA City", "Zomba city",
                                                                  if_else(agegroup == "ZA Rural", "Zomba rural", site))))))),
         
         vaxdose = as_factor(Doses),
         
         symp2wks = as_factor(if_else(Symptoms_2wks == "1. Yes", "Yes",
                                      if_else(Symptoms_2wks == "2. No", "No", NA_character_))),
         
         sympany = as_factor(if_else(Symptoms_any == "1. Yes", "Yes",
                                     if_else(Symptoms_any == "2. No", "No", NA_character_))),
         
         workpl = as_factor(if_else(workplace == "1. Outdoors", "Outdoor",
                                    if_else(workplace == "2. Indoors", "Indoor", NA_character_))),
         
         edu = as_factor(if_else(schoollevel == "1. Primary", "Primary",
                                 if_else(schoollevel == "2. Secondary", "Secondary",
                                         if_else(schoollevel == "3. Tertiary", "Tertiary",
                                                 if_else(schoollevel == "None", "None", NA_character_))))),
         
         occ = as_factor(if_else(occupation == "1. Unwaged", "Unwaged",
                                 if_else(occupation == "2. Irregular", "Irregular",
                                         if_else(occupation == "3. Regular", "Regular", NA_character_)))),
         
         ctrl = as_factor(if_else(publchlthmeas == "1. Facemask", "Facemask",
                                  if_else(publchlthmeas == "2. Hand hygiene", "Hygiene",
                                          if_else(publchlthmeas == "3. Social distancing", "Distance",
                                                  if_else(publchlthmeas == "None", "None", NA_character_))))),
         
         serostat = as_factor(if_else(Serostatus == "Positive", "Pos",
                                      if_else(Serostatus == "Negative", "Neg", NA_character_))),
         
         hivstat = as_factor(if_else(hivstatus == 0, "Neg",
                                     if_else(Serostatus == 1, "Pos-ART",
                                             if_else(Serostatus == 2, "Pos+ART", NA_character_)))),
         
         comob = as_factor(if_else(comorbidities == "Yes", "Yes",
                                   if_else(comorbidities == "No", "No", NA_character_)))) %>%
  
#=======================================================================
#select final variables
select(sid, serostat, agegp, sexgp, loc, site, vaxdose, symp2wks, sympany, workpl, edu, occ, ctrl, hivstat, comob) %>%
  left_join()

#=======================================================================

#cleaning/recoding variables for descriptive analysis
pcvpa.des <- pcvpa

#survey number
pcvpa.des$surv <- as.integer(pcvpa.des$surv)

#serogroup
pcvpa.des$serogroup <- if_else(pcvpa.des$serotype == "NoCarriage" | pcvpa.des$serotype == "Pending", "None", 
                               if_else(pcvpa.des$serotype == "NVT", "NVT", "VT"))

pcvpa.des <- select(pcvpa.des, pid, labid, date, surv, serotype, serogroup, age, sex, artdate, artreg, ctx, cd4cnt, nochild5, sescat)

#serotype
pcvpa.des$serotype <- if_else(pcvpa.des$serotype == "NoCarriage" | pcvpa.des$serotype == "Pending", NA_character_, pcvpa.des$serotype)

#age
pcvpa.des$age <- as.integer(pcvpa.des$age)

#sex
pcvpa.des$sex <- if_else(pcvpa.des$sex == 0, "Female", 
                         if_else(pcvpa.des$sex == 1, "Male", 
                                 if_else(pcvpa.des$sex == 2, "Female", NA_character_)))

#ART duration
pcvpa.des <- rename(pcvpa.des %>% 
                      mutate(artdate = as.integer((date-artdate)/365.25)), artdur = artdate)

#ART regimen
pcvpa.des$artreg <- if_else(pcvpa.des$artreg == "TDF/3TC/EFV" | pcvpa.des$artreg == "TDF/3TC/NVP" | pcvpa.des$artreg == "AZT/3TC/EFV", "First line",
                            if_else(pcvpa.des$artreg == "TDF/3TC/LPVr", "Second line", NA_character_))

#Cotrimoxozole
pcvpa.des$ctx <- if_else(pcvpa.des$ctx == 0, "No", 
                         if_else(pcvpa.des$ctx == 1,"Yes", NA_character_))

#CD4+ count
pcvpa.des$cd4cnt <- as.integer(if_else(pcvpa.des$cd4cnt > 1500, NA_integer_, pcvpa.des$cd4cnt))

#living with <5 years-old children
pcvpa.des$nochild5 <- as.integer(pcvpa.des$nochild5)

#social economic status
pcvpa.des$sescat <- if_else(pcvpa.des$sescat == 1, "Low", 
                            if_else(pcvpa.des$sescat == 2, "Middle", 
                                    if_else(pcvpa.des$sescat == 3, "High", NA_character_)))


#=======================================================================

#cleaning/recoding variables for modeling
pcvpa.mod <- pcvpa.des

#overall carriage
pcvpa.mod$carr <- if_else(pcvpa.mod$serogroup == "VT" | pcvpa.mod$serogroup == "NVT", 1L, 
                          if_else(pcvpa.mod$serogroup == "None", 0L, NA_integer_ ))

#vaccine types
pcvpa.mod$vtcarr <- if_else(pcvpa.mod$serogroup == "VT", 1L, 
                            if_else(pcvpa.mod$serogroup == "None" | pcvpa.mod$serogroup == "NVT", 0L, NA_integer_ ))

#non-vaccine types
pcvpa.mod$nvtcarr <- if_else(pcvpa.mod$serogroup == "NVT", 1L, 
                             if_else(pcvpa.mod$serogroup == "None" | pcvpa.mod$serogroup == "VT", 0L, NA_integer_))

#sensitivity of ST3
pcvpa.mod$vtcarr1 <- if_else(pcvpa.mod$serotype != "NVT" & !is.na(pcvpa.mod$serotype) & pcvpa.mod$serotype != 3, 1L, 0L)

pcvpa.mod$nvtcarr1 <- if_else(pcvpa.mod$serotype == "NVT" | pcvpa.mod$serotype == "3", 1L, 0L)
pcvpa.mod$nvtcarr1[is.na(pcvpa.mod$nvtcarr1)] <- 0L

#seasonal variables
pcvpa.mod$seas <- as.integer(month(ymd(pcvpa.mod$date)))
pcvpa.mod$year <- as.integer(year(ymd(pcvpa.mod$date)))

#sex
pcvpa.mod$sex <- if_else(pcvpa.mod$sex == "Female", 1L, 
                         if_else(pcvpa.mod$sex == "Male", 2L, NA_integer_))

#ART duration
pcvpa.mod$artdur <- if_else(pcvpa.mod$artdur <3, 1L, 
                            if_else(pcvpa.mod$artdur >=3 & pcvpa.mod$artdur <20, 2L, NA_integer_))

pcvpa.mod$nochild5 <- if_else(pcvpa.mod$nochild5 == 0, 1L, 
                              if_else(pcvpa.mod$nochild5 >=1 & pcvpa.mod$nochild5 <5, 2L, NA_integer_))

#social economic status
pcvpa.mod$sescat <- if_else(pcvpa.mod$sescat == "Low", 1L, 
                            if_else(pcvpa.mod$sescat == "Middle", 2L, 
                                    if_else(pcvpa.mod$sescat == "High", 2L, NA_integer_)))

#final dataset
pcvpa.mod <- select(pcvpa.mod, pid, labid, carr, nvtcarr, vtcarr, nvtcarr1, vtcarr1, year, age, seas, sex, sescat, artdur, nochild5)

#=======================================================================

#multiple imputation on 
set.seed(1988)
pcvpa.mod1 <- pcvpa.mod %>% 
  select(carr, vtcarr, vtcarr1, nvtcarr, nvtcarr1, year, age, seas, sex, artdur, nochild5, sescat) %>% 
  mutate(artdur = as_factor(artdur),
         nochild5 = as_factor(nochild5),
         sescat = as_factor(sescat))

pcvpa.mod2 <- missForest(pcvpa.mod1)
pcvpa.mod2 <- pcvpa.mod2$ximp
pcvpa.mod2 <- cbind(select(pcvpa.mod, pid, labid), pcvpa.mod2)

pcvpa.mod <- pcvpa.mod2 %>% 
  select(pid, labid, carr, vtcarr, vtcarr1, nvtcarr, nvtcarr1, year, age, seas, sex, artdur, nochild5, sescat) %>% 
  mutate(carr = as.integer(carr),
         vtcarr = as.integer(vtcarr),
         vtcarr1 = as.integer(vtcarr1),
         nvtcarr = as.integer(nvtcarr),
         nvtcarr1 = as.integer(nvtcarr1),
         year = as.integer(year),
         age = as.integer(age),
         seas = as.integer(seas),
         sex = as.integer(sex),
         artdur = as.integer(artdur),
         nochild5 = as.integer(nochild5),
         sescat = as.integer(sescat)) 

rm(pcvpa.mod1, pcvpa.mod2)

#=======================================================================

#descriptive of study population (figure 1)
source(here("script/Fig1_study_descriptive.R"))

#overall and VT carriage dynamics (figure 2)
source(here("script/Fig2_Overall_VT_prev_crude.R"))

#Risk factor VT carriage dynamics (figure 3)
source(here("script/Fig3a_Overall_VT_prev_risk_factors.R"))

#Risk factor VT carriage dynamics (figure 3)
source(here("script/Fig3b_Overall_VT_prev_risk_factors.R"))

#carriage heterogeneity (figure S1)
source(here("script/FigS1_age_time_heterogeneity.R"))

#multiple carriage (figure S2)
source(here("script/FigS2_multiple_carriage.R"))

#model autocorrelation (figure S3)
source(here("scriptFigS3_overall_carriage_ACF/.R"))

#model selection (figure S4)
source(here("script/FigS4_VT_carriage_ACF.R"))

#model selection (figure S5)
source(here("script/FigS5_model_selection.R"))

#model selection (figure S6)
source(here("script/FigS6_age_ART_duration.R"))




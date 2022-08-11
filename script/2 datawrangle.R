#Written by Deus Thindwa
#Estimating SARS-CoV-2 prevalence in Malawi
#Generalized additive model.
#01/8/2022 - 31/12/2022

#=======================================================================

#load PHIM dataset
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


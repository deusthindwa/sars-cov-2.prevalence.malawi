#written by Deus
#01/08/2022
#pneumococcal carriage ad serotype dynamics in HIV-infected adults in the infant PCV era

#====================================================================

#import spn patient level baseline and follow up data
spn_overall <- 
  import(here("data", "pneumo.xlsx")) %>%
  rowwise() %>%
  mutate(ses = sum(watch1, radio1, bank1, ironic1, sew1, mattress1, bed1, bike1, moto1, car1, mobil1, cd1, fanelec1, netyn1, tv1, fridge1), na.rm = TRUE)

#baseline data manipulation
spn_baseline <- 
  spn_overall %>%
  select(ParticipantID, visit_day0, visit_date0, ART_start_date0, HIV_status0, serotype0, sex0, no_und0, density1, cd4_count1, ses, age0, Takenabx1) %>% 
  mutate(visit = date(visit_date0), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit - ART_start_date0)/365.25), 
         carr = 1L) %>%
  ungroup() %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density1", 
         "cd4"= "cd4_count1", 
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype0",
         "age" = "age0",
         "date" = "visit_date0",
         "abx" = "Takenabx1") %>%
  select(pid, dens, sex, age, nochild, ses, hiv, artdur, cd4, serotype, abx)

#group baseline serotypes into VT and NVT
spn_baseline <- 
  spn_baseline %>% 
  ungroup() %>%
  mutate(serotype = if_else(serotype == 99, NA_character_, 
                            if_else(!is.na(serotype), serotype, NA_character_)),
         serogroup = if_else(serotype %in% c("1", "3", "4", "5", "6A", "6B", "7F ", "9V", "9V  n 3", "14", "18C", "19A", "19F", "23F") == TRUE, "VT", 
                             if_else(!is.na(serotype), "NVT", "None")))

#====================================================================

#baseline demographics
#pneumob %>%
#  mutate(nochild = if_else(nochild == 1, "1 child", "2+ children"),
#         agegp = if_else(age >= 18 & age <= 25, "18-25", 
#                if_else(age > 25 & age <= 35, "26-35", 
#                        if_else(age > 35, "36-45", NA_character_))))%>%
#  select(dens, sex, agegp, nochild, ses, hiv, cd4, serogroup) %>%
#  rename("density (CFU/ml)" = "dens", 
#         "number of children" = "nochild",
#         "social economic status" = "ses",
#         "cd4 count (cells/µl)" = "cd4",
#         "age group" = "agegp") %>%
#  tbl_summary(by = hiv, missing = "no") %>% 
#  add_p() %>% 
#  add_overall() %>% 
#  #add_ci(pattern = "{stat} ({ci})") %>% 
#  add_n() %>% 
#  bold_labels()

#Regression Models for density

#====================================================================

#follow up sampling data manipulation
v0 <- 
  spn_overall %>%
  select(ParticipantID, visit_day0, visit_date0, ART_start_date0, HIV_status0, serotype0, sex0, no_und0, ses, age0, lab_ID0) %>% 
  mutate(density0 = NA,
         Takenabx0 = NA,
         cd4_count0 = NA, 
         visit = date(visit_date0), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit - ART_start_date0)/365.25), 
         carr = 1L, 
         day = 0L,
         vday = 0L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density0", 
         "cd4"= "cd4_count0", 
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype0",
         "age" = "age0",
         "date" = "visit",
         "abx" = "Takenabx0",
         "Lab" = "lab_ID0") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)

#====================================================================

v1 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date1, ART_start_date0, HIV_status0, serotype1, density1, cd4_count1, sex0, no_und0, ses, age0, lab_ID1, Takenabx1) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date1 = date(visit_date1), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date1 - ART_start_date0)/365.25), 
         carr = if_else(serotype1 == "99", 0L, 1L, NA_integer_), 
         day = as.integer(visit_date1 - visit_date0),
         vday = 1L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density1", 
         "cd4"= "cd4_count1", 
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype1",
         "age" = "age0",
         "date" = "visit_date1",
         "abx" = "Takenabx1",
         "Lab" = "lab_ID1") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx) %>%
  filter(!is.na(carr))

#====================================================================

v2 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date2, ART_start_date0, HIV_status0, serotype2, density2, sex0, no_und0, ses, age0, lab_ID2, Takenabx2) %>% 
  mutate(cd4 = NA,
         visit_date0 = date(visit_date0), 
         visit_date2 = date(visit_date2), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date2 - ART_start_date0)/365.25), 
         carr = if_else(serotype2 == "99", 0L, 1L, NA_integer_), 
         day = as.integer(visit_date2 - visit_date0),
         vday = 2L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density2", 
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype2",
         "age" = "age0",
         "date" = "visit_date2",
         "abx" = "Takenabx2",
         "Lab" = "lab_ID2") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v3 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date3, ART_start_date0, HIV_status0, serotype3, density3, sex0, no_und0, ses, age0, lab_ID3, Takenabx3) %>% 
  mutate(cd4 = NA,
         visit_date0 = date(visit_date0), 
         visit_date3 = date(visit_date3), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date3 - ART_start_date0)/365.25), 
         carr = if_else(serotype3 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date3 - visit_date0),
         vday = 3L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density3", 
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype3",
         "age" = "age0",
         "date" = "visit_date3",
         "abx" = "Takenabx3",
         "Lab" = "lab_ID3") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v4 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date4, ART_start_date0, HIV_status0, serotype4, density4, sex0, no_und0, ses, age0, lab_ID4, Takenabx4) %>% 
  mutate(cd4 = NA,
         visit_date0 = date(visit_date0), 
         visit_date4 = date(visit_date4), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date4 - ART_start_date0)/365.25), 
         carr = if_else(serotype4 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date4 - visit_date0),
         vday = 4L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density4", 
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype4",
         "age" = "age0",
         "date" = "visit_date4",
         "abx" = "Takenabx4",
         "Lab" = "lab_ID4") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v5 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date5, ART_start_date0, HIV_status0, serotype5, density5, cd4_count5, sex0, no_und0, ses, age0, lab_ID5, Takenabx5) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date5 = date(visit_date5), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date5 - ART_start_date0)/365.25), 
         carr = if_else(serotype5 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date5 - visit_date0),
         vday = 5L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density5",
         "cd4" = "cd4_count5",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype5",
         "age" = "age0",
         "date" = "visit_date5",
         "abx" = "Takenabx5",
         "Lab" = "lab_ID5") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v6 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date6, ART_start_date0, HIV_status0, serotype6, density6, cd4_count6, sex0, no_und0, ses, age0, lab_ID6, Takenabx6) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date6 = date(visit_date6), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date6 - ART_start_date0)/365.25), 
         carr = if_else(serotype6 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date6 - visit_date0),
         vday = 6L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density6",
         "cd4" = "cd4_count6",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype6",
         "age" = "age0",
         "date" = "visit_date6",
         "abx" = "Takenabx6",
         "Lab" = "lab_ID6") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v7 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date7, ART_start_date0, HIV_status0, serotype7, density7, cd4_count7, sex0, no_und0, ses, age0, lab_ID7, Takenabx7) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date7 = date(visit_date7), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date7 - ART_start_date0)/365.25), 
         carr = if_else(serotype7 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date7 - visit_date0),
         vday = 7L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density7",
         "cd4" = "cd4_count7",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype7",
         "age" = "age0",
         "date" = "visit_date7",
         "abx" = "Takenabx7",
         "Lab" = "lab_ID7") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v8 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date8, ART_start_date0, HIV_status0, serotype8, density8, cd4_count8, sex0, no_und0, ses, age0, lab_ID8, Takenabx8) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date8 = date(visit_date8), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date8 - ART_start_date0)/365.25), 
         carr = if_else(serotype8 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date8 - visit_date0),
         vday = 8L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density8",
         "cd4" = "cd4_count8",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype8",
         "age" = "age0",
         "date" = "visit_date8",
         "abx" = "Takenabx8",
         "Lab" = "lab_ID8") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v9 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date9, ART_start_date0, HIV_status0, serotype9, density9, cd4_count9, sex0, no_und0, ses, age0, lab_ID9, Takenabx9) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date9 = date(visit_date9), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date9 - ART_start_date0)/365.25), 
         carr = if_else(serotype9 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date9 - visit_date0),
         vday = 9L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density9",
         "cd4" = "cd4_count9",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype9",
         "age" = "age0",
         "date" = "visit_date9",
         "abx" = "Takenabx9",
         "Lab" = "lab_ID9") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v10 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date10, ART_start_date0, HIV_status0, serotype10, density10, cd4_count10, sex0, no_und0, ses, age0, lab_ID10, Takenabx10) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date10 = date(visit_date10), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date10 - ART_start_date0)/365.25), 
         carr = if_else(serotype10 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date10 - visit_date0),
         vday = 10L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density10",
         "cd4" = "cd4_count10",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype10",
         "age" = "age0",
         "date" = "visit_date10",
         "abx" = "Takenabx10",
         "Lab" = "lab_ID10") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v11 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date11, ART_start_date0, HIV_status0, serotype11, density11, cd4_count11, sex0, no_und0, ses, age0, lab_ID11, Takenabx11) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date11 = date(visit_date11), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date11 - ART_start_date0)/365.25), 
         carr = if_else(serotype11 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date11 - visit_date0),
         vday = 11L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density11",
         "cd4" = "cd4_count11",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype11",
         "age" = "age0",
         "date" = "visit_date11",
         "abx" = "Takenabx11",
         "Lab" = "lab_ID11") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v12 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date12, ART_start_date0, HIV_status0, serotype12, density12, cd4_count12, sex0, no_und0, ses, age0, lab_ID12, Takenabx12) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date12 = date(visit_date12), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date12 - ART_start_date0)/365.25), 
         carr = if_else(serotype12 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date12 - visit_date0),
         vday = 12L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density12",
         "cd4" = "cd4_count12",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype12",
         "age" = "age0",
         "date" = "visit_date12",
         "abx" = "Takenabx12",
         "Lab" = "lab_ID12") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v13 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date13, ART_start_date0, HIV_status0, serotype13, density13, cd4_count13, sex0, no_und0, ses, age0, lab_ID13, Takenabx13) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date13 = date(visit_date13), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date13 - ART_start_date0)/365.25), 
         carr = if_else(serotype13 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date13 - visit_date0),
         vday = 13L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density13",
         "cd4" = "cd4_count13",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype13",
         "age" = "age0",
         "date" = "visit_date13",
         "abx" = "Takenabx13",
         "Lab" = "lab_ID13") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v14 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date14, ART_start_date0, HIV_status0, serotype14, density14, cd4_count14, sex0, no_und0, ses, age0, lab_ID14, Takenabx14) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date14 = date(visit_date14), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date14 - ART_start_date0)/365.25), 
         carr = if_else(serotype14 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date14 - visit_date0),
         vday = 14L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density14",
         "cd4" = "cd4_count14",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype14",
         "age" = "age0",
         "date" = "visit_date14",
         "abx" = "Takenabx14",
         "Lab" = "lab_ID14") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v15 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date15, ART_start_date0, HIV_status0, serotype15, density15, cd4_count15, sex0, no_und0, ses, age0, lab_ID15, Takenabx15) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date15 = date(visit_date15), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date15 - ART_start_date0)/365.25), 
         carr = if_else(serotype15 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date15 - visit_date0),
         vday = 15L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density15",
         "cd4" = "cd4_count15",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype15",
         "age" = "age0",
         "date" = "visit_date15",
         "abx" = "Takenabx15",
         "Lab" = "lab_ID15") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

v16 <- 
  spn_overall %>%
  select(ParticipantID, visit_date0, visit_date16, ART_start_date0, HIV_status0, serotype16, density16, cd4_count16, sex0, no_und0, ses, age0, lab_ID16, Takenabx16) %>% 
  mutate(visit_date0 = date(visit_date0), 
         visit_date16 = date(visit_date16), 
         ART_start_date0 = date(ART_start_date0), 
         artdur = as.numeric((visit_date16 - ART_start_date0)/365.25), 
         carr = if_else(serotype16 == "99", 0L,  1L, NA_integer_), 
         day = as.integer(visit_date16 - visit_date0),
         vday = 16L) %>%
  rename("pid" = "ParticipantID", 
         "dens" = "density16",
         "cd4" = "cd4_count16",
         "hiv"= "HIV_status0",
         "nochild" = "no_und0",
         "sex" = "sex0",
         "serotype" = "serotype16",
         "age" = "age0",
         "date" = "visit_date16",
         "abx" = "Takenabx16",
         "Lab" = "lab_ID16") %>%
  select(pid, vday, day, date, carr, dens, sex, age,  nochild, ses, hiv, artdur, cd4, serotype, abx)  %>%
  filter(!is.na(carr))

#====================================================================

#combine all spn sampling visits
spn_fup <- rbind(v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16)

#group visit serotypes into VT and NVT
spn_fup <- 
  spn_fup %>% 
  ungroup() %>%
  mutate(serotype = if_else(serotype == 99, NA_character_, 
                          if_else(!is.na(serotype), serotype, NA_character_)),
         serogroup = if_else(serotype %in% c("1", "3", "4", "5", "6A", "6B", "7F ", "9V", "9V  n 3", "14", "18C", "19A", "19F", "23F") == TRUE, "VT", 
                             if_else(!is.na(serotype), "NVT", "None")))

#delete individual fup spn datasets to remain with a merged fup dataset
rm(list= ls()[! (ls() %in% c( "spn_baseline", "spn_fup", "spn_overall"))])

#====================================================================

#formulate a new dataset for SIS markov modelling
spn_model <- 
  spn_fup %>%
  mutate(pid = pid,
         
         dys = day,
         
         vday = as.integer(vday),
         
         state = as.integer(if_else(serogroup == "None", 1L, 
                                    if_else(serogroup == "VT", 2L, 3L))),
         
         sex = as.factor(sex),
         
         agegp = as.factor(if_else(age >= 18 & age <= 24, "18-24", 
                                   if_else(age > 24 & age <= 34, "25-34", 
                                           if_else(age > 34 & age <= 45, "35-45", NA_character_)))),
         
         hiv = as.factor(hiv),
         
         artdur = as.factor(if_else(artdur <=3, "short", 
                                    if_else(artdur > 3 & artdur <= 25, "long", NA_character_))),
         
         cd4 = as.factor(if_else(cd4 <= 200, "low", #categorised similarly for HIVpos and HIVneg?
                                 if_else(cd4 > 200 & cd4 <= 500, "medium", 
                                         if_else(cd4 > 500, "high", NA_character_)))),
         
         dens = as.factor(if_else(carr == 0, "none",
                                  if_else(dens <= 1675, "low", #1st quantile
                                          if_else(dens > 1675 & dens <= 1206000000, "high", NA_character_)))), #2nd quantile

         
         season = as.factor(if_else(month(date) >= 5 & month(date) <= 10, "cooldry", 
                                    if_else(month(date) <= 4, "hotwet", 
                                            if_else(month(date) > 10, "hotwet", NA_character_)))),
         
         abx = as.factor(if_else(abx == 1, "taken", 
                                 if_else(abx == 0, "nottaken", NA_character_))),
         
         ses = as.factor(if_else(ses <= 3, "low", 
                                 if_else(ses > 3 & ses <=15, "high", NA_character_))),
         
         nochild = as.factor(if_else(nochild == 1, "1child", "2+child")),
         
         sstate = if_else(is.na(serotype), 1L,
                          if_else(serotype == "3", 2L,
                                  if_else(serotype == "4", 3L,
                                  if_else(serotype == "6A" | serotype == "6B", 4L,
                                          if_else(serotype == "7B/C", 5L,
                                                  if_else(serotype == "8", 6L,
                                                          if_else(serotype == "9V", 7L,
                                                                  if_else(serotype == "11", 8L,
                                                                          if_else(serotype == "15", 9L,
                                                                                  if_else(serotype == "17", 10L,
                                                                                          if_else(serotype == "19B/C", 11L,
                                                                                                  if_else(serotype == "19A" | serotype == "19F", 12L,
                                                                                                          if_else(serotype == "23A/B", 13L,
                                                                                                                  if_else(serotype %in% c("5", "7F", "14", "18C", "23F"), 14L, 15L))))))))))))))) %>%
  select(pid, dys, vday, state, sstate, sex, agegp, hiv, artdur, cd4, dens, season, abx, ses, nochild)

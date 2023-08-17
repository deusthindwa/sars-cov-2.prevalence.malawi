#written by Deus
#01/08/2022
#pneumococcal carriage ad serotype dynamics in HIV-infected adults in the infant PCV era

#====================================================================

#import spn patient level baseline and follow up data
spn_overall <- 
  rio::import(here("data", "nasomuneData.xlsx")) %>%
  dplyr::select(pid, Status, Date, Visit_no, carriage, PNS_serotype, serogroup, Density, sex, no_und5, age_inc, locna, ctx, abx_date, possession_comfort) %>%
  dplyr::rename( "hiv" = "Status", "sdate" = "Date", "vno" = "Visit_no", "state" = "carriage", "serotype" = "PNS_serotype", "dens" = "Density", "sex" = "sex", "nochild" = "no_und5", "agegp" = "age_inc", "loc" = "locna", "ses" = "possession_comfort") %>%
  dplyr::mutate(hiv = if_else(hiv == "HIV-", "HIV-",
                              if_else(hiv == "HIV <3m", "ART <3 month",
                              if_else(hiv == "HIV>1year", "ART >1 year", NA_character_)))) %>%
  dplyr::arrange(pid, sdate) %>%
  dplyr::mutate(sdate = date(sdate)) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(dys = as.integer(sdate - lag(sdate, default = first(sdate))),
                dys = cumsum(dys)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(vday = as.integer(vno),
                state = as.integer(state),
                serotype = if_else(serotype == "NVT[C+]" | serotype == "NVT[D+]" | serotype == "NVT[F+]" | serotype == "NVT[G+]" | serotype == "NVT[H+]" | serotype == "NVT[I+]", "NVT", serotype),
                serogroup = as.factor(serogroup),
                dens2 = as.factor(if_else(state == 0, "none",
                                         if_else(dens < 18420, "low", #less than median value
                                                 if_else(dens >= 18420 & dens <= 5.695e+09, "high", NA_character_)))),
                sex = as.factor(sex),
                nochild = as.factor(if_else(nochild == 1, "one child", 
                                            if_else(nochild > 1, "2+ children", NA_character_))),
                age = as.integer(agegp),
                agegp = as.factor(if_else(age >= 18 & age <= 34, "18-24y", 
                                                  if_else(age > 34, "35+y", NA_character_))),
                abx = as.factor(if_else(!is.na(abx_date), "yes", 
                              if_else(is.na(abx_date), ctx, NA_character_))),
                ses = as.factor(if_else(ses == "yes", "high",
                              if_else(ses == "no", "low", NA_character_)))) %>%
  dplyr::select(pid, dys, vday, state, serotype, serogroup, dens, dens2, sex, age, agegp, hiv, abx, nochild, loc)

#====================================================================

#formulate a new dataset for baseline
spn_baseline <- 
  spn_overall %>%
  dplyr::filter(dys == 0)
               
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

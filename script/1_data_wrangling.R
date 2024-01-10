#written by Deus
#01/08/2022
#pneumococcal carriage ad serotype dynamics in HIV-infected adults in the infant PCV era

#====================================================================

#import location data of blantyre areas and map where samples were collected
spn_loc <- 
  rio::import(here("data", "mapSample.xlsx"))

spn_btmap <- 
  rio::import(here("data", "spn_btmap.csv"))

#import spn patient level baseline and follow up data
spn_overall <- 
  rio::import(here("data", "nasomuneData.xlsx")) %>%
  dplyr::select(pid, status, data_date, visit_no, sex, age_inc, locna, no_und5, ctx, pns_carriage, serotype, serogroup, density, watch:netyn, mattress:otherelec) %>%
  dplyr::rename( "hiv" = "status", "sdate" = "data_date", "vno" = "visit_no", "state" = "pns_carriage", "dens" = "density", "agegp" = "age_inc",  "loc" = "locna", "nochild" = "no_und5") %>%
  dplyr::filter(!is.na(state)) %>%
  mutate(across(watch:otherelec, ~ if_else(.x == "Yes", 1, if_else(.x == "No", 0, NA_integer_)))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ses = sum(watch:otherelec, na.rm = TRUE))

  dplyr::mutate(hiv = if_else(hiv == "HIV-", "HIV neg",
                              if_else(hiv == "HIV <3m", "ART<3m",
                              if_else(hiv == "HIV>1year", "ART>1y", NA_character_))),
                hiv = factor(hiv, levels(factor(hiv))[c(3,1,2)])) %>%
  dplyr::arrange(pid, sdate) %>%
  dplyr::mutate(sdate = date(sdate)) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(dys = as.integer(sdate - lag(sdate, default = first(sdate))),
                dys = cumsum(dys)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(vday = as.integer(vno),
                state = as.integer(state),
                serotype = if_else(serotype == "NVT[C+]" | serotype == "NVT[D+]" | serotype == "NVT[F+]" | serotype == "NVT[G+]" | serotype == "NVT[H+]" | serotype == "NVT[I+]", "NVT", serotype),
                serotype = if_else(serotype == "NVT", "uNVT", serotype),
                serogroup = as.factor(serogroup),
                dens2 = as.factor(if_else(state == 0, "none",
                                         if_else(dens < 18420, "few", #less than median value
                                                 if_else(dens >= 18420 & dens <= 5.695e+09, "high", NA_character_)))),
                sex = factor(sex, levels(factor(sex))[c(2,1)]),
                nochild = as.factor(if_else(nochild == 1, "1child", 
                                            if_else(nochild > 1, "2+child", NA_character_))),
                age = as.integer(agegp),
                agegp = as.factor(if_else(age >= 18 & age <= 34, "18-34y", 
                                                  if_else(age > 34, "35-44y", NA_character_))))

#====================================================================
#====================================================================

#formulate a new dataset for baseline
spn_baseline <- 
  spn_overall %>%
  dplyr::filter(dys == 0) %>%
  group_by(pid) %>%
  rowwise() %>% 
  mutate(abx = sum(amox, fluclox, cef, genta, aug, doxy, metro, cipro, erythro, cot, na.rm = TRUE), 
         ses = sum(watch, radio, bank, iron, sew, mobile, cd, fanelec, netyn, netnum, mattress, bed, bike, motor, car, tv, fridge, otherelec, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(abxc = as.factor(if_else(abx == 0, "no", 
                                 if_else(abx >= 1, "yes", NA_character_))),
         sesc = as.factor(if_else(ses == 0, "no", 
                                  if_else(ses >= 1, "yes", NA_character_)))) %>%

dplyr::select(pid, dys, vday, state, serotype, sex, agegp, hiv, dens, nochild, serogroup, loc, abx, ses)


#====================================================================
#====================================================================

#formulate a new dataset for follow up
spn_model <- 
  spn_overall %>%
  dplyr::select(pid, sdate, dys, vday, state, serotype, serogroup, dens, dens2, sex, age, agegp, hiv, nochild, loc) %>%
  mutate(pid = pid,
         dys = as.integer(dys),
         vday = as.integer(vday),
         state = as.integer(if_else(serogroup == "None", 1L, 
                                    if_else(serogroup == "VT", 2L, 3L))),
         wstate = as.integer(if_else(state == 1L, 1L,
                                     if_else(state == 2L | state == 3L, 2L, NA_integer_))),
         dens0 = dens,
         dens = as.factor(dens2),
         sex = as.factor(sex),
         age = as.integer(age),
         agegp = as.factor(agegp),
         hiv = as.factor(hiv),
         season = as.factor(if_else(month(sdate) >= 5 & month(sdate) <= 10, "cooldry", 
                                    if_else(month(sdate) <= 4, "hotwet", 
                                            if_else(month(sdate) > 10, "hotwet", NA_character_)))),
         nochild = as.factor(nochild)) %>%
  dplyr::filter(pid != "NAS0057" & pid != "NAS1188") %>% #remove obs with single observation as adds no information to likelihood
  dplyr::mutate(serogroup = if_else(state == 2L, "VT",
                                    if_else(state == 3L, "NVT", 
                                            if_else(state == 1L, "None", NA_character_))),
                hivst = if_else(hiv == "HIV-" & serogroup == "VT", "VT,HIV-",
                                if_else(hiv == "HIV-" & serogroup == "NVT", "NVT,HIV-",
                                        if_else(hiv == "ART<3m" & serogroup == "VT", "VT,ART<3m",
                                                if_else(hiv == "ART<3m" & serogroup == "NVT", "NVT,ART<3m",
                                                        if_else(hiv == "ART>1y" & serogroup == "VT", "VT,ART>1y",
                                                                if_else(hiv == "ART>1y" & serogroup == "NVT", "NVT,ART>1y", NA_character_)))))),
                sstate = if_else(is.na(serotype), 1L, #pcv13 serotypes (1, 3, 4, 5, 6A, 6B, 7F, 9V, 14, 19A, 19F, 18C, and 23F)
                                                                 if_else(serotype == "11A/B/C/D/F", 2L,
                                                                         if_else(serotype == "7A/B/C", 3L,
                                                                                 if_else(serotype == "3", 4L,
                                                                                                 if_else(serotype == "19F", 5L,
                                                                                                         if_else(serotype == "17A/F", 6L,
                                                                                                                 if_else(serotype == "15A/B/C/F", 7L,
                                                                                                                         if_else(serotype %in% c("6D", "6C", "6C/D"), 8L,
                                                                                                                                 if_else(serotype == "23A/B", 9L,
                                                                                                                                         if_else(serotype == "10A/B/C/F", 10L,
                                                                                                                                                 if_else(serotype %in% c("6A", "6A/B"), 11L,
                                                                                                                                                                 if_else(serotype %in% c("1", "4", "14", "18C", "19A", "23F"), 12L,
                                                                                                                                                                         if_else(serotype == "uNVT", 13L, 14L)))))))))))))) %>%
                                                                                                                         
  dplyr::select(pid, dys, vday, wstate, state, sstate, serotype, sex, agegp, hiv, dens0, dens, nochild, serogroup, hivst, loc)


#written by Deus
#13/02/2024
#pneumococcal carriage and serotype dynamics by adult HIV status a mature PCV program

#====================================================================

#import location data of Blantyre areas and map where samples were collected
spn_loc <- 
  rio::import(here("data", "mapSample.xlsx"))

spn_map <- 
  rio::import(here("data", "spn_btmap.csv"))

#import spn patient level baseline and follow up data
spn_all <- 
  rio::import(here("data", "nasomuneData.xlsx")) %>%
  dplyr::select(pid, status, data_date, visit_no, sex, age_inc, locna, no_und5, ctx, pns_carriage, serotype, serogroup, density, watch:netyn, mattress:otherelec, abx_new_amo:abx_new_cot) %>%
  dplyr::rename( "hiv" = "status", "sdate" = "data_date", "vno" = "visit_no", "state" = "pns_carriage", "dens" = "density", "agegp" = "age_inc",  "loc" = "locna", "nochild" = "no_und5") %>%
  dplyr::filter(!is.na(state)) %>%

#antibiotic and social economic status
  dplyr::mutate(across(c(watch:otherelec), ~ if_else(. == "Yes", 1L, if_else(. == "No", 0L, NA_integer_)))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ses = sum(c(watch, radio, bank, iron, sew, mobile, cd, fanelec, netyn, mattress, bed, bike, motor, car, tv, fridge, otherelec))) %>% 
  dplyr::mutate(abx = sum(c(abx_new_amo, abx_new_flu, abx_new_cef, abx_new_gen, abx_new_aug, abx_new_dox, abx_new_met, abx_new_cip, abx_new_ery, abx_new_cot))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ses = factor(if_else(ses <=7, "low", if_else(ses >7, "high", NA_character_))),
                abx = factor(if_else(abx ==0, "no", if_else(abx >=1, "yes", NA_character_)))) %>%
  dplyr::select(everything(), -c(watch:otherelec), -c(abx_new_amo:abx_new_cot)) %>%

#hiv and art status (main exposure variable)
  dplyr::mutate(hiv = if_else(hiv == "HIV neg", "hiv_neg",
                              if_else(hiv == "HIV <3m", "art_3m",
                              if_else(hiv == "HIV>1year", "art_1y", NA_character_))),
                hiv = factor(hiv, levels(factor(hiv))[c(3,1,2)])) %>% #3-hivNeg,1-hiv1yPos, 2-hiv3mPos

#follow-up details 
  dplyr::arrange(pid, sdate) %>%
  dplyr::mutate(sdate = date(sdate)) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(dys = cumsum(as.integer(sdate - lag(sdate, default = first(sdate))))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(vday = as.integer(vno), state = as.integer(state)) %>%

#serotype and serotype group dynamics       
  dplyr::mutate(serotype = if_else(serotype == "NVT[C+]" | serotype == "NVT[D+]" | serotype == "NVT[F+]" | serotype == "NVT[G+]" | serotype == "NVT[H+]" | serotype == "NVT[I+]", "NVT", serotype),
                serotype = if_else(serotype == "6[B+Q+](b+c-d+)" | serotype == "6[B+Q+]b-c-d-)" | serotype == "19F and NVT[H+]", "NVT", serotype),
                serotype = if_else(serotype == "1 and 11A/B/C/D/F" | serotype == "3 and 11A/B/C/D/F", "11A/B/C/D/F", serotype),
                serogroup = as.factor(serogroup),

#spn carriage density
                densgp = if_else(is.na(dens), "none",
                                           if_else(dens <= 10720 , "low", "high")), #below median value
                densgp = factor(densgp, levels(factor(densgp))[c(2,1,3)]),

#other covariates
                sex = factor(sex, levels(factor(sex))[c(2,1)]),
                nochild = as.factor(if_else(nochild == 1, "1child", 
                                            if_else(nochild > 1, "2+child", NA_character_))),
                age = as.integer(agegp),
                agegp = as.factor(if_else(age >= 18 & age <= 33, "18-33y", #below median age
                                                  if_else(age > 33, "34-44y", NA_character_)))) %>% #above median age
  
  dplyr::select(pid, sdate, dys, vday, state, serotype, serogroup, hiv, sex, age, agegp, dens, densgp, nochild, ses, ctx, abx, loc)

#====================================================================

#baseline dataset
spn_base <- 
  spn_all %>%
  dplyr::filter(dys == 0) %>%
  dplyr::select(pid, dys, vday, state, serotype, serogroup, hiv, sex, age, agegp, dens, densgp, nochild, ses, ctx, loc)

#====================================================================

#follow-up dataset
spn_fup <- 

#configure a follow-up dataset with repeated base values
  left_join(
    spn_all %>%
      dplyr::select(pid, sdate, dys, vday, state, serotype, serogroup, hiv, sex, dens, densgp, abx),
    spn_all %>%
      dplyr::filter(vday ==1) %>% dplyr::select(pid, age, agegp, nochild, ses, ctx, loc)) %>%
  dplyr::select(pid, sdate, dys, vday, state, serotype, serogroup, hiv, sex, age, agegp, dens, densgp, nochild, ses, ctx, abx, loc) %>%

#configure data types for markov modelling
  dplyr::mutate(dys = as.integer(dys),
                vday = as.integer(vday),
                state = as.integer(if_else(serogroup == "None", 1L, if_else(serogroup == "VT", 2L, 3L))),
                wstate = as.integer(if_else(state == 1L, 1L, if_else(state == 2L | state == 3L, 2L, NA_integer_))),
                hiv = as.factor(hiv),
                sex = as.factor(sex),
                age = as.integer(age),
                agegp = as.factor(agegp),
                densgp = as.factor(densgp),
                nochild = as.factor(nochild),
                ses = as.factor(ses), 
                ctx = as.factor(ctx),
                abx = as.factor(abx),
                seas = as.factor(if_else(month(sdate) >= 5 & month(sdate) <= 10, "cooldry", if_else(month(sdate) <= 4, "hotwet", if_else(month(sdate) > 10, "hotwet", NA_character_)))),
  
#configure hiv and carriage status
                hivst = if_else(hiv == "hiv_neg" & serogroup == "VT", "VT,HIV-",
                                if_else(hiv == "hiv_neg" & serogroup == "NVT", "NVT,HIV-",
                                        if_else(hiv == "art_3m" & serogroup == "VT", "VT,ART<3m",
                                                if_else(hiv == "art_3m" & serogroup == "NVT", "NVT,ART<3m",
                                                        if_else(hiv == "art_1y" & serogroup == "VT", "VT,ART>1y",
                                                                if_else(hiv == "art_1y" & serogroup == "NVT", "NVT,ART>1y", NA_character_))))))) %>%

#remove obs with single observation as adds no information to likelihood
  dplyr::filter(pid != "NAS0057" & pid != "NAS1188" & pid != "NAS1535") %>%
  
#pcv13 serotypes (1, 3, 4, 5, 6A, 6B, 7F, 9V, 14, 19A, 19F, 18C, and 23F)
  dplyr::mutate(sstate = if_else(is.na(serotype), 1L,
                           if_else(serotype == "7A/B/C", 2L,
                                   if_else(serotype == "15A/B/C/F", 3L,
                                           if_else(serotype == "3", 4L,
                                                   if_else(serotype == "11A/B/C/D/F", 5L, 
                                                           if_else(serotype == "23A/B", 6L,
                                                                   if_else(serotype == "17A/F", 7L,
                                                                           if_else(serotype == "19F", 8L,
                                                                                   if_else(serotype == "10A/B/C/F", 9L, 
                                                                                           if_else(serotype == "20", 10L,
                                                                                                   if_else(serotype == "6C", 11L,
                                                                                                           if_else(serotype == "9A/L/N", 12L,
                                                                                                                   if_else(serotype == "19A", 13L,
                                                                                                                           if_else(serotype == "6A", 14L, 
                                                                                                                                   if_else(serotype %in% c("1", "4", "9V", "14", "18C", "23F"), 15L,
                                                                                                                                           if_else(serotype == "NVT", 16L, 17L)))))))))))))))))

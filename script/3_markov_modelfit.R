#written by Deus
#01/08/2022
#pneumococcal carriage ad serotype dynamics in HIV-infected adults in the infant PCV era

#====================================================================
# CONTINUOUS-TIME TIME-NONHOMOGENEOUS MARKOV MODEL - OVERALL
#====================================================================

# show transition frequency
spn_model <- arrange(spn_model, pid, dys)
statetable.msm(state, pid, data = spn_model)

#initiate transition intensity matrix Q
spn_Qmatrix <- rbind(c(0.03, 0.12, 0.26),
                     c(0.26, 0.74, 0.00), 
                     c(0.15, 0.00, 0.35))

rownames(spn_Qmatrix) <- c("Clear", "VT_carry", "NVT_carry")
colnames(spn_Qmatrix) <- c("Clear", "VT_carry", "NVT_carry")
spn_Qmatrix

#run the Markov model
spn_modelfit <- msm(state ~ dys, subject = pid, data = spn_model,
                qmatrix = spn_Qmatrix,
                covariates = list("1-2" = ~ hiv + agegp + sex + nochild + ses, "2-1" =~ hiv + agegp + sex + dens + abx, 
                                  "1-3" = ~ hiv + agegp + sex + nochild + ses, "3-1" =~ hiv + agegp + sex + dens + abx),
                pci = c(60, 120, 180, 240),
                opt.method = "bobyqa", control = list(maxfun = 10000000))

#====================================================================
# CONTINUOUS-TIME TIME-NONHOMOGENEOUS MARKOV MODEL - ART DURATION
#====================================================================

# show transition frequency
spn_modeldur <- filter(spn_model, hiv == "HIV+ART+")
spn_modeldur <- arrange(spn_modeldur, pid, dys)
statetable.msm(state, pid, data = spn_modeldur)

#initiate transition intensity matrix Q
spn_Qmatrix <- rbind(c(0.03, 0.12, 0.26),
                     c(0.26, 0.74, 0.00),
                     c(0.15, 0.00, 0.35))

rownames(spn_Qmatrix) <- c("Clear", "VT_carry", "NVT_carry")
colnames(spn_Qmatrix) <- c("Clear", "VT_carry", "NVT_carry")
spn_Qmatrix

#run the Markov model
spn_modelfitdur <- msm(state ~ dys, subject = pid, data = spn_modeldur,
                    qmatrix = spn_Qmatrix,
                    covariates = list("2-1" = ~ agegp + sex + dens + abx + artdur,
                                      "3-1" = ~ agegp + sex + dens + abx + artdur),
                    pci = c(60, 120, 180, 240),
                    opt.method = "bobyqa", control = list(maxfun = 10000000))

#====================================================================
# CONTINUOUS-TIME TIME-NONHOMOGENEOUS MARKOV MODEL - HIV POSITIVE
#====================================================================

# show transition frequency
spn_modelpos <- spn_model %>% filter(hiv == "HIV+ART+")
spn_modelpos <- arrange(spn_modelpos, pid, dys)
statetable.msm(state, pid, data = spn_modelpos)

#initiate transition intensity matrix Q
spn_Qmatrix <- rbind(c(0.03, 0.12, 0.26),
                     c(0.26, 0.74, 0.00), 
                     c(0.15, 0.00, 0.35))

rownames(spn_Qmatrix) <- c("Clear", "VT_carry", "NVT_carry")
colnames(spn_Qmatrix) <- c("Clear", "VT_carry", "NVT_carry")
spn_Qmatrix

#run the Markov model
spn_modelfitpos <- msm(state ~ dys, subject = pid, data = spn_modelpos,
                    qmatrix = spn_Qmatrix,
                    covariates = list("1-2" = ~ agegp + sex + nochild + ses, "2-1" =~ agegp + sex + dens + abx + artdur, 
                                      "1-3" = ~ agegp + sex + nochild + ses, "3-1" =~ agegp + sex + dens + abx + artdur),
                    pci = c(60, 120, 180, 240),
                    opt.method = "bobyqa", control = list(maxfun = 10000000))

#====================================================================
# CONTINUOUS-TIME TIME-NONHOMOGENEOUS MARKOV MODEL - HIV NEGATIVE
#====================================================================

# show transition frequency
spn_modelneg <- spn_model %>% filter(hiv == "HIV-")
spn_modelneg <- arrange(spn_modelneg, pid, dys)
statetable.msm(state, pid, data = spn_modelneg)

#initiate transition intensity matrix Q
spn_Qmatrix <- rbind(c(0.03, 0.12, 0.26),
                     c(0.26, 0.74, 0.00), 
                     c(0.15, 0.00, 0.35))

rownames(spn_Qmatrix) <- c("Clear", "VT_carry", "NVT_carry")
colnames(spn_Qmatrix) <- c("Clear", "VT_carry", "NVT_carry")
spn_Qmatrix

#run the Markov model
spn_modelfitneg <- msm(state ~ dys, subject = pid, data = spn_modelneg,
                       qmatrix = spn_Qmatrix,
                       covariates = list("1-2" = ~ agegp + sex + nochild + ses, "2-1" =~ agegp + sex + dens + abx, 
                                         "1-3" = ~ agegp + sex + nochild + ses, "3-1" =~ agegp + sex + dens + abx),
                       pci = c(60, 120, 180, 240),
                       opt.method = "bobyqa", control = list(maxfun = 10000000))

#====================================================================
# CONTINUOUS-TIME TIME-NONHOMOGENEOUS MARKOV MODEL (SEROTYPE-SPECIFIC)
#====================================================================

# show transition frequency
spn_model <- arrange(spn_model, pid, dys)
statetable.msm(sstate, pid, data = spn_model)

#initiate transition intensity matrix Q
spn_Qmatrix <-  rbind(c(0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26),
                      c(0.26, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26)
)

rownames(spn_Qmatrix) <- c("clear", "3", "4", "6A/6B", "7B/C", "8", "9V", "11", "15", "17", "19B/C", "19A/19F", "23A/B", "oVT", "oNVT")
colnames(spn_Qmatrix) <- c("clear", "3", "4", "6A/6B", "7B/C", "8", "9V", "11", "15", "17", "19B/C", "19A/19F", "23A/B", "oVT", "oNVT")
spn_Qmatrix

# #run the Markov model
spn_modelfitS <- msm(sstate ~ dys, subject = pid, data = spn_model,
                     qmatrix = spn_Qmatrix,
                     #pci = c(60, 120, 180, 240),
                     opt.method = "bobyqa", control = list(maxfun = 100000000))

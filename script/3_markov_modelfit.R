#written by Deus
#01/08/2022
#pneumococcal carriage ad serotype dynamics in HIV-infected adults in the infant PCV era

#====================================================================
# CONTINUOUS-TIME MARKOV MODEL - WHOLE CARRIAGE
#====================================================================

# show transition frequency
spn_model <- arrange(spn_model, pid, dys)
statetable.msm(wstate, pid, data = spn_model)

#initiate transition intensity matrix Q
spn_Qmatrix <- rbind(c(0.03, 0.12),
                     c(0.26, 0.74))

rownames(spn_Qmatrix) <- c("Clear", "Carry")
colnames(spn_Qmatrix) <- c("Clear", "Carry")
spn_Qmatrix

#run the Markov model
spn_modelfit1 <- msm(wstate ~ dys, subject = pid, data = spn_model,
                    qmatrix = spn_Qmatrix,
                    gen.inits = TRUE,
                    covariates = list("1-2" = ~ hiv + agegp + sex + nochild, "2-1" =~ hiv + agegp + sex + dens),
                    opt.method = "bobyqa", control = list(maxfun = 10000000))

#====================================================================
# CONTINUOUS-TIME MARKOV MODEL - SEROTYPE GROUP
#====================================================================

#show transition frequency
spn_model <- arrange(spn_model, pid, dys)
statetable.msm(state, pid, data = spn_model)

#initiate transition intensity matrix Q
spn_Qmatrix <- rbind(c(0.03, 0.12, 0.26),
                     c(0.26, 0.74, 0.00), 
                     c(0.15, 0.00, 0.35))

rownames(spn_Qmatrix) <- c("Clear", "VTcarry", "NVTcarry")
colnames(spn_Qmatrix) <- c("Clear", "VTcarry", "NVTcarry")
spn_Qmatrix

#run the Markov model
spn_modelfit2 <- msm(state ~ dys, subject = pid, data = spn_model,
                qmatrix = spn_Qmatrix,
                gen.inits = TRUE,
                covariates = list("1-2" = ~ hiv + agegp + sex + nochild, "2-1" =~ hiv + agegp + sex + dens, 
                                  "1-3" = ~ hiv + agegp + sex + nochild, "3-1" =~ hiv + agegp + sex + dens),
                opt.method = "bobyqa", control = list(maxfun = 10000000))

#====================================================================
# CONTINUOUS-TIME TIME-NONHOMOGENEOUS MARKOV MODEL (SEROTYPE-SPECIFIC)
#====================================================================

# show transition frequency
spn_model <- arrange(spn_model, pid, dys)
statetable.msm(sstate, pid, data = spn_model)

#initiate transition intensity matrix Q
spn_Qmatrix <-  rbind(c(0.26, 0.26, 0.26, 0.26, 0.26, 0.36, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26),
                      c(0.26, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.36, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26)
)

rownames(spn_Qmatrix) <- c("clear", "11A/B/C/D/F", "7A/B/C", "3", "19F", "17A/F", "15A/B/C/F", "6C/D", "23A/B", "10A/B/C/F", "6A/B", "oVT", "uNVT", "kNVT")
colnames(spn_Qmatrix) <- c("clear", "11A/B/C/D/F", "7A/B/C", "3", "19F", "17A/F", "15A/B/C/F", "6C/D", "23A/B", "10A/B/C/F", "6A/B", "oVT", "uNVT", "kNVT")
spn_Qmatrix

# #run the Markov model
spn_modelfit3 <- msm(sstate ~ dys, subject = pid, data = spn_model,
                     gen.inits = TRUE,
                     qmatrix = spn_Qmatrix,
                     opt.method = "bobyqa", control = list(maxfun = 100000000))

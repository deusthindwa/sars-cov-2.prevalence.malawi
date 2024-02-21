#written by Deus
#13/02/2024
#pneumococcal carriage and serotype dynamics by adult HIV status a mature PCV program

#====================================================================
# CONTINUOUS-TIME MARKOV MODEL - WHOLE CARRIAGE
#====================================================================

# show transition frequency
spn_fup <- arrange(spn_fup, pid, dys)
statetable.msm(wstate, pid, data = spn_fup)

#initiate transition intensity matrix Q
spn_Qmatrix <- rbind(c(0.05, 0.05),
                     c(0.05, 0.05))

rownames(spn_Qmatrix) <- c("Clear", "Carry")
colnames(spn_Qmatrix) <- c("Clear", "Carry")
spn_Qmatrix

#run the Markov model
spn_modelfit1 <- msm(wstate ~ dys, subject = pid, data = spn_fup,
                    qmatrix = spn_Qmatrix,
                    covariates = list("1-2" = ~ hiv + agegp + sex + nochild + ses + seas, "2-1" = ~ hiv + agegp + sex),
                    opt.method = "bobyqa", control = list(maxfun = 10000000))

#====================================================================
# CONTINUOUS-TIME MARKOV MODEL - PCV13 SEROTYPE GROUP
#====================================================================

#show transition frequency
spn_fup <- arrange(spn_fup, pid, dys)
statetable.msm(state, pid, data = spn_fup)

#initiate transition intensity matrix Q
spn_Qmatrix <- rbind(c(0.03, 0.12, 0.26),
                     c(0.26, 0.74, 0.0), 
                     c(0.15, 0.0, 0.35))

rownames(spn_Qmatrix) <- c("Clear", "VTcarry", "NVTcarry")
colnames(spn_Qmatrix) <- c("Clear", "VTcarry", "NVTcarry")
spn_Qmatrix

#run the Markov model
spn_modelfit2 <- msm(state ~ dys, subject = pid, data = spn_fup,
                qmatrix = spn_Qmatrix,
                covariates = list("1-2" = ~ hiv + agegp + sex + nochild + ses + seas, "2-1" = ~ hiv + agegp + sex,
                                  "1-3" = ~ hiv + agegp + sex + nochild + ses + seas, "3-1" = ~ hiv + agegp + sex),
                opt.method = "bobyqa", control = list(maxfun = 10000000))

#====================================================================
# CONTINUOUS-TIME TIME-NONHOMOGENEOUS MARKOV MODEL (SEROTYPE-SPECIFIC)
#====================================================================

# show transition frequency
spn_fup <- arrange(spn_fup, pid, dys)
statetable.msm(sstate, pid, data = spn_fup)

#initiate transition intensity matrix Q
spn_Qmatrix <-  rbind(c(0.26, 0.26, 0.26, 0.26, 0.26, 0.36, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26),
                      c(0.26, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.36, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26, 0.00),
                      c(0.26, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.26)
)

rownames(spn_Qmatrix) <- c("clear", "7A/B/C", "15A/B/C/F", "3", "11A/B/C/D/F", "23A/B", "17A/F", "19F", "10A/B/C/F", "20", "6C", "9A/L/N", "19A", "6A", "oVT", "uNVT", "kNVT")
colnames(spn_Qmatrix) <- c("clear", "7A/B/C", "15A/B/C/F", "3", "11A/B/C/D/F", "23A/B", "17A/F", "19F", "10A/B/C/F", "20", "6C", "9A/L/N", "19A", "6A", "oVT", "uNVT", "kNVT")
spn_Qmatrix

# #run the Markov model
spn_modelfit3 <- msm(sstate ~ dys, subject = pid, data = spn_fup,
                     gen.inits = TRUE,
                     qmatrix = spn_Qmatrix,
                     opt.method = "bobyqa", control = list(maxfun = 100000000))

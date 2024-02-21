#written by Deus
#13/02/2024
#pneumococcal carriage and serotype dynamics by adult HIV status a mature PCV program

#====================================================================

#other pneumococcal dynamics parameters

#estimate the number of episodes
#overall
envisits.msm(spn_modelfit1, fromt = 0, tot = 365.25, covariates = list(hiv = "hiv_neg"), ci = "normal", cl = 0.95)
envisits.msm(spn_modelfit1, fromt = 0, tot = 365.25, covariates = list(hiv = "art_3m"), ci = "normal", cl = 0.95)
envisits.msm(spn_modelfit1, fromt = 0, tot = 365.25, covariates = list(hiv = "art_1y"), ci = "normal", cl = 0.95)

#by serogroup
envisits.msm(spn_modelfit2, fromt = 0, tot = 365.25, ci = "normal", cl = 0.95)
envisits.msm(spn_modelfit2, fromt = 0, tot = 365.25, covariates = list(hiv = "hiv_neg"), ci = "normal", cl = 0.95)
envisits.msm(spn_modelfit2, fromt = 0, tot = 365.25, covariates = list(hiv = "art_3m"), ci = "normal", cl = 0.95)
envisits.msm(spn_modelfit2, fromt = 0, tot = 365.25, covariates = list(hiv = "art_1y"), ci = "normal", cl = 0.95)

#Probability that each state is next
pnext.msm(spn_modelfit2, ci = "normal", cl = 0.95)
pnext.msm(spn_modelfit2, covariates = list(hiv = "hiv_neg"), ci = "normal", cl = 0.95)
pnext.msm(spn_modelfit2, covariates = list(hiv = "art_3m"), ci = "normal", cl = 0.95)
pnext.msm(spn_modelfit2, covariates = list(hiv = "art_1y"), ci = "normal", cl = 0.95)

#For processes with successive periods of recovery and relapse, we may want to forecast the total time spent in each state
totlos.msm(spn_modelfit2, fromt = 0, tot = 365.25, ci = "normal", cl = 0.95)
totlos.msm(spn_modelfit2, fromt = 0, tot = 365.25, covariates = list(hiv = "hiv_neg"), ci = "normal", cl = 0.95)
totlos.msm(spn_modelfit2, fromt = 0, tot = 365.25, covariates = list(hiv = "art_3m"), ci = "normal", cl = 0.95)
totlos.msm(spn_modelfit2, fromt = 0, tot = 365.25, covariates = list(hiv = "art_1y"), ci = "normal", cl = 0.95)

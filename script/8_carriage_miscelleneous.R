#written by Deus
#01/08/2022
#pneumococcal carriage ad serotype dynamics in HIV-infected adults in the infant PCV era

#====================================================================

#other pneumococcal dynamics parameters

#estimate the number of episodes
envisits.msm(spn_modelfit, fromt = 0, tot = 365.25, covariates = list(hiv = "HIV+ART+"), ci = "normal", cl = 0.95)
envisits.msm(spn_modelfit, fromt = 0, tot = 365.25, covariates = list(hiv = "HIV-"), ci = "normal", cl = 0.95)

#Probability that each state is next
pnext.msm(spn_modelfit, covariates = list(hiv = "HIV+ART+"), ci = "normal", cl = 0.95)
pnext.msm(spn_modelfit, covariates = list(hiv = "HIV-"), ci = "normal", cl = 0.95)

#For processes with successive periods of recovery and relapse, we may want to forecast the total time spent healthy or diseased, before death.
totlos.msm(spn_modelfit, fromt = 0, tot = 365.25, covariates = list(hiv = "HIV+ART+"), ci = "normal", cl = 0.95)
totlos.msm(spn_modelfit, fromt = 0, tot = 365.25, covariates = list(hiv = "HIV-"), ci = "normal", cl = 0.95)








  



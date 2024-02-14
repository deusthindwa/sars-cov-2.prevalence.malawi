#written by Deus
#13/02/2024
#pneumococcal carriage and serotype dynamics by adult HIV status a mature PCV program

#====================================================================

#effects of covariates on whole carriage
hazard.msm(spn_modelfit1, hazard.scale = 1, cl = 0.95)

#effects of covariates on vaccine-serotype group carriage
hazard.msm(spn_modelfit2, hazard.scale = 1, cl = 0.95)

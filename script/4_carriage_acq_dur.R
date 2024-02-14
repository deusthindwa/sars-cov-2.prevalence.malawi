#written by Deus
#13/02/2024
#pneumococcal carriage and serotype dynamics by adult HIV status a mature PCV program

#====================================================================

#acquisition probability of whole carriage
modela <- pmatrix.msm(spn_modelfit1, t = 1, covariates = list(hiv = "hiv_neg"), ci = "normal", cl = 0.95)
modelb <- pmatrix.msm(spn_modelfit1, t = 1, covariates = list(hiv = "art_3m"), ci = "normal" , cl = 0.95)
modelc <- pmatrix.msm(spn_modelfit1, t = 1, covariates = list(hiv = "art_1y"), ci = "normal" , cl = 0.95)

pneumo0 <- data.frame("hivst" = c("HIV-", "ART<3m", "ART>1y"))
pneumo0$carry <- c(modela$estimates[1,2], modelb$estimates[1,2], modelc$estimates[1,2])
pneumo0$Lcarry <- c(modela$L[1,2], modelb$L[1,2], modelc$L[1,2])
pneumo0$Ucarry <- c(modela$U[1,2], modelb$U[1,2], modelc$U[1,2])

a <- 
  pneumo0 %>%
  ggplot() +
  geom_point(aes(hivst, carry, shape = hivst), color = "black", size = 3, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(hivst,  ymin = Lcarry, ymax = Ucarry), color = "black", width = 0, size = 1.2, position = position_dodge2(width = 0.5)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04), limits = c(0, 0.045)) + 
  labs(title = "(a)", x = "", y = "Daily carriage acquisition probability") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(shape = guide_legend(title = ""), color = "none") +
  theme(legend.text = element_text(size = 12), legend.position = "none") +
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#acquisition probability by HIV status overall
modela <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "hiv_neg"), ci = "normal", cl = 0.95)
modelb <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_3m"), ci = "normal" , cl = 0.95)
modelc <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_1y"), ci = "normal" , cl = 0.95)

pneumo0 <- data.frame("hivst" = c("VT,HIV-", "VT,ART<3m", "VT,ART>1y"))
pneumo0$carry <- c(modela$estimates[1,2], modelb$estimates[1,2], modelc$estimates[1,2])
pneumo0$Lcarry <- c(modela$L[1,2], modelb$L[1,2], modelc$L[1,2])
pneumo0$Ucarry <- c(modela$U[1,2], modelb$U[1,2], modelc$U[1,2])

pneumo1 <- data.frame("hivst" = c("NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y"))
pneumo1$carry <- c(modela$estimates[1,3], modelb$estimates[1,3], modelc$estimates[1,3])
pneumo1$Lcarry <- c(modela$L[1,3], modelb$L[1,3], modelc$L[1,3])
pneumo1$Ucarry <- c(modela$U[1,3], modelb$U[1,3], modelc$U[1,3])

b <- 
  bind_rows(pneumo0, pneumo1) %>%
  ggplot() +
  geom_point(aes(hivst, carry, color = hivst), size = 2, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(hivst,  ymin = Lcarry, ymax = Ucarry, color = hivst), width = 0, size = 1.2, position = position_dodge2(width = 0.5)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04), limits = c(0, 0.045)) + 
  labs(title = "(b)", x = "", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank()) + 
  guides(shape = guide_legend(title = ""), color = "none") +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.8)) +
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#acquisition probability by HIV status and sex
modela <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "hiv_neg", sex = "Male"), ci = "normal", cl = 0.95)
modelb <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_3m", sex = "Male"), ci = "normal", cl = 0.95)
modelc <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_1y", sex = "Male"), ci = "normal", cl = 0.95)
modeld <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "hiv_neg", sex = "Female"), ci = "normal", cl = 0.95)
modele <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_3m", sex = "Female"), ci = "normal", cl = 0.95)
modelf <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_1y", sex = "Female"), ci = "normal", cl = 0.95)

pneumo0 <- data.frame("hivst" = c("VT,HIV-", "VT,ART<3m", "VT,ART>1y", "VT,HIV-", "VT,ART<3m", "VT,ART>1y"), 
                      "sex" = c("Male", "Male", "Male","Female", "Female", "Female"))
pneumo0$carry <- c(modela$estimates[1,2], modelb$estimates[1,2], modelc$estimates[1,2], modeld$estimates[1,2], modele$estimates[1,2], modelf$estimates[1,2])
pneumo0$Lcarry <- c(modela$L[1,2], modelb$L[1,2], modelc$L[1,2], modeld$L[1,2], modele$L[1,2], modelf$L[1,2])
pneumo0$Ucarry <- c(modela$U[1,2], modelb$U[1,2], modelc$U[1,2], modeld$U[1,2], modele$U[1,2], modelf$U[1,2])

pneumo1 <- data.frame("hivst" = c("NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y", "NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y"), 
                      "sex" = c("Male", "Male", "Male","Female", "Female", "Female"))
pneumo1$carry <- c(modela$estimates[1,3], modelb$estimates[1,3], modelc$estimates[1,3], modeld$estimates[1,3], modele$estimates[1,3], modelf$estimates[1,3])
pneumo1$Lcarry <- c(modela$L[1,3], modelb$L[1,3], modelc$L[1,3], modeld$L[1,3], modele$L[1,3], modelf$L[1,3])
pneumo1$Ucarry <- c(modela$U[1,3], modelb$U[1,3], modelc$U[1,3], modeld$U[1,3], modele$U[1,3], modelf$U[1,3])

c <-
  bind_rows(pneumo0, pneumo1) %>%
  mutate(id = paste0(hivst, sex)) %>%
  ggplot() +
  geom_point(aes(id, carry, color = hivst, shape = sex), size = 2, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(id, carry, color = hivst, ymin = Lcarry, ymax = Ucarry), width = 0, size = 1.2, position = position_dodge2(width = 0.5)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04), limits = c(0, 0.045)) + 
  scale_shape_manual(values = c(4, 5)) +
  labs(title = "(c)", x = "", y = "") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank()) + 
  guides(shape = guide_legend(title = ""), color = "none") +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.8)) +
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#acquisition probability by HIV status and age groups
modela <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "hiv_neg", agegp = "18-33y"), ci = "normal", cl = 0.95)
modelb <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_3m", agegp = "18-33y"), ci = "normal", cl = 0.95)
modelc <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_1y", agegp = "18-33y"), ci = "normal", cl = 0.95)
modeld <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "hiv_neg", agegp = "34-44y"), ci = "normal", cl = 0.95)
modele <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_3m", agegp = "34-44y"), ci = "normal", cl = 0.95)
modelf <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_1y", agegp = "34-44y"), ci = "normal", cl = 0.95)

pneumo0 <- data.frame("hivst" = c("VT,HIV-", "VT,ART<3m", "VT,ART>1y", "VT,HIV-", "VT,ART<3m", "VT,ART>1y"), 
                      "agegp" = c("18-33y", "18-33y", "18-33y", "34-44y", "34-44y", "34-44y"))
pneumo0$carry <- c(modela$estimates[1,2], modelb$estimates[1,2], modelc$estimates[1,2], modeld$estimates[1,2], modele$estimates[1,2], modelf$estimates[1,2])
pneumo0$Lcarry <- c(modela$L[1,2], modelb$L[1,2], modelc$L[1,2], modeld$L[1,2],  modele$L[1,2], modelf$L[1,2]) 
pneumo0$Ucarry <- c(modela$U[1,2], modelb$U[1,2], modelc$U[1,2], modeld$U[1,2], modele$U[1,2], modelf$U[1,2])

pneumo1 <- data.frame("hivst" = c("NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y", "NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y"), 
                      "agegp" = c("18-33y", "18-33y", "18-33y", "34-44y", "34-44y", "34-44y"))
pneumo1$carry <- c(modela$estimates[1,3], modelb$estimates[1,3], modelc$estimates[1,3], modeld$estimates[1,3], modele$estimates[1,3], modelf$estimates[1,3])
pneumo1$Lcarry <- c(modela$L[1,3], modelb$L[1,3], modelc$L[1,3], modeld$L[1,3], modele$L[1,3], modelf$L[1,3]); 
pneumo1$Ucarry <- c(modela$U[1,3], modelb$U[1,3], modelc$U[1,3], modeld$U[1,3], modele$U[1,3], modelf$U[1,3])

d <-
  bind_rows(pneumo0, pneumo1) %>%
  mutate(id = paste0(hivst, agegp)) %>%
  ggplot() +
  geom_point(aes(id, carry, color = hivst, shape = agegp), size = 2, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(id, carry, color = hivst, ymin = Lcarry, ymax = Ucarry), width = 0, size = 1.2, position = position_dodge2(width = 0.5)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04), limits = c(0, 0.045)) + 
  scale_shape_manual(values = c(4, 5)) +
  labs(title = "(d)", x = "", y = "") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank()) + 
  guides(shape = guide_legend(title = ""), color = "none") +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.8)) +
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#acquisition probability by HIV status and number of under 5 children
modela <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "hiv_neg", nochild = "1child"), ci = "normal", cl = 0.95)
modelb <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_3m", nochild = "1child"), ci = "normal", cl = 0.95)
modelc <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_1y", nochild = "1child"), ci = "normal", cl = 0.95)
modeld <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "hiv_neg", nochild = "2+child"), ci = "normal", cl = 0.95)
modele <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_3m", nochild = "2+child"), ci = "normal", cl = 0.95)
modelf <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_1y", nochild = "2+child"), ci = "normal", cl = 0.95)

pneumo0 <- data.frame("hivst" = c("VT,HIV-", "VT,ART<3m", "VT,ART>1y", "VT,HIV-", "VT,ART<3m", "VT,ART>1y"), 
                      "nochild" = c("1 child", "1 child", "1 child", "2+ children", "2+ children", "2+ children"))
pneumo0$carry <- c(modela$estimates[1,2], modelb$estimates[1,2], modelc$estimates[1,2], modeld$estimates[1,2], modele$estimates[1,2], modelf$estimates[1,2])
pneumo0$Lcarry <- c(modela$L[1,2], modelb$L[1,2], modelc$L[1,2], modeld$L[1,2], modele$L[1,2], modelf$L[1,2]); 
pneumo0$Ucarry <- c(modela$U[1,2], modelb$U[1,2], modelc$U[1,2], modeld$U[1,2], modele$U[1,2], modelf$U[1,2])

pneumo1 <- data.frame("hivst" = c("NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y", "NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y"), 
                      "nochild" = c("1 child", "1 child", "1 child", "2+ children", "2+ children", "2+ children"))
pneumo1$carry <- c(modela$estimates[1,3], modelb$estimates[1,3], modelc$estimates[1,3], modeld$estimates[1,3], modele$estimates[1,3], modelf$estimates[1,3])
pneumo1$Lcarry <- c(modela$L[1,3], modelb$L[1,3], modelc$L[1,3], modeld$L[1,3], modele$L[1,3], modelf$L[1,3])
pneumo1$Ucarry <- c(modela$U[1,3], modelb$U[1,3], modelc$U[1,3], modeld$U[1,3], modele$U[1,3], modelf$U[1,3])

e <-
  bind_rows(pneumo0, pneumo1) %>%
  mutate(id = paste0(hivst, nochild)) %>%
  ggplot() + 
  geom_point(aes(id, carry, color = hivst, shape = nochild), size = 2, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(id, carry, color = hivst, ymin = Lcarry, ymax = Ucarry), width = 0, size = 1.2, position = position_dodge2(width = 0.5)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04), limits = c(0, 0.045)) + 
  scale_shape_manual(values = c(4, 5)) +
  labs(title = "(e)", x = "", y = "") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank()) + 
  guides(shape = guide_legend(title = ""), color = "none") + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.8)) + 
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm")) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#acquisition probability by HIV status and social economic status
modela <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "hiv_neg", ses = "low"), ci = "normal", cl = 0.95)
modelb <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_3m", ses = "low"), ci = "normal", cl = 0.95)
modelc <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_1y", ses = "low"), ci = "normal", cl = 0.95)
modeld <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "hiv_neg", ses = "high"), ci = "normal", cl = 0.95)
modele <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_3m", ses = "high"), ci = "normal", cl = 0.95)
modelf <- pmatrix.msm(spn_modelfit2, t = 1, covariates = list(hiv = "art_1y", ses = "high"), ci = "normal", cl = 0.95)

pneumo0 <- data.frame("hivst" = c("VT,HIV-", "VT,ART<3m", "VT,ART>1y", "VT,HIV-", "VT,ART<3m", "VT,ART>1y"), 
                      "ses" = c("Low SES", "Low SES", "Low SES", "High SES", "High SES", "High SES"))
pneumo0$carry <- c(modela$estimates[1,2], modelb$estimates[1,2], modelc$estimates[1,2], modeld$estimates[1,2], modele$estimates[1,2], modelf$estimates[1,2])
pneumo0$Lcarry <- c(modela$L[1,2], modelb$L[1,2], modelc$L[1,2], modeld$L[1,2], modele$L[1,2], modelf$L[1,2]); 
pneumo0$Ucarry <- c(modela$U[1,2], modelb$U[1,2], modelc$U[1,2], modeld$U[1,2], modele$U[1,2], modelf$U[1,2])

pneumo1 <- data.frame("hivst" = c("NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y", "NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y"), 
                      "ses" = c("Low SES", "Low SES", "Low SES", "High SES", "High SES", "High SES"))
pneumo1$carry <- c(modela$estimates[1,3], modelb$estimates[1,3], modelc$estimates[1,3], modeld$estimates[1,3], modele$estimates[1,3], modelf$estimates[1,3])
pneumo1$Lcarry <- c(modela$L[1,3], modelb$L[1,3], modelc$L[1,3], modeld$L[1,3], modele$L[1,3], modelf$L[1,3])
pneumo1$Ucarry <- c(modela$U[1,3], modelb$U[1,3], modelc$U[1,3], modeld$U[1,3], modele$U[1,3], modelf$U[1,3])

f <-
  bind_rows(pneumo0, pneumo1) %>%
  mutate(id = paste0(hivst, ses)) %>%
  ggplot() + 
  geom_point(aes(id, carry, color = hivst, shape = ses), size = 2, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(id, carry, color = hivst, ymin = Lcarry, ymax = Ucarry), width = 0, size = 1.2, position = position_dodge2(width = 0.5)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04), limits = c(0, 0.045)) + 
  scale_shape_manual(values = c(4, 5)) +
  labs(title = "(f)", x = "", y = "") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank()) + 
  guides(shape = guide_legend(title = ""), color = "none") + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.8)) + 
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm")) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#====================================================================
#====================================================================

#carriage duration of whole carriage (inverse clearance rate) by HIV status 
modela <- qmatrix.msm(spn_modelfit1, covariates = list(hiv = "hiv_neg"), ci = "normal", cl = 0.95)
modelb <- qmatrix.msm(spn_modelfit1, covariates = list(hiv = "art_3m"), ci = "normal" , cl = 0.95)
modelc <- qmatrix.msm(spn_modelfit1, covariates = list(hiv = "art_1y"), ci = "normal" , cl = 0.95)

pneumo0 <- data.frame("hivst" = c("HIV-", "ART<3m", "ART>1y"))
pneumo0$clear <- c(modela$estimates[2,1], modelb$estimates[2,1], modelc$estimates[2,1])
pneumo0$Lclear <- c(modela$L[2,1], modelb$L[2,1], modelc$L[2,1])
pneumo0$Uclear <- c(modela$U[2,1], modelb$U[2,1], modelc$U[2,1])

g <- 
  pneumo0 %>%
  mutate(clearx = 1/clear, Lclearx = 1/Uclear, Uclearx = 1/Lclear) %>%
  ggplot() +
  geom_point(aes(hivst, clearx, shape = hivst), color = "black", size = 3, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(hivst,  ymin = Lclearx, ymax = Uclearx), color = "black", width = 0, size = 1.2, position = position_dodge2(width = 0.5)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  scale_y_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30, 35, 40)) + 
  coord_cartesian(ylim = c(0, 40)) +
  labs(title = "(g)", x = "", y = "Average carriage duration (days)") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(shape = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 0), legend.position = "none") +
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#carriage duration (inverse clearance rate) by HIV status overall
modela <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "hiv_neg"), ci = "normal", cl = 0.95)
modelb <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "art_3m"), ci = "normal" , cl = 0.95)
modelc <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "art_1y"), ci = "normal" , cl = 0.95)

pneumo0 <- data.frame("hivst" = c("VT,HIV-", "VT,ART<3m", "VT,ART>1y"))
pneumo0$clear <- c(modela$estimates[2,1], modelb$estimates[2,1], modelc$estimates[2,1])
pneumo0$Lclear <- c(modela$L[2,1], modelb$L[2,1], modelc$L[2,1])
pneumo0$Uclear <- c(modela$U[2,1], modelb$U[2,1], modelc$U[2,1])

pneumo1 <- data.frame("hivst" = c("NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y"))
pneumo1$clear <- c(modela$estimates[3,1], modelb$estimates[3,1], modelc$estimates[3,1])
pneumo1$Lclear <- c(modela$L[3,1], modelb$L[3,1], modelc$L[3,1])
pneumo1$Uclear <- c(modela$U[3,1], modelb$U[3,1], modelc$U[3,1])

h <- 
  bind_rows(pneumo0, pneumo1) %>%
  mutate(clearx = 1/clear, Lclearx = 1/Uclear, Uclearx = 1/Lclear) %>%
  ggplot() +
  geom_point(aes(hivst, clearx, color = hivst), size = 2, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(hivst,  ymin = Lclearx, ymax = Uclearx, color = hivst), width = 0, size = 1.2, position = position_dodge2(width = 0.5)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  scale_y_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30, 35, 40)) + 
  coord_cartesian(ylim = c(0, 40)) +
  labs(title = "(h)", x = "", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank()) + 
  guides(shape = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 0), legend.position = "none") +
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#carriage duration (inverse clearance rate) by HIV status and sex
modela <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "hiv_neg", sex = "Male"), ci = "normal", cl = 0.95)
modelb <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "art_3m", sex = "Male"), ci = "normal", cl = 0.95)
modelc <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "art_1y", sex = "Male"), ci = "normal", cl = 0.95)
modeld <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "hiv_neg", sex = "Female"), ci = "normal", cl = 0.95)
modele <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "art_3m", sex = "Female"), ci = "normal", cl = 0.95)
modelf <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "art_1y", sex = "Female"), ci = "normal", cl = 0.95)

pneumo0 <- data.frame("hivst" = c("VT,HIV-", "VT,ART<3m", "VT,ART>1y", "VT,HIV-", "VT,ART<3m", "VT,ART>1y"), 
                      "sex" = c("Male", "Male", "Male", "Female", "Female", "Female"))
pneumo0$clear <- c(modela$estimates[2,1], modelb$estimates[2,1], modelc$estimates[2,1], modeld$estimates[2,1], modele$estimates[2,1], modelf$estimates[2,1])
pneumo0$Lclear <- c(modela$L[2,1], modelb$L[2,1], modelc$L[2,1], modeld$L[1,2], modele$L[2,1], modelf$L[1,2])
pneumo0$Uclear <- c(modela$U[2,1], modelb$U[2,1], modelc$U[2,1], modeld$U[2,1], modele$U[2,1], modelf$U[2,1])

pneumo1 <- data.frame("hivst" = c("NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y", "NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y"),
                      "sex" = c("Male", "Male", "Male", "Female", "Female", "Female"))
pneumo1$clear <- c(modela$estimates[3,1], modelb$estimates[3,1], modelc$estimates[3,1], modeld$estimates[3,1], modele$estimates[3,1], modelf$estimates[3,1])
pneumo1$Lclear <- c(modela$L[3,1], modelb$L[3,1], modelc$L[3,1], modeld$L[3,1], modele$L[3,1], modelf$L[3,1])
pneumo1$Uclear <- c(modela$U[3,1], modelb$U[3,1], modelc$U[3,1], modeld$U[3,1], modele$U[3,1], modelf$U[3,1])

i <-
  bind_rows(pneumo0, pneumo1) %>%
  mutate(clearx = 1/clear, Lclearx = 1/Uclear, Uclearx = 1/Lclear) %>%
  mutate(id = paste0(hivst, sex)) %>%
  ggplot() +
  geom_point(aes(id, clearx, color = hivst, shape = sex), size = 2, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(id, clearx, color = hivst, ymin = Lclearx, ymax = Uclearx), width = 0, size = 1.2, position = position_dodge2(width = 0.5)) +
  scale_shape_manual(values = c(4, 5)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  scale_y_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30, 35, 40)) + 
  coord_cartesian(ylim = c(0, 40)) +
  labs(title = "(i)", x = "", y = "") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank()) + 
  guides(shape = guide_legend(title = ""), color = "none") +
  theme(legend.text = element_text(size = 12), legend.position = c(0.3, 0.8)) +
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#carriage duration (inverse clearance rate) by HIV status and age groups
modela <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "hiv_neg", agegp = "18-33y"), ci = "normal", cl = 0.95)
modelb <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "art_3m", agegp = "18-33y"), ci = "normal", cl = 0.95)
modelc <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "art_1y", agegp = "18-33y"), ci = "normal", cl = 0.95)
modeld <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "hiv_neg", agegp = "34-44y"), ci = "normal", cl = 0.95)
modele <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "art_3m", agegp = "34-44y"), ci = "normal", cl = 0.95)
modelf <- qmatrix.msm(spn_modelfit2, covariates = list(hiv = "art_1y", agegp = "34-44y"), ci = "normal", cl = 0.95)

pneumo0 <- data.frame("hivst" = c("VT,HIV-", "VT,ART<3m", "VT,ART>1y", "VT,HIV-", "VT,ART<3m", "VT,ART>1y"), 
                      "agegp" = c("18-33y", "18-33y", "18-33y", "34-44y", "34-44y", "34-44y"))
pneumo0$clear <- c(modela$estimates[2,1], modelb$estimates[2,1], modelc$estimates[2,1], modeld$estimates[2,1], modele$estimates[2,1], modelf$estimates[2,1])
pneumo0$Lclear <- c(modela$L[2,1], modelb$L[2,1], modelc$L[2,1], modeld$L[2,1],  modele$L[2,1], modelf$L[2,1]); 
pneumo0$Uclear <- c(modela$U[2,1], modelb$U[2,1], modelc$U[2,1], modeld$U[2,1], modele$U[2,1], modelf$U[2,1])

pneumo1 <- data.frame("hivst" = c("NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y", "NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y"), 
                      "agegp" = c("18-33y", "18-33y", "18-33y", "34-44y", "34-44y", "34-44y"))
pneumo1$clear <- c(modela$estimates[3,1], modelb$estimates[3,1], modelc$estimates[3,1], modeld$estimates[3,1], modele$estimates[3,1], modelf$estimates[3,1])
pneumo1$Lclear <- c(modela$L[3,1], modelb$L[3,1], modelc$L[3,1], modeld$L[3,1], modele$L[3,1], modelf$L[3,1]); 
pneumo1$Uclear <- c(modela$U[3,1], modelb$U[3,1], modelc$U[3,1], modeld$U[3,1], modele$U[3,1], modelf$U[3,1])

j <-
  bind_rows(pneumo0, pneumo1) %>%
  mutate(clearx = 1/clear, Lclearx = 1/Uclear, Uclearx = 1/Lclear) %>%
  mutate(id = paste0(hivst, agegp)) %>%
  ggplot() +
  geom_point(aes(id, clearx, color = hivst, shape = agegp), size = 2, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(id, clearx, color = hivst, ymin = Lclearx, ymax = Uclearx), width = 0, size = 1.2, position = position_dodge2(width = 0.5)) +
  scale_shape_manual(values = c(4, 5, 1)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') +
  scale_y_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30, 35, 40)) + 
  coord_cartesian(ylim = c(0, 40)) +
  labs(title = "(j)", x = "", y = "") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank()) + 
  guides(shape = guide_legend(title = ""), color = "none") +
  theme(legend.text = element_text(size = 12), legend.position = c(0.3, 0.8)) +
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm")) + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#====================================================================
#====================================================================

#plot a blank plot with legend for HIV status
modela <- pmatrix.msm(spn_modelfit1, t = 1, ci = "normal", cl = 0.95)
modelb <- pmatrix.msm(spn_modelfit1, t = 1, ci = "normal" , cl = 0.95)
modelc <- pmatrix.msm(spn_modelfit1, t = 1, ci = "normal" , cl = 0.95)

pneumo0 <- data.frame("hivst" = c("HIV-; adults living\nwithout HIV", "ART<3m; adults with HIV\non ART <3 months", "ART>1y; adults with HIV\non ART >1 year"))
pneumo0$carry <- c(modela$estimates[1,2], modelb$estimates[1,2], modelc$estimates[1,2])
pneumo0$Lcarry <- c(modela$L[1,2], modelb$L[1,2], modelc$L[1,2])
pneumo0$Ucarry <- c(modela$U[1,2], modelb$U[1,2], modelc$U[1,2])

k <- 
  pneumo0 %>%
  ggplot() +
  geom_point(aes(hivst, carry, shape = hivst), color = "black", size = 3, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(hivst,  ymin = Lcarry, ymax = Ucarry), color = "black", width = 0, size = 1.2, position = position_dodge2(width = 0.5)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  ylim(0.8, 1) +
  labs(title = "", x = "", y = "") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.5, 0.5), legend.key.size = unit(2.5, 'lines'), legend.title = element_blank()) +
  theme(plot.margin = grid::unit(c(0.5,0,0,0), "mm")) +
  theme(panel.border = element_rect(colour = "white", fill = NA, size = 1))

#plot a blank plot with legend for vaccine-type and HIV status
modela <- pmatrix.msm(spn_modelfit2, t = 1, ci = "normal", cl = 0.95)
modelb <- pmatrix.msm(spn_modelfit2, t = 1, ci = "normal" , cl = 0.95)
modelc <- pmatrix.msm(spn_modelfit2, t = 1, ci = "normal" , cl = 0.95)

pneumo0 <- data.frame("hivst" = c("VT,HIV-","VT,ART<3m", "VT,ART>1y"))
pneumo0$carry <- c(modela$estimates[1,2], modelb$estimates[1,2], modelc$estimates[1,2])
pneumo0$Lcarry <- c(modela$L[1,2], modelb$L[1,2], modelc$L[1,2])
pneumo0$Ucarry <- c(modela$U[1,2], modelb$U[1,2], modelc$U[1,2])

pneumo1 <- data.frame("hivst" = c("NVT,HIV-", "NVT,ART<3m", "NVT,ART>1y"))
pneumo1$carry <- c(modela$estimates[1,3], modelb$estimates[1,3], modelc$estimates[1,3])
pneumo1$Lcarry <- c(modela$L[1,3], modelb$L[1,3], modelc$L[1,3])
pneumo1$Ucarry <- c(modela$U[1,3], modelb$U[1,3], modelc$U[1,3])

l <- 
  rbind(pneumo0, pneumo1) %>%
  ggplot() +
  geom_line(aes(hivst, carry, color = hivst), size = 2, stroke = 2, position = position_dodge2(width = 0), stat = "identity") +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  ylim(0.8, 1) +
  labs(title = "", x = "", y = "") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.5, 0.5), legend.key.size = unit(2.5, 'lines'), legend.title = element_blank()) +
  theme(plot.margin = grid::unit(c(0.5,0,0,0), "mm")) +
  theme(panel.border = element_rect(colour = "white", fill = NA, size = 1))

#====================================================================

ggsave(here("output", "Fig2_acq_dur.png"),
       plot = ((a | b | c | d | e | f | plot_layout(ncol = 6, width = c(0.5,0.5,1,1,1,1))) / 
                 (g | h | i | j | k | l | plot_layout(ncol = 6, width = c(0.5,0.5,1,1,1,1)))),
       width = 17, height = 8, unit="in", dpi = 300)

#rm(list = grep("spn", ls(), value = TRUE, invert = TRUE))

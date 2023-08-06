#written by Deus
#01/08/2022
#pneumococcal carriage ad serotype dynamics in HIV-infected adults in the infant PCV era


#====================================================================

#overview of hazard ratios for pneumococcal acquisition rates
#spn_modelfit_pos <- qmatrix.msm(spn_modelfitpos, covariates = list(hiv = "HIV+ART+", agegp = "18-24"))

H <- hazard.msm(spn_modelfitpos, hazard.scale = 1, cl = 0.95)

#create an empty acquisition table for forest plot
spn_acqratio <- tibble(
  index = c(1:10),
  carr = c(rep(c("VT", "NVT"), 5)),
  cova = c(rep("agegp", 4), rep("sex", 2), rep("nochild", 2), rep("ses", 2)),
  label = c(rep("Age group \n (18-24y) vs 25-34y", 2), rep("Age group \n (18-24y) vs 35-45y", 2), rep("Sex \n (female) vs male", 2), rep("Number of children \n (one) vs two+ ", 2), rep("Social economic status \n (high) vs low", 2)),
  HR = c(rep(0, 10)),
  HRlci = c(rep(0, 10)),
  HRuci = c(rep(0, 10))
)

#insert HR from acquisition values in the fitted model
l = 1; m = 1; n = 1
for(i in c("agegp25-34", "agegp35-45", "sexMale", "nochild2+child", "seslow")){
  
  for(j in 1:2){
    spn_acqratio[l, 5] = H[[i]][j]
    l = l+1
  }
  
  for(j in 5:6){
    spn_acqratio[m, 6] = H[[i]][j]
    m = m+1
  }
  
  for(j in 9:10){
    spn_acqratio[n, 7] = H[[i]][j]
    n = n+1
  }
  
}

#plot the acquisition hazard ratios
A <- 
  spn_acqratio %>%
  ggplot(aes(y = reorder(label, index), x = HR, color = carr)) +
  geom_errorbarh(aes(xmin = HRlci, xmax = HRuci), height = 0.5, size = 1.2, position = position_dodge2(width = 0.5)) +
  geom_point(shape = 5, size = 2,  stroke = 2, position = position_dodge2(width = 0.5)) +  
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", cex = 0.6, alpha = 0.8) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 7, 9, 10, 11, 12), limits = c(0, 12.5)) + 
  labs(title = "(a) ALWHIV", x = "", y = "") + 
  theme_bw(base_size = 14, base_family = 'American Typewriter') +
  guides(color = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.85)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  scale_color_manual(values=c("#7CAE00", "#C77CFF"))

#====================================================================

#overview of hazard ratios for pneumococcal acquisition rates
H <- hazard.msm(spn_modelfitneg, hazard.scale = 1, cl = 0.95)

#create an empty acquisition table for forest plot
spn_acqratio <- tibble(
  index = c(1:10),
  carr = c(rep(c("VT", "NVT"), 5)),
  cova = c(rep("agegp", 4), rep("sex", 2), rep("nochild", 2), rep("ses", 2)),
  label = c(rep("Age group \n (18-24y) vs 25-34y", 2), rep("Age group \n (18-24y) vs 35-45y", 2), rep("Sex \n (female) vs male", 2), rep("Number of children \n (one) vs two+ ", 2), rep("Social economic status \n (high) vs low", 2)),
  HR = c(rep(0, 10)),
  HRlci = c(rep(0, 10)),
  HRuci = c(rep(0, 10))
)

#insert HR from acquisition values in the fitted model
l = 1; m = 1; n = 1
for(i in c("agegp25-34", "agegp35-45", "sexMale", "nochild2+child", "seslow")){
  
  for(j in 1:2){
    spn_acqratio[l, 5] = H[[i]][j]
    l = l+1
  }
  
  for(j in 5:6){
    spn_acqratio[m, 6] = H[[i]][j]
    m = m+1
  }
  
  for(j in 9:10){
    spn_acqratio[n, 7] = H[[i]][j]
    n = n+1
  }
  
}

#plot the acquisition hazard ratios
B <- 
  spn_acqratio %>%
  ggplot(aes(y = reorder(label, index), x = HR, color = carr)) +
  geom_errorbarh(aes(xmin = HRlci, xmax = HRuci), height = 0.5, size = 1.2, position = position_dodge2(width = 0.5)) +
  geom_point(shape = 5, size = 2,  stroke = 2, position = position_dodge2(width = 0.5)) +  
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", cex = 0.6, alpha = 0.8) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 7, 9, 10, 11, 12), limits = c(0, 12.5)) + 
  labs(title = "(b) Adults without HIV", x = "Carriage incidence ratio (95% CI)", y = "") + 
  theme_bw(base_size = 14, base_family = 'American Typewriter') +
  guides(color = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.85), axis.text.y = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  scale_color_manual(values=c("#F8766D", "#00BFC4"))

#====================================================================

#overview of hazard ratios for pneumococcal acquisition rates
H <- hazard.msm(spn_modelfit, hazard.scale = 1, cl = 0.95)

#create an empty acquisition table for forest plot
spn_acqratio <- tibble(
  index = c(1:10),
  carr = c(rep(c("VT", "NVT"), 5)),
  cova = c(rep("agegp", 4), rep("sex", 2), rep("nochild", 2), rep("ses", 2)),
  label = c(rep("Age group \n (18-24y) vs 25-34y", 2), rep("Age group \n (18-24y) vs 35-45y", 2), rep("Sex \n (female) vs male", 2), rep("Number of children \n (one) vs two+ ", 2), rep("Social economic status \n (high) vs low", 2)),
  HR = c(rep(0, 10)),
  HRlci = c(rep(0, 10)),
  HRuci = c(rep(0, 10))
)

#insert HR from acquisition values in the fitted model
l = 1; m = 1; n = 1
for(i in c("agegp25-34", "agegp35-45", "sexMale", "nochild2+child", "seslow")){
  
  for(j in 1:2){
    spn_acqratio[l, 5] = H[[i]][j]
    l = l+1
  }
  
  for(j in 5:6){
    spn_acqratio[m, 6] = H[[i]][j]
    m = m+1
  }
  
  for(j in 9:10){
    spn_acqratio[n, 7] = H[[i]][j]
    n = n+1
  }
  
}

#plot the acquisition hazard ratios
C <- 
  spn_acqratio %>%
  ggplot(aes(y = reorder(label, index), x = HR, color = carr)) +
  geom_errorbarh(aes(xmin = HRlci, xmax = HRuci), height = 0.5, size = 1.2, position = position_dodge2(width = 0.5)) +
  geom_point(shape = 5, size = 2,  stroke = 2, position = position_dodge2(width = 0.5)) +  
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", cex = 0.6, alpha = 0.8) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 7, 9, 10, 11, 12), limits = c(0, 12.5)) + 
  labs(title = "(c) All adults", x = "", y = "") + 
  theme_bw(base_size = 14, base_family = 'American Typewriter') +
  guides(color = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.85), axis.text.y = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  scale_color_manual(values=c("#A9A9A9", "#36454F"))


#====================================================================
#====================================================================

#overview of hazard ratios for pneumococcal clearance rates
H <- hazard.msm(spn_modelfitpos, hazard.scale = 1, cl = 0.95)

#create an empty clearance table for forest plot
spn_cleratio <- tibble(
  index = c(1:12),
  carr = c(rep(c("VT", "NVT"), 6)),
  cova = c(rep("agegp", 4), rep("sex", 2), rep("dens", 2), rep("abx", 2), rep("artdur", 2)),
  label = c(rep("Age group \n (18-24y) vs 25-34y", 2), rep("Age group \n (18-24y) vs 35-45y", 2), rep("Sex \n (female) vs male", 2), rep("Carriage density \n (high) vs Low", 2), rep("Antibiotic use \n (no) vs yes", 2), rep("ART duration \n (long) vs short", 2)),
  HR = c(rep(0, 12)),
  HRlci = c(rep(0, 12)),
  HRuci = c(rep(0, 12))
)

#insert HR from clearance values in the fitted model
l = 1; m = 1; n = 1
for(i in c("agegp25-34", "agegp35-45", "sexMale", "denslow", "abxtaken", "artdurshort")){
  
  for(j in 3:4){
    spn_cleratio[l, 5] = H[[i]][j]
    l = l+1
  }
  
  for(j in 7:8){
    spn_cleratio[m, 6] = H[[i]][j]
    m = m+1
  }
  
  for(j in 11:12){
    spn_cleratio[n, 7] = H[[i]][j]
    n = n+1
  }
  
}

#plot the clearance hazard ratios
D <- 
  spn_cleratio %>%
  #mutate(HRuci = if_else(HRuci >20, 12, HRuci)) %>% #
  ggplot(aes(y = reorder(label, index), x = HR, color = carr)) +
  geom_errorbarh(aes(xmin = HRlci, xmax = HRuci), height = 0.5, size = 1.2, position = position_dodge2(width = 0.5)) +
  geom_point(shape = 5, size = 2,  stroke = 2, position = position_dodge2(width = 0.5)) +  
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", cex = 0.6, alpha = 0.8) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 7, 9, 10, 11, 12, 40, 44), limits = c(0,45)) + 
  coord_cartesian(xlim = c(0, 12.5)) +
  labs(title = "(d) ALWHIV", x = "", y = "") + 
  theme_bw(base_size = 14, base_family = 'American Typewriter') +
  guides(color = guide_legend(title = "")) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  scale_color_manual(values=c("#7CAE00", "#C77CFF"))


#====================================================================

#overview of hazard ratios for pneumococcal clearance rates
H <- hazard.msm(spn_modelfitneg, hazard.scale = 1, cl = 0.95)

#create an empty clearance table for forest plot
spn_cleratio <- tibble(
  index = c(1:12),
  carr = c(rep(c("VT", "NVT"), 6)),
  cova = c(rep("agegp", 4), rep("sex", 2), rep("dens", 2), rep("abx", 2), rep("artdur", 2)),
  label = c(rep("Age group \n (18-24y) vs 25-34y", 2), rep("Age group \n (18-24y) vs 35-45y", 2), rep("Sex \n (female) vs male", 2), rep("Carriage density \n (high) vs Low", 2), rep("Antibiotic use \n (no) vs yes", 2), rep("ART duration \n (long) vs short", 2)),
  HR = c(rep(0, 12)),
  HRlci = c(rep(0, 12)),
  HRuci = c(rep(0, 12))
)

#insert HR from clearance values in the fitted model
l = 1; m = 1; n = 1
for(i in c("agegp25-34", "agegp35-45", "sexMale", "denslow", "abxtaken")){
  
  for(j in 3:4){
    spn_cleratio[l, 5] = H[[i]][j]
    l = l+1
  }
  
  for(j in 7:8){
    spn_cleratio[m, 6] = H[[i]][j]
    m = m+1
  }
  
  for(j in 11:12){
    spn_cleratio[n, 7] = H[[i]][j]
    n = n+1
  }
  
}

#plot the clearance hazard ratios
E <- 
  spn_cleratio %>%
  mutate(HR = if_else(HR == 0, NA_real_, HR),
         HRlci = if_else(HR == 0, NA_real_, HRlci),
         HRuci = if_else(HR == 0, NA_real_, HRuci)) %>%
  ggplot(aes(y = reorder(label, index), x = HR, color = carr)) +
  geom_errorbarh(aes(xmin = HRlci, xmax = HRuci), height = 0.5, size = 1.2, position = position_dodge2(width = 0.5)) +
  geom_point(shape = 5, size = 2,  stroke = 2, position = position_dodge2(width = 0.5)) +  
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", cex = 0.6, alpha = 0.8) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 7, 9, 10, 11, 12), limits = c(0, 12.5)) + 
  labs(title = "(e) Adults without HIV", x = "Carriage clearance ratio (95% CI)", y = "") + 
  theme_bw(base_size = 14, base_family = 'American Typewriter') +
  guides(color = guide_legend(title = "")) +
  theme(legend.position = "none", axis.text.y = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  scale_color_manual(values=c("#F8766D", "#00BFC4"))

#====================================================================

#overview of hazard ratios for pneumococcal clearance rates
H <- hazard.msm(spn_modelfit, hazard.scale = 1, cl = 0.95)

#create an empty clearance table for forest plot
spn_cleratio <- tibble(
  index = c(1:12),
  carr = c(rep(c("VT", "NVT"), 6)),
  cova = c(rep("agegp", 4), rep("sex", 2), rep("dens", 2), rep("abx", 2), rep("artdur", 2)),
  label = c(rep("Age group \n (18-24y) vs 25-34y", 2), rep("Age group \n (18-24y) vs 35-45y", 2), rep("Sex \n (female) vs male", 2), rep("Carriage density \n (high) vs Low", 2), rep("Antibiotic use \n (no) vs yes", 2), rep("ART duration \n (long) vs short", 2)),
  HR = c(rep(0, 12)),
  HRlci = c(rep(0, 12)),
  HRuci = c(rep(0, 12))
)

#insert HR from clearance values in the fitted model
l = 1; m = 1; n = 1
for(i in c("agegp25-34", "agegp35-45", "sexMale", "denslow", "abxtaken")){
  
  for(j in 3:4){
    spn_cleratio[l, 5] = H[[i]][j]
    l = l+1
  }
  
  for(j in 7:8){
    spn_cleratio[m, 6] = H[[i]][j]
    m = m+1
  }
  
  for(j in 11:12){
    spn_cleratio[n, 7] = H[[i]][j]
    n = n+1
  }
  
}

#plot the clearance hazard ratios
F <- 
  spn_cleratio %>%
  mutate(HR = if_else(HR == 0, NA_real_, HR),
         HRlci = if_else(HR == 0, NA_real_, HRlci),
         HRuci = if_else(HR == 0, NA_real_, HRuci)) %>%
  ggplot(aes(y = reorder(label, index), x = HR, color = carr)) +
  geom_errorbarh(aes(xmin = HRlci, xmax = HRuci), height = 0.5, size = 1.2, position = position_dodge2(width = 0.5)) +
  geom_point(shape = 5, size = 2,  stroke = 2, position = position_dodge2(width = 0.5)) +  
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", cex = 0.6, alpha = 0.8) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 7, 9, 10, 11, 12), limits = c(0, 12.5)) + 
  labs(title = "(f) All adults", x = "", y = "") + 
  theme_bw(base_size = 14, base_family = 'American Typewriter') +
  guides(color = guide_legend(title = "")) +
  theme(legend.position = "none", axis.text.y = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  + 
  scale_color_manual(values=c("#A9A9A9", "#36454F"))


#====================================================================

ggsave(here("output", "Fig3_hazardratios.png"),
       plot = ((A | B | C) / (D | E | F)),
       width = 18, height = 10, unit="in", dpi = 300)


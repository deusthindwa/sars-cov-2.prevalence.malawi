#written by Deus
#01/08/2022
#pneumococcal carriage ad serotype dynamics in HIV-infected adults in the infant PCV era

#====================================================================

# #run multiple chains to assess convergence of the Markov model
# j = 0.001; k = 0.038; l = 0.01; m = 0.48
# 
# for(i in 1:5){
# 
#   sink("/Users/lsh1703394/Rproject/Pneumodude/data/spn_converge.txt", append = TRUE)
# 
#   spn_converge <- msm(state ~ dys, subject = pid, data = spn_model,
#                   qmatrix = rbind(c(0.0, j, l), c(k, 0.0, 0.0), c(m, 0.0, 0.0)),
#                   covariates = list("1-2" = ~ hiv + agegp + sex + nochild + ses, "2-1" =~ hiv + agegp + sex + dens + abx,
#                                     "1-3" = ~ hiv + agegp + sex + nochild + ses, "3-1" =~ hiv + agegp + sex + dens + abx),
#                   pci = c(60, 120, 180, 240),
#                   control = list(fnscale = 1000, maxit = 100000, trace = 1, REPORT = 1))
#   
#   sink()
# 
#   j = j + 0.015; k = k - 0.009; l = l + 0.03; m = m - 0.008
# 
# }

#====================================================================

#retain the convergence dataset for plotting with added variables
spn_convplot <- filter(import(here("data", "spn_converge.xlsx")), chain <=5)

spn_convplot <- 
  spn_convplot %>%
  select(chain,iter,likelihood) %>%
  mutate(chaincat = if_else(chain == 1, "1 (q12=0.001, q21=0.038) \n (q13=0.010, q31=0.048)",
                            if_else(chain == 2,"2 (q12=0.016, q21=0.029) \n (q13=0.040, q31=0.040)",
                                    if_else(chain == 3,"3 (q12=0.031, q21=0.022) \n (q13=0.070, q31=0.032)",
                                            if_else(chain == 4,"4 (q12=0.046, q21=0.011) \n (q13=0.100, q31=0.024)", "5 (q12=0.061, q21=0.002) \n (q13=0.130, q31=0.016)")))),
         likelihood = likelihood*1000)

#====================================================================


spn_convplot %>%
ggplot(aes(x = as.integer(likelihood), y = iter, color = chain), position = position_dodge(width = 0.8)) +
  geom_line(size = 10, position = position_dodge(width = 1.8)) + 
  coord_polar(clip = "off") +
  theme_bw(base_size = 12, base_family = "Lato", base_line_size = 1) +
  theme(axis.text.x = element_text(angle = 40, vjust = 0.5, hjust = 0.3)) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  labs(title = "(D) Model steady states", x = "", y = "") + 
  theme(legend.position = "right") +
  guides(color = guide_legend(title = "Markov chain"), strip.background = element_rect(fill = "light blue")) +
  theme(aspect.ratio = 1)

#====================================================================

D <-
spn_convplot %>%
  ggplot(aes(iter, likelihood, color = chaincat)) + 
  geom_line(size = 1) + 
  labs(title = "(D) Model steady states", x = "Number of iterations", y = "-2Log-likelihood") + 
  coord_cartesian(xlim = c(0, 153), ylim = c(1860, 1950)) +
  theme_bw(base_size = 14, base_family = 'Lato') +
  theme(axis.text.x = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(legend.position = c(0.6,0.7), legend.key.height = unit(1, "line"),legend.key.width = unit(1, "line"), legend.text = element_text(size = 10)) +  
  guides(color = guide_legend(title = "Chain # (initial intensities)")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#====================================================================
#obtain observed versus predicted prevalence of fitted model
spn_obsexp <- 
  tk_tbl(prevalence.msm(spn_modelfit, times = seq(0,338,30)), preserve_index = TRUE, rename_index = "Time", ci = "normal") %>%
  rename(c("obs.nS" = "Observed.State.1", 
           "obs.nVT" = "Observed.State.2", 
           "obs.nNVT" = "Observed.State.3", 
           "obs.tot" = "Observed.Total",
           "obs.pS" = "Observed.percentages.State.1",
           "obs.pVT" = "Observed.percentages.State.2", 
           "obs.pNVT" = "Observed.percentages.State.3",
           "exp.nS" = "Expected.Clear",
           "exp.nVT" = "Expected.VT_carry", 
           "exp.nNVT" = "Expected.NVT_carry", 
           "exp.tot" = "Expected.Total",
           "exp.pS" = "Expected.percentages.Clear", 
           "exp.pVT" = "Expected.percentages.VT_carry",
           "exp.pNVT" = "Expected.percentages.NVT_carry"))

#obtain 95% confidence intervals of the OBSERVED prevalence
spn_obsexp %<>% nest(data = c(obs.nS, obs.tot)) %>% mutate(CI = map(.x = data, ~exactci(.x$obs.nS, .x$obs.tot, conf.level = 0.95)) %>%
           map('conf.int') %>% map(~data.frame(obs.pS_lci = .x[1]*100, obs.pS_uci = .x[2]*100))) %>% unnest_wider(CI) %>% unnest_wider(data)

spn_obsexp %<>% nest(data = c(obs.nVT, obs.tot)) %>% mutate(CI = map(.x = data, ~exactci(.x$obs.nVT, .x$obs.tot, conf.level = 0.95)) %>%
           map('conf.int') %>% map(~data.frame(obs.pVT_lci = .x[1]*100, obs.pVT_uci = .x[2]*100))) %>% unnest_wider(CI) %>% unnest_wider(data)

spn_obsexp %<>% nest(data = c(obs.nNVT, obs.tot)) %>% mutate(CI = map(.x = data, ~exactci(.x$obs.nNVT, .x$obs.tot, conf.level = 0.95)) %>%
           map('conf.int') %>% map(~data.frame(obs.pNVT_lci = .x[1]*100, obs.pNVT_uci = .x[2]*100))) %>% unnest_wider(CI) %>% unnest_wider(data)

#obtain 95% confidence intervals of the FITTED prevalence
spn_obsexp %<>% nest(data = c(exp.nS, exp.tot)) %>% mutate(CI = map(.x = data, ~exactci(.x$exp.nS, .x$exp.tot, conf.level = 0.95)) %>%
           map('conf.int') %>% map(~data.frame(exp.pS_lci = .x[1]*100, exp.pS_uci = .x[2]*100))) %>% unnest_wider(CI) %>% unnest_wider(data)

spn_obsexp %<>% nest(data = c(exp.nVT, exp.tot)) %>% mutate(CI = map(.x = data, ~exactci(.x$exp.nVT, .x$exp.tot, conf.level = 0.95)) %>%
           map('conf.int') %>% map(~data.frame(exp.pVT_lci = .x[1]*100, exp.pVT_uci = .x[2]*100))) %>% unnest_wider(CI) %>% unnest_wider(data)

spn_obsexp %<>% nest(data = c(exp.nNVT, exp.tot)) %>% mutate(CI = map(.x = data, ~exactci(.x$exp.nNVT, .x$exp.tot, conf.level = 0.95)) %>%
           map('conf.int') %>% map(~data.frame(exp.pNVT_lci = .x[1]*100, exp.pNVT_uci = .x[2]*100))) %>% unnest_wider(CI) %>% unnest_wider(data)

spn_obsexp <- 
  spn_obsexp %>%
  select(Time, obs.pS, obs.pS_lci, obs.pS_uci, obs.pVT, obs.pVT_lci, obs.pVT_uci, obs.pNVT, obs.pNVT_lci, obs.pNVT_uci,
         exp.pS, exp.pS_lci, exp.pS_uci, exp.pVT, exp.pVT_lci, exp.pVT_uci, exp.pNVT, exp.pNVT_lci, exp.pNVT_uci) %>%
  slice(-1)

spn_obsexp <- spn_obsexp*0.01
spn_obsexp$Time <- spn_obsexp$Time*100
#====================================================================

#plot observed and predicted carriage
cols <- c("Observed vs predicted clearance" = "#0000FF", "Observed vs predicted carriage" = "#FF0000")

E <-
spn_obsexp %>%
  ggplot(aes(Time)) + 
  geom_point(aes(Time+3, obs.pS, color = "Observed vs predicted clearance"), size = 2, shape = 5, stroke = 2) + 
  geom_errorbar(aes(Time+3, ymin = obs.pS_lci, ymax = obs.pS_uci, color = "Observed vs predicted clearance"), width = 0, size = 1) + 
  geom_line(aes(Time, exp.pS, color = "Observed vs predicted clearance"), size = 1) + 
  geom_ribbon(aes(ymin = exp.pS_lci, ymax = exp.pS_uci, color = "Observed vs predicted clearance"), alpha = 0.2, size = 0.1) +
  geom_point(aes(Time, obs.pVT, color = "Observed vs predicted carriage"), size = 2, shape = 5, stroke = 2, position = position_dodge2(width = 0.1)) + 
  geom_errorbar(aes(Time, ymin = obs.pVT_lci, ymax = obs.pVT_uci, color = "Observed vs predicted carriage"), width = 0, size = 1) + 
  geom_line(aes(Time, exp.pVT, color = "Observed vs predicted carriage"), size = 1) + 
  geom_ribbon(aes(ymin = exp.pVT_lci, ymax = exp.pVT_uci, color = "Observed vs predicted carriage"), alpha = 0.2, size = 0.1) +
  labs(title = "(E) VT carriage fit", x = "Days",y = "Carriage prevalence") + 
  scale_x_continuous(breaks=c(0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330)) + 
  scale_y_continuous(limit = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) + 
  theme_bw(base_size = 14, base_family = 'Lato') +
  theme(axis.text.x=element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(legend.position = c(0.4, 0.85), legend.text = element_text(size = 12),) + 
  guides(color = guide_legend(title = "")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

F <-
spn_obsexp %>%
  ggplot(aes(Time)) + 
  geom_point(aes(Time+3, obs.pS, color = "Observed vs predicted clearance"), size = 2, shape = 5, stroke = 2) + 
  geom_errorbar(aes(Time+3, ymin = obs.pS_lci, ymax = obs.pS_uci, color = "Observed vs predicted clearance"), width = 0, size = 1) + 
  geom_line(aes(Time, exp.pS, color = "Observed vs predicted clearance"), size = 1) + 
  geom_ribbon(aes(ymin = exp.pS_lci, ymax = exp.pS_uci, color = "Observed vs predicted clearance"), alpha = 0.2, size = 0.1) +
  geom_point(aes(Time, obs.pNVT, color = "Observed vs predicted carriage"), size = 2, shape = 5, stroke = 2) + 
  geom_errorbar(aes(Time, ymin = obs.pNVT_lci, ymax = obs.pNVT_uci, color = "Observed vs predicted carriage"), width = 0, size = 1) + 
  geom_line(aes(Time, exp.pNVT, color = "Observed vs predicted carriage"), size = 1) + 
  geom_ribbon(aes(ymin = exp.pNVT_lci, ymax = exp.pNVT_uci, color = "Observed vs predicted carriage"), alpha = 0.2, size = 0.1) +
  labs(title = "(F) NVT carriage fit", x = "Days", y = "Carriage prevalence") + 
  scale_x_continuous(breaks=c(0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330)) + 
  scale_y_continuous(limit = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) + 
  theme_bw(base_size = 14, base_family = 'Lato') +
  theme(axis.text.x=element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(legend.position = c(0.4, 0.85), legend.text = element_text(size = 12)) + 
  guides(color = guide_legend(title = "")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#====================================================================

ggsave(here("reference", "Fig1a_convergencefit.png"),
       plot = (D | E | F),
       width = 18, height = 5, unit="in", dpi = 300)

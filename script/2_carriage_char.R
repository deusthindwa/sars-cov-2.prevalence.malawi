#written by Deus
#01/08/2022
#pneumococcal carriage ad serotype dynamics in HIV-infected adults in the infant PCV era

#====================================================================

A <- 
  spn_fup %>%
  ungroup() %>%
  group_by(vday, hiv, serogroup) %>%
  tally() %>%
  mutate(prev = n/sum(n), N = sum(n)) %>%
  filter(vday != 0, vday != 16) %>%
  rowwise() %>% 
  mutate(obs_lci = exactci(n, N, 0.95) [["conf.int"]][[1]], obs_uci = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%
  filter(serogroup != "None") %>%
  mutate(hivst = if_else(hiv=="HIV-" & serogroup == "VT", "VT, HIV-",
                         if_else(hiv=="HIV-" & serogroup == "NVT", "NVT, HIV-",
                                 if_else(hiv=="HIV+ART+" & serogroup == "VT", "VT, HIV+ART", "NVT, HIV+ART")))) %>%
  
  ggplot() + 
  geom_point(aes(x = vday, y = prev, color = hivst, size = n), shape = 1, stroke = 2) +
  geom_line(aes(x = vday, y = prev, color = hivst), size = 1.5) + 
  #geom_ribbon(aes(x = vday, y = prev, group = hivst, fill = hivst, color = hivst, ymin = obs_lci, ymax = obs_uci), alpha = 0.2, size = 0.1) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_y_continuous(limit = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) + 
  scale_x_continuous(limit = c(1, 15), breaks = seq(1, 15, 2)) + 
  labs(title = "(a)", x = "Visit number", y = "Pneumococcal carriage prevalence") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(legend.text=element_text(size = 12), legend.title = element_text(size = 12)) +
  guides(fill = "none", color = guide_legend(title = "Serotype, HIV status"), size = guide_legend(title = "Sample size")) +
  theme(legend.position = "right") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) 


X <-
  spn_fup %>%
  ungroup() %>%
  group_by(hiv, serogroup) %>%
  tally() %>%
  mutate(prev = n/sum(n), N = sum(n)) %>%
  rowwise() %>% 
  mutate(obs_lci = exactci(n, N, 0.95) [["conf.int"]][[1]], obs_uci = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%
  #filter(serogroup != "None") %>%
  mutate(hivst = if_else(hiv=="HIV-" & serogroup == "VT", "VT, HIV-",
                         if_else(hiv=="HIV-" & serogroup == "NVT", "NVT, HIV-",
                                 if_else(hiv=="HIV+ART+" & serogroup == "VT", "VT, HIV+ART", "NVT, HIV+ART")))) %>%
  
  ggplot(aes(x = prev, y = hivst, color = hivst)) +  
  geom_bar(stat = "identity", position = "dodge", fill = "white", width = 0.8, size = 1) +
  geom_text(aes(label = n), position = position_dodge(0.9), size = 5, hjust = -1.8) +
  geom_errorbar(aes(xmin = obs_lci, xmax = obs_uci), width = 0.2, position = position_dodge(0.9), size = 0.8) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_x_continuous(limit = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) + 
  labs(title="", x = "Total carriage prevalence", y = "") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) 


B <- 
  spn_fup %>%
  ungroup() %>%
  mutate(nochildx = if_else(nochild == 1, "1 child", "2+ children")) %>%
  group_by(hiv, nochildx, serogroup) %>%
  tally() %>%
  mutate(prev = n/sum(n), N = sum(n)) %>%
  rowwise() %>% 
  mutate(obs_lci = exactci(n, N, 0.95) [["conf.int"]][[1]], obs_uci = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%
  filter(serogroup != "None") %>%
  mutate(hivst = if_else(hiv=="HIV-" & serogroup == "VT", "VT, HIV-",
                         if_else(hiv=="HIV-" & serogroup == "NVT", "NVT, HIV-",
                                 if_else(hiv=="HIV+ART+" & serogroup == "VT", "VT, HIV+ART", "NVT, HIV+ART")))) %>%
  
  ggplot(aes(x = hivst, y = prev, color = hivst)) +  
  geom_bar(stat = "identity", position = position_dodge(0.9), fill = "white", size = 1.5) +
  geom_text(aes(label = n), position = position_dodge(0.9), size = 5, vjust = -4.5) +
  geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci), width = 0.2, position = position_dodge(0.9), size = 0.8) +
  facet_grid(.~nochildx) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_y_continuous(limit = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) + 
  labs(title = "(b)", x = "Number of children in the household", y = "") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), strip.background = element_rect(fill = "white"))


C <- 
  spn_fup %>%
  ungroup() %>%
  group_by(hiv, sex, serogroup) %>%
  tally() %>%
  mutate(prev = n/sum(n), N = sum(n)) %>%
  rowwise() %>% 
  mutate(obs_lci = exactci(n, N, 0.95) [["conf.int"]][[1]], obs_uci = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%
  filter(serogroup != "None") %>%
  mutate(hivst = if_else(hiv=="HIV-" & serogroup == "VT", "VT, HIV-",
                         if_else(hiv=="HIV-" & serogroup == "NVT", "NVT, HIV-",
                                 if_else(hiv=="HIV+ART+" & serogroup == "VT", "VT, HIV+ART", "NVT, HIV+ART")))) %>%
  
  ggplot(aes(x = hivst, y = prev, color = hivst)) +  
  geom_bar(stat = "identity", position = position_dodge(0.9), fill = "white", size = 1.5) +
  geom_text(aes(label = n), position = position_dodge(0.9), size = 5, vjust = -4.5) +
  geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci), width = 0.2, position = position_dodge(0.9), size = 0.8) +
  facet_grid(.~sex) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_y_continuous(limit = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) + 
  labs(title = "(c)", x = "Sex", y = "") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), strip.background = element_rect(fill = "white"))

D <- 
  spn_fup %>%
  ungroup() %>%
  filter(!is.na(ses)) %>%
  mutate(ses = as.factor(if_else(ses <= 3, "Low",
                                 if_else(ses > 3 & ses <=15, "High", NA_character_)))) %>%
  group_by(hiv, ses, serogroup) %>%
  tally() %>%
  mutate(prev = n/sum(n), N = sum(n)) %>%
  rowwise() %>% 
  mutate(obs_lci = exactci(n, N, 0.95) [["conf.int"]][[1]], obs_uci = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%
  filter(serogroup != "None") %>%
  mutate(hivst = if_else(hiv=="HIV-" & serogroup == "VT", "VT, HIV-",
                         if_else(hiv=="HIV-" & serogroup == "NVT", "NVT, HIV-",
                                 if_else(hiv=="HIV+ART+" & serogroup == "VT", "VT, HIV+ART", "NVT, HIV+ART")))) %>%
  
  ggplot(aes(x = hivst, y = prev, color = hivst)) +  
  geom_bar(stat = "identity", position = position_dodge(0.9), fill = "white", size = 1.5) +
  geom_text(aes(label = n), position = position_dodge(0.9), size = 5, vjust = -3.5) +
  geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci), width = 0.2, position = position_dodge(0.9), size = 0.8) +
  facet_grid(.~ses) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_y_continuous(limit = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) + 
  labs(title = "(d)", x = "Socioeconomic status", y = "Pneumococcal carriage prevalence") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), strip.background = element_rect(fill = "white"))

E <- 
  spn_fup %>%
  ungroup() %>%
  filter(!is.na(age), serogroup != "None") %>%
  mutate(hivst = if_else(hiv=="HIV-" & serogroup == "VT", "VT\nHIV-",
                         if_else(hiv=="HIV-" & serogroup == "NVT", "NVT\nHIV-",
                                 if_else(hiv=="HIV+ART+" & serogroup == "VT", "VT\nHIV+ART", "NVT\nHIV+ART")))) %>%
  
  ggplot(aes(y = age, x = hivst, color = hivst)) + 
  geom_boxplot(notch = TRUE, size = 1.6) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_y_continuous(limit = c(18, 45), breaks = seq(18, 45, 5)) + 
  labs(title = "(e)", x="", y = "Age (years)") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11)) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))


F <- 
  spn_fup %>%
  ungroup() %>%
  filter(!is.na(dens), serogroup != "None") %>%
  mutate(hivst = if_else(hiv=="HIV-" & serogroup == "VT", "VT\nHIV-",
                         if_else(hiv=="HIV-" & serogroup == "NVT", "NVT\nHIV-",
                                 if_else(hiv=="HIV+ART+" & serogroup == "VT", "VT\nHIV+ART", "NVT\nHIV+ART")))) %>%
  
  ggplot(aes(y = dens, x = hivst, color = hivst)) + 
  geom_boxplot(notch = TRUE, size = 1.6) +
  scale_y_continuous(trans = log10_trans()) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(f)", x="", y = "Pneumococcal carriage density (logCFU/ml)") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11)) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))


G <- 
  spn_fup %>%
  ungroup() %>%
  filter(!is.na(artdur), serogroup != "None") %>%
  group_by(serogroup, artdur)  %>%
  
  ggplot() + 
  geom_density(aes(x = artdur, fill = serogroup), color = "black", alpha = 0.5, size = 1.5) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(g)", x = "ART duration (Years) among HIV+ART", y = "Probability density") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text=element_text(size = 11), legend.title = element_text(size = 11)) + 
  guides(fill=guide_legend(title="Carriage")) +
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  scale_fill_manual("", values = c("NVT" = "#7CAE00", "VT" = "#C77CFF"))


Y <- 
  spn_fup %>%
  ungroup() %>%
  filter(!is.na(artdur), serogroup != "None") %>%
  
  ggplot(aes(y = artdur, x = serogroup, fill = serogroup)) + 
  geom_boxplot(notch = TRUE, color = "black", size = 1) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_y_continuous(trans = log10_trans()) +
  labs(title = "", x="", y = "ART duration (years)") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 12, angle = 90), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11)) + 
  theme(legend.position = "none") + 
  scale_fill_manual("", values = c("NVT" = "#7CAE00", "VT" = "#C77CFF"))


H <- 
  spn_fup %>%
  ungroup() %>%
  filter(!is.na(cd4), serogroup != "None", hiv == "HIV+ART+") %>%
  group_by(hiv, serogroup, cd4)  %>%
  
  ggplot() + 
  geom_density(aes(x = cd4, fill = serogroup), color = "black", alpha = 0.5, size = 1.5) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(h)", x = "CD4 count among HIV+ART", y = "Probability density") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(legend.text=element_text(size = 11), legend.title = element_text(size = 11)) + 
  guides(fill=guide_legend(title="Carriage")) +
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  scale_fill_manual("", values = c("NVT" = "#7CAE00", "VT" = "#C77CFF"))

#combined plots
ggsave(here("output", "Fig1_carriage_char.png"),
       plot = ((A | inset_element(X, right = 0.9, left = 0.3, bottom = 0.66, top = 0.99) | B | C | plot_layout(ncol = 3, width = c(2,1,1)))) / 
         (( D | E | F | G | inset_element(Y, right = 0.95, left = 0.50, bottom = 0.4, top = 0.99)) | plot_layout(ncol = 4, width = c(2,2,2,2,2))),
       width = 21, height = 14, unit="in", dpi = 300)

#delete individual plots after saving above
rm(A, B, C, D, E, F, G, H, X, Y)

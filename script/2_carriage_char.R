#written by Deus
#01/08/2022
#pneumococcal carriage ad serotype dynamics in HIV-infected adults in the infant PCV era

#====================================================================

#BASELINE SAMPLE DESCRIPTION

#total HIV
spn_baseline %>%
  group_by(hiv) %>%
  tally() %>%
  mutate(p = n/sum(n))

#hiv by serotype group
spn_baseline %>%
  group_by(hiv, serogroup) %>%
  tally() %>%
  mutate(p = n/sum(n))

#hiv by sex group
spn_baseline %>%
  group_by(hiv, sex) %>%
  tally() %>%
  mutate(p = n/sum(n))

#hiv by agegp group
spn_baseline %>%
  filter(!is.na(age)) %>%
  group_by(hiv) %>%
  summarise(mqage = quantile(age, 0.50),
            fqage = quantile(age, 0.25),
            tqage = quantile(age, 0.75))

spn_baseline %>%
  group_by(hiv, agegp) %>%
  tally() %>%
  mutate(p = n/sum(n))

#hiv by number of children in the household
spn_baseline %>%
  group_by(hiv, nochild) %>%
  tally() %>%
  mutate(p = n/sum(n))

#hiv by pneumococcal density
spn_baseline %>%
  group_by(hiv) %>%
  filter(!is.na(dens)) %>%
  summarise(mqdens = quantile(dens, 0.50),
            fqdens = quantile(dens, 0.25),
            tqdens = quantile(dens, 0.75))

spn_baseline %>%
  group_by(hiv, dens2) %>%
  tally() %>%
  mutate(p = n/sum(n))

#====================================================================

#FOLLOW UP SAMPLE DESCRIPTION

A <- 
  spn_model %>%
  ungroup() %>%
  group_by(vday, hivst) %>%
  tally() %>%
  mutate(prev = n/sum(n), N = sum(n)) %>%
  ungroup() %>%
  rowwise() %>% 
  mutate(obs_lci = exactci(n, N, 0.95) [["conf.int"]][[1]], obs_uci = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%
  ungroup() %>%
  filter(!is.na(hivst)) %>%
  
  ggplot() + 
  geom_point(aes(x = vday, y = prev, color = hivst, size = n), shape = 1, stroke = 2) +
  geom_line(aes(x = vday, y = prev, color = hivst), size = 1.5) + 
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_y_continuous(limit = c(0, 0.3), breaks = seq(0, 0.3, 0.05), labels = scales::percent_format(accuracy = 1)) + 
  scale_x_continuous(limit = c(1, 17), breaks = seq(1, 17, 2)) + 
  labs(title = "(a)", x = "Visit number", y = "Pneumococcal carriage prevalence") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(legend.text=element_text(size = 12), legend.title = element_text(size = 12)) +
  guides(fill = "none", color = guide_legend(title = "Serotype, HIV status"), size = guide_legend(title = "Sample size")) +
  theme(legend.position = "right") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) 


X <-
  spn_model %>%
  ungroup() %>%
  group_by(hivst) %>%
  tally() %>%
  mutate(prev = n/sum(n), N = sum(n)) %>%
  ungroup() %>%
  rowwise() %>% 
  mutate(obs_lci = exactci(n, N, 0.95) [["conf.int"]][[1]], obs_uci = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%
  ungroup() %>%
  filter(!is.na(hivst)) %>%
  
  ggplot(aes(x = prev, y = hivst, color = hivst)) +  
  geom_bar(stat = "identity", position = "dodge", fill = "white", width = 0.8, size = 1) +
  geom_text(aes(label = n), position = position_dodge(0.9), size = 5, hjust = -1.8) +
  geom_errorbar(aes(xmin = obs_lci, xmax = obs_uci), width = 0.2, position = position_dodge(0.9), size = 0.8) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_x_continuous(limit = c(0, 0.2), breaks = seq(0, 0.2, 0.05), labels = scales::percent_format(accuracy = 1)) + 
  labs(title="", x = "Total carriage prevalence", y = "") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) 


B <- 
  spn_model %>%
  filter(!is.na(hivst)) %>%
  
  mutate(hivst0 = if_else(hivst == "VT,HIV-", "VT\nHIV-",
                          if_else(hivst == "NVT,HIV-", "NVT\nHIV-",
                                  if_else(hivst == "VT,ART<3m", "VT\nART<3m",
                                          if_else(hivst == "NVT,ART<3m", "NVT\nART<3m",
                                                  if_else(hivst == "VT,ART>1y", "VT\nART>1y",
                                                          if_else(hivst == "NVT,ART>1y", "NVT\nART>1y", NA_character_))))))) %>%
  
  ggplot(aes(y = dens0, x = hivst0, color = hivst0)) + 
  geom_boxplot(notch = TRUE, size = 1.6) +
  scale_y_continuous(trans = log10_trans()) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(b)", x="", y = "Pneumococcal carriage density (log_CFU/ml)") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 11)) + 
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))


C <- 
  spn_model %>%
  filter(!is.na(hivst)) %>%
  group_by(nochild, hivst) %>%
  tally() %>%
  mutate(prev = n/sum(n), N = sum(n)) %>%
  ungroup() %>%
  
  ggplot(mapping = aes(x = nochild, y = prev, color = hivst, fill = hivst)) + 
  geom_bar(stat = "identity", color = "black", size = 0.7) +
  geom_text(aes(label = n, fontface = 2), size = 5, color = "black", position = position_stack(vjust = 0.5)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(c)", x = "Number of children in the house", y = "Share of total carriage samples") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), strip.background = element_rect(fill = "white"))

  
D <- 
  spn_model %>%
  filter(!is.na(hivst)) %>%
  group_by(sex, hivst) %>%
  tally() %>%
  mutate(prev = n/sum(n), N = sum(n)) %>%
  ungroup() %>%
  
  ggplot(mapping = aes(x = sex, y = prev, color = hivst, fill = hivst)) + 
  geom_bar(stat = "identity", color = "black", size = 0.7) +
  geom_text(aes(label = n, fontface = 2), size = 5, color = "black", position = position_stack(vjust = 0.5)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(d)", x = "Sex", y = "") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 0)) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), strip.background = element_rect(fill = "white"))

  
E <- 
  spn_model %>%
  filter(!is.na(hivst)) %>%
  group_by(agegp, hivst) %>%
  tally() %>%
  mutate(prev = n/sum(n), N = sum(n)) %>%
  ungroup() %>%
  
  ggplot(mapping = aes(x = agegp, y = prev, color = hivst, fill = hivst)) + 
  geom_bar(stat = "identity", color = "black", size = 0.7) +
  geom_text(aes(label = n, fontface = 2), size = 5, color = "black", position = position_stack(vjust = 0.5)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(e)", x = "Age group", y = "") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 0)) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), strip.background = element_rect(fill = "white"))


#load shape file of Malawi map and subset for Blantyre
# nasotmp <- tempfile()
# download.file("https://raw.githubusercontent.com/deusthindwa/dlnm.typhoid.nts.climate.blantyre.malawi/master/data/malawi_map.zip", destfile = nasotmp)
# unzip(nasotmp, exdir = ".")
# mw.map <- rgdal::readOGR(".","malawi_map")
# bt1.map <- mw.map@data$OBJECTID >289 & mw.map@data$OBJECTID <297 #id from 290 to 296 
# bt2.map <- mw.map@data$OBJECTID >308 & mw.map@data$OBJECTID <311 #id from 309 to 310
# bt3.map <- mw.map@data$OBJECTID >342  #id fom 243
# bt.map <- bind_rows(fortify(mw.map[bt1.map,]), fortify(mw.map[bt2.map,]), fortify(mw.map[bt3.map,]))
# readr::write_csv(x = bt.map, file = here("data", "spn_btmap.csv"))

F <-
  spn_btmap %>%
  ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "gray50", size = 1) + 
  geom_point(data = spn_loc, aes(x = lon, y = lat, size = nsample), color = "black", shape = 1, stroke = 3) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(f)", x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(face = "bold", size = 14), axis.text.y = element_text(face = "bold", size = 14)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), strip.background = element_rect(fill = "white")) + 
  guides(color = FALSE, fill = FALSE, size = guide_legend(title = "Sample size")) +
  theme(legend.position = c(0.15,0.85))

  
#combined plots
ggsave(here("output", "Fig1_carriage_char.png"),
       plot = ((A | inset_element(X, right = 0.9, left = 0.3, bottom = 0.46, top = 0.99) | B | plot_layout(ncol = 2, width = c(2,2)))) / 
         (( C | D | E | F ) | plot_layout(ncol = 4, width = c(1.1,1.1,1.1,3))),
       width = 21, height = 14, unit="in", dpi = 300)

#delete individual plots after saving above
rm(A, B, C, D, E, F, X)

#written by Deus
#01/08/2022
#pneumococcal carriage ad serotype dynamics in HIV-infected adults in the infant PCV era

#====================================================================

# serotype prevalence
A <- 
  spn_model %>% 
  filter(!is.na(serotype)) %>%
  group_by(serotype) %>%
  tally() %>%
  mutate(sp = n/sum(n)) %>%
  ggplot(aes(x = reorder(serotype, sp)), fill = serotype) + 
  geom_bar(aes(y = sp), color = "black", stat = 'identity', size = 1, width = 0.6) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(a)", x = "Occurence of detected serotype during all visits combined", y = "Share of total serotypes") +
  scale_fill_grey(start = 0.9, end = 0.2) +
  scale_y_continuous(limit = c(0, 0.55), breaks = seq(0, 0.55, 0.1), labels = scales::percent_format(accuracy = 1)) + 
  theme(axis.text.x = element_text(face = "bold", size = 10, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(face = "bold", size = 10)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 10), axis.title.y=element_text(face = "bold", size = 10)) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm"))

X <- 
spn_model %>%
  mutate(stnew = if_else(sstate == 1L, NA_character_,
                         if_else(sstate == 2L, "11A/B/C/D/F",
                                 if_else(sstate == 3L, "7A/B/C",
                                         if_else(sstate == 4L, "3",
                                                 if_else(sstate == 5L, "19F",
                                                         if_else(sstate == 6L, "17A/F",
                                                                 if_else(sstate == 7L, "15A/B/C/F",
                                                                         if_else(sstate == 8L, "6C/D",
                                                                                 if_else(sstate == 9L, "23A/B",
                                                                                         if_else(sstate == 10L, "10A/B/C/F",
                                                                                                 if_else(sstate == 11L, "6A/B",
                                                                                                         if_else(sstate == 12L, "oVT",
                                                                                                                 if_else(sstate == 13L, "uNVT", "kNVT")))))))))))))) %>%
  
  filter(!is.na(stnew)) %>%
  group_by(stnew) %>%
  tally() %>%
  mutate(N = sum(n),
         sp = n/N) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(lsp = exactci(n, N, 0.95) [["conf.int"]][[1]],
         usp = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%

  ggplot(aes(x = reorder(stnew, sp))) + 
  geom_bar(aes(y = sp, fill = stnew, color = stnew), stat = 'identity', size = 1, width = 0.5) + 
  geom_errorbar(aes(reorder(stnew, sp), sp, ymin = lsp, ymax = usp), color = "black", width = 0, size = 1.2) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "", x = "", y = "Share of total serotypes") +
  scale_y_continuous(limit = c(0, 0.60), breaks = seq(0, 0.55, 0.10), labels = scales::percent_format(accuracy = 1)) + 
  theme(axis.text.x = element_text(face = "bold", size = 10, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(face = "bold",size = 10)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(face = "bold", size = 10), axis.title.y=element_text(face = "bold", size = 10)) +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) 

# circle network graph data creation
spn_model <- arrange(spn_model, pid, dys)
spn_wiw <- 
  as.data.frame(statetable.msm(sstate, pid, data = spn_model)) %>% 
  mutate(Freq = Freq+1) %>%
  dplyr::rename("weight" = Freq) %>%
  filter(weight>1 & from != to & to != 1 & from !=1) %>%
  mutate(from = if_else(from == 2L, "11A/B/C/D/F",
                        if_else(from == 3L, "7A/B/C",
                                if_else(from == 4L, "3",
                                        if_else(from == 5L, "19F",
                                                if_else(from == 6L, "17A/F",
                                                        if_else(from == 7L, "15A/B/C/F",
                                                                if_else(from == 8L, "6C/D",
                                                                        if_else(from == 9L, "23A/B",
                                                                                if_else(from == 10L, "10A/B/C/F",
                                                                                        if_else(from == 11L, "6A/B",
                                                                                                if_else(from == 12L, "oVT",
                                                                                                        if_else(from == 13L, "uNVT", "kNVT")))))))))))),
         to = if_else(to == 2L, "11A/B/C/D/F",
                      if_else(to == 3L, "7A/B/C",
                              if_else(to == 4L, "3",
                                      if_else(to == 5L, "19F",
                                              if_else(to == 6L, "17A/F",
                                                      if_else(to == 7L, "15A/B/C/F",
                                                              if_else(to == 8L, "6C/D",
                                                                      if_else(to == 9L, "23A/B",
                                                                              if_else(to == 10L, "10A/B/C/F",
                                                                                      if_else(to == 11L, "6A/B",
                                                                                              if_else(to == 12L, "oVT",
                                                                                                      if_else(to == 13L, "uNVT", "kNVT"))))))))))))) %>%
  graph.data.frame(directed = FALSE)

# circle network graph plotting
B <- as.ggplot(
~plot(spn_wiw, 
     edge.width = log(E(spn_wiw)$weight)*3.5,
     #edge.label = (E(spn_wiw)$weight),
     edge.arrow.size = 0.8,
     edge.curved = 0.3,
     edge.color = "lightblue", 
     vertex.color = "white", 
     vertex.frame.color = "black", 
     vertex.size = 40,
     vertex.label.color = "darkred",
     vertex.label.font = 2,                          
     vertex.label.cex = 0.8,
     layout = layout.circle,
     cex.axis = 1,
     main = "",
     ), 
scale = 1.15
)

B <- B + 
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(b)", x = "Network graph of pneumococcal serotype acquisition chains") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#====================================================================

#serotype specific acquisition probability
modela <- pmatrix.msm(spn_modelfit3, t = 1, ci = "normal", cl = 0.95)

pneumo0 <- data.frame("serotype" = c("11A/B/C/D/F", "7A/B/C", "3", "19F", "17A/F", "15A/B/C/F", "6C/D", "23A/B", "10A/B/C/F", "6A/B", "oVT", "uNVT", "kNVT"))
pneumo0$carry <- c(modela$estimates[1,2], modela$estimates[1,3], modela$estimates[1,4], modela$estimates[1,5], modela$estimates[1,6], modela$estimates[1,7], modela$estimates[1,8], modela$estimates[1,9], modela$estimates[1,10], modela$estimates[1,11], modela$estimates[1,12], modela$estimates[1,13], modela$estimates[1,14])
pneumo0$Lcarry <- c(modela$L[1,2], modela$L[1,3], modela$L[1,4], modela$L[1,5], modela$L[1,6], modela$L[1,7], modela$L[1,8], modela$L[1,9], modela$L[1,10], modela$L[1,11], modela$L[1,12], modela$L[1,13], modela$L[1,14]) 
pneumo0$Ucarry <- c(modela$U[1,2], modela$U[1,3], modela$U[1,4], modela$U[1,5], modela$U[1,6], modela$U[1,7], modela$U[1,8], modela$U[1,9], modela$U[1,10], modela$U[1,11], modela$U[1,12], modela$U[1,13], modela$U[1,14]) 

C <- 
  pneumo0 %>%
  mutate(acq = if_else(serotype == "uNVT", "high", "low")) %>%
  ggplot() +
  geom_point(aes(reorder(serotype, carry), carry, color = serotype), size = 1, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(serotype,  ymin = Lcarry, ymax = Ucarry, color = serotype), width = 0, size = 1, position = position_dodge2(width = 0.5)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(c)", x = "Pneumococcal serotype", y = "Daily carriage acquisition probability") + 
  facet_grid(acq ~., scales = "free_y" ) +
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(shape = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = "none", legend.title = element_text(size = 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(strip.background = element_blank(),strip.text.y = element_blank())

#====================================================================

#serotype specific carriage duration (inverse clearance rate)
modelb <- qmatrix.msm(spn_modelfit3, ci = "normal", cl = 0.95)

pneumo1 <- data.frame("serotype" = c("11A/B/C/D/F", "7A/B/C", "3", "19F", "17A/F", "15A/B/C/F", "6C/D", "23A/B", "10A/B/C/F", "6A/B", "oVT", "uNVT", "kNVT"))
pneumo1$clear <- c(modelb$estimates[2,1], modelb$estimates[3,1], modelb$estimates[4,1], modelb$estimates[5,1], modelb$estimates[6,1], modelb$estimates[7,1], modelb$estimates[8,1], modelb$estimates[9,1], modelb$estimates[10,1], modelb$estimates[11,1], modelb$estimates[12,1], modelb$estimates[13,1], modelb$estimates[14,1])
pneumo1$Lclear <- c(modelb$L[2,1], modelb$L[3,1], modelb$L[4,1], modelb$L[5,1], modelb$L[6,1], modelb$L[7,1], modelb$L[8,1], modelb$L[9,1], modelb$L[10,1], modelb$L[11,1], modelb$L[12,1], modelb$L[13,1], modelb$L[14,1])
pneumo1$Uclear <- c(modelb$U[2,1], modelb$U[3,1], modelb$U[4,1], modelb$U[5,1], modelb$U[6,1], modelb$U[7,1], modelb$U[8,1], modelb$U[9,1], modelb$U[10,1], modelb$U[11,1], modelb$U[12,1], modelb$U[13,1], modelb$U[14,1])

D <- 
  pneumo1 %>%
  mutate(clearx = 1/clear, Lclearx = 1/Uclear, Uclearx = 1/Lclear) %>%
  ggplot() +
  geom_point(aes(reorder(serotype, clearx), clearx, color = serotype), size = 2, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_errorbar(aes(serotype,  ymin = Lclearx, ymax = Uclearx, color = serotype), width = 0, size = 1, position = position_dodge2(width = 0.5)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_y_continuous(breaks=c(10, 20, 30, 60, 75, 90), limits = c(0, 95)) + 
  labs(title = "(d)", x = "Pneumococcal serotype", y = "Average carriage duration (days)") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  guides(shape = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 10), legend.position = "none", legend.title = element_text(size = 10)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#====================================================================

#combined plots
ggsave(here::here("output", "Fig4_st_dynamics.png"),
       plot = ((A | inset_element(X, right = 0.8, left = 0.2, bottom = 0.12, top = 0.98))/(B | C | D | plot_layout(ncol = 3, width = c(2,1,1)))),
       width = 15, height = 10, unit = "in", dpi = 300)


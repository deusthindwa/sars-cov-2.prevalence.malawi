#written by Deus
#13/02/2024
#pneumococcal carriage and serotype dynamics by adult HIV status a mature PCV program

#====================================================================
X <-
spn_base %>%
  dplyr::group_by(serotype) %>%
  dplyr::tally() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(serotype = if_else(serotype == "10A/B/C/F", "10A/10B/10C/10F", serotype),
                serotype = if_else(serotype == "11A/B/C/D/F", "11A/11B/11C/11D/11F", serotype),
                serotype = if_else(serotype == "15A/B/C/F", "15A/15B/15C/15F", serotype),
                serotype = if_else(serotype == "17A/F", "17A/17F", serotype),
                serotype = if_else(serotype == "18A/B/C/F", "18A/18B/18C/18F", serotype),
                serotype = if_else(serotype == "19B/C", "19B/19C", serotype),
                serotype = if_else(serotype == "22A/F", "22A/22F", serotype),
                serotype = if_else(serotype == "23A/B", "23A/23B", serotype),
                serotype = if_else(serotype == "7A/B/C", "7A/7B/7C", serotype),
                serotype = if_else(serotype == "9A/L/N", "9A/9L/9N", serotype)) %>%
  dplyr::mutate(n = if_else(is.na(serotype), n, if_else(str_length(serotype)>3, n/(lengths(gregexpr("/", serotype)) + 1), n))) %>%
  tidyr::separate_rows(., serotype) %>% #tidy serotypes to atomic values
  dplyr::group_by(serotype) %>%
  dplyr::tally(n) %>%
  dplyr::ungroup() %>%
  
  dplyr::mutate(pcv7pfz = if_else(grepl("\\b(4|6B|9V|14|18C|19F|23F)\\b", serotype) == TRUE, "PCV7", if_else(is.na(serotype), NA_character_, "NVT")),
                pcv10sii = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", serotype) == TRUE, "PCV10-sii", if_else(is.na(serotype), NA_character_, "NVT")),
                pcv10gsk = if_else(grepl("\\b(1|4|5|6B|7F|9V|14|18C|19F|23F)\\b", serotype) == TRUE, "PCV10-gsk", if_else(is.na(serotype), NA_character_, "NVT")),
                pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", serotype) == TRUE, "PCV13", if_else(is.na(serotype), NA_character_, "NVT")),
                pcv15mek = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", serotype) == TRUE, "PCV15", if_else(is.na(serotype), NA_character_, "NVT")),
                pcv20pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", serotype) == TRUE, "PCV20", if_else(is.na(serotype), NA_character_, "NVT")))

X %>%
  dplyr::group_by(pcv7pfz) %>%
  dplyr::tally(n) %>%
  dplyr::mutate(p = n/sum(n))

X %>%
  dplyr::group_by(pcv10sii) %>%
  dplyr::tally(n) %>%
  dplyr::mutate(p = n/sum(n))

X %>%
  dplyr::group_by(pcv10gsk) %>%
  dplyr::tally(n) %>%
  dplyr::mutate(p = n/sum(n))

X %>%
  dplyr::group_by(pcv13pfz) %>%
  dplyr::tally(n) %>%
  dplyr::mutate(p = n/sum(n))

X %>%
  dplyr::group_by(pcv15mek) %>%
  dplyr::tally(n) %>%
  dplyr::mutate(p = n/sum(n))

X %>%
  dplyr::group_by(pcv20pfz) %>%
  dplyr::tally(n) %>%
  dplyr::mutate(p = n/sum(n))


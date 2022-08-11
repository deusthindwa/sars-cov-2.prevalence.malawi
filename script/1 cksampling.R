#Written by Deus Thindwa
#Estimating SARS-CoV-2 prevalence in Malawi
#Generalized additive model.
#01/8/2022 - 31/12/2022

#=======================================================================

#load PHIM dataset
ckcens <- import(here("data", "chikwawacensus.csv")) %>%
  
#=======================================================================

#data wrangling
select(data_date, cluster, ta, ea_code, village_name, household_id, household_name, child_id, child_name) %>% 
  mutate(cluster = as_factor(cluster), ta = as_factor(ta), village_name = as_factor(village_name)) 

#=======================================================================

# perform stratified sampling
ckcens_samp <- ckcens %>% filter(ta !="NULL" & child_id != "C7055446") #cleaning
set.seed(1988) #reproducibility
ckcens_samp <- stratified(ckcens_samp, c("cluster", "ta", "ea_code"), 4) #sampling
ckcens_samp <- arrange(ckcens_samp, cluster, ta, ea_code, village_name, household_id)
write.csv(ckcens_samp, here("output", "ckcensus_sampled.csv"))

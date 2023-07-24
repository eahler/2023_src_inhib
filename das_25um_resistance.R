# Classify DAS 25um variants as "resistant" based on 2x Std Dev. Compare with previous activity classifications
# from 2019 Ahler et. al. Relates to Figure S2A.

library(tidyverse)

das_25um_variants <- read_tsv("data/DAS_25um_Activity.tsv")
das_25um_syn <- read_tsv("data/DAS_25um_synonymous_scores.tsv")
activity_scores <- read_tsv("data/2019_mol_cell_activity_annotated.tsv") %>%
  select(variant, Classification) #This data is from a previous publication (2019 Ahler et al)


#Calculate the standard deviation and mean of the synonymous variants
syn_std_dev <- sd(das_25um_syn$Activity, na.rm = TRUE) #0.3274449
syn_mean <- mean(das_25um_syn$Activity, na.rm = TRUE) #0.09431948

#Determine the activity score boundary for calling a variant as "resistant"
resistant_activity_threshold <- syn_mean + (2*syn_std_dev) #0.74920935

# Define "resistant" variants in DAS 25um
resistant_df <- das_25um_variants %>%
  filter(Activity >= resistant_activity_threshold)

#Final dataframe
final_df <- resistant_df %>%
  left_join(activity_scores, by = "variant") %>%
  drop_na() 

# Calculate percent of each "Classification" group for DAS 25um resistance variants
final_df %>%
  group_by(Classification) %>%
  dplyr::summarise(percentage = n() / nrow(final_df) * 100 )



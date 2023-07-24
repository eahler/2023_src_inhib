# Generating base plots for Figure 1D/1E

library(dplyr)
library(readr)
library(ggplot2)

##### DAS 25um
das_25um_syn_df <- read_tsv("data/DAS_25um_synonymous_scores.tsv") %>% 
  mutate(label='synonymous')

das_25um_nonsyn_df <- read_tsv("data/DAS_25um_Activity.tsv") %>% 
  mutate(label = ifelse(grepl('Ter',variant)==TRUE, 'nonsense','missense')) #label missense and nonsense

das_25um_merge <- rbind(das_25um_syn_df,das_25um_nonsyn_df) %>% na.omit() #remove any variants found in one dasatinib concentration but not other


das_25um_plot <- ggplot(das_25um_merge, aes(Activity, fill = label, colour = label)) +
  geom_density(alpha = 0.1)+ 
  theme( axis.title = element_text(color="black",size=20), axis.text = element_text(color="black",size=15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) 

das_25um_plot

ggsave("DAS_25um_score_distribution.pdf", das_25um_plot,  height = 6, width =6, units = "in", dpi=300)

##### DAS 100um
das_100um_syn_df <- read_tsv("data/DAS_100um_synonymous_scores.tsv") %>% 
  mutate(label='synonymous')

das_100um_nonsyn_df <- read_tsv("data/DAS_100um_Activity.tsv") %>% 
  mutate(label = ifelse(grepl('Ter',variant)==TRUE, 'nonsense','missense')) #label missense and nonsense

das_100um_merge <- rbind(das_100um_syn_df,das_100um_nonsyn_df) %>% na.omit() #remove any variants found in one dasatinib concentration but not other


das_100um_plot <- ggplot(das_100um_merge, aes(Activity, fill = label, colour = label)) +
  geom_density(alpha = 0.1)+ 
  theme( axis.title = element_text(color="black",size=20), axis.text = element_text(color="black",size=15), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) 

das_100um_plot

ggsave("DAS_100um_score_distribution.pdf", das_100um_plot,  height = 6, width =6, units = "in", dpi=300)

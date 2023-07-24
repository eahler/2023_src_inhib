# To quantile normalize scores 
# Inhibitors 1, 2, and 4 from manuscript are denoted as 1, 15, and 2 in these data, respectively
# Python 2.7

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Read in scores for Inhibitors 1, 2, and 4
syn_df = pd.read_csv('Data/A1-Inhibitor_Syn', sep='\t')
nonsyn_df = pd.read_csv('Data/A1-Inhibitor_Nonsyn', sep='\t')

casted_nonsyn = nonsyn_df.pivot(index='variant',columns='variable')['score']
casted_syn = syn_df.pivot(index='variant',columns='variable')['score']

########### Nonsynonymous scores
rank_mean_nonsyn = casted_nonsyn.stack().groupby(casted_nonsyn.rank(method='first').stack().astype(int)).mean()
quant_normalized_nonsyn = casted_nonsyn.rank(method='min').stack().astype(int).map(rank_mean_nonsyn).unstack()

# Check nonsyn normalization 
quant_normalized_nonsyn.columns = [' '.join(col).strip() for col in quant_normalized_nonsyn.columns.values] #Collapse the column multiindex into a concatenation
quant_normalized_nonsyn.head()

# Remove WT and Synonymous score
quant_normalized_nonsyn = quant_normalized_nonsyn.drop(labels=['_sy','_wt'])

# Fix column names and remove index
quant_normalized_nonsyn.rename(columns = {list(quant_normalized_nonsyn)[0]: '15_High Score'},inplace = True) 
quant_normalized_nonsyn.rename(columns = {list(quant_normalized_nonsyn)[1]: '1_High Score'},inplace = True) 
quant_normalized_nonsyn.rename(columns = {list(quant_normalized_nonsyn)[2]: '2_High Score'},inplace = True) 
quant_normalized_nonsyn['variant'] = quant_normalized_nonsyn.index
quant_normalized_nonsyn = quant_normalized_nonsyn.reset_index(drop=True)


# Inverse scores to "Activity" to make more intuitive sense, save each inhibitor as separate data frame 
inhibitor_1 = quant_normalized_nonsyn[['1_High Score','variant']]
inhibitor_1['Activity'] = -(inhibitor_1['1_High Score'])
inhibitor_1 = inhibitor_1[['variant','Activity']]
inhibitor_1.to_csv("Data/Inhibitor_1_Activity_Scores", sep='\t',index=False)

inhibitor_2 = quant_normalized_nonsyn[['15_High Score','variant']]
inhibitor_2['Activity'] = -(inhibitor_2['15_High Score'])
inhibitor_2 = inhibitor_2[['variant','Activity']]
inhibitor_2.to_csv("Data/Inhibitor_2_Activity_Scores", sep='\t',index=False)

inhibitor_4 = quant_normalized_nonsyn[['2_High Score','variant']]
inhibitor_4['Activity'] = -(inhibitor_4['2_High Score'])
inhibitor_4 = inhibitor_4[['variant','Activity']]
inhibitor_4.to_csv("Data/Inhibitor_4_Activity_Scores", sep='\t',index=False)


########### Synonymous scores
rank_mean_syn = casted_syn.stack().groupby(casted_syn.rank(method='first').stack().astype(int)).mean()
quant_normalized_syn = casted_syn.rank(method='min').stack().astype(int).map(rank_mean_syn).unstack()

# Check synonymous normalization 
quant_normalized_syn.columns = [' '.join(col).strip() for col in quant_normalized_syn.columns.values] #Collapse the column multiindex into a concatenation
quant_normalized_syn.head()

# Fix column names and remove index
quant_normalized_syn.rename(columns = {list(quant_normalized_syn)[0]: '15_High Score'},inplace = True) 
quant_normalized_syn.rename(columns = {list(quant_normalized_syn)[1]: '1_High Score'},inplace = True) 
quant_normalized_syn.rename(columns = {list(quant_normalized_syn)[2]: '2_High Score'},inplace = True) 
quant_normalized_syn['variant'] = quant_normalized_syn.index
quant_normalized_syn = quant_normalized_syn.reset_index(drop=True)

# Inverse scores to "Activity" to make more intuitive sense, save each inhibitor as separate data frame 
inhibitor_1_syn = quant_normalized_syn[['1_High Score','variant']]
inhibitor_1_syn['Activity'] = -(inhibitor_1_syn['1_High Score'])
inhibitor_1_syn = inhibitor_1_syn[['variant','Activity']]
inhibitor_1_syn.to_csv("Data/Inhibitor_1_Activity_Scores_Synonymous", sep='\t',index=False)

inhibitor_2_syn = quant_normalized_syn[['15_High Score','variant']]
inhibitor_2_syn['Activity'] = -(inhibitor_2_syn['15_High Score'])
inhibitor_2_syn = inhibitor_2_syn[['variant','Activity']]
inhibitor_2_syn.to_csv("Data/Inhibitor_2_Activity_Scores_Synonymous", sep='\t',index=False)

inhibitor_4_syn = quant_normalized_syn[['2_High Score','variant']]
inhibitor_4_syn['Activity'] = -(inhibitor_4_syn['2_High Score'])
inhibitor_4_syn = inhibitor_4_syn[['variant','Activity']]
inhibitor_4_syn.to_csv("Data/Inhibitor_4_Activity_Scores_Synonymous", sep='\t',index=False)


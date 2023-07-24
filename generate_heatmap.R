# Hacky script to generate heatmap displaying activity scores for each mutant for a given condition in desired visual format 
# Generates separate positive, negative, and missing heatmaps, merge them together in Inkscape or Illustrator afterwards and add residue/amino acid labels manually
# Example data is for DAS 25um (Figure 2A), change the input score dataframe to generate maps for Inhibitors 1/2/4 (Figure S3B)


library(dplyr)
library(plyr)
library(readr)
library(ggplot2)
library(tidyr)
library(stringr)
library(gridExtra)
library(RColorBrewer)
Classified_Variants <-read_tsv("data/DAS_25um_Activity.tsv") #For Inhibitors 1/2/4, change this to respecgtive dataset
wt_pattern1 <- "p\\.[:upper:][:lower:][:lower:]"
Classified_Variants <- mutate(Classified_Variants, wt1 = str_extract(variant, wt_pattern1)) #Grabs first ".pXXX"
wt_pattern2 <- "[:upper:][:lower:][:lower:]"
Classified_Variants <- mutate(Classified_Variants, wt = str_extract(wt1, wt_pattern2)) #Grabs only the XXX (aka the WT)
Classified_Variants[3] <- NULL
residue_pattern <- "[:digit:]+"
Classified_Variants <- mutate(Classified_Variants, residue = str_extract(variant, residue_pattern)) #Grabs residue 
mutant_pattern <- "[:upper:][:lower:][:lower:]$"
Classified_Variants <- mutate(Classified_Variants, mutant = str_extract(variant, mutant_pattern)) #Grabs mutant residue


#For single AA, need to change all three letter codes to single letter, uses plyr's mapvalues function
Classified_Variants$wt <- mapvalues(Classified_Variants$wt, from=c("His", "Lys", "Arg", "Asp", "Glu", "Cys", "Met", "Asn", "Gln", "Ser", "Thr", "Ala", "Gly", "Ile", "Leu", "Pro", "Val", "Phe", "Trp", "Tyr", "Ter"), to=c("H","K","R","D","E","C","M","N","Q","S","T","A","G","I","L","P","V","F","W","Y","*"))
Classified_Variants$mutant <- mapvalues(Classified_Variants$mutant, from=c("His", "Lys", "Arg", "Asp", "Glu", "Cys", "Met", "Asn", "Gln", "Ser", "Thr", "Ala", "Gly", "Ile", "Leu", "Pro", "Val", "Phe", "Trp", "Tyr", "Ter"), to=c("H","K","R","D","E","C","M","N","Q","S","T","A","G","I","L","P","V","F","W","Y","*"))

WTSEQ <- "LRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGETGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIEDNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMPCPPECPESLHDLMCQCWRKEPEERPTFEYLQAFL"
wtTEST <- sapply(seq(from=1, to=nchar(WTSEQ), by=1), function(i) substr(WTSEQ, i, i))


ClassTemp <- Classified_Variants

ClassTemp$residue <- as.numeric(ClassTemp$residue)
ClassTemp <- ClassTemp[with(ClassTemp, order(residue)), ] #Orders by ascending residue number

#Creates vector of all possible amino acids

Amino_Acids <- c("H","K","R","D","E","C","M","N","Q","S","T","A","G","I","L","P","V","F","W","Y","*")


#Creates vector of all residues
All_Residues <- c(270:520)

mutant <- character()
for (i in 1:250) {
  for (j in 1:21) {
    mutant <- append(mutant, Amino_Acids[j])
  }
}
mutant
residue <- character()
for (i in 1:250) {
  for (j in 1:21) {
    residue <- append(residue, All_Residues[i])
  }
}




Empty_Res <- cbind(residue, mutant)
Empty_Res <- as.data.frame(Empty_Res)
Empty_Res <- tbl_df(Empty_Res)
ClassTemp$residue = as.character(as.numeric(ClassTemp$residue))


Empty_Res$residue = as.character(as.factor(Empty_Res$residue))
Empty_Res$mutant = as.character(as.factor(Empty_Res$mutant))


Empty_Res <- mutate(Empty_Res, wt = wtTEST[as.numeric(Empty_Res$residue)-269]) #Adds the wt AA to each position


ClassTemp$joinby <- paste(ClassTemp$residue, ClassTemp$mutant) #joinby is just used for merging purposes
Empty_Res$joinby <- paste(Empty_Res$residue, Empty_Res$mutant)

Missing_Residues <- anti_join(Empty_Res, ClassTemp, by = "joinby")
Missing_Residues <- mutate(Missing_Residues, Activity = NA)

ClassTemp[1] <- NULL #Gets rid of variant
Merge_Classification <- rbind(ClassTemp, Missing_Residues)
Merge_Classification <- Merge_Classification[with(Merge_Classification, order(residue)), ]
str(Merge_Classification)

Merge_Classification$Activity[Merge_Classification$mutant == Merge_Classification$wt] <- 0 #Changes "NA" to "0" for any WT residu
Missing_Residues <- filter(Missing_Residues, mutant != wt)






ClassTest <- Merge_Classification
ClassTest[5] <- NULL #Get rid of joinby column
ClassTest$residue <- as.numeric(ClassTest$residue)
ClassTest <- ClassTest[with(ClassTest, order(residue)), ]
ClassTest$Activity <- as.numeric(ClassTest$Activity)
str(ClassTest)

#This gets the missing data in the same format as the variant dataset
MissingTest <- Missing_Residues
MissingTest[4] <- NULL #Get rid of joinby column
MissingTest$residue <- as.numeric(MissingTest$residue)
MissingTest <- MissingTest[with(MissingTest, order(residue)), ]

#Generates just the positive activity score heatmap
positive <- filter(ClassTest, Activity >= 0)

#Generates just the negative activity score heatmap
negative <- filter(ClassTest, Activity < 0)


all_frames = ClassTest[ClassTest$wt == ClassTest$mutant, c("mutant", "residue")] #This gets the WT residue for the point

pos_color <- c("#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b")
neg_color <- c("#0571b0", "#92c5de", "#f7f7f7")


#Make positive, negative, and NA for nlobe positions
all_pos_graph <- ggplot(positive, aes(mutant, residue)) + 
  geom_tile(aes(fill = Activity, colour="black"))+ 
  coord_equal() +
  scale_y_reverse() + 
  scale_fill_gradientn(colours=pos_color,limits = c(0,4),na.value="#67001f") +
  scale_color_manual(values="black")+ theme(panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(),
                                            legend.position="none",axis.title=element_blank(),
                                            axis.text=element_blank(),
                                            axis.ticks=element_blank(),panel.background = element_blank(), line = element_blank())                                                                                                                                                                                                                                                            
all_pos_final <- all_pos_graph + 
  geom_point(data=all_frames, aes(x=mutant, y=residue), size =.05) + 
  scale_x_discrete(limits=Amino_Acids) 
all_pos_final
ggsave("DAS_25um_pos.pdf", all_pos_final, path = "heatmaps/", height = 30, width = 5, units = "cm", dpi=300)

all_neg_graph <- ggplot(negative, aes(mutant, residue)) + 
  geom_tile(aes(fill = Activity, color="black"))+ 
  coord_equal() + scale_y_reverse() + 
  scale_fill_gradientn(colours=neg_color) + 
  scale_color_manual(values="black") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none",axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),panel.background = element_blank(), line = element_blank()) + 
  scale_x_discrete(limits=Amino_Acids) 
                                                                                                                                                                                                                                  
all_neg_graph
ggsave("DAS_25um_neg.pdf", all_neg_graph, path = "heatmaps/", height = 30, width = 5, units = "cm", dpi=300)





all_missing_graph <- ggplot(MissingTest, aes(mutant, residue)) + 
  geom_tile(aes(fill = Activity, colour="black"))+ 
  coord_equal() + 
  scale_y_reverse() + 
  scale_fill_manual(values = "white", na.value = "grey70")+ 
  scale_color_manual(values="black") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position="none",axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),panel.background = element_blank(),line = element_blank()) + 
  scale_x_discrete(limits=Amino_Acids) 

all_missing_graph
ggsave("DAS_25um_missing.pdf", all_missing_graph, path = "heatmaps/", height = 30, width = 5, units = "cm", dpi=300)

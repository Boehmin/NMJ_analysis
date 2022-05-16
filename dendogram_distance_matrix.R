############################################################
# This script loads a .csv of your NMJ data from NMJ-morph #
# with additional columns for species, muscle and animal.  #
# Animal is for example Mouse 1 L, Mouse 1 R (L = left,    #
# R = right), or Cat 1 - Cat 9 etc...                      #
# This allows comparison between species, muscles and      #
# animals. Each row represents values for one NMJ.         #
############################################################
#########################
### LOAD DEPENDENCIES ###
#########################
# Ensure these are installed first with install.packages()
library(igraph)
library(tidyverse)
library(usedist)
library(vegan)
library(magrittr) # use pipes etc
library(dendextend)

rm(list = ls())

#########################
###    Main Code      ###
#########################

# 1. Load Spreadsheets
data <- read.csv("~/Filepath/PCA_muscavr.csv", header = TRUE, dec = ",")

# THE FOLLOWING STEP MIGHT NOT BE NECESSARY FOR YOU. CHECK NUMERICAL VALUES IN YOUR DATA FRAME
# 2. Change data frame from character form into numeric form
options(digits=5) # preserve 9 digits after comma, necessary for as.numeric function
char_columns <- sapply(data[4:23], is.character) # Identify character columns
data_chars_as_num <- data[4:23]    # Replicate data
data_chars_as_num[ , char_columns] <- as.data.frame(   # Recode characters as numeric
  apply(data_chars_as_num[ , char_columns], 2, as.numeric))
sapply(data_chars_as_num, mode)                       # Print classes of all colums
df <- data.frame(data[1:3], data_chars_as_num) # Replicate data in new dataframe merged with first few columns still as character columns

# 3. Calculate pairwise distances of the dataset in Euclidean space of Species
mrpp(df[4:23], df$Species, permutations = 999, distance ='euclidean', weight.type = 1)

# 3.1 Plot a Dendrogram using euclidean Distance of all variables across all animals
dm <- vegdist(df[4:23], method='euclidian')
meandm <- meandist(dm, df$Animal)

tiff("~/Filepath/Cluster_Dendrogram_allvariables.tiff", units="cm", width=10, height=15, res=600, compression = 'lzw')

plot(hclust(as.dist(meandm)))
dev.off()

# 3.2 Plot just a Dendrogram of pre-synaptic variables (minus axon diameter) across species
dpre <- vegdist(df[5:23], method='euclidian')
meandpre <- meandist(dpre, df$Species)

tiff("~/Filepath/Cluster_Dendrogram_pre.tiff", units="cm", width=10, height=15, res=600, compression = 'lzw')

plot(hclust(as.dist(meandpre)))
dev.off()

# 3.3 Plot just a Dendrogram of post-synaptic variables (minus fibre diameter) across species
dpost <- vegdist(df[12:22], method='euclidian')
meandpost <- meandist(dpost, df$Species)

tiff("~/Filepath/Cluster_Dendrogram_post.tiff", units="cm", width=10, height=15, res=600, compression = 'lzw')

plot(hclust(as.dist(meandpost)))
dev.off()

# 4. Comparison of pre and post-synapse - TANGLEGRAM
d1 <- meandpost %>% 
  dist() %>% 
  hclust(method="average") %>% 
  as.dendrogram()

d2 <- meandpre %>% 
  dist() %>% 
  hclust(method="complete") %>% 
  as.dendrogram()

# 4.1 Customise these, and place them in a list
# you can choose your colours here for your tanglegram
dl <- dendlist(
  d1 %>% 
    set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3),
  d2 %>% 
    set("labels_col", value = c("skyblue", "orange", "grey"), k=3) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3)
)

# 4.2 Plot your tanglegram
tanglegram(dl, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
           margin_inner=7,
           lwd=2
)
dend_list <- dendlist(d1,d2)

tiff("~/Filepath/tanglegram_post_pre.tiff", units="cm", width=15, height=15, res=600, compression = 'lzw')

tanglegram(d1,d2,
           highlight_distinct_edges = T, # Turn-off dashed lines
           common_subtrees_color_lines = T, # Turn-off line colors
           common_subtrees_color_branches = T, # Color common branches 
)
dev.off()
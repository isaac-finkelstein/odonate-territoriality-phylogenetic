#phylogenetic anaylsis of odonate territoriality

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(dplyr)

my_data<- read.csv("data/species_formatted_master.csv")
load(file="data/Odo.tree.Waller.Svensson.2017.rda") #odonate tree extracated from Waller and Svensson 2017
#str(tree) #plot(tree, no.margin=TRUE) #tree$Nnode #tree$tip.label -- to check the structure of the tree

#Put data into a readable dataframe -- one species per row
terr_table<- ftable(my_data$Formatted_species, my_data$Territorial) #this has 3 rows
prop_terr<- round(terr_table[,3]/(terr_table[,2]+ terr_table[,3]),3) # so this calculates the percentage "yes" for territorial
sn<- attr(terr_table,"row.vars")[[1]]
terr_data_with_na<- data.frame(sn,prop_terr) #so this is a dataframe with the proporitons of territorial "yes" for each species
terr_data<- terr_data_with_na[complete.cases(terr_data_with_na), ] #removed NA values 

#prune tree to match data
chk<-name.check(tree, terr_data, data.names=as.character(terr_data$sn))
summary(chk)
odonate_tree<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree
# Identify species to drop from terr_data
species_to_drop <- chk$data_not_tree
odonate_terr_data <- terr_data[!(terr_data$sn %in% species_to_drop), ] #droopped data_not_tree species from dataset
name.check(odonate_tree, odonate_terr_data, data.names=as.character(odonate_terr_data$sn))

#let's plot the tree
plot(odonate_tree, type="fan", cex=0.5, fsize=0.6, ftype="i") #it's very busy because it has 384 tips

#just anisoptera -- node 386 
#collapseTree(odonate_tree)  -- use this to find the node numbers

anisoptera_tree_extract<-unname(ape::extract.clade(odonate_tree, node = 386)$tip.label)
anisoptera_tree<-drop.tip(odonate_tree, anisoptera_tree_extract)
plot(anisoptera_tree, type="fan", cex=0.5, ftype="i")

#just zygoptera -- node 584
zygoptera_tree_extract<-unname(ape::extract.clade(odonate_tree, node = 584)$tip.label)
zygoptera_tree<-drop.tip(odonate_tree, zygoptera_tree_extract)
plot(zygoptera_tree, type="fan", cex=0.5, ftype="i")

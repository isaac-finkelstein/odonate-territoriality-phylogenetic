#clean anisoptera only analysis

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(dplyr)
library(corHMM)

my_data<- read.csv("data/data_vfinal.csv")

load(file="data/Odo.tree.Waller.Svensson.2017.rda") #odonate tree extracated from Waller and Svensson 2017
#str(tree) #plot(tree, no.margin=TRUE) #tree$Nnode #tree$tip.label #check the structure of the tree

#tree from Rocha-Ortega et al., 2020
#same structure as Waller and Svensson 2017 but has different species coverage.
#tree<-read.nexus(file="data/Rocha_ortega_tree.nex")

#my dataset has data that uses the same original source
#remove instances where data is recorded more than once and uses the same original source
rows_to_remove <- c(2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 27, 30, 32, 33, 34, 36, 37, 38, 41, 42, 43, 44, 46, 47, 48, 49, 50, 52, 53, 55, 56, 57, 58, 59, 62, 64, 66, 67, 69, 71, 72, 73, 75, 76, 77, 78, 79, 289, 317, 322, 324, 326, 345, 356, 360, 371, 374, 390, 393, 488, 1068, 1089, 8639, 8658, 8659, 8661, 8662)
my_data <- my_data[-rows_to_remove, ]

#Put data into a readable dataframe - one species per row
terr_table<- ftable(my_data$Formatted_species, my_data$Territorial) #3 columns, the first column is for NAs
prop_terr<- round(terr_table[,3]/(terr_table[,2]+ terr_table[,3]),2) #calculates the percentage "yes" for territorial
#set a 3:1 threshhold = 75% threshold
sp_terr <- ifelse(prop_terr >= 0.75, 1, ifelse(prop_terr <= 0.25, 0, NA))

sn<- attr(terr_table,"row.vars")[[1]]
terr_data_with_na<- data.frame(sn,sp_terr) #dataframe with territorial (1/0) for each species
terr_data<- terr_data_with_na[complete.cases(terr_data_with_na), ] #remove NA values 

#prune tree to match data
chk<-name.check(tree, terr_data, data.names=as.character(terr_data$sn))
summary(chk)
odonate_tree<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree
# Identify species to drop from terr_data
species_to_drop <- chk$data_not_tree
odonate_terr_data <- terr_data[!(terr_data$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_terr_data, data.names=as.character(odonate_terr_data$sn))

terr_mode<-setNames(odonate_terr_data$sp_terr, odonate_terr_data$sn)

#just anisoptera -- node 568 
#plot(odonate_tree, show.node.label = TRUE)  # Plot tree with tip labels
#nodelabels()  # Adds node numbers to the plot


anisoptera_tree_extract<-unname(ape::extract.clade(odonate_tree, node = 568)$tip.label)
anisoptera_tree<-drop.tip(odonate_tree, anisoptera_tree_extract)
plot(anisoptera_tree, type="fan", cex=0.5, ftype="i") 
dev.off()


#prune tree to match data
chk_anis<-name.check(anisoptera_tree, odonate_terr_data, data.names=as.character(odonate_terr_data$sn))
summary(chk_anis)
anis_tree<-drop.tip(anisoptera_tree, chk_anis$tree_not_data) #dropped tree_not_data species
anis_tree
# Identify species to drop from terr_data - using the factor "Yes or No"
species_to_drop_anis <- chk_anis$data_not_tree
anis_terr_data <- odonate_terr_data[!(odonate_terr_data$sn %in% species_to_drop_anis), ] #dropped data_not_tree species from dataset
name.check(anis_tree, anis_terr_data, data.names=as.character(anis_terr_data$sn))

anis_terr_mode <- setNames(anis_terr_data$sp_terr, anis_terr_data$sn)


#marginal ancestral state reconstruction
fit_marginal_anis<- corHMM(anis_tree, anis_terr_data, node.states = "marginal",
                           rate.cat=1, rate.mat=NULL, model="ARD")
fit_marginal_anis
head(fit_marginal_anis$states)       
#interpret this matrix as the posterior probabilities that each state is in each node

anis_terr_mode <- factor(ifelse(anis_terr_data$sp_terr == 1, "yes", "no"))
names(anis_terr_mode) <- anis_terr_data$sn

#plot this
cols<-setNames(c("turquoise", "brown"), levels(anis_terr_mode)) #brown =yes territorail, turquoise = not territorial
plotTree.datamatrix(anis_tree, as.data.frame(anis_terr_mode),
                    colors=list(cols), header=FALSE)
legend("topright", legend=levels(anis_terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
nodelabels(pie=fit_marginal_anis$states, piecol=cols, cex=0.3)
tiplabels(pie = to.matrix(anis_terr_mode, sort(unique(anis_terr_mode))), piecol=cols, cex=0.3)

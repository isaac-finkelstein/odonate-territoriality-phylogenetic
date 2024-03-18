#Phylo analysis for only Anisoptera
#phylogenetic anaylsis of odonate territoriality

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(dplyr)
library(corHMM)

my_data<- read.csv("data/data_v4.csv") #this dataset (3rd version) switches "tandem" for "contact" in De Recende's data.
#For some reason, they use both terms. But since they mean the same thing, I changed them all to "Contact"
#I also emailed to ask, but did not get a response. 
load(file="data/Odo.tree.Waller.Svensson.2017.rda") #odonate tree extracated from Waller and Svensson 2017
#str(tree) #plot(tree, no.margin=TRUE) #tree$Nnode #tree$tip.label -- to check the structure of the tree

#Put data into a readable dataframe -- one species per row
terr_table<- ftable(my_data$Formatted_species, my_data$Territorial) #this has 3 columns, the first column is for NAs
prop_terr<- round(terr_table[,3]/(terr_table[,2]+ terr_table[,3]),2) # so this calculates the percentage "yes" for territorial
#the 2 at the end rounds to 2 decimal 
#set a 3:1 threshhold = 75% threshold
sp_terr <- ifelse(prop_terr >= 0.75, 1, ifelse(prop_terr <= 0.25, 0, NA))

sn<- attr(terr_table,"row.vars")[[1]]
terr_data_with_na<- data.frame(sn,sp_terr) #so this is a dataframe with territorial (1/0) for each species
terr_data<- terr_data_with_na[complete.cases(terr_data_with_na), ] #removed NA values 

#to make the territorial variable binary (yes or no) make it a factor
prop_terr_factor <- factor(ifelse(sp_terr == 1, "yes", "no"))
terr_data_with_na_factor <- data.frame(sn, prop_terr_factor = prop_terr_factor) # Create a data frame
terr_data_factor <- terr_data_with_na_factor[complete.cases(terr_data_with_na_factor), ] #Removes NA values

#prune tree to match data
chk<-name.check(tree, terr_data, data.names=as.character(terr_data$sn))
summary(chk)
odonate_tree<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree
# Identify species to drop from terr_data
species_to_drop <- chk$data_not_tree
odonate_terr_data <- terr_data[!(terr_data$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_terr_data, data.names=as.character(odonate_terr_data$sn))

#prune data to match tree
odonate_tree_factor<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree_factor
species_to_drop <- chk$data_not_tree
odonate_terr_data_factor <- terr_data_factor[!(terr_data_factor$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_terr_data_factor, data.names=as.character(odonate_terr_data_factor$sn))

terr_mode<-setNames(odonate_terr_data_factor$prop_terr, odonate_terr_data_factor$sn)

#just anisoptera -- node 411 
#EVERYTIME I UPDATE THE DATASET, THE NODE NUMBER WILL CHANGE!! SO HAVE TO UPDATE THESE EVERYTIME!
#collapseTree(odonate_tree)  #-- use this to find the node numbers
anisoptera_tree_extract<-unname(ape::extract.clade(odonate_tree, node = 566)$tip.label)
anisoptera_tree<-drop.tip(odonate_tree, anisoptera_tree_extract)
plot(anisoptera_tree, type="fan", cex=0.5, ftype="i") 


#Plotting territoriality on Anisoptera tree
#prune tree to match data
chk_anis<-name.check(anisoptera_tree, terr_data_factor, data.names=as.character(terr_data_factor$sn))
summary(chk_anis)
anis_tree<-drop.tip(anisoptera_tree, chk_anis$tree_not_data) #dropped tree_not_data species
anis_tree
# Identify species to drop from terr_data - using the factor "Yes or No"
species_to_drop_anis <- chk_anis$data_not_tree
anis_terr_data_factor <- terr_data_factor[!(terr_data_factor$sn %in% species_to_drop_anis), ] #dropped data_not_tree species from dataset
name.check(anis_tree, anis_terr_data_factor, data.names=as.character(anis_terr_data_factor$sn))

#plot anisoptera tree
anis_terr_mode <- setNames(anis_terr_data_factor$prop_terr, anis_terr_data_factor$sn)
cols<-setNames(c("blue", "goldenrod1"), levels(anis_terr_mode))
plotTree.datamatrix(anis_tree, as.data.frame(anis_terr_mode),
                    colors= list(cols), header=FALSE)
legend("topright", legend=levels(anis_terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
#this works
#the solution to make it work was changing the setNames function line
#it's so crowded that I might want to remove names and then provide a large list so people can look up whatever species they want
#Can't figure out how to remove the species labels!

#choose a character model
#4 possible models for discreet response
#equal rates (ER)
#all rates different (ARD)
# irreversible modes (2) - going from 0 -> 1 or from 1 -> 0 but not vice versa)

fit_er_anis<-fitMk(anis_tree, terr_mode, model = "ER")
fit_ard_anis<-fitMk(anis_tree, terr_mode, model = "ARD")
#fit the no -> yes model (transition from non-territorial to territorial but not back)
fit_01_anis<-fitMk(anis_tree, terr_mode, model=matrix(c(0,1,0,0),2,2,byrow=TRUE))
#fit the yes to no model (transition from territorial to non-territorail but not back)
fit_10_anis<-fitMk(anis_tree, terr_mode, model=matrix(c(0,0,1,0),2,2, byrow=TRUE))

#extract AIC values from the models
aic_model_anis<-c(AIC(fit_er_anis), AIC(fit_ard_anis), AIC(fit_01_anis), AIC(fit_10_anis))
data.frame(model=c("ER", "ARD", "no -> yes", "yes->no"),
           logL=c(logLik(fit_er_anis), logLik(fit_ard_anis),
                  logLik(fit_01_anis), logLik(fit_10_anis)),
           AIC=aic_model_anis,delta.AIC=aic_model_anis-min(aic_model_anis))             
#0kay, so it seems that the all-rates-different is the way to go   
#ARD is quite strongly favoured according to AIC values 

#we can estimate q (the transition rate)
#Q is the expected number of transition given a particular amount of time 
print(fit_er_anis)
print(fit_ard_anis)
print(fit_01_anis)
print(fit_10_anis)



#joint ancestral state reconstruction
#character data must be numerical interger (0,1)
anis_data<-data.frame(Genus_sp=names(anis_terr_mode), anis_terr_mode=as.numeric(anis_terr_mode)-1)
head(anis_data, n=10)
#fit the joint model
fit_joint_anis<-corHMM(anis_tree, anis_data, node.states="joint", rate.cat=1, model="ARD",rate.mat=NULL)
fit_joint_anis
#now plot
plotTree.datamatrix(anis_tree, as.data.frame(anis_terr_mode),
                    colors=list(cols), header=FALSE)
legend("topright", legend=levels(anis_terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
nodelabels(pie=to.matrix(levels(anis_terr_mode)[fit_joint_anis$phy$node.label],
                         levels(anis_terr_mode)), piecol=cols, cex=0.3)

#marginal ancestral state reconstruction
fit_marginal_anis<- corHMM(anis_tree, anis_data, node.states = "marginal",
                           rate.cat=1, rate.mat=NULL, model="ARD")
fit_marginal_anis
head(fit_marginal_anis$states)       
#interpret this matrix as the posterior probabilities that each state is in each node

#plot this
cols<-setNames(c("turquoise", "brown"), levels(anis_terr_mode)) #brown =yes territorail, turquoise = not territorial
plotTree.datamatrix(anis_tree, as.data.frame(anis_terr_mode),
                    colors=list(cols), header=FALSE)
legend("topright", legend=levels(anis_terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
nodelabels(pie=fit_marginal_anis$states, piecol=cols, cex=0.3)
tiplabels(pie = to.matrix(anis_terr_mode, sort(unique(anis_terr_mode))), piecol=cols, cex=0.3)
#likely ancestral state = non-territorial! - that's interesting because it differs from odonates

#plot this with a fan shape
cols<-setNames(c("turquoise", "brown"), levels(anis_terr_mode)) #brown =yes territorail, turquoise = not territorial
anis_fit_ARD_again<-ace(anis_terr_mode, anis_tree, model="ARD", type="discrete") #this gives us the marginal ancestral state - also called "empirical Bayesian posterior probabilities"
round(anis_fit_ARD_again$lik.anc, 3)
plotTree(anis_tree, type="fan", fsize=0.6, ftype="i")
nodelabels(node=1:anis_tree$Nnode+Ntip(anis_tree),
           pie=anis_fit_ARD_again$lik.anc, piecol = cols, cex=0.3)
tiplabels(pie=to.matrix(anis_terr_mode, sort(unique(anis_terr_mode))), piecol=cols, cex=0.3)
legend("topright", legend=levels(anis_terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
#why does it look differnet from when I plotted marginal ancestral state reconstruction?
#I think it's because I used Ancestral Character Estimation using the function ace from package Ape
#whereas before I estimate ancestral state using corHMM
#the corHMM uses a markov model
#Idk why these would be different ugghhh
#or maybe it's just because there's less data when I only use anisoptera?
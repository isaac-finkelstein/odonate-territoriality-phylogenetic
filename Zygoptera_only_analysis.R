#Phylo analysis for only Zygoptera
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

#just zygoptera -- node 374
#EVERYTIME I UPDATE THE DATASET, THE NODE NUMBER WILL CHANGE!! SO HAVE TO UPDATE THESE EVERYTIME!
#collapseTree(odonate_tree)  #-- use this to find the node numbers
zygoptera_tree_extract<-unname(ape::extract.clade(odonate_tree, node = 374)$tip.label)
zygoptera_tree<-drop.tip(odonate_tree, zygoptera_tree_extract)
plot(zygoptera_tree, type="fan", cex=0.5, ftype="i")

#prune tree to match data
chk_zygo<-name.check(zygoptera_tree, terr_data_factor, data.names=as.character(terr_data_factor$sn))
summary(chk_zygo)
zygo_tree<-drop.tip(zygoptera_tree, chk_zygo$tree_not_data) #dropped tree_not_data species
zygo_tree
# Identify species to drop from terr_data - using the factor "Yes or No"
species_to_drop_zygo <- chk_zygo$data_not_tree
zygo_terr_data_factor <- terr_data_factor[!(terr_data_factor$sn %in% species_to_drop_zygo), ] #dropped data_not_tree species from dataset
name.check(zygo_tree, zygo_terr_data_factor, data.names=as.character(zygo_terr_data_factor$sn))

#plot zygoptera tree
zygo_terr_mode <- setNames(zygo_terr_data_factor$prop_terr, zygo_terr_data_factor$sn)
cols<-setNames(c("forestgreen", "goldenrod1"), levels(zygo_terr_mode))
plotTree.datamatrix(zygo_tree, as.data.frame(zygo_terr_mode),
                    colours= list(cols), header=FALSE)
legend("topright", legend=levels(zygo_terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
#Can't figure out how to remove the species labels!



#choose a character model
#4 possible models for discreet response
#equal rates (ER)
#all rates different (ARD)
# irreversible modes (2) - going from 0 -> 1 or from 1 -> 0 but not vice versa)
fit_er_zygo<-fitMk(zygo_tree, terr_mode, model = "ER")
fit_ard_zygo<-fitMk(zygo_tree, terr_mode, model = "ARD")
#fit the no -> yes model (transition from non-territorial to territorial but not back)
fit_01_zygo<-fitMk(zygo_tree, terr_mode, model=matrix(c(0,1,0,0),2,2,byrow=TRUE))
#fit the yes to no model (transition from territorial to non-territorail but not back)
fit_10_zygo<-fitMk(zygo_tree, terr_mode, model=matrix(c(0,0,1,0),2,2, byrow=TRUE))

#extract AIC values from the models
aic_model_zygo<-c(AIC(fit_er_zygo), AIC(fit_ard_zygo), AIC(fit_01_zygo), AIC(fit_10_zygo))
data.frame(model=c("ER", "ARD", "no -> yes", "yes->no"),
           logL=c(logLik(fit_er_zygo), logLik(fit_ard_zygo),
                  logLik(fit_01_zygo), logLik(fit_10_zygo)),
           AIC=aic_model_zygo,delta.AIC=aic_model_zygo-min(aic_model_zygo))             
#0kay, so the no->yes model has the lowest AIC, but the ARD model is close.  


#we can estimate q (the transition rate)
#Q is the expected number of transition given a particular amount of time 
print(fit_er_zygo)
print(fit_ard_zygo)
print(fit_01_zygo)
print(fit_10_zygo)


#Next: Ancestral state reconstruction - plots are very busy so may want to re-run these for just anisoptera/zygoptera
#I have to decide if I want to do joint or marginal acenstral state reconstruction
#marginal is more popular in biology (Revell and Harmon, 2022)
#marginal measures uncertainty about the specific values for ancestral states - so preferred
#do both, then compare results (joint can go in supplementary material)
#character data must be numerical interger (0,1)
zygo_data<-data.frame(Genus_sp=names(zygo_terr_mode), zygo_terr_mode=as.numeric(zygo_terr_mode)-1)
head(zygo_data, n=10)
#fit the joint model
fit_joint_zygo<-corHMM(zygo_tree, zygo_data, node.states="joint", rate.cat=1, model="no -> yes",rate.mat=NULL)
fit_joint_zygo
#now plot
plotTree.datamatrix(zygo_tree, as.data.frame(zygo_terr_mode),
                    colors=list(cols), header=FALSE)
legend("topright", legend=levels(zygo_terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
nodelabels(pie=to.matrix(levels(zygo_terr_mode)[fit_joint_zygo$phy$node.label],
                         levels(zygo_terr_mode)), piecol=cols, cex=0.3)

#marginal ancestral state reconstruction
fit_marginal_zygo<- corHMM(zygo_tree, zygo_data, node.states = "marginal",
                           rate.cat=1, rate.mat=NULL)
fit_marginal_zygo
head(fit_marginal_zygo$states)       
#interpret this matrix as the posterior probabilities that each state is in each node

#plot this
plotTree.datamatrix(zygo_tree, as.data.frame(zygo_terr_mode),
                    colors=list(cols), header=FALSE)
legend("topright", legend=levels(zygo_terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
nodelabels(pie=fit_marginal_zygo$states, piecol=cols, cex=0.3)
tiplabels(pie = to.matrix(zygo_terr_mode, sort(unique(zygo_terr_mode))), piecol=cols, cex=0.3)
#likely ancestral state = non-territorial! - that's interesting because it differs from Anisoptera



#plot this with a fan shape
cols<-setNames(c("turquoise", "brown"), levels(zygo_terr_mode))
fit_ARD_again<-ace(zygo_terr_mode, zygo_tree, model="ARD", type="discrete")
round(fit_ARD_again$lik.anc, 3)
plotTree(zygo_tree, type="fan", fsize=0.6, ftype="i")
nodelabels(node=1:zygo_tree$Nnode+Ntip(zygo_tree),
           pie=fit_ARD_again$lik.anc, piecol = cols, cex=0.3)
tiplabels(pie=to.matrix(zygo_terr_mode, sort(unique(zygo_terr_mode))), piecol=cols, cex=0.3)
legend("topright", legend=levels(zygo_terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)

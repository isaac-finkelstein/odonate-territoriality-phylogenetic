# Final R script

#evolutionary models, ancestral state reconstruction, phylogenetic signal
#note that I dropped courtship from the manuscript because it had a small sample size which limited any potential inferences.


library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(dplyr)
library(corHMM)
library(caper)
library(tidyr)
library(ggplot2)

my_data<- read.csv("data/data_vfinal.csv")

tree <- read.tree("data/iqtree_odonate_output/odonatetree.nex.treefile")
str(tree) #plot(tree, no.margin=TRUE) #tree$Nnode #tree$tip.label #check the structure of the tree

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

#plot the tree
plot(odonate_tree, type="fan", cex=0.5, fsize=0.6, ftype="i")
str(odonate_tree) #371 species/tips

#making territoriality a factor = yes/no
prop_terr_factor <- factor(ifelse(sp_terr == 1, "yes", "no"))
terr_data_with_na_factor <- data.frame(sn, prop_terr_factor = prop_terr_factor) # Create a data frame
terr_data_factor <- terr_data_with_na_factor[complete.cases(terr_data_with_na_factor), ] #Removes NA values
odonate_tree_factor<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree_factor
species_to_drop <- chk$data_not_tree
odonate_terr_data_factor <- terr_data_factor[!(terr_data_factor$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_terr_data_factor, data.names=as.character(odonate_terr_data_factor$sn))
terr_mode<-setNames(odonate_terr_data_factor$prop_terr, odonate_terr_data_factor$sn)


#Choose a discrete character model (rate-of-transition model)

#four possible models for discreet binary response:
#equal rates (ER)
#all rates different (ARD)
#two irreversible models: going from 0 -> 1,or from 1 -> 0, but not vice versa
terr_mode<-setNames(terr_data_factor$prop_terr, terr_data_factor$sn)

fit_er<-fitMk(odonate_tree, terr_mode, model = "ER")
fit_ard<-fitMk(odonate_tree, terr_mode, model = "ARD")
#fit the no -> yes model (transition from non-territorial to territorial but not back)
fit_01<-fitMk(odonate_tree, terr_mode, model=matrix(c(0,1,0,0),2,2,byrow=TRUE))
#fit the yes to no model (transition from territorial to non-territorail but not back)
fit_10<-fitMk(odonate_tree, terr_mode, model=matrix(c(0,0,1,0),2,2, byrow=TRUE))


#extract AIC values from the models
aic_model<-c(AIC(fit_er), AIC(fit_ard), AIC(fit_01), AIC(fit_10))
data.frame(model=c("ER", "ARD", "no -> yes", "yes->no"),
           logL=c(logLik(fit_er), logLik(fit_ard),
                  logLik(fit_01), logLik(fit_10)),
           AIC=aic_model,delta.AIC=aic_model-min(aic_model))   
#the all-rates-different model has the lowest AIC

#graph it
rates <- fit_ard$rates  # get transition rates
print(rates)


max_rate <- max(rates, na.rm = TRUE)
scaled_lwd <- rates / max_rate * 3

legend_rates <- pretty(range(rates, na.rm=TRUE), n=3)
legend_lwd <- legend_rates / max_rate * 3

plot(fit_ard, show.zeros=FALSE, mar=rep(0,4), signif=5, lwd=scaled_lwd)
legend("topright", legend=sprintf("%.3f", legend_rates), 
       lwd=legend_lwd, col="black", title="Transition Rate")
#the line thickness is flipped, and I can't figure out why. 


#estimate q (the transition rate)
#Q is the expected number of transition given a particular amount of time 
print(fit_er)
print(fit_ard)
print(fit_01)
print(fit_10)

#graph it:
plot(fit_ard, show.zeros=FALSE, mar=rep(0,4), signif=5)
#the ARD (all rates different) model, which had the lowest AIC, shows that it is
#easier/more common to go No -> Yes than to go Yes -> No

#Ancestral state reconstruction
#using marginal ancestral state reconstruction as this is more popular
#and provides uncertainty about specific values

#create a special data frame of species names and trait data
#character data must be numerical interger (0,1)
odonate_data<-data.frame(Genus_sp=names(terr_mode), terr_mode=as.numeric(terr_mode)-1)
head(odonate_data, n=10)

#prune tree to match data
chk<-name.check(tree, odonate_data, data.names=as.character(odonate_data$Genus_sp))
summary(chk)
odonate_tree<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree
# Identify species to drop
species_to_drop <- chk$data_not_tree
odonate_data <- odonate_data[!(odonate_data$Genus_sp%in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_data, data.names=as.character(odonate_data$Genus_sp))


#fit model
fit_marginal<- corHMM(odonate_tree, odonate_data, node.states = "marginal",
                      rate.cat=1, rate.mat=NULL)
fit_marginal
head(fit_marginal$states) 

#the node for the common ancestor of all odonates = 375, all Anisoptera = 376, all zygoptera = 568
plot(odonate_tree, show.node.label = TRUE)  # Plot tree with tip labels
nodelabels() 

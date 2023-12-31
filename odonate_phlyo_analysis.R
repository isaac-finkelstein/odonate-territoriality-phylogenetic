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
prop_terr<- round(terr_table[,3]/(terr_table[,2]+ terr_table[,3]),0) # so this calculates the percentage "yes" for territorial
#note that I am rounding to 1 whole number- so 0.5-> 1. Figure out how to set this to a 66% threshold later
sn<- attr(terr_table,"row.vars")[[1]]
terr_data_with_na<- data.frame(sn,prop_terr) #so this is a dataframe with the proporitons of territorial "yes" for each species
terr_data<- terr_data_with_na[complete.cases(terr_data_with_na), ] #removed NA values 

#to make the territorial variable binary (yes or no) make it a factor
prop_terr_factor <- factor(ifelse(prop_terr == 1, "yes", "no"))
terr_data_with_na_factor <- data.frame(sn, prop_terr_factor = prop_terr_factor)
terr_data_factor <- terr_data_with_na_factor[complete.cases(terr_data_with_na_factor), ]

# Create a data frame
terr_data_with_na <- data.frame(sn, prop_terr = prop_terr_factor)

# Remove NA values
terr_data_factor <- terr_data_with_na[complete.cases(terr_data_with_na), ]

#prune tree to match data
chk<-name.check(tree, terr_data, data.names=as.character(terr_data$sn))
summary(chk)
odonate_tree<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree
# Identify species to drop from terr_data
species_to_drop <- chk$data_not_tree
odonate_terr_data <- terr_data[!(terr_data$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
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

#let's plot territoriality on my tree
#I want to use the factor (binary yes vs no)
odonate_tree_factor<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree_factor
species_to_drop <- chk$data_not_tree
odonate_terr_data_factor <- terr_data_factor[!(terr_data_factor$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_terr_data_factor, data.names=as.character(odonate_terr_data_factor$sn))


#let's plot territoriality on my tree
terr_mode<-setNames(odonate_terr_data_factor$prop_terr, odonate_terr_data_factor$sn)
cols<-setNames(c("forestgreen", "goldenrod1"), levels(terr_mode))
plotTree.datamatrix(odonate_tree_factor, as.data.frame(terr_mode),
                    colours= list(cols), header=FALSE, fsize=0.45)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
#it's so crowded that I might want to remove names and then provide a large list so people can look up whatever species they want

#Just Anisoptera
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
cols<-setNames(c("forestgreen", "goldenrod1"), levels(anis_terr_mode))
plotTree.datamatrix(anis_tree, as.data.frame(anis_terr_mode),
                    colours= list(cols), header=FALSE, fsize=0.45)
legend("topright", legend=levels(anis_terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)

#this seems to be working but the colours are wrong
#the solution to make it work was changing the setNames function line
#it's so crowded that I might want to remove names and then provide a large list so people can look up whatever species they want



#choose a character model
#4 possible models for discreet repsonse
#equal rates (ER)
#all rates different (ARD)
# irreversible modes (2) - going from 0 -> 1 or from 1 -> 0 but not vice versa)

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
#0kay, so it seems that the all-rates-different is the way to go   
#ARD is quite strongly favoured with a differnce of 22 with the next lowest AIC. 

#we can estimate q (the transition rate)
#Q ia the expecrted number of transition given a particular amount of time 
print(fit_er)
print(fit_ard)
print(fit_01)
print(fit_10)

#graph it:
plot(fit_ard, show.zeros=FALSE, mar=rep(0,4), signif=5)
plot(fit_er)
#the ARD (all rates different) model, which had the lowest AIC -- shows that it is
#easier/more common to go No -> Yes than to go Yes -> No


#Next: Ancestral state reconstruction - plots are very busy so may want to re-run these for just anisoptera/zygoptera
#I have to decide if I want to do joint or marginal acenstral state reconstruction
#marginal is more popular in biology (Revell and Harmon, 2022)
#marginal measures uncertainty about the specific values for ancestral states - so preferred
#do both, then compare results (joint can go in supplementary material)
library(corHMM)
#create a special data frame of species names and trait data
#character data must be numerical interger (0,1,2)
odonate_data<-data.frame(Genus_sp=names(terr_mode), terr_mode=as.numeric(terr_mode)-1)
head(odonate_data, n=10)


#fit the joint model
fit_joint<-corHMM(odonate_tree, odonate_data, node.states="joint", rate.cat=1, model="ARD",rate.mat=NULL)
fit_joint
#now plot
plotTree.datamatrix(odonate_tree, as.data.frame(terr_mode),
                    colors=list(cols), header=FALSE, fsize=0.45)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
nodelabels(pie=to.matrix(levels(terr_mode)[fit_joint$phy$node.label],
                         levels(terr_mode)), piecol=cols, cex=0.1)
#it's very small and hard to see -- adjust cex to make the pie charts larger

#let's do marginal ancestral state reconstruction
fit_marginal<- corHMM(odonate_tree, odonate_data, node.states = "marginal",
                      rate.cat=1, rate.mat=NULL)
fit_marginal
head(fit_marginal$states)       
#interpret this matrix as the posterior probabilities that each state is in each node

#plot this
plotTree.datamatrix(odonate_tree, as.data.frame(terr_mode),
                    colors=list(cols), header=FALSE, fsize=0.45)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
nodelabels(pie=fit_marginal$states, piecol=cols, cex=0.3)
#likely ancestral state = territorial!

#one more method using an MCMC approach: Stochasitic character mapping
#so you take the whole distribution from many sample stochastic maps and sample
#then you use Bayesian MCMC to sample q from its posterior distribution
#this tells us the number of changes, the proportion of time spent in each state,
#and the posterior probabilities of each node being a certain state
#this tells us the same thing as the marginal ancestral states
#but the transition rate between states (q) is the maximum likelihood value
#we sampled q from it's posterior distribution
#see Huelsenbeck et al., 2003


#generate 1000 stochastic character maps in which the transition rate is sampled from its posterior distribution
#this may work, but it took over an hour to run and didn't finished
#mtrees<-make.simmap(odonate_tree, terr_mode, model="ARD", nsim=1000, Q="mcmc", vQ=0.01,
#                    prior=list(use.empirical=TRUE), samplefreq=10)


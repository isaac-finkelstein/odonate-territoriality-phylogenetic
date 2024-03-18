#evolutionary models, ancestral state reconstruction, phylogenetic signal

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

#I can also try the tree from Rocha-Ortega et al., 2020 in Proceedings
#it uses the same structure as Waller and Svensson 2017 but has different species coverage.
#tree<-read.nexus(file="data/Rocha_ortega_tree.nex")

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

#let's plot the tree
plot(odonate_tree, type="fan", cex=0.5, fsize=0.6, ftype="i") #it's very busy because it has 372 tips


#EVERYTIME I UPDATE THE DATASET, THE NODE NUMBER WILL CHANGE!! SO HAVE TO UPDATE THESE EVERYTIME!
#collapseTree(odonate_tree)  #-- use this to find the node numbers

#just anisoptera -- node 411 
anisoptera_tree_extract<-unname(ape::extract.clade(odonate_tree, node = 560)$tip.label)
anisoptera_tree<-drop.tip(odonate_tree, anisoptera_tree_extract)
plot(anisoptera_tree, type="fan", cex=0.5, ftype="i") 

#just zygoptera -- node 366
zygoptera_tree_extract<-unname(ape::extract.clade(odonate_tree, node = 371)$tip.label)
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
                    colours= list(cols), header=FALSE)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
#it's so crowded that I might want to remove names and then provide a large list so people can look up whatever species they want
terr_mode
#Plot this in a fan tree format
plotTree(odonate_tree_factor, type="fan", fsie=0.8, ftype="i")
cols<-setNames(palette()[1:length(unique(terr_mode))],sort(unique(terr_mode)))
tiplabels(pie=to.matrix(terr_mode, sort(unique(terr_mode))), piecol=cols, cex=0.2)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
#garbage - just delete this


#choose a character model
#4 possible models for discreet response
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
#Q is the expected number of transition given a particular amount of time 
print(fit_er)
print(fit_ard)
print(fit_01)
print(fit_10)

#graph it:
plot(fit_ard, show.zeros=FALSE, mar=rep(0,4), signif=5)
plot(fit_er)
#the ARD (all rates different) model, which had the lowest AIC -- shows that it is
#easier/more common to go No -> Yes than to go Yes -> No

#One more model:
#ordered evolution model
#the evolution of territoriality should be ordered: strong->weak->no territoraility (or no->weak->strong)
#my data set for this delineation between strong and weak territoriality is more limited

#make a single column depicting territorial: no, weak, strong = ordered categorical variable
#strength of territoriality
strength_terr<-ftable(my_data$Formatted_species, my_data$Strength.of.territoriality..strong.vs.weak.)
prop_strength<-round(strength_terr[,3]/(strength_terr[,2]+strength_terr[,3]),2) #this is proportion weak
#so 1=strong, 2 = weak??? FIGURE THIS OUT
#the 2 at the end rounds to 2 decimal 
#set a 75% threshold = 3:1 threshold
sp_strength_terr<-ifelse(prop_strength >= 0.75, 2, ifelse(prop_strength <=0.25, 3, NA)) # so 2 = weak, 3 = strong
sn<-attr(strength_terr, "row.vars")[[1]]
strength_data_with_na<-data.frame(sn, sp_strength_terr, stringsAsFactors=TRUE)
strength_data<- strength_data_with_na[complete.cases(strength_data_with_na), ] #removed NA values 
#this worked. 


#next combine into a new column so I have three levels: "no, weak, strong"

strength_of_terr_data<- merge(terr_data, strength_data, by="sn", all= TRUE)
strength_of_terr_data$combined_strength<-with(strength_of_terr_data, ifelse(sp_terr=="0", "0", sp_strength_terr))
strength_of_terr_data <- strength_of_terr_data[, !(names(strength_of_terr_data) %in% c("sp_terr", "sp_strength_terr"))]
strength_of_terr_data <- na.omit(strength_of_terr_data)
#to not get confused with territorial? No = 0, yes = 1, for strength I set weak = 2, strong = 3
#This is just trying something: I'm making non-territorial to be 1, weak to be 2, strong to be 3
#let's turn it back to 0,1,2 because that just makes sense
strength_of_terr_data <- strength_of_terr_data %>%
  mutate(
    combined_strength = as.numeric(combined_strength),
    combined_strength = case_when(
      combined_strength == 0 ~ 1,
      TRUE ~ combined_strength
    )
  )
#note that there is more "no" here because it was easier to find binary data - yes/no, than it was to find data with strong/weak territoriality
#so "no" is over-represented in the data.

rownames(strength_of_terr_data) <- strength_of_terr_data$sn
strength_of_terr_data$sn <- NULL  # Remove the "sn" column

#I need to prune a tree to the data.
chk_strength<-name.check(tree, strength_of_terr_data)
summary(chk_strength)
strength_tree<-drop.tip(tree, chk_strength$tree_not_data) #dropped tree_not_data species
strength_tree
# Identify species to drop from terr_data
#species_to_drop_strength <- chk_strength$data_not_tree
strength_terr_data_pruned <- strength_of_terr_data[!(rownames(strength_of_terr_data) %in% 
                                                         chk_strength$data_not_tree),,drop=FALSE] #dropped data_not_tree species from dataset
strength_data_pruned<-setNames(as.factor(strength_terr_data_pruned[,"combined_strength"]), rownames(strength_terr_data_pruned))
name.check(strength_tree, strength_terr_data_pruned)

#ordered model
ordered_model<-matrix(c(
  0,1,0,
  2,0,3,
  0,4,0), 3,3, byrow=TRUE,
  dimnames=list(1:3, 1:3))
ordered_model

fit_ordered<-fitDiscrete(strength_tree, strength_terr_data_pruned, model=ordered_model, suppressWarnings=TRUE)
#not working
#could be the model? Could be a coding error? I'm not sure. 
#could also plot it.





#Next: Ancestral state reconstruction - plots are very busy so may want to re-run these for just anisoptera/zygoptera
#I have to decide if I want to do joint or marginal acenstral state reconstruction
#marginal is more popular in biology (Revell and Harmon, 2022)
#marginal measures uncertainty about the specific values for ancestral states - so preferred
#do both, then compare results (joint can go in supplementary material)

#create a special data frame of species names and trait data
#character data must be numerical interger (0,1)
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
                         levels(terr_mode)), piecol=cols, cex=0.2)
#adjust cex to make the pie charts larger
#if I can figure out how to remove species labels, it might be interesting to present the whole odonate tree for ASR
#because you can see that zygoptera and Anisoptera have different ancestral states!

#let's do marginal ancestral state reconstruction
fit_marginal<- corHMM(odonate_tree, odonate_data, node.states = "marginal",
                      rate.cat=1, rate.mat=NULL)
fit_marginal
head(fit_marginal$states)       
#interpret this matrix as the posterior probabilities that each state is in each node

cols<-setNames(c("turquoise", "brown"), levels(terr_mode))
#plot this
plotTree.datamatrix(odonate_tree, as.data.frame(terr_mode),
                    colors=list(cols), header=FALSE, fsize=0.45)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
nodelabels(pie=fit_marginal$states, piecol=cols, cex=0.3)
tiplabels(pie=to.matrix(terr_mode, sort(unique(terr_mode))), piecol=cols, cex=0.3)
#likely ancestral state = territorial!
#if I can figure out how to remove species labels, it might be interesting to present the whole odonate tree for ASR
#because you can see that zygoptera and Anisoptera have different ancestral states!

#plot this with a fan shape
fit_ARD_again<-ace(terr_mode, odonate_tree, model="ARD", type="discrete")
round(fit_ARD_again$lik.anc, 3)
plotTree(odonate_tree, type="fan", fsize=0.5, ftype="i")
nodelabels(node=1:odonate_tree$Nnode+Ntip(odonate_tree),
           pie=fit_ARD_again$lik.anc, piecol = cols, cex=0.3)
tiplabels(pie=to.matrix(terr_mode, sort(unique(terr_mode))), piecol=cols, cex=0.3)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
#this is great! still too busy but if I remove the species names it looks okay


#http://www.phytools.org/eqg2015/asr.html
#NOTICE that the results differ a bit. The estimate ancestral state changes between ace and corHMM.



#calculating delta, a measure of phylogenetic signal
#from Borges et al., 2019 github
#https://github.com/mrborges23/delta_statistic
source("code.R") ## This is taken from mrborges23's github https://github.com/mrborges23/delta_statistic/tree/master

#all branches must be positive
#DONT UNDERSTAND THIS STEP BUT COPYING FROM THE README FILE.
odonate_tree$edge.length[odonate_tree$edge.length==0] <-quantile(odonate_tree$edge.length, 0.1)*0.1

#next define the trait vector
#odonate_tree$tip.label #shows the species order, note it's alphabetical, split into Zygoptera and Anisoptera
order_in_terr_data <- match(odonate_tree$tip.label, odonate_terr_data$sn) #must be same order as tree
odonate_terr_data_reordered <- odonate_terr_data[order_in_terr_data, ] #now tree and data are in the same order
trait_with_na<-odonate_terr_data_reordered$sp_terr #this is a vector of territoriality, with 1= yes, 0 = no, at a 3:1 threshold in my dataset
trait<- trait_with_na[!is.na(trait_with_na)] #removed NA values 
#set at 3:1 threshold (75% concordance in the dataset)

#now calculate delta
deltaA<-delta(trait, odonate_tree, 0.1, 0.0589, 10000, 10, 100)
#in the example in the README, he uses 0.0589 for the SE, but in example.R file
#he uses 0.5, unless I have misunderstood this.
# WHICH TO USE? -- when I try both, I get basically the same result (but slightly different)
print(deltaA)
#I get 8.2995. What does this mean?
#Borges et al 2019: "delta can be any positive real number:
#"the higher the delta-value, the higher the degree of phylogenetic signal between a given trait and the phylogeny"
#so this means that 0 = no phylogenetic signal? in their example, they had 1 = no phylogenetic signal.

#note that delta should be used with at least 20 species - i have way more so I am good here.

#finally, we can calculate a p-value for this
#First, create a vector of random deltas that will be our null hypothesis
random_delta<-rep(NA,100)
for (i in 1:100){
  rtrait <- sample(trait)
  random_delta[i] <- delta(rtrait, odonate_tree, 0.1, 0.0589, 10000, 10, 100)
}
p_value<-sum(random_delta>deltaA)/length(random_delta) #results in a p-value of 0...did I do something wrong?
#random_delta>deltaA
#there are no instances of random_delta>deltaA
boxplot(random_delta)
abline(h=deltaA, col="red")


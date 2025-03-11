#evolutionary models, ancestral state reconstruction, phylogenetic signal

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(dplyr)
library(corHMM)
library(caper)
library(tidyr)

my_data<- read.csv("data/data_v6.csv")

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

#_________________________________________________________________________________________________
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
rates <- fit_ard$rates  # get transition rates
print(rates)

max_rate <- max(rates, na.rm = TRUE)
scaled_lwd <- rates / max_rate * 3


legend_rates <- pretty(range(rates, na.rm=TRUE), n=4)
legend_lwd <- legend_rates / max_rate * 3

plot(fit_ard, show.zeros=FALSE, mar=rep(0,4), signif=5, lwd=scaled_lwd)

legend("topright", legend=sprintf("%.4f", legend_rates), 
       lwd=legend_lwd, col="black", title="Transition Rate")
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




#______________________________________________________________________________________________________________________________
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
figure1<-plotTree.datamatrix(odonate_tree, as.data.frame(terr_mode),
                    colors=list(cols), header=FALSE, fsize=0.45)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
nodelabels(pie=fit_marginal$states, piecol=cols, cex=0.3)
tiplabels(pie=to.matrix(terr_mode, sort(unique(terr_mode))), piecol=cols, cex=0.3)

#likely ancestral state = territorial!
#if I can figure out how to remove species labels, it might be interesting to present the whole odonate tree for ASR
#because you can see that zygoptera and Anisoptera have different ancestral states!


#Adding trait data to this figure
#I created the trait_data dataframe in "Testing correlations between variables" file - go there (at the very end)

#trying to just make the matrix, no phylogeny
trait_data <- trait_data %>%
  mutate(Species = as.character(Species))
trait_data <- trait_data %>%
  mutate(across(c(Prop_Territorial, Prop_Mate_Guard, Prop_Flier_vs_Percher, 
                  Prop_Oviposition, Prop_Courtship, Mate_guarding_cat, 
                  lentic_lotic_size, Prop_lo_len), as.character))
trait_data_long <- trait_data %>%
  pivot_longer(cols = c(Prop_Territorial, Prop_Mate_Guard, Prop_Flier_vs_Percher, Prop_Oviposition, Prop_Courtship, Mate_guarding_cat, lentic_lotic_size, Prop_lo_len),
               names_to = "Trait",
               values_to = "Value") %>%
  mutate(Value = as.character(Value))

# it needs to be in the same order as the figure
species_order_from_fig1 <- odonate_tree$tip.label

# remove any Species that do not appear in the figure
trait_data_long <- trait_data_long %>%
  filter(Species %in% species_order_from_fig1)

# convert Species column to a factor with levels matching the tree order
trait_data_long$Species <- factor(trait_data_long$Species, levels = species_order_from_fig1)

#It needs to appear in the same order as the figure 1
# Extract the correct plotting order
tip_order <- odonate_tree$tip.label[odonate_tree$edge[odonate_tree$edge[,2] <= length(odonate_tree$tip.label), 2]]

#check order
print(tip_order) #This is the order for figure 1, with the last species at the top of the figure and the first species at the bottom of the figure


#match order
trait_data_long$Species <- factor(trait_data_long$Species, levels = tip_order)


#Plot heatmap

#predictor traits
trait_data_filtered <- trait_data_long %>%
  filter(Trait %in% c("lentic_lotic_size", "Mate_guarding_cat", "Prop_Flier_vs_Percher", "Prop_Oviposition", "Prop_Courtship", "Prop_lo_len")) 

# Convert NA values to a string ("NA") for ggplot to recognize them
trait_data_filtered <- trait_data_filtered %>%
  mutate(Value = ifelse(is.na(Value), "NA", Value))  # Convert NA to string


trait_data_filtered <- trait_data_filtered %>%
  mutate(Value = ifelse(grepl("^.*_", Value), Value, paste(Trait, Value, sep = "_")))

trait_data_filtered <- trait_data_filtered %>%
  mutate(Value = ifelse(Trait == "lentic_lotic_size" & !grepl("^lentic_lotic_size_", Value), 
                        paste("lentic_lotic_size", Value, sep = "_"), 
                        Value))

#order of presenting the traits
trait_order <- c("Prop_Flier_vs_Percher", "Prop_Oviposition", 
                 "Prop_lo_len", "Prop_Courtship", "Mate_guarding_cat", "lentic_lotic_size")
trait_data_filtered$Trait <- factor(trait_data_filtered$Trait, levels = trait_order)

lentic_colors <- c(
  "lentic_lotic_size_lentic_medium" = "darkred",
  "lentic_lotic_size_lentic_large" = "darkorange",
  "lentic_lotic_size_lentic_small" = "darkgreen",
  "lentic_lotic_size_stream" = "blue",
  "lentic_lotic_size_river" = "lightblue",
  "lentic_lotic_size_generalist" = "black",
  "lentic_lotic_size_both" = "aquamarine2"
)


guarding_colors <- c(
  "Mate_guarding_cat_No" = "cyan3",  
  "Mate_guarding_cat_Contact" = "firebrick1",
  "Mate_guarding_cat_Non-contact" = "mediumblue"
)

# Colours for Flier vs Percher
flying_colors <- c(
  "Prop_Flier_vs_Percher_0" = "pink",  # Flier
  "Prop_Flier_vs_Percher_1" = "purple" # Percher
)

# Colours for Oviposition
ovi_colors <- c(
  "Prop_Oviposition_0" = "gray34",  # Endophytic
  "Prop_Oviposition_1" = "coral4"    # Exophytic
)

# Colours for Courtship
court_colors <- c(
  "Prop_Courtship_0" = "darkslateblue",  # No
  "Prop_Courtship_1" = "khaki3"    # Yes
)

# Colours for lotic vs lentic
lo_len_colors <- c(
  "Prop_lo_len_0" = "seagreen",  # lentic
  "Prop_lo_len_1" = "coral"    # lotic
)

# Combine color schemes
all_trait_colors <- c(lentic_colors, guarding_colors, flying_colors, ovi_colors, court_colors, lo_len_colors, "NA" = "white")  # Assign white for missing values

#plot
ggplot(trait_data_filtered, aes(x = Trait, y = Species, fill = Value)) +
  geom_tile(color = "white") +  
  scale_fill_manual(values = all_trait_colors, na.value = "white") +
  scale_x_discrete(labels = c(
    "Prop_Flier_vs_Percher" = "Active behaviour",
    "Prop_Courtship" = "Courtship",
    "Prop_Oviposition" = "Oviposition method",
    "Prop_lo_len" = "Oviposition habitat type",
    "Mate_guarding_cat" = "Mate guarding",
    "lentic_lotic_size" = "Oviposition habitat size"
  )) + 
  theme_minimal() +
  labs(title = "Predictor trait Heatmap",
       x = "Trait",
       y = "Species",
       fill = "Trait Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 6))





#try using the ace function from the package "ape"
fit_ARD_again<-ace(terr_mode, odonate_tree, model="ARD", type="discrete", marginal = TRUE)
round(fit_ARD_again$lik.anc, 3)
plotTree(odonate_tree, fsize=0.5, ftype="i")
nodelabels(node=1:odonate_tree$Nnode+Ntip(odonate_tree),
           pie=fit_ARD_again$lik.anc, piecol = cols, cex=0.3)
tiplabels(pie=to.matrix(terr_mode, sort(unique(terr_mode))), piecol=cols, cex=0.3)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)

#plot this with a fan shape
plotTree(odonate_tree, type="fan", fsize=0.5, ftype="i")
nodelabels(node=1:odonate_tree$Nnode+Ntip(odonate_tree),
           pie=fit_ARD_again$lik.anc, piecol = cols, cex=0.3)
tiplabels(pie=to.matrix(terr_mode, sort(unique(terr_mode))), piecol=cols, cex=0.3)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
#this is great! still too busy but if I remove the species names it looks okay

#http://www.phytools.org/eqg2015/asr.html
#NOTICE that the results differ a bit. The estimate ancestral state changes between ace and corHMM.





#using stochastic character mapping to reveal the estimated history frequency of the transitions in the tree (the number of transitions that have occures)
terr_mode <- setNames(as.character(terr_mode), names(terr_mode))

odonate_simmaps <- make.simmap(odonate_tree, terr_mode, model="ARD", nsim=100)

summary(odonate_simmaps)

#________________________________________________________________________________________________________
#calculating phylogenetic signal using the D statistic from Fritz & Purvis 2010
#terr_data_factor$names <- terr_data_factor$sn
#terr_data_factor <- terr_data_factor[, -which(names(terr_data_factor) == "sn")]
signal_terr<-phylo.d(odonate_terr_data, odonate_tree, names.col=sn, binvar=sp_terr)

#do this for my other traits
#mate guarding 
#remove "No" and "Both". Below, I include "no" as a third category. this compares mate guarding as a binary variable contact vs non-contact
#there is only one instance of "both"
#note that this variable is only for species that exhibit mate guarding 
my_data$Mate.guarding <- trimws(my_data$Mate.guarding)
filtered_mate_guarding <- subset(my_data, Mate.guarding %in% c("Contact", "Non-contact"))
mate_guard_var<-ftable(filtered_mate_guarding$Formatted_species, filtered_mate_guarding$Mate.guarding)
prop_mate_guard<-round(mate_guard_var[,1]/(mate_guard_var[,1]+mate_guard_var[,2]),2) #this is the proportion that is contact
#so 1 = contact, 0 = non-contact. 
#the 2 at the end rounds to 2 decimal 
#set a 75% threshold = 3:1 threshold
sp_binary_mate_guard<-ifelse(prop_mate_guard >= 0.75, 1, ifelse(prop_mate_guard <=0.25, 0, NA))
sn<-attr(mate_guard_var, "row.vars")[[1]]
binary_mate_guard_df<-data.frame(sn, sp_binary_mate_guard, stringsAsFactors = TRUE)
#this worked, I checked manually.

#Flier vs percher
binary_fly_v_perch<-ftable(my_data$Formatted_species, my_data$Flier.vs.percher)
prop_fly_v_perch<-round(binary_fly_v_perch[,3]/(binary_fly_v_perch[,2]+binary_fly_v_perch[,3]),2) #this is proportion Percher
#so 1=percher, 0 = flier
#the 2 at the end rounds to 2 decimal 
#set a 75% threshold = 3:1 threshold
sp_fly_v_perch<-ifelse(prop_fly_v_perch >= 0.75, 1, ifelse(prop_fly_v_perch <=0.25, 0, NA))
sn<-attr(binary_fly_v_perch, "row.vars")[[1]]
binary_fly_v_perch_df<-data.frame(sn, sp_fly_v_perch, stringsAsFactors=TRUE)
#this worked, but it includes all species - so there are NAs. 

#courtship
binary_courtship<-ftable(my_data$Formatted_species, my_data$Courtship)
prop_court<-round(binary_courtship[,3]/(binary_courtship[,2]+binary_courtship[,3]),2) #this is proportion "Yes"
#so 1=yes, 0 = no
#set a 75% threshold = 3:1 threshold
sp_courtship<-ifelse(prop_court>=0.75,1, ifelse(prop_court<=0.25, 0, NA))
sn<-attr(binary_courtship, "row.vars")[[1]]
binary_court_df<-data.frame(sn, sp_courtship, stringsAsFactors=TRUE) #includes NAs

#oviposition (endophytic vs exophytic)
#I need to make all epiphytic = exophytic - to make it binary. 
my_data_mutated <- my_data %>%
  mutate(Oviposition.type..endophytic.vs.exophytic. = ifelse(Oviposition.type..endophytic.vs.exophytic. %in% c("Epiphytic", "Exophytic"), "Exophytic", Oviposition.type..endophytic.vs.exophytic.))
# Convert the variable to a factor with specified levels
my_data_mutated$Oviposition.type..endophytic.vs.exophytic. <- factor(
  my_data_mutated$Oviposition.type..endophytic.vs.exophytic.,
  levels = c("Endophytic", "Exophytic")
)
binary_ovi<-ftable(my_data_mutated$Formatted_species, my_data_mutated$Oviposition.type..endophytic.vs.exophytic.)
prop_ovi<-round(binary_ovi[,2]/(binary_ovi[,1]+binary_ovi[,2]),2) #This is proportion exophytic
#so 1=exophytic, 0=endophytic
#the 2 at the end rounds to 2 decimal 
#set a 75% threshold = 3:1 threshold
sp_ovi<-ifelse(prop_ovi >= 0.75, 1, ifelse(prop_ovi <=0.25, 0, NA))
sn<-attr(binary_ovi, "row.vars")[[1]]
binary_ovi_df<-data.frame(sn, sp_ovi, stringsAsFactors=TRUE)
#this worked, but it includes all species - so there are NAs. 

filtered_lo_len <- subset(my_data, Lotic.vs.lentic..breeding.habitat. %in% c("Lotic", "Lentic"))
lo_len_var<-ftable(filtered_lo_len$Formatted_species, filtered_lo_len$Lotic.vs.lentic..breeding.habitat.)
prop_lo_len<-round(lo_len_var[,2]/(lo_len_var[,1]+lo_len_var[,2]),2) #this is the proportion that is lotic
#the 2 at the end rounds to 2 decimal 
#set a 75% threshold = 3:1 threshold
sp_lo_len<-ifelse(prop_lo_len >= 0.75, 1, ifelse(prop_lo_len <=0.25, 0, NA))
sn<-attr(lo_len_var, "row.vars")[[1]]
binary_lo_len<-data.frame(sn, sp_lo_len, stringsAsFactors = TRUE)

#now prune them to match the tree
# Identify species to drop from mate_guard_terr_data_old
chk_mate_guard<-name.check(tree, binary_mate_guard_df, data.names=binary_mate_guard_df$sn)
tree_mate_guard <- drop.tip(tree, chk_mate_guard$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
mate_guard_species_to_drop<-chk_mate_guard$data_not_tree
mate_guard<-binary_mate_guard_df[!(binary_mate_guard_df$sn %in% mate_guard_species_to_drop),] #dropped data_not_tree species from dataset
name.check(tree_mate_guard, mate_guard, data.names=as.character(mate_guard$sn))

chk_fly_v_perch<-name.check(tree, binary_fly_v_perch_df, data.names=as.character(binary_fly_v_perch_df$sn))
tree_fly_v_perch <- drop.tip(tree, chk_fly_v_perch$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
fly_v_perch_species_to_drop<-chk_fly_v_perch$data_not_tree
fly_v_perch_data<-binary_fly_v_perch_df[!(binary_fly_v_perch_df$sn %in% fly_v_perch_species_to_drop),] #dropped data_not_tree species from dataset
name.check(tree_fly_v_perch, fly_v_perch_data, data.names=as.character(fly_v_perch_data$sn))

chk_court<-name.check(tree, binary_court_df, data.names=as.character(binary_court_df$sn))
tree_court <- drop.tip(tree, chk_court$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
court_species_to_drop<-chk_court$data_not_tree
court_data<-binary_court_df[!(binary_court_df$sn %in% court_species_to_drop),] #dropped data_not_tree species from dataset
name.check(tree_court, court_data, data.names=as.character(court_data$sn))

chk_ovi<-name.check(tree, binary_ovi_df, data.names=as.character(binary_ovi_df$sn))
tree_ovi <- drop.tip(tree, chk_ovi$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
ovi_species_to_drop<-chk_ovi$data_not_tree
ovi_data<-binary_ovi_df[!(binary_ovi_df$sn %in% ovi_species_to_drop),] #dropped data_not_tree species from dataset
name.check(tree_ovi,ovi_data, data.names=as.character(ovi_data$sn))

chk_lo_len<-name.check(tree, binary_lo_len, data.names=as.character(binary_lo_len$sn))
tree_lo_len <- drop.tip(tree, chk_lo_len$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
lo_len_species_to_drop<-chk_lo_len$data_not_tree
lo_len_data<-binary_lo_len[!(binary_lo_len$sn %in% lo_len_species_to_drop),] #dropped data_not_tree species from dataset
name.check(tree_lo_len, lo_len_data, data.names=as.character(lo_len_data$sn))



#Check the phylogenetic signal of these binary traits
signal_mate_guard<-phylo.d(mate_guard, tree_mate_guard, names.col=sn, binvar=sp_binary_mate_guard)
signal_fly_v_perch<-phylo.d(fly_v_perch_data, tree_fly_v_perch, names.col=sn, binvar=sp_fly_v_perch)
signal_courtship<-phylo.d(court_data, tree_court, names.col=sn, binvar=sp_courtship)
signal_oviposition<-phylo.d(ovi_data, tree_ovi, names.col=sn, binvar = sp_ovi)
signal_lo_len<-phylo.d(lo_len_data, tree_lo_len, names.col = sn, binvar = sp_lo_len)



#Using another calculation of phlyogenetic signal
#calculating delta, a measure of phylogenetic signal
#from Borges et al., 2019 github
#https://github.com/mrborges23/delta_statistic
source("code.R") ## This is taken from mrborges23's github https://github.com/mrborges23/delta_statistic/tree/master

#all branches must be positive
#DONT UNDERSTAND THIS STEP BUT COPYING FROM THE README FILE.
odonate_tree$edge.length[odonate_tree$edge.length==0] <-quantile(odonate_tree$edge.length, 0.1)*0.1

#define the trait vector
#odonate_tree$tip.label #shows the species order, note it's alphabetical, split into Zygoptera and Anisoptera
order_in_terr_data <- match(odonate_tree$tip.label, odonate_terr_data$sn) #must be same order as tree
odonate_terr_data_reordered <- odonate_terr_data[order_in_terr_data, ] #now tree and data are in the same order
trait_with_na<-odonate_terr_data_reordered$sp_terr #this is a vector of territoriality, with 1= yes, 0 = no, at a 3:1 threshold in my dataset
trait<- trait_with_na[!is.na(trait_with_na)] #removed NA values 

#now calculate delta
deltaA<-delta(trait, odonate_tree, 0.1, 0.0589, 10000, 10, 100)
#in the example in the README, he uses 0.0589 for the SE, but in example.R file
#he uses 0.5, unless I have misunderstood this.
# WHICH TO USE? -- when I try both, I get basically the same result (but slightly different)
print(deltaA)
#I get 8.197. What does this mean?
#Borges et al 2019: "delta can be any positive real number:
#"the higher the delta-value, the higher the degree of phylogenetic signal between a given trait and the phylogeny"
#so this means that 0 = no phylogenetic signal? in their example, they had 1 = no phylogenetic signal.

#note that delta should be used with at least 20 species sample size - I have way more so I am good here.

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


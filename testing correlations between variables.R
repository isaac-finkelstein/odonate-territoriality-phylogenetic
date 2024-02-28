#Testing correlations between variables

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(dplyr)

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

#Pagel94 plot
#this is to see if there are evolutionary relationships between traits. Are certain combinations of traits more likely to evolve? 
#is there correlated evolution?
#the traits I am interested in are: mate guarding ~territoriality, flier vs percher ~ territoriality, oviposition (endo vs exo) ~ territoriality

#I have to make a new dataframe
#so I make individual dataframes for each trait I want, then stitch them together
#territoriality
binary_terr<-ftable(my_data$Formatted_species, my_data$Territorial)
prop_binary_terr<-round(binary_terr[,3]/(binary_terr[,2]+ binary_terr[,3]),2)
#the 2 at the end rounds to 2 decimal 
#set a 75% threshold = 3:1 threshold
sp_binary_terr<-ifelse(prop_binary_terr >= 0.75, 1, ifelse(prop_binary_terr <=0.25, 0, NA))
sn<-attr(binary_terr, "row.vars")[[1]]
binary_terr_df<-data.frame(sn,sp_binary_terr, stringsAsFactors = TRUE)
#so 1 = territorial, 0 = non-territorial. 

#mate guarding 
#remove "No" and "Both". I am only interested in comparing the binary variable contact vs non-contact
#there is only one instance of "both"
#so note that this is variable is only for species that exhibit mate guarding 
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

#now I need to stitch these together in a single dataframe
binary_data <- merge(binary_terr_df, binary_mate_guard_df, by = "sn", all = TRUE)
binary_data <- merge(binary_data, binary_fly_v_perch_df, by = "sn", all = TRUE)
binary_data <- merge(binary_data, binary_ovi_df, by = "sn", all = TRUE)

colnames(binary_data) <- c("Species", "Prop_Territorial", "Prop_Mate_Guard", "Prop_Flier_vs_Percher", "Prop_Oviposition")
#This worked! I checked manually.

#one more step: I need to remove NA values, but if I do that to the entire dataframe, I will have no rows left
#because almost no rows have no missing values

#So I make dataframes of each pair of traits:
mate_guard_terr_data_with_na<-merge(binary_terr_df, binary_mate_guard_df, by = "sn", all = TRUE)
fly_v_perch_terr_data_with_na<-merge(binary_terr_df, binary_fly_v_perch_df, by = "sn", all = TRUE)
ovi_terr_data_with_na<-merge(binary_terr_df, binary_ovi_df, by = "sn", all = TRUE)
#now remove NAs
mate_guard_terr_data_old<- mate_guard_terr_data_with_na[complete.cases(mate_guard_terr_data_with_na), ] 
fly_v_perch_terr_data_old<-fly_v_perch_terr_data_with_na[complete.cases(fly_v_perch_terr_data_with_na), ] 
ovi_terr_data_old<-ovi_terr_data_with_na[complete.cases(ovi_terr_data_with_na), ] 


#fitting the pagel (1994) model

#Mate guarding and territoriality
# Identify species to drop from mate_guard_terr_data
chk_mate_guard<-name.check(tree, mate_guard_terr_data_old, data.names=as.character(mate_guard_terr_data_old$sn))
summary(chk_mate_guard)
tree_mate_guard <- drop.tip(tree, chk_mate_guard$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
mate_guard_species_to_drop<-chk_mate_guard$data_not_tree
mate_guard_terr_data_old_dropped<-mate_guard_terr_data_old[!(mate_guard_terr_data_old$sn %in% mate_guard_species_to_drop),] #dropped data_not_tree species from dataset
name.check(tree_mate_guard, mate_guard_terr_data_old_dropped, data.names=as.character(mate_guard_terr_data_old_dropped$sn))
#these have to be in the right format
row_names_mate_guard <- mate_guard_terr_data_old_dropped$sn
mate_guard_terr_data<-data.frame(sp_binary_terr = ifelse(mate_guard_terr_data_old_dropped$sp_binary_terr == 1, "territorial", "non-territorial"),
                                 sp_binary_mate_guard = ifelse(mate_guard_terr_data_old_dropped$sp_binary_mate_guard == 1, "contact", "non-contact"))
rownames(mate_guard_terr_data) <- row_names_mate_guard

#run pagel 94 model
terr_mode_pagel_mate_guard<-setNames(mate_guard_terr_data[,1],
                                     rownames(mate_guard_terr_data))
mate_mode_pagel_mate_guard<-setNames(mate_guard_terr_data[,2],
                                     rownames(mate_guard_terr_data))
mate_guard_fit<-fitPagel(tree_mate_guard, terr_mode_pagel_mate_guard, mate_mode_pagel_mate_guard)
mate_guard_fit
#the dependent model has lower AIC! And a very small p-value. 
#plot this
plot(mate_guard_fit, signif=2, cex.main=1, cex.sub=0.8, cex.traits=0.7, cex.rates=0.7, lwd=1)

# Pagel 94 model for flier vs percher and territoriality
# Identify species to drop from fly_v_perch_terr_data
chk_fly_v_perch<-name.check(tree, fly_v_perch_terr_data_old, data.names=as.character(fly_v_perch_terr_data_old$sn))
summary(chk_fly_v_perch)
tree_fly_v_perch <- drop.tip(tree, chk_fly_v_perch$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
fly_v_perch_species_to_drop<-chk_fly_v_perch$data_not_tree
fly_v_perch_terr_data_old_dropped<-fly_v_perch_terr_data_old[!(fly_v_perch_terr_data_old$sn %in% fly_v_perch_species_to_drop),] #dropped data_not_tree species from dataset
name.check(tree_fly_v_perch, fly_v_perch_terr_data_old_dropped, data.names=as.character(fly_v_perch_terr_data_old_dropped$sn))
#these have to be in the right format:
row_names_fly_v_perch <- fly_v_perch_terr_data_old_dropped$sn
fly_v_perch_terr_data<-data.frame(sp_binary_terr = ifelse(fly_v_perch_terr_data_old_dropped$sp_binary_terr == 1, "territorial", "non-territorial"),
                                  sp_fly_v_perch = ifelse(fly_v_perch_terr_data_old_dropped$sp_fly_v_perch == 1, "percher", "flier"))
rownames(fly_v_perch_terr_data) <- row_names_fly_v_perch

#run pagel 94 model
terr_mode_pagel_fly_v_perch<-setNames(fly_v_perch_terr_data[,1],
                                      rownames(fly_v_perch_terr_data))
fly_mode_pagel_fly_v_perch<-setNames(fly_v_perch_terr_data[,2],
                                     rownames(fly_v_perch_terr_data))
fly_v_perch_fit<-fitPagel(tree_fly_v_perch, terr_mode_pagel_fly_v_perch, fly_mode_pagel_fly_v_perch)
fly_v_perch_fit
#independent model has lower AIC and p-value is insignificant. 
#plot this
plot(fly_v_perch_fit, signif=2, cex.main=1, cex.sub=0.8, cex.traits=0.7, cex.rates=0.7, lwd=1)

# Pagel 94 model for oviposition (endophytic vs exophytic) and territoriality
# Identify species to drop from ovi_terr_data
chk_ovi<-name.check(tree, ovi_terr_data_old, data.names=as.character(ovi_terr_data_old$sn))
summary(chk_ovi)
tree_ovi <- drop.tip(tree, chk_ovi$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
ovi_species_to_drop<-chk_ovi$data_not_tree
ovi_terr_data_old_dropped<-ovi_terr_data_old[!(ovi_terr_data_old$sn %in% ovi_species_to_drop),] #dropped data_not_tree species from dataset
name.check(tree_ovi,ovi_terr_data_old_dropped, data.names=as.character(ovi_terr_data_old_dropped$sn))
#these have to be in the right format:
row_names_ovi_terr <- ovi_terr_data_old_dropped$sn
ovi_terr_data<-data.frame(sp_binary_terr = ifelse(ovi_terr_data_old_dropped$sp_binary_terr == 1, "territorial", "non-territorial"),
                          sp_ovi = ifelse(ovi_terr_data_old_dropped$sp_ovi == 1, "exophytic", "endophytic"))
rownames(ovi_terr_data) <- row_names_ovi_terr

#run pagel 94 model
terr_mode_pagel_ovi<-setNames(ovi_terr_data[,1],
                              rownames(ovi_terr_data))
ovi_mode_pagel_ovi<-setNames(ovi_terr_data[,2],
                             rownames(ovi_terr_data))
ovi_fit<-fitPagel(tree_ovi, terr_mode_pagel_ovi, ovi_mode_pagel_ovi)
ovi_fit
#lower AIC is the independent model
#p-value is insignificant
#plot this
plot(ovi_fit, signif=2, cex.main=1, cex.sub=0.8, cex.traits=0.7, cex.rates=0.7, lwd=1)

#we can also plot the trees to visually display the data:
#not really working
#object<-plotTree.datamatrix(tree_mate_guard, mate_guard_terr_data, fsize=0.5, yexp=1, header=FALSE, xexp=1.45, palettes=c("YlOrRd", "PuBuGn"))
#leg<-legend(x="topright", names(object$colors$prop_binary_terr), cex=0.7, pch=22, pt.bg=object$colors$prop_binary_terr, pt.cex=1.5, bty="n", title="territory mode")
#second legend
#leg<- legend (x=leg$rect$left+4.7, y=leg$rect$top-leg$rect$h,
#             names(object$colors$prop_pate_guard), cex=0.7, pch=22, pt.bg=object$colors$prop_mate_guard, pt.cex=1.5, bty="n", title="mate guarding")
#could make a tree for all of them

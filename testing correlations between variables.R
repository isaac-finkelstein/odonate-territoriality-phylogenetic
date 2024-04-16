#Testing correlations between variables

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(dplyr)
library("phylolm")
library(ggplot2)

my_data<- read.csv("data/data_v5.csv") #this dataset (3rd version) switches "tandem" for "contact" in De Recende's data.
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


#now I stitch these together in a single dataframe
binary_data <- merge(binary_terr_df, binary_mate_guard_df, by = "sn", all = TRUE)
binary_data <- merge(binary_data, binary_fly_v_perch_df, by = "sn", all = TRUE)
binary_data<-merge(binary_data, binary_court_df, by="sn", all=TRUE)
binary_data <- merge(binary_data, binary_ovi_df, by = "sn", all = TRUE)


colnames(binary_data) <- c("Species", "Prop_Territorial", "Prop_Mate_Guard", "Prop_Flier_vs_Percher", "Prop_Oviposition", "Prop_Courtship")
#This worked! I checked manually.

#one more step: I need to remove NA values, but if I do that to the entire dataframe, I will have no rows left
#because almost no rows have no missing values

#So I make dataframes of each pair of traits:
mate_guard_terr_data_with_na<-merge(binary_terr_df, binary_mate_guard_df, by = "sn", all = TRUE)
fly_v_perch_terr_data_with_na<-merge(binary_terr_df, binary_fly_v_perch_df, by = "sn", all = TRUE)
ovi_terr_data_with_na<-merge(binary_terr_df, binary_ovi_df, by = "sn", all = TRUE)
court_terr_data_with_na<-merge(binary_terr_df, binary_court_df, by="sn", all=TRUE)
#now remove NAs
mate_guard_terr_data_old<- mate_guard_terr_data_with_na[complete.cases(mate_guard_terr_data_with_na), ] 
fly_v_perch_terr_data_old<-fly_v_perch_terr_data_with_na[complete.cases(fly_v_perch_terr_data_with_na), ]
court_terr_data_old<-court_terr_data_with_na[complete.cases(court_terr_data_with_na), ]
ovi_terr_data_old<-ovi_terr_data_with_na[complete.cases(ovi_terr_data_with_na), ] 



#fitting the pagel (1994) model

#Mate guarding and territoriality
# Identify species to drop from mate_guard_terr_data_old
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

#plot this -- how many instances of territorial/contact, territorial/non-contact etc. 
max_obs <- nrow(mate_guard_terr_data)
breaks <- seq(0, max_obs, by = 10)
ggplot(mate_guard_terr_data, aes(x = sp_binary_terr, fill = sp_binary_mate_guard)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Territorial", y = "Number of species", fill = "Mate guarding") +
  scale_fill_manual(values = c("contact" = "lightblue", "non-contact" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(breaks = breaks)
#let's check to make sure that's correct:
count_non_contact_territorial <- sum(mate_guard_terr_data$sp_binary_mate_guard == "non-contact" & mate_guard_terr_data$sp_binary_terr == "territorial")
#looks good


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
fly_v_perch_terr_data<-data.frame(sp_binary_terr = ifelse(fly_v_perch_terr_data_old_dropped$sp_binary_terr == 1, "Territorial", "Non-territorial"),
                                  sp_fly_v_perch = ifelse(fly_v_perch_terr_data_old_dropped$sp_fly_v_perch == 1, "Percher", "Flier"))
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

#plot this
max_obs_fly <- nrow(fly_v_perch_terr_data)
breaks_fly <- seq(0, max_obs_fly, by = 20)
ggplot(fly_v_perch_terr_data, aes(x = sp_binary_terr, fill = sp_fly_v_perch)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Territorial", y = "Number of species", fill = "Flier or percher") +
  scale_fill_manual(values = c("Flier" = "lightblue", "Percher" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(breaks = breaks_fly)

#pagel 94 model for courtship and territoriality
#identify species to drop
chk_court<-name.check(tree, court_terr_data_old, data.names=as.character(court_terr_data_old$sn))
summary(chk_court)
tree_court <- drop.tip(tree, chk_court$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
court_species_to_drop<-chk_court$data_not_tree
court_terr_data_old_dropped<-court_terr_data_old[!(court_terr_data_old$sn %in% court_species_to_drop),] #dropped data_not_tree species from dataset
name.check(tree_court, court_terr_data_old_dropped, data.names=as.character(court_terr_data_old_dropped$sn))
#these have to be in the right format:
row_names_court <- court_terr_data_old_dropped$sn
court_terr_data<-data.frame(sp_binary_terr = ifelse(court_terr_data_old_dropped$sp_binary_terr == 1, "Territorial", "Non-territorial"),
                                  sp_courtship = ifelse(court_terr_data_old_dropped$sp_courtship == 1, "Yes", "No"))
rownames(court_terr_data) <- row_names_court

#run pagel 94 model
terr_mode_pagel_court<-setNames(court_terr_data[,1],
                                      rownames(court_terr_data))
court_pagel_court<-setNames(court_terr_data[,2],
                                     rownames(court_terr_data))
court_fit<-fitPagel(tree_court, terr_mode_pagel_court, court_pagel_court)
court_fit
#independent model has lower AIC and p-value is insignificant. 
#plot this
plot(court_fit, signif=2, cex.main=1, cex.sub=0.8, cex.traits=0.7, cex.rates=0.7, lwd=1)

#plot this
max_obs_court <- nrow(court_terr_data)
breaks_court<- seq(0, max_obs_court, by = 5)
ggplot(court_terr_data, aes(x = sp_binary_terr, fill = sp_courtship)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Territorial", y = "Number of species", fill = "Courtship") +
  scale_fill_manual(values = c("Yes" = "lightblue", "No" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(breaks = breaks_court)

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
ovi_terr_data<-data.frame(sp_binary_terr = ifelse(ovi_terr_data_old_dropped$sp_binary_terr == 1, "Territorial", "Non-territorial"),
                          sp_ovi = ifelse(ovi_terr_data_old_dropped$sp_ovi == 1, "Exophytic", "Endophytic"))
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

#plot this:
max_obs_ovi <- nrow(ovi_terr_data)
breaks_ovi <- seq(0, max_obs_ovi, by = 10)
ggplot(ovi_terr_data, aes(x = sp_binary_terr, fill = sp_ovi)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Territorial", y = "Number of species", fill = "Oviposition") +
  scale_fill_manual(values = c("Endophytic" = "lightblue", "Exophytic" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(breaks = breaks_ovi)
#we can also plot the trees to visually display the data:
#not really working
#object<-plotTree.datamatrix(tree_mate_guard, mate_guard_terr_data, fsize=0.5, yexp=1, header=FALSE, xexp=1.45, palettes=c("YlOrRd", "PuBuGn"))
#leg<-legend(x="topright", names(object$colors$prop_binary_terr), cex=0.7, pch=22, pt.bg=object$colors$prop_binary_terr, pt.cex=1.5, bty="n", title="territory mode")
#second legend
#leg<- legend (x=leg$rect$left+4.7, y=leg$rect$top-leg$rect$h,
#             names(object$colors$prop_pate_guard), cex=0.7, pch=22, pt.bg=object$colors$prop_mate_guard, pt.cex=1.5, bty="n", title="mate guarding")
#could make a tree for all of them


#Testing variables that I expect will have a causal relationship: Lotic vs lentic breeding sites, size of breeding site (pond, lake, puddle), temporary vs permanent oviposition site
#lotic vs lentic and temporary vs permanent are binary variables, so do like above. For the size of breeding site, I probably have to do a logistic regression.


#lotic vs lentic
filtered_lo_len <- subset(my_data, Lotic.vs.lentic..breeding.habitat. %in% c("Lotic", "Lentic"))
lo_len_var<-ftable(filtered_lo_len$Formatted_species, filtered_lo_len$Lotic.vs.lentic..breeding.habitat.)
prop_lo_len<-round(lo_len_var[,2]/(lo_len_var[,1]+lo_len_var[,2]),2) #this is the proportion that is lotic
#the 2 at the end rounds to 2 decimal 
#set a 75% threshold = 3:1 threshold
sp_lo_len<-ifelse(prop_lo_len >= 0.75, 1, ifelse(prop_lo_len <=0.25, 0, NA))
sn<-attr(lo_len_var, "row.vars")[[1]]
binary_lo_len<-data.frame(sn, sp_lo_len, stringsAsFactors = TRUE)
#so, 1 = lotic, 0 = lentic


lo_len_terr_data_with_na<-merge(binary_terr_df, binary_lo_len, by = "sn", all = TRUE)
#now remove NAs
lo_len_terr_data_old<- lo_len_terr_data_with_na[complete.cases(lo_len_terr_data_with_na), ] 


#fitting the pagel (1994) model
#Lotic vs Lentic and territoriality
# Identify species to drop from lo_len_terr_data_old
chk_lo_len<-name.check(tree, lo_len_terr_data_old, data.names=as.character(lo_len_terr_data_old$sn))
summary(chk_lo_len)
tree_lo_len <- drop.tip(tree, chk_lo_len$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
lo_len_species_to_drop<-chk_lo_len$data_not_tree
lo_len_terr_data_old_dropped<-lo_len_terr_data_old[!(lo_len_terr_data_old$sn %in% lo_len_species_to_drop),] #dropped data_not_tree species from dataset
name.check(tree_lo_len, lo_len_terr_data_old_dropped, data.names=as.character(lo_len_terr_data_old_dropped$sn))
#these have to be in the right format
row_names_lo_len <- lo_len_terr_data_old_dropped$sn
lo_len_terr_data<-data.frame(sp_binary_terr = ifelse(lo_len_terr_data_old_dropped$sp_binary_terr == 1, "Territorial", "Non-territorial"),
                                 sp_lo_len = ifelse(lo_len_terr_data_old_dropped$sp_lo_len == 1, "Lotic", "Lentic"))
rownames(lo_len_terr_data) <- row_names_lo_len

#run pagel 94 model
terr_mode_pagel_lo_len<-setNames(lo_len_terr_data[,1],
                                     rownames(lo_len_terr_data))
lo_len_pagel_lo_len<-setNames(lo_len_terr_data[,2],
                                     rownames(lo_len_terr_data))
lo_len_fit<-fitPagel(tree_lo_len, terr_mode_pagel_lo_len, lo_len_pagel_lo_len)
lo_len_fit
#the independent model has lower AIC. And insignificant p-value. 
#plot this
plot(lo_len_fit, signif=2, cex.main=1, cex.sub=0.8, cex.traits=0.7, cex.rates=0.7, lwd=1)

#plot this
max_obs_lo_len <- nrow(lo_len_terr_data)
breaks_lo_len <- seq(0, max_obs_lo_len, by = 20)
ggplot(lo_len_terr_data, aes(x = sp_binary_terr, fill = sp_lo_len)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Territorial", y = "Number of species", fill = "Lotic vs lentic") +
  scale_fill_manual(values = c("Lotic" = "lightblue", "Lentic" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(breaks = breaks_lo_len)

#temporary vs permanent oviposition sites
#I use only data from Renner et al., 2020 to maintain consistency in categorization and observation opportunity
temp_v_perm_table<-ftable(my_data$Formatted_species, my_data$Temporary.vs.permanent.oviposition.site)

prop_generalist<-round(temp_v_perm_table[,2]/(temp_v_perm_table[,3]+temp_v_perm_table[,4]+temp_v_perm_table[,2]), 2)
prop_permanent<-round(temp_v_perm_table[,3]/(temp_v_perm_table[,3]+temp_v_perm_table[,4]+temp_v_perm_table[,2]), 2)
prop_temporary<-round(temp_v_perm_table[,4]/(temp_v_perm_table[,3]+temp_v_perm_table[,4]+temp_v_perm_table[,2]), 2)
prop_oviposition_temp_perm_gen<-data.frame(prop_generalist=prop_generalist,
                                           prop_permanent=prop_permanent,
                                           prop_temporary=prop_temporary)
prop_temp_perm_df<-data.frame(prop_generalist=ifelse(prop_generalist>0.75, 1, NA),
                              prop_permanent=ifelse(prop_permanent>0.75, 1, NA),
                              prop_temporary=ifelse(prop_temporary>0.75, 1, NA))

sn<-attr(temp_v_perm_table, "row.vars")[[1]]
temp_v_perm_data_na<-data.frame(sn, prop_temp_perm_df)
temp_v_perm_data<-temp_v_perm_data_na %>%
  filter(!is.na(prop_generalist | !is.na(prop_permanent) | !is.na(prop_temporary)))
oviposition_site_cat<-temp_v_perm_data %>%
  mutate(oviposition_cat = case_when(
    prop_generalist==1 ~"generalist",
    prop_permanent==1 ~"permanent",
    TRUE ~"temporary" #default to temporary if not generalist or permanent
  ))
oviposition_site_cat<-oviposition_site_cat %>%
  select(sn, oviposition_cat)

#make dataset
data_oviposition_site_old<-merge(binary_terr_df, oviposition_site_cat, by="sn", all=TRUE)
colnames(data_oviposition_site_old)<-c("Species", "Territorial", "Oviposition_cat")
data_oviposition_site_old<-data_oviposition_site_old[complete.cases(data_oviposition_site_old), ]

#species to drop
chk_ovi_site<-name.check(tree, data_oviposition_site_old, data.names=as.character(data_oviposition_site_old$Species))
summary(chk_ovi_site)
tree_ovi_site<-drop.tip(tree, chk_ovi_site$tree_not_data)
#drop species from data
ovi_site_sp_to_drop<-chk_ovi_site$data_not_tree
data_oviposition_site<-na.omit(data_oviposition_site_old[!(data_oviposition_site_old$Species %in% ovi_site_sp_to_drop),])
rownames(data_oviposition_site)<-data_oviposition_site$Species
name.check(tree_ovi_site, data_oviposition_site, data.names= as.character(data_oviposition_site$Species))
#sadly, this is only 13 observations
#and it is entirely territorial species so can't do a regression test
#test
#mod_ovi_site<-phyloglm(Territorial~Oviposition_cat, data=data_oviposition_site, phy=tree_ovi_site, boot=1000, method='logistic_MPLE', btol=10)


#Let's test if territorial species are more likely to have conservation concerns
#7 ordered categories: Least concern, Near threatened, Vulnerable, Endangered, Critically Endangered, Extinct
#unfortunately, this results in very little power. So, I need to combine everything (NT, V, E, CE) vs Least Concern (LC).
#There is no odonate species that is extinct in the wild
#remove "Data Deficient" and "No Data"
redlist_category <- subset(my_data, Redlist.category %in% c("Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Critically Endangered", "Extinct"))
redlist_cat<-ftable(redlist_category$Formatted_species, redlist_category$Redlist.category)

#For some reason the order of the columns is ordered alphabetically
lc<-round(redlist_cat[,4]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)
nt<-round(redlist_cat[,5]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)
vul<-round(redlist_cat[,6]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)
end<-round(redlist_cat[,2]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)
cr_end<-round(redlist_cat[,1]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)
ext<-round(redlist_cat[,3]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)

#here I make (c) concern (all categories except lc) versus (lc) least concern categories
redlist_df<-data.frame(
  lc=ifelse(lc>0.75,1,NA),
  c=ifelse(nt>0.75,1,NA),
  c=ifelse(vul>0.75,1,NA),
  c=ifelse(end>0.75,1,NA),
  c=ifelse(cr_end>0.75,1,NA),
  c=ifelse(ext>0.75,1,NA))

sn<-attr(redlist_cat, "row.vars")[[1]]
redlist_w_na<-data.frame(sn, redlist_df)
redlist<-redlist_w_na %>%
  filter(!is.na(lc) | !is.na(c))
sp_redlist_cat<-redlist %>%
  mutate(redlist_cat_final = case_when(
    lc==1 ~"least concern",
    c==1 ~"concern"))

sp_redlist_cat_ordered<-sp_redlist_cat %>%
  select(sn, redlist_cat_final)

#make dataset
data_redlist_terr_old<-merge(binary_terr_df, sp_redlist_cat_ordered, by="sn", all=TRUE)
colnames(data_redlist_terr_old)<-c("Species", "Territorial", "Redlist_category")
data_redlist_terr_old<-data_redlist_terr_old[complete.cases(data_redlist_terr_old), ]

#species to drop
chk_redlist<-name.check(tree, data_redlist_terr_old, data.names=as.character(data_redlist_terr_old$Species))
summary(chk_redlist)
tree_redlist<-drop.tip(tree, chk_redlist$tree_not_data) #dropped tree_not_species
#drop from data
drop_redlist<-chk_redlist$data_not_tree
data_redlist_terr<-na.omit(data_redlist_terr_old[!(data_redlist_terr_old$Species %in% drop_redlist),])
rownames(data_redlist_terr)<-data_redlist_terr$Species
name.check(tree_redlist, data_redlist_terr, data.names=as.character(data_redlist_terr$Species))
data_redlist_terr<-data_redlist_terr[complete.cases(data_redlist_terr), ]


#these are ordered variables from most to least concern
mod_redlist<-phyloglm(Territorial~Redlist_category, data= data_redlist_terr, phy=tree_redlist, btol=20, boot=1000, method='logistic_MPLE')


#plot
max_obs_iucn<-nrow(data_redlist_terr)
breaks<-seq(0, max_obs_iucn, by =5)
ggplot(data_redlist_terr, aes(x = Territorial, fill = Redlist_category)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Territorial", y = "Number of species", fill = "Redlist category") +
  scale_fill_manual(values = c("least concern" = "lightblue", "concern" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(breaks = breaks)

#Finally, we can test the size of the water body used as breeding habitat

# For this I am doing a phylgenetic anova with 5 categories:
#lotic+small, lotic+medium, lotic+big, lentic+small, lentic+big

#first Lentic:
#create an ordered categorical variable
#3 categories: small, medium, large
#specify that these are ordinal
lentic_size_table<- ftable(my_data$Formatted_species, my_data$Lentic.size)
prop_large<-round(lentic_size_table[,2]/(lentic_size_table[,3]+ lentic_size_table[,4]+ lentic_size_table[,2]),2)
prop_medium<-round(lentic_size_table[,3]/(lentic_size_table[,2]+ lentic_size_table[,3]+ lentic_size_table[,4]),2)
prop_small<-round(lentic_size_table[,4]/(lentic_size_table[,2]+ lentic_size_table[,3]+ lentic_size_table[,4]),2)
prop_data<-data.frame(prop_large=prop_large,
                      prop_medium=prop_medium,
                      prop_small=prop_small)
prop_data <- data.frame(
  prop_large = ifelse(prop_large > 0.75, 1, NA),
  prop_medium = ifelse(prop_medium > 0.75, 1, NA),
  prop_small = ifelse(prop_small > 0.75, 1, NA)
) #this makes a 3:1 cutoff - so if my data disagrees, I only include it if there is at least 3:1 sources saying the size - this is rarely met so basically I exclude species where there is disagreement
sn<- attr(lentic_size_table,"row.vars")[[1]]
size_data_with_na<-data.frame(sn,prop_data)
size_data <- size_data_with_na %>%
  filter(!is.na(prop_large) | !is.na(prop_medium) | !is.na(prop_small))

lentic_size_data <- size_data %>%
  mutate(lentic_size = case_when(
    prop_large == 1 ~ "lentic_large",
    prop_medium == 1 ~ "lentic_medium",
    TRUE ~ "lentic_small"  # Default to "small" if neither prop_large nor prop_medium is 1
  ))%>%
  select(-prop_large, -prop_medium, -prop_small)
#size_data_ordered <- size_data_ordered %>%
#  select(sn, lentic_size)
#size_data_ordered$lentic_size <- factor(size_data_ordered$lentic_size, ordered = TRUE, levels = c("lentic_small", "lentic_medium", "lentic_large")) # make it an ordered categorical variable

#make dataset
data_lentic_size_territoriality_old <- merge(binary_terr_df, lentic_size_data, by = "sn", all = TRUE)
colnames(data_lentic_size_territoriality_old) <- c("Species", "Prop_Territorial", "Lentic_size")
data_lentic_size_territoriality_old<- data_lentic_size_territoriality_old[complete.cases(data_lentic_size_territoriality_old), ] 


#lotic size
my_data$Description.of.lotic.oviposition..river..stream. <- gsub("Stream, River", "River, Stream", my_data$Description.of.lotic.oviposition..river..stream.)
#3 categories: Stream, Both (River,Stream), and River
#also this is ordinal
lotic_size_table<-ftable(my_data$Formatted_species, my_data$Description.of.lotic.oviposition..river..stream.)

prop_stream<-round(lotic_size_table[,4]/(lotic_size_table[,3]+lotic_size_table[,4]+ lotic_size_table[,2]), 2)
prop_both<-round(lotic_size_table[,3]/(lotic_size_table[,3]+lotic_size_table[,4]+ lotic_size_table[,2]), 2)
prop_river<-round(lotic_size_table[,2]/(lotic_size_table[,3]+lotic_size_table[,4]+ lotic_size_table[,2]), 2)
prop_lotic_data<-data.frame(prop_stream=prop_stream,
                            prop_both=prop_both,
                            prop_river=prop_river)
prop_lotic_data<-data.frame(
  prop_stream=ifelse(prop_stream>0.75,1,NA),
  prop_both=ifelse(prop_both>0.75,1, NA),
  prop_river=ifelse(prop_river>0.75,1, NA))
#this makes a 3:1 cutoff - so if my data disagrees, I only include it if there is at least 3:1 sources saying the size - this is rarely met so basically I exclude species where there is disagreement
sn<- attr(lotic_size_table,"row.vars")[[1]]
lotic_data_with_na<-data.frame(sn,prop_lotic_data)
lotic_size_data <- lotic_data_with_na %>%
  filter(!is.na(prop_stream) | !is.na(prop_both) | !is.na(prop_river))

lotic_size<-lotic_size_data %>%
  mutate(lotic_size = case_when(
    prop_stream == 1 ~"stream",
    prop_both == 1 ~ "both",
    TRUE ~ "river" #default to river if neither prop_stream or prop_both
  ))%>%
  select(-prop_stream, -prop_both, -prop_river)
#lotic_size_ordered<-lotic_size_ordered %>%
#  select(sn, lotic_size)
#lotic_size_ordered$lotic_size<- factor(lotic_size_ordered$lotic_size, ordered = TRUE, levels = c("stream", "both", "river")) #make it ordered categorical

combined_data <- full_join(data_lentic_size_territoriality_old, lotic_size, by = c("Species" = "sn"))
combined_data$lentic_lotic_size <- coalesce(combined_data$Lentic_size, combined_data$lotic_size)
combined_data <- combined_data[, c("Species", "Prop_Territorial", "lentic_lotic_size")]
#remove three species that oviposit in both lentic and lotic locations (these would make new groups but they are not common so just removing)
ovi_size_data_old <- combined_data %>%
  filter(Species != "Orthetrum_cancellatum" & 
           Species != "Cordulia_aenea" & 
           Species != "Lestes_viridis")
ovi_size_data_old<-na.omit(ovi_size_data_old)


#now let's run an phylogenetic anova
# Identify species to drop
chk_ovi_size<-name.check(tree, ovi_size_data_old, data.names=as.character(ovi_size_data_old$Species))
summary(chk_ovi_size)
tree_ovi_size <- drop.tip(tree, chk_ovi_size$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
ovi_size_species_to_drop<-chk_ovi_size$data_not_tree
ovi_size_data<-ovi_size_data_old[!(ovi_size_data_old$Species %in% ovi_size_species_to_drop),] #dropped data_not_tree species from dataset
rownames(ovi_size_data)<-ovi_size_data$Species
name.check(tree_ovi_size, ovi_size_data, data.names=as.character(ovi_size_data$Species))

library(nlme)
spp<- rownames(ovi_size_data)
corBM<-corBrownian(phy=tree_ovi_size, form=~spp)
corBM
anova<-gls(Prop_Territorial~lentic_lotic_size, data=ovi_size_data, correlation=corBM)
anova(anova)
#Did I do this right?


#plot this
ovi_size_data <- ovi_size_data %>%
  mutate(Prop_Territorial = ifelse(Prop_Territorial == 1, "Territorial", "Non-territorial"))

ggplot(ovi_size_data, aes(x = Prop_Territorial, fill = lentic_lotic_size)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = after_stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Territorial", y = "Number of species", fill = "Oviposition size") +
  scale_fill_manual(values = c("lentic_small" ="darkblue", "lentic_medium" = "darkorange", "lentic_large" = "darkred", "stream"="lightblue", "both" = "purple", river="turquoise")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))


#Including a "no" category for mate guarding

#this results in 52 non-contact data instead of 54 when it was binary
#when I coded the binary variable, I had a line similar to this:
#my_data_filtered<-subset(my_data, Mate.guarding %in% c("No", "Contact", "Non-contact"))
#however, when I run that line instead, I get an error and I can't figure out why. 
my_data$Mate.guarding <- trimws(my_data$Mate.guarding)
my_data_filtered <- my_data[!(my_data$Mate.guarding == "Both"), ] #There's only 1 "Both" so I'm removing it.
mate_guarding_no_table<-ftable(my_data_filtered$Formatted_species, my_data_filtered$Mate.guarding)

prop_no<-round(mate_guarding_no_table[,3]/(mate_guarding_no_table[,3]+mate_guarding_no_table[,4]+mate_guarding_no_table[,2]), 2)
prop_contact_no<-round(mate_guarding_no_table[,2]/(mate_guarding_no_table[,3]+mate_guarding_no_table[,4]+mate_guarding_no_table[,2]), 2)
prop_non_contact_no<-round(mate_guarding_no_table[,4]/(mate_guarding_no_table[,3]+mate_guarding_no_table[,4]+mate_guarding_no_table[,2]), 2)
prop_mate_guard_no_df<-data.frame(prop_no=prop_no,
                                  prop_contact_no=prop_contact_no,
                                  prop_non_contact_no=prop_non_contact_no)
prop_mate_guard_no_df<-data.frame(
  prop_no=ifelse(prop_no>0.75,1,NA),
  prop_contact_no=ifelse(prop_contact_no>0.75,1,NA),
  prop_non_contact_no=ifelse(prop_non_contact_no>0.75, 1, NA))
#this makes a 3:1 cutoff - so if my data disagrees, I only include it if there is at least 3:1 sources saying the size - this is rarely met so basically I exclude species where there is disagreement

sn<-attr(mate_guarding_no_table, "row.vars")[[1]]
mate_guard_no_data_na<-data.frame(sn, prop_mate_guard_no_df)
mate_guard_no_data<-mate_guard_no_data_na %>%
  filter(!is.na(prop_no) | !is.na(prop_contact_no) | !is.na(prop_non_contact_no))
mate_guarding_three_cat<-mate_guard_no_data %>%
  mutate(Mate_guarding_cat = case_when(
    prop_no==1 ~"No",
  prop_contact_no == 1 ~"Contact",
  TRUE ~"Non-contact" #default to non-contact if neither contact or no
  ))
mate_guarding_three_cat<-mate_guarding_three_cat %>%
  select(sn, Mate_guarding_cat)
#make dataset
data_mate_guard_no_terr_old<-merge(binary_terr_df, mate_guarding_three_cat, by="sn", all=TRUE)
colnames(data_mate_guard_no_terr_old)<-c("Species", "Prop_territorial", "Mate_guarding_cat")
data_mate_guard_no_terr_old<-data_mate_guard_no_terr_old[complete.cases(data_mate_guard_no_terr_old), ]

#identify species to drop
chk_mate_guard_no<-name.check(tree, data_mate_guard_no_terr_old, data.names=as.character(data_mate_guard_no_terr_old$Species))
summary(chk_mate_guard_no)
tree_mate_guard_no<-drop.tip(tree, chk_mate_guard_no$tree_not_data) #dropped tree_not_data species
#identify species to drop from data
guard_no_species_to_drop<-chk_mate_guard_no$data_not_tree
data_mate_guard_no_terr<-na.omit(data_mate_guard_no_terr_old[!(data_mate_guard_no_terr_old$Species %in% guard_no_species_to_drop),]) #dropped data_not_tree species
rownames(data_mate_guard_no_terr)<-data_mate_guard_no_terr$Species
name.check(tree_mate_guard_no, data_mate_guard_no_terr, data.names = as.character(data_mate_guard_no_terr$Species))

#test if mate guarding with 3 categories (includes "no") predicts territoriality
mod_mate_guard<-phyloglm(Prop_territorial~Mate_guarding_cat, data=data_mate_guard_no_terr, phy=tree_mate_guard_no, boot=1000, method='logistic_MPLE', btol=10)
summary(mod_mate_guard)

#plot this
custom_levels <- c("No", "Contact", "Non-contact")
ggplot(data_mate_guard_no_terr, aes(x = factor(Prop_territorial), fill = factor(Mate_guarding_cat, levels = custom_levels))) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Territorial", y = "Number of species", fill = "Mate guarding") +
  scale_x_discrete(labels = c("0" = "Non-Territorial", "1" = "Territorial")) +
  scale_fill_manual(values = c("No" = "darkblue", "Contact" = "darkorange", "Non-contact" = "darkred")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

# Calculate counts of observations for each level of "Mate_guarding_cat"
count_data <- data_mate_guard_no_terr %>%
  group_by(Mate_guarding_cat) %>%
  summarise(count = n())

# Plot
ggplot(data_mate_guard_no_terr, aes(x = factor(Prop_territorial), fill = Mate_guarding_cat)) +
  geom_bar(position = "dodge") +
  geom_text(data = count_data, aes(label = count), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Territorial", y = "Number of species", fill = "Mate guarding") +
  scale_x_discrete(labels = c("0" = "Non-Territorial", "1" = "Territorial")) +
  scale_fill_manual(values = c("no" ="darkblue", "contact" = "darkorange", "non-contact" = "darkred")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))


#Doing a regression of all my variables
#Make a single dataframe
big_data <- merge(binary_terr_df, mate_guarding_three_cat, by = "sn", all = TRUE)
big_data <- merge(big_data, binary_fly_v_perch_df, by = "sn", all = TRUE)
big_data <- merge(big_data, binary_ovi_df, by = "sn", all = TRUE)

colnames(big_data) <- c("Species", "territorial", "mate_guard", "flier_vs_percher", "oviposition")
big_data<-big_data[complete.cases(big_data), ]
#identify species to drop
chk_big_data<-name.check(tree, big_data, data.names=as.character(big_data$Species))
summary(chk_big_data)
tree_big_data<-drop.tip(tree, chk_big_data$tree_not_data) #dropped tree_not_data species
#identify species to drop from data
big_data_species_to_drop<-chk_big_data$data_not_tree
data_big_new<-na.omit(big_data[!(big_data$Species %in% big_data_species_to_drop),]) #dropped data_not_tree species
rownames(data_big_new)<-data_big_new$Species
name.check(tree_big_data, data_big_new, data.names = as.character(data_big_new$Species))

mod_all_var<-phyloglm(territorial~mate_guard + flier_vs_percher + oviposition, data=data_big_new, phy=tree_big_data, boot=1000, method='logistic_MPLE', btol=10)
summary(mod_all_var)
#but small sample sizes:
View(data_big_new)

#Plot each variable:
#mate guard
ggplot(data_big_new, aes(x = territorial, fill = mate_guard)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 2) +
  labs(x = "Territorial", y = "Count", fill = "Mate guarding") +
  scale_fill_manual(values = c("no" = "gray", "contact" = "lightblue", "non-contact" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
#not sure why the territorial/non-territorial labels on the x-axis are not showing

data_big_new$flier_vs_percher <- factor(data_big_new$flier_vs_percher, levels = c(0, 1), labels = c("flier", "percher"))

# Plot fly_v_perch
ggplot(data_big_new, aes(x = factor(territorial), fill = flier_vs_percher)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 2) +
  labs(x = "Territorial", y = "Count", fill = "Flier or percher") +
  scale_fill_manual(values = c("flier" = "lightblue", "percher" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

#oviposition
data_big_new$oviposition <- factor(data_big_new$oviposition, levels = c(0, 1), labels = c("endophytic", "exophytic"))

# Plot fly_v_perch
ggplot(data_big_new, aes(x = factor(territorial), fill = oviposition)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 2) +
  labs(x = "Territorial", y = "Count", fill = "oviposition") +
  scale_fill_manual(values = c("endophytic" = "lightblue", "exophytic" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))



#Testing if abundance predicts territoriality
#Using data from De Resende et al., 2021 for abundance data
abundance_with_na<-my_data$Abundance
sn<-my_data$Formatted_species
abundance <- na.omit(abundance_with_na)
abundance_df <- data.frame(sn = sn[!is.na(abundance_with_na)], Abundance = abundance)


#make dataset
data_abundance_terr_old<-merge(binary_terr_df, abundance_df, by = "sn", all=TRUE)
colnames(data_abundance_terr_old)<- c("Species", "Territorial", "Abundance")
data_abundance_terr_old<-data_abundance_terr_old[complete.cases(data_abundance_terr_old), ]


#identify species to drop
chk_abun<-name.check(tree, data_abundance_terr_old, data.names=as.character(data_abundance_terr_old$Species))
summary(chk_abun)
tree_abundance<-drop.tip(tree, chk_abun$tree_not_data)
#identify species to drop from data
abund_species_to_drop<-chk_abun$data_not_tree
data_abundance_terr<-data_abundance_terr_old[!(data_abundance_terr_old$Species %in% abund_species_to_drop), ]
rownames(data_abundance_terr)<-data_abundance_terr$Species
name.check(tree_abundance, data_abundance_terr, data.names=as.character(data_abundance_terr$Species))
#sadly, this only leaves 20 species

#Test if abundance predicts territoriality
mod_abund<-phyloglm(Territorial~Abundance, data=data_abundance_terr, phy=tree_abundance, boot=1000, method= 'logistic_MPLE', btol=10)
summary(mod_abund)

#plot this
#plot this
#unfortunately I have to rework the dataset so that territoriality is coded as "territorial"/"non-territorial", instead of 0/1.
#there is probably a faster way to do this but I just make the dataset again. 
terr_df_len<-data.frame(
  sn=binary_terr_df$sn,
  binary_terr_df = ifelse(binary_terr_df$sp_binary_terr == 1, "Territorial (n=60)", "Non-territorial (n=35)"),
  stringsAsFactors = TRUE)
#make dataset
abundance_terr_with_na <- merge(terr_df_len, abundance_df, by = "sn", all = TRUE)
colnames(abundance_terr) <- c("Species", "Territorial", "Abundance")
abundance_terr<- abundance_terr_with_na[complete.cases(abundance_terr_with_na), ]

colnames(abundance_terr) <- c("sn", "Territorial", "Abundance")

ggplot(abundance_terr, aes(x = Territorial, y = Abundance, fill = Territorial)) +
  geom_boxplot() +
  labs(x = "Territorial", y = "Abundance") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(breaks = seq(0, max(abundance_terr$Abundance), by = 25))

#means are very close, but territorial species seem a bit more inclined to be more abundance
#could also be that they tend to be easier to spot?




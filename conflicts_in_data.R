#finding conflicts in the data
library(phytools)
library(dplyr)
library(geiger)

my_data<- read.csv("data/data_v4.csv") #this dataset (3rd version) switches "tandem" for "contact" in De Recende's data.
#For some reason, they use both terms. But since they mean the same thing, I changed them all to "Contact"
#I also emailed to ask, but did not get a response. 
load(file="data/Odo.tree.Waller.Svensson.2017.rda") #odonate tree extracated from Waller and Svensson 2017
#str(tree) #plot(tree, no.margin=TRUE) #tree$Nnode #tree$tip.label -- to check the structure of the tree

#tree from Rocha-Ortega et al., 2020 in Proceedings
#it uses the same structure as Waller and Svensson 2017 but has different species coverage.
rocha_tree<-read.nexus(file="data/Rocha_ortega_tree.nex")
#I compare the trees starting at line 167

#total number of species in dataset
total_species <- length(unique(my_data$Formatted_species))
#total species = 1276

#total number of species with territorial data in dataset
terr_table<- ftable(my_data$Formatted_species, my_data$Territorial) #this has 3 columns, the first column is for NAs
prop_terr<- round(terr_table[,3]/(terr_table[,2]+ terr_table[,3]),2) # so this calculates the percentage "yes" for territorial
sn<- attr(terr_table,"row.vars")[[1]]
terr_data_with_na<- data.frame(sn,prop_terr) #so this is a dataframe with territorial (1/0) for each species
terr_data<- terr_data_with_na[complete.cases(terr_data_with_na), ] #removed NA values 
str(terr_data)
#total number of species = 643



#Put data into a readable dataframe -- one species per row
terr_table<- ftable(my_data$Formatted_species, my_data$Territorial) #this has 3 columns, the first column is for NAs
prop_terr<- round(terr_table[,3]/(terr_table[,2]+ terr_table[,3]),2) # so this calculates the percentage "yes" for territorial


#find conflicts in territoriality
conflict_terr<-ifelse(prop_terr <1 & prop_terr >0,2,NA) #any conlficts will be denoted as "2"
conflicts_sn<-attr(terr_table, "row.vars")[[1]]
conflict_terr<-data.frame(conflicts_sn, conflict_terr)
colnames(conflict_terr)<-c("Species", "Conflict")
conflicts_territoriality<-conflict_terr[conflict_terr$Conflict == 2, ]
conflicts_territoriality<-na.omit(conflicts_territoriality)
str(conflicts_territoriality) #so there are 36 conflicts
#If you check this with my dataset, you will notice that most involve from Waller et al., 2019. 



#conflicts in mate guarding data
filtered_mate_guarding <- subset(my_data, Mate.guarding %in% c("Contact", "Non-contact"))
mate_guard_var<-ftable(filtered_mate_guarding$Formatted_species, filtered_mate_guarding$Mate.guarding)
prop_mate_guard<-round(mate_guard_var[,1]/(mate_guard_var[,1]+mate_guard_var[,2]),2) #this is the proportion that is contact


#find conflicts
conflict_mate_guard<-ifelse(prop_mate_guard<1 & prop_mate_guard>0, 2, NA) #conflicts set to 2
conflicts_sn_mate_guard<-attr(mate_guard_var, "row.vars")[[1]]
conflict_mate_guard<-data.frame(conflicts_sn_mate_guard, conflict_mate_guard)
colnames(conflict_mate_guard)<-c("Species", "Conflict")
conflicts_mate_guarding<-conflict_mate_guard[conflict_mate_guard$Conflict ==2,]
conflicts_mate_guarding<-na.omit(conflicts_mate_guarding)
str(conflicts_mate_guarding) #so there are 5 conflicts




#conflicts in flier vs percher data:
binary_fly_v_perch<-ftable(my_data$Formatted_species, my_data$Flier.vs.percher)
prop_fly_v_perch<-round(binary_fly_v_perch[,3]/(binary_fly_v_perch[,2]+binary_fly_v_perch[,3]),2) #this is proportion Percher


#find conflicts
conflict_fly_v_perch<-ifelse(prop_fly_v_perch<1 & prop_fly_v_perch>0, 2, NA) #conflicts set to 2
conflicts_sn_fly_v_perch<-attr(binary_fly_v_perch, "row.vars")[[1]]
conflict_fly_v_perch<-data.frame(conflicts_sn_fly_v_perch, conflict_fly_v_perch)
colnames(conflict_fly_v_perch)<-c("Species", "Conflict")
conflicts_flying_vs_perching<-conflict_fly_v_perch[conflict_fly_v_perch$Conflict==2,]
conflicts_flying_vs_perching<-na.omit(conflicts_flying_vs_perching)
str(conflicts_flying_vs_perching) #8 conflicts



#conflicts in oviposition (endophytic vs exophytic)
my_data_mutated <- my_data %>%
  mutate(Oviposition.type..endophytic.vs.exophytic. = ifelse(Oviposition.type..endophytic.vs.exophytic. %in% c("Epiphytic", "Exophytic"), "Exophytic", Oviposition.type..endophytic.vs.exophytic.))
# Convert the variable to a factor with specified levels
my_data_mutated$Oviposition.type..endophytic.vs.exophytic. <- factor(
  my_data_mutated$Oviposition.type..endophytic.vs.exophytic.,
  levels = c("Endophytic", "Exophytic")
)
binary_ovi<-ftable(my_data_mutated$Formatted_species, my_data_mutated$Oviposition.type..endophytic.vs.exophytic.)
prop_ovi<-round(binary_ovi[,2]/(binary_ovi[,1]+binary_ovi[,2]),2) #This is proportion exophytic


#find conflicts
conflict_ovi<- ifelse(prop_ovi<1 & prop_ovi>0, 2, NA) #conflicts set to 2
conflicts_sn_ovi<-attr(binary_ovi, "row.vars")[[1]]
conflict_ovi<-data.frame(conflicts_sn_ovi, conflict_ovi)
colnames(conflict_ovi)<-c("Species", "Conflict")
conflicts_oviposition<-conflict_ovi[conflict_ovi$Conflict==2,]
conflicts_oviposition<-na.omit(conflicts_oviposition)
str(conflicts_oviposition) #4 conflicts



#conflicts in lotic vs lentic breeding habitat
filtered_lo_len <- subset(my_data, Lotic.vs.lentic..breeding.habitat. %in% c("Lotic", "Lentic"))
lo_len_var<-ftable(filtered_lo_len$Formatted_species, filtered_lo_len$Lotic.vs.lentic..breeding.habitat.)
prop_lo_len<-round(lo_len_var[,2]/(lo_len_var[,1]+lo_len_var[,2]),2) #this is the proportion that is lotic


#find conflicts
conflict_lo_v_len<- ifelse(prop_lo_len<1 & prop_lo_len>0, 2, NA) #conflict set to 2
conflict_sn_lo_len<- attr(lo_len_var, "row.vars")[[1]]
conflict_lo_v_len<-data.frame(conflict_sn_lo_len, conflict_lo_v_len)
colnames(conflict_lo_v_len)<-c("Species", "Conflict")
conflicts_lotic_vs_lentic<-conflict_lo_v_len[conflict_lo_v_len$Conflict==2,]
conflicts_lotic_vs_lentic<-na.omit(conflicts_lotic_vs_lentic)
str(conflicts_lotic_vs_lentic) # 41 conflicts




#conflicts in lentic size
lentic_size_table<- ftable(my_data$Formatted_species, my_data$Lentic.size)
prop_large<-round(lentic_size_table[,2]/(lentic_size_table[,3]+ lentic_size_table[,4]+ lentic_size_table[,2]),2)
prop_medium<-round(lentic_size_table[,3]/(lentic_size_table[,2]+ lentic_size_table[,3]+ lentic_size_table[,4]),2)
prop_small<-round(lentic_size_table[,4]/(lentic_size_table[,2]+ lentic_size_table[,3]+ lentic_size_table[,4]),2)
prop_data_original<-data.frame(prop_large=prop_large,
                      prop_medium=prop_medium,
                      prop_small=prop_small)

#find conflicts
conflict_len_size<-ifelse(prop_data_original<1 & prop_data_original>0, 5, NA) #conflicts set to 5
conflict_sn_len_size<-attr(lentic_size_table, "row.vars") [[1]]
conflict_len_size<-data.frame(conflict_sn_len_size, conflict_len_size)
colnames(conflict_len_size)<-c("Species", "Conflict")
conflict_lentic_size<-conflict_len_size[conflict_len_size$Conflict==5,]
conflict_lentic_size<-na.omit(conflict_lentic_size)
str(conflict_lentic_size) #no conflicts. 


#conflicts in lotic size:
my_data$Description.of.lotic.oviposition..river..stream. <- gsub("Stream, River", "River, Stream", my_data$Description.of.lotic.oviposition..river..stream.)
#3 categories: Stream, Both (River,Stream), and River
#also this is ordinal
lotic_size_table<-ftable(my_data$Formatted_species, my_data$Description.of.lotic.oviposition..river..stream.)

prop_stream<-round(lotic_size_table[,4]/(lotic_size_table[,3]+lotic_size_table[,4]+ lotic_size_table[,2]), 2)
prop_both<-round(lotic_size_table[,3]/(lotic_size_table[,3]+lotic_size_table[,4]+ lotic_size_table[,2]), 2)
prop_river<-round(lotic_size_table[,2]/(lotic_size_table[,3]+lotic_size_table[,4]+ lotic_size_table[,2]), 2)
prop_lotic_data_original<-data.frame(prop_stream=prop_stream,
                            prop_both=prop_both,
                            prop_river=prop_river)

#find conflicts
conflict_lot_size<-ifelse(prop_lotic_data_original<1 & prop_lotic_data_original>0, 5, NA) #conflicts set to 5
conflict_sn_lot_size<-attr(lotic_size_table, "row.vars")[[1]]
conflict_lot_size<-data.frame(conflict_sn_lot_size, conflict_lot_size)
colnames(conflict_lot_size)<-c("Species", "Conflict")
conflict_lotic_size<-conflict_lot_size[conflict_lot_size$Conflict==5,]
conflict_lotic_size<-na.omit(conflict_lotic_size)
str(conflict_lotic_size) #1 conflict




#comparing trees

#start with Waller and Svensson tree
#prune tree to match data
chk<-name.check(tree, terr_data, data.names=as.character(terr_data$sn))
summary(chk)
odonate_tree<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree
# Identify species to drop from terr_data
species_to_drop <- chk$data_not_tree
odonate_terr_data <- terr_data[!(terr_data$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_terr_data, data.names=as.character(odonate_terr_data$sn))

#let's plot the entire tree - busy
plot(odonate_tree, type="fan", cex=0.5, fsize=0.6, ftype="i") #it's very busy because it has 372 tips
#let's plot territoriality on my tree


#plot territoriality on tree
#to make the territorial variable binary (yes or no) make it a factor
prop_terr_factor <- factor(ifelse(prop_terr == 1, "yes", "no"))
terr_data_with_na_factor <- data.frame(sn, prop_terr_factor = prop_terr_factor) # Create a data frame
terr_data_factor <- terr_data_with_na_factor[complete.cases(terr_data_with_na_factor), ] #Removes NA values

odonate_tree_factor<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree_factor
species_to_drop <- chk$data_not_tree
odonate_terr_data_factor <- terr_data_factor[!(terr_data_factor$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_terr_data_factor, data.names=as.character(odonate_terr_data_factor$sn))


terr_mode<-setNames(odonate_terr_data_factor$prop_terr, odonate_terr_data_factor$sn)
cols<-setNames(c("forestgreen", "goldenrod1"), levels(terr_mode))
plotTree.datamatrix(odonate_tree_factor, as.data.frame(terr_mode),
                    colours= list(cols), header=FALSE)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)


#Now let's do the same with the Rocha-Ortega tree
#prune tree to match data
chk_rocha<-name.check(rocha_tree, terr_data, data.names=as.character(terr_data$sn))
summary(chk_rocha)
odonate_tree_rocha<-drop.tip(rocha_tree, chk_rocha$tree_not_data) #dropped tree_not_data species
odonate_tree_rocha
# Identify species to drop from terr_data
species_to_drop_rocha <- chk_rocha$data_not_tree
odonate_terr_data_rocha <- terr_data[!(terr_data$sn %in% species_to_drop_rocha), ] #dropped data_not_tree species from dataset
name.check(odonate_tree_rocha, odonate_terr_data_rocha, data.names=as.character(odonate_terr_data_rocha$sn))

#let's plot the entire tree - busy
plot(odonate_tree_rocha, type="fan", cex=0.5, fsize=0.6, ftype="i") #it's very busy because it has 372 tips
#let's plot territoriality on my tree


#plot territoriality on tree
#to make the territorial variable binary (yes or no) make it a factor
prop_terr_factor <- factor(ifelse(prop_terr == 1, "yes", "no"))
terr_data_with_na_factor <- data.frame(sn, prop_terr_factor = prop_terr_factor) # Create a data frame
terr_data_factor <- terr_data_with_na_factor[complete.cases(terr_data_with_na_factor), ] #Removes NA values

odonate_tree_factor_rocha<-drop.tip(rocha_tree, chk_rocha$tree_not_data) #dropped tree_not_data species
odonate_tree_factor_rocha
species_to_drop_rocha <- chk_rocha$data_not_tree
odonate_terr_data_factor_rocha <- terr_data_factor[!(terr_data_factor$sn %in% species_to_drop_rocha), ] #dropped data_not_tree species from dataset
name.check(odonate_tree_rocha, odonate_terr_data_factor_rocha, data.names=as.character(odonate_terr_data_factor_rocha$sn))


terr_mode<-setNames(odonate_terr_data_factor_rocha$prop_terr, odonate_terr_data_factor_rocha$sn)
cols<-setNames(c("forestgreen", "goldenrod1"), levels(terr_mode))
plotTree.datamatrix(odonate_tree_factor_rocha, as.data.frame(terr_mode),
                    colours= list(cols), header=FALSE)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)

str(odonate_terr_data_factor_rocha)
#vs
str(odonate_terr_data_factor)
#much fewer species on the Rocha_tree

#what species differ?
species_rocha <- odonate_terr_data_factor_rocha$sn
species_sven <- odonate_terr_data_factor$sn

#unique species to each dataset
species_unique_to_rocha <- setdiff(species_rocha, species_sven)
species_unique_to_sven <- setdiff(species_sven, species_rocha)

str(species_unique_to_rocha)
str(species_unique_to_sven)

#I can do the same thing for the specific traits wherer the results differ when using the different trees


#finding conflicts in the data
library(tidyverse)


my_data<- read.csv("data/data_v4.csv") #this dataset (3rd version) switches "tandem" for "contact" in De Recende's data.
#For some reason, they use both terms. But since they mean the same thing, I changed them all to "Contact"
#I also emailed to ask, but did not get a response. 
load(file="data/Odo.tree.Waller.Svensson.2017.rda") #odonate tree extracated from Waller and Svensson 2017
#str(tree) #plot(tree, no.margin=TRUE) #tree$Nnode #tree$tip.label -- to check the structure of the tree

#total number of species in dataset
total_species <- length(unique(my_data$Formatted_species))
#total species = 1276

#total number of species with territorial data in dataset
terr_table<- ftable(my_data$Formatted_species, my_data$Territorial) #this has 3 columns, the first column is for NAs
prop_terr<- round(terr_table[,3]/(terr_table[,2]+ terr_table[,3]),2) # so this calculates the percentage "yes" for territorial
sn<- attr(terr_table,"row.vars")[[1]]
terr_data_with_na<- data.frame(sn,sp_terr) #so this is a dataframe with territorial (1/0) for each species
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

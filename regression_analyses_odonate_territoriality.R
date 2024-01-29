#regression analyses for odonate territoriailtiy
wing_data<- read.csv("data/odo_data_wing_pigment.csv") #same dataset as old one (V3) but added wing pigment and climate data


load(file="data/Odo.tree.Waller.Svensson.2017.rda") #odonate tree extracted from Waller and Svensson 2017
#str(tree) #plot(tree, no.margin=TRUE) #tree$Nnode #tree$tip.label -- to check the structure of the tree


#make it into a readable dataframe with only one species per row
#the variables I want are territoriality, wing pigmentation (male and female), wing dichromatism, male/female wing dichromatism, and climate

#territoriality
terr_table<- ftable(wing_data$Formatted_species, wing_data$Territorial) #this has 3 columns, the first column is for NAs
prop_terr<- round(terr_table[,3]/(terr_table[,2]+ terr_table[,3]),2) # so this calculates the percentage "yes" for territorial
#the 2 at the end rounds to 2 decimal 
#set a 3:1 threshhold = 75% threshold
sp_terr <- ifelse(prop_terr >= 0.75, 1, ifelse(prop_terr <= 0.25, 0, NA))
sn<- attr(terr_table,"row.vars")[[1]]
terr_data_with_na<- data.frame(sn,sp_terr) #so this is a dataframe with territorial (1/0) for each species
terr_data<- terr_data_with_na[complete.cases(terr_data_with_na), ] #removed NA values
colnames(terr_data) <- c("sn", "terr")

#wing pigmentation
filtered_male_wing_pig<- subset(wing_data, Male.wing.pigment %in% c("Yes", "No"))
male_wing_pig_var <- ftable(filtered_male_wing_pig$Formatted_species, filtered_male_wing_pig$Male.wing.pigment)
male_wing_pig_var_df <- as.data.frame(male_wing_pig_var)
male_wing_pig_var_df <- male_wing_pig_var_df[, -3]
colnames(male_wing_pig_var_df) <- c("sn", "male_wing_pigment")
#I do not calculate the proportion of my data here because there are no duplicated species rows for this variable
filtered_female_wing_pig<-subset(wing_data, Female.wing.pigment %in% c("Yes", "No"))
female_wing_pig_var<-ftable(filtered_female_wing_pig$Formatted_species, filtered_female_wing_pig$Female.wing.pigment)
female_wing_pig_var_df <- as.data.frame(female_wing_pig_var)
female_wing_pig_var_df <- female_wing_pig_var_df[, -3]
colnames(female_wing_pig_var_df) <- c("sn", "female_wing_pigment")

#wing dichromatism (yes/no)
filtered_wing_dichro<-subset(wing_data, Sexual.wing.dichromatism %in% c("Yes", "No"))
wing_dichro_var<-ftable(filtered_wing_dichro$Formatted_species, filtered_wing_dichro$Sexual.wing.dichromatism)
wing_dichro_var_df <- as.data.frame(wing_dichro_var)
wing_dichro_var_df <- wing_dichro_var_df[, -3]
colnames(wing_dichro_var_df) <- c("sn", "wing_dichromatism")
#I do not calculate the proportion of my data here because there are no duplicated species rows for this variable

#male or female biased sexual wing dichromatism
filtered_biased_wing_dichro<-subset(wing_data, Male.or.female.biased.sexual.wing.dichromatism %in% c("Male", "Female"))
biased_dichro_var<-ftable(filtered_biased_wing_dichro$Formatted_species, filtered_biased_wing_dichro$Male.or.female.biased.sexual.wing.dichromatism)
biased_dichro_var_df <- as.data.frame(biased_dichro_var)
biased_dichro_var_df <- biased_dichro_var_df[, -3]
colnames(biased_dichro_var_df) <- c("sn", "male_or_female_biased_wing_dichromatism")
#I do not calculate the proportion of my data here because there are no duplicated species rows for this variable

#climate
filtered_climate<-subset(wing_data, !is.na(Mean.annual.temperature))
climate_var<-ftable(filtered_climate$Formatted_species, filtered_climate$Mean.annual.temperature)
climate_var_df <- as.data.frame(climate_var)
climate_var_df <- climate_var_df[, -3]
colnames(climate_var_df) <- c("sn", "mean_annual_temperature")

#now stitch all variables together into a single dataframe
reg_data<-merge(terr_data, male_wing_pig_var_df, by = "sn", all=TRUE)
reg_data<-merge(reg_data, female_wing_pig_var_df, by = "sn", all=TRUE)
reg_data<-merge(reg_data, wing_dichro_var_df, by = "sn", all=TRUE)
reg_data<-merge(reg_data, biased_dichro_var_df, by="sn", all=TRUE)
reg_data<-merge(reg_data, climate_var_df, by="sn", all=TRUE)



#rest of example
#I still need to make my data readable
#then I need to prune it to the tree
#then do the regressions. 
colnames(binary_data) <- c("Species", "Prop_Territorial", "Prop_Mate_Guard", "Prop_Flier_vs_Percher", "Prop_Oviposition")
#This worked! I checked manually.

#now remove NAs
mate_guard_terr_data_old<- mate_guard_terr_data_with_na[complete.cases(mate_guard_terr_data_with_na), ] 

#next prune tree to match data
chk<-name.check(tree, terr_data, data.names=as.character(terr_data$sn))
summary(chk)
odonate_tree<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree
# Identify species to drop from terr_data
species_to_drop <- chk$data_not_tree
odonate_terr_data <- terr_data[!(terr_data$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_terr_data, data.names=as.character(odonate_terr_data$sn))

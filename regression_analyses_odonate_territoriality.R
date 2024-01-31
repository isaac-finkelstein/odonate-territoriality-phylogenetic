library("geiger")
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
#male
filtered_male_wing_pig<- subset(wing_data, Male.wing.pigment %in% c("Yes", "No"))
male_wing_pig_var <- ftable(filtered_male_wing_pig$Formatted_species, filtered_male_wing_pig$Male.wing.pigment)
prop_wing_pig_var <- round(male_wing_pig_var[,2]/(male_wing_pig_var[,1]+male_wing_pig_var[,2]),0) #rounds to 0 or 1, with 1 = Yes, 0 = No
male_wing_sn<-attr(male_wing_pig_var, "row.vars")[[1]]
male_wing_pig_df<-data.frame(male_wing_sn, prop_wing_pig_var)
colnames(male_wing_pig_df) <- c("sn", "male_wing_pigment")
#This worked, I checked. 1=pigmented wings, 0 = non-pigmented wings. 
#female
filtered_female_wing_pig<-subset(wing_data, Female.wing.pigment %in% c("Yes", "No"))
female_wing_pig_var<-ftable(filtered_female_wing_pig$Formatted_species, filtered_female_wing_pig$Female.wing.pigment)
prop_female_wing_pig_var<-round(female_wing_pig_var[,2]/(female_wing_pig_var[,1]+female_wing_pig_var[,2]),0)
female_wing_sn<-attr(female_wing_pig_var, "row.vars")[[1]]
female_wing_pig_df <- data.frame(female_wing_sn, prop_female_wing_pig_var)
colnames(female_wing_pig_df) <- c("sn", "female_wing_pigment")
#This worked, I checked. 1=pigmented wings, 0 = non-pigmented wings.

#wing dichromatism (yes/no)
filtered_wing_dichro<-subset(wing_data, Sexual.wing.dichromatism %in% c("Yes", "No"))
wing_dichro_var<-ftable(filtered_wing_dichro$Formatted_species, filtered_wing_dichro$Sexual.wing.dichromatism)
prop_wing_dichro<-round(wing_dichro_var[,2]/(wing_dichro_var[,1]+wing_dichro_var[,2]),0)
wing_dichro_sn<- attr(wing_dichro_var, "row.vars")[[1]]
wing_dichro_var_df <- data.frame(wing_dichro_sn, prop_wing_dichro)
colnames(wing_dichro_var_df) <- c("sn", "wing_dichromatism")
#this worked

#male or female biased sexual wing dichromatism
filtered_biased_wing_dichro<-subset(wing_data, Male.or.female.biased.sexual.wing.dichromatism %in% c("Male", "Female"))
biased_dichro_var<-ftable(filtered_biased_wing_dichro$Formatted_species, filtered_biased_wing_dichro$Male.or.female.biased.sexual.wing.dichromatism)
prop_biased_dichro<-round(biased_dichro_var[,2]/(biased_dichro_var[,1]+biased_dichro_var[,2]),0)
biased_dichro_sn<-attr(biased_dichro_var, "row.vars")[[1]]
biased_dichro_var_df <- data.frame(biased_dichro_sn, prop_biased_dichro)
colnames(biased_dichro_var_df) <- c("sn", "male_or_female_biased_wing_dichromatism")
#1=male biased, 0= female biased. 

#climate
filtered_climate<-subset(wing_data, !is.na(Mean.annual.temperature))
climate_var<-table(filtered_climate$Formatted_species, filtered_climate$Mean.annual.temperature)
climate_summary <- tapply(filtered_climate$Mean.annual.temperature, filtered_climate$Formatted_species, FUN = mean)

# Create a data frame with the summarized information
climate_var_df <- data.frame(sn = names(climate_summary), mean_annual_temperature = climate_summary)

#now stitch all variables together into a single dataframe
reg_data<-merge(terr_data, male_wing_pig_df, by = "sn", all=TRUE)
reg_data<-merge(reg_data, female_wing_pig_df, by = "sn", all=TRUE)
reg_data<-merge(reg_data, wing_dichro_var_df, by = "sn", all=TRUE)
reg_data<-merge(reg_data, biased_dichro_var_df, by="sn", all=TRUE)
reg_data<-merge(reg_data, climate_var_df, by="sn", all=TRUE)
#so there's 921 species total because that's how many rows I have data for at least one variable
#remember, I have territorial data for about 650 species. So About 250 species, I don't have territorial data for -- there may be sample size issues down the road. 

#next prune tree to match data
chk<-name.check(tree, reg_data, data.names=as.character(reg_data$sn))
summary(chk)
odonate_tree<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree
# Identify species to drop from terr_data
species_to_drop <- chk$data_not_tree
odonate_reg_data <- reg_data[!(reg_data$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_reg_data, data.names=as.character(odonate_reg_data$sn))

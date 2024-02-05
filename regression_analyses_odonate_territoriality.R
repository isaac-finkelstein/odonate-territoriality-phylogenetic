library("geiger")
library("phylolm")
#wing pigmentation analyses for odonate territoriailtiy
wing_data<- read.csv("data/data_v4.csv") #same dataset as old one (V3) but added wing pigment and climate data


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
prop_wing_pig_var <- round(male_wing_pig_var[,2]/(male_wing_pig_var[,1]+male_wing_pig_var[,2]),2) #rounds to 0 or 1, with 1 = Yes, 0 = No

#find any conflicts:
con_wing<-ifelse(prop_wing_pig_var <1 & prop_wing_pig_var>0, 2, NA) #see if there is any conflicts
conflicts_sn<-attr(male_wing_pig_var, "row.vars")[[1]]
conflicting_df_male<-data.frame(conflicts_sn, con_wing)
colnames(conflicting_df_male)<-c("sn", "con")
conflicts_male <- conflicting_df_male[conflicting_df_male$con == 2, ]
conflicts_male <- na.omit(conflicts_male)
#so there is a few conflicts -- so I should deal with these. It could be because of the different definitions?

#for now, I am just removing any disagreements in the data
sp_male_wing_pig<-ifelse(prop_wing_pig_var == 1 | prop_wing_pig_var == 0, prop_wing_pig_var, NA)
male_wing_sn<-attr(male_wing_pig_var, "row.vars")[[1]]
male_wing_pig_df<-data.frame(male_wing_sn, sp_male_wing_pig)
male_wing_pig_df<- male_wing_pig_df[complete.cases(male_wing_pig_df), ]#remove NAs (conflicting values)
colnames(male_wing_pig_df) <- c("sn", "male_wing_pigment")
#This worked, I checked. 1=pigmented wings, 0 = non-pigmented wings. 

#female wing pigment
#no conflicts becuase only one entry per species 
filtered_female_wing_pig<-subset(wing_data, Female.wing.pigment %in% c("Yes", "No"))
female_wing_pig_var<-ftable(filtered_female_wing_pig$Formatted_species, filtered_female_wing_pig$Female.wing.pigment)
prop_female_wing_pig_var<-round(female_wing_pig_var[,2]/(female_wing_pig_var[,1]+female_wing_pig_var[,2]),0)
female_wing_sn<-attr(female_wing_pig_var, "row.vars")[[1]]
female_wing_pig_df <- data.frame(female_wing_sn, prop_female_wing_pig_var)
colnames(female_wing_pig_df) <- c("sn", "female_wing_pigment")
#This worked, I checked. 1=pigmented wings, 0 = non-pigmented wings.

#wing dichromatism (yes/no)
#no conflicts because only one entry per species
filtered_wing_dichro<-subset(wing_data, Sexual.wing.dichromatism %in% c("Yes", "No"))
wing_dichro_var<-ftable(filtered_wing_dichro$Formatted_species, filtered_wing_dichro$Sexual.wing.dichromatism)
prop_wing_dichro<-round(wing_dichro_var[,2]/(wing_dichro_var[,1]+wing_dichro_var[,2]),0)
wing_dichro_sn<- attr(wing_dichro_var, "row.vars")[[1]]
wing_dichro_var_df <- data.frame(wing_dichro_sn, prop_wing_dichro)
colnames(wing_dichro_var_df) <- c("sn", "wing_dichromatism")
#this worked

#male or female biased sexual wing dichromatism
#no conflicts because only one entry per species 
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
#so there's 1111 species total because that's how many rows I have data for at least one variable
#remember, I have territorial data for about 450 species. So there's many species I don't have territorial data for -- there may be sample size issues down the road. 

#next prune tree to match data
chk<-name.check(tree, reg_data, data.names=as.character(reg_data$sn))
summary(chk)
odonate_tree<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree
# Identify species to drop from terr_data
species_to_drop <- chk$data_not_tree
odonate_reg_data <- reg_data[!(reg_data$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_reg_data, data.names=as.character(odonate_reg_data$sn))

#what species I would need to collect climate data for
terr_and_male_wing<-merge(terr_data, male_wing_pig_df, by="sn", all=TRUE)
terr_and_male_wing <- na.omit(terr_and_male_wing)
common_species <- intersect(terr_and_male_wing$sn, climate_var_df$sn) #74 species
#so there's about 110 species that I don't have climate data for


#testing a correlation between territoriality and wing pigmentation in males
#description of the method from Waller and Svensson 2013:
#"We used the concentrated-changes test (Pagel 1994) as implemented in Mesquite
#v2.74 to test for a correlation between male wing pigment and latitudinal range
#limits. Pagel’s method controls for the lack of independence be- tween closely 
#related clades and uses a Markov model of evolutionary change, where two models
#are compared. One model assumes independent evolution and the other correlated 
#evolution between the two traits."
#Pagel’s method requires that the traits that are compared must be coded as two
#binary characters"

#first prune the tree again
chk_pagel<-name.check(tree, terr_and_male_wing, data.names=as.character(terr_and_male_wing$sn))
summary(chk_pagel)
pagel_tree<-drop.tip(tree, chk_pagel$tree_not_data) #dropped tree_not_data species
pagel_tree
# Identify species to drop from terr_data
species_to_drop_pagel <- chk_pagel$data_not_tree
pagel_data <- terr_and_male_wing[!(terr_and_male_wing$sn %in% species_to_drop_pagel), ] #dropped data_not_tree species from dataset
name.check(pagel_tree, pagel_data, data.names=as.character(pagel_data$sn))
row_names_pagel <- pagel_data$sn
pagel_data<-data.frame(terr_species = ifelse(pagel_data$terr == 1, "territorial", "non-territorial"),
                                 wing_pig = ifelse(pagel_data$male_wing_pigment == 1, "pigmented", "non-pigmented"))
rownames(pagel_data) <- row_names_pagel
#run pagel 94 model
pagel_terr<-setNames(pagel_data[,1],
                                     rownames(pagel_data))
pagel_wing<-setNames(pagel_data[,2],
                                     rownames(pagel_data))
terr_wing_fit<-fitPagel(pagel_tree, pagel_terr, pagel_wing)
#the independent model has lower AIC! But insignificant p-value. #sample size is 150, so fairly large. 
#plot this
plot(terr_wing_fit, signif=2, cex.main=1, cex.sub=0.8, cex.traits=0.7, cex.rates=0.7, lwd=1)




#sample regression
#note that right now, the climate data is only for the Moore et al., 2021 data
#also Waller et al., 2019 may have some more wing pigment data i could add. 
#prune tree to the data for this regression only
data_for_male_wing_reg<-merge(male_wing_pig_df, climate_var_df, by="sn", all=TRUE)
data_for_male_wing_reg <- na.omit(data_for_male_wing_reg)
chk_first_reg<-name.check(odonate_tree, data_for_male_wing_reg, data.names = as.character(data_for_male_wing_reg$sn))
summary(chk_first_reg)
odonate_tree_first_reg<-drop.tip(odonate_tree, chk_first_reg$tree_not_data)
species_to_drop_first_reg<-chk_first_reg$data_not_tree
odonate_data_for_male_wing_reg<-data_for_male_wing_reg[!(data_for_male_wing_reg$sn %in% species_to_drop_first_reg), ]
name.check(odonate_tree_first_reg, odonate_data_for_male_wing_reg, data.names=as.character(odonate_data_for_male_wing_reg$sn))

mod_male_wing_pig<-phyloglm(male_wing_pigment~mean_annual_temperature, data = odonate_data_for_male_wing_reg, phy=odonate_tree_first_reg, boot=1000)

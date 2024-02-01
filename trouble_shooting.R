#trouble shooting
#the problem: I took the data from Moore et al., 2021.
#there's 319 rows of data
#Moore et al., 2021 apparently have a tree the fully covers the 319 rows of data
#however my tree only covers 134 rows
#We should be using the same tree -- both of us took the tree from from Svennson et al., 2017
#as a rusult of my smaller sample size (I think) when I re-run the regression from Moore et al., 2021
#my result is insiginificant (p=0.4379), but the same regression from Moore et al., 2021 R script is siginificant (p= 0.001008)
#but we are using the same data!

#so I can't figure out how they got 319 rows of coverage.
library("geiger")
library("phylolm")


#the tree from Moore et al., 2021 - I can't find any code for this tree. But in methods they say "taken from Svensson et al., 2017"
drags <- read.newick('drag.phylo.for.sub.tre')
#str(drags$tip.label)

#the data from Moore et al., 2021
dat <- read.csv('moore.inter.data.csv')



#my tree
load(file="data/Odo.tree.Waller.Svensson.2017.rda") #odonate tree extracted from Waller and Svensson 2017

#my data
wing_data<- read.csv("data/odo_data_wing_pigment.csv")
#I need to extract the variables I want from my larger dataset (the data for these variables are extracted from Moore et al., 2021)
#wing pigmentation
#male
filtered_male_wing_pig<- subset(wing_data, Male.wing.pigment %in% c("Yes", "No"))
male_wing_pig_var <- ftable(filtered_male_wing_pig$Formatted_species, filtered_male_wing_pig$Male.wing.pigment)
prop_wing_pig_var <- round(male_wing_pig_var[,2]/(male_wing_pig_var[,1]+male_wing_pig_var[,2]),0) #rounds to 0 or 1, with 1 = Yes, 0 = No
male_wing_sn<-attr(male_wing_pig_var, "row.vars")[[1]]
male_wing_pig_df<-data.frame(male_wing_sn, prop_wing_pig_var)
colnames(male_wing_pig_df) <- c("sn", "male_wing_pigment")

#climate
filtered_climate<-subset(wing_data, !is.na(Mean.annual.temperature))
climate_var<-table(filtered_climate$Formatted_species, filtered_climate$Mean.annual.temperature)
climate_summary <- tapply(filtered_climate$Mean.annual.temperature, filtered_climate$Formatted_species, FUN = mean)
# Create a data frame with the summarized information
climate_var_df <- data.frame(sn = names(climate_summary), mean_annual_temperature = climate_summary)

#combine these into a single dataframe
data_sample_reg<-merge(male_wing_pig_df, climate_var_df, by="sn", all=TRUE)
data_sample_reg <- na.omit(data_sample_reg)
str(data_sample_reg)
#this has 319 observations, same as Moore et al., 2021

#pruning my tree to match my data
#this is where I lose my coverage and end up with a sample size of 134
chk<-name.check(tree, data_sample_reg, data.names = as.character(data_sample_reg$sn))
summary(chk)
tree_first_reg<-drop.tip(tree, chk$tree_not_data)
species_to_drop<-chk$data_not_tree
pruned_data_reg<-data_sample_reg[!(data_sample_reg$sn %in% species_to_drop), ]
rownames(pruned_data_reg)<-pruned_data_reg$sn
name.check(tree_first_reg, pruned_data_reg, data.names=as.character(pruned_data_reg$sn))
str(pruned_data_reg)
#see, 134 observations!


#here is the regression that Moore et al., 2021 perform
#the regression
sample_reg<-phyloglm(male_wing_pigment~mean_annual_temperature, data = pruned_data_reg, phy=tree_first_reg, boot=1000, method = 'logistic_MPLE', btol = 35)
summary(sample_reg)
#not significant (but it should be!)
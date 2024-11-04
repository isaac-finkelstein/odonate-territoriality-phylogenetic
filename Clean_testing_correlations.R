#cleaner version
#testing for associations between variables

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(dplyr)
library(phylolm)
library(ggplot2)

my_data<- read.csv("data/data_vfinal.csv")

load(file="data/Odo.tree.Waller.Svensson.2017.rda") #odonate tree extracated from Waller and Svensson 2017
#str(tree) #plot(tree, no.margin=TRUE) #tree$Nnode #tree$tip.label #check the structure of the tree

#tree from Rocha-Ortega et al., 2020
#same structure as Waller and Svensson 2017 but has different species coverage.
#tree<-read.nexus(file="data/Rocha_ortega_tree.nex")

#my dataset has data that uses the same original source
#remove instances where data is recorded more than once and uses the same original source
rows_to_remove <- c(2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 27, 30, 32, 33, 34, 36, 37, 38, 41, 42, 43, 44, 46, 47, 48, 49, 50, 52, 53, 55, 56, 57, 58, 59, 62, 64, 66, 67, 69, 71, 72, 73, 75, 76, 77, 78, 79, 289, 317, 322, 324, 326, 345, 356, 360, 371, 374, 390, 393, 488, 1068, 1089, 8639, 8658, 8659, 8661, 8662)
my_data <- my_data[-rows_to_remove, ]

#Put data into a readable dataframe - one species per row
terr_table<- ftable(my_data$Formatted_species, my_data$Territorial) #3 columns, the first column is for NAs
prop_terr<- round(terr_table[,3]/(terr_table[,2]+ terr_table[,3]),2) #calculates the percentage "yes" for territorial
#set a 3:1 threshhold = 75% threshold
sp_terr <- ifelse(prop_terr >= 0.75, 1, ifelse(prop_terr <= 0.25, 0, NA))

sn<- attr(terr_table,"row.vars")[[1]]
terr_data_with_na<- data.frame(sn,sp_terr) #dataframe with territorial (1/0) for each species
terr_data<- terr_data_with_na[complete.cases(terr_data_with_na), ] #remove NA values 

#testing correlations between binary traits
#territoriality
binary_terr<-ftable(my_data$Formatted_species, my_data$Territorial)
prop_binary_terr<-round(binary_terr[,3]/(binary_terr[,2]+ binary_terr[,3]),2)
#set a 75% threshold = 3:1 threshold
sp_binary_terr<-ifelse(prop_binary_terr >= 0.75, 1, ifelse(prop_binary_terr <=0.25, 0, NA))
sn<-attr(binary_terr, "row.vars")[[1]]
binary_terr_df<-data.frame(sn,sp_binary_terr, stringsAsFactors = TRUE)


#active behaviour
binary_fly_v_perch<-ftable(my_data$Formatted_species, my_data$Flier.vs.percher)
prop_fly_v_perch<-round(binary_fly_v_perch[,3]/(binary_fly_v_perch[,2]+binary_fly_v_perch[,3]),2) #this is proportion Percher
#1=percher, 0 = flier
#set a 75% threshold = 3:1 threshold
sp_fly_v_perch<-ifelse(prop_fly_v_perch >= 0.75, 1, ifelse(prop_fly_v_perch <=0.25, 0, NA))
sn<-attr(binary_fly_v_perch, "row.vars")[[1]]
binary_fly_v_perch_df<-data.frame(sn, sp_fly_v_perch, stringsAsFactors=TRUE)

#courtship
binary_courtship<-ftable(my_data$Formatted_species, my_data$Courtship)
prop_court<-round(binary_courtship[,3]/(binary_courtship[,2]+binary_courtship[,3]),2) #this is proportion "Yes"
# 1=yes, 0 = no
#set a 75% threshold = 3:1 threshold
sp_courtship<-ifelse(prop_court>=0.75,1, ifelse(prop_court<=0.25, 0, NA))
sn<-attr(binary_courtship, "row.vars")[[1]]
binary_court_df<-data.frame(sn, sp_courtship, stringsAsFactors=TRUE)

#oviposition method (endophytic vs exophytic)
#I need to make all epiphytic = exophytic - to make it binary. 
my_data_mutated <- my_data %>%
  mutate(Oviposition.type..endophytic.vs.exophytic. = ifelse(Oviposition.type..endophytic.vs.exophytic. %in% c("Epiphytic", "Exophytic"), "Exophytic", Oviposition.type..endophytic.vs.exophytic.))
#Convert the variable to a factor with specified levels
my_data_mutated$Oviposition.type..endophytic.vs.exophytic. <- factor(
  my_data_mutated$Oviposition.type..endophytic.vs.exophytic.,
  levels = c("Endophytic", "Exophytic")
)
binary_ovi<-ftable(my_data_mutated$Formatted_species, my_data_mutated$Oviposition.type..endophytic.vs.exophytic.)
prop_ovi<-round(binary_ovi[,2]/(binary_ovi[,1]+binary_ovi[,2]),2) #This is proportion exophytic
#1=exophytic, 0=endophytic
#set a 75% threshold = 3:1 threshold
sp_ovi<-ifelse(prop_ovi >= 0.75, 1, ifelse(prop_ovi <=0.25, 0, NA))
sn<-attr(binary_ovi, "row.vars")[[1]]
binary_ovi_df<-data.frame(sn, sp_ovi, stringsAsFactors=TRUE)

#oviposition habitat (lotic vs lentic)
filtered_lo_len <- subset(my_data, Lotic.vs.lentic..breeding.habitat. %in% c("Lotic", "Lentic"))
lo_len_var<-ftable(filtered_lo_len$Formatted_species, filtered_lo_len$Lotic.vs.lentic..breeding.habitat.)
prop_lo_len<-round(lo_len_var[,2]/(lo_len_var[,1]+lo_len_var[,2]),2) #this is the proportion that is lotic
#set a 75% threshold = 3:1 threshold
sp_lo_len<-ifelse(prop_lo_len >= 0.75, 1, ifelse(prop_lo_len <=0.25, 0, NA))
sn<-attr(lo_len_var, "row.vars")[[1]]
binary_lo_len<-data.frame(sn, sp_lo_len, stringsAsFactors = TRUE)
#1 = lotic, 0 = lentic


#one more step: I need to remove NA values, but if I do that to the entire dataframe, I will have no rows left
#because almost no rows have no missing values

#make dataframes of each pair of traits:
fly_v_perch_terr_data_with_na<-merge(binary_terr_df, binary_fly_v_perch_df, by = "sn", all = TRUE)
ovi_terr_data_with_na<-merge(binary_terr_df, binary_ovi_df, by = "sn", all = TRUE)
court_terr_data_with_na<-merge(binary_terr_df, binary_court_df, by="sn", all=TRUE)
lo_len_terr_data_with_na<-merge(binary_terr_df, binary_lo_len, by = "sn", all = TRUE)
#now remove NAs
fly_v_perch_terr_data_old<-fly_v_perch_terr_data_with_na[complete.cases(fly_v_perch_terr_data_with_na), ]
court_terr_data_old<-court_terr_data_with_na[complete.cases(court_terr_data_with_na), ]
ovi_terr_data_old<-ovi_terr_data_with_na[complete.cases(ovi_terr_data_with_na), ] 
lo_len_terr_data_old<- lo_len_terr_data_with_na[complete.cases(lo_len_terr_data_with_na), ] 
#fitting the pagel (1994) model

#active behaviour + territoriality
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

#plot
plot(fly_v_perch_fit, signif=2, cex.main=1, cex.sub=0.8, cex.traits=0.7, cex.rates=0.7, lwd=1)

max_obs_fly <- nrow(fly_v_perch_terr_data)
breaks_fly <- seq(0, max_obs_fly, by = 20)
ggplot(fly_v_perch_terr_data, aes(x = sp_binary_terr, fill = sp_fly_v_perch)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = NULL, y = "Number of species", fill = "Active behaviour") +
  scale_fill_manual(values = c("Flier" = "lightblue", "Percher" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(breaks = breaks_fly)


#courtship + territoriality
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

#plot
plot(court_fit, signif=2, cex.main=1, cex.sub=0.8, cex.traits=0.7, cex.rates=0.7, lwd=1)

max_obs_court <- nrow(court_terr_data)
breaks_court<- seq(0, max_obs_court, by = 5)
ggplot(court_terr_data, aes(x = sp_binary_terr, fill = sp_courtship)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = NULL, y = "Number of species", fill = "Courtship") +
  scale_fill_manual(values = c("Yes" = "lightblue", "No" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(breaks = breaks_court)


#oviposition method
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

#plot
plot(ovi_fit, signif=2, cex.main=1, cex.sub=0.8, cex.traits=0.7, cex.rates=0.7, lwd=1)

max_obs_ovi <- nrow(ovi_terr_data)
breaks_ovi <- seq(0, max_obs_ovi, by = 10)
ggplot(ovi_terr_data, aes(x = sp_binary_terr, fill = sp_ovi)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = NULL, y = "Number of species", fill = "Oviposition method") +
  scale_fill_manual(values = c("Endophytic" = "lightblue", "Exophytic" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(breaks = breaks_ovi)


#oviposition habitat
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

#plot
plot(lo_len_fit, signif=2, cex.main=1, cex.sub=0.8, cex.traits=0.7, cex.rates=0.7, lwd=1)

#plot this
max_obs_lo_len <- nrow(lo_len_terr_data)
breaks_lo_len <- seq(0, max_obs_lo_len, by = 20)
ggplot(lo_len_terr_data, aes(x = sp_binary_terr, fill = sp_lo_len)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = NULL, y = "Number of species", fill = "Oviposition habitat") +
  scale_fill_manual(values = c("Lotic" = "lightblue", "Lentic" = "darkorange")) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_y_continuous(breaks = breaks_lo_len)

#test if territorial species are more likely to have conservation concerns
#7 ordered categories: Least concern, Near threatened, Vulnerable, Endangered, Critically Endangered, Extinct

#remove "Data Deficient" and "No Data"
redlist_category <- subset(my_data, Redlist.category %in% c("Least Concern", "Near Threatened", "Vulnerable", "Endangered", "Critically Endangered", "Extinct"))
redlist_cat<-ftable(redlist_category$Formatted_species, redlist_category$Redlist.category)

lc<-round(redlist_cat[,4]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)
nt<-round(redlist_cat[,5]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)
vul<-round(redlist_cat[,6]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)
end<-round(redlist_cat[,2]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)
cr_end<-round(redlist_cat[,1]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)
ext<-round(redlist_cat[,3]/(redlist_cat[,1]+redlist_cat[,2]+redlist_cat[,3]+ redlist_cat[,4]+ redlist_cat[,5]+redlist_cat[,6]), 2)

#make c (concern) = all categories except lc, (lc) least concern
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

least_concern_count <- sum(data_redlist_terr$Redlist_category == "least concern")
concern_count <- sum(data_redlist_terr$Redlist_category == "concern")
print(least_concern_count)
print(concern_count)
#I only end up with 5 instances of "concern". 


#test the size of the water body used as breeding habitat
#6 categories:
#Generalist, lotic+small, lotic+medium, lotic+big, lentic+small, lentic+big


#first the lentic categories
new_lentic_lotic_data <- my_data[my_data$Lotic.vs.lentic..breeding.habitat. != "Both", ] #Remove species that have oviposition sites in both lentic and lotic locations
lentic_size_table<- ftable(new_lentic_lotic_data$Formatted_species, new_lentic_lotic_data$Lentic.size)
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
)
sn<- attr(lentic_size_table,"row.vars")[[1]]
size_data_with_na<-data.frame(sn,prop_data)
size_data <- size_data_with_na %>%
  filter(!is.na(prop_large) | !is.na(prop_medium) | !is.na(prop_small))

lentic_size_data <- size_data %>%
  mutate(lentic_size = case_when(
    prop_large == 1 ~ "lentic_large",
    prop_medium == 1 ~ "lentic_medium",
    TRUE ~ "lentic_small"  # Default to "small" if neither prop_large nor prop_medium is 1
  ))
lentic_size <- lentic_size_data[, !(names(lentic_size_data) %in% c("prop_large", "prop_medium", "prop_small"))]


#next lotic categories
new_lentic_lotic_data$Description.of.lotic.oviposition..river..stream. <- gsub("Stream, River", "River, Stream", new_lentic_lotic_data$Description.of.lotic.oviposition..river..stream.)
lotic_size_table<-ftable(new_lentic_lotic_data$Formatted_species, new_lentic_lotic_data$Description.of.lotic.oviposition..river..stream.)

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

sn<- attr(lotic_size_table,"row.vars")[[1]]
lotic_data_with_na<-data.frame(sn,prop_lotic_data)
lotic_size_data <- lotic_data_with_na %>%
  filter(!is.na(prop_stream) | !is.na(prop_both) | !is.na(prop_river))

lotic_size<-lotic_size_data %>%
  mutate(lotic_size = case_when(
    prop_stream == 1 ~"stream",
    prop_both == 1 ~ "both",
    TRUE ~ "river" #default to river if neither prop_stream or prop_both
  ))
lotic_size <- lotic_size[, !(names(lotic_size) %in% c("prop_stream", "prop_both", "prop_river"))]

#add species described as ovipositing in both lotic and lentic habitats (generalists)
paulson_data <- my_data %>%
  filter(grepl("\\(Paulson, 2012\\)|\\(Paulson, 2009\\)", Reference))
generalist_data_table<-ftable(paulson_data$Formatted_species, paulson_data$Lotic.vs.lentic..breeding.habitat.)
prop_generalist<-round(generalist_data_table[,2]/(generalist_data_table[,2]), 2)
sp_generalist<-ifelse(prop_generalist>0.99, 1, ifelse(prop_generalist<0.98, 0, NA))
sn<-attr(generalist_data_table,"row.vars")[[1]]
sp_generalist_with_na<- data.frame(sn,sp_generalist)
sp_generalist<- sp_generalist_with_na[complete.cases(sp_generalist_with_na), ]

colnames(sp_generalist)[colnames(sp_generalist) == "sp_generalist"] <- "generalist"
sp_generalist$generalist[sp_generalist$generalist == 1] <- "generalist"



#combine these

#make dataset
data_lentic_size_territoriality_old <- merge(binary_terr_df, lentic_size, by = "sn", all = TRUE)
colnames(data_lentic_size_territoriality_old) <- c("Species", "Territorial", "Habitat_size")
data_lentic_size_territoriality_old<- data_lentic_size_territoriality_old[complete.cases(data_lentic_size_territoriality_old), ] 

lotic_size_terr_data <- merge(binary_terr_df, lotic_size, by = "sn", all = TRUE)
colnames(lotic_size_terr_data) <- c("Species", "Territorial", "Habitat_size")
lotic_size_terr_data<-lotic_size_terr_data[complete.cases(lotic_size_terr_data), ] 

generalist_terr_data<-merge(binary_terr_df, sp_generalist, by = "sn", all=TRUE)
colnames(generalist_terr_data) <- c("Species", "Territorial", "Habitat_size")
generalist_terr_data<- generalist_terr_data[complete.cases(generalist_terr_data), ]

combined_data <- rbind(
  data_lentic_size_territoriality_old[, c("Species", "Territorial", "Habitat_size")],
  lotic_size_terr_data[, c("Species", "Territorial", "Habitat_size")],
  generalist_terr_data[, c("Species", "Territorial", "Habitat_size")]
)


#find conflicts
common_sn_generalist_lentic <- lentic_size %>%
  inner_join(sp_generalist, by = "sn")
common_sn_generalist_lotic <- lotic_size %>%
  inner_join(sp_generalist, by = "sn")
common_sn_lentic_lotic <- lotic_size %>%
  inner_join(lentic_size, by = "sn")
#remove conflicts:
common_sn <- unique(c(common_sn_generalist_lentic$sn, 
                      common_sn_generalist_lotic$sn, 
                      common_sn_lentic_lotic$sn))
combined_data <- combined_data[!combined_data$Species %in% common_sn, ]

chk_ovi_size<-name.check(tree, combined_data, data.names=as.character(combined_data$Species))
summary(chk_ovi_size)
tree_ovi_size <- drop.tip(tree, chk_ovi_size$tree_not_data)#dropped tree_not_data species
#identify species to drop from data
ovi_size_species_to_drop<-chk_ovi_size$data_not_tree
ovi_size_data<-combined_data[!(combined_data$Species %in% ovi_size_species_to_drop),] #dropped data_not_tree species from dataset
rownames(ovi_size_data)<-ovi_size_data$Species
name.check(tree_ovi_size, ovi_size_data, data.names=as.character(ovi_size_data$Species))


ovi_size_data$Habitat_size <- as.factor(ovi_size_data$Habitat_size)
#set "generalist" as the reference level for Habitat_size
ovi_size_data$Habitat_size <- relevel(ovi_size_data$Habitat_size, ref = "generalist")


#phylogenetic logistic regression
ovi_size_log_reg<-phyloglm(Territorial~Habitat_size, data=ovi_size_data, phy=tree_ovi_size, boot=1000, method='logistic_MPLE', btol= 30)
summary(ovi_size_log_reg)

oddrat_ovi <- exp(coef(ovi_size_log_reg))


#plot
ovi_size_data <- ovi_size_data %>%
  mutate(Territorial = ifelse(Territorial == 1, "Territorial", "Non-Territorial"))

ggplot(ovi_size_data, aes(x = factor(Territorial), fill = factor(Habitat_size, levels = c("lentic_small", "lentic_medium", "lentic_large", "stream", "both", "river", "generalist" )))) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = after_stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = NULL, y = "Number of species", fill = "Oviposition habitat size") +
  scale_fill_manual(values = c("lentic_small" = "darkblue", "lentic_medium" = "darkorange", "lentic_large" = "darkred", "stream" = "lightblue", "both" = "purple", "river" = "turquoise", "generalist" = "darkgrey"),
                    labels = c("lentic_small" = "Small Lentic", "lentic_medium" = "Medium Lentic", "lentic_large" = "Large Lentic", "stream" = "Stream", "both" = "Stream and River", "river" = "River", "generalist" = "Generalist")) +
  scale_x_discrete(labels = c("0" = "Non-Territorial", "1" = "Territorial")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

count_ovi <- ovi_size_data %>%
  group_by(Habitat_size) %>%
  summarise(count = n())




#mate guarding
my_data$Mate.guarding <- trimws(my_data$Mate.guarding)
my_data_filtered <- my_data[!(my_data$Mate.guarding == "Both"), ] #There's only 1 "Both" so I'm removing it.
mate_guarding_table<-ftable(my_data_filtered$Formatted_species, my_data_filtered$Mate.guarding)

prop_no<-round(mate_guarding_table[,3]/(mate_guarding_table[,3]+mate_guarding_table[,4]+mate_guarding_table[,2]), 2)
prop_contact<-round(mate_guarding_table[,2]/(mate_guarding_table[,3]+mate_guarding_table[,4]+mate_guarding_table[,2]), 2)
prop_non_contact<-round(mate_guarding_table[,4]/(mate_guarding_table[,3]+mate_guarding_table[,4]+mate_guarding_table[,2]), 2)
prop_mate_guard_df<-data.frame(prop_no=prop_no,
                                  prop_contact=prop_contact,
                                  prop_non_contact=prop_non_contact)
prop_mate_guard_df<-data.frame(
  prop_no=ifelse(prop_no>0.75,1,NA),
  prop_contact=ifelse(prop_contact>0.75,1,NA),
  prop_non_contact=ifelse(prop_non_contact>0.75, 1, NA))
#3:1 cutoff

sn<-attr(mate_guarding_table, "row.vars")[[1]]
mate_guard_data_na<-data.frame(sn, prop_mate_guard_df)
mate_guard_data<-mate_guard_data_na %>%
  filter(!is.na(prop_no) | !is.na(prop_contact) | !is.na(prop_non_contact))
mate_guarding_three_cat<-mate_guard_data %>%
  mutate(Mate_guarding_cat = case_when(
    prop_no==1 ~"No",
    prop_contact == 1 ~"Contact",
    TRUE ~"Non-contact" #default to non-contact if neither contact or no
  ))
mate_guarding_three_cat<-mate_guarding_three_cat %>%
  select(sn, Mate_guarding_cat)
#make dataset
data_mate_guard_terr_old<-merge(binary_terr_df, mate_guarding_three_cat, by="sn", all=TRUE)
colnames(data_mate_guard_terr_old)<-c("Species", "Prop_territorial", "Mate_guarding_cat")
data_mate_guard_terr_old<-data_mate_guard_terr_old[complete.cases(data_mate_guard_terr_old), ]

#identify species to drop
chk_mate_guard<-name.check(tree, data_mate_guard_terr_old, data.names=as.character(data_mate_guard_terr_old$Species))
summary(chk_mate_guard)
tree_mate_guard<-drop.tip(tree, chk_mate_guard$tree_not_data) #dropped tree_not_data species
#identify species to drop from data
guard_species_to_drop<-chk_mate_guard$data_not_tree
data_mate_guard_terr<-na.omit(data_mate_guard_terr_old[!(data_mate_guard_terr_old$Species %in% guard_species_to_drop),]) #dropped data_not_tree species
rownames(data_mate_guard_terr)<-data_mate_guard_terr$Species
name.check(tree_mate_guard, data_mate_guard_terr, data.names = as.character(data_mate_guard_terr$Species))

data_mate_guard_terr$Mate_guarding_cat <- as.factor(data_mate_guard_terr$Mate_guarding_cat)

#phylogenetic logistic regression
mate_guard_reg<-phyloglm(Prop_territorial~Mate_guarding_cat, data=data_mate_guard_terr, phy=tree_mate_guard, boot=1000, method='logistic_MPLE', btol=10)
summary(mate_guard_reg)

oddrat_mate_guard <- exp(coef(mate_guard_reg))

#plot
custom_levels <- c("No", "Contact", "Non-contact")
ggplot(data_mate_guard_terr, aes(x = factor(Prop_territorial), fill = factor(Mate_guarding_cat, levels = custom_levels))) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", aes(label = stat(count)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = NULL, y = "Number of species", fill = "Mate guarding") +
  scale_x_discrete(labels = c("0" = "Non-Territorial", "1" = "Territorial")) +
  scale_fill_manual(values = c("No" = "darkblue", "Contact" = "darkorange", "Non-contact" = "darkred")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 12),
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
#This only leaves 20 species, and only two of these are non-territorial. 

#Test if abundance predicts territoriality
mod_abund<-phyloglm(Territorial~Abundance, data=data_abundance_terr, phy=tree_abundance, boot=1000, method= 'logistic_MPLE', btol=10)
summary(mod_abund)

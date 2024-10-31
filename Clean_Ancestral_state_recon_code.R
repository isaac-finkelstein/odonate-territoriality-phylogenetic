#cleaner version
#evolutionary models, ancestral state reconstruction, phylogenetic signal

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(dplyr)
library(corHMM)
library(caper)

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

#plot the tree
plot(odonate_tree, type="fan", cex=0.5, fsize=0.6, ftype="i")
str(odonate_tree) #371 species/tips


#Choose a discrete character model (rate-of-transition model)

#four possible models for discreet binary response:
#equal rates (ER)
#all rates different (ARD)
#two irreversible models: going from 0 -> 1,or from 1 -> 0, but not vice versa
terr_mode<-setNames(terr_data_factor$prop_terr, terr_data_factor$sn)

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
#the all-rates-different model has the lowest AIC

#estimate q (the transition rate)
#Q is the expected number of transition given a particular amount of time 
print(fit_er)
print(fit_ard)
print(fit_01)
print(fit_10)

#graph it:
plot(fit_ard, show.zeros=FALSE, mar=rep(0,4), signif=5)
#the ARD (all rates different) model, which had the lowest AIC, shows that it is
#easier/more common to go No -> Yes than to go Yes -> No


#Ancestral state reconstruction
#using marginal ancestral state reconstruction as this is more popular
#and provides uncertainty about specific values

#create a special data frame of species names and trait data
#character data must be numerical interger (0,1)
odonate_data<-data.frame(Genus_sp=names(terr_mode), terr_mode=as.numeric(terr_mode)-1)
head(odonate_data, n=10)

#prune tree to match data
chk<-name.check(tree, odonate_data, data.names=as.character(odonate_data$Genus_sp))
summary(chk)
odonate_tree<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree
# Identify species to drop
species_to_drop <- chk$data_not_tree
odonate_data <- odonate_data[!(odonate_data$Genus_sp%in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, odonate_data, data.names=as.character(odonate_data$Genus_sp))


#fit model
fit_marginal<- corHMM(odonate_tree, odonate_data, node.states = "marginal",
                      rate.cat=1, rate.mat=NULL)
fit_marginal
head(fit_marginal$states) 

#plot
cols<-setNames(c("turquoise", "brown"), levels(terr_mode))
#plot this
plotTree.datamatrix(odonate_tree, as.data.frame(terr_mode),
                    colors=list(cols), header=FALSE, fsize=0.45)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)
nodelabels(pie=fit_marginal$states, piecol=cols, cex=0.3)
tiplabels(pie=to.matrix(terr_mode, sort(unique(terr_mode))), piecol=cols, cex=0.3)

#Ancestral state reconstruction using the ace function from the package "ape"
fit_ARD_again<-ace(odonate_data$terr_mode, odonate_tree, model="ARD", type="discrete", marginal = TRUE)
round(fit_ARD_again$lik.anc, 3)
head(fit_ARD_again$lik.anc)
plotTree(odonate_tree, fsize=0.5, ftype="i")
nodelabels(node=1:odonate_tree$Nnode+Ntip(odonate_tree),
           pie=fit_ARD_again$lik.anc, piecol = cols, cex=0.3)
tiplabels(pie=to.matrix(terr_mode, sort(unique(terr_mode))), piecol=cols, cex=0.3)
legend("topright", legend=levels(terr_mode), pch=22, pt.cex=1.5, pt.bg=cols, bty="n", cex=0.8)



#Calculate phylogenetic signal

#using the D statistic from Fritz & Purvis 2010

#territoriality
signal_terr<-phylo.d(odonate_terr_data, odonate_tree, names.col=sn, binvar=sp_terr)


#for other binary traits I need to create a dataframe for each trait
#active behaviour (flier/percher)
binary_fly_v_perch<-ftable(my_data$Formatted_species, my_data$Flier.vs.percher)
prop_fly_v_perch<-round(binary_fly_v_perch[,3]/(binary_fly_v_perch[,2]+binary_fly_v_perch[,3]),2) #this is proportion Percher
#1=percher, 0 = flier
#set a 75% threshold = 3:1 threshold
sp_fly_v_perch<-ifelse(prop_fly_v_perch >= 0.75, 1, ifelse(prop_fly_v_perch <=0.25, 0, NA))
sn<-attr(binary_fly_v_perch, "row.vars")[[1]]
binary_fly_v_perch_df<-data.frame(sn, sp_fly_v_perch, stringsAsFactors=TRUE) #includes NA

#courtship
binary_courtship<-ftable(my_data$Formatted_species, my_data$Courtship)
prop_court<-round(binary_courtship[,3]/(binary_courtship[,2]+binary_courtship[,3]),2) #this is proportion "Yes"
#1=yes, 0 = no
#set a 75% threshold = 3:1 threshold
sp_courtship<-ifelse(prop_court>=0.75,1, ifelse(prop_court<=0.25, 0, NA))
sn<-attr(binary_courtship, "row.vars")[[1]]
binary_court_df<-data.frame(sn, sp_courtship, stringsAsFactors=TRUE) #includes NA


#oviposition method (endophytic vs exophytic)
#make all epiphytic = exophytic - to make it binary. 
my_data_mutated <- my_data %>%
  mutate(Oviposition.type..endophytic.vs.exophytic. = ifelse(Oviposition.type..endophytic.vs.exophytic. %in% c("Epiphytic", "Exophytic"), "Exophytic", Oviposition.type..endophytic.vs.exophytic.))
# Convert the variable to a factor with specified levels
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
binary_ovi_df<-data.frame(sn, sp_ovi, stringsAsFactors=TRUE) #includes NA

#oviposition habitat (lotic/lentic)
filtered_lo_len <- subset(my_data, Lotic.vs.lentic..breeding.habitat. %in% c("Lotic", "Lentic"))
lo_len_var<-ftable(filtered_lo_len$Formatted_species, filtered_lo_len$Lotic.vs.lentic..breeding.habitat.)
prop_lo_len<-round(lo_len_var[,2]/(lo_len_var[,1]+lo_len_var[,2]),2) #this is the proportion that is lotic
#set a 75% threshold = 3:1 threshold
sp_lo_len<-ifelse(prop_lo_len >= 0.75, 1, ifelse(prop_lo_len <=0.25, 0, NA))
sn<-attr(lo_len_var, "row.vars")[[1]]
binary_lo_len<-data.frame(sn, sp_lo_len, stringsAsFactors = TRUE) #includes NA


#prune binary variables to match the tree
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


#phylogenetic signal of binary traits
signal_terr<-phylo.d(odonate_terr_data, odonate_tree, names.col=sn, binvar=sp_terr)
signal_fly_v_perch<-phylo.d(fly_v_perch_data, tree_fly_v_perch, names.col=sn, binvar=sp_fly_v_perch)
signal_courtship<-phylo.d(court_data, tree_court, names.col=sn, binvar=sp_courtship)
signal_oviposition<-phylo.d(ovi_data, tree_ovi, names.col=sn, binvar = sp_ovi)
signal_lo_len<-phylo.d(lo_len_data, tree_lo_len, names.col = sn, binvar = sp_lo_len)

#phylogenetic signal for non-binary traits
#calculating delta, a measure of phylogenetic signal

#from Borges et al., 2019 github
#https://github.com/mrborges23/delta_statistic
source("code.R") #taken from https://github.com/mrborges23/delta_statistic/tree/master

#all branches must be positive
odonate_tree$edge.length[odonate_tree$edge.length==0] <-quantile(odonate_tree$edge.length, 0.1)*0.1


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

prop_mate_guard_df <- prop_mate_guard_df %>%
  mutate(mate_guard = case_when(
    prop_no == 1 ~ 0,         # 0 = no mate guarding
    prop_contact == 1 ~ 2,     # 2 = contact mate guarding
    prop_non_contact == 1 ~ 1, # 1 = non-contact mate guarding
    TRUE ~ NA_real_            # NA if none of the conditions are met
  ))
prop_mate_guard_df <- prop_mate_guard_df[, !(names(prop_mate_guard_df) %in% c("prop_no", "prop_contact", "prop_non_contact"))]
sn<-attr(mate_guarding_table, "row.vars")[[1]]
mate_guard_df<-data.frame(sn, prop_mate_guard_df, stringsAsFactors=TRUE)
mate_guard_df<- mate_guard_df[complete.cases(mate_guard_df), ] #remove NA values 

#prune tree to match data
chk<-name.check(tree, mate_guard_df, data.names=as.character(mate_guard_df$sn))
summary(chk)
odonate_tree<-drop.tip(tree, chk$tree_not_data) #dropped tree_not_data species
odonate_tree
# Identify species to drop from terr_data
species_to_drop <- chk$data_not_tree
mate_guard_data <- mate_guard_df[!(mate_guard_df$sn %in% species_to_drop), ] #dropped data_not_tree species from dataset
name.check(odonate_tree, mate_guard_data, data.names=as.character(mate_guard_data$sn))

#define the trait vector
mate_guard_data$sn <- droplevels(mate_guard_data$sn) #drop unused levels in the species (sn) column
order_in_mate_guard <- match(odonate_tree$tip.label, mate_guard_data$sn) #must be same order as tree
mate_guard_data_reordered <- mate_guard_data[order_in_mate_guard, ] #tree and data are in the same order
trait_with_na<-mate_guard_data_reordered$prop_mate_guard_df
trait<- trait_with_na[!is.na(trait_with_na)] #removed NA values 

#calculate delta
deltaA<-delta(trait, odonate_tree, 0.1, 0.0589, 10000, 10, 100)

random_delta<-rep(NA,100)
for (i in 1:100){
  rtrait <- sample(trait)
  random_delta[i] <- delta(rtrait, odonate_tree, 0.1, 0.0589, 10000, 10, 100)
}
p_value<-sum(random_delta>deltaA)/length(random_delta)

#finally, let's do this for the size of oviposition habitat
#make 7 categories: lentic small, lentic medium, lentic large, stream, stream and river, river, generalist
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
generalist <- sp_generalist %>%
  rename(habitat_size = generalist)
lotic <- lotic_size %>%
  rename(habitat_size = lotic_size)
lentic <- lentic_size %>%
  rename(habitat_size = lentic_size)

combined_data <- bind_rows(
  generalist[, c("sn", "habitat_size")],  # Select columns by index
  lotic[, c("sn", "habitat_size")],
  lentic[, c("sn", "habitat_size")]
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
combined_data <- combined_data[!combined_data$sn %in% common_sn, ]


#data needs to be numbers
combined_data$habitat_size[combined_data$habitat_size == "generalist"] <- 0
combined_data$habitat_size[combined_data$habitat_size == "lentic_small"] <- 1
combined_data$habitat_size[combined_data$habitat_size == "lentic_medium"] <- 2
combined_data$habitat_size[combined_data$habitat_size == "lentic_large"] <- 3
combined_data$habitat_size[combined_data$habitat_size == "stream"] <- 4
combined_data$habitat_size[combined_data$habitat_size == "both"] <- 5
combined_data$habitat_size[combined_data$habitat_size == "river"] <- 6


#prune tree to match data
chk_hab_size<-name.check(tree, combined_data, data.names=as.character(combined_data$sn))
summary(chk_hab_size)
tree_hab_size <- drop.tip(tree, chk_hab_size$tree_not_data)
#identify species to drop from data
hab_size_species_to_drop<-chk_hab_size$data_not_tree
hab_size_data<-combined_data[!(combined_data$sn %in% hab_size_species_to_drop),] #dropped data_not_tree species from dataset
rownames(hab_size_data)<-hab_size_data$sn
name.check(tree_hab_size, hab_size_data, data.names=as.character(hab_size_data$sn))

hab_size_data$sn <- as.factor(hab_size_data$sn)
hab_size_data$habitat_size <- as.numeric(hab_size_data$habitat_size)
#define the trait vector
order_hab_size <- match(tree_hab_size$tip.label, hab_size_data$sn) #must be same order as tree
hab_size_data_reordered <- hab_size_data[order_hab_size, ] #tree and data are in the same order
hab_trait_with_na<-hab_size_data_reordered$habitat_size
hab_trait<- hab_trait_with_na[!is.na(hab_trait_with_na)] #removed NA values 

#calculate delta
deltaB<-delta(hab_trait, tree_hab_size, 0.1, 0.0589, 10000, 10, 100)

random_delta<-rep(NA,100)
for (i in 1:100){
  rtrait <- sample(trait)
  random_delta[i] <- delta(rtrait, odonate_tree, 0.1, 0.0589, 10000, 10, 100)
}
p_value<-sum(random_delta>deltaB)/length(random_delta)

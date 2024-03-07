#finding conflicts in the data


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
#so 1 = contact, 0 = non-contact. 
#the 2 at the end rounds to 2 decimal 
#set a 75% threshold = 3:1 threshold
sp_binary_mate_guard<-ifelse(prop_mate_guard >= 0.75, 1, ifelse(prop_mate_guard <=0.25, 0, NA))
sn<-attr(mate_guard_var, "row.vars")[[1]]
binary_mate_guard_df<-data.frame(sn, sp_binary_mate_guard, stringsAsFactors = TRUE)

#find conflicts
conflict_mate_guard<-ifelse(prop_mate_guard<1 & prop_mate_guard>0, 2, NA) #conflicts set to 2
conflicts_sn_mate_guard<-attr(mate_guard_var, "row.vars")[[1]]
conflict_mate_guard<-data.frame(conflicts_sn_mate_guard, conflict_mate_guard)
colnames(conflict_mate_guard)<-c("Species", "Conflict")
conflicts_mate_guarding<-conflict_mate_guard[conflict_mate_guard$Conflict ==2,]
conflicts_mate_guarding<-na.omit(conflicts_mate_guarding)
str(conflicts_mate_guarding) #so there are 5 conflicts




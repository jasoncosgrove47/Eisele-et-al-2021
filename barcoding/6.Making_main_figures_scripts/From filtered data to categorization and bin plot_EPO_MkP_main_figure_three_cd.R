library("ggplot2")
library("plyr")
library("reshape2")
library("dplyr")
library("gdata")
library("data.table")


##### 0 Import data
x3<-read.table("/Users/almut/Dropbox/Eisele_et_al_EPO/Rscripts_and_data_by_Almut/Filtered_barcoding_data/EPO_1000ngml_DC_part_I_all_filtered_.txt", sep="\t", header=TRUE)
x2<-read.table("/Users/almut/Dropbox/Eisele_et_al_EPO/Rscripts_and_data_by_Almut/Filtered_barcoding_data/EPO_1000ngml_DC_part_II_all_filtered_.txt", sep="\t", header=TRUE)

x1 <- merge(x3,x2,all=T)

##### 1 Table creation

# make a file per mouse and list
# adapt mouse number according to experiment
H1 <- x1[,grep("_2_", colnames(x1))]
H1<- cbind(x1[,1],H1)
colnames(H1)[1]<-"tag"
H2 <- x1[,grep("_3_", colnames(x1))]
H2<- cbind(x1[,1],H2)
colnames(H2)[1]<-"tag"
H3 <- x1[,grep("_4_", colnames(x1))]
H3<- cbind(x1[,1],H3)
colnames(H3)[1]<-"tag"
H4 <- x1[,grep("_5_", colnames(x1))]
H4<- cbind(x1[,1],H4)
colnames(H4)[1]<-"tag"
Hmouse<-list(H1, H2, H3, H4)


# Take the mean of the replicates and remove rest of values
# Adapt grep names according to experiments
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Bmean"=rowMeans(x[grep("_Bone.B_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Emean"=rowMeans(x[grep("_Bone.E_", names(x))])))
## Merging M samples part
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Monocytesmean"=rowMeans(x[grep("_Bone.Monocytes_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Macrophagesmean"=rowMeans(x[grep("_Bone.Macrophages_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Eosinophilsmean"=rowMeans(x[grep("_Bone.Eosinophils_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Neutrophilsmean"=rowMeans(x[grep("_Bone.Neutrophils_", names(x))])))

# CHosing weights by which is merged
weight = c(0.075, 0.15, 0.075, 0.70)
Hmouse <- lapply(Hmouse, function(x) cbind(x, "Mmean"= apply(x[,c((ncol(x)-3):(ncol(x)))], 1, weighted.mean, weight)))

# Remove the replicate values
Hmouse <- lapply(Hmouse, function(x) x[,-c(2:(ncol(x)-7))])
Hmouse <- lapply(Hmouse, function(x) x[,-c(4:(ncol(x)-1))])

#Split list again
for (i in 1:length(Hmouse)) {
  assign(paste0("HBi", i), as.data.frame(Hmouse[[i]]))
}
HBimouse<-list(HBi1, HBi2, HBi3, HBi4)


# Renormalize per barcode
for (i in 1:length(Hmouse)) {
  HBiX<-as.data.frame(prop.table(as.matrix(HBimouse[[i]][,-c(1)]),margin=1))
  HBimouse[[i]]<-cbind(HBiX, HBimouse[[i]])
  assign(paste0("HBii", i), as.data.frame(HBimouse[[i]]))
}

# make new list and rename columns
HBiimouse<-list(HBii1, HBii2, HBii3, HBii4)
HBiimouse <- lapply(HBiimouse, function(x) setnames(x, c("BmeanR","EmeanR", "MmeanR", "tag", "Bmean", "Emean", "Mmean")))


##### 2 Assign lineage bias

# Function to run on tables to assign classes based on M, E, and B lineage.

allbinarizationwDC <- function( a, y){
  
  y$progenitor<-0
  y$progenitor<-with(y, ifelse((y[,1]>a) & (y[,2]<=a) & (y[,3]<=a), "B1", y$progenitor))
  y$progenitor<-with(y, ifelse((y[,3]>a) & (y[,1]<=a) & (y[,2]<=a), "M1", y$progenitor))
  y$progenitor<-with(y, ifelse((y[,2]>a) & (y[,1]<=a) & (y[,3]<=a), "E1", y$progenitor))
  y$progenitor<-with(y, ifelse((y[,3]>a) & (y[,1]>a) & (y[,2]>a), "all1", y$progenitor))
  y$progenitor<-with(y, ifelse((y[,1]>a) & (y[,3]>a) & (y[,2]<=a), "MB1", y$progenitor))
  y$progenitor<-with(y, ifelse((y[,2]>a) & (y[,3]>a) & (y[,1]<=a), "ME1", y$progenitor))
  y$progenitor<-with(y, ifelse((y[,1]>a) & (y[,2]>a) & (y[,3]<=a), "BE1", y$progenitor))
  y$progenitor<-with(y, ifelse((y[,1]<=a) & (y[,2]<=a) & (y[,3]<=a), "low1", y$progenitor))
  
  y<-y[,-c(1:3)]
  y<-y[,-c(2:4)]
  return(y)
}


# run the function on a list of dataframes
X<- lapply(HBiimouse, function(x) allbinarizationwDC(0.1,x))


# Get dataframes with DC reads to merge binarization with it

Hmouse<-list(H1, H2, H3, H4)


# Take the mean of the replicates and remove rest of values
# Adapt grep names according to experiments
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Bmean"=rowMeans(x[grep("_Bone.B_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"MkPmean"=rowMeans(x[grep("_Bone..MkP_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Emean"=rowMeans(x[grep("_Bone.E_", names(x))])))
## Merging M samples part
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Monocytesmean"=rowMeans(x[grep("_Bone.Monocytes_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Macrophagesmean"=rowMeans(x[grep("_Bone.Macrophages_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Eosinophilsmean"=rowMeans(x[grep("_Bone.Eosinophils_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Neutrophilsmean"=rowMeans(x[grep("_Bone.Neutrophils_", names(x))])))

# CHosing weights by which is merged
weight = c(0.075, 0.15, 0.075, 0.70)
Hmouse <- lapply(Hmouse, function(x) cbind(x, "Mmean"= apply(x[,c((ncol(x)-3):(ncol(x)))], 1, weighted.mean, weight)))

# Remove the replicate values
Hmouse <- lapply(Hmouse, function(x) x[,-c(2:(ncol(x)-8))])
Hmouse <- lapply(Hmouse, function(x) x[,-c(5:(ncol(x)-1))])

#Split list again
for (i in 1:length(Hmouse)) {
  assign(paste0("HBi", i), as.data.frame(Hmouse[[i]]))
}
HBimouse<-list(HBi1, HBi2, HBi3, HBi4)

# Merge result with initial dataframe
for (i in 1:length(HBiimouse)) {
  M<-merge(X[[i]],HBimouse[[i]], all=T)
  assign(paste0("H", i), as.data.frame(M))
}


# Make one big table out of it and assign conditions
# Adapt number of mice according to experiment
# Adapt Condition values according to experiment
NORM<-gdata::combine(H1, H2, H3, H4)

NORM$Condition<- 0
NORM<- within(NORM, Condition[source == "H2"]<- "EPO")
NORM<- within(NORM, Condition[Condition == "0"]<- "No EPO")

####### The resulting table can be used to make plots 
####### The plots in the paper are made in Prism, but here are versions of it

####### contribution of progenitor classes to MkP lineage by Condition

NC<-NORM

NC<- melt(NC[,c("progenitor", "MkPmean", "Condition", "source")], id=c("progenitor","Condition", "source"))
ggplot(NC, aes(y = value, x = variable, colour = as.factor(progenitor))) +
  geom_bar(stat="identity", position="fill")+
  facet_wrap(~Condition)+ 
  scale_colour_manual(values = c("BE1" = "darkred", "ME1" = "forestgreen", "MB1" = "royalblue3", "all1" = "gold", "E1" = "lightcoral", "M1"= "palegreen", "B1"="lightskyblue"))


###### percentage progenitors producing MkP






















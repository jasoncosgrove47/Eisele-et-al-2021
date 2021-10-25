library("ggplot2")
library("plyr")
library("reshape2")
library("dplyr")
library("gdata")
library("data.table")


##### 0 Import data

# x1 will be EPO high
# x2 will be EPO low
x1<-read.table("/Users/almut/Dropbox/Eisele_et_al_EPO/Rscripts_and_data_by_Almut/Filtered_barcoding_data/EPO_1000ngml_HSPC_all_filtered_EPO_9_10_11_Ctrl_13_14
               .txt", sep="\t", header=TRUE)
x2<-read.table("/Users/almut/Dropbox/Eisele_et_al_EPO/Rscripts_and_data_by_Almut/Filtered_barcoding_data/EPO_160ngml_HSPC_all_filtered_EPO_5_6_Ctrl_1_2_3.txt", sep="\t", header=TRUE)


##### 1 Table creation

# make a file per mouse and list
# adapt mouse number according to experiment
H1 <- x1[,grep("_9_", colnames(x1))]
H1<- cbind(x1[,1],H1)
H2 <- x1[,grep("_10_", colnames(x1))]
H2 <- cbind(x1[,1],H2)
H3 <- x1[,grep("_11_", colnames(x1))]
H3<- cbind(x1[,1],H3)
H4 <- x2[,grep("_1_", colnames(x2))]
H4<- cbind(x2[,1],H4)
H5 <- x2[,grep("_2_", colnames(x2))]
H5 <- cbind(x2[,1],H5)
H6 <- x2[,grep("_3_", colnames(x2))]
H6<- cbind(x2[,1],H6)
H7 <- x2[,grep("_6_", colnames(x2))]
H7 <- cbind(x2[,1],H7)
H8 <- x2[,grep("_7_", colnames(x2))]
H8 <- cbind(x2[,1],H8)
H9 <- x1[,grep("_13_", colnames(x1))]
H9<- cbind(x1[,1],H9)
H10 <- x1[,grep("_14_", colnames(x1))]
H10 <- cbind(x1[,1],H10)


Hmouse<-list(H1, H2, H3, H4, H5, H6, H7, H8, H9, H10)


# Take the mean of the replicates and remove rest of values
# Adapt grep names according to experiments
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Bmean"=rowMeans(x[grep("_B_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Emean"=rowMeans(x[grep("_E.P7_|_E_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Mmean"=rowMeans(x[grep("_M_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) x[,-c(2:(ncol(x)-3))])

#Split list again
for (i in 1:length(Hmouse)) {
  assign(paste0("HBi", i), as.data.frame(Hmouse[[i]]))
}
HBimouse<-list(HBi1, HBi2, HBi3, HBi4, HBi5, HBi6, HBi7, HBi8, HBi9, HBi10)


# Renormalize per barcode
for (i in 1:length(Hmouse)) {
  HBiX<-as.data.frame(prop.table(as.matrix(HBimouse[[i]][,-c(1)]),margin=1))
  HBimouse[[i]]<-cbind(HBiX, HBimouse[[i]])
  assign(paste0("HBii", i), as.data.frame(HBimouse[[i]]))
}

# make new list and rename columns
HBiimouse<-list(HBii1, HBii2, HBii3, HBii4, HBii5, HBii6, HBii7, HBii8, HBii9, HBii10)
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

# Import again to get HSC values too

# x4 will be EPO high
# x5 will be EPO low
x4<-read.table("/Users/almut/Dropbox/Eisele_et_al_EPO/Rscripts_and_data_by_Almut/Filtered_barcoding_data/EPO_1000ngml_HSPC_all_filtered_EPO_9_10_11.txt", sep="\t", header=TRUE)
x5<-read.table("/Users/almut/Dropbox/Eisele_et_al_EPO/Rscripts_and_data_by_Almut/Filtered_barcoding_data/EPO_160ngml_HSPC_all_filtered_EPO_5_6_Ctrl_1_2_3.txt", sep="\t", header=TRUE)


H1x2 <- x4[,grep("_9_", colnames(x4))]
H1x2<- cbind(x4[,1],H1x2)
H2x2 <- x4[,grep("_10_", colnames(x4))]
H2x2<- cbind(x4[,1],H2x2)
H3x2 <- x4[,grep("_11_", colnames(x4))]
H3x2<- cbind(x4[,1],H3x2)
H4x2 <- x5[,grep("_1_", colnames(x5))]
H4x2<- cbind(x5[,1],H4x2)
H5x2 <- x5[,grep("_2_", colnames(x5))]
H5x2<- cbind(x5[,1],H5x2)
H6x2 <- x5[,grep("_3_", colnames(x5))]
H6x2<- cbind(x5[,1],H6x2)
H7x2 <- x5[,grep("_6_", colnames(x5))]
H7x2<- cbind(x4[,1],H7x2)
H8x2 <- x5[,grep("_7_", colnames(x5))]
H8x2<- cbind(x5[,1],H8x2)
H9x2 <- x4[,grep("_13_", colnames(x4))]
H9x2<- cbind(x4[,1],H9x2)
H10x2 <- x4[,grep("_14_", colnames(x4))]
H10x2<- cbind(x4[,1],H10x2)

Hmousex2<-list(H1x2,H2x2, H3x2, H4x2, H5x2, H6x2, H7x2, H8x2, H9x2, H10x2)

for (i in 1:length(Hmousex2)) {
  colnames(Hmousex2[[i]])<-sub("var.Totalreads.AE1_HSC_", "", colnames(Hmousex2[[i]]))
  assign(paste0("HB", i), as.data.frame(Hmousex2[[i]]))}

Hmousex2<-list( HB1, HB2, HB3, HB4, HB5, HB6, HB7, HB8, HB9, HB10)

Hmousex2 <- lapply(Hmousex2, function(x) cbind(x,"HSCmean"=rowMeans(x[grep("_HSC_", names(x))])))
Hmousex2 <- lapply(Hmousex2, function(x) cbind(x,"Bmean"=rowMeans(x[grep("_B_", names(x))])))
Hmousex2 <- lapply(Hmousex2, function(x) cbind(x,"Emean"=rowMeans(x[grep("_E.P7_|_E_", names(x))])))
Hmousex2 <- lapply(Hmousex2, function(x) cbind(x,"Mmean"=rowMeans(x[grep("_M_", names(x))])))

Hmouse <- lapply(Hmousex2, function(x) x[,-c(2:(ncol(x)-4))])
Hmouse <- lapply(Hmouse, function(x){colnames(x)[1] <- "tag"; x}) 

for (i in 1:length(Hmouse)) {
  M<-merge(X[[i]],Hmouse[[i]], all=T)
  assign(paste0("H", i), as.data.frame(M))
}


# Make one big table out of it and assign conditions
# Adapt number of mice according to experiment
# Adapt Condition values according to experiment
NORM<-gdata::combine(H1, H2, H3, H4, H5, H6, H7, H8, H9, H10)

NORM$Condition<- 0
NORM<- within(NORM, Condition[source == "H1"]<- "EPO high")
NORM<- within(NORM, Condition[source == "H2"]<- "EPO high")
NORM<- within(NORM, Condition[source == "H3"]<- "EPO high")
NORM<- within(NORM, Condition[source == "H7"]<- "EPO low")
NORM<- within(NORM, Condition[source == "H8"]<- "EPO low")
NORM<- within(NORM, Condition[Condition == "0"]<- "No EPO")

####### The resulting table can be used to make plots 
####### The plots in the paper are made in prism, 
####### but here are versions of it

##### Making bin plot coloured by overlap mature HSPCs

NC<-NORM
NC <- filter(NC, !progenitor %in% c(NA))
NC$mature<-rowSums(NC[,c(4,5,6)]) 
#NC$bonemarrow<-rowSums(NC[,c(7,8,9)]) #P1 taken out or not
NC$bonemarrow<-NC[,c(3)]
NC$overlap<-NC$mature * NC$bonemarrow 
NC$overlap[which(NC$overlap >0)] = "overlap"
NC$overlap[which(NC$overlap ==0)] = "no overlap"


NC<- melt(NC[,c("progenitor","Mmean","Bmean","Emean","Condition", "source", "overlap")], id=c("progenitor","Condition", "source", "overlap"))
ggplot(NC, aes(y = value, x = variable, fill=overlap, group=overlap, alpha=overlap)) +
  geom_bar(stat="identity", position="fill", colour="black")+
  facet_wrap(~Condition)+
  scale_fill_hue(l=40)+
  theme_grey(base_size = 20)+
  scale_alpha_manual(values=c(0, 0.5, 1))

##### contribution of progenitor classes to HSPC reads













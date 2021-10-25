library("ggplot2")
library("plyr")
library("reshape2")
library("dplyr")
library("gdata")
library("data.table")


##### 0 Import data
x1<-read.table("/Users/almut/Dropbox/Eisele_et_al_EPO/Rscripts_and_data_by_Almut/5.Filtered_barcoding_data/EPO_1000ngml_exp1_all_filtered_Ctrl_1_2_3_4_5_EPO_6_7.txt", sep="\t", header=TRUE)


##### 1 Table creation

# make a file per mouse and list
# adapt mouse number according to experiment
H1 <- x1[,grep("_1_", colnames(x1))]
H1<- cbind(x1[,1],H1)
H2 <- x1[,grep("_2_", colnames(x1))]
H2<- cbind(x1[,1],H2)
H3 <- x1[,grep("_3_", colnames(x1))]
H3<- cbind(x1[,1],H3)
H4 <- x1[,grep("_4_", colnames(x1))]
H4<- cbind(x1[,1],H4)
H5 <- x1[,grep("_5_", colnames(x1))]
H5<- cbind(x1[,1],H5)
H6 <- x1[,grep("_6_", colnames(x1))]
H6<- cbind(x1[,1],H6)
H7 <- x1[,grep("_7_", colnames(x1))]
H7<- cbind(x1[,1],H7)
Hmouse<-list(H1, H2, H3, H4, H5, H6, H7)


# Take the mean of the replicates and remove rest of values
# Adapt grep names according to experiments
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Bmean"=rowMeans(x[grep("_B_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Emean"=rowMeans(x[grep("_E.P7_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) cbind(x,"Mmean"=rowMeans(x[grep("_M_", names(x))])))
Hmouse <- lapply(Hmouse, function(x) x[,-c(2:(ncol(x)-3))])

#Split list again
for (i in 1:length(Hmouse)) {
  assign(paste0("HBi", i), as.data.frame(Hmouse[[i]]))
}
HBimouse<-list(HBi1, HBi2, HBi3, HBi4, HBi5, HBi6, HBi7)


# Renormalize per barcode
for (i in 1:length(Hmouse)) {
  HBiX<-as.data.frame(prop.table(as.matrix(HBimouse[[i]][,-c(1)]),margin=1))
  HBimouse[[i]]<-cbind(HBiX, HBimouse[[i]])
  assign(paste0("HBii", i), as.data.frame(HBimouse[[i]]))
}

# make new list and rename columns
HBiimouse<-list(HBii1, HBii2, HBii3, HBii4, HBii5, HBii6, HBii7)
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

# Merge result with initial dataframe
for (i in 1:length(HBiimouse)) {
  M<-merge(X[[i]],HBiimouse[[i]], all=T)
  assign(paste0("H", i), as.data.frame(M))
}

# Make one big table out of it and assign conditions
# Adapt number of mice according to experiment
# Adapt Condition values according to experiment
NORM<-gdata::combine(H1, H2, H3, H4, H5, H6, H7)

NORM$Condition<- 0
NORM<- within(NORM, Condition[source == "H6"]<- "EPO")
NORM<- within(NORM, Condition[source == "H7"]<- "EPO")
NORM<- within(NORM, Condition[Condition == "0"]<- "No EPO")

####### The resulting table can be used to make plots 

##### Making triangle plots

NORM$allreadsBEM<-rowSums(NORM[,c(6,7,8)]) 

ggtern(data=NORM,aes(Bmean,Emean,Mmean, size=allreadsBEM)) + 
  geom_mask() + 
  geom_point(shape=1)+
  theme_bw(base_size = 10)+
  labs(x="B",y="E",z="M",title="100ngml exp1")+
  facet_wrap(~Condition)


####### The plots in the paper are made in Prism, but here are versions of it
####### contribution of progenitor classes to lineages

NC<-NORM

NC<- melt(NC[,c("progenitor","Mmean","Bmean","Emean","Condition", "source")], id=c("progenitor","Condition", "source"))
ggplot(NC, aes(y = value, x = variable, colour = as.factor(progenitor))) +
  geom_bar(stat="identity", position="fill")+
  facet_wrap(~Condition)+ 
  scale_colour_manual(values = c("BE1" = "darkred", "ME1" = "forestgreen", "MB1" = "royalblue3", "all1" = "gold", "E1" = "lightcoral", "M1"= "palegreen", "B1"="lightskyblue"))


#### Make plot of progenitors percentage per mouse 

NC<-NORM
NC[is.na(NC)] <- 0

NC<-NC[which((NC$progenitor!="0")),]

detach("package:plyr", unload=TRUE)
NC<- NC %>% arrange(source, progenitor) %>%
  group_by(source) %>% 
  mutate(progenitornumber = length(source))

NC<- NC %>% arrange(source, progenitor) %>%
  group_by(source,progenitor) %>% 
  mutate(progenitorpercent = (length(progenitor)*100/progenitornumber))

ggplot(NC, aes(y = progenitorpercent, x = progenitor, colour = source)) +
  geom_bar(stat="identity", position="dodge")+
  theme_grey(base_size = 20)


#### Get values for binarization plots for prism

NC<-NORM

NC<- melt(NC[,c("progenitor","Mmean","Bmean","Emean","Condition", "source")], id=c("progenitor","Condition", "source"))

detach("package:ggtern", unload=TRUE)
detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)

NC<- NC %>% arrange(Condition, source, progenitor, variable) %>%
  group_by(Condition, source, progenitor, variable) %>% 
  mutate(bin = sum(value))

NC<-NC[,-5]
NC<-unique(NC)

NC<- NC %>% arrange(Condition, source, variable) %>%
  group_by(Condition, source, variable) %>% 
  mutate(all = sum(bin))

NC<- NC %>% arrange(Condition, source, progenitor, variable) %>%
  group_by(Condition, source, progenitor, variable) %>% 
  mutate(bin =((bin*100)/all))

NC<-NC[,-c(2,6)]

#### Progenitor numbers

NC<-NORM

NC<- melt(NC[,c("progenitor","allreadsBEM","Condition", "source")], id=c("progenitor","Condition", "source"))
NC<-NC[which((NC$value!="0")),]
NC<-NC[,-c(4,5)]

detach("package:ggtern", unload=TRUE)
detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)

NC<- NC %>% arrange(Condition, source, progenitor) %>%
  group_by(Condition, source, progenitor) %>% 
  mutate(prognumber = length(progenitor))

NC<-unique(NC)


#### Barcode numbers

NC<-NORM

NC<- melt(NC[,c("progenitor","allreadsBEM","Condition", "source")], id=c("progenitor","Condition", "source"))
NC<-NC[which((NC$value!="0")),]
NC<-NC[,-c(4,5)]

detach("package:ggtern", unload=TRUE)
detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)

NC<- NC %>% arrange(Condition, source) %>%
  group_by(Condition, source) %>% 
  mutate(prognumber = length(source))

NC<-unique(NC)














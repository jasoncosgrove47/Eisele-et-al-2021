
# Import the data filtered for quantile script
x1<-read.table("/Users/Almut/Dropbox/Institut_Curie/Desktop/Sequencing_data_new_library/Analysis AE3:AE5:AE4 before standardizing protocol for new library/Analysis AE3 High/AE03High filtering/AE03High 16_all filtered.txt", sep="\t", header=TRUE)

# add column with merged reads by mouse
x1$"_1_" <- rowSums(x1[,grep("_1_", names(x1))])
x1$"_2_" <- rowSums(x1[,grep("_2_", names(x1))])
x1$"_3_" <- rowSums(x1[,grep("_3_", names(x1))])
x1$"_4_" <- rowSums(x1[,grep("_4_", names(x1))])
x1$"_5_" <- rowSums(x1[,grep("_5_", names(x1))])
x1$"_6_" <- rowSums(x1[,grep("_6_", names(x1))])
x1$"_7_" <- rowSums(x1[,grep("_7_", names(x1))])
x1$"_8_" <- rowSums(x1[,grep("_8_", names(x1))])
x1$"_9_" <- rowSums(x1[,grep("_9_", names(x1))])
x1$"_10_" <- rowSums(x1[,grep("_10_", names(x1))])

# remove all other columns
rownames(x1)<- x1[,1]
x1<-x1[,-c(1:77)]

# binarize on highest reads for a sepcific barcode in a specific mouse
# values are set to one for a mouse if the barcode has the highest reads in this mice compared to all other mice
# all other values are set to zero
x3<- x1
x3<-t(apply(x3, 1, function(x) replace(x, x== max(x), 1000000)))
x3<-t(apply(x3, 1, function(x) replace(x, x!= 1000000, 0)))
x3<-t(apply(x3, 1, function(x) replace(x, x== 1000000, 1))) 
x3<-as.data.frame(x3)
x3 <- cbind(tag = rownames(x3), x3)
rownames(x3)<-NULL

# Import the file from data filtered for quantile script in which plate indices are combined
# And merge it with the table summarizing highest presence of the barcode over the mice
x1<- AE03High.4_plate.index.combined.txt
x1[,1]<-NULL
x1<- merge(x1,x3, all=TRUE)
x1[is.na(x1)] <- 1 

# Extract all data belonging to one mouse
H1 <- x1[,grep("_1_", colnames(x1))]
H1<- cbind(x1[,1],H1)
dimnames(H1)[[2]][1]<- "tag"
H2 <- x1[,grep("_2_", colnames(x1))]
H2<- cbind(x1[,1],H2)
dimnames(H2)[[2]][1]<- "tag"
H3 <- x1[,grep("_3_", colnames(x1))]
H3<- cbind(x1[,1],H3)
dimnames(H3)[[2]][1]<- "tag"
H4 <- x1[,grep("_4_", colnames(x1))]
H4<- cbind(x1[,1],H4)
dimnames(H4)[[2]][1]<- "tag"
H5 <- x1[,grep("_5_", colnames(x1))]
H5<- cbind(x1[,1],H5)
dimnames(H5)[[2]][1]<- "tag"
H6 <- x1[,grep("_6_", colnames(x1))]
H6<- cbind(x1[,1],H6)
dimnames(H6)[[2]][1]<- "tag"
H7 <- x1[,grep("_7_", colnames(x1))]
H7<- cbind(x1[,1],H7)
dimnames(H7)[[2]][1]<- "tag"
H8 <- x1[,grep("_8_", colnames(x1))]
H8<- cbind(x1[,1],H8)
dimnames(H8)[[2]][1]<- "tag"
H9 <- x1[,grep("_9_", colnames(x1))]
H9<- cbind(x1[,1],H9)
dimnames(H9)[[2]][1]<- "tag"
H10 <- x1[,grep("_10_", colnames(x1))]
H10<- cbind(x1[,1],H10)
dimnames(H10)[[2]][1]<- "tag"

# extract the barcodes with "noise level" (shared) or highest reads (unique) per mouse
H1shared = split(H1,H1$"_1_")[['0']]
H1unique = split(H1,H1$"_1_")[['1']]
H2shared = split(H2,H2$"_2_")[['0']]
H2unique = split(H2,H2$"_2_")[['1']]
H3shared = split(H3,H3$"_3_")[['0']]
H3unique = split(H3,H3$"_3_")[['1']]
H4shared = split(H4,H4$"_4_")[['0']]
H4unique = split(H4,H4$"_4_")[['1']]
H5shared = split(H5,H5$"_5_")[['0']]
H5unique = split(H5,H5$"_5_")[['1']]
H6shared = split(H6,H6$"_6_")[['0']]
H6unique = split(H6,H6$"_6_")[['1']]
H7shared = split(H7,H7$"_7_")[['0']]
H7unique = split(H7,H7$"_7_")[['1']]
H8shared = split(H8,H8$"_8_")[['0']]
H8unique = split(H8,H8$"_8_")[['1']]
H9shared = split(H9,H9$"_9_")[['0']]
H9unique = split(H9,H9$"_9_")[['1']]
H10shared = split(H10,H10$"_10_")[['0']]
H10unique = split(H10,H10$"_10_")[['1']]

# Merge all samples where barcodes have noise level ("shared")
Hshared<-list(H1shared, H2shared, H3shared, H4shared, H5shared, H6shared, H7shared, H8shared, H9shared, H10shared)
Hshared <- Reduce(function(...) merge(..., by="tag", all=TRUE), Hshared)
#Hshared<-Hshared[,-1]

# Removing the binarization identifier columns
Hshared<-Hshared[,!grepl("_1_$", colnames(Hshared))]
Hshared<-Hshared[,!grepl("_2_$", colnames(Hshared))]
Hshared<-Hshared[,!grepl("_3_$", colnames(Hshared))]
Hshared<-Hshared[,!grepl("_4_$", colnames(Hshared))]
Hshared<-Hshared[,!grepl("_5_$", colnames(Hshared))]
Hshared<-Hshared[,!grepl("_6_$", colnames(Hshared))]
Hshared<-Hshared[,!grepl("_7_$", colnames(Hshared))]
Hshared<-Hshared[,!grepl("_8_$", colnames(Hshared))]
Hshared<-Hshared[,!grepl("_9_$", colnames(Hshared))]
Hshared<-Hshared[,!grepl("_10_$", colnames(Hshared))]


# Get different quantiles of the "noise" (shared) barcode reads
Hshared[Hshared == 0] <- NA
quantile(unlist(Hshared),probs=c(0,5,90,95,98,100)/100, na.rm=TRUE)

# The 95th quantile value can be used in the next script to generate the final filtered output.





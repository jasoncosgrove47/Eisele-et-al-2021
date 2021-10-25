
##______________________________________________________##
## Import data after merging per sample (all plate indices)
##______________________________________________________##

m5<-AE03High.4_plate.index.combined.txt
m5 <- m5[,-1]
rownames(m5) <- m5[,1]
m5[,1] <- NULL

##______________________________________________________##
## Filter based on 95th quantile value
##______________________________________________________##

m5 <- apply(m5, 2, function(x) ifelse(x < 186, 0, x))
m5<-as.data.frame(m5)

##______________________________________________________##
## Normalisation reads to 100000 per column (sample)
##______________________________________________________##

norm.data <- apply(m5,2, function(x) (x/sum(x))*100000)
# Export data
write.table(norm.data, "Desktop/AE03High 5_norm 100K-95th.txt",row.names=TRUE,col.names=NA, sep="\t")

##______________________________________________________##
## Get Correlation between replicates of a same sample
##______________________________________________________##

d<- norm.data
#d<- d[,-1]

#transform wild to long
d <- melt(d, id.vars="X") # convert from wide to long format
d2 <- cbind( d, ldply(strsplit(as.character(d[,2]), "_"), identity) ) # split the identifier column
#d2 <- d2[,-8]#delete the empty column due to __
names(d2) <- c("tag","id","var", "exp","prog","mouse","type", "rep")  # rename all columns

#delete the zeros
d3<- d2[which(d2$var>0),]

#make matrix with replicates per samples
#a vs b
da <- d3[which(d3$rep == 'a'),]#take lines of replicate a
da <- da[,!colnames(da)=="rep"]#eliminate the other lines
da <- da[,!colnames(da)=="id"]
row.names(da) <- NULL #élimine la première colonne avec les numeros de ligne Eliminate first column with line numbers
dimnames(da) [[2]][2]<- "vara"#renommer la colonne var en vara

db <- d3[which(d3$rep == 'b'),]#idem avec b
row.names(db) <- NULL
db <- db[,!colnames(db)=="rep"]
db <- db[,!colnames(db)=="id"]
dimnames(db) [[2]][2]<- "varb"

dab <- merge(da,db,all=T)#merge da et db
dab[is.na(dab)]<-0
write.table( dab, "Desktop/AE03High 6_norm 100K a vs b-95th.txt", sep="\t", row.names=F )


#plot self/self, replicate against each other
#dev.copy(png,'selfself exp1 before cor filtering.png') # to save directly as a pdf if too many plots are generated
#x <- dab[which(dab$exp=="MPP"),] # if want to chose a certain subset

theme_set(theme_gray(base_size = 5)) #to change font
qplot(asinh(vara), asinh(varb), data=dab) + facet_wrap(~exp~prog~mouse~type)
#dev.off()
#ggsave(file="self e7 w4 before corr filtering.pdf") # or via export in plot window
#save the figure

##______________________________________________________##
## Filter out samples with low replicate correlation
##______________________________________________________##

x <- ddply(dab, c("exp","prog","mouse","type"), summarize, cor=cor(vara,varb,use="na"))
x[is.na(x)]<-0
hist(x$cor, breaks=50) # plot the distribution of correlation and save figure
write.table( x, "Desktop/AE03High 9_norm 100K a vs b corr linear result-95th.txt", sep="\t", row.names=F ) #save the result of correlation per sample

y <- merge(dab,x, all=TRUE)
z<- y[which(y$cor>0.004),]
Z<- y[which(y$cor<0.004),]
write.table( z, "Desktop/AE03High 10_norm 100K a vs b samples removed by corr filtering-95th.txt", sep="\t", row.names=F ) # save the table with low corraletion samples removed
write.table( Z, "Desktop/AE03High 11_norm 100K a vs b corfiltered-95th.txt", sep="\t", row.names=F ) #save the rows that have been removed 

#plot new self/self after corr filtering
qplot(asinh(vara), asinh(varb), data=z) + facet_wrap(~exp~prog~mouse~type)
#save the figure


##______________________________________________________##
#eliminate bc that are not present in one half of a sample
##______________________________________________________##

z1 <- z[which(z$vara>0 & z$varb>0),]
filt3 <- z[which((z$vara>0 & z$varb==0) | (z$vara==0 & z$varb>0)),]
write.table(filt3, "Desktop/AE03High 13_norm 100K a vs b samples removed by axis filtering-95th.txt",quote=F, sep="\t", row.names=F ) #save the rows that have been removed
write.table(z1, "Desktop/AE03High 14_norm 100K a vs b axis filtered-95th.txt",quote=F, sep="\t", row.names=F )#save the table with barcodes on the axis removed
qplot(asinh(vara), asinh(varb), data=z1) + facet_wrap(~exp~prog~mouse~type)
#save figure

##______________________________________________________##
##transform format back to wide format and save
##______________________________________________________##

z1[,8] <-paste(z1$exp,z1$prog,z1$mouse,z1$type, sep="_")
x1 <- z1[,-c(1,2,3,4)]
xa <- x1[,-3]
names(xa) <- c("tag","var", "id")
xa$id <-paste(xa$id,"a", sep="_")
xb <- x1[,-2]
names(xb) <- c("tag","var", "id")
xb$id <-paste(xb$id,"b", sep="_")
x1<- rbind(xa,xb)

x1 <- reshape(x1,direction="wide", timevar="id", idvar="tag" )
x1[is.na(x1)] <- 0


write.table( x1, "Desktop/AE03High 16_all filtered-95th.txt",quote=F, sep="\t", row.names=F )
















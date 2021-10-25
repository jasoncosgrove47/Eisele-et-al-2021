##______________________________________________________##
## Import data after first filtering steps outside of R
##______________________________________________________##

M<- AE03.High.3_withnames_without_nohits

##______________________________________________________##
## Merge the columns belonging to the same sample (plate indices)
##______________________________________________________##

m <- melt(M, id.vars="tag")
m2 <- cbind( m, ldply(strsplit(as.character(m[,2]), "_"), identity) ) # split the identifier column
m2 <- m2[,-8] # delete the empty column due to __
names(m2) <- c("tag","id","var", "exp","prog","mouse","type", "rep", "plate_index", "sample_index")  # rename all columns
m3<- data.table(m2)
m4<- m3[ , .(Totalreads = sum(var)), by = .(exp, prog, mouse, type, rep, sample_index, tag), drop=FALSE]

# sample_index column can be deleted now
m4[,6]<- NULL
m4$id<-paste(m4$exp,m4$prog,m4$mouse,m4$type, m4$rep, m4$sample_index, sep="_")
m5 <- m4[,-c(1,2,3,4,5)]
m5 <- reshape(m5,direction="wide", timevar="id", idvar="tag" )

# Export data
write.table(m5, "Desktop/AE03High 4_plate index combined.txt.txt",row.names=TRUE,col.names=NA, sep="\t")

##______________________________________________________##
## Make heatmap to check overlap in differents mice
##______________________________________________________##

#heatmap function
mydistM = function(z) dist(z,'manhattan') 
#mydistK = function(z) as.dist(cor(z,method='pearson'))
myhclu = function(d) hclust(d,method='complete')
# make sure that your input is a matrix, else it won't work
library("gplots")

#heatmap with adjustable scale
d <- asinh(m5[,-1]/1000)# the data to be plotted without the tag column transformed as asin values so that zero=zero
d <- d[which(rowSums(d)>0),] # exclude from table the barcodes with no values (should not be any, should there?)
pairs.breaks <- seq(round(min(d)), round(max(d))+1, by=0.5) # usually by 0,5
mycol <- colorpanel(n=length(pairs.breaks)-1,low="black",mid="green",high="red")
mymar <- c(10, 4) + 0.1
myhm = function(y) heatmap.2(y, dendrogram="col", distfun=mydistM, hclustfun=myhclu,scale="none", breaks=pairs.breaks,
                             col=mycol, margins=mymar,
                             labCol = NULL, cexCol=0.9, density.info="none", trace="none", 
                             key = TRUE, keysize = 0.5, key.title=NA, # no title
                             key.xlab=NA)  # no xlabel of color key
m <-as.matrix(d)
x <- myhm(m)

##______________________________________________________##
## Normalisation reads to 100000 per column (sample)
##______________________________________________________##

m5<-as.data.frame(m5)
rownames(m5) <- m5[,1]
data <- m5[,-1]
data <- data[-1295,]
norm.data <- apply(data,2, function(x) (x/sum(x))*100000)

# Export data
write.table(norm.data, "Desktop/AE03High 5_norm 100K.txt.txt",row.names=TRUE,col.names=NA, sep="\t")

##______________________________________________________##
## Get Correlation between replicates of a same sample
##______________________________________________________##

d<- norm.data
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
row.names(da) <- NULL #éEliminate first column with line numbers
dimnames(da) [[2]][2]<- "vara"#rname column var in vara

db <- d3[which(d3$rep == 'b'),]#idem with b
row.names(db) <- NULL
db <- db[,!colnames(db)=="rep"]
db <- db[,!colnames(db)=="id"]
dimnames(db) [[2]][2]<- "varb"

dab <- merge(da,db,all=T)#merge da et db
dab[is.na(dab)]<-0
write.table( dab, "Desktop/AE03High 6_norm 100K a vs b.txt", sep="\t", row.names=F )

theme_set(theme_gray(base_size = 5))
qplot(asinh(vara), asinh(varb), data=dab) + facet_wrap(~exp~prog~mouse~type)

##______________________________________________________##
## Filter out samples with low replicate correlation
##______________________________________________________##

#correlation on self/self
x <- ddply(dab, c("exp","prog","mouse","type"), summarize, cor=cor(vara,varb,use="na"))
x[is.na(x)]<-0
hist(x$cor, breaks=50) # plot the distribution of correlation and save figure
write.table( x, "Desktop/AE03High 9_norm 100K a vs b corr linear result.txt", sep="\t", row.names=F ) #save the result of correlation per sample

y <- merge(dab,x, all=TRUE)
z<- y[which(y$cor>0.004),]
Z<- y[which(y$cor<0.004),]
write.table( z, "Desktop/AE03High 10_norm 100K a vs b samples removed by corr filtering.txt", sep="\t", row.names=F ) # save the table with low corraletion samples removed
write.table( Z, "Desktop/AE03High 11_norm 100K a vs b corfiltered.txt", sep="\t", row.names=F ) #save the rows that have been removed 

#plot new self/self after corr filtering
qplot(asinh(vara), asinh(varb), data=z) + facet_wrap(~exp~prog~mouse~type)

##______________________________________________________##
#eliminate bc that are not present in one half of a sample
##______________________________________________________##

z1 <- z[which(z$vara>0 & z$varb>0),]
filt3 <- z[which((z$vara>0 & z$varb==0) | (z$vara==0 & z$varb>0)),]
write.table(filt3, "Desktop/AE03High 13_norm 100K a vs b samples removed by axis filtering.txt",quote=F, sep="\t", row.names=F ) #save the rows that have been removed
write.table(z1, "Desktop/AE03High 14_norm 100K a vs b axis filtered.txt",quote=F, sep="\t", row.names=F )#save the table with barcodes on the axis removed
qplot(asinh(vara), asinh(varb), data=z1) + facet_wrap(~exp~prog~mouse~type)

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

write.table( x1, "Desktop/AE03High 16_all filtered.txt",quote=F, sep="\t", row.names=F )



#### From here on, move to making your heatmaps or binarization
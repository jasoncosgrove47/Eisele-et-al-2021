###### Process the output of xcalibr for a barcoding experiment  #######
#
# @author: Jason Cosgrove (jason.cosgrove@curie.fr)
# @date: 06/07/2018
#
########################################################################

#input parameters
args <- commandArgs(trailingOnly=TRUE)
path_to_results = as.character("/Users/Almut/Desktop/D")  #root directory for the xcalibr outputs
output_path = as.character("/Users/Almut/Desktop/D") #where you want to send the final matrix to
sample_prefix = as.character("A1166R") #the sample prefix for the xcalibr outputs



# iterate over the samples 
# TODO should add this in as a parameter
for(i in 10:96){
 
  if(i != 11){
  print(i)
  #set the sample code

  sample <-  paste(sample_prefix,i,sep="")

  #read in the processed matrix
  barcode.data <- read.table(paste(path_to_results,"/", sample,"/",sample,"-res-all.txt",sep=""),skip=1)
  
  
  
  #update the rownames so that it takes the plate index and sample index into account
  colnames(barcode.data) <- c("name",paste(sample_prefix,i,"_AACAACCG",sep=""),
                              paste(sample_prefix,i,"_ACGGAATG",sep=""),
                              paste(sample_prefix,i,"_CTAACTCC",sep=""),
                              paste(sample_prefix,i,"_GATGGTCA",sep=""),
                              paste(sample_prefix,i,"_TGGCAGAA",sep=""),
                              paste(sample_prefix,i,"_TGTGACGT",sep=""),
                              paste(sample_prefix,i,"_nohit_cols",sep=""))
  
  

  #combine matrices
  if(i ==10){ processed.data <- barcode.data}                           
  else{ processed.data <- cbind(processed.data,barcode.data[,2:7])}
  }
}

#get rid of the barcodes that dont have any counts
processed.data <- processed.data[rowSums(processed.data[,2:ncol(processed.data)]) > 0,]
#send the processed matrix to file
write.csv(processed.data,paste(output_path,"/barcode_data_processed.csv",sep=""),row.names=FALSE)

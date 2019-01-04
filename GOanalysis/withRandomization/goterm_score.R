

goterm_score <- function(barseq,samples){

library(tidyverse)

######################################################
######################################################

#load GO annotation
GO <- read.delim(file="/Volumes/GoogleDrive/My Drive/genome_data/Yeast/Datasets/go_slim_mapping.tab",sep="\t",header=FALSE)
colnames(GO)[1] = "systematic"
colnames(GO)[2] = "gene"
colnames(GO)[5] = "goterm"

#add systematic gene names to data
barseq2 <- AddSystematic(barseq)

#list of GO terms
GOtib <- as.tibble(GO)
GOlist <- as.vector(GOtib %>% group_by(goterm) %>% summarise())
GOlist <- as.vector(GOlist$goterm)


#dataframe for putting scores
MyMeans <-data.frame(GO = GOlist) #this is where final data will go
uniquesamples <- unique(samples)
for(i in 1:length(uniquesamples)){
  MyMeans$placeholder_name = ""
  names(MyMeans)[names(MyMeans) == "placeholder_name"] <- toString(uniquesamples[i])
}

for(m in 1:length(GOlist)){
######get scoring for each GO term in GOlist
  GOgenes <- as.vector(GO[GO$goterm %in% GOlist[m],]$systematic) #genes in the m'th GOterm
  barseq2_GOgenes <- barseq2[barseq2$systematic %in% GOgenes,] #get data associated with GOterm
  numgenes <- dim(barseq2_GOgenes) #number of genes in dataset from Goterm

  if(numgenes[1] == 0) next

#avg_data: dataframe for putting log2 means of data for each gene
  avg_data <- data.frame(gene = barseq2_GOgenes$gene)
  for(i in 1:length(uniquesamples)){
    avg_data$placeholder_name = ""
    names(avg_data)[names(avg_data) == "placeholder_name"] <- toString(uniquesamples[i])
  }
  
  for(i in 1:dim(avg_data)[1]){
    genedata <- as.vector(as.matrix(barseq2_GOgenes[,1:length(samples)+1][i,])) #data from i'th gene

    for(j in 1:length(uniquesamples)){ #compute mean for each gene from each sample
      grab <- samples == uniquesamples[j]
      avg_data[i,names(avg_data)==uniquesamples[j]] = mean(genedata[grab],na.rm=TRUE)
    }
  }

#avg each sample separately
  for(k in 1:length(uniquesamples)){
    MyMeans[m,toString(uniquesamples[k])] = mean(as.numeric(avg_data[,toString(uniquesamples[k])]),na.rm=TRUE)
  }
}


return(MyMeans)

##get genes from particular GO category
View(barseq2[barseq2$systematic %in% as.vector(GO[GO$goterm %in% "cellular respiration",]$systematic),])


}




#functions
######################################################
######################################################

AddSystematic <- function(data){
library('org.Sc.sgd.db')
xx= as.list(org.Sc.sgdCOMMON2ORF)
yy = sapply(xx, function(x){x[1]})
common=names(xx)
orf = as.vector(yy)
names(common)=orf

#add systematic names to dataset
data$systematic <- as.vector(yy[as.vector(data$gene)])
sysmissing_positions <- is.na(data$systematic)
data$systematic[sysmissing_positions] = as.vector(data$gene[sysmissing_positions])
rmlist <- data$systematic[duplicated(data$systematic)] #get dups
data <-data[!(data$systematic %in% rmlist),] #remove dups

return(data)

}


library(tidyverse)

######################################################
######################################################

#load GO annotation
GO <- read.delim(file="/Volumes/GoogleDrive/My Drive/genome_data/Yeast/Datasets/go_slim_mapping.tab",sep="\t",header=FALSE)
colnames(GO)[1] = "systematic"
colnames(GO)[2] = "gene"
colnames(GO)[5] = "goterm"

#load data
barseq <- read.delim(file="/Volumes/GoogleDrive/My Drive/Projects/BARSeq/GOtermAnalysis/allBarSeq_20181112.txt",sep="\t",header = TRUE)

#annotate samples
samples <- c(1,1,1,2,3,4,5) #label replicates here

#add systematic gene names to data
barseq2 <- AddSystematic(barseq)

#list of GO terms
GOtib <- as.tibble(GO)
GOlist <- as.vector(GOtib %>% group_by(goterm) %>% summarise())
GOlist <- as.vector(GOlist$goterm)


#Atlas:
#MyMeans: For each sample, the average score for each GO term
#avg_data: The sample average of each gene in a specific GO term

#MyMeans_Pairwise: take differences of genes for each GO term and then calculate average coherence of that matrix


#dataframe for putting scores for each sample
MyMeans <-data.frame(GO = GOlist) #this is where final data will go
uniquesamples <- unique(samples)
for(i in 1:length(uniquesamples)){
  MyMeans$placeholder_name = ""
  names(MyMeans)[names(MyMeans) == "placeholder_name"] <- toString(uniquesamples[i])
}

################################################
################################################
#store_pairs: dataframe for putting pairwise scores
#note: if there are N samples, then there are N*(N-1) pairs to go through

store_pairs <- c()
MyMeans_Pairwise <- data.frame(GO = GOlist) 

for(i in 1:(length(uniquesamples)-1)){
  temp <- paste0(uniquesamples[i],uniquesamples[i+1:(length(uniquesamples)-i)])
  store_pairs <- c(store_pairs,temp)
}


for(i in 1:length(store_pairs)){
  MyMeans_Pairwise$placeholder_name = ""
  names(MyMeans_Pairwise)[names(MyMeans_Pairwise) == "placeholder_name"] <- store_pairs[i]
}

TotalPairs = length(store_pairs)

################################################
################################################

for(m in 1:length(GOlist)){
  
######get scoring for each GO term in GOlist
  GOgenes <- as.vector(GO[GO$goterm %in% GOlist[m],]$systematic) #genes in the m'th GOterm
  barseq2_GOgenes <- barseq2[barseq2$systematic %in% GOgenes,] #get data associated with GOterm
  numgenes <- dim(barseq2_GOgenes) #number of genes in dataset and in Goterm

  if(numgenes[1] == 0) next

#avg_data: dataframe for putting log2 means of data for each gene
  

  ################################################################
  ################do each sample 1 at a time
  ################################################################
  
  #create avg_data dataframe
  avg_data <- data.frame(gene = barseq2_GOgenes$gene)
  for(i in 1:length(uniquesamples)){
    avg_data$placeholder_name = ""
    names(avg_data)[names(avg_data) == "placeholder_name"] <- toString(uniquesamples[i])
  }

  #populate avg_data dataframe for this particular GO term
  for(i in 1:dim(avg_data)[1]){
    genedata <- as.vector(as.matrix(barseq2_GOgenes[,1:length(samples)+1][i,])) #data from i'th gene
    for(j in 1:length(uniquesamples)){ #compute mean for each gene from each sample
      grab <- samples == uniquesamples[j]
      avg_data[i,names(avg_data)==uniquesamples[j]] = mean(genedata[grab],na.rm=TRUE)
    }
  }

  #calculate average coherence for each sample --> store in MyMeans
  for(k in 1:length(uniquesamples)){
    MyMeans[m,toString(uniquesamples[k])] = mean(as.numeric(avg_data[,toString(uniquesamples[k])]),na.rm=TRUE)
  } 

  ################################################################
  ################do pairwise comparisons of datasets
  ################################################################
  
#create avg_data_pairwise dataframe
  avg_data_pairwise <- data.frame(gene = barseq2_GOgenes$gene)
  for(i in 1:length(store_pairs)){
    avg_data_pairwise$placeholder_name = ""
    names(avg_data_pairwise)[names(avg_data_pairwise) == "placeholder_name"] <- store_pairs[i]
  }

  counter = 1
#populate avg_data_pairwise dataframe
  for(i in 1:(length(uniquesamples)-1)){k
    for(j in (i+1):(length(uniquesamples))){
      avg_data_pairwise[,paste0(i,j)] = as.numeric(avg_data[,toString(i)]) - as.numeric(avg_data[,toString(j)])
      counter = counter + 1
    }
    #
  }

##calculate average coherence for the DIFFERENCE of each sample --> store in MyMeans_Pairwise
#remember: TotalPairs = number of samples to compare
  for(k in 1:TotalPairs){
    MyMeans_Pairwise[m,colnames(avg_data_pairwise)[k+1]] = mean(as.numeric(avg_data_pairwise[,colnames(avg_data_pairwise)[k+1]]),na.rm=TRUE)
  } 
  


}


write.table(MyMeans,file="/Volumes/GoogleDrive/My Drive/Projects/BARSeq/GOtermAnalysis/GOscores.txt",col.names=NA,sep="\t")
write.table(MyMeans_Pairwise,file="/Volumes/GoogleDrive/My Drive/Projects/BARSeq/GOtermAnalysis/GOscores_Differences.txt",col.names=NA,sep="\t")




#Make histograms
a <- MyMeans_Pairwise[,-1]
a <- gather(a)
a$value <- as.numeric(a$value)
ggplot(gather(a),aes(value)) + 
  geom_histogram(bins=12) +
  facet_wrap(~key,scales="free_x")




##get genes from particular GO category
View(barseq2[barseq2$systematic %in% as.vector(GO[GO$goterm %in% "transposition",]$systematic),])





#get longformat MyMeans_Pairwise and convert to tidy format
pairwise_tidy <- gather(MyMeans_Pairwise,sample,value,-GO)
pairwise_tidy$value <- as.vector(sapply(pairwise_tidy$value, as.numeric))
pairwise_tidy <- as.tibble(pairwise_tidy)

pairwise_tidyfilter <- pairwise_tidy %>% filter(abs(value) > 0.8) 

d <- pairwise_tidyfilter %>% count(GO) %>% 
  arrange(desc(n))
  
p <- ggplot(d,aes(x=reorder(GO,-n),y=n)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p























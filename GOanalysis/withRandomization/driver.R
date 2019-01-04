setwd("/Volumes/GoogleDrive/My Drive/Projects/Calico/BARSeq/GOtermAnalysis/withRandomization")

source("AddSystematic.R")
source("goterm_score.R")

#load data
ourdata <- read.delim(file="/Volumes/GoogleDrive/My Drive/Projects/Calico/BARSeq/GOtermAnalysis/allBarSeq_20181112.txt",sep="\t",header = TRUE)

#sample info
samples <- c(1,1,1,2,3,4,5) #label replicates here

#non-randomized data
scored_data <- goterm_score(ourdata,samples)



####################################
########################################################################
########################################################################
########################################################################
####################################
#strategy: 100 randomized trials
# 1. create empty list
# start for loop: 1:100 [number of randomizations]
# save randomizations to mybiglist
# save goterm resultes in gotermlist


NumTrials = 60
mybiglist <- list()
randscoreslist <- list()
#randomize
for(i in 1:NumTrials){
  ourdata_randomize = ourdata
  for(j in 2:dim(ourdata)[2]){
    ourdata_randomize[,j] <- sample(ourdata_randomize[,j]) #randomize columns
  }
  mybiglist[[i]] <- ourdata_randomize
}


# run goterm_score on randomized data in mybiglist! [THIS IS SLOW]
for(i in 1:NumTrials){
  randscoreslist[[i]] <- goterm_score(mybiglist[[i]],samples)
}

saveRDS(randscoreslist,'gotermscores_ofrandomized_data.rds')
write.table(scored_data,file="/Volumes/GoogleDrive/My Drive/Projects/Calico/BARSeq/GOtermAnalysis/GOscores_20190102.txt",col.names=NA,sep="\t")

####################################
########################################################################
########################################################################
########################################################################
####################################
#########new section
####################

#load randomized_data into df_asnumeric

df <- readRDS('gotermscores_ofrandomized_data.rds')
df_asnumeric <- df

numsamples <- dim(df[[1]])[2]
#convert columns to numeric

for(i in 1:length(df)){
  
  for(j in 2:numsamples){
    
    df_asnumeric[[i]][,j] <- as.numeric(df_asnumeric[[i]][,j])
    
  }
  
}

####################################
########################################################################
########################################################################
########################################################################
####################################


#compare actual results [scored_data] to randomized_data [df_asnumeric]
#scoring > 0.8 means that the value of the actual coherence score is more extreme
#than the value from randomized data

SignificanceMatrix <- scored_data
SignificanceMatrix[,2:dim(SignificanceMatrix)[2]] <- 0

DimX <- dim(scored_data)[1]
DimY <- dim(scored_data)[2]

for(i in 1:DimX){
  
  for(j in 2:DimY){
    
    random_scores <- unlist(lapply(df_asnumeric,function(x) x[i,j]))
    
    actual_score = as.numeric(scored_data[i,j])
    
    print(random_scores)
    print(actual_score)
    
    if(is.na(as.numeric(scored_data[i,j]))){next} #deals with the case that the value is NA
    
    if(actual_score < 0){
      
      SignificanceMatrix[i,j] <- sum(actual_score < random_scores)/length(random_scores)
      
    } else {
      
      SignificanceMatrix[i,j] <- sum(actual_score > random_scores)/length(random_scores)
      
    }
  }
}

####################################
########################################################################
########################################################################
########################################################################
####################################

#once we have the signifiance matrix, we want to only keep data in scored_data if it's significant

Threshold = 0.75
scored_data_thresholded <- scored_data

for(i in 1:DimX){
  for(j in 2:DimY){
    
    if(is.na(SignificanceMatrix[i,j] == TRUE)){
      
      print(c(i,j))
      SignificanceMatrix[i,j] <- 0
    }
    
    if(SignificanceMatrix[i,j] > Threshold){
      
      scored_data_thresholded[i,j] <- scored_data[i,j]
      
    } else {
      
      scored_data_thresholded[i,j] <- 0
      
    }
    
  }
}


write.table(scored_data_thresholded,file="/Volumes/GoogleDrive/My Drive/Projects/Calico/BARSeq/GOtermAnalysis/GOscores_20190103_thresholded.txt",col.names=NA,sep="\t")


####################################
########################################################################
########################################################################
####################################join the thresholed and non-thresholded data

left_join(scored_data,scored_data_thresholded,by="GO") %>% View()


####################################
########################################################################
########################################################################
####################################

D <- ourdata[,2:8]

colMeans(ourdata[,2:8],na.rm=TRUE) #means are all negative

apply(Tgiven, 2, FUN = median)

colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(ourdata[,2],na.rm=TRUE),col=colramp[1],lwd=3,ylim=c(0,1),xlim = c(-2,2))
for(i in 3:8){lines(density(ourdata[,i],na.rm=TRUE),lwd=3,col=colramp[i])}












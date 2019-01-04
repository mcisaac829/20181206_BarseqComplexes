#Hi!
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
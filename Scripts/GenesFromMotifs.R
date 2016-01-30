ID <- c("A","A","A",'B','B','C',"C",'C')
seqnames <- c("a","b",'c','a','c','a','b','c')
dfMotifs <- data.frame(ID, seqnames)
found <- c()
subset(dfMotifs, ID==(dfMotifs$seqnames))
for (i in seq(1, (length(unique(dfMotifs$seqnames))))){
  setMotifs <- subset(dfMotifs, seqnames==(unique(dfMotifs$seqnames)[i]))
  uniek <- unique(setMotifs$ID)
  if (length(uniek) == length(unique(dfMotifs$ID))){
    found[i] <- (unique(dfMotifs$seqnames)[i])
    found[[i]] <- dfMotifs$seqnames[i,]
    found <- found[!is.na(found)]
  }
}  




fullFrame <- unique(fullFrame)
test <- removePathwayGenes(fullFrame)
testvar<-paste(genevector[1], genevector[1], sep=".")
fullFrame[rownames(fullFrame)==testvar,]
write.table(fullFrame, "AllGenes.csv")

removePathwayGenes <- function(dfGenes){
  dfGenes <- unlist(dfGenes)
  indices <- c()
  pathwayGenes <- read.table("PathwayInfo.csv")
  pathwayGenes <- as.vector(t(pathwayGenes))
  pathwayGenes <- pathwayGenes[-1]
  geneVector <- as.vector(fullFrame)
  outFrame <- dfGenes
  for (i in seq(1, (length(pathwayGenes)))){
    
    string <- paste(pathwayGenes[i], pathwayGenes[i], sep=".")
    index<- which(dfGenes == string)
    if(length(index) != 0){
    indices[i] <- index[i]

    dfGenes <- dfGenes[-index+1, drop=FALSE]
    }else{
      print("NF")
    }
    
  }
  
  return(dfGenes)
}

which(fullFrame == "4282.4282")

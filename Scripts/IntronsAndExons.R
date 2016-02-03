#Intron en Exonlengtes
library(Biostrings)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library("hgu95av2.db")
library("BSgenome.Hsapiens.UCSC.hg38")

setwd("C:/Users/jesse/Documents/Bio-informatica/Jaar 3/Periode 2/Bapgc/")
#Gene database en genoom worden gedefinieerd
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genome <- BSgenome.Hsapiens.UCSC.hg38
load("AllCoregulatedGenes.csv")

geneset <- fullFrame[c(1:20),]
geneset <- getEntrez(geneset)
lengths.ma <- matrix(ncol=2)
i<-1

for(i in seq(1, length(geneset))){
  
  #lengths <- tryCatch({
  exon <- exonsBy(txdb, by="gene")[geneset[i]]
  exon.seq <- unlist(getSeq(genome, exon))
  exon.length <- sum(width(exon.seq))
  transcripts <- transcriptsBy(txdb, by="gene")[geneset[i]]
  transcripts.seq <- unlist(getSeq(genome, transcripts))
  intron.length <- sum(width(transcripts.seq))-exon.length
  print(intron.length)
  #print(paste(exon.length, intron.length))
  #print("YOOOO")
  length <- c(exon.length, intron.length)
  print(length)
  lengths.ma <- rbind(lengths.ma, length)
  #return(length)
  #},
 # error = function(err){
  #print("NAh")
 # return(c(NA,NA))
#})

}
dflengths <- as.data.frame(lengths.ma[-1,])
plot(dflengths$V2, dflengths$V1, xlab="Intronlengte", ylab="Exonlengte",
     main="Scatterplot van 20 gecorelugeerde genen")
getEntrez <- function(dfGenes){
  #De functie splitst de ID's op een punt en de unieke lijst
  #Wordt teruggegeven
  entrez <- unlist(strsplit(dfGenes, "[.]"))
  entrez <- unique(entrez)
  return(entrez)
}
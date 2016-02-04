#Intron en Exonlengtes
library(Biostrings)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library("hgu95av2.db")
library("BSgenome.Hsapiens.UCSC.hg38")
library(png)
IntroExons <- function(){
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
    length <- c(exon.length, intron.length)
    lengths.ma <- rbind(lengths.ma, length)
  }
  dflengths <- as.data.frame(lengths.ma[-1,])
  png("IntronsExons.png")
  plot(dflengths$V2, dflengths$V1, xlab="Intronlengte", ylab="Exonlengte",
     main="Scatterplot van 20 gecorelugeerde genen")
  dev.off()
	 
}
getEntrez <- function(dfGenes){
  #De functie splitst de ID's op een punt en de unieke lijst
  #Wordt teruggegeven
  entrez <- unlist(strsplit(dfGenes, "[.]"))
  entrez <- unique(entrez)
  return(entrez)
}

main <- function(args)
{
  if (length(args) > 0)
  {
    if (args[1] == "-h" |  args[1] == "-help" | args[1] == "--h" | args[1] == "--help")
    {
      showUsageInformation()
    }
    else
    {
      setwd(args[1])
      IntroExons()
    }
  }
  
}


main(commandArgs(T))

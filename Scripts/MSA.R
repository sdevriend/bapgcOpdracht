#!/usr/bin/Rscript
library("msa")
library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg38")
library(Biostrings)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(KEGGREST)

MSA <- function(){
  #Gene database en genoom worden gedefinieerd
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genome <- BSgenome.Hsapiens.UCSC.hg38

  #De lijst met gecoreguleerde genen wordt geladen
  load("AllCoregulatedGenes.csv")

  # De bovenste tien hits worden hiervan genomen (moet nog naar 50) en de ID's 
  # Worden omgezet naar Entrez Gene ID's
  msaData <- fullFrame[c(1:50),]
  msaData <- getEntrez(msaData)

  # Alle transcripts van de genen worden geladen,
  # hiervan wordt de sequentie geladen
  x <- 6
  sequences <- DNAStringSetList()
  widths <- c()
  for(x in seq(1,length(msaData))){
    #print(msaData[[x]])
    gene.info <- keggGet(paste("hsa:",msaData[[x]], sep=""))
    sequences[[x]] <- gene.info[[1]]$NTSEQ
    widths[[x]] <- width(gene.info[[1]]$NTSEQ)
  }
  #alignment<- msa(unlist(sequences))
# Een DNASeqSet wordt gemaak en uniek gemaakt.

  alignment <- msa(unlist(sequences))
# Het langste gen wordt opgeslagen voor later gebruik
  longestGene <- which(widths == max(widths))
  longestGene.name <- msaData[longestGene]
  longest.seq <- sequences[longestGene]



  

  #longest.gene <- longestGene # TESTVAR
  
  perc <- getConservedPercentage(longest.gene, alignment)
  print(perc)
}
getEntrez <- function(dfGenes){
  entrez <- unlist(strsplit(dfGenes, "[.]"))
  entrez <- unique(entrez)
  return(entrez)
}
getConservedPercentage <- function(longest.gene, alignment){
  
  alignment.matrix <- consensusMatrix(unmasked(alignment))
  alignment.unmasked <- unmasked(alignment)
  consensus <- gsub("[?]","", consensusString(alignment.matrix))
  consensus.vector <- unlist(strsplit(consensus, split=""))
  
  NT.count<-length(consensus.vector)
  longest.seq <- alignment.unmasked[longest.gene]
  longest <- gsub("[?]","", longest.seq)
  longest.vector <- unlist(strsplit(longest, split=""))
  longest.count <- length(longest.vector)
  count <- 0
  for(y in seq(1, length(consensus.vector))){
    if(consensus.vector[[y]] == longest.vector[[y]]){
      count <- count + 1
    }
  }
  percentage <- count/NT.count*100
  return(percentage)
}

printSplitString <- function(x, width=getOption("width") - 1)
{
  starts <- seq(from=1, to=nchar(x), by=width)
  for (i in 1:length(starts))
    cat(substr(x, starts[i], starts[i] + width - 1), "\n")
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
      MSA()
    }
  }
  
}


main(commandArgs(T))
#!/usr/bin/Rscript
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(JASPAR2014)
library(TFBSTools)
library(Biostrings)

showUsageInformation <- function()
{
  print("") 
  print("Dit script geeft als output een bestand met alle 
        promotorsequenties voor een lijst met alle promotorsequenties 
        van een lijst met genen als input.")
  quit()
}


setMotifOnGenome <- function(){
  
  
  lstMatchPos <- list(c(1, 2500), c(2501,5000), c(5001, 7500), c(7501, 10000),
				      c(10001, 12500), c(12501, 15000), c(15001, 17500), 
				      c(17501, 20000), c(20001, 22500), c(22501, 23779))
  foundMotifsTable <- read.csv("FoundMotifs.csv")
  foundMotifs <- foundMotifsTable$x
  amount <- length(foundMotifs)
  dfTenth1 <- getMatches(1, 2500)
  save(dfTenth1, file='motifspart1.RDATA')
  load('motifspart1.RDATA')
  geneset11 <- getGenes(dfTenth1, amount)

  dfTenth2 <- getMatches(2501,5000)
  save(dfTenth2, file='motifspart2.RDATA')
  load('motifspart2.RDATA')
  geneset12 <- getGenes(dfTenth2, amount)

  dfTenth3 <- getMatches(5001, 7500)
  save(dfTenth3, file='motifspart3.RDATA')
  geneset3 <- getGenes(dfTenth3, amount)

  dfTenth4 <- getMatches(7501, 10000)
  save(dfTenth4, file='motifspart4.RDATA')
  geneset4 <- getGenes(dfTenth4, amount)

  dfTenth5 <- getMatches(10001, 12500)
  save(dfTenth5, file='motifspart5.RDATA')
  geneset5 <- getGenes(dfTenth5, amount)

  dfTenth6 <- getMatches(12501, 15000)
  save(dfTenth6, file='motifspart6.RDATA')
  geneset6 <- getGenes(dfTenth6, amount)

  dfTenth7 <- getMatches(15001, 17500)
  save(dfTenth7, file='motifspart7.RDATA')
  geneset7 <- getGenes(dfTenth7, amount)

  dfTenth8 <- getMatches(17501, 20000)
  save(dfTenth8, file='motifspart8.RDATA')
  geneset8 <- getGenes(dfTenth8, amount)

  dfTenth9 <- getMatches(20001, 22500)
  save(dfTenth9, file='motifspart9.RDATA')
  geneset9 <- getGenes(dfTenth9, amount)

  dfTenth10 <- getMatches(22501, 23779)
  save(dfTenth10, file='motifspart10.RDATA')
  geneset10 <- getGenes(dfTenth10, amount)

  fullFrame <- data.frame(c(geneset1,geneset10, geneset2, geneset3, 
                          geneset4, geneset5, geneset6, geneset7, geneset8, geneset9))
  save(fullFrame, file="AllCoregulatedGenes.csv")

}

getMatches <- function(left, right){
  foundMotifsTable <- read.csv("FoundMotifs.csv")
  foundMotifs <- foundMotifsTable$x
  genome <- BSgenome.Hsapiens.UCSC.hg38
  Txdb<- TxDb.Hsapiens.UCSC.hg38.knownGene
  allGenes <- genes(Txdb)
  allTranscripts <- transcriptsBy(Txdb, by="gene")[allGenes$gene_id]
  Transcripts <- allTranscripts[IRanges(start=left, end=right)]
  PromoterSeqs <- getPromoterSeq(Transcripts, genome, 3000, 100)
  opts = list()
  opts[["ID"]] <- foundMotifs
  opts[["matrixtype"]] <- "PWM"
  
  # Maakt motifList m.b.v JASPAR2014. De opties voor mens en PWM worden
  # Meegegeven
  motifListSubSet <- getMatrixSet(JASPAR2014, opts)
  
  AllSSeqs <- searchSeq(motifListSubSet, unlist(PromoterSeqs), min.score="80%")
  dfSseqs<- as(AllSSeqs, "DataFrame")
  return(dfSseqs)
  }

getGenes <- function(dfMotifs, amount){
  found <- c()
  for (i in seq(1, (length(unique(dfMotifs$seqnames))))){
    setMotifs <- subset(dfMotifs, seqnames==(unique(dfMotifs$seqnames)[i]) & relScore > 0.97)
    uniek <- unique(setMotifs$ID)
    if (length(uniek) == amount){
      found[[i]] <- (unique(dfMotifs$seqnames)[i])
      found <- found[!is.na(found)]
    }
  }  
  return(found)
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
      setMotifOnGenome()
    }
  }
  
}


main(commandArgs(T))

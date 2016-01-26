#!/usr/bin/Rscript
setwd("C:/Users/jesse/Documents/Bio-informatica/Jaar 3/Periode 2/Bapgc/")
# Global variables and libraries
biocLite("GenomicRanges")
biocLite("BSgenome.Hsapiens.UCSC.hg38")

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

showUsageInformation <- function()
{
  print("") 
  print("Dit script geeft als output een bestand met alle 
        promotorsequenties voor een lijst met alle promotorsequenties 
        van een lijst met genen als input.")
  quit()
}

getPromotorSequences <- function()
{
  # Gebruikte libraries:
  library(Biostrings)
  library(org.Hs.eg.db)
  library (BSgenome.Hsapiens.UCSC.hg38)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  # Gene Database (HG38)
  txdb<- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  # Het Hsa genoom
  genome <- BSgenome.Hsapiens.UCSC.hg38
  
  # De Entrezgene codes van de pathway inladen
  genes <- read.table("PathwayInfo.csv")
  
  # Het dataframe omzetten in een vector
  genevector <- as.vector(t(genes))
  genevector <- genevector[-1]
  
  # De transcripts van de lijst met genen zoeken
  transcripts.grL <- transcriptsBy(txdb, by="gene")[genevector]
  
  # De namen van de transcripts opslaan (misschien later nodig)
  transcript.names <- mcols(unlist(transcripts.grL))$tx_name  
  
  # De promotersequenties worden opgehaald met 3000 upstream 300 
  promoseqs  = getPromoterSeq(transcripts.grL, genome, upstream=3000, downstream = 300)
  # Schrijf naar rdata file.
  save(promoseqs, "promoseqs.RData")
  # Schrijf de sequenties naar een csv bestand (als tussenproduct)
  #write.table(promoseqs, "PromoterSequences.csv", row.names=FALSE)
  
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
      getPromotorSequences()

    }
  }
  else
  {
    getPromotorSequences()
  }
  
}

main(commandArgs(T))

# Additional information:
# =======================
#
# Remarks about the Skeleton Script R itself.
# Description how it works.
# Description which improvements can be done to improve the Skeleton Script R itself.
#


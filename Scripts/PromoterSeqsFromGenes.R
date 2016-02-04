#!/usr/bin/Rscript


library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


showUsageInformation <- function()
{
  print("") 
  print("Dit script geeft als output een bestand met alle 
        promotorsequenties voor een lijst met alle promotorsequenties 
        van een lijst met genen als input.")
  print("Het script kan aangeroepen worden met:")
  print("./PromoterSeqsFromGenes Werkmap")
  print("Let er wel op dat bestanden van eerdere stappen")
  print(" nodig zijn voor een succesvolle run voor dit script.")
  quit()
}

getPromotorSequences <- function()
{
  # De functie maakt gebruik van packages om promotor sequenties
  # op te halen van de genen die bij een pathway horen.
  # De promotor sequenties worden gezocht in de regio's 3000 upstream
  # en 300 downstream voor optimale resultaten. De promotor
  # sequenties die gevonden worden, worden opgeslagen in de variable
  # promoseqs. Deze wordt opgeslagen als R.Datafile om in latere
  # stappen gebruikt te worden.
  # Gene Database (HG38)
  txdb<- TxDb.Hsapiens.UCSC.hg38.knownGene
  # Het Hsa genoom
  genome <- BSgenome.Hsapiens.UCSC.hg38
  # De Entrezgene codes van de pathway inladen
  dfGenes <- read.table("PathwayInfo.csv")
  # Het dataframe omzetten in een vector
  vctGenevector <- as.vector(t(dfGenes))
  vctGenevector <- vctGenevector[-1]
  # De transcripts van de lijst met genen zoeken
  transcripts.grL <- transcriptsBy(txdb, by="gene")[vctGenevector]
  # De promotersequenties worden opgehaald met 3000 upstream 300 
  promoseqs  = getPromoterSeq(transcripts.grL, genome, upstream=3000, downstream = 300)
  # Schrijf naar rdata file.
  save(promoseqs, file="promoseqs.RData")
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
}

main(commandArgs(T))

# Additional information:
# =======================
#
# Dit script is geschreven door Jesse Kerkvliet en
# kan gebruikt worden bij de opdracht van bapgc door
# Sebastiaan de Vriend.
# Het script zoekt promotorsequenties van genen van een ingeladen
# pathway. Het script kan als volgt aangeroepen worden:
# ./PromoterSeqsFromGenes Werkmap
# Als output wordt een Rdata file opgeslagen in de werkmap.
# Deze kan gebruikt worden voor vervolgstappen.



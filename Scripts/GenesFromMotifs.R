#!/usr/bin/Rscript


library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(JASPAR2014)
library(TFBSTools)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


showUsageInformation <- function()
{
  print("") 
  print("Dit script laad gevonden motifs uit FoundMotifs.csv en ")
  print("zoekt naar genen die bij de motifs horen. Dit gebeurt door de")
  print(" genen op te verdelen in tien stuks. De gevonden genen worden")
  print(" later samengevoegd en opgeslagen in AllCoregulatedGenes.csv.")
  print(" Het script kan als volgt worden aangeroepen:")
  print(" ./GenesFromMotifs.R Werkmap")
  quit()
}

setMotifOnGenome <- function(){
  # De functie is een wat langere functie. Dit is omdat de input data
  # in tien stukken wordt gehakt om het zo mogelijk te maken om
  # het script op een computer te draaien zonder dat er memmory
  # errors opkomen.
  # De functie leest als eerste FoundMotifs.csv in en zet deze in
  # de variable foundMotifs. Daarna worden 2500 genen geladen met
  # de functie getMatches en worden deze opgeslagen. Daarna worden de
  # genen opgehaald met de functie getGenes. Als laaste wordt het
  # totale resultaat samengevoegd en weggeschreven naar
  # AllCoregulatedGenes.csv
  foundMotifsTable <- read.csv("FoundMotifs.csv")
  foundMotifs <- foundMotifsTable$x
  amount <- length(foundMotifs)
  dfTenth1 <- getMatches(1, 2500)
  geneset1 <- getGenes(dfTenth1, amount)

  dfTenth2 <- getMatches(2501,5000)
  geneset2 <- getGenes(dfTenth2, amount)

  dfTenth3 <- getMatches(5001, 7500)
  geneset3 <- getGenes(dfTenth3, amount)

  dfTenth4 <- getMatches(7501, 10000)
  geneset4 <- getGenes(dfTenth4, amount)

  dfTenth5 <- getMatches(10001, 12500)
  geneset5 <- getGenes(dfTenth5, amount)

  dfTenth6 <- getMatches(12501, 15000)
  geneset6 <- getGenes(dfTenth6, amount)

  dfTenth7 <- getMatches(15001, 17500)
  geneset7 <- getGenes(dfTenth7, amount)

  dfTenth8 <- getMatches(17501, 20000)
  geneset8 <- getGenes(dfTenth8, amount)

  dfTenth9 <- getMatches(20001, 22500)
  geneset9 <- getGenes(dfTenth9, amount)

  dfTenth10 <- getMatches(22501, 23779)
  geneset10 <- getGenes(dfTenth10, amount)

  fullFrame <- data.frame(c(geneset1,geneset10, geneset2, geneset3, 
                          geneset4, geneset5, geneset6, geneset7, 
						  geneset8, geneset9))
  save(fullFrame, file="AllCoregulatedGenes.csv")
}

getMatches <- function(left, right){
  # input: 2
  # left: nummer van gen om links mee te beginnen
  # right: nummer van gen om mee te eindigen.
  #
  # De functie haalt de data van FoundMotifs op en gebruikt
  # deze samen met genome om alle transcripten op te zoeken
  # die voor de range van left en right gelden. Daarna worden motifs
  # opgezocht voor de mens die uit FoundMotifs komen en daar worden
  # matchende genen voor gezocht door searchSeq. Als laaste
  # wordt er een dataframe terug gegeven met de gevonden genen en
  # daarbij de scores.
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
  # Input: 2
  # dfMotifs: Dataframe met Motifs en de daarbij behorende score en genen.
  # amount: : integer met lengte van aantal motifs van pathway.
  #
  # De functie kijkt per subset motifs of de score hoger is dan 97%
  # daarna worden alle unieke motifs toegevoegd aan de vector found.
  # Daarna worden alle unieke motifs terug gegeven.
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
# Additional information:
# =======================
# Dit script is geschreven door Jesse Kerkvliet en
# kan gebruikt worden bij de opdracht van bapgc door
# Sebastiaan de Vriend.
# Het script zoekt naar genen die matchen voor motifs uit het bestand
# FoundMotifs.csv en slaat deze op. 
# Het script kan als volgt worden aangeroepen:
# ./GenesFromMotifs.R Werkmap
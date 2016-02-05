#!/usr/bin/Rscript


library(Biostrings)
library("hgu95av2.db")
library("BSgenome.Hsapiens.UCSC.hg38")
library(org.Hs.eg.db)
library(png)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


showUsageInformation <- function()
{
  print("") 
  print("Genereert een scatterplot van de lengten van intronen en exonen")
  print("aan de hand van de twintig beste medegereguleerde genen.")
  print("Het script kan aangeroepen worden door IntronsAndExons.R")
  print("")
  quit()
}

IntroExons <- function(){
  #Gene database en genoom worden gedefinieerd
  #De medegereguleerde genen worden ingeladen 
  #en de 20 beste hits worden meegenomen en omgezet naar Entrez ID's
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genome <- BSgenome.Hsapiens.UCSC.hg38
  load("AllCoregulatedGenes.csv")
  geneset <- fullFrame[c(1:20),]
  geneset <- getEntrez(geneset)
  
  #Voor ieder gen worden de exonen en de transcripten opgehaald.
  #De lengte hiervan wordt berekend en aan een matrix toegevoegd.
  mt.lengths <- matrix(ncol=2)
  for(i in seq(1, length(geneset))){
    exon <- exonsBy(txdb, by="gene")[geneset[i]]
    exon.seq <- unlist(getSeq(genome, exon))
    exon.length <- sum(width(exon.seq))
    transcripts <- transcriptsBy(txdb, by="gene")[geneset[i]]
    transcripts.seq <- unlist(getSeq(genome, transcripts))
    intron.length <- sum(width(transcripts.seq))-exon.length
    length <- c(exon.length, intron.length)
    mt.lengths <- rbind(mt.lengths, length)
  }
  #De matrix met lengten wordt omgezet naar een dataframe en een plot
  #Wordt weggeschreven naar een PNG bestand.
  dflengths <- as.data.frame(mt.lengths[-1,])
  png("IntronsExons.png")
  plot(dflengths$V2, dflengths$V1, xlab="Intronlengte", ylab="Exonlengte",
     main="Scatterplot van 20 gecorelugeerde genen")
  dev.off() 
}

getEntrez <- function(vctGenes){
  #De functie splitst de genen op een punt en geeft de 
  #Entrez ID's terug
  vctEntrez <- unlist(strsplit(as.character(vctGenes), "[.]"))
  vctEntrez <- unique(vctEntrez)
  return(vctEntrez)
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

# Additional information:
# =======================
# Dit script is geschreven door Jesse Kerkvliet en
# kan gebruikt worden bij de opdracht van bapgc door
# Sebastiaan de Vriend.
#
# Een scatterplot van de lengten van de intronen en exonen
# wordt gegenereerd voor de beste 20 medegereguleerde genen.

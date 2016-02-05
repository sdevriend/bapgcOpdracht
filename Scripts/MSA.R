#!/usr/bin/Rscript
library("msa")
library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg38")
library(Biostrings)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(KEGGREST)

showUsageInformation <- function()
{
  print("") 
  print("Berkent het percentage conservatie van het langste gen")
  print("aan de hand van de vijftig beste medegereguleerde genen.")
  print("Dit wordt gedaan met behulp van eem Multiple Sequence Alignment (MSA)")
  print("Het script kan aangeroepen worden door MSA.R")
  print("")
  quit()
}

MSA <- function(){
  #Gene database en genoom worden gedefinieerd
  #De lijst met gecoreguleerde genen wordt geladen
  # De bovenste vijftig hits worden hiervan genomen en de ID's 
  # Worden omgezet naar Entrez Gene ID's
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genome <- BSgenome.Hsapiens.UCSC.hg38
  load("AllCoregulatedGenes.csv")
  msaData <- fullFrame[c(1:50),]
  msaData <- getEntrez(msaData)

  # De sequenties en de lengten van de genen worden in een DNAStringSetList
  # en een vector opgeslagen
  sequences <- DNAStringSetList()
  widths <- c()
  for(x in seq(1,length(msaData))){
    gene.info <- keggGet(paste("hsa:",msaData[[x]], sep=""))
    sequences[[x]] <- gene.info[[1]]$NTSEQ
    widths[[x]] <- width(gene.info[[1]]$NTSEQ)
  }

  # De alignment wordt uitgevoerd
  alignment <- msa(unlist(sequences))

  # De naam en sequentie van het langste gen worden opgeslagen
  longestGene <- which(widths == max(widths))
  longestGene.name <- msaData[longestGene]
  longest.seq <- sequences[longestGene]
  
  # Het percentage conservering van het langste gen wordt berekend 
  # En de resultaten worden weggeschreven
  perc <- getConservedPercentage(longestGene, alignment)
  maPerc <- as.matrix(perc)
  colnames(maPerc) <- c(longestGene.name)
  taWriteRes <- as.table(maPerc)
  write.table(taWriteRes, row.names=FALSE, file="MSARESULT.txt")
}

getEntrez <- function(dfGenes){
  #De functie splitst de genen op een punt en geeft de 
  #Entrez ID's terug
  entrez <- unlist(strsplit(dfGenes, "[.]"))
  entrez <- unique(entrez)
  return(entrez)
}

getConservedPercentage <- function(longest.gene, alignment){
  #De consensus van de alignment wordt genomen en
  #De 
  alignment.matrix <- consensusMatrix(unmasked(alignment))
  alignment.unmasked <- unmasked(alignment)
  consensus <- gsub("[?]","", consensusString(alignment.matrix))
  consensus.vector <- unlist(strsplit(consensus, split=""))
  
  # De lengte van de sequentie van het langste gen 
  # en van de alignment worden geteld. Als twee tekens overeen komen,
  # worden deze opgeslagen. Hieruit wordt een percentage berekend.
  # Het percentage wordt teruggegeven.
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

# Additional information:
# =======================
# Dit script is geschreven door Jesse Kerkvliet en
# kan gebruikt worden bij de opdracht van bapgc door
# Sebastiaan de Vriend.
#
# Genereert een MSA en berekent het percentage conservering van 
# Het langste gen
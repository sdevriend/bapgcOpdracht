#!/usr/bin/Rscript


library(Biostrings)
library("BSgenome.Hsapiens.UCSC.hg38")
library("GenomicRanges")
library(KEGGREST)
library("msa")
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


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
  fctMsaData <- fullFrame[c(1:50),]
  fctMsaData <- getEntrez(fctMsaData)

  # De sequenties en de lengten van de genen worden in een DNAStringSetList
  # en een vector opgeslagen
  sequences <- DNAStringSetList()
  vctWidths <- c()
  for(x in seq(1,length(fctMsaData))){
    gene.info <- keggGet(paste("hsa:",fctMsaData[[x]], sep=""))
    sequences[[x]] <- gene.info[[1]]$NTSEQ
    vctWidths[[x]] <- width(gene.info[[1]]$NTSEQ)
  }

  # De alignment wordt uitgevoerd
  alignment <- msa(unlist(sequences))

  # De naam en sequentie van het langste gen worden opgeslagen
  longestGene <- which(vctWidths == max(vctWidths))
  longestGene.name <- fctMsaData[longestGene]
  longest.seq <- sequences[longestGene]
  
  # Het percentage conservering van het langste gen wordt berekend 
  # En de resultaten worden weggeschreven
  perc <- getConservedPercentage(longestGene, alignment)
  maPerc <- as.matrix(perc)
  colnames(maPerc) <- c(longestGene.name)
  taWriteRes <- as.table(maPerc)
  write.table(taWriteRes, row.names=FALSE, file="MSARESULT.txt")
}

getEntrez <- function(vctGenes){
  #De functie splitst de genen op een punt en geeft de 
  #Entrez ID's terug
  vctEntrez <- unlist(strsplit(as.character(vctGenes), "[.]"))
  vctEntrez <- unique(vctEntrez)
  return(vctEntrez)
}

getConservedPercentage <- function(longest.gene, alignment){
  #De consensus van de alignment wordt genomen en
  #De ?s worden verwijderd
  mtAlignment <- consensusMatrix(unmasked(alignment))
  alignment.unmasked <- unmasked(alignment)
  consensus <- gsub("[?]","", consensusString(mtAlignment))
  vctConsensus <- unlist(strsplit(consensus, split=""))
  
  # De lengte van de sequentie van het langste gen 
  # en van de alignment worden geteld. Als twee tekens overeen komen,
  # worden deze opgeslagen. Hieruit wordt een percentage berekend.
  # Het percentage wordt teruggegeven.
  NT.count<-length(vctConsensus)
  longest.seq <- alignment.unmasked[longest.gene]
  longest <- gsub("[?]","", longest.seq)
  vctLongest <- unlist(strsplit(longest, split=""))
  longest.count <- length(vctLongest)
  count <- 0
  for(y in seq(1, length(vctConsensus))){
    if(vctConsensus[[y]] == vctLongest[[y]]){
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
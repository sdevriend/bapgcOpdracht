#!/usr/bin/Rscript


library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


showUsageInformation <- function(){
  print("") 
  print("Genereert een barplot voor de hydrofobiciteit van de 10 beste")
  print("medegereguleerde genen. De aminozuren worden verdeeld in")
  print("de categorieen: hydrofoob, hydrofiel en neutraal.")
  print("Het script kan aangeroepen worden door Hydrofobiciteit.R")
  print("")
  quit()
}

PreferenceRatio <- function(){
  #Gene database en genoom worden gedefinieerd
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genome <- BSgenome.Hsapiens.UCSC.hg38

  #Alle gevonden genen worden geladen en de eerste 10 worden meegenomen
  #De ID's worden omgezet naar Entrez Gene ID's 
  load("AllCoregulatedGenes.csv")
  vctProtein.data <- fullFrame[c(1:10),]
  vctProtein.data <- getEntrez(vctProtein.data)

  # Alle exons van de genen worden geladen,
  # hiervan wordt de sequentie geladen
  # De sequenties worden getransleerd naar AA sequenties
  exons <- exonsBy(txdb, by="gene")[vctProtein.data]
  protein.dna <- getSeq(genome, exons)
  Protein.AA <- translate(unlist(protein.dna))

  #De hydrofobiciteit wordt geteld en de barplot wordt gegenereerd
  dfAmounts <- getRatios(Protein.AA)
  generateBarplot(Protein.AA, dfAmounts)
}

getEntrez <- function(vctGenes){
  #De functie splitst de genen op een punt en geeft de 
  #Entrez ID's terug
  vctEntrez <- unlist(strsplit(as.character(vctGenes), "[.]"))
  vctEntrez <- unique(vctEntrez)
  return(vctEntrez)
}

getRatios <- function(AAStringSet){
  #Een set van unieke ID's wordt gegenereerd en een lege matrix wordt gemaakt
  IDs.Unique <- unique(names(AAStringSet))
  dfAmounts <- matrix(ncol=3)
  #Voor ieder ID in de unieke ID list wordt een set van eiwitsequenties
  #gemaakt, horend bij dat ID.
  #Ieder aminozuur wordt geteld. De set van hydrofobische AA's wordt geteld
  #en de set van hydrofiele AA's wordt geteld.
  #Een volledige som van het eiwit wordt berekend. Alle sets worden
  #in de matrix opgeslagen. Deze wordt teruggegeven.
  for(x in seq(1, length(IDs.Unique))){
    ID <- IDs.Unique[x]
    ID.indices <- which((names(AAStringSet) == ID)== T)
    ID.AASet <- AAStringSet[ID.indices]
    ID.counts <- alphabetFrequency(ID.AASet)
    df.counts <- as.data.frame(ID.counts)
    hydrophobic <- sum(df.counts$A, df.counts$C, df.counts$F, df.counts$G,
                       df.counts$H, df.counts$I, df.counts$K, df.counts$L,
                       df.counts$M, df.counts$S, df.counts$T, df.counts$V,
                       df.counts$W, df.counts$Y)
    hydrophilic <- sum(df.counts$D, df.counts$E, df.counts$K, df.counts$R)
    total <- sum(df.counts)
    hydrophilic <- hydrophilic/total*100
    hydrophobic <- hydrophobic/total*100
    total <- 100
    dfAmounts <- rbind(dfAmounts, c(hydrophobic, hydrophilic, 
		(total-hydrophobic-hydrophilic)))
    }
  return(dfAmounts[-1,])
}

generateBarplot <- function(Protein.AA, dfAmounts){
  # De aantallen matrix wordt omgedraaid (rijen worden kolommen en v.v)
  # Een figuur wordt gemaakt, en een plot wordt gegenereerd en weggeschreven
  dfAmounts.t <- t(dfAmounts)
  png("Hydrofobiciteit.png")
  rownames(dfAmounts.t) <- c("Hydrophoob","Hydrofiel","Totaal")
  colnames(dfAmounts.t) <- unique(names(Protein.AA))
  barplot(dfAmounts.t, 
	main="Verdeling hydrofobiciteit in gecoreguleerde genen",
    xlab="Entrez ID", ylab="Frequentie",
	col=c("darkblue","blue","lightblue"),
    legend = rownames(dfAmounts.t), beside=FALSE)
  dev.off()
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
      PreferenceRatio()
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
# Een barplot van de hydrofobiciteit van de beste 10 medegereguleerde 
# genen wordt gegenereerd en opgeslagen in een afbeelding.

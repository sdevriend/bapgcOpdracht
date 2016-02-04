#!/usr/bin/Rscript
library("hgu95av2.db")
library(limma)
Enzymes <- function{
  load("AllCoregulatedGenes.csv")

  # Neemt de honderd beste hits van een genenset (in dit geval een testvar)
  HundredBestHits <- fullFrame[c(1:100),]
  geneset <- getEntrez(HundredBestHits)
  # Genereert annotaties (Entrez, Symbol, Ensembl, functie en EC-nummer)
  annotations <- getEnzymes(unlist(HundredBestHits))
  # Haalt alle annotaties met een EC-nummer op
  enzymes <- annotations[!is.na(annotations$ENZYME),]
  # Zoekt naar hypotheticals, unknowns en uncharacterizeds.
  funcless <- getEnzymesWithKnownFunctions(annotations)
  # Genereert venn diagram
  generateVennDiagram(enzymes)
}

getEnzymes <- function(geneset){
  #De entrezID's worden genomen door de set te splitsen op een punt
  entrez <- unlist(strsplit(geneset, "[.]"))
  entrez <- unique(entrez)
  #De duplicaten worden verwijderd
  #De annotaties worden geselecterd
  annotation<- select(hgu95av2.db, key=entrez, keytype="ENTREZID", 
                      columns=c("SYMBOL", "ENSEMBL","GENENAME","ENZYME", 
                                "ENTREZID"))
  #Geeft de annotatie terug
  return (annotation)
}

getEnzymesWithKnownFunctions<-function(enzymeset){
    #Haalt alle enzymen zonder functie op
    noFuncs <- grepl("uncharacterized | hypothetical | ^unknown$", enzymeset$GENENAME)
    #Er wordt gekeken naar de functieloze enzymen. Deze worden opgeslagen in een
    #lijst. De lijst wordt ontdaan van NA's en teruggegeven.    
    noFuncs2 <- c()
    for(i in seq(1,(length(noFuncs)))){
    
    if(noFuncs[i] == TRUE){
      print(enzymeset$SYMBOL[i])
      noFuncs2[[i]] <- enzymeset$SYMBOL[i]
    }
  }
  noFuncs2 <- noFuncs2[!is.na(noFuncs2)]
  return(unique(noFuncs2))
  
}


generateVennDiagram <- function(enzymeset){
  #De functieloze enzymen worden gezocht
  noFuncs <- grepl("uncharacterized | hypothetical | ^unknown$", enzymeset$GENENAME)
  
  #Een logical object wordt gemaakt met alle enzymen die een EC-nummer hebben
  #(dat zijn ze allemaal)
    totalEnzyms <- (!is.na(enzymeset$ENZYME))
  #De twee logische klassen worden vergeleken en een venndiagram wordt gemaakt
  c3 <- cbind(totalEnzyms, noFuncs)
  count <- vennCounts(c3)
  vennDiagram(count, mar=c(1,1,1,1),
              names=c("Enzymen","Enzymes zonder bekende functie"), 
              main="Enzymen zonder functie in de set",
              lwd=2, circle.col="blue", cex=1)
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
      Enzymes()
    }
  }
  
}


main(commandArgs(T))

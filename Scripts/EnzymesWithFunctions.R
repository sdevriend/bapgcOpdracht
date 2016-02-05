#!/usr/bin/Rscript


library("hgu95av2.db")
library(limma)
library(png)


showUsageInformation <- function()
{
  print("") 
  print("Het script genereert een venndiagram")
  print("Aan de hand van de honderd beste medegereguleerde genen.")
  print("Het diagram toont de hoeveelheid enzymen in de set")
  print("en de hoeveelheid enzymen zonder bekende functie")
  print("Het script kan aangeroepen worden door EnzymesWithFunctions.R")
  print("")
  quit()
}

Enzymes <- function(){
	# Alle gevonden medegereguleerde genen worden ingeladen.
	# De honderd bovenste hits van de lijst worden genomen en de 
	# ID's worden omgezet naar Entrez-gene ID's
	# De annotaties worden opgehaald. Hierin zitten de Entrez-gene ID,
	# functie en EC-nummer. Alle annotaties met een EC-nummer worden 
	# opgehaald. De hypothetical, unknowns en uncharaterizeds worden
	# gezocht en verwijdert en de functie om een venndiagram te genereren
	# wordt aangeroepen.
	load("AllCoregulatedGenes.csv")
	vctHundredBestHits <- fullFrame[c(1:100),]
	vctGeneset <- getEntrez(vctHundredBestHits)
	dfAnnotations<- select(hgu95av2.db, key=vctGeneset, keytype="ENTREZID", 
                      columns=c("GENENAME","ENZYME", 
                                "ENTREZID"))
	dfEnzymes <- dfAnnotations[!is.na(dfAnnotations$ENZYME),]
	generateVennDiagram(dfEnzymes)
}

getEntrez <- function(dfGenes){
    #De functie splitst de ID's op een punt en de unieke lijst
    #Wordt teruggegeven
    vctEntrez <- unlist(strsplit(as.character(dfGenes), '\\.'))
    vctEntrez <- unique(vctEntrez)
    return(vctEntrez)
}

generateVennDiagram <- function(dfEnzymeset){
    #De functieloze enzymen worden gezocht,
	#Een logical object wordt gemaakt met alle enzymen die een EC-nummer hebben
	#(dat zijn ze allemaal, het typeert een enzym)
	#De twee logische klassen worden vergeleken en een venndiagram wordt gemaakt
	dfNoFuncs <- grepl("uncharacterized | hypothetical | ^unknown$",
        dfEnzymeset$GENENAME)
	dfTotalEnzyms <- (!is.na(dfEnzymeset$ENZYME))

	mtDelen <- cbind(dfTotalEnzyms, dfNoFuncs)
	count <- vennCounts(mtDelen)
	png("VennDiagram.png")
	vennDiagram(count, mar=c(1,1,1,1),
        names=c("Enzymen","onbekende functie"), 
        main="Enzymen zonder functie in de set",
        lwd=2, circle.col="blue", cex=1)

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
      Enzymes()
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
# Het script genereert een venndiagram met enzymen zonder bekende functie.
# en slaat deze op als afbeelding

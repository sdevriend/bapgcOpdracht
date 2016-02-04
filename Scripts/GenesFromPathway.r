#!/usr/bin/Rscript


library("KEGGREST")
library("png")


showUsageInformation <- function()
{
  print("") 
  print("Het script download informatie van een pathway.")
  print("De informatie wordt in een csv file gezet met naam pathway")
  print(" en de genen die op een pathway zitten.")
  print("Daarnaast haalt het script de afbeelding van de meegegeven")
  print(" pathwaynaam op van Kegg.")
  print("")
  print("Het script kan aangeroepen worden door ./GenesFromPathway pathwaynaam.")
  print("")
  quit()
}

infoFromPathway <- function(strPWNAME){
  # Functie haalt informatie van pathway op en zet deze in een csv file.
  #
  # lstPathwayInfo: Pathway object van kegg
  # strPathwayName: String met pathway naam
  # strGeneIDs: Char string. 1e element wordt genomen. 
  #     2e element wordt overgeslagen.
  # dfGenes: Omgezette geneid's.
  # png: Variable met pathway afbeelding.
  #
  # Na het omzetten van de geneid's worden deze in een dataframe gestopt
  # zodat de pathwaynaam als kolom meegegeven kan worden.
  # Het bestand wordt opgeslagen als csv bestand.
  # Daarnaast wordt de afbeelding van de pathway gedownload en
  # opgeslagen.
  lstPathwayInfo <- keggGet(strPWNAME)
  strPathwayName <- lstPathwayInfo[[1]]$NAME
  strGeneIDs <- lstPathwayInfo[[1]]$GENE[c(TRUE, FALSE)]
  png <- keggGet(strPWNAME, "image")
  dfGenes <- as.data.frame.character(strGeneIDs)
  colnames(dfGenes) <- strPathwayName
  write.table(dfGenes[1], "PathwayInfo.csv", row.names=FALSE)
  writePNG(png, paste(strPWNAME, ".png", sep=""))
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
      infoFromPathway(args[2])
    }
  }
}

main(commandArgs(T))

# Additional information:
# =======================
# Dit script is geschreven door Sebastiaan de Vriend en
# kan gebruikt worden bij de opdracht van bapgc door
# Jesse Kerkvliet
#
# Het script download de gegevens van een opgegven pathway
# en zet deze in de opgegeven directory in args 1.

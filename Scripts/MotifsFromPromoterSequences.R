#!/usr/bin/Rscript


source("https://bioconductor.org/biocLite.R")
biocLite()# Installation bioclite
library(JASPAR2014)
library(TFBSTools)
library(Biostrings)
showUsageInformation <- function()
{
  print("") 
  print("Het script download informatie van een pathway.")
  print("De informatie wordt in een csv file gezet met naam pathway")
  print(" en de genen die op een pathway zitten.")
  print("")
  print("Het script kan aangeroepen worden door GenesFromPathway pathwaynaam.")
  print("")
  quit()
}

getMotifsFromPromoterSeqs <- function(){
  # Maakt de parameters voor de volgende stap
  opts = list()
  opts[["species"]] <- 9606
  opts[["matrixtype"]] <- "PWM"
  
  # Maakt motifList m.b.v JASPAR2014. De opties voor mens en PWM worden
  # Meegegeven
  motifList <- getMatrixSet(JASPAR2014, opts)
  
  # De motifs worden in de promotorsequenties opgezocht
  seqs1 <- unlist(promoseqs)
  test<- searchSeq(motifList, seqs1, min.score="90%") # DUURT LANG!
  ?searchSeq
  # Schrijft tussenstap weg naar csv file
  write.table(test, "MotifsFromPromoterSeqs.csv", row.names=FALSE)
  dfTest <- as(test, "DataFrame") # DUURT LANG!
  write.table(dfTest, "Boekenlegger-25-1.csv", row.names=FALSE)
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
      print ("NOT FINISHED! DO NOT EXECUTE!")
      infoFromPathway()
    }
  }
  else
  {
    infoFromPathway() #test var
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
# en zet deze in een nog default script plek.

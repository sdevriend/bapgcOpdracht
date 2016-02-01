#!/usr/bin/Rscript

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
  load("promoseqs.RData")
  amount <- length(promoseqs)
  seqs1 <- unlist(promoseqs)
  # De genen worden langs de motiflijst gehaald. Threshold 80%
  sequences<- searchSeq(motifList, seqs1, min.score="80%") # DUURT LANG!

  # Maakt een dataframe van de searchSeq output
  dfSequences <- as(sequences, "DataFrame") # DUURT LANG!
  
  # Haalt de namen van de gevonden motifs uit het dataframe
  foundMotifs <- getMotifs(dfSequences, amount)
  write.csv(foundMotifs, "FoundMotifs.csv")
  return(foundMotifs)
}

getMotifs <- function(dfMotifs, amount){
  found = c()
  # De motifs die een relatieve score van 95% hebben en voor ieder gen gelden
  # Worden opgeslagen in een vector. Deze wordt teruggegeven.
  for (i in seq(1, (length(unique(dfMotifs$ID))))){
    setMotifs <- subset(dfMotifs, ID==(unique(dfMotifs$ID)[i]) & relScore > 0.95)
    uniek <- unique(setMotifs$seqnames)
    if (length(uniek) == amount){
      found[[i]] <- (unique(dfMotifs$ID)[i])
      
    }
  }
  found <- found[!is.na(found)]
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
      getMotifsFromPromoterSeqs()
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
# en zet deze in een nog default script plek.

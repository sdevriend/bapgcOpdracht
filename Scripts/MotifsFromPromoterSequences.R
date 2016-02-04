#!/usr/bin/Rscript


library(Biostrings)
library(JASPAR2014)
library(TFBSTools)


showUsageInformation <- function()
{
  print("") 
  print("Het script zoekt motifs bij bijbehorende promotorsequenties die")
  print(" geladen worden uit de R datafile promoseqs. De motifs moeten")
  print(" minimaal een score van 80% hebben om opgenomen te worden")
  print(" in de motiflist. Van de motiflist wordt een dataframe ")
  print(" gemaakt en deze wordt gefilterd op unieke motifs.")
  print(" De gevonden motifs worden opgeslagen in FoundMotifs.csv.")
  print(" Het script kan als volgt worden aangeroepen: ")
  print("./MotifsFromPromoterSequences.R WerkMap")
  print("")
  quit()
}

getMotifsFromPromoterSeqs <- function(){
  # De functie laad een optielijst in waar PWM en de JASPAR2014
  # nummer van de mens in gestopt worden. De lijst wordt gebruikt
  # om een lijst met alle motifs op te halen. 
  # Ook wordt de variable promoseqs geladen en daarvan wordt 
  # bepaald hoeveel sequenties erin zitten. Daarna wordt er een
  # unlist uitgevoerd die gebruikt wordt door searchSeq waar
  # motifs eruit komen met een score van minimaal 80%. Dit is gedaan
  # omdat er veel resultaten zijn, en er nog zat gefilterd kan worden
  # op deze manier op het resultaat.
  #
  # De sequences en motifs worden vervolgens in een dataframe gezet en
  # wordt de dataframe doorgestuurd naar getMotifs. Het resultaat
  # daarvan wordt opgeslagen in FoundMotifs.csv.
  # Maakt de parameters voor de volgende stap
  lstOpts = list()
  lstOpts[["species"]] <- 9606 
  lstOpts[["matrixtype"]] <- "PWM"
  # Maakt motifList m.b.v JASPAR2014. De opties voor mens en PWM worden
  # Meegegeven.
  motifList <- getMatrixSet(JASPAR2014, lstOpts)
  # De motifs worden in de promotorsequenties opgezocht.
  load("promoseqs.RData")
  intAmount <- length(promoseqs)
  seqs1 <- unlist(promoseqs)
  # De genen worden langs de motiflijst gehaald. Threshold 80%
  sequences <- searchSeq(motifList, seqs1, min.score="80%")
  # Maakt een dataframe van de searchSeq output
  dfSequences <- as(sequences, "DataFrame")
  # Haalt de namen van de gevonden motifs uit het dataframe
  foundMotifs <- getMotifs(dfSequences, intAmount)
  write.csv(foundMotifs, "FoundMotifs.csv")
  return(foundMotifs)
}

getMotifs <- function(dfMotifs, intAmount){
  # dfMotifs: Dataframe met gevonden motieven en scores per gevonden hit.
  # intAmount: Integer: Lengte van aantal PromoterSequenties.
  # De functie kijkt per subset motifs of de score hoger is dan 95%
  # daarna worden alle unieke motifs toegevoegd aan de vector found.
  # Daarna worden alle unieke motifs terug gegeven.
  found = c()
  # De motifs die een relatieve score van 95% hebben en voor ieder gen gelden
  # Worden opgeslagen in een vector. Deze wordt teruggegeven.
  for (i in seq(1, (length(unique(dfMotifs$ID))))){
    setMotifs <- subset(dfMotifs, ID==(unique(dfMotifs$ID)[i]) & relScore > 0.95)
    uniek <- unique(setMotifs$seqnames)
    if (length(uniek) == intAmount){
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
# Dit script is geschreven door Jesse Kerkvliet en
# kan gebruikt worden bij de opdracht van bapgc door
# Sebastiaan de Vriend.
#
# Het script zoekt motifs bij bijbehorende promotorsequenties die 
# geladen worden uit de R datafile promoseqs.
# Unieke motifs worden opgeslagen in: FoundMotifs.csv.
# Het script kan als volgt worden aangeroepen:
# ./MotifsFromPromoterSequences.R WerkMap

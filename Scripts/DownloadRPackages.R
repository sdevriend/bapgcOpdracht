#!/usr/bin/Rscript

showUsageInformation <- function(){
  print("")
  print("Dit script installeerd alle packages die gebruikt worden voor de bapgc opdracht.")
  print("Hiervoor is een internet verbinding nodig.")
  quit()
}

installdata <- function(){
  source("https://bioconductor.org/biocLite.R")
  biocLite()# Installation bioclite
  biocLite("KEGGREST")# installation keggrest?
  biocLite("GenomicRanges")
  biocLite("BSgenome.Hsapiens.UCSC.hg38")
  biocLite("org.Hs.eg.db")
  biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
  biocLite("JASPAR2014")
  #biocLite("gsl") # nodig voor TFBSTOOLS?
  biocLite("TFBSTools")
  biocLite("hgu95av2.db")
  biocLite("limma")
  #install.packages("foreach")
  #install.packages("doParallel")
  #install.packages("plyr")


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
      
    }
    
  }
  else
  {
    installdata()
  }
}


main(commandArgs(T))
# Additional information:
# =======================
# Dit script is geschreven door Sebastiaan de Vriend en
# kan gebruikt worden bij de opdracht van bapgc door
# Jesse Kerkvliet
#
# Het script download de benodigde packages die gebruikt
# worden bij de bapgc opdract.

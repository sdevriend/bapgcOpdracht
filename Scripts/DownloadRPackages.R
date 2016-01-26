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

#!/usr/bin/Rscript


showUsageInformation <- function(){
  print("")
  print("Dit script installeerd alle packages die gebruikt worden voor de ")
  print("bapgc opdracht.")
  print("Hiervoor is een internet verbinding nodig samen met root rechten.")
  quit()
}

installPackages <- function(){
  # Functie download alle onderstaande packages. Dit is handig voor bij
  # een schone installatie. De functie heeft root rechten nodig
  # om alle files te installeren.
  source("https://bioconductor.org/biocLite.R")
  biocLite()
  biocLite("KEGGREST")
  biocLite("GenomicRanges")
  biocLite("BSgenome.Hsapiens.UCSC.hg38")
  biocLite("org.Hs.eg.db")
  biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
  biocLite("JASPAR2014")
  biocLite("TFBSTools")
  biocLite("hgu95av2.db")
  biocLite("limma")
  biocLite("msa")
  biocLite("biomaRt")
  biocLite("seqinr")
  biocLite("ape")
  biocLite("ggtree")
}

main <- function(args)
{
  if (length(args) > 0)
  {
    if (args[1] == "-h" |  args[1] == "-help" | args[1] == "--h" | args[1] == "--help")
    {
      showUsageInformation()
    }
  }
  else
  {
    installPackages()
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

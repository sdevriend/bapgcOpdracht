#!/usr/bin/Rscript


library(ape)
library(ggtree)
library(png)


showUsageInformation <- function()
{
  print("") 
  print("Het script laad een phylogenetische boom in en visualiseerd ")
  print(" deze. Het script kan aangeroepen worden met: ")
  print("./generateTree.R Werkmap")
  print("")
  quit()
}
generateTree <- function(){
  # De functie leest tree.ph in en visualiseerd deze met ggtree
  # De visualisatie wordt opgeslagen als tree.png
  phy <- read.tree("tree.ph")
  ggtree(phy, branch.length=20) + geom_tiplab(size=11)  
  ggsave(filename = "tree.png", width=35, height=20)
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
      generateTree()
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
# Het script maakt een visualisatie van tree.ph en slaat
# deze op als tree.png. Het script kan aangeroepen worden met:
# ./generateTree.R Werkmap

#!/usr/bin/Rscript
library(ape)
library(ggtree)
library(png)
generateTree <- function(){
  phy <- read.tree("tree.ph")
  
  #png('tree.png')
  ggtree(phy, branch.length=20) + geom_tiplab(size=11)  
 # dev.off()
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
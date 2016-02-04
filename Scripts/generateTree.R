#!/usr/bin/Rscript
setwd("/home/bapgc/bapgcOpdracht/temphsa04916/")
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
generateTree()

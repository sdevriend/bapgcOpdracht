library(ape)
library(ggtree)
library(png)
generateTree <- function(){
  #phy <- read.tree("tree.ph")
  
  png('tree.png')
  ggtree(phy, layout="circular") + geom_tiplab()  
  dev.off()
  ggsave(filename = "tree.png")
}
generateTree()

source("https://bioconductor.org/biocLite.R")
biocLite()# Installation bioclite
biocLite("KEGGREST")# installation keggrest?
library("KEGGREST")

png<- keggGet("hsa04920", "image")

 library(png)
 writePNG(png,"test.png")
 
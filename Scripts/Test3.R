source("http://bioconductor.org/biocLite.R")
biocLite("pathview")
library(pathview)
data(gse16873.d)
1+1
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa", out.suffix = "gse16873")
str(pv.out)
head(pv.out$plot.data.gene)
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                      species = "hsa", out.suffix = "gse16873.2layer", kegg.native = T,
                      same.layer = F)
 pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                      species = "hsa", out.suffix = "gse16873", kegg.native = F,
                      sign.pos = demo.paths$spos[i])
 #pv.out remains the same
   dim(pv.out$plot.data.gene)

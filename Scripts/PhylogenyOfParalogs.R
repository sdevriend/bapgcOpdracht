#Paralogs
library(biomaRt)
library("msa")
library(Biostrings)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library("hgu95av2.db")
library(KEGGREST)
library(seqinr)
library(ape)
getSequences <- function(id.names){
  
  annotation<- select(hgu95av2.db, key=c(id.names), keytype="ENSEMBL", 
                      columns=c("ENTREZID"))  
  annotation<-na.omit(annotation)
  entrez.ids <- annotation$ENTREZID
  
  
  blacklist <- c()
  for( i in seq(1,(length(entrez.ids)))){
    
    #blacklist[[i]] 
    black<- tryCatch({
      transcriptsBy(txdb, by="gene")[entrez.ids[i]]},
      error = function(err){
        return(entrez.ids[[i]])
      })
    if(class(black) != "GRangesList"){
      blacklist[[i]] <- black
    }
  }
  blacklist <- na.omit(blacklist)
  
  if(length(blacklist) > 0){
    for(x in seq(1,length(blacklist))){
      index <- which(entrez.ids == blacklist[x] )
      entrez.ids <- entrez.ids[-index]
    }
  }
  sequences <- DNAStringSetList()
  for(gene in seq(1,length(entrez.ids))){
    gene.info <- keggGet(paste("hsa:",entrez.ids[[gene]], sep=""))
    sequences[[gene]] <- gene.info[[1]]$NTSEQ
    
  }
  
  return(sequences)
  
}


filterParalogs <- function(paralogs){
  ensembl.IDs <- unique(paralogs$ensembl_gene_id)
  paralogs.matrix <- matrix(ncol=2)
  for(x in seq(1, length(ensembl.IDs))){
    ID <- ensembl.IDs[x]
    ID.indices <- which((paralogs$ensembl_gene_id == ID)== T)
    print(length(ID.indices))
    if(length(ID.indices) > 4){
      paralogs.filterd <-  as(paralogs[ID.indices[c(1:4)],],"matrix")
      paralogs.matrix <- rbind(paralogs.matrix, paralogs.filterd)
    }else{
      paralogs.filtered.small <- as(paralogs[ID.indices,], "matrix")
      
      paralogs.matrix <- rbind(paralogs.matrix, paralogs.filtered.small)
    }
  }
  return(paralogs.matrix[-1,])
}

getEntrez <- function(dfGenes){
  #De functie splitst de ID's op een punt en de unieke lijst
  #Wordt teruggegeven
  entrez <- unlist(strsplit(dfGenes, "[.]"))
  entrez <- unique(entrez)
  return(entrez)
}
PhyloTree <- function(){
setwd("/home/bapgc/bapgcOpdracht/temphsa04916/")
#Gene database en genoom worden gedefinieerd
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genome <- BSgenome.Hsapiens.UCSC.hg38
load("AllCoregulatedGenes.csv")
genes <- fullFrame[c(1:10),]
genes <- getEntrez(genes)

g_ensembldatabase = useMart("ensembl", dataset="hsapiens_gene_ensembl")
paralogs <- getBM(c("ensembl_gene_id", "hsapiens_paralog_ensembl_gene"),
                  filters="entrezgene", values=genes, g_ensembldatabase)




four.paralogs <- filterParalogs(paralogs)
four.paralogs[four.paralogs==""] <- NA
four.paralogs <- na.omit(four.paralogs)
names <- (as.data.frame(four.paralogs)$ensemble_gene_id)
dfParalogs <- as.data.frame(four.paralogs)
names <- unique(dfParalogs$ensembl_gene_id)
names <- as(names, 'character')


paralog.names <- unique(dfParalogs$hsapiens_paralog_ensembl_gene)
paralog.names <- as(paralog.names, 'character')
nameset <- c(paralog.names, names)
sequences <- getSequences(nameset)


makeTree <- function(nameset, sequences){
  alignment <- seqinr::as.alignment(length(nameset),nameset,unlist(sequences))
  distance <- dist.alignment(alignment, matrix=c("similarity","identity"))
  tree <- bionj(distance)
  return(tree)

}

tree <- makeTree(nameset, sequences)
plot(tree)
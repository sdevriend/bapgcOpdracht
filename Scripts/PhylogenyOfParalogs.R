#Paralogs
biocLite("TxDb.Hsapiens.UCSC.hg38.ensGene")
library(biomaRt)
library("msa")
library(Biostrings)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library("hgu95av2.db")

setwd("C:/Users/jesse/Documents/Bio-informatica/Jaar 3/Periode 2/Bapgc/")
#Gene database en genoom worden gedefinieerd
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genome <- BSgenome.Hsapiens.UCSC.hg38
paralogs <- getBM(c("ensembl_gene_id", "hsapiens_paralog_ensembl_gene"),
                                                       filters="entrezgene", values=genes, g_ensembldatabase)
g_ensembldatabase = useMart("ensembl", dataset="hsapiens_gene_ensembl")
load("AllCoregulatedGenes.csv")


genes <- fullFrame[c(1:10),]
genes <- getEntrez(genes)
four.paralogs <- filterParalogs(paralogs)
four.paralogs[four.paralogs==""] <- NA
four.paralogs <- na.omit(four.paralogs)
names <- (as.data.frame(four.paralogs)$ensemble_gene_id)
dfParalogs <- as.data.frame(four.paralogs)
names <- unique(dfParalogs$ensembl_gene_id)
names <- as(names, 'character')


paralog.names <- unique(dfParalogs$hsapiens_paralog_ensembl_gene)
paralog.names <- as(paralog.names, 'character')

sequences <- getSequences(paralog.names)
alignment <- msa(unlist(sequences))

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
      print(paste("Error in",entrez.ids[i]))
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

  transcripts <- transcriptsBy(txdb, by="gene")[entrez.ids]
  transcripts.dna <-getSeq(genome, transcripts)
  return(transcripts.dna)
}
getAvailableGenes <- function(geneset){
  
  for(i in seq(1,(length(geneset)))){
    
  }
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


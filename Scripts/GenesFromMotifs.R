library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

genome <- BSgenome.Hsapiens.UCSC.hg38
Txdb<- TxDb.Hsapiens.UCSC.hg38.knownGene
allGenes <- genes(Txdb)
allTranscripts <- transcriptsBy(Txdb, by="gene")[allGenes$gene_id]

foundMotifsTable <- read.csv("FoundMotifs.csv")
foundMotifs <- foundMotifsTable$x

dfTenth1 <- getMatches(1, 2500)
save(dfTenth1, file='motifspart1.RDATA')
load('motifspart1.RDATA')
geneset11 <- getGenes(dfTenth1)

dfTenth2 <- getMatches(2501,5000)
save(dfTenth2, file='motifspart2.RDATA')
load('motifspart2.RDATA')
geneset12 <- getGenes(dfTenth2)

dfTenth3 <- getMatches(5001, 7500)
save(dfTenth3, file='motifspart3.RDATA')
geneset3 <- getGenes(dfTenth3)

dfTenth4 <- getMatches(7501, 10000)
save(dfTenth4, file='motifspart4.RDATA')
geneset4 <- getGenes(dfTenth4)

dfTenth5 <- getMatches(10001, 12500)
save(dfTenth5, file='motifspart5.RDATA')
geneset5 <- getGenes(dfTenth5)

dfTenth6 <- getMatches(12501, 15000)
save(dfTenth6, file='motifspart6.RDATA')
geneset6 <- getGenes(dfTenth6)

dfTenth7 <- getMatches(15001, 17500)
save(dfTenth7, file='motifspart7.RDATA')
geneset7 <- getGenes(dfTenth7)

dfTenth8 <- getMatches(17501, 20000)
save(dfTenth8, file='motifspart8.RDATA')
geneset8 <- getGenes(dfTenth8)

dfTenth9 <- getMatches(20001, 22500)
save(dfTenth9, file='motifspart9.RDATA')
geneset9 <- getGenes(dfTenth9)

dfTenth10 <- getMatches(22501, 23779)
save(dfTenth10, file='motifspart10.RDATA')
geneset10 <- getGenes(dfTenth10)

fullFrame <- data.frame(c(geneset1,geneset10, geneset2, geneset3, 
                          geneset4, geneset5, geneset6, geneset7, geneset8, geneset9))
save(fullFrame, file="AllCoregulatedGenes.csv")

getMatches <- function(left, right){
  Transcripts <- allTranscripts[IRanges(start=left, end=right)]
  PromoterSeqs <- getPromoterSeq(Transcripts, genome, 3000, 100)
  opts = list()
  opts[["ID"]] <- foundMotifs
  opts[["matrixtype"]] <- "PWM"
  
  # Maakt motifList m.b.v JASPAR2014. De opties voor mens en PWM worden
  # Meegegeven
  motifListSubSet <- getMatrixSet(JASPAR2014, opts)
  
  AllSSeqs <- searchSeq(motifListSubSet, unlist(PromoterSeqs), min.score="80%")
  dfSseqs<- as(AllSSeqs, "DataFrame")
  return(dfSseqs)
}

getGenes <- function(dfMotifs){
  found <- c()
  for (i in seq(1, (length(unique(dfMotifs$seqnames))))){
    setMotifs <- subset(dfMotifs, seqnames==(unique(dfMotifs$seqnames)[i]) & relScore > 0.97)
    uniek <- unique(setMotifs$ID)
    if (length(uniek) == 8){
      found[[i]] <- (unique(dfMotifs$seqnames)[i])
      found <- found[!is.na(found)]
      
    }
  }  
  return(found)
}

biocLite("msa")
library("msa")
library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg38")
library(Biostrings)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


setwd("C:/Users/jesse/Documents/Bio-informatica/Jaar 3/Periode 2/Bapgc/")
#Gene database en genoom worden gedefinieerd
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genome <- BSgenome.Hsapiens.UCSC.hg38

#De lijst met gecoreguleerde genen wordt geladen
load("AllCoregulatedGenes.csv")

# De bovenste tien hits worden hiervan genomen (moet nog naar 50) en de ID's 
# Worden omgezet naar Entrez Gene ID's
msaData <- fullFrame[c(1:10),]
msaData <- getEntrez(msaData)

# Alle transcripts van de genen worden geladen,
# hiervan wordt de sequentie geladen
transcripts <- transcriptsBy(txdb, by="gene")[msaData]
transcripts.dna <- getSeq(genome, transcripts)

# Een DNASeqSet wordt gemaak en uniek gemaakt.
geneSeqSet <- unlist(transcripts.dna)
geneSeqSet.u <- unique(geneSeqSet)

# Het langste gen wordt opgeslagen voor later gebruik
longestGene <- geneSeqSet[[which(width(geneSeqSet) == max(width(geneSeqSet)))]]

# Een msa van de genen wordt uitgevoerd en opgeslagen als stringset
msa <- msa(geneSeqSet.u) # DUURT LANG!
myMaskedAlignment <- msa
rowM <- IRanges(start=1, end=2)
rowmask(myMaskedAlignment) <- rowM
msa.stringset <- unmasked(myMaskedAlignment)


getEntrez <- function(dfGenes){
  entrez <- unlist(strsplit(dfGenes, "[.]"))
  entrez <- unique(entrez)
  return(entrez)
}
testing <- unlist(transcripts.dna)

class(testing)
letterFrequency(testing[[1]], letters="ACGT", OR=0)
testing[[1]]


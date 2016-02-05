#!/usr/bin/Rscript


library(ape)
library(biomaRt)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library("hgu95av2.db")
library(KEGGREST)
library("msa")
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(seqinr)


showUsageInformation <- function()
{
  print("") 
  print("Het script genereert fasta bestanden")
  print("Aan de hand van de twintig beste medegereguleerde genen en")
  print("per gen hooguit vier paralogen")
  print("Het script kan aangeroepen worden door PhologenyOfParalogs.R")
  print("")
  quit()
}

getSequences <- function(id.names){
  #De functie haalt alle sequenties op uit de knownGene database.
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene  
  annotation<- select(hgu95av2.db, key=c(id.names), keytype="ENSEMBL", 
                      columns=c("ENTREZID"))  
  annotation<-na.omit(annotation)
  entrez.ids <- annotation$ENTREZID
  
  # Er wordt gekeken of ieder gen een bijbehorende sequentie heeft.
  # Is dit niet het geval, wordt het gen in de blacklist gestopt.
  # De genen in de blacklist worden verwijderd uit de set.
  blacklist <- c()
  for( i in seq(1,(length(entrez.ids)))){
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
  
  # Als er genen zonder sequentie zijn, worden deze verwijderd.
  if(length(blacklist) > 0){
    for(x in seq(1,length(blacklist))){
      index <- which(entrez.ids == blacklist[x] )
      entrez.ids <- entrez.ids[-index]
    }
  }
  # De sequenties worden voor ieder gen opgehaald
  sequences <- DNAStringSetList()
  for(gene in seq(1,length(entrez.ids))){
    gene.info <- keggGet(paste("hsa:",entrez.ids[[gene]], sep=""))
    sequences[[gene]] <- gene.info[[1]]$NTSEQ
    
  }
  # De sequenties worden teruggegeven
  return(sequences)
}

filterParalogs <- function(paralogs){
  #Voor ieder gen wordt gekeken hoe veel paralogen
  #aanwezig zijn. Als er meer dan vier paralogen worden
  #gevonden, worden alleen de eerste vier meegegeven.
  #Een matrix met hooguit vier paralogen per gen wordt teruggegeven.
  ensembl.IDs <- unique(paralogs$ensembl_gene_id)
  paralogs.matrix <- matrix(ncol=2)
  for(x in seq(1, length(ensembl.IDs))){
    ID <- ensembl.IDs[x]
    ID.indices <- which((paralogs$ensembl_gene_id == ID)== T)
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
  vctEntrez <- unlist(strsplit(as.character(dfGenes), '\\.'))
  vctEntrez <- unique(vctEntrez)
  return(vctEntrez)
}

getGenesforPhylo <- function(){
  #Gene database en genoom worden gedefinieerd
  #De beste 10 beste medegereguleerde worden geladen
  #En omgezet naar Entrez-gene ID's
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genome <- BSgenome.Hsapiens.UCSC.hg38
  load("AllCoregulatedGenes.csv")
  genes <- fullFrame[c(1:10),]
  genes <- getEntrez(genes)
  
  # Een biomart met de ensembl genen wordt aangemaakt.
  # Voor alle genen worden de 
  g_ensembldatabase = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  paralogs <- getBM(c("ensembl_gene_id", "hsapiens_paralog_ensembl_gene"),
                  filters="entrezgene", values=genes, g_ensembldatabase)
  # Voor ieder gen worden hooguit vier paralogen opgehaald.
  # De genen zonder paralogen worden verwijderd.
  # De namen worden opgeslagen van de oorspronkelijke genen
  # En de paralogen. De namen worden aan een vector met alleen
  # namen toegevoegd en de sequenties worden opgehaald.

  mtFourParalogs <- filterParalogs(paralogs)
  mtFourParalogs[mtFourParalogs==""] <- NA
  mtFourParalogs <- na.omit(mtFourParalogs)
  vctNames <- (as.data.frame(mtFourParalogs)$ensemble_gene_id)
  dfParalogs <- as.data.frame(mtFourParalogs)
  vctNames <- unique(dfParalogs$ensembl_gene_id)
  vctNames <- as(vctNames, 'character')
  paralog.names <- unique(dfParalogs$hsapiens_paralog_ensembl_gene)
  paralog.names <- as(paralog.names, 'character')
  vctNameset <- c(paralog.names, vctNames)
  sequences <- getSequences(vctNameset)
  
  #Iedere sequentie wordt weggeschreven naar een fasta bestand.
  for(i in seq(1, length(sequences))){
  outname = paste(vctNameset[[i]],".fasta", sep="")
  write.fasta(unlist(sequences)[[i]], names = vctNameset[[i]], file.out = outname)
  }
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
      getGenesforPhylo()
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
# Het script genereert fasta bestanden voor de twintig best medegereguleerde
# genen en tot vier paralogen.

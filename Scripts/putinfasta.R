for(i in seq(1, length(sequences))){
  outname = paste(nameset[[i]],".fasta", sep="")
  write.fasta(unlist(sequences)[[i]], names = nameset[[i]], file.out = outname)
}
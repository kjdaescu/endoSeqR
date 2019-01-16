# Hello, endoSeqR!
# This is a way to test that the package installed
# which prints 'Hello, endoSeqR!'.


hello_endoSeqR <- function() {
  print("Hello, endoSeqR!")
}

#' Package_Load
#'
#' Loads dependencies
#' Package_Load()

Package_Load<-function(){
  library(openxlsx)
  library(plyr)
  library(Biostrings)
}

#' generate_fasta
#'
#' This function generates custom, library-specific IDs for a list of sequences from a library.
#' @param x - file of small RNA reads. It can have a "fastq/fq suffix", a "fasta/fa suffix", or a "txt" suffix with a list of sequences. Any file can be gzipped.
#' @param desired name of file. Output file will be a fasta file, so the suffix ".fasta" is automatticaly added (example: "Desired_name.fasta")
#' @export correctly formatted fasta file
#' generate_fasta()

generate_fasta<-function(x,library){
  y<-paste0(library)
  file_type<-substring(paste0(x),nchar(x),nchar(x))
  if (file_type %in% "z"){
    gunzip(paste0(x),overwrite=TRUE,remove=FALSE)
    file_type<-substring(paste0(x),nchar(x)-3,nchar(x)-3)
  }
  x<-readLines(x)
  x<-c(as.character((x)))
  file_spaces<-ifelse(file_type %in% "q", 4,
                      ifelse(file_type %in% "a", 2,
                             ifelse(file_type %in% "t", 1,NA)))
  start<-ifelse(file_type %in% "q", 2,
                ifelse(file_type %in% "a", 2,
                       ifelse(file_type %in% "t", 1,
                              ifelse(file_type %in% "t", 1,NA))))
  x<-x[seq(start, length(x), file_spaces)]
  x<-as.data.frame(x)
  x$sequence<-as.character(x[,1])
  x$sequence<-gsub("U","T",x$sequence)
  x$number<-seq(1,nrow(x))
  x$length<-nchar(as.character(x[,1]))
  print("Read length distribution:")
  print(count(x$length))
  x$ID<-paste0(">",library,"_",x$number)
  x<-as.vector(rbind(x$ID,as.character(x[,1])))
  print("Starting reads:")
  print(head(x))
  writeLines(x, paste0(y,".fasta"))
}

#' Map_endosiRNA
#'
#' This is internal function to get the bed files containing siRNA and it's genomic coordinates and ensuring length is 18-28 nts.
#' @param fasta_file_name is the fasta file
#' @param index is the species prefix to direct bowtie to the correct index
#' @export linux generates 2 files: file.siRNA.bed and file.siRNA.fasta.table
#' Map_endosiRNA()


Map_endosiRNA<-function(index, fasta_file_name){
      system(paste0("endo_MAP ",index," ",fasta_file_name),show.output.on.console=TRUE,intern=FALSE)
    }

#' Import_endosiRNA
#'
#' This function produces the output xlsx and RData files that contains all endo-siRNA information
#' @param fasta_file_name
#' @export xlsx bed file, detailed data frame, and RData detailed endo-siRNA file
#' Import_endosiRNA()

Import_endosiRNA<-function(fasta_file_name){
  fa<-paste0(fasta_file_name,".siRNA.fasta.table")
  fa<-read.table(fa)
  colnames(fa)<-c("Read_ID","Sequence")
  fa$Length<-nchar(as.character(fa$Sequence))
  fa<-fa[as.numeric(fa$Length) <= 28,]
  fa<-fa[as.numeric(fa$Length) >= 18,]
  print("Length siRNA")
  print(count(fa$Length))
  bed<-paste0(fasta_file_name,".siRNA.bed")
  bed<-read.table(paste0(bed))
  bed<-bed[,c(1:4,6)]
  colnames(bed)<-c("Chr","Start","Stop","Read_ID","Strand")
  fa$Read_ID<-substring(as.character(fa$Read_ID),2,nchar(as.character(fa$Read_ID)))
  bed<-bed[order(bed$Read_ID),]
  fa<-fa[order(fa$Read_ID),]
  df<-merge.data.frame(fa,bed)
  sequence<-DNAStringSet(as.character(df$Sequence))
  sequence<-reverseComplement(sequence)
  sequence2<-data.frame(sequence)
  sequence2<-unlist(sequence2[,1])
  df$target.mRNA<-sequence2
  save(df,file=paste0(fasta_file_name,"_endosiRNA.RData"))
  write.xlsx(df,file=paste0(fasta_file_name,"_detailed_endosiRNA.xlsx"),colNames=TRUE,rowNames=FALSE)
  write.xlsx(bed,file=paste0(fasta_file_name,"_endosiRNA_bed.xlsx"),colNames=TRUE,rowNames=FALSE)
  system(paste0("clean_up",""),show.output.on.console=TRUE,intern=FALSE)
  return(df)
}


#' PrepsiRNAntPlot
#'
#' This internal function generates leading base pair information
#' @param x - bed file with endosiRNA and sequences
#' PrepsiRNAntPlot()

PrepsiRNAntPlot<-function(x){
  x<-x[,c('Chr','Start','Stop','Sequence','Strand')]
  x<-unique(x)
  x$Leading.base<-substring(x$Sequence, 1, 1)
  y<-cbind.data.frame(nrow(x[x$Leading.base %in% "N",]),
                      nrow(x[x$Leading.base %in% "A",]),
                      nrow(x[x$Leading.base %in% "T",]),
                      nrow(x[x$Leading.base %in% "C",]),
                      nrow(x[x$Leading.base %in% "G",]))
  colnames(y)<-c("N", "A","U", "C", "G")
  #print("leading base distribution")
  #print(y)
  return(y)
}

#' HistBase1Plot
#'
#' This internal function generates leading base pair information for each size of endo-siRNA (18-28 base pairs)
#' @param x - bed file with endosiRNA and sequences
#' @export linux generates information for the leading base/length histogram and a data frame with the same information
#' HistBase1Plot()

HistBase1Plot<-function(x){
siRNA_lengths<-c(18:28)
x_size<-data.frame()
x_all<-data.frame()
for (i in 1:length(siRNA_lengths)){
  print(paste("endo-siRNA length: ", siRNA_lengths[i]))
  x_size<-x[nchar(x$Sequence) == siRNA_lengths[i],]
  x_size<-PrepsiRNAntPlot(x_size)
  x_all<-rbind.data.frame(x_all, x_size)
}
  z<-as.matrix(x_all)
  rownames(z)<-c(18:28)
  #print(count(z))
  return(z)
}

#' PlotNTLetter
#'
#' This function generates a histogram of endo-siRNA lengths and divides distribution by leading base
#' @param x - bed file with endosiRNA and sequences
#' PlotNTLetter()

PlotNTLetter<-function(x){
  x$Sequence<-as.character(x$Sequence)
  x<-HistBase1Plot(x)
  print(t(x))
  barplot(t(x),
          col=c("firebrick1","darkorchid3","darkturquoise","blue","darkseagreen1"),
          names=rownames(x),
          xlab="Length of small RNA (nucleotides)",
          main=paste0("endo-siRNA leading base distribution"),
          ylab="Leading Nucleotide")

  legend("topright",
         legend=c("G","C","U", "A"),
         col=c("darkseagreen1","blue","darkturquoise","darkorchid3"), pch=15)
  return(t(x))
}









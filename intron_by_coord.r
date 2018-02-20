##############################################################
# Extract introns for a list of genes using exon coordinates #
##############################################################
#if biomaRt and the corresponding genome packages are not downloaded yet
#install.packages('Rcpp')
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#biocLite('BSgenome.Mmusculus.UCSC.mm10')
#biocLite("BSgenome.Hsapiens.NCBI.GRCh38")
#or
#biocLite("BSgenome.Hsapiens.UCSC.hg38")

#the csv file with the genes of interest
#it should cointain a column with the ensembl gene ids and the corresponding header should be: Ensembl Gene ID
#the species: m for mouse, h for human
setwd("/Users/GG/Desktop/Kata/pipeline_python_part/pou5f1")
lst = 'pou5f1.csv'
sp = 'h'

library(biomaRt)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome)
library(Biostrings)
if(sp == 'm'){
  library(BSgenome.Mmusculus.UCSC.mm10)
  mart<-useEnsembl("ensembl","mmusculus_gene_ensembl")
}else if(sp == 'h'){
  library(BSgenome.Hsapiens.UCSC.hg38)
  mart<-useEnsembl("ensembl","hsapiens_gene_ensembl")
}

#import the gene list of interest
gene_list=read.delim(lst,sep=',',header = T)
gene_ids = gene_list$Ensembl.Gene.ID

#get the exon coordinates
exoncoord=getBM(attributes = c('ensembl_gene_id','ensembl_exon_id','chromosome_name','exon_chrom_start','exon_chrom_end'),
                filters='ensembl_gene_id',values=gene_ids,mart=mart,bmHeader = TRUE)
#sort by chromosome name and exon coordinate
g_ec=exoncoord[order(exoncoord$`Chromosome/scaffold name`),]
e_ec=g_ec[order(g_ec$`Exon region start (bp)`),]

#get the start and end cordinate for the genes
generange=getBM(attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position"),
                filters = "ensembl_gene_id",values = gene_ids,mart = mart,bmHeader = T)

#group by exon coordinate - to have overlapping exons
gname=e_ec$`Gene stable ID`[1]
e_strt = e_ec$`Exon region start (bp)`[1]
lst_end = e_ec$`Exon region end (bp)`[1]
for(gene in 1:dim(e_ec)[1]){
  #find the exon with the max end coordinate among overlapping ones
  e_end=e_ec$`Exon region end (bp)`[gene]
  if(gname==e_ec$`Gene stable ID`[gene]){
    e_strt =e_ec$`Exon region start (bp)`[gene]
    if(e_strt < lst_end){
      if(e_end > lst_end){
        lst_end = e_end
      }
    }else{
      while(e_strt > lst_end){
        #get intron seq with getSeq function by coordinates
        iend = e_ec$`Exon region start (bp)`[gene]
        intron_seq=getSeq(Hsapiens,paste("chr",e_ec$`Chromosome/scaffold name`[gene],sep = ""),
                          start=lst_end, end=iend)
        write(paste(paste(">",e_ec$`Gene stable ID`[gene],"_chr",e_ec$`Chromosome/scaffold name`[gene],":",
                          lst_end,"-",iend,sep=""),
                    as.character(intron_seq),sep="\n"),
              file = paste(as.character(gname),'.fa',sep=""),append = T,sep = "\n")
        write(paste(paste("chr",e_ec$`Chromosome/scaffold name`[gene],sep=""),lst_end,iend,e_ec$`Gene stable ID`[gene],sep="\t"),
              file="bed_version.txt",append=T,sep="\n")
        e_strt =e_ec$`Exon region start (bp)`[gene]
        lst_end=e_ec$`Exon region end (bp)`[gene]
        rm(intron_seq)
      }
    }
  }else{
    gname=e_ec$`Gene stable ID`[gene]
    e_strt = e_ec$`Exon region start (bp)`[gene]
    lst_end = e_ec$`Exon region end (bp)`[gene]
  }
}

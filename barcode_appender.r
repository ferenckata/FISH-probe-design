#barcode appender
fbc=read.csv('/Volumes/Macintosh HD 2/Common space/Kata/eleni_genes/mouse_barocdes/fwd_barcodes.csv',header = F)
rbc=read.csv('/Volumes/Macintosh HD 2/Common space/Kata/eleni_genes/mouse_barocdes/rev_barcodes.csv',header = F)
cbc=read.csv('/Volumes/Macintosh HD 2/Common space/Kata/eleni_genes/mouse_barocdes/col_barcodes.csv',header = F)
T7 = 'CGATTGAGGCCGGTAATACGACTCACTATAGGG'
f=1
r=1
k=1
m=1
setwd('/Volumes/Macintosh HD 2/Common space/Kata/eleni_genes/FINAL/all')
to=read.delim('40genes_to_order.fa', as.is=TRUE, header = FALSE)
ogn=''
frt=matrix(data=NA,nrow=31,ncol=31)
frt[1,2:31]=paste('F',seq(1,30),sep="")
frt[2:31,1]=paste('R',seq(1,30),sep="")
for(line in 1:dim(to)[1]){
  if(grepl('>',as.character(to[line,1]))){
    gi=unlist(strsplit(as.character(to[line,1]),'_',fixed = TRUE))
   if(line>1&&gi[2]!=ogn){
      r=r+1
      if(r==dim(rbc)[1]+1){
        f=f+1
        r=1
      }
      frt[r+1,f+1]=substr(gi[2],2,nchar(gi[2]))
      name=paste(fbc[f,1],'_T_',rbc[r,1],sep = '')
      seq=paste(fbc[f,3],to[line+1,1],rbc[r,3],sep='')
      note=paste('F:',fbc[f,2],'_T:',substr(gi[2],2,nchar(gi[2])),'_cds:',gi[10],'_R:',rbc[r,2],sep='')
      write(paste(name,seq,note,sep='\t'),'barcoded_40.tsv',append = TRUE)
    }else{
      frt[r+1,f+1]=substr(gi[2],2,nchar(gi[2]))
      name=paste(fbc[f,1],'_T_',rbc[r,1],sep = '')
      seq=paste(fbc[f,3],to[line+1,1],rbc[r,3],sep='')
      note=paste('F:',fbc[f,2],'_T:',substr(gi[2],2,nchar(gi[2])),'_cds:',gi[10],'_R:',rbc[r,2],sep='')
      write(paste(name,seq,note,sep='\t'),'barcoded_40.tsv',append = TRUE)
    }
    #store the gene name from the previous line
    ogn=gi[2]
  }
}
write.table(frt,'fwd_rev_table_for_40genes.tsv',sep = '\t',na='0',row.names = F,col.names = F,quote = F)
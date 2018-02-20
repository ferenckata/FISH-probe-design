#setwd("..")
f<-list.files(path=".",pattern="_co.fa")
max_thr=.6
min_thr=.4

#########################################################################
## calculating ##
#################

l<-length(f)
m<-matrix(data=NA,nrow=l,ncol=3)
for (file in 1:l){
  # saving the gene name #
  gn=unlist(strsplit(f[file],'.',fixed=TRUE))
  # total number of oligos with homopolymer stretch greater than 4 nt #
  hp=0
  data<-read.delim(f[file],header=FALSE,sep=":")
  
  hp=sum(data[,5],na.rm = T)
  cat(sprintf("In %s gene there are %d oligos in total.\nThere are %d oligos with homopolymer stretch greater than 4 nt.\n\n",
          gn[1], dim(data)[1]/2, hp))
  # saving the median of the GC and Tm values for each gene #
  m[file,1]=gn[1]
  m[file,2]<-median(data[,3],na.rm = T)
  m[file,3]<-median(data[,4],na.rm=T)
}

###########################################################################
## plotting ##
##############

# for multiple genes #
if(file>1){
  ## GC content ##
  png(filename='GC_median.png')
  barplot(as.numeric(m[,2]))
  abline(h=(summary(as.numeric(m[,2]))),col="green")
  abline(h=c(.4,.6),col="red")
  dev.off()
  cat("The quartiles of the GC content of the genes are:
  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. \n",
  summary(as.numeric(m[,2])),'\n\n')
  cat(sprintf("There are %d genes which have the median GC content below the threshold.\n",
          length(which(as.numeric(m[,2])<min_thr))))
  ab=m[which(as.numeric(m[,2])>.6),1]
  if(length(ab)>0){
    cat(ab)
  }
  cat(sprintf("There are %d genes which have the median GC content above the threshold.\n",
          length(which(as.numeric(m[,2])>max_thr))))
  ab=m[which(as.numeric(m[,2])<.4),1]
  if(length(ab)>0){
    cat(ab)
  }
  ## melting temperature ##
  png(filename='Tm_median.png')
  barplot(as.numeric(m[,3]))
  abline(h=(summary(as.numeric(m[,3]))),col="green")
  dev.off()
  cat("The quartiles of the Tm content of the genes are:
  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. \n",
  summary(as.numeric(m[,3])))
# for a single gene - in oligo resolution #
}else{
  ## GC content ##
  png(filename='GC_median.png')
  barplot(as.numeric(data[,3]))
  abline(h=summary(as.numeric(m[,2])),col="green")
  abline(h=c(.4,.6),col="red")
  dev.off()
  cat(sprintf("%d oligos have higher GC content than 0.6. \n%d oligos have optimal GC content. \n%d oligos have lower GC content than 0.4.\n\n",
          length(which(as.numeric(data[,3])>max_thr))/2,
          (length(which(as.numeric(data[,3])>=min_thr))+length(which(as.numeric(data[,3])<=max_thr))-1)/2,
          (length(which(as.numeric(data[,3])<min_thr))-1)/2))
  cat("The quartiles of the GC content of the oligos are:
  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. \n",
  summary(data[which(is.na(as.numeric(data[,3]))==F),3]),'\n\n')
  ## melting temperature ##
  png(filename='Tm_median.png')
  barplot(as.numeric(data[,4]))
  abline(h=summary(as.numeric(m[,3])),col="green")
  dev.off()
  cat("The quartiles of the Tm content of the oligos are:
  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. \n",
  summary(data[which(is.na(as.numeric(data[,4]))==F),4]))
}

#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Kata Ferenc
# Email: katagyik@gmail.com
# Version: 1.0
# Description: Retrieve cDNA or CDS sequences from biomart
# ------------------------------------------------------------------------------


# DEPENDENCIES =================================================================

#install.packages('argparser')
library(argparser)
library(biomaRt)

# INPUT ========================================================================

# Create argument parser
parser = arg_parser('cds_cdna_biomart', name = 'cds_cdna.r')

# Define mandatory arguments
parser = add_argument(parser, arg = 'din',
                      help = 'Directory where you have your gene list.')
parser = add_argument(parser, arg = 'f',
                      help = 'CSV file that contains the gene (Ensembl Gene or Transcript ID) list of interest.')
parser = add_argument(parser, arg = 'sq',
                      help = 'The type of sequence you want to download: cdna for cDNA, cds for CDS, b for both.')

# Define elective arguments
parser = add_argument(parser, arg = 'tp',
                      help = 'The ID type you would like to use: g for Ensembl Gene ID (default), t for Ensembl Transcript ID.',
                      default = 'g')
parser = add_argument(parser, arg = 'sp',
                      help = 'The species, h for human (default), m for mouse.',
                      default = 'h')
parser = add_argument(parser, arg = 'dout',
                      help = 'Directory where you would like to save your sequences.',
                      default = './')

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# RUN ==========================================================================

setwd(din)

 # Read gene list 
  gene_list<-read.csv(f,header=TRUE,sep=",")

  # Get the mart of the species
  if (sp == 'h'){
    mart<-useEnsembl("ensembl","hsapiens_gene_ensembl")
  }else if (sp == 'm'){
    mart<-useEnsembl("ensembl","mmusculus_gene_ensembl")
  }

 # Get sequences from Ensembl Biomart

 #----------------------------------------------  gene ID -----------------------------------------------------

  if (tp == 'g'){
    
    # Get the IDs
    genes<-as.character(gene_list$Ensembl.Gene.ID)
    
    for (i in 1:length(genes)){
      GOI<-genes[i]
      
      # Get all transcript information for the gene of interest
      atr<-getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id",
                                "transcript_length","cds_length"),
                 filters = "ensembl_gene_id",values = GOI,mart=mart,bmHeader = TRUE)
      
      # CDNA ### Get the sequence of the longest transcript (including CDS and UTRs)
      if (sq=='b'||sq=='cdna'){
        TOI<-atr$`Transcript stable ID`[which.max(atr$`Transcript length (including UTRs and CDS)`)]
        ltr<-getBM(attributes = c("hgnc_symbol","ensembl_gene_id","ensembl_transcript_id",
                                  "transcript_length","percentage_gene_gc_content",
                                  "cdna"),
                   filters="ensembl_transcript_id",values=TOI,mart=mart,bmHeader=TRUE)
        
        # Save the sequence
        write(paste(paste(paste(">",gene_list$Gene.Symbol[i],sep=""),ltr$'Gene stable ID'[1],ltr$'Transcript stable ID'[1],
                          "cdna",sep="_"),
              as.character(ltr$`cDNA sequences`),sep="\n"),
              file = paste(dout,gene_list$Gene.Symbol[i],'_cdna','.fa',sep=""),append = T,sep = "\n")
      }
      
      # CDS ### Get the sequence of the longest transcript (CDS without UTRs)
      if (sq=='b'||sq=='cds'){
        COI<-atr$`Transcript stable ID`[which.max(atr$`CDS`)]
        if (length(COI)==0){
          COI=TOI
        }
        lcds<-getBM(attributes = c("hgnc_symbol","ensembl_gene_id","ensembl_transcript_id",
                                   "percentage_gene_gc_content","coding"),
                    filters="ensembl_transcript_id",values=COI,mart=mart,bmHeader=TRUE)
        
        # Save the sequence
        write(paste(paste(paste(">",gene_list$Gene.Symbol[i],sep=""),lcds$'Gene stable ID'[1],lcds$'Transcript stable ID'[1],
                          "cds",sep="_"),
                    as.character(lcds$`Coding sequence`),sep="\n"),
              file = paste(dout,gene_list$Gene.Symbol[i],'_cds','.fa',sep=""),append = T,sep = "\n")
      }
      
      #clear variables
      rm(GOI,atr,TOI,ltr,COI,lcds)
    }
    
  # ---------------------------------------------------- transcript ID ----------------------------------------------

  }else if (tp == 't'){
    
      # Get the IDs
      genes<-as.character(gene_list$Ensembl.Transcript.ID)
      
      for (i in 1:length(genes)){
        TOI<-genes[i]
        
        # CDNA ### Get the sequence of the longest transcript (including CDS and UTRs)
        if (sq=='b'||sq=='cdna'){
          cdna<-getBM(attributes = c("hgnc_symbol","ensembl_gene_id","ensembl_transcript_id",
                                     "transcript_length","percentage_gene_gc_content",
                                     "cdna"),
                      filters="ensembl_transcript_id",values=TOI,mart=mart,bmHeader=TRUE)
          
          # Save the sequence
          f2<-paste(dout,paste(gene_list$Gene.Symbol[i],cdna$'Gene stable ID'[1],cdna$'Transcript stable ID'[1],
                               "cdna.csv",sep="_"),sep="")
          write.csv(cdna,f2,row.names = FALSE)
        }
        
        # CDS ### Get the sequence of the longest transcript (CDS without UTRs)
        if (sq=='b'||sq=='cds'){
          cds<-getBM(attributes = c("hgnc_symbol","ensembl_gene_id","ensembl_transcript_id",
                                    "percentage_gene_gc_content","coding"),
                     filters="ensembl_transcript_id",values=TOI,mart=mart,bmHeader=TRUE)
          
          # Save the sequence
          f3<-paste(dout,paste(tr$`HGNC symbol`[1],cds$'Gene stable ID'[1],
                                   cds$`Transcript stable ID`[1],"cds.csv",sep="_"),sep="")
          write.csv(cds,f3,row.names = FALSE)
        }
        
        #clear variables
        rm(TOI,cdna,f2,lcds,f3)
      }
    
    }
      

# END --------------------------------------------------------------------------

################################################################################
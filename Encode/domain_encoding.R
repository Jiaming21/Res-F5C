library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(AnnotationDbi)
library(m6ALogisticModel)
txdb <- loadDb("/Users/jiaming/Desktop/f5c/S288C.sqlite")
pos_C <- readRDS("/Users/jiaming/Desktop/f5c/f5c_positive_all.rds")

same_transcript <- subsetByOverlaps(transcripts(txdb),pos_C)
neg_C <- m6ALogisticModel::sample_sequence("C",same_transcript,Scerevisiae,Fixed = T)
neg_C <- subsetByOverlaps(neg_C,exons(txdb))
olp <- findOverlaps(neg_C,pos_C,ignore.strand= T)
neg_C <- neg_C[-queryHits(olp)]
neg_C <- neg_C[sample(seq_along(neg_C),10*length(pos_C))]
rm(olp,same_transcript)

GFgenreation_f5c <- function(data){
  analysis_data <- data
  matureSE <- SummarizedExperiment()
  rowRanges(matureSE) <- analysis_data 
  data_standardized <- predictors_annot(se = matureSE,
                                        txdb = txdb,
                                        bsgnm = Scerevisiae,
                                        motif_clustering = "C",
                                        isoform_ambiguity_method = "longest_tx",
                                        genes_ambiguity_method = "average",
                                        annot_clustering = matureSE,
                                        standardization = T) 
  GF <- mcols(data_standardized)
  return(GF)
}

pos_encoding <- GFgenreation_f5c(pos_C)
neg_encoding <- GFgenreation_f5c(neg_C)

transform <- function(data){
  
  data$UTR5[data$UTR5==FALSE] <- 0
  data$UTR5[data$UTR5==TRUE] <- 1
  
  data$UTR3[data$UTR3==FALSE] <- 0
  data$UTR3[data$UTR3==TRUE] <- 1
  
  data$cds[data$UTR3==FALSE] <- 0
  data$cds[data$UTR3==TRUE] <- 1
  
  data$Stop_codons[data$Stop_codons==FALSE] <- 0
  data$Stop_codons[data$Stop_codons==TRUE] <- 1
  
  data$Start_codons[data$Start_codons==FALSE] <- 0
  data$Start_codons[data$Start_codons==TRUE] <- 1
  
  data$TSS[data$TSS==FALSE] <- 0
  data$TSS[data$TSS==TRUE] <- 1
  
  data$TSS_A[data$TSS_A==FALSE] <- 0
  data$TSS_A[data$TSS_A==TRUE] <- 1
  
  data$exon_stop[data$exon_stop==FALSE] <- 0
  data$exon_stop[data$exon_stop==TRUE] <- 1
  
  data$TSS[data$TSS==FALSE] <- 0
  data$TSS[data$TSS==TRUE] <- 1
  
  data$alternative_exon[data$alternative_exon==FALSE] <- 0
  data$alternative_exon[data$alternative_exon==TRUE] <- 1
  
  data$constitutive_exon[data$constitutive_exon==FALSE] <- 0
  data$constitutive_exon[data$constitutive_exon==TRUE] <- 1
  
  data$internal_exon[data$internal_exon==FALSE] <- 0
  data$internal_exon[data$internal_exon==TRUE] <- 1
  
  data$long_exon[data$long_exon==FALSE] <- 0
  data$long_exon[data$long_exon==TRUE] <- 1
  
  data$last_exon[data$last_exon==FALSE] <- 0
  data$last_exon[data$last_exon==TRUE] <- 1
  
  data$last_exon_400bp[data$last_exon_400bp==FALSE] <- 0
  data$last_exon_400bp[data$last_exon_400bp==TRUE] <- 1
  
  data$last_exon_sc400[data$last_exon_sc400==FALSE] <- 0
  data$last_exon_sc400[data$last_exon_sc400==TRUE] <- 1
  
  data$intron[data$intron==FALSE] <- 0
  data$intron[data$intron==TRUE] <- 1
  
  data$sncRNA[data$long_exon==FALSE] <- 0
  data$sncRNA[data$long_exon==TRUE] <- 1
  
  data$lncRNA[data$long_exon==FALSE] <- 0
  data$lncRNA[data$long_exon==TRUE] <- 1
  
  return(data)
}

pos_encoding <- as.data.frame(transform(pos_encoding))
neg_encoding <- as.data.frame(transform(neg_encoding))

pos_encoding <- pos_encoding[,-c(17:18, 30:41, 52)]
neg_encoding <- neg_encoding[,-c(17:18, 30:41, 52)]

pos_encoding <- pos_encoding[,-c(38)]
neg_encoding <- neg_encoding[,-c(38)]

#pos_encoding <- pos_encoding[,-c(19,20,30,31)]
#neg_encoding <- neg_encoding[,-c(19,20,30,31)]

saveRDS(pos_encoding,"/Users/jiaming/Desktop/f5c/pos_domain_encoding.rds")
saveRDS(neg_encoding,"/Users/jiaming/Desktop/f5c/neg_domain_encoding.rds")

write.csv(pos_encoding,"/Users/jiaming/Desktop/f5c/pos_domain_encoding.csv")
write.csv(neg_encoding,"/Users/jiaming/Desktop/f5c/neg_domain_encoding.csv")


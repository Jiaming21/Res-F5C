library(seqinr)
library("m6ALogisticModel")
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(AnnotationDbi)
library(stringr)

#得到pos的序列
pos_C <- readRDS('/Users/jiaming/Desktop/f5c/f5c_positive_all.rds')
pos_C_seq <- as.data.frame(pos_C$refseq)

#得到neg的序列
txdb <- loadDb("/Users/jiaming/Desktop/f5c/S288C.sqlite")
pos_C <- pos_C + 20 
same_transcript <- subsetByOverlaps(transcripts(txdb), pos_C) #得到位点信息
neg_C <- m6ALogisticModel::sample_sequence("C", same_transcript, Scerevisiae, Fixed = T)
neg_C <- subsetByOverlaps(neg_C,exons(txdb))
olp <- findOverlaps(neg_C,pos_C,ignore.strand= T)
neg_C <- neg_C[-queryHits(olp)]
neg_C <- neg_C[sample(seq_along(neg_C),10*length(pos_C))]
rm(olp,same_transcript)
neg_C <- neg_C + 20
neg_C_seq <- as.data.frame(DNAStringSet(Views(Scerevisiae,neg_C)))

set.seed(1)
index <- c(sample(1:18920,1892))
neg_C_seq <- as.data.frame(neg_C_seq[index,])

#编码
#1. NCP & ELLP
encoding_NCP_ELLP <- function(seqs){ #传入dataframe
  num_rows <- nrow(seqs)
  num_cols <- 4*41 #NCP(3)+ELLP(1)
  na_matrix <- matrix(NA, nrow = num_rows, ncol = num_cols)
  encoding_result <- as.data.frame(na_matrix)
  for (i in 1:num_rows) { #锁定序列
    for(j in 1:41){ #锁定碱基
      if(substr(seqs[i,], start = j, stop = j) == "A"){
        encoding_result[i,((j-1)*4+1):((j-1)*4+4)] <- c(1, 1, 1, 0.1260)
      }
      if(substr(seqs[i,], start = j, stop = j) == "C"){
        encoding_result[i,((j-1)*4+1):((j-1)*4+4)] <- c(0, 1, 0, 0.1340)
      }
      if(substr(seqs[i,], start = j, stop = j) == "G"){
        encoding_result[i,((j-1)*4+1):((j-1)*4+4)] <- c(1, 0, 0, 0.0806)
      }
      if(substr(seqs[i,], start = j, stop = j) == "T"){
        encoding_result[i,((j-1)*4+1):((j-1)*4+4)] <- c(0, 0, 1, 0.1335)
      }
    }
  }
  return(encoding_result)
}

pos_encoding_NCP_ELLP <- encoding_NCP_ELLP(pos_C_seq)
neg_encoding_NCP_ELLP <- encoding_NCP_ELLP(neg_C_seq)

write.csv(pos_encoding_NCP_ELLP,'/Users/jiaming/Desktop/f5c/pos_encoding_NCP_ELLP.csv')
write.csv(neg_encoding_NCP_ELLP,'/Users/jiaming/Desktop/f5c/neg_encoding_NCP_ELLP.csv')

#2. OH & ELLP
encoding_OH_ELLP <- function(seqs){ #传入dataframe
  num_rows <- nrow(seqs)
  num_cols <- 5*41 #NCP(3)+ELLP(1)
  na_matrix <- matrix(NA, nrow = num_rows, ncol = num_cols)
  encoding_result <- as.data.frame(na_matrix)
  for (i in 1:num_rows) { #锁定序列
    for(j in 1:41){ #锁定碱基
      if(substr(seqs[i,], start = j, stop = j) == "A"){
        encoding_result[i,((j-1)*5+1):((j-1)*5+5)] <- c(1, 0, 0, 0, 0.1260)
      }
      if(substr(seqs[i,], start = j, stop = j) == "C"){
        encoding_result[i,((j-1)*5+1):((j-1)*5+5)] <- c(0, 1, 0, 0, 0.1340)
      }
      if(substr(seqs[i,], start = j, stop = j) == "G"){
        encoding_result[i,((j-1)*5+1):((j-1)*5+5)] <- c(0, 0, 1, 0, 0.0806)
      }
      if(substr(seqs[i,], start = j, stop = j) == "T"){
        encoding_result[i,((j-1)*5+1):((j-1)*5+5)] <- c(0, 0, 0, 1, 0.1335)
      }
    }
  }
  return(encoding_result)
}

pos_encoding_OH_ELLP <- encoding_OH_ELLP(pos_C_seq)
neg_encoding_OH_ELLP <- encoding_OH_ELLP(neg_C_seq)

write.csv(pos_encoding_OH_ELLP,'/Users/jiaming/Desktop/f5c/pos_encoding_OH_ELLP.csv')
write.csv(neg_encoding_OH_ELLP,'/Users/jiaming/Desktop/f5c/neg_encoding_OH_ELLP.csv')

#3. NCP & ND
encoding_NCP_ND <- function(seqs){ #传入dataframe
  num_rows <- nrow(seqs)
  num_cols <- 4*41 #NCP(3)+ELLP(1)
  na_matrix <- matrix(NA, nrow = num_rows, ncol = num_cols)
  encoding_result <- as.data.frame(na_matrix)
  for (i in 1:num_rows) { #锁定序列
    for(j in 1:41){ #锁定碱基
      ND <- str_count(substr(seqs[i,], start = 1, stop = j), substr(seqs[i,], start = j, stop = j)) / j
      if(substr(seqs[i,], start = j, stop = j) == "A"){
        encoding_result[i,((j-1)*4+1):((j-1)*4+4)] <- c(1, 1, 1, ND)
      }
      if(substr(seqs[i,], start = j, stop = j) == "C"){
        encoding_result[i,((j-1)*4+1):((j-1)*4+4)] <- c(0, 1, 0, ND)
      }
      if(substr(seqs[i,], start = j, stop = j) == "G"){
        encoding_result[i,((j-1)*4+1):((j-1)*4+4)] <- c(1, 0, 0, ND)
      }
      if(substr(seqs[i,], start = j, stop = j) == "T"){
        encoding_result[i,((j-1)*4+1):((j-1)*4+4)] <- c(0, 0, 1, ND)
      }
    }
  }
  return(encoding_result)
}

pos_encoding_NCP_ND <- encoding_NCP_ND(pos_C_seq)
neg_encoding_NCP_ND <- encoding_NCP_ND(neg_C_seq)

write.csv(pos_encoding_NCP_ND,'/Users/jiaming/Desktop/f5c/pos_encoding_NCP_ND.csv')
write.csv(neg_encoding_NCP_ND,'/Users/jiaming/Desktop/f5c/neg_encoding_NCP_ND.csv')

#4. OH & ND
encoding_OH_ND <- function(seqs){ #传入dataframe
  num_rows <- nrow(seqs)
  num_cols <- 5*41 #NCP(3)+ELLP(1)
  na_matrix <- matrix(NA, nrow = num_rows, ncol = num_cols)
  encoding_result <- as.data.frame(na_matrix)
  for (i in 1:num_rows) { #锁定序列
    for(j in 1:41){ #锁定碱基
      ND <- str_count(substr(seqs[i,], start = 1, stop = j), substr(seqs[i,], start = j, stop = j)) / j
      if(substr(seqs[i,], start = j, stop = j) == "A"){
        encoding_result[i,((j-1)*5+1):((j-1)*5+5)] <- c(1, 0, 0, 0, ND)
      }
      if(substr(seqs[i,], start = j, stop = j) == "C"){
        encoding_result[i,((j-1)*5+1):((j-1)*5+5)] <- c(0, 1, 0, 0, ND)
      }
      if(substr(seqs[i,], start = j, stop = j) == "G"){
        encoding_result[i,((j-1)*5+1):((j-1)*5+5)] <- c(0, 0, 1, 0, ND)
      }
      if(substr(seqs[i,], start = j, stop = j) == "T"){
        encoding_result[i,((j-1)*5+1):((j-1)*5+5)] <- c(0, 0, 0, 1, ND)
      }
    }
  }
  return(encoding_result)
}

pos_encoding_OH_ND <- encoding_OH_ND(pos_C_seq)
neg_encoding_OH_ND <- encoding_OH_ND(neg_C_seq)

write.csv(pos_encoding_OH_ND,'/Users/jiaming/Desktop/f5c/pos_encoding_OH_ND.csv')
write.csv(neg_encoding_OH_ND,'/Users/jiaming/Desktop/f5c/neg_encoding_OH_ND.csv')


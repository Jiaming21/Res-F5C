library(jsonlite)

input_json <- commandArgs(trailingOnly = T)
jobID <- input_json[1]

a <- as.data.frame(fromJSON(paste0('/var/www/html/Deep-f5C/job/',jobID,'/',jobID,'_para.json')))
file <- as.character(a$file)
jobID <- as.character(a$jobID)

target_dir1 <- file 
target_dir2 <- '/home/jiaming/webserver/webserver_home'
write.table(target_dir1, paste0(target_dir2,'/target_dir1.txt'))
target_dir1 <- read.table(paste0(target_dir2,'/target_dir1.txt'))[1,1]

################################################################################
system("R -f /home/jiaming/webserver/webserver_step1_sequence_encoding_41bp_4types.R")

system("R -f /home/jiaming/webserver/webserver_step2_domain_encoding.R")

system("python3 /home/jiaming/webserver/webserver_step3.py")
################################################################################
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(AnnotationDbi)
library(seqinr)
library("m6ALogisticModel")
library(stringr)

data <- read.table(target_dir1, header = TRUE)
print(data)
data <- GRanges(seqnames = data$seqnames,
                ranges = IRanges(start = data$position,
                                 end = data$position,
                                 width = 1),
                strand = data$strand)
webserver_input_example <- data
df <- as.data.frame(webserver_input_example)

preds <- read.table('/home/jiaming/webserver/webserver_home/preds.csv', header = TRUE)
print(preds)

is_f5C <- data.frame(is_f5C = rep(FALSE, length(data)))
is_f5C$is_f5C[which(preds$preds>=0.5)] = "True"
is_f5C$is_f5C[which(preds$preds<0.5)] = "False"

LR <- (preds$preds)/(1-preds$preds)
indx99 <- which(LR > 99)
if(length(indx99) != 0){
  LR[indx99] <- 99
}

CL <- LR
threshold_1 <- 99
threshold_2 <- 80
if (CL >= threshold_1){
  CL <- 'high'
}else if (CL >= threshold_2){
  CL <- 'medium'
}else{
  CL <- 'low'
}

results <- data.frame(
  seqnames=as.character(df$seqnames),
  start = as.character(df$start),
  end = as.character(df$end),
  width = as.character(df$width),
  strand = as.character(df$strand),
  preds = as.character(preds$preds),
  is_f5C = as.character(is_f5C$is_f5C),
  'Likelihood Ratio' = as.character(LR$x),
  'Confidence Level' = as.character(CL$x)
  )

################################################################################
write_json(results, paste0('/var/www/html/Deep-f5C/job/',jobID,'/result/results.json'))
write.csv(results, paste0('/var/www/html/Deep-f5C/job/',jobID,'/result/results.csv'), row.names = FALSE)

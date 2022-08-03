# Copyright © 2022 CDI Laboratories Inc. | All Rights Reserved |

version <- "C0.3"

library(preprocessCore)
#
# This is the code to install preprocessCore
#
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("preprocessCore")

#
# Input and Output Folder arguments: must be two full pathnames or nothing
#

input_fpath <- "input"
output_fpath <- "output"

args <- commandArgs(trailingOnly=TRUE)

if(length(args) == 2){
  input_fpath <- args[1]
  output_fpath <- args[2]
} else {
  if(length(args) != 0){
    quit(save = "no", status = 0)
  }
}

CONTROL_PATTERN <- "Control"

SUPPRESSED_LABELS <- c("TRIM21", 'control', 'blank', 'empty', 'nd', 'n.d', 'n.d.', 'none', 'bsa', 'buffer', 'subclass', "igha1", "igha2", "ighd", "ighg1", "ighg3", "ighg4", "ighv5-78", "ighv7-81", "igip", "igk", "igkc", "igkv1or2-108", "igl@", "iglc2", "igll1", "igll5", "iglon5", "iglv6-57")

channel_path = paste(input_fpath, "/channels.txt", sep="")
if(!file.exists(channel_path)){
  print(channel_path)
  quit(save = "no", status = 1)
}

channels <- read.delim(channel_path,sep="\t", stringsAsFactors=FALSE)

drop_gpr_suffix <- function(aName){
  substr(aName,1,nchar(aName)-4)
}

info_column_count <- 5  # Block, Column, Row, Name, ID
load_gpr <- function(filename, keep_info_columns=FALSE, channel="red"){
  x = readLines(filename, 2)
  if(strsplit(trimws(x[1]), '\t')[[1]][1] != "ATF"){
    stop("Expected an Axon Text File!!!")
  }
  if(as.double(strsplit(trimws(x[1]), '\t')[[1]][2]) != 1.0){
    stop("Expected an Axon Text File ver=1.0!!!")
  }
  to_skip <- strtoi(strsplit(x[2],"\t")[[1]][1]) + 2
  d <- read.delim(filename, skip=to_skip, stringsAsFactors=FALSE)
  if(channel=="red"){
    channel_data <- data.frame(d$F635.Median)
  }
  if(channel=="green"){
    channel_data <- data.frame(d$F532.Median)
  }
  colnames(channel_data) <- c(drop_gpr_suffix(filename))
  if(keep_info_columns){
    # The next line is why info_column_count == 5...
    info_columns <- data.frame(Block=d$Block,Row=d$Row,Column=d$Column,Name=d$Name,ID=d$ID,stringsAsFactors=FALSE)
    info_columns$Name <- trimws(info_columns$Name)  # This is necessary because '\n' was found in the Names...
    info_columns$ID <- trimws(info_columns$ID)
    return(cbind(info_columns, channel_data))
  } else {
    return(channel_data)
  }
}


do_A_vs_B_comparisons <- function(dcols, info_cols, comparison_fname, output_name){
  comparisons <- read.delim(comparison_fname, sep="\t", stringsAsFactors=FALSE)
  comparisons$Sample <- drop_gpr_suffix(comparisons$Sample)
  
  data <- as.matrix(dcols)
  data[data == 0] <- 1
  data <- log2(data)
  stopifnot(colnames(data) == comparisons$Sample)
  
  complete <- info_cols
  for(col in 2:ncol(comparisons)){
    labels <- levels(factor(comparisons[,col], exclude=""))
    if(tolower(labels[1]) == "control"){
      control <- labels[1]
      treatment <- labels[2]
    } else {
      control <- labels[1]
      treatment <- labels[2]
    }
    title <- colnames(comparisons)[col]
    AvsB_pval <- c()
    AvsB_ratio <- c()
    for(i in 1:nrow(data)){
      A <- comparisons[,col]==treatment
      B <- comparisons[,col]==control
      AvsB_ratio <- c(AvsB_ratio,mean(data[i,A])-mean(data[i,B]))
      if(sd(data[i,A])>0 && sd(data[i,B])>0){
        AvsB_pval <- c(AvsB_pval,t.test(data[i,A],data[i,B])$p.value)
      } else {
        AvsB_pval <- c(AvsB_pval,1)
      }
    }
    
    output <- cbind( 
      AvsB_ratio, 
      AvsB_pval, 
      p.adjust(AvsB_pval, method="bonferroni"), 
      p.adjust(AvsB_pval, method="BH"),
      p.adjust(AvsB_pval, method="BY"))
    colnames(output)[1] <- paste(title, ": Ratio (Log2)", sep="")
    colnames(output)[2] <- paste(title, ": P-value (T-test)", sep="")
    colnames(output)[3] <- paste(title, ": Bonferroni", sep="")
    colnames(output)[4] <- paste(title, ": FDR (BH)", sep="")
    colnames(output)[5] <- paste(title, ": FDR (BY)", sep="")
    complete <- cbind(complete, output)
  }
  write.table(complete, file=output_name, sep="\t", row.names = FALSE)
}

run_pipeline <- function(input_fpath, output_fpath, channel="red", secondary=""){

  startup_path <- getwd()
  setwd(input_fpath)

  files <- list.files(".", pattern="\\.gpr$")
  for(i in 1:length(files)){
    if(i == 1){
      out <- load_gpr(files[i], keep_info_columns=TRUE, channel=channel)
    } else {
      out <- cbind(out, load_gpr(files[i], channel=channel))
    }
  }
  
  setwd(startup_path)

  spot_num <- nrow(out)
  max_col <- max(out$Column)
  full_rows <- spot_num / max_col
  m <- matrix(1:spot_num, max_col, full_rows)
  mevens <- m[,seq(2, full_rows, 2)]
  modds <- m[,seq(1, full_rows-1, 2)]
  
  controls <- out$ID[modds] == CONTROL_PATTERN
  
  average <- (out[mevens, (-1):(-info_column_count)] + out[modds, (-1):(-info_column_count)])/2

  suppressed <- (tolower(out$Name)[mevens] %in% SUPPRESSED_LABELS) | controls
  dcols <- average
  info_cols <- out[mevens, 1:info_column_count]
  info_cols$Caveat <- suppressed

  dcolsn <- normalize.quantiles(as.matrix(dcols))
  colnames(dcolsn) <- colnames(dcols)

  if(channel == "red"){
    if(secondary != ""){
      channel_name = paste(secondary, "_635nm_Red", sep="")
    } else{
      channel_name = "635nm_Red"
    }
  } else {
    if(secondary != ""){
      channel_name = paste(secondary, "_532nm_Green", sep="")
    } else{
      channel_name = "532nm_Green"
    }
  }
  
  A_vs_B_comparison_fname = paste(input_fpath, "/comparisons.tsv", sep="")
  if(file.exists(A_vs_B_comparison_fname)){
    do_A_vs_B_comparisons(dcols, info_cols, A_vs_B_comparison_fname, paste(output_fpath, "/Case_vs_Control_Raw_PairAvg_", channel_name, ".tsv",sep=""))
    do_A_vs_B_comparisons(dcolsn, info_cols, A_vs_B_comparison_fname, paste(output_fpath, "/Case_vs_Control_QuantileNorm_PairAvg_", channel_name, ".tsv",sep=""))
  }
  
  combined_gprsn <-cbind(info_cols, dcolsn)
  write.table(combined_gprsn, file=paste(output_fpath, "/Combined_QuantileNorm_PairAvg_", channel_name, ".tsv", sep=""), sep="\t", row.names=FALSE)
  
  combined_gprs <-cbind(info_cols, dcols)
  write.table(combined_gprs, file=paste(output_fpath, "/Combined_Raw_PairAvg_", channel_name, ".tsv", sep=""), sep="\t", row.names=FALSE)
  
  combined_ctrls <-cbind(out[mevens, 1:info_column_count][controls, ], average[controls,])
  write.table(combined_ctrls, file=paste(output_fpath, "/Controls_Raw_PairAvg_", channel_name, ".tsv", sep=""), sep="\t", row.names=FALSE)
}

for(i in 1:nrow(channels)){
  if(channels$Label[i] == "red"){
    run_pipeline(input_fpath, output_fpath, channel="red", secondary=channels$Secondary[i])
  }
  if(channels$Label[i] == "green"){
    run_pipeline(input_fpath, output_fpath, channel="green", secondary=channels$Secondary[i])
  }
}

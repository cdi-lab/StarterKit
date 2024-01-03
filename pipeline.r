# Copyright 2022 CDI Laboratories Inc. | All Rights Reserved |

version <- "C0.1"

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

run_pipeline <- function(input_fpath, output_fpath, channel="red"){

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
    channel_name = "635nm_Red"
  } else {
    channel_name = "532nm_Green"
  }
  
  combined_gprsn <-cbind(info_cols, dcolsn)
  write.table(combined_gprsn, file=paste(output_fpath, "/Combined_QuantileNorm_PairAvg_", channel_name, ".tsv", sep=""), sep="\t", row.names=FALSE)
  
  combined_gprs <-cbind(info_cols, dcols)
  write.table(combined_gprs, file=paste(output_fpath, "/Combined_Raw_PairAvg_", channel_name, ".tsv", sep=""), sep="\t", row.names=FALSE)
  
  combined_ctrls <-cbind(out[mevens, 1:info_column_count][controls, ], average[controls,])
  write.table(combined_ctrls, file=paste(output_fpath, "/Controls_Raw_PairAvg_", channel_name, ".tsv", sep=""), sep="\t", row.names=FALSE)
}

run_pipeline(input_fpath, output_fpath, channel="red")
run_pipeline(input_fpath, output_fpath, channel="green")

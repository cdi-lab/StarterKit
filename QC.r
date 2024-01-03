# Copyright 2023 CDI Laboratories Inc. | All Rights Reserved |

version <- "C0.1"

library(tidyverse)
library(limma)

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


#
# QC Parameters:
#

BAD_RATIO_PROPORTION = 0.02
MIN_SPOT_PAIR_CORRELATION = 0.95

CONTROL_PATTERN <- "Control"

targets_path <- paste(input_fpath, "/targets.txt", sep="")
if(file.exists(targets_path)){
  targets <- readTargets(targets_path)
} else {
  files <- list.files(input_fpath,pattern="\\.gpr$")
  targets <- data.frame(FileName=files)
  rownames(targets) <- substr(files, 1, nchar(files)-4)
}

startup_path <- getwd()
setwd(input_fpath)

RG <- read.maimages(
  targets, 
  source="genepix.median",
  annotation=c("Block", "Row","Column", "ID", "Name", "X")
)

setwd(startup_path)

if(is.na(RG$printer$ngrid.c)){
  nblocks <- max(RG$genes$Block, na.rm=TRUE)
  blocksize <- RG$printer$nspot.r * RG$printer$nspot.c
  i <- (1:(nblocks - 1)) * blocksize
  RG$printer$ngrid.c <- which(RG$genes$X[i] > RG$genes$X[i + 1])[1]
  RG$printer$ngrid.r <- nblocks / RG$printer$ngrid.c
}

RG$genes$Control <- 0
RG$genes$Control[RG$genes$ID == CONTROL_PATTERN] <- 1

spot_names <- paste("B",RG$genes$Block,"R",RG$genes$Row,"C",RG$genes$Column,":",RG$genes$ID,sep="")

RG$R[RG$R <= 0] <- 1
RG$Rb[RG$Rb <= 0] <- 1
RG$G[RG$G <= 0] <- 1
RG$Gb[RG$Gb <= 0] <- 1

r <- log2(RG$R)
g <- log2(RG$G)

rownames(r) <- spot_names
rownames(g) <- spot_names

spot_num <- dim(r)[1]
full_rows <- spot_num / RG$printer$nspot.c
m <- matrix(1:spot_num, RG$printer$nspot.c, full_rows)
mevens <- m[,seq(2, full_rows, 2)]
modds <- m[,seq(1, full_rows-1, 2)]

critical_value = BAD_RATIO_PROPORTION * length(mevens)


r_cors <- c()
r_ccors <- c()
r_outs <- c()
r_couts <- c()

for(i in 1:dim(r)[2]){
  r_cors <- c(r_cors, cor(r[mevens,i],r[modds,i]))
  r_ccors <- c(r_ccors, cor(r[mevens[RG$genes$Control[mevens]==1],i],r[modds[RG$genes$Control[mevens]==1],i]))
  r_outs <- c(r_outs, sum(abs(r[mevens,i]-r[modds,i])>1))
  r_couts <- c(r_couts, sum(abs(r[mevens[RG$genes$Control[mevens]==1],i]-r[modds[RG$genes$Control[mevens]==1],i])>1))
}

targets$r_cors <- r_cors
targets$r_ccors <- r_ccors
targets$r_outs <- r_outs
targets$r_couts <- r_couts

png(filename = paste(output_fpath, "/red_cors.png", sep=""), width = 900, height = 900, res=150)
plot(r_ccors,type="l",col="red",  ylim=c(0,1.0), main="replicate correlations",ylab="correlation coefficient", xlab="Slide#", sub="black: total, red: controls, purple: QC-threshold")
lines(r_cors,col="black")
abline(h=MIN_SPOT_PAIR_CORRELATION, col="purple", lty=3, lwd=3)
dev.off()

png(filename = paste(output_fpath, "/red_outs.png", sep=""), width = 900, height = 900, res=150)
plot(r_outs,type="l",col="black", ylim=c(0,max(max(r_outs),2*critical_value)), main="replicate-spots having ratios greater than 2-fold",ylab="abnormal spot-pairs", xlab="Slide#", sub="black: total, red: controls, purple: QC-threshold")
lines(r_couts,col="red")
abline(h=critical_value, col="purple", lty=3, lwd=3)
dev.off()


g_cors <- c()
g_ccors <- c()
g_outs <- c()
g_couts <- c()

for(i in 1:dim(g)[2]){
  g_cors <- c(g_cors, cor(g[mevens,i],g[modds,i]))
  g_ccors <- c(g_ccors, cor(g[mevens[RG$genes$Control[mevens]==1],i],g[modds[RG$genes$Control[mevens]==1],i]))
  g_outs <- c(g_outs, sum(abs(g[mevens,i]-g[modds,i])>1))
  g_couts <- c(g_couts, sum(abs(g[mevens[RG$genes$Control[mevens]==1],i]-g[modds[RG$genes$Control[mevens]==1],i])>1))
}

png(filename = paste(output_fpath, "/green_cors.png", sep=""), width = 900, height = 900, res=150)
plot(g_ccors,type="l",col="red",  ylim=c(0,1.0), main="replicate correlations",ylab="correlation coefficient", xlab="Slide#", sub="black: total, red: controls, purple: QC-threshold")
lines(g_cors,col="black")
abline(h=MIN_SPOT_PAIR_CORRELATION, col="purple", lty=3, lwd=3)
dev.off()

png(filename = paste(output_fpath, "/green_outs.png", sep=""), width = 900, height = 900, res=150)
plot(g_outs,type="l",col="black", ylim=c(0,max(max(g_outs),2*critical_value)), main="replicate-spots having ratios greater than 2-fold",ylab="abnormal spot-pairs", xlab="Slide#", sub="black: total, red: controls, purple: QC-threshold")
lines(g_couts,col="red")
abline(h=critical_value, col="purple", lty=3, lwd=3)
dev.off()

targets$g_cors <- g_cors
targets$g_ccors <- g_ccors
targets$g_outs <- g_outs
targets$g_couts <- g_couts


write.table(targets, file=paste(output_fpath, "/reports.tsv", sep=""), sep="\t", row.names=FALSE)

#
# Control Tables
#

red_avg <- cbind(RG$genes[modds,],(r[mevens,] + r[modds,])/2)
red_controls <- red_avg[red_avg$Control==1, -which(names(red_avg) %in% c("Control", "X"))]
write.table(red_controls, file=paste(output_fpath, "/red_controls.tsv", sep=""),sep="\t",row.names = FALSE)

green_avg <- cbind(RG$genes[modds,],(g[mevens,] + g[modds,])/2)
green_controls <- green_avg[green_avg$Control==1, -which(names(green_avg) %in% c("Control", "X"))]
write.table(green_controls, file=paste(output_fpath, "/green_controls.tsv", sep=""),sep="\t",row.names = FALSE)

#
# IgG Control Plot
#

pdf(paste(output_fpath, "/Human_IgG Controls_Red.pdf", sep=""))
plot_obj <- red_controls %>% 
  filter(startsWith(Name,"human IgG ")) %>% 
  mutate(Row=paste("Row", floor((Block-1) / RG$printer$ngrid.c) + 1)) %>%
  mutate(Column=paste("Column", Block - RG$printer$ngrid.c*floor((Block-1) / RG$printer$ngrid.c))) %>%
  mutate(Conc=as.numeric(str_extract(Name, "\\-*\\d+\\.*\\d*"))) %>%
  pivot_longer(6:ncol(red_controls),"Sample") %>%
  ggplot(aes(Conc,value, hwy, colour = Sample, group=Sample)) + 
  geom_point() + geom_line() + facet_grid(Row ~ Column)+
  scale_x_continuous(trans = 'log10') +
  ggtitle("Human IgG Controls by Block (Red Channel)") +
  xlab("Concentration [ng/ul]") + ylab("Log Intensity")
  print(plot_obj + theme(legend.position = "none"))
dev.off()



#
# Anti-IgG Control Plot
#
pdf(paste(output_fpath, "/Anti_Human_IgG Controls_Red.pdf", sep=""))
plot_obj <- red_controls %>% 
  filter(startsWith(Name,"Anti-human IgG ")) %>% 
  mutate(Row=paste("Row", floor((Block-1) / RG$printer$ngrid.c) + 1)) %>%
  mutate(Column=paste("Column", Block - RG$printer$ngrid.c*floor((Block-1) / RG$printer$ngrid.c))) %>%
  mutate(Conc=as.numeric(str_extract(Name, "\\-*\\d+\\.*\\d*"))) %>%
  pivot_longer(6:ncol(red_controls),"Sample") %>%
  ggplot(aes(Conc,value, hwy, colour = Sample, group=Sample)) + 
  geom_point() + geom_line() + facet_grid(Row ~ Column)+
  scale_x_continuous(trans = 'log10') +
  ggtitle("Anti-Human IgG Controls by Block (Red Channel)") +
  xlab("Concentration [ng/ul]") + ylab("Log Intensity")
  print(plot_obj + theme(legend.position = "none"))
dev.off()


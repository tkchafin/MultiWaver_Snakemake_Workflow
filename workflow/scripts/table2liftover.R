############################################################
#   R script for converting ancestry breakpoints to centimorgans
#   & creating input files for 
############################################################

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

#setwd("~/projects/red_wolf/liftover_coords")

args = commandArgs(trailingOnly=TRUE)

#if (length(args) < 3){
#  stop("Wrong number of inputs")
#}else{
#  oname <- args[1]
#  liftover <- args[2]
#  fpath <- args[3]
#}
fpath<-"/Users/tyler/projects/red_wolf/liftover_coords"
oname<-"data/rep_1.seg"
liftover<-"data/coords.tsv.liftover"
#oname<-paste0("rep_",rep,".seg")

print(paste0("Output file will be called ",oname))

#read in linkage map data
print("Reading liftover file")
lo<-read.table(liftover, sep="\t", header=T)

lo$chr<-replace(as.character(lo$chr), lo$chr == "X", "39")

#initialize data frame
wolf_dat <- data.frame(matrix(ncol = 25, nrow = 0))

#grab names of Rdata files (from trestles)
print(paste0("Loading Rdata files from ", fpath, "... (this might take a while"))
files<-list.files(path=fpath, pattern="*.Rdata", full.names=T)
print(paste0("Found files: ", files))
#info <- file.info(files)
#info$size_mb <- info$size/(1024 * 1024)
#print(subset(info, select=c("size_mb")))

#read in the RDS and append to existing df (...this might take a while)
datalist = list()
i<-1
for (f in files){
  datalist[[i]] <- readRDS(f)
  i <- i+1
}
wolf_dat <- rbind.fill(datalist)

wolf_dat[,c(4:25)]<-sapply(wolf_dat[,c(4:25)], as.numeric)

colnames(wolf_dat) <- c("chrom", "id", "sister","mrca_time", 
                        "cg", "rg", "cr", "sister_branch_rate", 
                        "D", "Dsig", "Nabba", "Nbaba", "RND", "Dxy", 
                        "start", "end", "aln_len", "nucdiv", "R-logL", 
                        "R-deltaL", "R-bp.RELL", 
                        "G-logL", "G-deltaL", "G-bp.RELL", 
                        "int.het")

nrow(wolf_dat)
ncol(wolf_dat)

chroms<-unique(wolf_dat$chrom)

subs<-c("CMW", "YSW")

getHetRand <- function(sis, het, threshold){
  subs<-c("CMW", "YSW")
  if (!is.na(het) && het > threshold || sis=="NONE"){
    return(sample(subs, 1))
  }else{
    return(sis)
  }
}
getHetExclude <- function(sis, het, threshold){
  subs<-c("CMW", "YSW")
  if (!is.na(het) && het > threshold || sis=="NONE"){
    return("NONE")
  }else{
    return(sis)
  }
}

getChromId <- function(chrom){
  sub<-substr(chrom, 7, 8)
  sub <- gsub("(^|[^0-9])0+", "\\1", sub, perl = TRUE)
  return(sub)
}


for (rep in 1:100){

  oname<-paste0("rep_",rep,".seg")
  
print("Computing ancestry assignments for each segment using a int.het cutoff of 0.5")
wolf_dat["ASSIGN"] <- mapply(getHetRand, wolf_dat$sister, wolf_dat$int.het, 0.5)
wolf_dat["ASSIGN2"] <- mapply(getHetExclude, wolf_dat$sister, wolf_dat$int.het, 0.5)
wolf_dat["CHROM_ID"] <- mapply(getChromId, wolf_dat$chrom)

chroms<-unique(wolf_dat$CHROM_ID)

library(pracma)
library(plyr)
library(dplyr)

otable<-data.frame("start"=numeric(), end=numeric(), ancestry=numeric())
print("Interpolating cM positions")
for (c in chroms){
  print(c)
  
  first_row<-data.frame(c,1, 0.0000000, 1)
  names(first_row)<-c("chr","bp", "cM", "liftover.bp")
  
  sub_lo <- rbind(first_row, lo[lo$chr==c,])
  sub_lo <- sub_lo[order(sub_lo$liftover.bp),]
  sub_df <- wolf_dat[wolf_dat$CHROM_ID==c,]
  sub_df <- sub_df[order(sub_df$start),]
  #plot(sub_lo$liftover.bp, sub_lo$cM)
  #xs<-seq(1, max(sub_lo$liftover.bp), by = 10000)
  #ys<- cubicspline(sub_lo$liftover.bp, sub_lo$cM, xs)
  ynew<-cubicspline(sub_lo$liftover.bp, sub_lo$cM, sub_df$start)
  #plot(sub_df$start, ynew, cex=0.2, 
  #     xlab ="Physical position (bp)",
  #     ylab="Genetic position (cM)")
  #points(sub_lo$liftover.bp, sub_lo$cM, col="darkred", pch=1, cex=2)
  last<-cubicspline(sub_lo$liftover.bp, sub_lo$cM, tail(sub_df, 1)$end)
  sub_df["interp.cM"] <- ynew
  #print(sub_df)
  t<-mutate(sub_df, new = lead(interp.cM, default = "NA"))
  t[nrow(t),]$new <- last
  t[t$ASSIGN=="CMW",]$ASSIGN<-1
  t[t$ASSIGN=="YSW",]$ASSIGN<-2
  ovars<-c("interp.cM", "new", "ASSIGN")
  o<-t[ovars]
  otable<-rbind(otable, o)
}

print("Writing output")
write.table(otable, file=oname, append=F, quote=F, col.names=F,
            row.names=F, sep = "\t")

}


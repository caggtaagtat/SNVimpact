args <- commandArgs()
## Get hap-file name
filename <- args[6]

setwd("/pathTo/haploSNV/newVCF/haps/csv/")

##load respective hap-vcf as table
tt <- read.table(filename, quote="\"", comment.char="")

## Restructure table for downstream analysis
chr <- as.character(lapply(tt$V1 , function(x){  strsplit(x, "\\:")[[1]][[1]] }))
tt2 <- data.frame(CHROM=chr)

tt2$CHROM <- substr(tt2$CHROM, 4, nchar(tt2$CHROM))
tt2$POS <- tt$V3
tt2$ID <- "."
tt2$REF <- tt$V4
tt2$ALT <- tt$V5
tt2$QUAL <- 100
tt2$FILTER <- "*"
tt2$INFO <- ""
tt2$H1 <- tt$V6
tt2$H2 <- tt$V7

tt2 <- tt2[tt2$CHROM %in% c(as.character(1:24), "X", "Y", "M"),]
tt2$CHROM[tt2$CHROM == "M"] <- "MT"

tt2$nREF <- nchar(tt2$REF)
tt2$nALT <- nchar(tt2$ALT)

tt2$INFO[tt2$nREF == tt2$nALT] <- "TYPE=SUBSTITUTE"
tt2$INFO[tt2$nREF < tt2$nALT] <- "TYPE=INSERT"
tt2$INFO[tt2$nREF > tt2$nALT] <- "TYPE=DELETE"

tt2$nREF <- NULL
tt2$nALT <- NULL

## Separate by haplotype
tt2_h1 <- tt2[tt2$H1 == 1,]
tt2_h2 <- tt2[tt2$H2 == 1,]

tt2_h1$H1 <- NULL
tt2_h1$H2 <- NULL
tt2_h2$H1 <- NULL
tt2_h2$H2 <- NULL

filename_short <- gsub("hprc-v1.1-mc-chm13.GRCh38.", "",basename(filename))
filename_short <- gsub(".hap.vcf.hap","", filename_short)

## Save variations as table per haplotype
write.csv2(tt2_h1, file=paste0(filename_short,"_h1_SD_variations.csv"), row.names = F)
write.csv2(tt2_h2, file=paste0(filename_short,"_h2_SD_variations.csv"), row.names = F)
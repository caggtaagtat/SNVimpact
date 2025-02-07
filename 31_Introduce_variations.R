args <- commandArgs()
## Load csv of haplotyped variations
filename <- args[6]

## Load required R packages
library("VarCon")
library("stringr")
library("ModCon")
library("Biostrings")

## Load data
all_data <- read.csv2(filename)
names(all_data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

## Load reference data
load("genes_ens105")
load("GRCH38_fasta")
exonList <- read.csv("exonList_ens105.txt")
Transcript_cords_chr_ss_seqs <- read.csv("Transcript_cords_chr_ss_seqs_ens105.txt")

## Load custom functions
load("32_SDvariation.R")
load("33_GetReferenceSequence.R")

all_data$type[str_detect(all_data$INFO, "SUBSTITUTE")] <- "SNV"
all_data$type[str_detect(all_data$INFO, "INSERT")] <- "INS"
all_data$type[str_detect(all_data$INFO, "DELETE")] <- "DEL"
all_data$ID <- 1:nrow(all_data)

#first use fast Granges function for annotation
refGR <- makeGRangesFromDataFrame(genes)

all_data$start <- all_data$POS
all_data$end <- all_data$POS
all_data <- all_data[all_data$CHROM %in% c(1:22, "X", "Y"),]
inputGR <- makeGRangesFromDataFrame(all_data)

## find hits
hits <- findOverlaps(refGR, inputGR)
hits <- as.data.frame(hits)
hits$transcript <- genes$gene_name[match(hits$queryHits, genes$id)]

hits$ID    <- all_data$ID[match(hits$subjectHits, all_data$ID)]
hits$CHROM <- all_data$CHROM[match(hits$subjectHits, all_data$ID)]
hits$POS   <- all_data$POS[match(hits$subjectHits, all_data$ID)]
hits$REF   <- all_data$REF[match(hits$subjectHits, all_data$ID)]
hits$ALT   <- all_data$ALT[match(hits$subjectHits, all_data$ID)]
hits$type  <- all_data$type[match(hits$subjectHits, all_data$ID)]
hits$queryHits <- NULL
hits$subjectHits <- NULL
hits$txID <- hits$transcript
hits$transcript <- NULL
hits$strand <- Transcript_cords_chr_ss_seqs$Strand[match(hits$txID, Transcript_cords_chr_ss_seqs$Transcript.stable.ID )]
hits$geneid <- Transcript_cords_chr_ss_seqs$Gene.stable.ID[match(hits$txID, Transcript_cords_chr_ss_seqs$Transcript.stable.ID )]

all_data <- hits

## Summarize per SNP
all_data$ID <- 1:nrow(all_data)
snp <- data.frame(ID=unique(all_data$ID))
snp$gene <- ""
snp$chromosome <- all_data$CHROM[match(snp$ID, all_data$ID)]
snp$genome_position <- all_data$POS[match(snp$ID, all_data$ID)]
snp <- na.omit(snp)
snp$type <-  all_data$type[match(snp$ID, all_data$ID)]
snp$ref <-  all_data$REF[match(snp$ID, all_data$ID)]
snp$alt <-  all_data$ALT[match(snp$ID, all_data$ID)]
snp$ref_len <- nchar(snp$ref)
snp$cord1 <- snp$genome_position
snp$cord2 <- snp$genome_position + snp$ref_len - 1
snp$genome_position[snp$type == "DEL"] <- snp$genome_position[snp$type == "DEL"]+1
snp$cord1 <- snp$genome_position

snp$c.DNA <- paste0("g.",all_data$POS[match(snp$ID, all_data$ID)],
                    all_data$REF[match(snp$ID, all_data$ID)], ">",
                    all_data$ALT[match(snp$ID, all_data$ID)] )
snp$c.DNA[snp$type == "DEL"] <- paste0("g.",snp$cord1[snp$type == "DEL"],"_",snp$cord2[snp$type == "DEL"],
                                       "del", substring(snp$ref[snp$type == "DEL"], 2, nchar(snp$ref[snp$type == "DEL"])) )
snp$c.DNA[snp$type == "INS"] <- paste0("g.",snp$cord1[snp$type == "INS"],"_",snp$cord2[snp$type == "INS"],
                                       "ins", snp$alt[snp$type == "INS"] )

snp$strand <-  Transcript_cords_chr_ss_seqs$Strand[match(snp$ID, all_data$ID)]
snp$ENSEMBL_transcriptID <- all_data$txID[match(snp$ID, all_data$ID)]
snp$gene <- exonList$Gene.name[match(snp$ENSEMBL_transcriptID, exonList$Transcript.stable.ID)]
snp <- snp[snp$gene != "",]

snps <- snp[,c(-1, -2, -12, -13)]
snps <- unique(snps)

##Delete INDELS
snps <- snps[!(snps$type == "INS" & nchar(snps$ref) > 1),]


load("sd_ens105")

sd$ID <- 1:nrow(sd)



for(i in c(1:nrow(sd))){

print(i)
sd$sequences[i] <- SDvariation(sd[i,])

}


sd$ref_seq <- as.character(lapply(sd$sequences, function(x){
  
  strsplit(x, ",")[[1]][[1]]
  
}))

sd$alt_seq <- as.character(lapply(sd$sequences, function(x){
  
  strsplit(x, ",")[[1]][[2]]
  
}))



sd$n_SNV <- as.numeric(lapply(sd$sequences, function(x){
  
  strsplit(x, ",")[[1]][[length( strsplit(x, ",")[[1]])-2]]
  
}))

sd$n_DEL <- as.numeric(lapply(sd$sequences, function(x){
  
  strsplit(x, ",")[[1]][[length( strsplit(x, ",")[[1]])-1]]
  
}))

sd$n_INS <- as.numeric(lapply(sd$sequences, function(x){
  
  strsplit(x, ",")[[1]][[length( strsplit(x, ",")[[1]])]]
  
}))


sd$alt_seq[sd$alt_seq == 0] <- sd$ref_seq[sd$alt_seq == 0]

sd$donor_seq <- substr(sd$ref_seq, 61, 71)
sd$donor_seq_alt <- substr(sd$alt_seq, 61, 71)

sd$refup60 <- substr(sd$ref_seq, 1, 60)
sd$refdown60 <- substr(sd$ref_seq, 72, 131)

sd$altup60 <- substr(sd$alt_seq, 1, 60)
sd$altdown60 <- substr(sd$alt_seq, 72, 131)

library(stringr)

sd$Ns <- str_detect(sd$sequences, "N")
sd <- sd[!sd$Ns,]

library(ModCon)
cl <- makeCluster(10)
clusterExport(cl=cl, list("sd", "calculateHZEIint" ))

sd$refup60_HEZI <- as.numeric(parLapply(cl, sd$refup60, calculateHZEIint))
sd$refdown60_HEZI <- as.numeric(parLapply(cl, sd$refdown60, calculateHZEIint))
sd$altup60_HEZI <- as.numeric(parLapply(cl, sd$altup60, calculateHZEIint))
sd$altdown60_HEZI <- as.numeric(parLapply(cl, sd$altdown60, calculateHZEIint))

stopCluster(cl)


sd$ref_SSHW <- sd$refup60_HEZI - sd$refdown60_HEZI
sd$alt_SSHW <- sd$altup60_HEZI - sd$altdown60_HEZI

sd$ref_HBS <- hbg$hbs[match(sd$donor_seq, hbg$seq)]
sd$alt_HBS <- hbg$hbs[match(sd$donor_seq_alt, hbg$seq)]
sd$HBS_diff <- sd$alt_HBS-sd$ref_HBS 
sd$HBS_diff[sd$ref_HBS] <- 999

sd$SSHW_diff <- sd$alt_SSHW - sd$ref_SSHW

sd$n_variations <- sd$n_DEL + sd$n_INS + sd$n_SNV


write.csv2(sd, file=paste0("newVCFout/",basename(filename),"_SD_variations.csv"))








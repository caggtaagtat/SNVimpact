## Load exon data
ex <- read.csv("exons105.txt")
ex$c <- 1

## Find last exon per transcript
maxExonsTX <- tapply( ex$c, ex$Transcript.stable.ID, sum)
maxExonsTX <- data.frame(maxExonsTX)
maxExonsTX$tx <- row.names(maxExonsTX)
ex$maxExonRank <- maxExonsTX$maxExonsTX[match(ex$Transcript.stable.ID, maxExonsTX$tx)]

SDs <- ex[ex$maxExonRank != 1,]
SDs$SDcor <- SDs$Exon.region.end..bp.
SDs$SDcor[SDs$Strand == -1] <- SDs$Exon.region.start..bp.[SDs$Strand == -1]
SDs <- SDs[SDs$Exon.rank.in.transcript != SDs$maxExonRank,]

## Load transcript data information TSL
tsl105 <- read.csv("transcriptTSL105.txt")
library(stringr)
tsl105$anyTSL1 <- str_detect(tsl105$Transcript.support.level..TSL., "tsl1")
SDs$transcriptTSL <- tsl105$Transcript.support.level..TSL.[match(SDs$Transcript.stable.ID, tsl105$Transcript.stable.ID)]
SDs$TSL1 <- str_detect(SDs$transcriptTSL, "tsl1")
SDs$ensembl_canonical <- tsl105$Ensembl.Canonical[match(SDs$Transcript.stable.ID, tsl105$Transcript.stable.ID)]
SDs$sdID <- paste(SDs$Chromosome.scaffold.name, SDs$Strand, SDs$SDcor)


## Load haplotype SD variation data
library("VarCon")
library("stringr")

res <- list.files(path = "/path/to/csvFiles")

## Load data from first genome as template
ad <- read.csv(paste0("/path/to/csvFiles/",res[1]), sep=";")
ad$genome <- "HG00438"
ad2 <- read.csv(paste0("/path/to/csvFiles/",res[2]), sep=";")

## Create new columns for the second allele
ad$donor_seq_alt_h2 <-  ad2$donor_seq_alt
ad$alt_SSHW_h2 <- ad2$alt_SSHW
ad$alt_HBS_h2 <- ad2$alt_HBS
ad$HBS_diff_h2 <- ad2$HBS_diff
ad$SSHW_diff_h2 <- ad2$SSHW_diff
ad$n_variations_h2 <- ad2$n_variations

## Create new lines for every genome haplotype
setwd("/path/to/csvFiles/")

for(i in res[-1]){
   
   #i <- res[3]
   fileName <- str_split(i, "_")[[1]][[1]]
   haplotype <- "h1"
   if(str_detect(i, "h2")) haplotype <- "h2"
   
   if(haplotype == "h1"){
      
      curTbl <- read.csv(i, sep=";")
      curTbl$genome <- fileName
      
      i2 <- gsub("h1", "h2", i)
      curTbl2 <- read.csv(i2, sep=";")
      
      curTbl$donor_seq_alt_h2 <-  curTbl2$donor_seq_alt
      curTbl$alt_SSHW_h2 <- curTbl2$alt_SSHW
      curTbl$alt_HBS_h2 <- curTbl2$alt_HBS
      curTbl$HBS_diff_h2 <- curTbl2$HBS_diff
      curTbl$SSHW_diff_h2 <- curTbl2$SSHW_diff
      curTbl$n_variations_h2 <- curTbl2$n_variations
      
      ad <- rbind(ad, curTbl)
   }
}

## Format new clomuns
ad$ref_HBS <- as.numeric(gsub(",","\\.", ad$ref_HBS))
ad$ref_SSHW <- as.numeric(gsub(",","\\.", ad$ref_SSHW))
ad$ref_SSHW <- round(ad$ref_SSHW, 1)
ad$alt_SSHW_h2 <- round(as.numeric(gsub(",","\\.", ad$alt_SSHW_h2)),1)
ad$alt_HBS_h2 <- as.numeric(gsub(",","\\.", ad$alt_HBS_h2))
ad$HBS_diff_h2 <- as.numeric(gsub(",","\\.", ad$HBS_diff_h2))
ad$SSHW_diff_h2 <- round(as.numeric(gsub(",","\\.", ad$SSHW_diff_h2)),1)
ad$alt_SSHW<- round(as.numeric(gsub(",","\\.", ad$alt_SSHW)),1)
ad$alt_HBS <- as.numeric(gsub(",","\\.", ad$alt_HBS))
ad$HBS_diff <- as.numeric(gsub(",","\\.", ad$HBS_diff))
ad$SSHW_diff <- round(as.numeric(gsub(",","\\.", ad$SSHW_diff)),1)

## Add lables to data
ad$sdID <- paste(ad$Chromosome.scaffold.name, ad$Strand, ad$SD_cor)

## Create cluster for parallel processing
library(parallel)
cl <- makeCluster(20)
clusterExport(cl=cl, varlist=c("ad", "SDs"), envir=environment())

ad$TSL1_site <- as.character(parLapply(cl, ad$sdID, function(x){
   
   #x <- sd$sdID[1]
   curentTX <- SDs[SDs$sdID == x,]
   
   if(any(curentTX$TSL1)) y <- "partly TSL1"
   if(all(curentTX$TSL1)) y <- "only TSL1"
   if(all(!curentTX$TSL1)) y <-  "never TSL1"
   
   y
   
}))

clusterExport(cl=cl, varlist=c("ad", "SDs"), envir=environment())

ad$protein_coding <- as.character(lapply( ad$sdID, function(x){
   
   #x <- sd$sdID[1]
   curentTX <- SDs[SDs$sdID == x,]
   
   if(any(curentTX$transcriptType == "protein_coding")) y <- "partly protein_coding"
   if(all(curentTX$transcriptType  == "protein_coding")) y <- "only protein_coding"
   if(all(!curentTX$transcriptType  == "protein_coding")) y <-  "never protein_coding"
   
   y
   
}))

## make all sequences upper case
ad$ref_seq <- toupper(ad$alt_seq)
ad$alt_seq <- toupper(ad$alt_seq)
ad$donor_seq <- toupper(ad$donor_seq)
ad$donor_seq_alt <- toupper(ad$donor_seq_alt)
ad$donor_seq_alt_h2 <- toupper(ad$donor_seq_alt_h2)


ad$homozygot_HBSdiff <- 0
ad$homozygot_HBSdiff[ad$HBS_diff == ad$HBS_diff_h2] <- 1

ad$homozygot_SDvariation <- 0
ad$homozygot_SDvariation[ad$donor_seq_alt != ad$donor_seq & ad$donor_seq_alt_h2 != ad$donor_seq] <- 1

ad$protein_coding[ad$protein_coding == "never protein_coding"] <- "not protein coding"
ad$protein_coding[ad$protein_coding == "only protein_coding"] <- "protein coding"
ad$protein_coding[ad$protein_coding == "partly protein_coding"] <- "protein coding"

ad$TSL1_site[ad$TSL1_site == "never TSL1"] <- " not TSL1"
ad$TSL1_site[ad$TSL1_site == "only TSL1"] <- "TSL1"
ad$TSL1_site[ad$TSL1_site == "partly TSL1"] <- "TSL1"

ad$protein_coding[ad$protein_coding == "not protein coding"] <- "Not protein coding"
ad$protein_coding[ad$protein_coding == "protein coding"] <- "Protein coding"
ad$TSL1_site[ad$TSL1_site == " not TSL1"] <- "Not TSL1"

### Check prevalence of SNV, DEL und INS in genomes around splice sites
table(ad$n_SNV)
table(ad$n_DEL)
table(ad$n_INS)


## Homozygot results
ho <- ad[ad$homozygot_HBSdiff == 1 & ad$homozygot_SDvariation == 1,]
length(unique(ho$sdID))

ho$c <- 1
test <- tapply(ho$c, ho$sdID, sum)
test <- data.frame(test)
test$sdID <- row.names(test)
names(test) <- c("n_indi", "sdID")
test2 <- test

sum(test2$n_indi!= 1)/nrow(test2)

test <- tapply(ho$HBS_diff, ho$sdID, sd)
test <- data.frame(test)
test$sdID <- row.names(test)
names(test) <- c("HBSdiff_sd", "sdID")
test <- test[test$sdID %in% test2$sdID[test2$n_indi > 1],]
test <- na.omit(test)

sum(test$HBSdiff_sd == 0)/nrow(test)

test2$HBSdiff_sd <- test$HBSdiff_sd[match(test2$sdID, test$sdID)]

mean(test2$n_indi)
table(test2$HBSdiff_sd)/sum(table(test2$HBSdiff_sd))

## homo SD group table#
test <- unique(data.frame(id=ho$sdID, prot=ho$protein_coding, tsl = ho$TSL1_site))
table(test$prot, test$tsl)

test <- unique(data.frame(id=ad$sdID, prot=ad$protein_coding, tsl = ad$TSL1_site))
table(test$prot, test$tsl)

## HBS diff per TSL1 and coding labels
library(ggplot2)
g <- ggplot(ho,aes(x=protein_coding, y=HBS_diff, fill=TSL1_site)) + geom_boxplot()+
   ylab("Homozygous HBS difference to the reference")+xlab("") + theme(text = element_text(size = 14)) 
g
ggsave(g, filename="Figure 1 A.png",width=5, height=5)

## Test distribution
NN <- ho$HBS_diff[ho$protein_coding == "Not protein coding" & ho$TSL1_site == "Not TSL1" ]
NT <- ho$HBS_diff[ho$protein_coding == "Not protein coding" & ho$TSL1_site == "TSL1" ]
PN <- ho$HBS_diff[ho$protein_coding == "Protein coding" & ho$TSL1_site == "Not TSL1" ]
PT <- ho$HBS_diff[ho$protein_coding == "Protein coding" & ho$TSL1_site == "TSL1" ]

## compare distributions
ks.test(PT,NN )
ks.test(PT,NT )
ks.test(PT,PN )
ks.test(PT,NN )



### Calcualte MaxEnt diff 
ho$ref_seq_MaxEnt <- substr(ho$donor_seq, 1, 9)
ho$alt_seq_MaxEnt <- substr(ho$donor_seq_alt, 1, 9)

ho_MaxEnt <- ho[ho$ref_seq_MaxEnt != ho$alt_seq_MaxEnt,]

library(ModCon)
ho_MaxEnt$ref_MaxEnt <- calculateMaxEntScanScore(ho_MaxEnt$ref_seq_MaxEnt, 5)
ho_MaxEnt$alt_MaxEnt <- calculateMaxEntScanScore(ho_MaxEnt$alt_seq_MaxEnt, 5)

ho_MaxEnt$ref_MaxEnt <- as.numeric(ho_MaxEnt$ref_MaxEnt)
ho_MaxEnt$alt_MaxEnt <- as.numeric(ho_MaxEnt$alt_MaxEnt)

ho_MaxEnt$MaxEnt_diff <- ho_MaxEnt$alt_MaxEnt - ho_MaxEnt$ref_MaxEnt

ho_MaxEnt_tsl1PC <- ho_MaxEnt[ho_MaxEnt$protein_coding == "Protein coding",]
ho_MaxEnt_tsl1PC <- ho_MaxEnt_tsl1PC[ho_MaxEnt_tsl1PC$TSL1_site ==  "TSL1",]

sum(ho_MaxEnt_tsl1PC$HBS_diff == 0)/nrow(ho_MaxEnt_tsl1PC)
sum(abs(ho_MaxEnt_tsl1PC$MaxEnt_diff) < 1)/nrow(ho_MaxEnt_tsl1PC)

g <- ggplot(ho_MaxEnt, aes(protein_coding, MaxEnt_diff, fill=TSL1_site ))+ geom_boxplot()+
   ylab("Homozygous HBS difference to the reference")+xlab("") +
   scale_fill_brewer(palette="Set2")+
   theme(text = element_text(size=15))
g

ggsave(g, filename="Supplementary Figure 5 B.png",width=5, height=5)






## ref HBS of 5'ss categories
allRef <- ad[ad$genome == "HG00438",]
allRef <- na.omit(allRef)

g <- ggplot(allRef,aes(x=protein_coding, y=ref_HBS, fill=TSL1_site)) + geom_boxplot()+
   ylab("HBS of reference 5'ss")+xlab("") + theme(text = element_text(size = 14)) 
g

ggsave(g, filename="Supplementary Figure 1.png", width = 5, height = 5)

NN <- allRef$ref_HBS[allRef$protein_coding == "Not protein coding" & allRef$TSL1_site == "Not TSL1" ]
NT <- allRef$ref_HBS[allRef$protein_coding == "Not protein coding" & allRef$TSL1_site == "TSL1" ]
PN <- allRef$ref_HBS[allRef$protein_coding == "Protein coding" & allRef$TSL1_site == "Not TSL1" ]
PT <- allRef$ref_HBS[allRef$protein_coding == "Protein coding" & allRef$TSL1_site == "TSL1" ]

shapiro.test(PT[sample(1:length(NN), 5000)])

##distributions are not normal => Mann-Whitney-U-Test
library(rstatix)
wilcox.test(NT,PT)




## Vizialize ref HBS for every exon type
allRef_cat <- ad[ad$genome == ad$genome[1],]

## add information about whether SD coordinates can be found in canonical transcripts of the gene
exons$cano <- 0
exons$cano[exons$Transcript.stable.ID %in% Transcript_startEnd_coord$Transcript.stable.ID[Transcript_startEnd_coord$Ensembl.Canonical == 1]] <- 1

exons$SDcor <- exons$Exon.region.end..bp.
exons$SDcor[exons$Strand == -1] <- exons$Exon.region.start..bp.[exons$Strand == -1]
exons$SD_id <- paste(exons$Chromosome.scaffold.name, exons$Strand, exons$SDcor)

exons$transcriptTSL <- tsl105$Transcript.support.level..TSL.[match(exons$Transcript.stable.ID, tsl105$Transcript.stable.ID)]
exons$TSL1 <- str_detect(exons$transcriptTSL, "tsl1")

exonsTSL1 <- exons[exons$TSL1 == 1,]

coSS <- c(exonsTSL1$SD_id[exonsTSL1$cano == 1] )
nocoSS <- c(exonsTSL1$SD_id[exonsTSL1$cano == 0] )

allRef_cat$canonicalSD <- "Not canonical 5'ss"
allRef_cat$canonicalSD[allRef_cat$sdID %in% coSS] <- "Canonical 5'ss"

## add information if SD coordinate is found exclusivley at the first exon of TSL1 transcripts than not

library(parallel)

cl <- makeCluster(8)

allRef_cat$sequences <- NULL
allRef_cat$ref_seq <- NULL
allRef_cat$alt_seq <- NULL
allRef_cat$donor_seq <- NULL
allRef_cat$donor_seq_alt <- NULL
allRef_cat$refup60 <- NULL
allRef_cat$refdown60 <- NULL
allRef_cat$altup60 <- NULL
allRef_cat$altdown60 <- NULL

allRef_cat <- allRef_cat[allRef_cat$protein_coding == "Protein coding" & allRef_cat$TSL1_site == "TSL1",]

clusterExport(cl, c("allRef_cat", "exonsTSL1"))

## Find 5'ss of micro exons
allRef_cat$micro_exon <- 0
clusterExport(cl, c("allRef_cat"))

allRef_cat$micro_exon <- parLapply(cl, allRef_cat$sdID, function(x){
   
   cur_exons <- exonsTSL1[exonsTSL1$SD_id == x,]
   cur_exons$length <- cur_exons$Exon.region.end..bp. - cur_exons$Exon.region.start..bp. +1
   cur_exons$length3 <- cur_exons$length <= 30
   if(all(cur_exons$length3 )) 1 else 0
   
})

allRef_cat$micro_exon <- unlist(allRef_cat$micro_exon)

allRef_cat$micro_exon[allRef_cat$micro_exon == 0] <- "Exon length > 30nt"
allRef_cat$micro_exon[allRef_cat$micro_exon == 1] <- "Exon length <= 30nt"

library(ggplot2)

g <- ggplot(allRef_cat, aes(micro_exon,   ref_HBS    ))+ geom_boxplot()+
   ylab("HBond score of reference 5'ss")+xlab("") +
   scale_fill_brewer(palette="Set2")+
   theme(text = element_text(size=15))
g

ggsave(g, filename="Supplementary Figure 3 A.png", width = 4, height = 4)




## HBS diff in coding TSL1 SDs
tsl1_ho <- ho[ho$TSL1_site == "TSL1" & ho$protein_coding == "Protein coding",]

tsl1_ho$HBS_diff_group <- cut(tsl1_ho$HBS_diff_h2, breaks=seq(-8,8,0.5) )   

tsl1_hos <- tsl1_ho[,c(38, 44)]
tsl1_hos$c <- 1

tsl1_ho_df <- data.frame(tapply(tsl1_hos$c, tsl1_hos$HBS_diff_group, sum))
tsl1_ho_df$lab <- row.names(tsl1_ho_df)

tsl1_hos <- unique(tsl1_hos)
tsl1_ho_df_unqiue <- data.frame(tapply(tsl1_hos$c, tsl1_hos$HBS_diff_group, sum))

tsl1_ho_df2 <- cbind(tsl1_ho_df, tsl1_ho_df_unqiue)
names(tsl1_ho_df2) <- c("all","lab","unique")

tab1 <- tsl1_ho_df2[,c(1,2)]
tab1$count <- "Total"

tab2 <- tsl1_ho_df2[,c(3,2)]
tab2$all <- tab2$unique
tab2$unique <- NULL
tab2 <- tab2[,c(2,1)]
tab2$count <- "Unique 5'ss"

data2 <- rbind(tab1, tab2)

data2$lab <- factor(data2$lab, levels= tab2$lab)

data2$logAll <- log10(data2$all)

data2$all[data2$all ==  1] <- 1.05


g <- ggplot(data2, aes(fill=count , y= all , x=lab)) + 
   geom_bar(position="identity", stat="identity") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   ylab("Number of 5'ss affected [log10]") + xlab("Group of homozygous HBS difference to the reference")+
   theme(text = element_text(size = 14))+
   scale_fill_manual(name= "lab", values = c("#A6D854", "#8DA0CB"))+
   scale_color_manual(name = "lab", values = c("#A6D854", "#8DA0CB"))+
   scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
   ) +
   annotation_logticks(sides = 'lr') 

g

ggsave(g, filename="Figure 1 B.png",width=8, height=5)

## amount of 5'ss with diff only 0
test <- data.frame(id=tsl1_ho$sdID , diff= tsl1_ho$HBS_diff_h2)
test <- unique(test)
sum(test$diff == 0)/nrow(test)

sum(tsl1_ho$HBS_diff == 0)/nrow(tsl1_ho)


## HBSdiff related to refHBS
lm_eqn <- function(df){
   m <- lm(y ~ x, df);
   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                    list(a = format(unname(coef(m)[1]), digits = 2),
                         b = format(unname(coef(m)[2]), digits = 2),
                         r2 = format(summary(m)$r.squared, digits = 3)))
   as.character(as.expression(eq));
}

tsl1_ho <- ho[ho$TSL1_site == "TSL1" & ho$protein_coding == "Protein coding",]

g <- ggplot(tsl1_ho[tsl1_ho$HBS_diff_h2 != 0,],aes(x=ref_HBS, y=HBS_diff)) + geom_point() +
   geom_smooth(method='lm', formula= y~x)+ylab("HBond score difference to the reference 5'ss")+
   xlab("HBond score of the reference 5'ss")

tsl1_ho$x <- tsl1_ho$ref_HBS
tsl1_ho$y <- tsl1_ho$HBS_diff
lab <- lm_eqn(tsl1_ho[tsl1_ho$HBS_diff != 0,])
lab

p1 <- g + geom_text(x = 17, y = 7.7, label = lab, parse = TRUE)
ggsave(p1, filename="Figure 2 B.png",width=4, height=4)



g <- ggplot(tsl1_ho[tsl1_ho$HBS_diff_h2 != 0,],aes(x=ref_HBS, y=alt_HBS)) + geom_point() +
   geom_smooth(method='lm', formula= y~x)+ylab("Alternative HBond score")+
   xlab("HBond score of the reference 5'ss")

tsl1_ho$x <- tsl1_ho$ref_HBS
tsl1_ho$y <- tsl1_ho$alt_HBS
lab <- lm_eqn(tsl1_ho[tsl1_ho$HBS_diff != 0,])
#lab <- lm_eqn(tsl1_ho)
lab

p1 <- g + geom_text(x = 10, y = 22, label = lab, parse = TRUE)
ggsave(p1, filename="Figure 2 A.png",width=4, height=4)



tsl1_ho_cat <- tsl1_ho

## add information about whether SD coordinates can be found in canonical transcripts of the gene
exons$SDcor <- exons$Exon.region.end..bp.
exons$SDcor[exons$Strand == -1] <- exons$Exon.region.start..bp.[exons$Strand == -1]
exons$SD_id <- paste(exons$Chromosome.scaffold.name, exons$Strand, exons$SDcor)

exons$transcriptTSL <- tsl105$Transcript.support.level..TSL.[match(exons$Transcript.stable.ID, tsl105$Transcript.stable.ID)]
exons$TSL1 <- str_detect(exons$transcriptTSL, "tsl1")

exonsTSL1 <- exons[exons$TSL1 == 1,]

## Find 5'ss of micro exons
tsl1_ho_cat$micro_exon <- 0

tsl1_ho_cat$micro_exon <- lapply(tsl1_ho_cat$sdID, function(x){
   
   cur_exons <- exonsTSL1[exonsTSL1$SD_id == x,]
   cur_exons$length <- cur_exons$Exon.region.end..bp. - cur_exons$Exon.region.start..bp. +1
   cur_exons$length3 <- cur_exons$length <= 30
   if(all(cur_exons$length3 )) 1 else 0
   
})

tsl1_ho_cat$micro_exon <- unlist(tsl1_ho_cat$micro_exon)

tsl1_ho_cat$micro_exon[tsl1_ho_cat$micro_exon == 0] <- "Exon length > 30nt"
tsl1_ho_cat$micro_exon[tsl1_ho_cat$micro_exon == 1] <- "Exon length <= 30nt"

g <- ggplot(tsl1_ho_cat, aes(micro_exon, HBS_diff))+ geom_boxplot()+
   ylab("Homozygous HBS difference to the reference")+xlab("") +
   scale_fill_brewer(palette="Set2")+
   theme(text = element_text(size=15))

tsl1_ho_cat$test <- 0
tsl1_ho_cat$test[tsl1_ho_cat$micro_exon == "Exon length <= 30nt" & 
                    tsl1_ho_cat$exon_position =="All associated exons NOT at first position" &
                    tsl1_ho_cat$canonicalSD == "Canonical 5'ss"] <- 1

unique(tsl1_ho_cat$sdID[tsl1_ho_cat$test == 1])

ggsave(g, filename="Supplementary Figure 4 B.png", width = 5, height = 5)


ks.test(tsl1_ho_cat$HBS_diff[tsl1_ho_cat$micro_exon == "Exon length <= 30nt"],
        tsl1_ho_cat$HBS_diff[tsl1_ho_cat$micro_exon == "Exon length > 30nt"])



## Distribution of affected SD positions

ho_positions <-  lapply(1:nrow(ho), function(x){
   
   ref <- strsplit(ho$donor_seq[x], "")[[1]]
   alt <- strsplit(ho$donor_seq_alt[x], "")[[1]]
   
   which(alt != ref)
   
})

ho_positions <- unlist(ho_positions)
ho_positions_tbl <- table(ho_positions)
ho_positions_tbl2 <- data.frame(ho_positions_tbl)
ho_positions_tbl2$ho_positions <- as.numeric(as.character(ho_positions_tbl2$ho_positions))
ho_positions_tbl2$ho_positions <- as.character(ho_positions_tbl2$ho_positions)

ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 1] <- "-3"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 2] <- "-2"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 3] <- "-1"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 6] <- "+3"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 7] <- "+4"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 8] <- "+5"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 9] <- "+6"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 10] <- "+7"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 11] <- "+8"

ho_positions_tbl3 <- ho_positions_tbl2[c(1:2),]
ho_positions_tbl3$ho_positions <- c("+1","+2")
ho_positions_tbl3$Freq <- 0

ho_positions_tbl2 <- rbind(ho_positions_tbl2, ho_positions_tbl3)

ho_positions_tbl2$ho_positions <- factor(ho_positions_tbl2$ho_positions, levels=c( "-3", "-2", "-1",
                                                                                   "+1", "+2" ,"+3",
                                                                                   "+4", "+5", "+6",
                                                                                   "+7", "+8"  ))

ho_positions_tbl2$percentage <- ho_positions_tbl2$Freq / sum(ho_positions_tbl2$Freq)
ho_positions_tbl2$percentage  <- ho_positions_tbl2$percentage *100
ho_positions_tbl2$percentage <- round(ho_positions_tbl2$percentage, 1)
ho_positions_tbl2$percentage <- paste0(ho_positions_tbl2$percentage, "%")
ho_positions_tbl2$percentage[ho_positions_tbl2$percentage == "0%"] <- ""

g <- ggplot(ho_positions_tbl2, aes(x=ho_positions, y=Freq))+geom_bar(stat="identity")+xlab("5'ss sequence position")+
   ylab("Frequency of nucleotide substitutions")+geom_text(aes(label = percentage), vjust = -0.5)
g
ggsave(g, filename="Supplementary Figure 5 A.png", width = 5, height = 4)


## Does it look similar when analysing only 5ss with only 1 SNV?
ho2 <- ho[ho$SD_ntDiff == 1,]

ho_positions <-  lapply(1:nrow(ho2), function(x){
   
   ref <- strsplit(ho2$donor_seq[x], "")[[1]]
   alt <- strsplit(ho2$donor_seq_alt[x], "")[[1]]
   
   which(alt != ref)
   
})

ho_positions <- unlist(ho_positions)
ho_positions_tbl <- table(ho_positions)
ho_positions_tbl2 <- data.frame(ho_positions_tbl)
ho_positions_tbl2$ho_positions <- as.numeric(as.character(ho_positions_tbl2$ho_positions))
ho_positions_tbl2$ho_positions <- as.character(ho_positions_tbl2$ho_positions)

ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 1] <- "-3"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 2] <- "-2"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 3] <- "-1"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 6] <- "+3"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 7] <- "+4"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 8] <- "+5"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 9] <- "+6"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 10] <- "+7"
ho_positions_tbl2$ho_positions[ho_positions_tbl2$ho_positions == 11] <- "+8"

ho_positions_tbl3 <- ho_positions_tbl2[c(1:2),]
ho_positions_tbl3$ho_positions <- c("+1","+2")
ho_positions_tbl3$Freq <- 0

ho_positions_tbl2 <- rbind(ho_positions_tbl2, ho_positions_tbl3)

ho_positions_tbl2$ho_positions <- factor(ho_positions_tbl2$ho_positions, levels=c( "-3", "-2", "-1",
                                                                                   "+1", "+2" ,"+3",
                                                                                   "+4", "+5", "+6",
                                                                                   "+7", "+8"  ))

ggplot(ho_positions_tbl2, aes(x=ho_positions, y=Freq))+geom_bar(stat="identity")+xlab("5'ss sequence position")+
   ylab("Frequency of nucleotide substitutions")




### Do Figure 2 B but with percentage of possible increase decrease for SDs with single SNV within 5ss

tsl1_ho_single <- data.frame(id=tsl1_ho$sdID, ref_SD_seq=tsl1_ho$donor_seq, 
                             HBSdiff=tsl1_ho$HBS_diff, ntDiff=tsl1_ho$SD_ntDiff,
                             refHBS=tsl1_ho$ref_HBS, altSDseq=tsl1_ho$donor_seq_alt)

tsl1_ho_single <- unique(tsl1_ho_single)
tsl1_ho_single$altHBS <- tsl1_ho_single$refHBS + tsl1_ho_single$HBSdiff

tsl1_ho_single <- tsl1_ho_single[tsl1_ho_single$ntDiff == 1,]

tsl1_ho_single$HBSgroup <- cut(tsl1_ho_single$refHBS, breaks=seq(1.8, 23.8, 1))

tsl1_ho_single$maxHBS <- 0
tsl1_ho_single$minHBS <- 0

for(single in 1:nrow(tsl1_ho_single)){
   
   # single <- 1
   
   refSD <- tsl1_ho_single$ref_SD_seq[single]
   
   curSD <- strsplit(refSD, "")[[1]]
   
   seqlist <- list()
   
   for(pos in c(1:11)){
      
      curSD1 <- curSD
      curSD2 <- curSD
      curSD3 <- curSD
      curSD4 <- curSD
      
      curSD1[pos] <- "A"
      curSD2[pos] <- "G"
      curSD3[pos] <- "C"
      curSD4[pos] <- "T"
      
      seqs <- c(paste(curSD1, collapse = ""),paste(curSD2, collapse = ""),
                paste(curSD3, collapse = ""),paste(curSD4, collapse = ""))
      
      seqs <- seqs[seqs != refSD]
      seqsHBS <- hbg$hbs[match(seqs, hbg$seq)]
      seqsHBS <- na.omit(seqsHBS)
      seqlist[[pos]] <- seqsHBS
      
   }   
   seqlist <- unlist(seqlist)   
   seqlist <- seqlist - tsl1_ho_single$refHBS[single] 
   
   tsl1_ho_single$maxHBS[single] <- max(seqlist)
   tsl1_ho_single$minHBS[single] <- min(seqlist)
   
}

tsl1_ho_single$percentual_change <- 0

for(single in 1:nrow(tsl1_ho_single)){
   
   # single <- 1
   
   maxHBSdiff <- tsl1_ho_single$maxHBS[single]
   minHBSdiff <- tsl1_ho_single$minHBS[single]
   
   if(tsl1_ho_single$HBSdiff[single] > 0) tsl1_ho_single$percentual_change[single] <- round((tsl1_ho_single$HBSdiff[single]/maxHBSdiff)*100)
   if(tsl1_ho_single$HBSdiff[single] < 0) tsl1_ho_single$percentual_change[single] <- round((tsl1_ho_single$HBSdiff[single]/minHBSdiff)*100)*-1
   
}


tsl1_ho_single <- tsl1_ho_single[tsl1_ho_single$refHBS > 9.8 & tsl1_ho_single$refHBS <= 20.8,]
g <- ggplot(tsl1_ho_single[tsl1_ho_single$HBSdiff != 0,], aes(y= percentual_change , x= HBSgroup)) + 
   geom_boxplot( ) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   ylab("Percentual HBS difference relative to possible maxima [%]") + xlab("Hbond score group of reference 5'ss")+
   theme(text = element_text(size = 14))


g

ggsave(g, filename="Figure 2 C.png", width = 5, height = 6.1)

cor(tsl1_ho_single$refHBS, tsl1_ho_single$percentual_change)

## spearson correlation test
corr <- cor.test(x=tsl1_ho_single$refHBS, y=tsl1_ho_single$percentual_change, method = 'spearman')




## Statistically test, if SD diff is significantly shifted from random distribution
tsl1_ho$SD_ntDiff <- as.numeric(lapply(1:nrow(tsl1_ho), function(x){
   
   ref <- strsplit(tsl1_ho$donor_seq[x], "")[[1]]
   alt <- strsplit(tsl1_ho$donor_seq_alt[x], "")[[1]]
   
   sum(alt != ref)
   
}))

tsl1_ho_single <- data.frame(id=tsl1_ho$sdID, ref_SD_seq=tsl1_ho$donor_seq, 
                             HBSdiff=tsl1_ho$HBS_diff, ntDiff=tsl1_ho$SD_ntDiff,
                             refHBS=tsl1_ho$ref_HBS, altSDseq=tsl1_ho$donor_seq_alt)

tsl1_ho_single <- unique(tsl1_ho_single)
tsl1_ho_single$altHBS <- tsl1_ho_single$refHBS + tsl1_ho_single$HBSdiff

tsl1_ho_single <- tsl1_ho_single[tsl1_ho_single$ntDiff == 1,]


tsl1_ho_single$HBSgroup <- cut(tsl1_ho_single$refHBS, breaks=seq(1.8, 23.8, 1))

g <- table(tsl1_ho_single$HBSgroup)

test <- tsl1_ho_single[tsl1_ho_single$HBSgroup %in% names(g[g>10]),]

allVars <- list()

for(group in unique(tsl1_ho_single$HBSgroup)){
   
   #group <- unique(tsl1_ho_single$HBSgroup)[1]
   
   varLIST <- list()
   
   cur <- tsl1_ho_single[tsl1_ho_single$HBSgroup == group,]
   
   for(i in 1:nrow(cur)){
      
      #i <- 1
      sdSeq <- cur$ref_SD_seq[i]
      
      curSD <- strsplit(sdSeq, "")[[1]]
      
      seqlist <- list()
      
      for(pos in c(1:11)){
         
         curSD1 <- curSD
         curSD2 <- curSD
         curSD3 <- curSD
         curSD4 <- curSD
         
         curSD1[pos] <- "A"
         curSD2[pos] <- "G"
         curSD3[pos] <- "C"
         curSD4[pos] <- "T"
         
         seqs <- c(paste(curSD1, collapse = ""),paste(curSD2, collapse = ""),
                   paste(curSD3, collapse = ""),paste(curSD4, collapse = ""))
         
         seqs <- seqs[seqs != sdSeq]
         seqsHBS <- hbg$hbs[match(seqs, hbg$seq)]
         seqsHBS <- na.omit(seqsHBS)
         seqsHBSdiff <- seqsHBS-cur$refHBS[i]
         
         seqlist[[pos]] <- seqsHBSdiff
         
      }
      varLIST[[paste(sdSeq, i, collapse = "_")]] <- unlist(seqlist)
   }
   
   x <- unlist(varLIST)
   names(x) <- NULL
   allVars[[group]] <- x
   
}

all <- unlist(allVars)

varSDs <- data.frame(hbs_diff= as.numeric(all), group = names(all))
varSDs$group <- as.character(lapply(varSDs$group, function(x){
   
   paste(strsplit(x, "]")[[1]][[1]], "]", collapse = "")
   
}))

varSDs$group <- gsub(" ", "", varSDs$group)


res <- list()

for(i in unique(varSDs$group)){
   
   ## i <- "(14.8,15.8]"
   
   background <- varSDs$hbs_diff[varSDs$group == i]
   real <- tsl1_ho_single$HBSdiff[tsl1_ho_single$HBSgroup == i]
   
   background <- data.frame(diff= as.numeric(background), diffgroup="", c=1)
   background$diffgroup <-  cut(background$diff, breaks=c(-23.8, -1, 1, 23.8))
   backgroundTbl <- data.frame(tapply(background$c, background$diffgroup, sum))
   
   real <- data.frame(diff= as.numeric(real), diffgroup="", c=1)
   real$diffgroup <-  cut(real$diff, breaks=c(-23.8, -1, 1, 23.8))
   realTbl <- data.frame(tapply(real$c, real$diffgroup, sum))
   
   realTbl$typ <- "real"
   backgroundTbl$type <- "background"
   
   dat <- data.frame(
      "real" = realTbl[,1],
      "background" = backgroundTbl[,1],
      row.names = row.names(realTbl),
      stringsAsFactors = FALSE
   )
   
   dat[is.na(dat)] <- 0
   test <- fisher.test(dat)
   
   mosaicplot(dat,
              main = i,
              color = TRUE
   )
   
   
   res[[i]] <- test$p.value
   
   
}

check <- data.frame(table(tsl1_ho_single$HBSgroup))

resDF <- data.frame(group=names(res), p=unlist(res))
check$pvalue <- resDF$p[match(check$Var1, resDF$group )]



tbl <- table(tsl1_ho_single$HBSgroup)


res <- list()

for(i in names(tbl[tbl>1])){
   
   ## i <- "(14.8,15.8]"
   
   background <- varSDs$hbs_diff[varSDs$group == i]
   real <- tsl1_ho_single$HBSdiff[tsl1_ho_single$HBSgroup == i]
   
   background <- data.frame(diff= as.numeric(background), diffgroup="", c=1)
   background$diffgroup <-  cut(background$diff, breaks=c(-23.8, -1, 23.8))
   backgroundTbl <- data.frame(tapply(background$c, background$diffgroup, sum))
   
   real <- data.frame(diff= as.numeric(real), diffgroup="", c=1)
   real$diffgroup <-  cut(real$diff, breaks=c(-23.8, -1, 23.8))
   realTbl <- data.frame(tapply(real$c, real$diffgroup, sum))
   
   realTbl$typ <- "real"
   backgroundTbl$type <- "background"
   
   dat <- data.frame(
      "real" = realTbl[,1],
      "background" = backgroundTbl[,1],
      row.names = row.names(realTbl),
      stringsAsFactors = FALSE
   )
   
   dat[is.na(dat)] <- 0
   test <- fisher.test(dat)
   
   mosaicplot(dat,
              main = i,
              color = TRUE
   )
   
   
   res[[i]] <- test$p.value
   
   
}

check <- data.frame(table(tsl1_ho_single$HBSgroup))

resDF <- data.frame(group=names(res), p=unlist(res))
check$pvalue <- resDF$p[match(check$Var1, resDF$group )]

check$percentage <- round(check$Freq/sum(check$Freq)*100,1)


plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)

file.copy(from=plots.png.paths, to="/saveToPath/haplotype genome SNVs")

write.csv2(resDF, file="FisherExact Pvales for suppl. table 1.csv", row.names=F)


## show example
group <-"(14.8,15.8]"

varLIST <- list()
varLIST2 <- list()

cur <- tsl1_ho_single[tsl1_ho_single$HBSgroup == group,]

for(i in 1:nrow(cur)){
   
   #i <- 164
   sdSeq <- cur$ref_SD_seq[i]
   
   curSD <- strsplit(sdSeq, "")[[1]]
   
   seqlist <- list()
   seqlist2 <- list()
   
   for(pos in c(1:11)){
      
      curSD1 <- curSD
      curSD2 <- curSD
      curSD3 <- curSD
      curSD4 <- curSD
      
      curSD1[pos] <- "A"
      curSD2[pos] <- "G"
      curSD3[pos] <- "C"
      curSD4[pos] <- "T"
      
      seqs <- c(paste(curSD1, collapse = ""),paste(curSD2, collapse = ""),
                paste(curSD3, collapse = ""),paste(curSD4, collapse = ""))
      
      seqs <- seqs[seqs != sdSeq]
      seqsHBS <- hbg$hbs[match(seqs, hbg$seq)]
      seqsHBS <- na.omit(seqsHBS)
      seqsHBSdiff <- seqsHBS-cur$refHBS[i]
      
      seqlist[[pos]] <- seqsHBSdiff
      seqlist2[[pos]] <- seqs[substr(seqs, 4, 5) == "GT"]
      
   }
   varLIST[[paste(sdSeq, i, collapse = "_")]] <- unlist(seqlist)
   varLIST2[[paste(sdSeq, cur$id[i], i, collapse = "_")]] <- unlist(seqlist2)
}


res_testDF <- list()

for(i in 1:length(varLIST2)){
   
   #i <- 1
   
   testDF <- data.frame(ref_5ss_chromosome="", ref_5ss_strand="", ref_5ss_coordinate="", 
                        ref_5ss_seq="", reference_5ss_HBS=0, random_5ss_SNV_variation=varLIST2[[i]])
   refSDseq <- str_split(names(varLIST2)[i], " ")[[1]][[1]]
   refSDchr <- str_split(names(varLIST2)[i], " ")[[1]][[2]]
   refSDcor <- str_split(names(varLIST2)[i], " ")[[1]][[4]]
   refSDstrand<- str_split(names(varLIST2)[i], " ")[[1]][[3]]
   
   testDF$ref_5ss_chromosome <-   refSDchr
   testDF$ref_5ss_strand <-   refSDstrand
   testDF$ref_5ss_coordinate <-   refSDcor
   testDF$ref_5ss_seq <-   refSDseq
   
   testDF$reference_5ss_HBS <- hbg$hbs[hbg$seq ==  refSDseq]
   testDF$random_5ss_SNV_variation_position <- rep(c(1,2,3,6,7,8,9,10,11), each=3)
   
   res_testDF[[i]] <- testDF
}

library(data.table)
res_testDF2 <- rbindlist(res_testDF)

res_testDF2$random_5ss_SNV_variation_HBS <- hbg$hbs[match(res_testDF2$random_5ss_SNV_variation, hbg$seq)]

res_testDF2$HBS_difference  <- res_testDF2$random_5ss_SNV_variation_HBS - res_testDF2$reference_5ss_HBS

res_testDF2$diffgroup <-  cut(res_testDF2$HBS_difference, breaks=c(-23.8, -1, 23.8))
res_testDF2$c <- 1
res_testDF2Tbl <- data.frame(tapply(res_testDF2$c, res_testDF2$diffgroup, sum))

# tapply.res_testDF2.c..res_testDF2.diffgroup..sum.
# (-23.8,-1]                                              4879
# (-1,23.8]                                               4247

## i <- "(14.8,15.8]"



realDFL <- list()

for(i in 1:nrow(tsl1_ho_single)){
   
   #i <- 1
   
   realDF <- data.frame(ref_5ss_chromosome="", ref_5ss_strand="", ref_5ss_coordinate="", 
                        ref_5ss_seq="", reference_5ss_HBS=0, observed_5ss_seq_variation= "")
   
   refSDseq <- tsl1_ho_single$ref_SD_seq[i]
   refSDchr <- str_split(tsl1_ho_single$id[i], " ")[[1]][[3]]
   refSDcor <- str_split(tsl1_ho_single$id[i], " ")[[1]][[1]]
   refSDstrand<- str_split(tsl1_ho_single$id[i], " ")[[1]][[2]]
   
   realDF$ref_5ss_chromosome <-   refSDchr
   realDF$ref_5ss_strand <-   refSDstrand
   realDF$ref_5ss_coordinate <-   refSDcor
   realDF$ref_5ss_seq <-   refSDseq
   
   realDF$observed_5ss_seq_variation <- tsl1_ho_single$altSDseq[i]
   
   realDF$reference_5ss_HBS <- hbg$hbs[hbg$seq ==  refSDseq]
   
   ref <- strsplit(refSDseq, "")[[1]]
   alt <- strsplit(tsl1_ho_single$altSDseq[i], "")[[1]]
   
   posi <- which(alt != ref)
   
   realDF$random_5ss_SNV_variation_position <- posi
   
   realDFL[[i]] <- realDF
   
}

realDFL2 <- rbindlist(realDFL)

realDFL2$altHBS <- hbg$hbs[match(realDFL2$observed_5ss_seq_variation, hbg$seq)]
realDFL2$HBS_diff <- realDFL2$altHBS - realDFL2$reference_5ss_HBS
realDFL2$c <- 1
realDFL2$diffgroup <-  cut(realDFL2$HBS_diff, breaks=c(-23.8, -1, 23.8))
realTbl <- data.frame(tapply(realDFL2$c, realDFL2$diffgroup, sum))

realDFL2 <- realDFL2[realDFL2$reference_5ss_HBS > 14.8 & realDFL2$reference_5ss_HBS <= 15.9,]


write.csv2(res_testDF2, file="Supplementary table 1 A.csv", row.names=F)
write.csv2(realDFL2, file="Supplementary table 1 B.csv", row.names=F)




## Heterozygos variations


select  <- ((ad$donor_seq_alt_h2 != ad$donor_seq) & (ad$donor_seq_alt != ad$donor_seq) & (ad$HBS_diff != ad$HBS_diff_h2))
select1 <- ((ad$donor_seq_alt_h2 == ad$donor_seq) & (ad$donor_seq_alt != ad$donor_seq))
select2 <- ((ad$donor_seq_alt_h2 != ad$donor_seq) & (ad$donor_seq_alt == ad$donor_seq))


ad$heterozygot_SDvariation <- 0
ad$heterozygot_SDvariation[select | select1 | select2 ] <- 1

ad$het_allele_diff <- ""
ad$het_allele_diff[select] <- "both"
ad$het_allele_diff[select1] <- "one"
ad$het_allele_diff[select2] <- "one"

he <- ad[ ad$heterozygot_SDvariation == 1,]

# table(he$protein_coding, he$TSL1_site)
length(unique(he$sdID))

he$c <- 1
test <- tapply(he$c, he$sdID, sum)
test <- data.frame(test)
test$sdID <- row.names(test)
names(test) <- c("n_indi", "sdID")
test2 <- test

test <- tapply(he$HBS_diff, he$sdID, sd)
test <- data.frame(test)
test$sdID <- row.names(test)
names(test) <- c("HBSdiff_sd", "sdID")

test2$HBSdiff_sd <- test$HBSdiff_sd[match(test2$sdID, test$sdID)]

mean(test2$n_indi)
table(test2$HBSdiff_sd)/sum(table(test2$HBSdiff_sd))

sum(test2$n_indi > 1)/nrow(test2)
test3 <- na.omit(test2)
sum(test3$HBSdiff_sd == 0)/nrow(test2)

## homo SD group table#
test <- unique(data.frame(id=he$sdID, prot=he$protein_coding, tsl = he$TSL1_site))
table(test$prot, test$tsl)

test <- unique(data.frame(id=ad$sdID, prot=ad$protein_coding, tsl = ad$TSL1_site))
table(test$prot, test$tsl)


6173  /83910
1001/11404
3082 /42920
7200/173199



## Same analysis for heterzyous changes as for homozygous

tsl1_he1 <- he[he$het_allele_diff == "one",]

IndexVector1 <- (tsl1_he1$donor_seq_alt != tsl1_he1$donor_seq & 
                    tsl1_he1$donor_seq_alt_h2 == tsl1_he1$donor_seq)
IndexVector2 <- (tsl1_he1$donor_seq_alt == tsl1_he1$donor_seq &
                    tsl1_he1$donor_seq_alt_h2 != tsl1_he1$donor_seq)



length(unique(tsl1_he1$sdID))

tsl1_he_df <- data.frame(HBS_diff= c(tsl1_he1$HBS_diff[IndexVector1],
                                     tsl1_he1$HBS_diff_h2[IndexVector2]),
                         protein_coding= c(tsl1_he1$protein_coding[IndexVector1],
                                           tsl1_he1$protein_coding[IndexVector2]),
                         TSL1_site= c(tsl1_he1$TSL1_site[IndexVector1],
                                      tsl1_he1$TSL1_site[IndexVector2]))

tsl1_he_df <- na.omit(tsl1_he_df)

g <- ggplot(tsl1_he_df,aes(x=protein_coding, y=HBS_diff, fill=TSL1_site)) + geom_boxplot()+
   ylab("Heterozygous HBS difference to the reference")+xlab("") + theme(text = element_text(size = 14))  +
   ylim(-20,12)
ggsave(g, filename="Figure 3 A.png",width=5, height=5)




tsl1_he1 <- he[he$het_allele_diff == "both",]

length(unique(tsl1_he1$sdID))

plot(density(tsl1_he1$HBS_diff - tsl1_he1$HBS_diff_h2))


tsl1_he_df <- data.frame(HBS_diff= c(tsl1_he1$HBS_diff,
                                     tsl1_he1$HBS_diff_h2),
                         protein_coding= c(tsl1_he1$protein_coding,
                                           tsl1_he1$protein_coding),
                         TSL1_site= c(tsl1_he1$TSL1_site,
                                      tsl1_he1$TSL1_site))


g <- ggplot(tsl1_he_df,aes(x=protein_coding, y=HBS_diff, fill=TSL1_site)) + geom_boxplot()+
   ylab("Heterozygous HBS difference to the reference")+xlab("") + theme(text = element_text(size = 14)) +
   ylim(-20,12)
g
ggsave(g, filename="Figure 3 B.png",width=5, height=5)




## Test distribution
NN <- tsl1_he_df$HBS_diff[tsl1_he_df$protein_coding == "Not protein coding" & tsl1_he_df$TSL1_site == "Not TSL1" ]
NT <- tsl1_he_df$HBS_diff[tsl1_he_df$protein_coding == "Not protein coding" & tsl1_he_df$TSL1_site == "TSL1" ]
PN <- tsl1_he_df$HBS_diff[tsl1_he_df$protein_coding == "Protein coding" & tsl1_he_df$TSL1_site == "Not TSL1" ]
PT <- tsl1_he_df$HBS_diff[tsl1_he_df$protein_coding == "Protein coding" & tsl1_he_df$TSL1_site == "TSL1" ]

## compare distributions
ks.test(PT,NN )
ks.test(PT,NT )
ks.test(PT,PN )
ks.test(PT,NN )

PN_he <- PN
PT_he <- PT

ks.test(PN_he,PN )
ks.test(PT_he,PT )






## HBSdiff related to refHBS
lm_eqn <- function(df){
   m <- lm(y ~ x, df);
   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                    list(a = format(unname(coef(m)[1]), digits = 2),
                         b = format(unname(coef(m)[2]), digits = 2),
                         r2 = format(summary(m)$r.squared, digits = 3)))
   as.character(as.expression(eq));
}




tsl1_he <- he[he$het_allele_diff == "one",]
tsl1_he <- tsl1_he[tsl1_he$TSL1_site == "TSL1" & tsl1_he$protein_coding == "Protein coding",]

tsl1_he1 <- tsl1_he

IndexVector2 <- (tsl1_he1$donor_seq_alt == tsl1_he1$donor_seq &
                    tsl1_he1$donor_seq_alt_h2 != tsl1_he1$donor_seq)
IndexVector1 <- (tsl1_he1$donor_seq_alt != tsl1_he1$donor_seq & 
                    tsl1_he1$donor_seq_alt_h2 == tsl1_he1$donor_seq)


length(unique(tsl1_he1$sdID))

tsl1_he_df <- data.frame(HBS_diff= c(tsl1_he1$HBS_diff[IndexVector1],
                                     tsl1_he1$HBS_diff_h2[IndexVector2]),
                         protein_coding= c(tsl1_he1$protein_coding[IndexVector1],
                                           tsl1_he1$protein_coding[IndexVector2]),
                         TSL1_site= c(tsl1_he1$TSL1_site[IndexVector1],
                                      tsl1_he1$TSL1_site[IndexVector2]),
                         refHBS= c(tsl1_he1$ref_HBS[IndexVector1],
                                   tsl1_he1$ref_HBS[IndexVector2]),
                         altSeq1=c(tsl1_he1$donor_seq_alt[IndexVector1],
                                   tsl1_he1$donor_seq_alt_h2[IndexVector2]),
                         refSeq=c(tsl1_he1$donor_seq[IndexVector1],
                                  tsl1_he1$donor_seq[IndexVector2]),
                         sdID= c(tsl1_he1$sdID[IndexVector1],
                                 tsl1_he1$sdID[IndexVector2]))

g <- ggplot(tsl1_he_df[tsl1_he_df$HBS_diff != 0,],aes(x=refHBS, y=HBS_diff)) + geom_point() +
   geom_smooth(method='lm', formula= y~x)+ylab("HBond score difference to the reference 5'ss")+
   xlab("HBond score of the reference 5'ss")+ylim(-15,8)+xlim(5,23.8)

tsl1_he_df$x <- tsl1_he_df$refHBS
tsl1_he_df$y <- tsl1_he_df$HBS_diff
lab <- lm_eqn(tsl1_he_df[tsl1_he_df$HBS_diff != 0,])
lab

p1 <- g + geom_text(x = 17, y = 8, label = lab, parse = TRUE)
ggsave(p1, filename="Figure 4 A.png",width=4, height=4)




tsl1_he <- he[he$het_allele_diff == "both",]
tsl1_he <- tsl1_he[tsl1_he$TSL1_site == "TSL1" & tsl1_he$protein_coding == "Protein coding",]


tsl1_he1 <- tsl1_he

IndexVector2 <- T
IndexVector1 <- T


length(unique(tsl1_he1$sdID))

tsl1_he_df <- data.frame(HBS_diff= c(tsl1_he1$HBS_diff[IndexVector1],
                                     tsl1_he1$HBS_diff_h2[IndexVector2]),
                         protein_coding= c(tsl1_he1$protein_coding[IndexVector1],
                                           tsl1_he1$protein_coding[IndexVector2]),
                         TSL1_site= c(tsl1_he1$TSL1_site[IndexVector1],
                                      tsl1_he1$TSL1_site[IndexVector2]),
                         refHBS= c(tsl1_he1$ref_HBS[IndexVector1],
                                   tsl1_he1$ref_HBS[IndexVector2]),
                         altSeq1=c(tsl1_he1$donor_seq_alt[IndexVector1],
                                   tsl1_he1$donor_seq_alt_h2[IndexVector2]),
                         refSeq=c(tsl1_he1$donor_seq[IndexVector1],
                                  tsl1_he1$donor_seq[IndexVector2]),
                         sdID= c(tsl1_he1$sdID[IndexVector1],
                                 tsl1_he1$sdID[IndexVector2]))

g <- ggplot(tsl1_he_df[tsl1_he_df$HBS_diff != 0,],aes(x=refHBS, y=HBS_diff)) + geom_point() +
   geom_smooth(method='lm', formula= y~x)+ylab("HBond score difference to the reference 5'ss")+
   xlab("HBond score of the reference 5'ss")+ylim(-15,8)+xlim(5,23.8)

tsl1_he_df$x <- tsl1_he_df$refHBS
tsl1_he_df$y <- tsl1_he_df$HBS_diff
lab <- lm_eqn(tsl1_he_df[tsl1_he_df$HBS_diff != 0,])
lab

p1 <- g + geom_text(x = 17, y = 8, label = lab, parse = TRUE)
ggsave(p1, filename="Figure 4 B.png",width=4, height=4)





## HBS diff in coding TSL1 SDs
tsl1_he_df <- tsl1_he_df[tsl1_he_df$TSL1_site == "TSL1" & tsl1_he_df$protein_coding == "Protein coding",]
tsl1_he_df_save <- tsl1_he_df

tsl1_he_df$HBS_diff_group <- cut(tsl1_he_df$HBS_diff, breaks=seq(-15,8,0.5) )   

tsl1_hes <- tsl1_he_df[,c(7, 8)]
tsl1_hes$c <- 1

tsl1_he_df <- data.frame(tapply(tsl1_hes$c, tsl1_hes$HBS_diff_group, sum))
tsl1_he_df$lab <- row.names(tsl1_he_df)

tsl1_hes <- unique(tsl1_hes)
tsl1_he_df_unqiue <- data.frame(tapply(tsl1_hes$c, tsl1_hes$HBS_diff_group, sum))

tsl1_he_df2 <- cbind(tsl1_he_df, tsl1_he_df_unqiue)
names(tsl1_he_df2) <- c("all","lab","unique")

tab1 <- tsl1_he_df2[,c(1,2)]
tab1$count <- "Total"

tab2 <- tsl1_he_df2[,c(3,2)]
tab2$all <- tab2$unique
tab2$unique <- NULL
tab2 <- tab2[,c(2,1)]
tab2$count <- "Unique 5'ss"

data2 <- rbind(tab1, tab2)

data2$lab <- factor(data2$lab, levels= tab2$lab)

data2$logAll <- log10(data2$all)

data2$all[data2$all ==  1] <- 1.05


g <- ggplot(data2, aes(fill=count , y= all , x=lab)) + 
   geom_bar(position="identity", stat="identity") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   ylab("Number of 5'ss affected [log10]") + xlab("Group of heterozygous HBS difference to the reference")+
   theme(text = element_text(size = 14))+
   scale_fill_manual(name= "lab", values = c("#A6D854", "#8DA0CB"))+
   scale_color_manual(name = "lab", values = c("#A6D854", "#8DA0CB"))+
   scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
   ) +
   annotation_logticks(sides = 'lr')  +
   theme(legend.title = element_blank())

ggsave(g, filename="Supplementary Figure 2.png", width = 8, height = 5)



ks.test(tsl1_he_df_save$HBS_diff, tsl1_ho$HBS_diff)

sd(tsl1_he_df_save$HBS_diff, na.rm=T)





## Check for compensating mutations

## Statistically test, if SD diff is significantly shifted from random distribution
tsl1_he_df$SD_ntDiff <- as.numeric(lapply(1:nrow(tsl1_he_df), function(x){
   
   ref <- strsplit(tsl1_he_df$refSeq[x], "")[[1]]
   alt <- strsplit(tsl1_he_df$altSeq1[x], "")[[1]]
   
   sum(alt != ref)
   
}))


tsl1_he_df_2 <- tsl1_he_df[tsl1_he_df$SD_ntDiff == 2,]
tsl1_ho_2 <- ho[ho$SD_ntDiff > 2,]
i <- 1

tsl1_ho_2$snv1_change <- 0

for(i in 1:nrow(tsl1_ho_2)){
   
   refSeq <- tsl1_ho_2$donor_seq[i]
   refSeq_split <- strsplit(refSeq, "")[[1]]
   altSeq <- tsl1_ho_2$donor_seq_alt[i]
   altSeq_split <- strsplit(altSeq, "")[[1]]
   
   changedPos <- which((altSeq_split != refSeq_split))
   changedPosDiff <- list()
   snv1 <- refSeq_split
   
   for(change in changedPos){
      
      snv2 <- snv1
      snv2[change] <- altSeq_split[change]
      snv2HBS <- hbg$hbs[hbg$seq == paste(snv2, collapse = "")]
      changedPosDiff[[change]] <- snv2HBS - refHBS
   }
   
   changedPosDiff <- unlist(changedPosDiff)
   
   
   tsl1_ho_2$snv1_change[i] <- changedPosDiff[which(abs(changedPosDiff) == max(abs(changedPosDiff)))]
   
}

tsl1_ho_2$id <- 1:nrow(tsl1_ho_2)

tsl1_ho_2$overallChange <- ""
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change < 0 & tsl1_ho_2$HBS_diff > tsl1_ho_2$snv1_change] <- "Reduction compensated"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change > 0 & tsl1_ho_2$HBS_diff < tsl1_ho_2$snv1_change] <- "Increase compensated"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change < 0 & tsl1_ho_2$HBS_diff >= 0] <- "Reduction compensated to 0 or greater"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change > 0 & tsl1_ho_2$HBS_diff <= 0] <- "Increase compensated to 0 or less"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change < 0 & tsl1_ho_2$HBS_diff < tsl1_ho_2$snv1_change] <- "Reduction enhanced"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change > 0 & tsl1_ho_2$HBS_diff > tsl1_ho_2$snv1_change] <- "Increase enhanced"
# tsl1_ho_2$overallChange[(tsl1_ho_2$snv1_change > 0 & tsl1_ho_2$HBS_diff < tsl1_ho_2$snv1_change)|
#                            (tsl1_ho_2$snv1_change < 0 & tsl1_ho_2$HBS_diff > 0)] <- "Increasing and Reducing"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change == 0 & tsl1_ho_2$HBS_diff == 0] <- "No change"

tsl1_ho_2_2 <- tsl1_ho_2[,c(14,38,47,49,28)]
tsl1_ho_2_2 <- unique(tsl1_ho_2_2)

tsl1_ho_2_2$id <- 1:nrow(tsl1_ho_2_2)

v1 <- data.frame(HBSdiff = 0, id=tsl1_ho_2_2$id, type="Ref HBS\nstart",overallChange=tsl1_ho_2_2$overallChange)
v2 <- data.frame(HBSdiff = tsl1_ho_2_2$snv1_change, id=tsl1_ho_2_2$id, type="Strongest SNV \nHBS change",overallChange=tsl1_ho_2_2$overallChange)
v3 <- data.frame(HBSdiff = tsl1_ho_2_2$HBS_diff, id=tsl1_ho_2_2$id, type="Total HBS\ndifference",overallChange=tsl1_ho_2_2$overallChange)

dataTbl <- rbind(v1, v2, v3)

dataTbl$type <- factor(dataTbl$type, levels=c("Ref HBS\nstart", "Strongest SNV \nHBS change", "Total HBS\ndifference"))

g <- ggplot(dataTbl, aes(x = type, y = HBSdiff, group = id)) + 
   geom_point() + geom_line() + xlab("") + ylab("")
g <- g + facet_grid(. ~ dataTbl$overallChange, scales="free", shrink=T)+theme(plot.title = element_text(hjust = 0.5))

ggsave(g, filename="Supplementary Figure 6 A.png", width = 17, height = 4)

outTBLho <- table(dataTbl$overallChange)
outTBLho <- data.frame(outTBLho)



tsl1_he_df_2 <- tsl1_he_df[tsl1_he_df$SD_ntDiff > 2,]
tsl1_ho_2 <-tsl1_he_df_2
tsl1_ho_2 <- na.omit(tsl1_ho_2)
i <- 1

tsl1_ho_2$snv1_change <- 0

for(i in 1:nrow(tsl1_ho_2)){
   
   refSeq <- tsl1_ho_2$refSeq[i]
   refSeq_split <- strsplit(refSeq, "")[[1]]
   altSeq <- tsl1_ho_2$altSeq1[i]
   altSeq_split <- strsplit(altSeq, "")[[1]]
   
   changedPos <- which((altSeq_split != refSeq_split))
   changedPosDiff <- list()
   snv1 <- refSeq_split
   
   for(change in changedPos){
      
      snv2 <- snv1
      snv2[change] <- altSeq_split[change]
      snv2HBS <- hbg$hbs[hbg$seq == paste(snv2, collapse = "")]
      changedPosDiff[[change]] <- snv2HBS - refHBS
   }
   
   changedPosDiff <- unlist(changedPosDiff)
   
   
   tsl1_ho_2$snv1_change[i] <- changedPosDiff[which(abs(changedPosDiff) == max(abs(changedPosDiff)))]
   
}

tsl1_ho_2$id <- 1:nrow(tsl1_ho_2)


tsl1_ho_2$overallChange <- ""
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change < 0 & tsl1_ho_2$HBS_diff > tsl1_ho_2$snv1_change] <- "Reduction compensated"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change > 0 & tsl1_ho_2$HBS_diff < tsl1_ho_2$snv1_change] <- "Increase compensated"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change < 0 & tsl1_ho_2$HBS_diff >= 0] <- "Reduction compensated to 0 or greater"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change > 0 & tsl1_ho_2$HBS_diff <= 0] <- "Increase compensated to 0 or less"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change < 0 & tsl1_ho_2$HBS_diff < tsl1_ho_2$snv1_change] <- "Reduction enhanced"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change > 0 & tsl1_ho_2$HBS_diff > tsl1_ho_2$snv1_change] <- "Increase enhanced"
# tsl1_ho_2$overallChange[(tsl1_ho_2$snv1_change > 0 & tsl1_ho_2$HBS_diff < tsl1_ho_2$snv1_change)|
#                            (tsl1_ho_2$snv1_change < 0 & tsl1_ho_2$HBS_diff > 0)] <- "Increasing and Reducing"
tsl1_ho_2$overallChange[tsl1_ho_2$snv1_change == 0 & tsl1_ho_2$HBS_diff == 0] <- "No change"

tsl1_ho_2_2 <- tsl1_ho_2[,c(5,1,7,9,11)]
tsl1_ho_2_2 <- unique(tsl1_ho_2_2)

tsl1_ho_2_2$id <- 1:nrow(tsl1_ho_2_2)

v1 <- data.frame(HBSdiff = 0, id=tsl1_ho_2_2$id, type="Ref HBS\nstart",overallChange=tsl1_ho_2_2$overallChange)
v2 <- data.frame(HBSdiff = tsl1_ho_2_2$snv1_change, id=tsl1_ho_2_2$id, type="Strongest SNV \nHBS change",overallChange=tsl1_ho_2_2$overallChange)
v3 <- data.frame(HBSdiff = tsl1_ho_2_2$HBS_diff, id=tsl1_ho_2_2$id, type="Total HBS\ndifference",overallChange=tsl1_ho_2_2$overallChange)

dataTbl <- rbind(v1, v2, v3)

dataTbl$type <- factor(dataTbl$type, levels=c("Ref HBS\nstart", "Strongest SNV \nHBS change", "Total HBS\ndifference"))

g <- ggplot(dataTbl, aes(x = type, y = HBSdiff, group = id)) + 
   geom_point() + geom_line() + xlab("") + ylab("")
g <- g + facet_grid(. ~ dataTbl$overallChange, scales="free", shrink=T)+theme(plot.title = element_text(hjust = 0.5))

ggsave(g, filename="Supplementary Figure 6 B.png", width = 17, height = 4)


outTBL <- table(dataTbl$overallChange)
outTBL <- data.frame(outTBL)



outTBLho$Homozygot <- outTBLho$Freq
outTBLho$Heterozygot <- outTBL$Freq

outTBLho$Freq <- NULL


library(ggplot2)

data_long <- reshape2::melt(outTBLho, id.vars = "Var1")

data_long$Var1 <- as.character(data_long$Var1)

data_long$Var1[data_long$Var1 == "Increase compensated"] <- "Increase\ncompensated"
data_long$Var1[data_long$Var1 == "Increase compensated to 0 or less"] <- "Increase compensated\nto 0 or less"
data_long$Var1[data_long$Var1 == "Increase enhanced"] <- "Increase\nenhanced"
data_long$Var1[data_long$Var1 == "Reduction compensated"] <- "Reduction\ncompensated"
data_long$Var1[data_long$Var1 == "Reduction compensated to 0 or greater"] <- "Reduction compensated\nto 0 or greater"
data_long$Var1[data_long$Var1 == "Reduction enhanced"] <- "Reduction\nenhanced"

data_long$total[data_long$variable == "Homozygot"] <- 342
data_long$total[data_long$variable == "Heterozygot"] <- 648

data_long$percentage <- round(data_long$value/data_long$total*100,1)
data_long$percentage <- paste0(data_long$percentage , "%")


g <- ggplot(data_long, aes(x = Var1, y = value, fill = variable)) +
   geom_bar(stat = "identity", position = "dodge") +
   theme(legend.title = element_blank()) +
   labs(title = "", x = "", y = "Frequency")+
   geom_text(aes(x = Var1, y = value, 
                 label = format(percentage, digits = 1)),
             size = 4.5,
             position = position_dodge(.9),
             vjust=-.1)


g
ggsave(g, filename="Supplementary Figure 6 C.png", width = 9, height = 4)


sum(data_long$value[data_long$variable == "Homozygot"])
sum(data_long$value[data_long$variable == "Heterozygot"])

81/342
210/648

chisq.test(outTBLho[, 2:3])
#X-squared = 78.593, df = 5, p-value = 1.652e-15




load("essential Genes")

tsl1_ho_dieter2 <- tsl1_ho[tsl1_ho$gene %in% x,]
tsl1_ho <- ho[ho$TSL1_site == "TSL1" & ho$protein_coding == "Protein coding",]


## Compare essential genes with non-essential

tsl1_all <- tsl1_ho[!tsl1_ho$gene %in% tsl1_ho_dieter2$gene, ]

tsl1_all$type <- "Others"
tsl1_ho_dieter2$type <- "5'ss for breast cancer\nrisk assesment" 

tsl1_ho_dieter2$HBS_diff_group <- NULL

tsl1_all <- rbind(tsl1_all,tsl1_ho_dieter2 )

tsl1_all_unique <- data.frame(SDid=tsl1_all$sdID, type=tsl1_all$type, HBS_diff=tsl1_all$HBS_diff,
                              const = tsl1_all$constitutiveSD )
tsl1_all_unique <- unique(tsl1_all_unique)

tsl1_all_unique <- tsl1_all_unique[ tsl1_all_unique$const == 1,]


# Check percentage of varying SDs
sum(tsl1_ho_dieter2$HBS_diff[tsl1_ho_dieter2$constitutiveSD == 1] == 0)/ nrow(tsl1_ho_dieter2[tsl1_ho_dieter2$constitutiveSD == 1,])
sum(tsl1_all$HBS_diff[tsl1_all$constitutiveSD == 1] == 0)/ nrow(tsl1_all[tsl1_all$constitutiveSD == 1,])

tsl1_all <- tsl1_ho[!tsl1_ho$gene %in% tsl1_ho_dieter2$gene, ]

ks.test(tsl1_ho_dieter2$HBS_diff[tsl1_ho_dieter2$constitutiveSD == 1], tsl1_all$HBS_diff[tsl1_all$constitutiveSD == 1])

sum(tsl1_all_unique$HBS_diff[tsl1_all_unique$type == "Others"] == 0)/ nrow(tsl1_all_unique[tsl1_all_unique$type == "Others",])

sum(tsl1_all_unique$HBS_diff[tsl1_all_unique$type == "5'ss for breast cancer\nrisk assesment"] == 0)/ nrow(tsl1_all_unique[tsl1_all_unique$type == "5'ss for breast cancer\nrisk assesment",])


tsl1_ho <- ho[ho$TSL1_site == "TSL1" & ho$protein_coding == "Protein coding",]


allSDids <- ad[ad$genome ==ad$genome[1],]
allSDids <- data.frame(sdID=allSDids$sdID, gene=allSDids$gene, cano=allSDids$constitutiveSD, coding=allSDids$protein_coding, tsl=allSDids$TSL1_site)
allSDids <- unique(allSDids)

allSDids <- allSDids[allSDids$coding == "Protein coding" & allSDids$tsl == "TSL1",]

length(unique(tsl1_ho_dieter2$sdID[tsl1_ho_dieter2$constitutiveSD == 1]))/length(unique(allSDids$sdID[allSDids$gene %in% x & allSDids$cano == 1]))
length(unique(tsl1_ho$sdID[tsl1_ho$constitutiveSD == 1]))/length(unique(allSDids$sdID[(!allSDids$gene %in% x) & allSDids$cano == 1]))


sum(tsl1_ho_dieter2$HBS_diff == 0)/nrow(tsl1_ho_dieter2)
sum(tsl1_ho$HBS_diff == 0)/nrow(tsl1_ho) 

length(unique(tsl1_ho$sdID))/length(unique(ad$sdID[ad$TSL1_site == "TSL1" & ad$protein_coding == "Protein coding" ]))
sum(tsl1_ho$HBS_diff < -1)/nrow(tsl1_ho) 






#Analyzing homozygous SSHW variations

ho_sshw <- ad
ho_sshw <- ho_sshw[ho_sshw$SSHW_diff ==  ho_sshw$SSHW_diff_h2,]
ho_sshw <- ho_sshw[ho_sshw$SSHW_diff != 0,]

ho_sshw$sdID <- paste(ho_sshw$Chromosome.scaffold.name, ho_sshw$Strand, ho_sshw$SD_cor)


tsl1_ho_sshw <- ho_sshw[ho_sshw$TSL1_site == "TSL1" & ho_sshw$protein_coding == "Protein coding",]

length(unique(ho_sshw$sdID ))

# length(unique(ho_sshw$sdID )) 59316


ho_sshw2 <- data.frame(id=ho_sshw$sdID, tsl1=ho_sshw$TSL1_site, prot=ho_sshw$protein_coding)
ho_sshw2 <- unique(ho_sshw2)
table(ho_sshw2$prot, ho_sshw2$tsl1)

#                     Not TSL1  TSL1
# Not protein coding    37115  5898
# Protein coding        22164  79797

ho_sshw3 <- data.frame(id=c(ho_sshw$sdID[ho_sshw$SSHW_diff != 0],ho_sshw$sdID[ho_sshw$SSHW_diff_h2 != 0]),
                       TSL1_site=c(ho_sshw$TSL1_site[ho_sshw$SSHW_diff != 0],ho_sshw$TSL1_site[ho_sshw$SSHW_diff_h2 != 0]),
                       protein_coding=c(ho_sshw$protein_coding[ho_sshw$SSHW_diff != 0],ho_sshw$protein_coding[ho_sshw$SSHW_diff_h2 != 0]),
                       SSHW_diff= c(ho_sshw$SSHW_diff[ho_sshw$SSHW_diff != 0],ho_sshw$SSHW_diff_h2[ho_sshw$SSHW_diff_h2 != 0]))


library(ggplot2)
g <- ggplot(ho_sshw3,aes(x=protein_coding, y=SSHW_diff, fill=TSL1_site)) + geom_boxplot()+xlab("")+
   ylab("SSHW difference to the reference 5'ss")+ylim(-1500, 1500)
ggsave(g, filename="Figure 5 A.png",width=5, height=5)

table(ho_sshw$protein_coding, ho_sshw$TSL1_site)

hist(tsl1_ho_sshw$SSHW_diff)
qqnorm(tsl1_ho_sshw$SSHW_diff)
qqline(tsl1_ho_sshw$SSHW_diff)





ho_sshw3 <- data.frame(id=c(ho_sshw$sdID[ho_sshw$SSHW_diff != 0],ho_sshw$sdID[ho_sshw$SSHW_diff_h2 != 0]),
                       TSL1_site=c(ho_sshw$TSL1_site[ho_sshw$SSHW_diff != 0],ho_sshw$TSL1_site[ho_sshw$SSHW_diff_h2 != 0]),
                       protein_coding=c(ho_sshw$protein_coding[ho_sshw$SSHW_diff != 0],ho_sshw$protein_coding[ho_sshw$SSHW_diff_h2 != 0]),
                       SSHW_diff= c(ho_sshw$SSHW_diff[ho_sshw$SSHW_diff != 0],ho_sshw$SSHW_diff_h2[ho_sshw$SSHW_diff_h2 != 0]))

tsl1_ho_sshw <- ho_sshw3[ho_sshw3$TSL1_site == "TSL1" & ho_sshw3$protein_coding == "Protein coding",]
tsl1_ho_sshw$SSHW_diff_group <- cut(tsl1_ho_sshw$SSHW_diff, breaks=seq(-250,250,10) )   

tsl1_ho_sshws <- tsl1_ho_sshw[,c(1, 5)]
tsl1_ho_sshws$c <- 1

tsl1_ho_sshws_df <- data.frame(tapply(tsl1_ho_sshws$c, tsl1_ho_sshws$SSHW_diff_group, sum))
tsl1_ho_sshws_df$lab <- row.names(tsl1_ho_sshws_df)

tsl1_ho_sshws <- unique(tsl1_ho_sshws)
tsl1_ho_sshws_df_unqiue <- data.frame(tapply(tsl1_ho_sshws$c, tsl1_ho_sshws$SSHW_diff_group, sum))

tsl1_ho_sshw_df2 <- cbind(tsl1_ho_sshws_df, tsl1_ho_sshws_df_unqiue)
names(tsl1_ho_sshw_df2) <- c("all","lab","unique")

tab1 <- tsl1_ho_sshw_df2[,c(1,2)]
tab1$count <- "Total"

tab2 <- tsl1_ho_sshw_df2[,c(3,2)]
tab2$all <- tab2$unique
tab2$unique <- NULL
tab2 <- tab2[,c(2,1)]
tab2$count <- "Unique 5'ss"

data2 <- rbind(tab1, tab2)

data2$lab <- factor(data2$lab, levels= tab2$lab)

library(ggplot2)
g <- ggplot(data2, aes(fill=count , y= all , x=lab)) + 
   geom_bar(position="identity", stat="identity") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   ylab("Number of 5'ss affected [log10]") + xlab("Group of homozygous SSHW difference to the reference")+
   theme(text = element_text(size = 14))+
   scale_fill_manual(name= "lab", values = c("#A6D854", "#8DA0CB"))+
   scale_color_manual(name = "lab", values = c("#A6D854", "#8DA0CB"))+
   scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
   ) +
   annotation_logticks(sides = 'lr') 




g
ggsave(g, filename="Figure 5 B.png",width=9, height=5)





ho_sshw3 <- data.frame(id=c(ho_sshw$sdID[ho_sshw$SSHW_diff != 0],ho_sshw$sdID[ho_sshw$SSHW_diff_h2 != 0]),
                       TSL1_site=c(ho_sshw$TSL1_site[ho_sshw$SSHW_diff != 0],ho_sshw$TSL1_site[ho_sshw$SSHW_diff_h2 != 0]),
                       protein_coding=c(ho_sshw$protein_coding[ho_sshw$SSHW_diff != 0],ho_sshw$protein_coding[ho_sshw$SSHW_diff_h2 != 0]),
                       SSHW_diff= c(ho_sshw$SSHW_diff[ho_sshw$SSHW_diff != 0],ho_sshw$SSHW_diff_h2[ho_sshw$SSHW_diff_h2 != 0]),
                       HBS_diff= c(ho_sshw$HBS_diff[ho_sshw$SSHW_diff != 0],ho_sshw$HBS_diff_h2[ho_sshw$SSHW_diff_h2 != 0]),
                       ref_HBS=c(ho_sshw$ref_HBS[ho_sshw$SSHW_diff != 0],ho_sshw$ref_HBS[ho_sshw$SSHW_diff_h2 != 0]),
                       ref_SSHW=c(ho_sshw$ref_SSHW[ho_sshw$SSHW_diff != 0],ho_sshw$ref_SSHW[ho_sshw$SSHW_diff_h2 != 0]),
                       alt_SSHW=c(ho_sshw$alt_SSHW[ho_sshw$SSHW_diff != 0],ho_sshw$alt_SSHW[ho_sshw$SSHW_diff_h2 != 0]))

tsl1_ho_sshw <- ho_sshw3[ho_sshw3$TSL1_site == "TSL1" & ho_sshw3$protein_coding == "Protein coding",]
tsl1_ho_sshw$SSHW_diff_group <- cut(tsl1_ho_sshw$SSHW_diff, breaks=seq(-250,250,10) )  


sum(tsl1_ho_sshw$SSHW_diff <= -55 | tsl1_ho_sshw$SSHW_diff >= 55 )/nrow(tsl1_ho_sshw)

tsl1_ho_sshw2 <- tsl1_ho_sshw[tsl1_ho_sshw$SSHW_diff <= -55 | tsl1_ho_sshw$SSHW_diff >= 55,]
length(unique(tsl1_ho_sshw2$id))/length(unique(tsl1_ho_sshw$id))

tsl1_ho2 <- data.frame(coding= tsl1_ho_sshw2$protein_coding, sdID=tsl1_ho_sshw2$id,
                       SSHW_diff=tsl1_ho_sshw2$SSHW_diff, 
                       ref_SSHW=tsl1_ho_sshw2$ref_SSHW,
                       alt_SSHW=tsl1_ho_sshw2$alt_SSHW)


tsl1_ho2$SSHW_diff_group <-  as.numeric(cut_number(tsl1_ho2$SSHW_diff,5)) 
table(tsl1_ho2$SSHW_diff_group )

tsl1_ho2$SSHW_diff_group <-  as.numeric(cut_number(tsl1_ho2$SSHW_diff,5)) 
table(tsl1_ho2$SSHW_diff_group )

test <- tapply(tsl1_ho2$SSHW_diff, tsl1_ho2$SSHW_diff_group, mean)

tsl1_ho2$SSHW_diff_group <- as.factor(tsl1_ho2$SSHW_diff_group)

#SSHW_diff
g <- ggplot(tsl1_ho2,aes(x=SSHW_diff_group, y=SSHW_diff)) + geom_boxplot()+xlab("")+
   ylab("SSHW difference")+xlab("Quantile SSHW difference group")+
   theme(text=element_text(size=18))+ theme(axis.text = element_text(size = 18))
g


ggsave(g, filename="Figure 6 A.png",width=5, height=6)



test <- tapply(tsl1_ho2$ref_SSHW, tsl1_ho2$SSHW_diff_group, mean)

tsl1_ho2$SSHW_diff_group <- as.factor(tsl1_ho2$SSHW_diff_group)



g <- ggplot(tsl1_ho2,aes(x=SSHW_diff_group, y=ref_SSHW)) + geom_boxplot()+xlab("")+
   ylab("SSHW of reference 5'ss")+xlab("Decentile SSHW difference group ")+
   theme(text=element_text(size=18))+ theme(axis.text = element_text(size = 18))+
   scale_y_continuous(limits = c(-700, 2000))

g


ggsave(g, filename="Figure 6 B.png",width=5, height=6)





tsl1_ho2 <- data.frame(coding= tsl1_ho_sshw2$protein_coding, id=tsl1_ho_sshw2$id, 
                       SSHW_diff=tsl1_ho_sshw2$SSHW_diff,
                       HBS_alt=tsl1_ho_sshw2$ref_HBS+tsl1_ho_sshw2$HBS_diff)
                       #HBS_alt=tsl1_ho_sshw2$ref_HBS)

tsl1_ho2 <- na.omit(tsl1_ho2)

tsl1_ho2$SSHW_diff_group <-  as.numeric(cut_number(tsl1_ho2$SSHW_diff,5)) 
table(tsl1_ho2$SSHW_diff_group )

test <- tapply(tsl1_ho2$HBS_alt, tsl1_ho2$SSHW_diff_group, mean)


tsl1_ho2$SSHW_diff_group <- as.factor(tsl1_ho2$SSHW_diff_group)

g <- ggplot(tsl1_ho2,aes(x=SSHW_diff_group, y=HBS_alt)) + geom_boxplot()+xlab("")+
   ylab("HBS of 5'ss")+xlab("Quantile SSHW difference group ")+
   theme(text=element_text(size=18))+ theme(axis.text = element_text(size = 18))
g


ggsave(g, filename="Figure 6 C.png",width=5, height=6)



ho_sshw3 <- data.frame(id=c(ho_sshw$sdID[ho_sshw$SSHW_diff != 0],ho_sshw$sdID[ho_sshw$SSHW_diff_h2 != 0]),
                       TSL1_site=c(ho_sshw$TSL1_site[ho_sshw$SSHW_diff != 0],ho_sshw$TSL1_site[ho_sshw$SSHW_diff_h2 != 0]),
                       protein_coding=c(ho_sshw$protein_coding[ho_sshw$SSHW_diff != 0],ho_sshw$protein_coding[ho_sshw$SSHW_diff_h2 != 0]),
                       SSHW_diff= c(ho_sshw$SSHW_diff[ho_sshw$SSHW_diff != 0],ho_sshw$SSHW_diff_h2[ho_sshw$SSHW_diff_h2 != 0]),
                       HBS_diff= c(ho_sshw$HBS_diff[ho_sshw$SSHW_diff != 0],ho_sshw$HBS_diff_h2[ho_sshw$SSHW_diff_h2 != 0]),
                       ref_HBS=c(ho_sshw$ref_HBS[ho_sshw$SSHW_diff != 0],ho_sshw$ref_HBS[ho_sshw$SSHW_diff_h2 != 0]),
                       ref_SSHW=c(ho_sshw$ref_SSHW[ho_sshw$SSHW_diff != 0],ho_sshw$ref_SSHW[ho_sshw$SSHW_diff_h2 != 0]),
                       alt_SSHW=c(ho_sshw$alt_SSHW[ho_sshw$SSHW_diff != 0],ho_sshw$alt_SSHW[ho_sshw$SSHW_diff_h2 != 0]))


tsl1_ho_sshw <- ho_sshw3[ho_sshw3$TSL1_site == "TSL1" & ho_sshw3$protein_coding == "Protein coding",]
tsl1_ho_sshw$SSHW_diff_group <- cut(tsl1_ho_sshw$SSHW_diff, breaks=seq(-250,250,10) )  


sum(tsl1_ho_sshw$SSHW_diff <= -55 | tsl1_ho_sshw$SSHW_diff >= 55 )/nrow(tsl1_ho_sshw)

tsl1_ho_sshw2 <- tsl1_ho_sshw[tsl1_ho_sshw$SSHW_diff <= -55 | tsl1_ho_sshw$SSHW_diff >= 55,]
length(unique(tsl1_ho_sshw2$id))/length(unique(tsl1_ho_sshw$id))

tsl1_ho2 <- data.frame(coding= tsl1_ho_sshw2$protein_coding, sdID=tsl1_ho_sshw2$id,
                       SSHW_diff=tsl1_ho_sshw2$SSHW_diff, 
                       HBS_diff=tsl1_ho_sshw2$HBS_diff)
tsl1_ho2 <- na.omit(tsl1_ho2)

tsl1_ho2$SSHW_diff_group <-  as.numeric(cut_number(tsl1_ho2$SSHW_diff,5)) 
table(tsl1_ho2$SSHW_diff_group )

test <- tapply(tsl1_ho2$SSHW_diff, tsl1_ho2$SSHW_diff_group, mean)
test

tsl1_ho2$SSHW_diff_group <- as.factor(tsl1_ho2$SSHW_diff_group)

g <- ggplot(tsl1_ho2[tsl1_ho2$HBS_diff != 0,],aes(x=SSHW_diff_group, y=HBS_diff)) + geom_boxplot()+xlab("")+
   ylab("HBS difference")+xlab("Quantile SSHW difference group")+
   theme(text=element_text(size=18))+ theme(axis.text = element_text(size = 18))+
   scale_y_continuous(limits = c(-15, 6))

g


ggsave(g, filename="Figure 6 D.png",width=5, height=6)





ho_sshw3 <- data.frame(id=c(ho_sshw$sdID[ho_sshw$SSHW_diff != 0],ho_sshw$sdID[ho_sshw$SSHW_diff_h2 != 0]),
                       TSL1_site=c(ho_sshw$TSL1_site[ho_sshw$SSHW_diff != 0],ho_sshw$TSL1_site[ho_sshw$SSHW_diff_h2 != 0]),
                       protein_coding=c(ho_sshw$protein_coding[ho_sshw$SSHW_diff != 0],ho_sshw$protein_coding[ho_sshw$SSHW_diff_h2 != 0]),
                       SSHW_diff= c(ho_sshw$SSHW_diff[ho_sshw$SSHW_diff != 0],ho_sshw$SSHW_diff_h2[ho_sshw$SSHW_diff_h2 != 0]),
                       HBS_diff= c(ho_sshw$HBS_diff[ho_sshw$SSHW_diff != 0],ho_sshw$HBS_diff_h2[ho_sshw$SSHW_diff_h2 != 0]),
                       ref_HBS=c(ho_sshw$ref_HBS[ho_sshw$SSHW_diff != 0],ho_sshw$ref_HBS[ho_sshw$SSHW_diff_h2 != 0]),
                       ref_SSHW=c(ho_sshw$ref_SSHW[ho_sshw$SSHW_diff != 0],ho_sshw$ref_SSHW[ho_sshw$SSHW_diff_h2 != 0]),
                       donor_seq=c(ho_sshw$donor_seq[ho_sshw$SSHW_diff != 0],ho_sshw$donor_seq[ho_sshw$SSHW_diff_h2 != 0]),
                       alt_donor_seq=c(ho_sshw$donor_seq_alt[ho_sshw$SSHW_diff != 0],ho_sshw$donor_seq_alt_h2[ho_sshw$SSHW_diff_h2 != 0] ),
                       gene=c(ho_sshw$gene[ho_sshw$SSHW_diff != 0],ho_sshw$gene[ho_sshw$SSHW_diff_h2 != 0] ) )               


# ho_sshw3 <- data.frame(id=ho_sshw$sdID,
#                        TSL1_site=ho_sshw$TSL1_site,
#                        protein_coding=ho_sshw$protein_coding,
#                        SSHW_diff= ho_sshw$SSHW_diff,
#                        HBS_diff= ho_sshw$HBS_diff,
#                        ref_HBS=ho_sshw$ref_HBS,
#                        ref_SSHW=ho_sshw$ref_SSHW,
#                        donor_seq=ho_sshw$donor_seq,
#                        alt_donor_seq=ho_sshw$donor_seq_alt
# )

ho_sshw3 <- ho_sshw3[ho_sshw3$TSL1_site == "TSL1" & ho_sshw3$protein_coding == "Protein coding",]


ho_sshw3$const <- ad$constitutiveSD[match(ho_sshw3$id, ad$sdID)]
ho_sshw3 <- ho_sshw3[ho_sshw3$const == 1,]

tsl1_ho_dieter <- ho_sshw3[ho_sshw3$gene %in% x,]
tsl1_all <- ho_sshw3[!ho_sshw3$gene %in% x,]

tsl1_ho_dieter$endSSHW <- tsl1_ho_dieter$ref_SSHW + tsl1_ho_dieter$SSHW_diff
tsl1_ho_dieter$altHBS <- tsl1_ho_dieter$ref_HBS+tsl1_ho_dieter$HBS_diff
tsl1_all$endSSHW <- tsl1_all$ref_SSHW + tsl1_all$SSHW_diff
tsl1_all$altHBS <- tsl1_all$ref_HBS+tsl1_all$HBS_diff


tsl1_all$type <- "Others"
tsl1_ho_dieter$type <- "5'ss of essential genes" 

tsl1_ho_dieter$HBS_diff_group <- NULL

tsl1_all <- rbind(tsl1_all,tsl1_ho_dieter )


dieterDF <- data.frame(altSSHW=tsl1_all$endSSHW, sdID=tsl1_all$id, altHBS=tsl1_all$altHBS,
                       refHBS=tsl1_all$ref_HBS, diffSSHW=tsl1_all$SSHW_diff, 
                       refSSHW = tsl1_all$ref_SSHW, const=tsl1_all$const, type=tsl1_all$type)

dieterDF$HBS_group <-  cut(dieterDF$altHBS, breaks=seq(1.8,23.8,2) ) 
table(dieterDF$HBS_group )

dieterDF <- dieterDF[dieterDF$const == 1,]

dieterDF$HBS_group <- as.factor(dieterDF$HBS_group)



dieterDF2 <- unique(dieterDF)
dieterDF2 <- na.omit(dieterDF2)


g <- ggplot(dieterDF2,aes(x=HBS_group, y=altSSHW, color=type)) + geom_boxplot(position= position_dodge(preserve = "single"))+xlab("")+
   ylab("Alternative 5'ss SSHW")+xlab("HBond score group")+
   theme(text=element_text(size=18))+ theme(axis.text.x = element_text(size = 18,
                                                                       angle = 90, vjust = 0.5, hjust=1)) +
   theme(legend.title=element_blank())     

g

ggsave(g, filename="Figure 7.png", width=10, height =5)









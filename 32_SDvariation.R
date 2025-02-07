SDvariation <- function(sdCUR){
  
  ## Get basic info
  SDcor <- sdCUR$SD_cor
  chrom <- sdCUR$Chromosome.scaffold.name
  strand <- sdCUR$Strand
  
  ## Convert strand annotation
  if(strand== "+") strand <- "1"
  if(strand== "-") strand <- "-1"
  
  listOfResults <- list()
  
  ## Get SD sequence coordinates range
  donor_cords <- c((SDcor-62) : (SDcor+68))
  if(strand == "-1") donor_cords <- rev(c((SDcor-68) : (SDcor+62)))
  
  ## Search for SNPs within the complete region from -60nt to +60nt around the SD
  rel_snps <- snps[ ((snps$cord1 < (SDcor-62) & snps$cord2 >= (SDcor-62)) | snps$cord2  %in% donor_cords)   &
                      snps$chromosome == chrom ,]
  
  ## if on negative strand, adjust SDcor start accordingly
  if(strand == "-1"){
    
    rel_snps <- snps[ ((snps$cord2 > (SDcor+62) & snps$cord1 <= (SDcor+62)) | snps$cord1  %in% donor_cords)   &
                        snps$chromosome == chrom ,]
    
  }
  
  ## Get Sequence of Donor surrounding with and without the variation
  upborder <- 62
  downborder <- 68
  sequence_range <- getReferenceSequence(chromosome = chrom, 
                                         indexCoordinate = SDcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, hg38)
  
  listOfResults$SDseq <- sequence_range
  
  if(nrow(rel_snps) > 0){ 
    
    ## calculate SNV start-location relative to SD
    rel_snps$distance <- rel_snps$cord1 - SDcor
    if(strand == "-1") rel_snps$distance <-  SDcor - rel_snps$cord2
    rel_snps$location <- "downstream"
    rel_snps$location[rel_snps$distance <= 0] <- "upstream"
    
    ## For every deletion increase SD sequence cords by
    ## number of nt which overlap with variation coordinates
    if(any(rel_snps$type == "DEL")){
      
      dels <- rel_snps[rel_snps$type == "DEL",]
      for(del in 1:nrow(dels)){
        
        #del <- 1
        ## check how many nucleotides are upstream of SD GT (del pos <64)
        varCoos <- c(dels$cord1[del]:(dels$cord2[del]))
        names(varCoos) <-  varCoos-SDcor
        if(strand == "-1") names(varCoos) <-  SDcor - varCoos
        
        n_upSDdel <- varCoos[names(varCoos) < 1]
        n_downSDdel <- varCoos[names(varCoos)> 0]
        
        ## If located upstream from GT dinucleotide, increase upstream border for sequence retrieval
        upborder <- upborder+length(n_upSDdel)
        
        ## If located downstream from GT dinucleotide, increase downstream border for sequence retrieval
        downborder <- downborder+length(n_downSDdel)
        
      }
      
      ## Get the potentially extended surrounding reference sequence
      sequence_range <- getReferenceSequence(chromosome = chrom, 
                                             indexCoordinate = SDcor,
                                             upRange = upborder, downRange = downborder, 
                                             strand = strand, hg38)
      
      ## Get SD sequence coordinates range
      donor_cords <- c((SDcor-upborder) : (SDcor+downborder))
      if(strand == "-1") donor_cords <- rev(c((SDcor-downborder) : (SDcor+upborder)))
    }
    
    ## Split sequence
    sdHBS_SNV <- strsplit(sequence_range, "")[[1]]
    
    ## Insert information of every SNV to generate alternative seq
    for(snv in seq_len(nrow(rel_snps))){
      
      ## Calcualte new HBS
      #snv <- 1
      pos <- which(donor_cords %in% rel_snps$cord1[snv])
      if(length(pos) == 0) next
      varCoos <- c(rel_snps$cord1[snv]:(rel_snps$cord2[snv]))
      varType <-  rel_snps$type[snv]
      altNuc <- rel_snps$alt[snv]
      
      if(strand == "-1"){
        altNuc <- gsub("A","t",altNuc)
        altNuc <- gsub("C","g",altNuc)
        altNuc <- gsub("T","a",altNuc)
        altNuc <- gsub("G","c",altNuc)
        altNuc <- toupper(altNuc)
        altNuc <- reverse(altNuc)
        
      }
      
      if(varType == "SNV" | varType == "DUP" ) sdHBS_SNV[pos] <- altNuc
      #if(varType == "DEL") sdHBS_SNV[pos:(pos+rel_snps$ref_len[snv]-1-1)] <- ""
      #if(varType == "DEL") donor_cords <- donor_cords[-(pos:(pos+rel_snps$ref_len[snv]-1-1))]
      if(varType == "INS") sdHBS_SNV[pos] <- altNuc
      
      dc <- data.frame(donor_cords)
      dc$rep <-  unlist(lapply(sdHBS_SNV, nchar))
      dc$nt <- sdHBS_SNV
      if(varType == "DEL") dc$rep[dc$donor_cords %in% varCoos] <- 0
      if(varType == "DEL") dc$nt[dc$donor_cords %in% varCoos] <- ""
      donor_cords <- unlist(lapply(seq_len(nrow(dc)),
                                   function(dcord){ rep(dc$donor_cords[dcord], 
                                                        dc$rep[dcord]) }))
      
      sdHBS_SNV <- unlist( strsplit(dc$nt, ""))
      
      if(varType == "INS" & 
         rel_snps$location[snv] == "upstream") donor_cords <- donor_cords[nchar(altNuc): length(sdHBS_SNV)]
      if(varType == "INS" & 
         rel_snps$location[snv] == "upstream") sdHBS_SNV <- sdHBS_SNV[nchar(altNuc): length(sdHBS_SNV)]
      
      if(varType == "INS" & 
         rel_snps$location[snv] == "downstream") donor_cords <- donor_cords[1:(length(donor_cords)-nchar(altNuc)+1)]
      if(varType == "INS" & 
         rel_snps$location[snv] == "downstream") sdHBS_SNV <- sdHBS_SNV[1:(length(sdHBS_SNV)-nchar(altNuc)+1)]
      # if(varType == "DUP" & 
      #      rel_snps$location == "upstream") sdHBS_SNV <- sdHBS_SNV[2: nchar(sdHBS_SNV)]
      # if(varType == "DUP" & 
      #      rel_snps$location == "downstream") sdHBS_SNV <- sdHBS_SNV[1:(nchar(sdHBS_SNV)-1)]
      
    }
    
    sdHBS_SNV <- paste(sdHBS_SNV, collapse = "")
    
    listOfResults$SDseq_alt <- sdHBS_SNV
    
    
  }
  
  listOfResults$n_SNV <- sum(rel_snps$type == "SNV")
  listOfResults$n_DEL <- sum(rel_snps$type == "DEL")
  listOfResults$n_INS <- sum(rel_snps$type == "INS")
  
  return(paste(listOfResults, collapse = ","))
}
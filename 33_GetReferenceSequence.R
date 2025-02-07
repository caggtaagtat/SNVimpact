getReferenceSequence <- function(chromosome = "1", indexCoordinate = 1000, upRange = 20,
                                 downRange = 2, strand = "1", dnaStringSet){
  
  ## Convert strand annotation
  if(strand== "+") strand <- "1"
  if(strand== "-") strand <- "-1"
  
  ## Get the respective genomic sequence
  sub.info <- data.frame(chrom=as.character(chromosome), start=indexCoordinate-upRange, end=indexCoordinate+downRange)
  if(strand == "-1") sub.info <- data.frame(chrom=as.character(chromosome), start=indexCoordinate-downRange,
                                            end=indexCoordinate+upRange)
  sequence_range <- as.character(getSeq(dnaStringSet, as(sub.info, "GRanges")))
  
  ## Minus strand sequecne manipulation
  if(strand == "-1")sequence_range <- gsub("A","t",sequence_range)
  if(strand == "-1")sequence_range <- gsub("C","g",sequence_range)
  if(strand == "-1")sequence_range <- gsub("T","a",sequence_range)
  if(strand == "-1")sequence_range <- gsub("G","c",sequence_range)
  if(strand == "-1")sequence_range <- toupper(sequence_range)
  if(strand == "-1")sequence_range <- reverse(sequence_range)
  
  return(sequence_range)
}
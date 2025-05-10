library("ArchR")
library("argparse")
library("tidyverse")
library("data.table")
library("GenomicRanges")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
# library("BSgenome.Mfascicularis.NCBI.5.0")
# library("BSgenome.Cjacchus.UCSC.calJac4")
# library("BSgenome.Mmusculus.UCSC.mm10")
# seqnames(BSgenome.Mfascicularis.NCBI.5.0)=gsub("MFA","chr",seqnames(BSgenome.Mfascicularis.NCBI.5.0))

### Input archr project
proj = readRDS("/cluster/home/huangyiming/cross_species/cross_species_int_2/human/human_1/Save_in_anno/Save-ArchR-Project.rds")
# proj = readRDS("/cluster/home/huangyiming/cross_species/cross_species_int_2/human/human_1/Save_in_anno/Save-ArchR-Project.rds")
# proj = readRDS("/cluster/home/huangyiming/cross_species/cross_species_int_2/human/human_1/Save_in_anno/Save-ArchR-Project.rds")
# proj = readRDS("/cluster/home/huangyiming/cross_species/cross_species_int_2/human/human_1/Save_in_anno/Save-ArchR-Project.rds")


### input parameter
args = commandArgs(T)
type_id = args[1]
chromLengths <- getChromLengths(proj)
chromLengths <- chromLengths [names(chromLengths) %ni% "chrY"]
blacklist = getBlacklist(proj)


# read bed to gr
read2gr <- function(bedF, label){
  df <- fread(bedF, sep="\t", header=F)
  colnames(df) <- c("chr", "start", "end", "name", "score")
  df$label <- label
  gr <- GRanges(
    df$chr,
    IRanges(df$start, df$end)
  );
  mcols(gr)$score <- df$score
  mcols(gr)$name <- df$name
  mcols(gr)$label <- df$label
  return(gr)
}

# extend summit to 500 bp
extendSummit <- function(gr, size=500){
  gr <- resize(gr, width=size, fix="center")
  return(gr)
}


# filter blacklist
filter4blacklist <- function(gr,  black_list.gr = blacklist){
  black_list = as.data.frame(blacklist);
  black_list.gr = GRanges(
    black_list[,1],
    IRanges(black_list[,2], black_list[,3])
  );
  
  idx = queryHits(
    findOverlaps(gr, black_list.gr)
  );
  if(length(idx) > 0){
    gr = gr[-idx]
  }
  return(gr)
}

# filter non-chromosome
filter4chrom <- function(gr, chromL=chromLengths) {
  chrom.gr <- GRanges(
    names(chromLengths),
    IRanges(0, chromLengths)
  )
  idx <- queryHits(
    findOverlaps(gr, chrom.gr, type="within")
  )
  if(length(idx) > 0) {
    gr <- gr[idx]
  }
  return(gr)
}


# filter N containing regions
filter4N <- function(gr, genome=BSgenome.Hsapiens.UCSC.hg38){
  genome <- getBSgenome(genome)
  nucFreq <- BSgenome::alphabetFrequency(getSeq(genome, gr))
  mcols(gr)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  mcols(gr)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  gr[which(mcols(gr)$N < 0.001)] #Remove N Containing Peaks
  return(gr)
}

nonOverlappingGR <- function(
    gr = NULL, 
    by = "score", 
    decreasing = TRUE, 
    verbose = FALSE
) {
  stopifnot(by %in% colnames(mcols(gr)))
  #-----------
  # Cluster GRanges into islands using reduce and then select based on input
  #-----------
  .clusterGRanges <- function(gr = NULL, filter = TRUE, by = "score", decreasing = TRUE) {
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r, ignore.strand = TRUE)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  
  if(verbose) { message("Converging", appendLF = FALSE) }
  i <-  0
  grConverge <- gr
  while(length(grConverge) > 0) {
    if(verbose){ message(".", appendLF = FALSE) }
    i <-  i + 1
    grSelect <- .clusterGRanges(
      gr = grConverge, 
      filter = TRUE, 
      by = by, 
      decreasing = decreasing)
    
    grConverge <- subsetByOverlaps(
      grConverge,
      grSelect, 
      invert=TRUE, 
      ignore.strand = TRUE) #blacklist selected gr
    
    if(i == 1){ #if i=1 then set gr_all to clustered
      grAll <- grSelect
      
    }else{
      grAll <- c(grAll, grSelect)
    } 
  }
  message(sprintf("Converged after %s iterations!", i))
  
  if(verbose){
    message("\nSelected ", length(grAll), " from ", length(gr))
  }
  grAll <- sort(sortSeqlevels(grAll))
  return(grAll)
}

# normlize to score per million
norm2spm <- function(gr, by = "score") {
  mlogp <- mcols(gr)[,by]
  normmlogp <- 10^6 * mlogp / sum(mlogp)
  mcols(gr)$spm <- normmlogp
  return(gr)
}


for (i in 1:length(a)) {
  type_id = a[i]
  print(type_id)
  path = paste0(type_id,"/",type_id,"_0.05","/",type_id,"_summits.bed")
  p.gr <- read2gr(path, label = type_id)
  p.gr <- extendSummit(p.gr, size = 501)
  p.gr <- filter4blacklist(p.gr, black_list.gr=blacklist)
  p.gr <- filter4chrom(p.gr)
  p.gr <- filter4N(p.gr)
  p.gr <- nonOverlappingGR(p.gr, by = "score", decreasing = TRUE)
  p.gr <- norm2spm(p.gr, by = "score")
  outPeak <- as.data.frame(p.gr)
  write.table(outPeak,paste0(type_id,"/",type_id,"_filter_peak.bed"),sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}


####  Conserve peak
type = read.table("cell_type")
type = type$V1
result = NULL
for (i in 1:length(type)) {
  data = fread(paste0(type[i],"/",type[i],"_filter_peak.bed"))
  result[[i]]  = data
}

total = do.call(rbind,result)
gr <- GRanges(
  total$V1,
  IRanges(total$V2, total$V3)
);
mcols(gr)$score <- total$V11
mcols(gr)$name <- total$V7
p.gr <- nonOverlappingGR(gr, by = "score", decreasing = TRUE)

write.table(p.gr,"consensus_regions2.bed")
library(ArchR)
library(data.table)

## HUMAN
data = readRDS("/cluster/home/huangyiming/public_data/human_7/Save-ProjHeme1_anno/Save-ArchR-Project.rds")

## MACACA
data =readRDS("/cluster/home/huangyiming/macaca_PFC_2/Save-ProjHeme2_anno/Save-ArchR-Project.rds")

## MARMOSET
data = readRDS("~/marmoset/marmoset_7/Save-ProjHeme2_anno/Save-ArchR-Project.rds")

## MOUSE
data =readRDS("~/public_data/mouse_3/fragment/Mouse_FPC_RemoveDoublet_anno/Save-ArchR-Project.rds")

### Read consuse_peaks
peak = fread("consensus_regions2.bed")

a <- GRanges(
  seqnames = peak$V1,
  ranges = IRanges(start = peak$V2, end = peak$V3),
  symbol = peak$V7,
  score = peak$V6
)


projHeme2 = addPeakSet(
  ArchRProj = data,
  peakSet = a,
  genomeAnnotation = getGenomeAnnotation(data),
  force =T
)

## Add peak matrix
projHeme4 = addPeakMatrix(
  ArchRProj = projHeme2,
  ceiling = 4,
  binarize = FALSE,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addPeakMatrix")
)

## Get peakmatrix
data = getMatrixFromProject(
  ArchRProj = projHeme4,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)


## Calculate peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = data, 
  useMatrix = "PeakMatrix", 
  groupBy = "anno",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
markerList

lapply(markerList,function(x) {print(x)})
for (i in 1:length(names(markerList))) {
  markerList[[i]]$group = names(markerList)[i]
}

total_DEG_peak = do.call(rbind,markerList)
total_DEG_peak$ID = paste(total_DEG_peak$seqnames,":",total_DEG_peak$start,"-",total_DEG_peak$end,sep = "")

saveRDS(total_DEG_peak,"Marker_DAG.rds")





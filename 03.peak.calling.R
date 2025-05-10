library(tidyverse)
library(ArchR)


args = commandArgs(T)
type_id = args[1]
print(Sys.time())
print(type_id)
print("start")

projHeme2 = readRDS("ArchR_project/Save-ArchR-Project.rds")
ID_check_mouse = projHeme2@cellColData
print(dim(ID_check_mouse))

bed = getFragmentsFromProject(
  ArchRProj = projHeme2,
  cellNames = row.names(ID_check_mouse[ID_check_mouse$anno %in% type_id,]),
  verbose = T,
  logFile = createLogFile("getFragmentsFromProject")
)

bed <- do.call(c, bed)
grList_new_strand <- lapply(bed, function(gr) {
  strand(gr) <- "+"
  return(gr)
})

grList_new_strand <- GRangesList(grList_new_strand)
all_beds <- unlist(grList_new_strand)
exchr = getBlacklist(projHeme2)
all_beds <- subsetByOverlaps(all_beds, exchr, invert=TRUE)
dir.create(type_id)

rtracklayer::export(all_beds,con=paste0(type_id,"/",type_id,"_blacklistrm.bed"), format="bed")
all_beds = as_tibble(all_beds) %>% select(seqnames,start,end,RG,strand)
fwrite(all_beds,paste0(type_id,"/",type_id,"_matrix.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Input path
input = paste0(type_id,"/",type_id,"_blacklistrm.bed")
outname = type_id
outdir = paste0(type_id,"/",type_id,"_0.05")

## 注意基因组长度
## human:3.2e9/macaca:2.95e9/marmoset:2.90e9/mouse:2.73e9
print("call peak")
system(paste0("macs2 callpeak -t ", input,
              " -f BED -n ", outname,
              " --outdir ", outdir,
              " -g 2.95e9 -q 0.05 -B --SPMR",
              " --nomodel --keep-dup all --call-summits --nolambda --shift 75 --extsize 150"))

print(Sys.time())
print("end")

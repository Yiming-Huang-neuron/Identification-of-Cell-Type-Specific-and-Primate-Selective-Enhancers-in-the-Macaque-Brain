#########Genome browser tracks plot
library(Gviz)
library(ggplotify)
library(patchwork)



track = function(ymax) {
  list(
    bw_L2 <- DataTrack(range = "L2.bw",genome="hg38",col.histogram = c("#f0e6c8"),background.title = "white",fill.histogram = c("#f0e6c8"),type="histogram",  ylim = c(0,ymax)),
    bw_L23 <- DataTrack(range = "L23.bw",genome="hg38",col.histogram = c("#f9ed69"),background.title = "white",fill.histogram = c("#f9ed69"),type="histogram",  ylim = c(0,ymax)),
    bw_L34 <- DataTrack(range = "L34.bw",genome="hg38",col.histogram = c("#fce825"),background.title = "white",fill.histogram = c("#fce825"),type="histogram",   ylim = c(0,ymax)),
    bw_L4 <- DataTrack(range = "L4.bw",genome="hg38",col.histogram = c("#fed26f"),background.title = "white",fill.histogram = c("#fed26f"),type="histogram",ylim = c(0,ymax)),
    bw_L45 <- DataTrack(range = "L45.bw",genome="hg38",col.histogram = c("#dfc621"),background.title = "white",fill.histogram = c("#dfc621"),type="histogram", ylim = c(0,ymax)),
    bw_L456 <- DataTrack(range = "L456.bw",genome="hg38",col.histogram = c("#f8ab5d"),background.title = "white",fill.histogram = c("#f8ab5d"),type="histogram",  ylim = c(0,ymax)),
    bw_L56 <- DataTrack(range = "L56.bw",genome="hg38",col.histogram = c("#fb8a2e"),background.title = "white",fill.histogram = c("#fb8a2e"),type="histogram",  ylim = c(0,ymax)),
    bw_L6 <- DataTrack(range = "L6.bw",genome="hg38",col.histogram = c("#f28863"),background.title = "white",fill.histogram = c("#f28863"),type="histogram",  ylim = c(0,ymax)),
    bw_PV_CHC <- DataTrack(range = "PV_CHC.bw",genome="hg38",col.histogram = c("#96cd8c"),background.title = "white",fill.histogram = c("#96cd8c"),type="histogram",  ylim = c(0,ymax)),
    bw_PVALB <- DataTrack(range = "PVALB.bw",genome="hg38",col.histogram = c("#5bbc87"),background.title = "white",fill.histogram = c("#5bbc87"),type="histogram",  ylim = c(0,ymax)),
    bw_SST <- DataTrack(range = "SST.bw",genome="hg38",col.histogram = c("#bd9ac6"),background.title = "white",fill.histogram = c("#bd9ac6"),type="histogram",  ylim = c(0,ymax)),
    bw_VIP <- DataTrack(range = "VIP.bw",genome="hg38",col.histogram = c("#6dc5dc"),background.title = "white",fill.histogram = c("#6dc5dc"),type="histogram",  ylim = c(0,ymax)),
    bw_RELN <- DataTrack(range = "VIP_RELN.bw",genome="hg38",col.histogram = c("#3b9bcc"),background.title = "white",fill.histogram = c("#3b9bcc"),type="histogram",  ylim = c(0,ymax)),
    bw_RELN <- DataTrack(range = "RELN.bw",genome="hg38",col.histogram = c("#79b9e7"),background.title = "white",fill.histogram = c("#79b9e7"),type="histogram",  ylim = c(0,ymax)),
    bw_LAMP5 <- DataTrack(range = "LAMP5.bw",genome="hg38",col.histogram = c("#bc9b6a"),background.title = "white",fill.histogram = c("#bc9b6a"),type="histogram",  ylim = c(0,ymax)),
    bw_ASC <- DataTrack(range = "ASC.bw",genome="hg38",col.histogram = c("#edb9d5"),background.title = "white",fill.histogram = c("#edb9d5"),type="histogram",  ylim = c(0,ymax)),
    bw_MG <- DataTrack(range = "MG.bw",genome="hg38",col.histogram = c("#a2ceed"),background.title = "white",fill.histogram = c("#a2ceed"),type="histogram",  ylim = c(0,ymax)),
    bw_OLG <- DataTrack(range = "OLG.bw",genome="hg38",col.histogram = c("#96588a"),background.title = "white",fill.histogram = c("#96588a"),type="histogram",  ylim = c(0,ymax)),
    bw_OPC <- DataTrack(range = "OPC.bw",genome="hg38",col.histogram = c("#939393"),background.title = "white",fill.histogram = c("#939393"),type="histogram",  ylim = c(0,ymax))
  )
}


file_location = read.table("Positive_13_enhancer.bed")
result = NULL
for (i in 1:length(file_location$V1)) {
  a = track(100)
  result[[i]] = as.ggplot(grid2grob(plotTracks(c(a),from = file_location[i,2]-50, to = file_location[i,3]+50,chromosome = file_location[i,1],add = T,type="hist",extend.right=0,extend.left=0,widow = 20,showTitle = F)))
}
plot_layout <- wrap_plots(result,nrow=1)
plot_layout
ggsave(paste0("PEAKS/Genome_browser_tracks.pdf"),plot_layout,width = 10,height = 10)

#########Re-annotaion genome browser tracks plot
library(Gviz)
library(ggplotify)
library(patchwork)

track = function(ymax) {
  list(
    bw_L2IT <- DataTrack(range = "L2.bw",genome="hg38",col.histogram = c("#f0e6c8"),background.title = "white",filsl.histogram = c("#f0e6c8"),type ="histogram",  ylim = c(0,ymax)),
    bw_L23IT <- DataTrack(range = "L23.bw",genome="hg38",col.histogram = c("#f9ed69"),background.title = "white",filsl.histogram = c("#f9ed69"),type="histogram",  ylim = c(0,ymax)),
    bw_L34IT <- DataTrack(range = "L34.bw",genome="hg38",col.histogram = c("#fce825"),background.title = "white",fill.histogram = c("#fce825"),type="histogram",   ylim = c(0,ymax)),
    bw_L4IT <- DataTrack(range = "L4.bw",genome="hg38",col.histogram = c("#fed26f"),background.title = "white",fill.histogram = c("#fed26f"),type="histogram",ylim = c(0,ymax)),
    bw_L45IT <- DataTrack(range = "L45.bw",genome="hg38",col.histogram = c("#dfc621"),background.title = "white",fill.histogram = c("#dfc621"),type="histogram", ylim = c(0,ymax)),
    bw_L5CT <- DataTrack(range = "L5.bw",genome="hg38",col.histogram = c("#f8ab5d"),background.title = "white",fill.histogram = c("#f8ab5d"),type="histogram",  ylim = c(0,ymax)),
    bw_L56NP <- DataTrack(range = "L56.bw",genome="hg38",col.histogram = c("#fb8a2e"),background.title = "white",fill.histogram = c("#fb8a2e"),type="histogram",  ylim = c(0,ymax)),
    bw_L6CAR3 <- DataTrack(range = "L6CAR3.bw",genome="hg38",col.histogram = c("#F47D2B"),background.title = "white",fill.histogram = c("#F47D2B"),type="histogram",  ylim = c(0,ymax)),
    bw_L6IT <- DataTrack(range = "L6_1.bw",genome="hg38",col.histogram = c("#FF5744"),background.title = "white",fill.histogram = c("#FF5744"),type="histogram",  ylim = c(0,ymax)),
    bw_L6CT <- DataTrack(range = "L6_2.bw",genome="hg38",col.histogram = c("#C70039"),background.title = "white",fill.histogram = c("#C70039"),type="histogram",  ylim = c(0,ymax)),
    bw_L6B <- DataTrack(range = "L6B.bw",genome="hg38",col.histogram = c("#900C3F"),background.title = "white",fill.histogram = c("#900C3F"),type="histogram",  ylim = c(0,ymax)),
    bw_LAMP5 <- DataTrack(range = "LAMP5.bw",genome="hg38",col.histogram = c("#bc9b6a"),background.title = "white",fill.histogram = c("#bc9b6a"),type="histogram",  ylim = c(0,ymax)),
    bw_LAMP5_LHX6 <- DataTrack(range = "LAMP5_LHX6.bw",genome="hg38",col.histogram = c("#bc9b1a"),background.title = "white",fill.histogram = c("#bc9b1a"),type="histogram",  ylim = c(0,ymax)),
    bw_LAMP5_RELN <- DataTrack(range = "LAMP5_RELN.bw",genome="hg38",col.histogram = c("#aea400"),background.title = "white",fill.histogram = c("#aea400"),type="histogram",  ylim = c(0,ymax)),
    bw_RELN <- DataTrack(range = "RELN.bw",genome="hg38",col.histogram = c("#79b9e7"),background.title = "white",fill.histogram = c("#79b9e7"),type="histogram",  ylim = c(0,ymax)),
    bw_PAX6 <- DataTrack(range = "PAX6.bw",genome="hg38",col.histogram = c("#90b1e7"),background.title = "white",fill.histogram = c("#90b1e7"),type="histogram",  ylim = c(0,ymax)),
    bw_VIP_RELN_1 <- DataTrack(range = "VIP_RELN_1.bw",genome="hg38",col.histogram = c("#74d2e7"),background.title = "white",fill.histogram = c("#74d2e7"),type="histogram",  ylim = c(0,ymax)),
    bw_VIP_RELN_2 <- DataTrack(range = "VIP_RELN_2.bw",genome="hg38",col.histogram = c("#48a9c5"),background.title = "white",fill.histogram = c("#48a9c5"),type="histogram",  ylim = c(0,ymax)),
    bw_VIP_RELN_3 <- DataTrack(range = "VIP_RELN_3.bw",genome="hg38",col.histogram = c("#0085ad"),background.title = "white",fill.histogram = c("#0085ad"),type="histogram",  ylim = c(0,ymax)),
    bw_VIP <- DataTrack(range = "VIP.bw", genome="hg38",col.histogram = c("#005670"),background.title = "white",fill.histogram = c("#005670"),type="histogram",  ylim = c(0,ymax)),
    bw_PV_CHC <- DataTrack(range = "PV_CHC.bw",genome="hg38",col.histogram = c("#96cd8c"),background.title = "white",fill.histogram = c("#96cd8c"),type="histogram",  ylim = c(0,ymax)),
    bw_PVALB_1 <- DataTrack(range = "PVALB_1.bw",genome="hg38",col.histogram = c("#5bbc17"),background.title = "white",fill.histogram = c("#5bbc17"),type="histogram",  ylim = c(0,ymax)),
    bw_PVALB_2 <- DataTrack(range = "PVALB_2.bw",genome="hg38",col.histogram = c("#5bbc87"),background.title = "white",fill.histogram = c("#5bbc87"),type="histogram",  ylim = c(0,ymax)),
    bw_PVALB_3 <- DataTrack(range = "PVALB_3.bw",genome="hg38",col.histogram = c("#5bbc97"),background.title = "white",fill.histogram = c("#5bbc97"),type="histogram",  ylim = c(0,ymax)),
    bw_SST_L23 <- DataTrack(range = "SST_L23.bw",genome="hg38",col.histogram = c("#bd5ac6"),background.title = "white",fill.histogram = c("#bd5ac6"),type="histogram",  ylim = c(0,ymax)),
    bw_SST_L23_27 <- DataTrack(range = "SST_L3.bw",genome="hg38",col.histogram = c("#bd6ac6"),background.title = "white",fill.histogram = c("#bd6ac6"),type="histogram",  ylim = c(0,ymax)),
    bw_SST_L35 <- DataTrack(range = "SST_L35.bw",genome="hg38",col.histogram = c("#bd7ac6"),background.title = "white",fill.histogram = c("#bd7ac6"),type="histogram",  ylim = c(0,ymax)),
    bw_SST_L56 <- DataTrack(range = "SST_L56.bw",genome="hg38",col.histogram = c("#bd8ac6"),background.title = "white",fill.histogram = c("#bd8ac6"),type="histogram",  ylim = c(0,ymax)),
    bw_ASC <- DataTrack(range = "ASC.bw",genome="hg38",col.histogram = c("#edb9d5"),background.title = "white",fill.histogram = c("#edb9d5"),type="histogram",  ylim = c(0,ymax)),
    bw_MG <- DataTrack(range = "MG.bw",genome="hg38",col.histogram = c("#a2ceed"),background.title = "white",fill.histogram = c("#a2ceed"),type="histogram",  ylim = c(0,ymax)),
    bw_OLG <- DataTrack(range = "OLG.bw",genome="hg38",col.histogram = c("#96588a"),background.title = "white",fill.histogram = c("#96588a"),type="histogram",  ylim = c(0,ymax)),
    bw_OPC <- DataTrack(range = "OPC.bw",genome="hg38",col.histogram = c("#939393"),background.title = "white",fill.histogram = c("#939393"),type="histogram",  ylim = c(0,ymax))
  )
}

file_location = read.table("Positive_13_enhancer.bed")
result = NULL
for (i in 1:length(file_location$V1)) {
  a = track(100)
  result[[i]] = as.ggplot(grid2grob(plotTracks(c(a),from = file_location[i,2]-50, to = file_location[i,3]+50,chromosome = file_location[i,1],add = T,type="hist",extend.right=0,extend.left=0,widow = 20,showTitle = F)))
}
plot_layout <- wrap_plots(result,nrow=1)
plot_layout
ggsave(paste0("PEAKS/Genome_browser_tracks.pdf"),plot_layout,width = 5,height = 10)

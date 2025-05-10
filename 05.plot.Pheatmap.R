library(tidyverse)
library(pheatmap)
library(ArchR)

data= readRDS("cor_RNA_ATAC.rds")

atac = c("L2.scATAC","L23.scATAC","L34.scATAC","L4.scATAC","L45.scATAC","L456.scATAC","L56.scATAC","L6.scATAC",
         "PV.CHC.scATAC","PVALB.scATAC","SST.scATAC","VIP.scATAC","VIP.RELN.scATAC","RELN.scATAC","LAMP5.scATAC",
         "ASC.scATAC","MG.scATAC","OLG.scATAC","OPC.scATAC")
atac = data.frame(group = atac)
row.names(atac) = atac$group

rna = c("L2.scRNA","L2.3.scRNA","L3.4.scRNA","L4.scRNA","L4.5.scRNA","L4.5.6.scRNA","L5.6.scRNA","L6.scRNA",
        "PV.CHC.scRNA","PVALB.scRNA","SST.scRNA","VIP.scRNA","VIP.RELN.scRNA","RELN.scRNA","LAMP5.scRNA",
        "ASC.scRNA","MG.scRNA","OLG.scRNA","OPC.scRNA")
rna = data.frame(group_1 = rna)
row.names(rna) = rna$group_1

color_scATAC = c("L2.scATAC" = "#f0e6c8","L23.scATAC" = "#f9ed69","L34.scATAC" ="#fce825","L4.scATAC" ="#fed26f",
                 "L45.scATAC" = "#dfc621","L456.scATAC" ="#f8ab5d","L56.scATAC" = "#fb8a2e","L6.scATAC" ="#f28863",
                 "ASC.scATAC" = "#edb9d5","OLG.scATAC" ="#96588a","OPC.scATAC" ="#939393","MG.scATAC" ="#a2ceed",
                 "PV.CHC.scATAC" = "#96cd8c","PVALB.scATAC" ="#5bbc87","RELN.scATAC" = "#79b9e7","SST.scATAC" = "#bd9ac6",
                 "VIP.scATAC" = "#6dc5dc","VIP.RELN.scATAC" ="#3b9bcc","LAMP5.scATAC" ="#bc9b6a")


color_scRNA = c("L2.scRNA" = "#ec6400","L2.3.scRNA" = "#ea4335","L3.4.scRNA" = "#d85697","L4.scRNA" = "#b6015c",
                "L4.5.scRNA" = "#ae2900","L4.5.6.scRNA" ="#831e38","L5.6.scRNA" = "#8f1010","L6.scRNA" ="#cd694a",
                "ASC.scRNA" = "#ad5f7d","OLG.scRNA" = "#7e29d3","OPC.scRNA" = "#818499","MG.scRNA" ="#2b86ee",
                "PV.CHC.scRNA" = "#009776","PVALB.scRNA" = "#027f44","RELN.scRNA" ="#44712e","SST.scRNA" ="#631d76",
                "VIP.scRNA" = "#005c3a","LAMP5.scRNA" = "#744f28","VIP.RELN.scRNA" ="#226d00")

bk <- c(seq(0,0.5,by=0.01),seq(0.51,1,by=0.01))
pheatmap(data[rna$group_1,atac$group],
         cluster_rows = F,cluster_cols = F,
         annotation_colors = list(group = color_scATAC,group_1 = color_scRNA),
         annotation_col = atac,annotation_row = rna,legend_breaks=seq(0,1,0.5),breaks=bk,
         color = c(colorRampPalette(colors = c("navy","#8db9ca","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#ffdc80","#af0809"))(length(bk)/2)))
library(tidyverse)


### Umap plot
color_scATAC = c("L2-scATAC" = "#f0e6c8","L2/3-scATAC" = "#f9ed69","L3/4-scATAC" ="#fce825","L4-scATAC" ="#fed26f",
                 "L4/5-scATAC" = "#dfc621","L4/5/6-scATAC" ="#f8ab5d","L56-scATAC" = "#fb8a2e","L6-scATAC" ="#f28863",
                 "ASC-scATAC" = "#edb9d5","OLG-scATAC" ="#96588a","OPC-scATAC" ="#939393","MG-scATAC" ="#a2ceed",
                 "PV_CHC-scATAC" = "#96cd8c","PVALB-scATAC" ="#5bbc87","RELN-scATAC" = "#79b9e7","SST-scATAC" = "#bd9ac6",
                 "VIP-scATAC" = "#6dc5dc","VIP_RELN-scATAC" ="#3b9bcc","LAMP5-scATAC" ="#bc9b6a")

data= readRDS("umap_plot.rds")
names(data) = c("UMAP_1","UMAP_2","int_2")
cell_type_med <- data %>%
  group_by(int_2) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

##########去掉outlier
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  pos1=which( x < (qnt[1] - H))
  pos2=which( x > (qnt[2] + H))
  c(pos1,pos2)
}


data_split = split(data,data$int_2)
sapply(data_split,nrow)
re=data.frame()
outlier_number=c()

for(i in 1:length(data_split)){
  # i = 1
  print(names(data_split[i]))
  a=data_split[[i]]
  median_X=median(a[,"UMAP_1"])
  median_Y=median(a[,"UMAP_2"])
  distance_original_point <- apply(a,1,function(x) 
    (as.numeric(x["UMAP_1"])-median_X)**2+(as.numeric(x["UMAP_2"])-median_Y)**2 )
  a$distance_original_point=distance_original_point
  a$distance_original_point1=sqrt(distance_original_point)
  pos=remove_outliers(a$distance_original_point1)
  length(pos)
  a[pos,"Sub"]="outlier"
  re=rbind(re,a)
  outlier_number=c(outlier_number,length(pos))
}

table(re$Sub)
re = re[!re$Sub %in% "outlier",]
head(re)

p = ggplot(re,aes(UMAP_1,UMAP_2,color = int_2)) + 
  geom_point(size = 0.1) + 
  scale_color_manual(values = c(color_scRNA,color_scATAC)) +
  theme_classic() + facet_wrap(~int) + 
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white")) 
p
ggsave("umap.pdf",p)

##### ArchR marker gene plot
### Get marker genes matrix
GS_matrix_markergenes = ArchR:::.getMatrixValues(ArchRProj = data, 
                                         name = c("SLC17A7","GAD1","CX3CR1","PDGFRA","PLP1","SLC1A2"), 
                                         matrixName = "GeneScoreMatrix", 
                                         log2Norm = T, 
                                         threads = getArchRThreads(),
                                         logFile = NULL)

gene_expression = as.matrix(t(GS_matrix_markergenes)) %>% as.data.frame() %>% rownames_to_column(var = "ID") %>% as_tibble 
umapharmoney = getEmbedding(ArchRProj = data, embedding = "UMAPHarmony", returnDF = TRUE)
marker_gene_point = cbind(gene_expression,umapharmoney)

## plot
for (i in 1:length(marker_gene)) {
  gene_ID = marker_gene[i]
  p = ggplot(total_plot,aes(UMAP_1,UMAP_2 )) + 
    geom_point(aes_string(color = gene_ID),size = 0.1) + 
    scale_color_gradientn(colours = c("#e4e0ee","#8f5db7","#8f499c","#59267c")) + 
    theme_classic() +  theme(
      legend.position = "none",
      axis.line = element_blank(),        
      axis.text.x = element_blank(),      
      axis.text.y = element_blank(),       
      axis.ticks = element_blank(),      
      axis.title.x = element_blank(),      
      axis.title.y = element_blank()      
    )
  p
  plot_slice_headline = total_plot[total_plot[gene_ID] > 2.5,]
  p3 = p+geom_point(data = plot_slice_headline,aes(UMAP_1,UMAP_2),size = 0.1,color = "#59267c") + 
    xlab("UMAP_1") + ylab("UMAP_2") +  theme(
      legend.position = "none",
      axis.line = element_blank(),        
      axis.text.x = element_blank(),      
      axis.text.y = element_blank(),       
      axis.ticks = element_blank(),      
      axis.title.x = element_blank(),      
      axis.title.y = element_blank()      
    )
  # p3
  p3
  ggsave(paste(gene_ID,"_m1.tiff",sep = ""),p3,width = 10,height = 10,dpi = 600,compression = "lzw",bg = "white")
}

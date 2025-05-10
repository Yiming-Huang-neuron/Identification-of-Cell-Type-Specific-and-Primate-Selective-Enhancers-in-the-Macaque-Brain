library(tidyverse)

data = read.table("annotation_meta.txt")
adjustedRandIndex(data$old_anno,data$new_anno)
meta6 = table(data$old_anno,data$new_anno) %>% melt()
head(meta6)

color_scRNA = c("L2" = "#f0e6c8","L23" = "#f9ed69","L34" ="#fce825","L4" ="#fed26f",
                "L45" = "#dfc621","L5" ="#f8ab5d","L56" = "#fb8a2e","L6CAR3" ="#F47D2B","L6_1" ="#FF5744",
                "L6_2" ="#C70039","L6B" ="#900C3F",
                "LAMP5_LHX6" ="#bc9b1a","LAMP5" ="#bc9b6a","LAMP5_RELN" ="#aea400","PAX6" = "#90b1e7","RELN" = "#79b9e7",
                "VIP_RELN_1" ="#74d2e7","VIP_RELN_2" ="#48a9c5","VIP_RELN_3" ="#0085ad","VIP" = "#005670",
                "PVALB_1" ="#5bbc17","PVALB_2" ="#5bbc87","PVALB_3" ="#5bbc97",
                "PV_CHC" = "#96cd8c",
                "SST_L23" = "#bd5ac6","SST_L3" = "#bd6ac6","SST_L35" = "#bd7ac6","SST_L56" = "#bd8ac6","SST_CHODL" = "#bd6ac6",
                "ASC" = "#edb9d5","MG" ="#a2ceed","OLG" ="#96588a","OPC" ="#939393")

### 
type1 = c("OPC","OLG","MG","ASC",
          "PAX6","SST_CHODL","VIP","SST","PVALB","PV_CHC","LAMP5","RELN",
          "L6CAR3","L6B","L6_2","L6_1","L6","L56","L5","L45","L4","L34","L23")

type2 = c(rev(names(color_scRNA)))

meta6$Var1 = factor(meta6$Var1,levels = type1)
meta6$Var2 = factor(meta6$Var2,levels = type2)

## filter cells
meta6 = meta6[meta6$value>10,]
table(meta6$Var1)

p = ggplot(meta6,
           aes(y =value,
               axis1 = Var1, axis2 = Var2))+
  theme_bw() +  
  geom_alluvium(aes(fill = Var2),width = 0.7, reverse = FALSE,discern = TRUE)+
  geom_stratum(width = 1/3, reverse = FALSE,discern = TRUE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),reverse = FALSE, size = 4,angle=0,discern = TRUE)+ 
  scale_x_continuous(breaks = 1:2, labels = c("human", "integration"))+
  theme(legend.position = "none") +
  scale_fill_manual(values = color_scRNA) 
p


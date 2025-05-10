library(tidyverse)
library(RColorBrewer)

boxPlot_color =c("CTCF-only,CTCF-bound"="#fbb040",
                 dELS  = "#4078c0","dELS,CTCF-bound" = "#00457c" ,
                 "DNase-H3K4me3" = "#d25238","DNase-H3K4me3,CTCF-bound" = "#ad182d",
                 "pELS" = "#8ba753","pELS,CTCF-bound" = "#30660f",         
                 "PLS"  = "#8a8acb" ,"PLS,CTCF-bound" ="#631d76",
                 "unconserve_region" = "#caccd1",
                 "CR" = "#629aa9")

#### human FigureD
data= readRDS("macaca_liftover_human_CRE.rds")
data = data[!grepl("PLS|PLS,CTCF-bound",data$Var1),]
data$Var1 = factor(data$Var1,levels = c("unconserve_region",
                                        "CR",
                                        "CTCF-only,CTCF-bound",
                                        "pELS,CTCF-bound",
                                        "pELS",
                                        "dELS,CTCF-bound",
                                        "dELS",
                                        "DNase-H3K4me3,CTCF-bound",
                                        "DNase-H3K4me3"))
# 绘制饼图
p = ggplot(data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  ggtitle("Liftover human and annotation") + 
  scale_fill_manual(values = c(boxPlot_color))
p 
ggsave("human_pie.pdf",p)


#### mouse FigureD
data= readRDS("macaca_liftover_mouse_CRE.rds")
boxPlot_color =c("CTCF-only,CTCF-bound"="#fbb040",
                 dELS  = "#4078c0","dELS,CTCF-bound" = "#00457c" ,
                 "DNase-H3K4me3" = "#d25238","DNase-H3K4me3,CTCF-bound" = "#ad182d",
                 "pELS" = "#8ba753","pELS,CTCF-bound" = "#30660f",         
                 "PLS"  = "#8a8acb" ,"PLS,CTCF-bound" ="#631d76",
                 "unconserve" = "#caccd1",
                 "regin_couserve" = "#629aa9")

data = data[!grepl("PLS|PLS,CTCF-bound",data$Var1),]
data$Var1 = factor(data$Var1,levels = c("unconserve",
                                        "regin_couserve",
                                        "CTCF-only,CTCF-bound",
                                        "pELS,CTCF-bound",
                                        "pELS",
                                        "dELS,CTCF-bound",
                                        "dELS",
                                        "DNase-H3K4me3,CTCF-bound",
                                        "DNase-H3K4me3"))
# 绘制饼图
p = ggplot(data, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  ggtitle("Liftover human and annotation") + 
  scale_fill_manual(values = c(boxPlot_color))
p


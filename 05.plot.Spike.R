library(tidyverse)
library(data.table)
library(reshape2)

# Inhibited
"#006934"
# Activated
"#5F1985"
# No change 
"#E83828"

######on
data = fread("Raster_P02_240803_195441_A_003_P_13b_ON.csv")
actual_rows <- nrow(data)
desired_rows <- 60
if (actual_rows < desired_rows) {
  rows_to_add <- desired_rows - actual_rows
  additional_rows <- data.frame(matrix(NA, nrow = rows_to_add, ncol = ncol(data)))
  colnames(additional_rows) <- colnames(data)
  data <- rbind(data, additional_rows)
}

data$group1 = rep(paste0("Trial",c(1:10)),6)
data$group2 = rep(c(1:10),6)
data$group2 = as.factor(data$group2)
data$group3 = c(rep("0°",10),rep("30°",10),rep("60°",10),
                rep("90°",10),rep("120°",10),rep("150°",10))

data1 = melt(data)
data1$group3 = factor(data1$group3,levels = c("0°","30°","60°","90°","120°","150°"))
data1$group2 = factor(data1$group2,levels = c(1:10))

pdf("Raster_on.pdf",width = 8,height = 1.7)
ggplot(data1,aes(value,group2)) + 
  geom_rect(aes(xmin = -600, xmax = 1400, ymin = -Inf, ymax = Inf), 
            fill = "#C9CACA", alpha = 0.5)  + 
  geom_point(size =1,color = "#006934",shape = "|", stroke =1) + 
  theme_classic() + facet_wrap(~group3,nrow = 1) + 
  scale_y_discrete(breaks = c(1,2,3,4,5,6,7,8,9,10))+ 
  geom_vline(xintercept=c(0,1000),lty=5,col="black",lwd=0.1) + 
  scale_x_continuous(limits = c(-1000,2000),breaks = c(-1000,-500,0,500,1000,1500, 2000)) + 
  guides(x = guide_axis(angle = 45)) 
dev.off()

######off
data = fread("Raster_P02_240803_195441_A_003_P_13b_off.csv")
actual_rows <- nrow(data)
desired_rows <- 60
if (actual_rows < desired_rows) {
  rows_to_add <- desired_rows - actual_rows
  additional_rows <- data.frame(matrix(NA, nrow = rows_to_add, ncol = ncol(data)))
  colnames(additional_rows) <- colnames(data)
  data <- rbind(data, additional_rows)
}

data$group1 = rep(paste0("Trial",c(1:10)),6)
data$group2 = rep(c(1:10),6)
data$group2 = as.factor(data$group2)
data$group3 = c(rep("0°",10),rep("30°",10),rep("60°",10),
                rep("90°",10),rep("120°",10),rep("150°",10))

data1 = melt(data)
data1$group3 = factor(data1$group3,levels = c("0°","30°","60°","90°","120°","150°"))
data1$group2 = factor(data1$group2,levels = c(1:10))

pdf("Raster_off.pdf",width = 8,height =1.7)
ggplot(data1,aes(value,group2)) + 
  geom_point(size =1,color = "#006934",shape = "|", stroke =1) + 
  theme_classic() + 
  facet_wrap(~group3,nrow = 1) + 
  scale_y_discrete(breaks = c(1,2,3,4,5,6,7,8,9,10))+ 
  geom_vline(xintercept=c(0,1000),lty=5,col="black",lwd=0.1) + 
  scale_x_continuous(limits = c(-1000,2000),breaks = c(-1000,-500,0,500,1000,1500, 2000)) + 
  guides(x = guide_axis(angle = 45)) 
dev.off()


#######  on
data = fread("PSTH_P01_240816_142114_A_022_P_23b_on.csv")
names(data) = c("value1","0°","30°","60°","90°","120°","150°")
data = melt(data,c(1))
max(data$value)

pdf("PSTH_on.pdf",width = 13,height = 2)
ggplot(data,aes(value1,value)) + 
  geom_rect(aes(xmin = -600, xmax = 1400, ymin = -Inf, ymax = Inf), 
            fill = "grey", alpha = 0.2) + ylim(0,60) + 
  geom_col(color = "#E83828") + 
  theme_classic() + facet_wrap(~ variable,nrow = 1) + 
  geom_vline(xintercept=c(0,1000),lty=3,col="black",lwd=0.5) + 
  scale_x_continuous(limits = c(-1000,2000),breaks = c(-1000,-500,0,500,1000,1500, 2000)) + 
  guides(x = guide_axis(angle = 45)) + 
  scale_y_continuous(limits = c(0,60))
dev.off()

#####  off
data = fread("PSTH_P01_240816_142114_A_022_P_23b_off.csv")
head(data)
names(data) = c("value1","0°","30°","60°","90°","120°","150°")
head(data)

data = melt(data,c(1))
max(data$value)
pdf("PSTH_off.pdf",width = 13,height = 2)
ggplot(data,aes(value1,value)) + 
  geom_col(color = "#E83828") + 
  theme_classic() + facet_wrap(~ variable,nrow = 1) + 
  geom_vline(xintercept=c(0,1000),lty=3,col="black",lwd=0.5) + 
  scale_x_continuous(limits = c(-1000,2000),breaks = c(-1000,-500,0,500,1000,1500, 2000)) + 
  scale_y_continuous(limits = c(0,60)) + 
  guides(x = guide_axis(angle = 45))
dev.off()

#### Line chart
data = fread("Oritun_P02_240803_195441_A_004_P_12b_line.csv",header = T)
head(data)
names(data) = c("0°","30°","60°","90°","120°","150°","180°","group")
data = data %>% column_to_rownames(var = "group")
data = as.data.frame(t(data))
data$time = row.names(data)
data1 = data[,c(1,2,5)]
data1$group = "on"
names(data1) = c("mean","sd","group1","group2")
data2 = data[,c(3,4,5)]
data2$group = "off"
names(data2) = c("mean","sd","group1","group2")
data = rbind(data1,data2)
head(data)
data$group1 = factor(data$group1,levels = c("0°","30°","60°","90°","120°","150°","180°"))

pdf("line.pdf",width = 3,height = 2.7)
ggplot(data, aes(x = group1, y = mean, group = group2,color = group2,fill = group2)) +
  geom_line() + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.1) +
  theme_classic() + 
  scale_color_manual(values = c("grey","#E83828")) + ylim(0,2.5)+ 
  guides(x = guide_axis(angle = 45)) 
dev.off()


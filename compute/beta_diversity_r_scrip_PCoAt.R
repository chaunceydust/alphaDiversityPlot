source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
#安装phyloseq
library("phyloseq")
library("ggplot2")
#加载必要的包
setwd("E:/lesson2")
#设置工作目录，即数据存储的目录
qiimedata <- import_qiime(otufilename = "feature-table.taxonomy.txt", mapfilename = "mapping_file.txt", treefilename = "tree.rooted.nwk", refseqfilename = "dna-sequences.fasta")

qiimedata <- import_qiime(otufilename = "../Total.gut.sub3W1.filter_4_noChimera_rmAlignFail.normalized_otu_table.txt", mapfilename = "mapping_file2.txt", treefilename = "../Total.gut.sub3W1_rep_set.filter_4_noChimera_tree.tre", refseqfilename = "../Total.gut.sub3W1_rep_set.filter_4_noChimera_aligned_pfiltered.fasta")
#读取数据,读取数据，参数都是文件名，注意加后缀
#otufilename指定out表格，mapfilename指定map文件（分组数据）， treefilename指定有根进化树文件，refseqfilename指定代表序列文件
otu<-qiimedata@otu_table@.Data
#从导入的数据中提取otu表格
#otu<-otu_table(qiimedata)#可选的otu表格提取方法

sum_of_otus<-colSums(t(otu))
#计算各个otu检测到的总序列数
selected_otu<-names(sum_of_otus)[sum_of_otus>0]

#获取总序列数大于10的otu id
sub_qiimedata <- prune_taxa(selected_otu, qiimedata)
#筛选总序列数大于10的otu phyloseq数据
#sub_qiimedata=subset_taxa(sub_qiimedata,Kingdom=="Bacteria")#根据注释分类进行筛选otu的方法

weighted_unifrac <- distance(sub_qiimedata, method='wunifrac')
#计算样本间加权UniFrac矩阵
unweighted_unifrac <- distance(sub_qiimedata, method='unifrac')
#计算样本间非加权UniFrac矩阵
bray_curtis <- distance(sub_qiimedata, method='bray')
#计算样本间Bray-Curtis距离矩阵
write.table(as.matrix(weighted_unifrac),"weighted_unifrac.txt",sep = '\t',quote = FALSE,col.names = NA)
write.table(as.matrix(unweighted_unifrac),"unweighted_unifrac.txt",sep = '\t',quote = FALSE,col.names = NA)
write.table(as.matrix(bray_curtis),"bray_curtis.txt",sep = '\t',quote = FALSE,col.names = NA)
#保存三个距离矩阵


nmds_of_bray_curtis<-ordinate(physeq=sub_qiimedata,distance = 'bray',method = "NMDS")
#基于Bray-Curtis距离矩阵的NMDS排序分析
p<-plot_ordination(sub_qiimedata, nmds_of_bray_curtis, type="samples", color="Group1") 
p
#将NMDS排序分析结果可视化
p<-p + geom_point(size=3) +ggtitle("NMDS of Bray-Curtis distance") + stat_ellipse()+theme(text = element_text(size = 15))
p
#对图片进行适当修饰， stat_ellipse()加椭圆， ggtitle()加标题
ggsave(plot = p,"nmds_of_bary_curtis.pdf",dpi = 300,width = 7,height = 6)
#保存图片
pcoa_of_unifrac <- ordinate(physeq=sub_qiimedata,distance = 'unifrac',method = "PCoA")
# unifrac, wunifrac
#基于Bray-Curtis距离矩阵的PCoA排序分析
pcoa <- pcoa_of_unifrac
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]

group <- read_tsv("../Total-color3.txt")
plotdata <- data.frame(
  SampleID = rownames(pcoa$vectors),
  PC1,PC2
) %>% 
  right_join(group)
colnames(plotdata)[1:4] <-c("sample","PC1","PC2","Group")

pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)
plotdata$Group <- factor(
  plotdata$Group,
  levels = c("Control", "L.fermentum", "L.animalis",
             "L.salvarius", "L.reuteri")
)

library(dplyr)
yf <- plotdata
yd1 <- yf %>% group_by(Group) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(Group) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1

fit1 <- aov(PC1~Group,data = plotdata)

library(multcomp)
tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
res1 <- cld(tuk1,alpah=0.05)

fit2 <- aov(PC2~Group,data = plotdata)
tuk2<-glht(fit2,linfct=mcp(Group="Tukey"))
res2 <- cld(tuk2,alpah=0.05)


test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                   yd1 = yd1$Max,yd2 = yd2$Max,Group = yd1$Group)
test$Group <- factor(
  test$Group,
  levels = c("Control", "L.fermentum", "L.animalis",
             "L.salvarius", "L.reuteri")
)

cbbPalette <- c(
  "Control" = "#C1CDCD", "L.fermentum" = "#329A54", "L.animalis" = "#2467C8" , "L.salvarius" = "#E4945E", "L.reuteri" = "#9A3254"
)

p1 <- ggplot(plotdata,aes(Group,PC1)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
            size = 7,color = "black",fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(colour='black',size=20,face = "bold"),
        axis.text.x=element_blank(),
        legend.position = "none")

p3 <- ggplot(plotdata,aes(Group,PC2)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
            size = 7,color = "black",fontface = "bold") +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=20,angle = 45,
                                 vjust = 1,hjust = 1,face = "bold"),
        axis.text.y=element_blank(),
        legend.position = "none")

p2<-ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(fill=Group),size=8,pch = 21)+
  scale_fill_manual(values=cbbPalette,name = "Group")+
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
  ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
  theme(text=element_text(size=30))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=34),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=34,vjust = 7),
        axis.title.y=element_text(colour='black', size=34,vjust = -2),
        axis.text=element_text(colour='black',size=28),
        legend.title=element_text(size = 24,face = "bold"),
        legend.text=element_text(size=20),
        legend.key=element_blank(),
        legend.position = c(0.92,0.83),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(1,"cm")) +
  guides(fill = guide_legend(ncol = 1))

library(vegan)
otu.adonis=adonis(unweighted_unifrac~Group,data = group,distance = "unifrac")

p4 <- ggplot(plotdata, aes(PC1, PC2)) +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",otu.adonis$aov.tab$Df[1],
                                               "\nR2 = ",round(otu.adonis$aov.tab$R2[1],4),
                                               "\np-value = ",otu.adonis$aov.tab$`Pr(>F)`[1],sep = "")),
            size = 7) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

library(patchwork)
p5 <- p1 + p4 + p2 + p3 + 
  plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)

pdf("PCoA_unifrac.pdf",height=12,width=15)
p5
dev.off()
png(filename="PCoA_unifrac.png",res=600,height=7000,width=9000)
p5
dev.off()


# weight ------------------------------------------------------------------

pcoa_of_wunifrac <- ordinate(physeq=sub_qiimedata,distance = 'wunifrac',method = "PCoA")
# unifrac, wunifrac
#基于Bray-Curtis距离矩阵的PCoA排序分析
pcoa <- pcoa_of_wunifrac
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]

# group <- read_tsv("../Total-color3.txt")
plotdata <- data.frame(
  SampleID = rownames(pcoa$vectors),
  PC1,PC2
) %>% 
  right_join(group)
colnames(plotdata)[1:4] <-c("sample","PC1","PC2","Group")

pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)
plotdata$Group <- factor(
  plotdata$Group,
  levels = c("Control", "L.fermentum", "L.animalis",
             "L.salvarius", "L.reuteri")
)

library(dplyr)
yf <- plotdata
yd1 <- yf %>% group_by(Group) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(Group) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1

fit1 <- aov(PC1~Group,data = plotdata)

library(multcomp)
tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
res1 <- cld(tuk1,alpah=0.05)

fit2 <- aov(PC2~Group,data = plotdata)
tuk2<-glht(fit2,linfct=mcp(Group="Tukey"))
res2 <- cld(tuk2,alpah=0.05)


test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                   yd1 = yd1$Max,yd2 = yd2$Max,Group = yd1$Group)
test$Group <- factor(
  test$Group,
  levels = c("Control", "L.fermentum", "L.animalis",
             "L.salvarius", "L.reuteri")
)

# cbbPalette <- c(
#   "Control" = "#C1CDCD", "L.fermentum" = "#329A54", "L.animalis" = "#2467C8" , "L.salvarious" = "#E4945E", "L.reuteri" = "#9A3254"
# )

p1 <- ggplot(plotdata,aes(Group,PC1)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd1,label = PC1),
            size = 7,color = "black",fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(colour='black',size=20,face = "bold"),
        axis.text.x=element_blank(),
        legend.position = "none")

p3 <- ggplot(plotdata,aes(Group,PC2)) +
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test,aes(x = Group,y = yd2,label = PC2),
            size = 7,color = "black",fontface = "bold") +
  scale_fill_manual(values=cbbPalette) +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=20,angle = 45,
                                 vjust = 1,hjust = 1,face = "bold"),
        axis.text.y=element_blank(),
        legend.position = "none")

p2<-ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(fill=Group),size=8,pch = 21)+
  scale_fill_manual(values=cbbPalette,name = "Group")+
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  xlim(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range) +
  ylim(ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range) +
  theme(text=element_text(size=30))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=34),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=34,vjust = 7),
        axis.title.y=element_text(colour='black', size=34,vjust = -2),
        axis.text=element_text(colour='black',size=28),
        legend.title=element_text(size = 24,face = "bold"),
        legend.text=element_text(size=20),
        legend.key=element_blank(),
        legend.position = c(0.92,0.83),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(1,"cm")) +
  guides(fill = guide_legend(ncol = 1))

library(vegan)
otu.adonis=adonis(weighted_unifrac~Group,data = group,distance = "wunifrac")

p4 <- ggplot(plotdata, aes(PC1, PC2)) +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",otu.adonis$aov.tab$Df[1],
                                               "\nR2 = ",round(otu.adonis$aov.tab$R2[1],4),
                                               "\np-value = ",otu.adonis$aov.tab$`Pr(>F)`[1],sep = "")),
            size = 7) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

library(patchwork)
p5 <- p1 + p4 + p2 + p3 + 
  plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)

pdf("PCoA_wunifrac.pdf",height=12,width=15)
p5
dev.off()
png(filename="PCoA_wunifrac.png",res=600,height=7000,width=9000)
p5
dev.off()

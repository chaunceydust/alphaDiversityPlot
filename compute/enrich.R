library(tidyverse)
library(vegan)
# library(ape)
# library(ggplot2)
# library(grid)
# library(dplyr)
# library(multcomp)
# library(patchwork)

library(tidyverse)

data1 <- read_tsv("../Total.Genus.absolute.abundance_full.txt")

group <- read_tsv("../Total-color4.txt")

data2 <- data1 %>% 
  select(c("ID", group$SampleID))

otu_table <- select(data2, -ID) %>% 
  as.data.frame()
rownames(otu_table) <- data2$ID

design <- read.table("../Total-color4.txt",
                     header = TRUE, row.names = 1,
                     sep = "\t")
design <- design[1]
colnames(design) <- "SampleType"
design$SampleType <- factor(
  design$SampleType,
  levels = c("Control", "L.fermentum", "L.animalis",
             "L.salvarius", "L.reuteri", "L.reuteriAPS", "APS")
)

idx =which(rownames(design) %in% colnames(otu_table))

#提取分组文件

sub_design <- design[idx, , drop = FALSE]

#提取OTUtable文件

count <- otu_table[, rownames(sub_design), drop = FALSE]

head(count)


library(edgeR)
# create DGElist

d =DGEList(counts=count, group=sub_design$SampleType)

d =calcNormFactors(d)#默认为TMM标准化

# 生成实验设计矩阵

#生成矩阵

design.mat =model.matrix(~ 0 + d$samples$group)

#矩阵命名

colnames(design.mat)=levels(design$SampleType)

d2 =estimateGLMCommonDisp(d, design.mat)

d2 =estimateGLMTagwiseDisp(d2, design.mat)

fit = glmFit(d2,design.mat)

# 设置比较组

BvsA <-makeContrasts(contrasts = "L.fermentum-Control",levels=design.mat)#注意是以GF1为对照做的比较

# 组间比较,统计Foldchange, Pvalue

lrt =glmLRT(fit,contrast=BvsA)

# FDR检验，控制假阳性率小于5%

de_lrt =decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05,lfc=0)#lfc=0这个是默认值

summary(de_lrt)
# 导出计算结果

x=lrt$table

x$sig=de_lrt

#这里加一列调整后的p值

x <- cbind(x,padj = p.adjust(x$PValue, method = "fdr"))

enriched =row.names(subset(x,sig==1))

depleted =row.names(subset(x,sig==-1))

x$level =as.factor(ifelse(x$sig==1, "enriched",ifelse(x$sig==-1,"depleted","nosig")))



# 2 -----------------------------------------------------------------------

# read data

readCount<-read.table(file="../Total.Genus.absolute.abundance_full.txt", 
                      header = T, row.names = 1)
conditions<-read.table(file="../Total-color4.txt", header = F, sep = "\t")
conditions<-factor(t(conditions))

readCount <- otu_table
conditions <- design$SampleType

# edger TMM normalize

# y <- DGEList(counts=readCount,group=conditions)
# ##Remove rows conssitently have zero or very low counts
# keep <- filterByExpr(y)
# y <- y[keep,keep.lib.sizes=FALSE]
# ##Perform TMM normalization and transfer to CPM (Counts Per Million)
# y <- calcNormFactors(y,method="TMM")
# count_norm=cpm(y)
# count_norm<-as.data.frame(count_norm)
count_norm <- otu_table

ddd <- rowMeans(otu_table) %>% 
  unname()
ccc <- data.frame(
  strain = rownames(otu_table),
  mm = ddd
)

top15 <- ccc %>% 
  arrange(desc(mm)) %>% 
  slice_head(n = 15) %>% 
  select(strain) %>% 
  as_vector()
# Run the Wilcoxon rank-sum test for each gene
levels(conditions)
ppp <- c()
for (m in 2:length(levels(conditions))) {
  idx <- which(design$SampleType %in% levels(conditions)[c(1, m)])
  count_norm <- otu_table[idx]
  row_idx <- which(rowMeans(count_norm) > 0)
  count_norm <- count_norm[row_idx,]
  conditions <- design$SampleType[idx]
  pvalues <- sapply(1:nrow(count_norm),function(i){
    data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
    p=wilcox.test(gene~conditions, data)$p.value
    return(p)
  })
  fdr=p.adjust(pvalues,method = "fdr")
  
  conditionsLevel<-levels(conditions)[c(1, m)]
  dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  mybaseMean = rowMeans(dataCon2)
  mylogP = -log10(pvalues)
  mylogFDR = -log10(fdr)
  
  outRst<-data.frame(
    log2foldChange=foldChanges, 
    pValues=pvalues, 
    logP = mylogP,
    baseMean = mybaseMean,
    FDR=fdr,
    logFDR = mylogFDR
  )
  rownames(outRst)=rownames(count_norm)
  outRst$log2foldChange[which(outRst$log2foldChange == Inf)] <- 1.1*max(abs(outRst$log2foldChange[which(is.finite(outRst$log2foldChange))]))
  outRst$log2foldChange[which(outRst$log2foldChange == -Inf)] <- -1.1*max(abs(outRst$log2foldChange[which(is.finite(outRst$log2foldChange))]))
  # outRst=na.omit(outRst)
  outRst$group <- "nosig"
  outRst$group[which(outRst$log2foldChange > 1 & outRst$logP > 1)] <- "enriched"
  outRst$group[which(outRst$log2foldChange < -1 & outRst$logP > 1)] <- "depleted"
  outRst$group <- factor(
    outRst$group,
    levels = c("enriched", "nosig", "depleted")
  )
  
  outRst$label[rownames(outRst) %in% top15] <- rownames(outRst)[rownames(outRst) %in% top15]
  outRst$bar <- 1
  
  fdrThres=0.05
  
  write.table(outRst, 
              paste0(paste0(conditionsLevel, collapse = "-"), ".WilcoxonTest.rst51.tsv"),
              row.names = TRUE, col.names = TRUE,
              sep = "\t",
              quote = FALSE)
  
  library(ggrepel)
  mycol <- c("red", "grey", "blue")
  a <- ggplot(outRst,
              aes(log2foldChange, logP)) +
    geom_vline(xintercept = c(-1, 1),
               linetype = 'dashed',
               color = "grey") +
    geom_hline(yintercept = c(1),
               linetype = 'dashed',
               color = "grey") +
    geom_point(
      aes(size  = baseMean,
          color = group)
    ) +
    geom_label_repel(aes(label = label,
                         color = group,
                         szie = baseMean)) +
    scale_size_continuous(range = c(1, 8)) +
    scale_color_manual(
      values = mycol
    ) +
    # guides(color)+
    labs(
      x = "Log2 fold change",
      y = "-Log10 P"
    ) +
    ggtitle(
      label = conditionsLevel[2],
      subtitle = paste0(
        "(total = ",
        nrow(outRst),
        " genus)"
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(
        hjust = .5,
        size = 18,
        face = 3
      ),
      plot.subtitle = element_text(
        hjust = .5,
        size = 12
      )
    )
  
  b <- ggplot(outRst) +
    geom_bar(
      aes(bar,
          # color = group,
          fill = group),
      position = "stack"
    ) +
    scale_fill_manual(
      values = mycol
    ) +
    ggtitle(
      label = " ",
      subtitle = paste0(
        paste0(names(table(outRst$group)), collapse = ":"),
        "\n(",
        paste0(table(outRst$group), collapse = ":"),
        ")"
      )
    ) +
    guides(
      fill = guide_legend(title = "change")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      plot.subtitle = element_text(hjust = .5, size = 12),
      panel.grid = element_blank(),
      
    ) + 
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))
  
  library(patchwork)
  g1 <- a + b +
    plot_layout(widths = c(4, 1))
  ggsave(
    file = paste0(paste0(conditionsLevel, collapse = "-"), ".enriched51.pdf"),
    width = 6, height = 4
  )
}





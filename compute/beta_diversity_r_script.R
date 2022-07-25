source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
#安装phyloseq
library("phyloseq")
library("ggplot2")
#加载必要的包
setwd("E:/lesson2")
#设置工作目录，即数据存储的目录
qiimedata <- import_qiime(otufilename = "feature-table.taxonomy.txt", mapfilename = "mapping_file.txt", treefilename = "tree.rooted.nwk", refseqfilename = "dna-sequences.fasta")
#读取数据,读取数据，参数都是文件名，注意加后缀
#otufilename指定out表格，mapfilename指定map文件（分组数据）， treefilename指定有根进化树文件，refseqfilename指定代表序列文件
otu<-qiimedata@otu_table@.Data
#从导入的数据中提取otu表格
#otu<-otu_table(qiimedata)#可选的otu表格提取方法

sum_of_otus<-colSums(t(otu))
#计算各个otu检测到的总序列数
selected_otu<-names(sum_of_otus)[sum_of_otus>10]

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
pcoa_of_bray_curtis<-ordinate(physeq=sub_qiimedata,distance = 'bray',method = "PCoA")
#基于Bray-Curtis距离矩阵的PCoA排序分析
p<-plot_ordination(sub_qiimedata, pcoa_of_bray_curtis, type="samples", color="Group1",shape = "Group1") 
p
#将PCoA排序分析结果可视化
p<-p+ scale_colour_manual(values=c("#DC143C","#808000","#00CED1")) + geom_point(size=2) +ggtitle("PCoA of Bray-Curtis distance")+theme(text = element_text(size = 15))
p
#对图片进行适当修饰
#用scale_colour_manual(values=c())自定义颜色，可查颜色的16进制对照表
ggsave(plot = p,"pcoa_of_bary_curtis.pdf",dpi = 300,width = 7,height = 6)
#保存图片

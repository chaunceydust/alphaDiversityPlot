source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
#��װphyloseq
library("phyloseq")
library("ggplot2")
#���ر�Ҫ�İ�
setwd("E:/lesson2")
#���ù���Ŀ¼�������ݴ洢��Ŀ¼
qiimedata <- import_qiime(otufilename = "feature-table.taxonomy.txt", mapfilename = "mapping_file.txt", treefilename = "tree.rooted.nwk", refseqfilename = "dna-sequences.fasta")
#��ȡ����,��ȡ���ݣ����������ļ�����ע��Ӻ�׺
#otufilenameָ��out����mapfilenameָ��map�ļ����������ݣ��� treefilenameָ���и��������ļ���refseqfilenameָ�����������ļ�
otu<-qiimedata@otu_table@.Data
#�ӵ������������ȡotu����
#otu<-otu_table(qiimedata)#��ѡ��otu������ȡ����

sum_of_otus<-colSums(t(otu))
#�������otu��⵽����������
selected_otu<-names(sum_of_otus)[sum_of_otus>10]

#��ȡ������������10��otu id
sub_qiimedata <- prune_taxa(selected_otu, qiimedata)
#ɸѡ������������10��otu phyloseq����
#sub_qiimedata=subset_taxa(sub_qiimedata,Kingdom=="Bacteria")#����ע�ͷ������ɸѡotu�ķ���

weighted_unifrac <- distance(sub_qiimedata, method='wunifrac')
#�����������ȨUniFrac����
unweighted_unifrac <- distance(sub_qiimedata, method='unifrac')
#����������Ǽ�ȨUniFrac����
bray_curtis <- distance(sub_qiimedata, method='bray')
#����������Bray-Curtis�������
write.table(as.matrix(weighted_unifrac),"weighted_unifrac.txt",sep = '\t',quote = FALSE,col.names = NA)
write.table(as.matrix(unweighted_unifrac),"unweighted_unifrac.txt",sep = '\t',quote = FALSE,col.names = NA)
write.table(as.matrix(bray_curtis),"bray_curtis.txt",sep = '\t',quote = FALSE,col.names = NA)
#���������������


nmds_of_bray_curtis<-ordinate(physeq=sub_qiimedata,distance = 'bray',method = "NMDS")
#����Bray-Curtis��������NMDS�������
p<-plot_ordination(sub_qiimedata, nmds_of_bray_curtis, type="samples", color="Group1") 
p
#��NMDS�������������ӻ�
p<-p + geom_point(size=3) +ggtitle("NMDS of Bray-Curtis distance") + stat_ellipse()+theme(text = element_text(size = 15))
p
#��ͼƬ�����ʵ����Σ� stat_ellipse()����Բ�� ggtitle()�ӱ���
ggsave(plot = p,"nmds_of_bary_curtis.pdf",dpi = 300,width = 7,height = 6)
#����ͼƬ
pcoa_of_bray_curtis<-ordinate(physeq=sub_qiimedata,distance = 'bray',method = "PCoA")
#����Bray-Curtis��������PCoA�������
p<-plot_ordination(sub_qiimedata, pcoa_of_bray_curtis, type="samples", color="Group1",shape = "Group1") 
p
#��PCoA�������������ӻ�
p<-p+ scale_colour_manual(values=c("#DC143C","#808000","#00CED1")) + geom_point(size=2) +ggtitle("PCoA of Bray-Curtis distance")+theme(text = element_text(size = 15))
p
#��ͼƬ�����ʵ�����
#��scale_colour_manual(values=c())�Զ�����ɫ���ɲ���ɫ��16���ƶ��ձ�
ggsave(plot = p,"pcoa_of_bary_curtis.pdf",dpi = 300,width = 7,height = 6)
#����ͼƬ
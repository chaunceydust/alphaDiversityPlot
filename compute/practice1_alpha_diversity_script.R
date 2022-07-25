
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("agricolae")
install.packages("vegan")
#��װ��Ҫ�İ������Ѱ�װ����ʡ����һ��
library(ggplot2)
library(ggpubr)
library(agricolae)
library(vegan)
#������Ҫ�ķ������������û�а�װ��Щ������ʹ��

otu=read.table('E:/lesson1/feature-table.taxonomy.txt',row.names = 1,skip=1,header=T,comment.char='',sep='\t')

#'E:/lesson1/feature-table.taxonomy.txt'Ϊ�ļ�·����ע��б�߷���
#row.names = 1ָ����һ��Ϊ����
#skip=1������һ�в���
#header=Tָ����һ����Ч��Ϊ����
#sep='\t'��ʾָ���Ʊ���Ϊ�ָ���
#comment.char=''��ʾ����ע�ͷ���Ϊ���ַ�����������#��������ݾͲ��ᱻʡ��

otu=otu[,-ncol(otu)]
#ɾ���������
otu=t(otu)
#ת��

shannon=diversity(otu,"shannon")
#������ũָ��
simpson=diversity(otu,"simpson")
#��������ɭָ��
alpha=data.frame(shannon,simpson)
#�ϲ�����
write.table(alpha,'E:/lesson1/alpha-summary.tsv',sep = '\t',quote=F)
#������


map<-read.table('E:/lesson1/mapping_file.txt',row.names = 1,header = T,sep='\t',comment.char='',check.names=F)
#��ȡ�������
#row.names = 1��ʾָ����һ��Ϊ����
#header = T��ʾָ����һ����Ч��Ϊ����
#sep='\t'��ʾָ���Ʊ���Ϊ�ָ���
#comment.char=''��ʾ����ע�ͷ���Ϊ���ַ�����������#��������ݾͲ��ᱻʡ��
#check.names=F��ʾ��ȡ�����в����������������κ��޸�
group<-map['Group1']
#��ȡ��Ҫ�ķ��飬'Group1'Ϊ���з�������
alpha<-read.table('E:/lesson1/alpha-summary.tsv',header = T,row.names = 1,sep = '\t')
#��ȡalpha�����Ա�
alpha<-alpha[match(rownames(group),rownames(alpha)),]
#����alpha���е�˳��ʹ����group������id��������˳��һ��
data<-data.frame(group,alpha)
#�ϲ���������


p=ggplot(data = data,aes(x=Group1,y=shannon))+geom_boxplot(fill=rainbow(7)[5])
p
#data = dataָ�����ݱ���
#x=Group1ָ����Ϊx�����������
#y=shannonָ����Ϊy�����������
#geom_boxplot()��ʾ������ͼ
mycompare=list(c('A','B'),c('A','C'),c('B','C'))
#ָ�����رȽϵķ����
p<-p+stat_compare_means(comparisons=mycompare,label = "p.signif",method = 'wilcox')
p
#�����Ǻţ����������Ա�ǵĵ�һ�ַ���,ʹ��wilcoxon�ǲ������鷽��

anova <- aov(shannon~Group1,data = data)
plotdata<-duncan.test(anova,"Group1",console = TRUE, alpha = 0.05)
plotdata<-data.frame(id=rownames(plotdata$groups),plotdata$groups)
p<-p+geom_text(data = plotdata,aes(x=id,y=5,label=groups))
p

#������������ĸ�����������Ա�ǵĵڶ��ַ�����ע��˷���Ϊ�������飬Ҫ��alpha������ָ��������̫�ֲ�

p=p+xlab("")+ylab("Shannon")
#���ĺ������������������
p=p+theme(text = element_text(size = 15,face = "bold"))
p
#�޸��ֺ�����
ggsave(plot = p,'E:/lesson1/shannon_boxplot2.png',height = 7,width = 6,dpi = 300)
#plot = pָ���ոյ������κõ�p
# ��E:/practice1/shannon_boxplot.pdf����ʾ�洢·�����ļ����֣�
#ע���ļ��ĺ�׺������������ͼƬ��ʽ
#���鴢��Ϊpdf��ʽ������ԭ�����ע����Ľ���



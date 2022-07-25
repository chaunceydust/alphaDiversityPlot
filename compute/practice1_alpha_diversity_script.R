
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("agricolae")
install.packages("vegan")
#安装需要的包，如已安装，可省略这一步
library(ggplot2)
library(ggpubr)
library(agricolae)
library(vegan)
#加载需要的分析包，如果还没有安装这些包，请使用

otu=read.table('E:/lesson1/feature-table.taxonomy.txt',row.names = 1,skip=1,header=T,comment.char='',sep='\t')

#'E:/lesson1/feature-table.taxonomy.txt'为文件路径，注意斜线方向
#row.names = 1指定第一列为行名
#skip=1跳过第一行不读
#header=T指定第一个有效行为列名
#sep='\t'表示指定制表符为分隔符
#comment.char=''表示设置注释符号为空字符‘’，这样#后面的内容就不会被省略

otu=otu[,-ncol(otu)]
#删除多余的列
otu=t(otu)
#转置

shannon=diversity(otu,"shannon")
#计算香农指数
simpson=diversity(otu,"simpson")
#计算辛普森指数
alpha=data.frame(shannon,simpson)
#合并数据
write.table(alpha,'E:/lesson1/alpha-summary.tsv',sep = '\t',quote=F)
#储存结果


map<-read.table('E:/lesson1/mapping_file.txt',row.names = 1,header = T,sep='\t',comment.char='',check.names=F)
#读取分组表格
#row.names = 1表示指定第一列为行名
#header = T表示指定第一个有效行为列名
#sep='\t'表示指定制表符为分隔符
#comment.char=''表示设置注释符号为空字符‘’，这样#后面的内容就不会被省略
#check.names=F表示读取过程中不对行名和列名做任何修改
group<-map['Group1']
#提取需要的分组，'Group1'为表中分组列名
alpha<-read.table('E:/lesson1/alpha-summary.tsv',header = T,row.names = 1,sep = '\t')
#读取alpha多样性表
alpha<-alpha[match(rownames(group),rownames(alpha)),]
#重排alpha的行的顺序，使其与group的样本id（行名）顺序一致
data<-data.frame(group,alpha)
#合并两个表格


p=ggplot(data = data,aes(x=Group1,y=shannon))+geom_boxplot(fill=rainbow(7)[5])
p
#data = data指定数据表格
#x=Group1指定作为x轴的数据列名
#y=shannon指定作为y轴的数据列名
#geom_boxplot()表示画箱线图
mycompare=list(c('A','B'),c('A','C'),c('B','C'))
#指定多重比较的分组对
p<-p+stat_compare_means(comparisons=mycompare,label = "p.signif",method = 'wilcox')
p
#添加星号，添加显著性标记的第一种方法,使用wilcoxon非参数检验方法

anova <- aov(shannon~Group1,data = data)
plotdata<-duncan.test(anova,"Group1",console = TRUE, alpha = 0.05)
plotdata<-data.frame(id=rownames(plotdata$groups),plotdata$groups)
p<-p+geom_text(data = plotdata,aes(x=id,y=5,label=groups))
p

#添加显著性字母，添加显著性标记的第二种方法，注意此方法为参数检验，要求alpha多样性指数服从正太分布

p=p+xlab("")+ylab("Shannon")
#更改横坐标标题和纵坐标标题
p=p+theme(text = element_text(size = 15,face = "bold"))
p
#修改字号字形
ggsave(plot = p,'E:/lesson1/shannon_boxplot2.png',height = 7,width = 6,dpi = 300)
#plot = p指定刚刚叠加修饰好的p
# ‘E:/practice1/shannon_boxplot.pdf’表示存储路径和文件名字，
#注意文件的后缀名，它将决定图片格式
#建议储存为pdf格式，具体原因，请关注后面的讲座




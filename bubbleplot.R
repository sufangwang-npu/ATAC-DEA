##KEGG pathway气泡图
library(ggplot2)
pathway <- read.csv("KEGG.csv",header=T,sep=",")
# 画图
jpeg(filename="Fig_kegg.jpeg",width=5400,height=5400,res=900)

pathway <- BiologicalPath
p <- ggplot(pathway,aes(pathway$count.InDataset,pathway$path.term))
p + geom_point()
# 改变点的大小
p + geom_point(aes(size=pathway$count.InDataset))
# 四维数据的展示
Pvalue <- pathway$pvalue
pbubble = p + geom_point(aes(size=pathway$count.InDataset,color=-1*log10(Pvalue)))
# 自定义渐变颜色
pbubble =pbubble+ scale_colour_gradient(low="blue",high="red")
# 绘制pathway富集散点图
pbubble + scale_colour_gradient(low="blue",high="red") + labs(color=expression(-log[10](Pvalue)),size="Gene number",x="Rich factor",y="Pathway name",title="KEGG pathway enrichment")
dev.off()

ggplot(pathway,aes(count.InDataset,path.term))+
  geom_point()+
  geom_point(aes(size=count.InDataset))+
  geom_point(aes(size=count.InDataset,color=-1*log10(pvalue)))+
  scale_colour_gradient(low="blue",high="red")+
  scale_colour_gradient(low="blue",high="red") +
  labs(color=expression(-log[10](pvalue)),size="Gene number",x="Rich factor",y="Pathway name",title="KEGG pathway enrichment")
dev.off()


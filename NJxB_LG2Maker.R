library('ggplot2')
library('reshape2')
library('patchwork')
setwd('F:/01.work/04.DS_HBL比较分析/01.Data/18.GenomeEvaluvationAndWholeGenomeStat/NBS')
data=read.table("NJxB_LG2Maker.txt",sep="\t",header=T,check.names=F,quote="")
ChrLen=read.table('LG2Len.txt',header = T)
data
ChrLen$Pdan_hap2
ggplot() + 
  geom_point(data = data, aes(x= Chromosome, y=Position,color=Type) , shape = 95, size = 10) + 
  geom_segment(data=ChrLen,aes(x=Chr,y=Start,xend=Chr,yend=LG2),size=1 )+
  labs( x = 'Chromosome', y = 'Chromosomal Position(Mb)')  + 
  theme_test() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + 
  guides(color = guide_legend(title = 'NBS Type'))
ggsave('NJxB_LG2Maker.pdf',width=4,height=8)  

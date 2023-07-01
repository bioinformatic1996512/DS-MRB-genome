library('ggplot2')
library('reshape2')
library('patchwork')
setwd('D:/01.work/04.DS_HBL/01.Data/18.GenomeEvaluvationAndWholeGenomeStat/NBS')
data=read.table("Pdan_hap1.NBS.bed",sep="\t",header=T,check.names=F,quote="")
ChrLen=read.table('GenomeChrLen.txt',header = T)

ChrLen$Pdan_hap2
ggplot() + 
  geom_point(data = data, aes(x= Chr, y=Start/1000000, color =Type) , shape = 95, size = 10) + 
  geom_segment(data=ChrLen,aes(x=Chr,y=Start,xend=Chr,yend=HBL_hap1/1000000),size=1 )+
  labs( x = 'Chromosome', y = 'Chromosomal Position(Mb)')  + 
  theme_test() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + 
  guides(color = guide_legend(title = 'NBS Type'))
ggsave('Pdan_hap1_NBSgeneDist.pdf',width=9,height=4)  

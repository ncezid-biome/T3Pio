library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)

df1 <- read.delim('BioNum_3_Core_Comp')
df2 <- read.delim('Sean_Allele_3_Col_Comp_Primers')

df3 <- merge(x=df1,y=df2,by=c('ISO1','ISO2'))

df4 <- df3 %>% select(ISO1,ISO2,BN,SL) %>% filter(ISO1 != ISO2)

df5 <- subset(df4, BN <= 25, select=c(ISO1,ISO2,BN,SL))

df6 <- subset(df4, SL <= 6, select=c(ISO1,ISO2,BN,SL))

cn <- df6 %>% group_by(BN,SL) %>% summarize(count=n())

df7 <- df5 %>% group_by(BN,SL) %>% summarize(count=n())

write.table(df4,file='CorePairwise',quote=FALSE,row.names=FALSE)

mydf <- read.table('CorePairwise',sep=' ',skip=1,header=FALSE,stringsAsFactors = FALSE,col.names=c('ISO1','ISO2','BN','SL'))

sero <- read.table('Sal_266_genome_serotypes.txt',header=TRUE,sep='\t',stringsAsFactors = FALSE)

mydf_sero <- left_join(mydf, sero, by = c("ISO1" = "Genome"))
mydf_sero <- select(mydf_sero, -Serotype)
mydf_sero <- rename(mydf_sero, Iso1_sero = Serotype_code)
mydf_sero <- left_join(mydf_sero, sero, by = c("ISO2" = "Genome"))
mydf_sero <- select(mydf_sero, -Serotype)
mydf_sero <- rename(mydf_sero, Iso2_sero = Serotype_code)

mydf_sero <- mutate(mydf_sero,Serotype_comp = if_else(Iso1_sero==Iso2_sero,'Same Serotype Comparison','Different Serotype Comparison'))

ggplot(cn, aes(BN,SL, size=count))+
  geom_point(alpha=0.4, shape=20,color='darkblue')+
  scale_size_continuous(breaks = c(1,2))+
  scale_y_continuous(breaks = seq(0,20,1))+
  scale_x_continuous(breaks=seq(0,35,5))+
  theme(axis.text = element_text(size=10,face='bold'))+
  theme(axis.title = element_text(size=14,face='bold'))+
  labs(x='BioNumerics Allele Difference', y = 'Pipeline Allele Difference',size='Frequency')+
  

ggsave('SLZoom95BNCore.png',width=10,height =9,dpi=300)

ggplot(df7, aes(BN,SL, size=count,fill=BN))+
  guides(fill=FALSE)+
  geom_point(alpha=0.4,shape=21) +
  scale_size_continuous(breaks = c(1,2,3))+
  scale_x_continuous(breaks = seq(0,25,5)) +
  scale_y_continuous(breaks = seq(0,40,1)) +
  labs(x='BioNumerics Allele Difference', y = 'Pipeline Allele Difference',size='Frequency')+
  theme(legend.background = element_rect(fill='lightblue',size=0.5,linetype ='solid'))


ggplot(df4, aes(BN,SL))+
  geom_point(size=3.0,stroke=0,shape='+',color='#0072B2')+
  labs(x='Bionumerics Allele Differences',y='Pipeline Allele Differences')+
  scale_y_continuous(breaks = seq(0,4000,100))+
  scale_x_continuous(breaks = seq(0,5000,500))+
  ggtitle('Pairwise Comparison of Allele Differences for \n BioNumerics wgMLST and SL Pipeline')+
  theme(plot.title = element_text(hjust = (0.5 )))

ggsave('PipelinevBNCore.png',width=8,height=10,dpi=300)

ggplot(mydf_sero, aes(x=BN,y=SL))+
    geom_point(size=3.0,stroke=0,shape='+',aes(color=Serotype_comp))+
    labs(x='Bionumerics Allele Differences',y='Pipeline Allele Differences',color='Serotype Comparison')+
    scale_y_continuous(breaks = seq(0,4000,100))+
    scale_x_continuous(breaks=seq(0,5000,500))+
    ggtitle('Pairwise Comparison with Comparisons Between \n Same and Different Serogroups')+
    theme(plot.title = element_text(hjust = (0.5)))+
    scale_color_manual(values=c('#0072B2','#D55E00'))+
   
  
ggsave('PipelineVBNCorewithSero.png',width=8,height=10,dpi=300)


write.table(df6,file='Pipeline17Under',quote=FALSE,row.names=FALSE)

underS <- read.table('Pipeline17Under',sep=' ',skip =1, stringsAsFactors = FALSE,header=FALSE,col.names=c('ISO1','ISO2','BN','SL'))

underS_sero <- left_join(underS, sero, by = c("ISO1" = "Genome"))
underS_sero <- select(underS_sero, -Serotype_code)
underS_sero <- rename(underS_sero, Iso1_sero = Serotype)
underS_sero <- left_join(underS_sero, sero, by = c("ISO2" = "Genome"))
underS_sero <- select(underS_sero, -Serotype_code)
underS_sero <- rename(underS_sero, Iso2_sero = Serotype)

subBN <- subset(underS_sero, BN >= 20,select = c(ISO1,ISO2,BN,SL,Iso1_sero,Iso2_sero))


ggplot(subBN,aes(x=BN,y=SL))+
  geom_point(size=5.0,stroke=0,shape='+',aes(color=Iso1_sero))+
  labs(x='Bionumerics Allele Differences',y='Pipeline Allele Differences',color='Serotype')+
  scale_y_continuous(breaks = seq(10,20,1))+
  scale_x_continuous(breaks=seq(20,40,2))+
  ggtitle('Pairwise Comparison with Comparisons Between \n Same and Different Serogroups')+
  theme(plot.title = element_text(hjust = (0.5)))+

ggsave('ProblematicIsolates.png',width=8,height=10,dpi=300)



ggplot(subSFAF, aes(x=BN.x,y=SL.x))+
  geom_point(size=5.0,stroke=0,shape='+',aes(color=Serotype_comp))+
  labs(x='Bionumerics Allele Differences',y='Pipeline Allele Differences',color='Serotype Comparison')+
  scale_y_continuous(breaks = seq(0,10,2))+
  scale_x_continuous(breaks=seq(0,70,10))+
  ggtitle('Pairwise Comparison with Comparisons Between \n Same and Different Serogroups')+
  theme(plot.title = element_text(hjust = (0.5)))+
  scale_color_manual(values=c('#D55E00','#0072B2'))+
  theme_bw()

subSFAF <- subset(df_with_SFAF_info, SL.x <= 6,select = c(ISO1,ISO2,BN.x,SL.x,Serotype_comp))

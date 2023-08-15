#unsupervised learning approaches for HCF drug data
install.packages('ggplot2')
install.packages('tidyverse')
install.packages('ggthemes')
install.packages('ggrepel')
install.packages('stats')
install.packages("factoextra")
install.packages("scales")
install.packages('ggpubr')
install.packages("vegan")
install.packages('cluster')
install.packages('dplyr')
remove.packages('rlang')
install.packages('rlang')
install.packages('foreign')
install.packages('rstatix')
install.packages('mclust')
install.packages('Rtsne')
install.packages('rgl')
install.packages('ppclust')
install.packages('NbClust')
install.packages('ClusterR')
#can skip the installation step if these are already installed
library(ggplot2)
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(stats)
library(factoextra)
library(scales)
library(ggpubr)
library (vegan)
library (cluster)
library(dplyr)
library(rlang)
library(mclust)
library(Rtsne)
library(rgl)
library(ppclust)
library(NbClust)
library(ClusterR)

setwd('Y:\\Anders\\Drug Paper\\Code')
########### Plot the K-means 9 clustering for loadings, figures for drug paper #######################

### Plot scores after applying PCA on the reduced feature set 8-5-2022 ###



df=read_csv('PCA_Scores.csv')
colnames(df)=c("X1","DG1","PC1","PC2","PC3","PC4","PC5","PC6")

drugs=df$DG1
drugsNew=list()
for (d in drugs){
  if(grepl('_1',d)){
    dose='Low'
  }else if(grepl('_2',d)){dose='Med'
  }else if (grepl('_3',d)){dose='High'
  }else{dose=''}
    
  drug=strsplit(d,'_')[[1]][1]
  
  
  if(grepl('TGFB',d) & !grepl('IL1',d)){
    cyto='TGFB'
  }else if(grepl('IL1',d) & !grepl('TGFB',d)){
    cyto='IL1'
  }else if (grepl('IL1',d) & grepl('TGFB',d)){
    cyto='TGFB+IL1'
  }else{cyto=''}
  
  if (cyto == ''){
  treat=paste(c(drug,'_',dose),collapse = '')}else{
    treat=paste(c(drug,'_',dose,'+',cyto),collapse = '')}
  
  if (grepl('Val_',treat)){treat=gsub('Val_','ValBNP_',treat)}
  
drugsNew=c(drugsNew,treat)
}

df$newDrugs=drugsNew


f_name=paste0('Fig3S_Scores_12.pdf')
pdf(f_name,width=18,height=16)
plt=ggplot(df)+geom_point(aes_string(x='PC1',y='PC2'),size=4,alpha=0.6)+
  geom_text_repel(aes_string(x='PC1',y='PC2',label='newDrugs'),max.overlaps=10000,max.iter=10000,point.size = 6, size=4, box.padding = 0.3,point.padding=0.5)+
  theme_bw()+xlab('PC1: 39% Variance Explained')+ylab('PC2: 20% Variance Explained')
print(plt)
dev.off()

f_name=paste0('Fig3S_Scores_13.pdf')
pdf(f_name,width=18,height=16)
plt=ggplot(df)+geom_point(aes_string(x='PC1',y='PC3'),size=4,alpha=0.6)+
  geom_text_repel(aes_string(x='PC1',y='PC3',label='newDrugs'),max.overlaps=10000,max.iter=10000,point.size = 6, size=4, box.padding = 0.3,point.padding=0.5)+
  theme_bw()+xlab('PC1: 39% Variance Explained')+ylab('PC3: 13% Variance Explained')
print(plt)
dev.off()


#Publication PDFs with larger points for the groups of interest#

treatFilter=c('control','IL1','TGFB','TGFB_IL1',
              'Fasudil_3_TGFB','Pirfenidone_3_TGFB',
              'SB203580_3_TGFB','HY12289A_3_TGFB',
              'WH4023_3_TGFB','Glutathione_3_TGFB')
dfb=df[df$DG1 %in% treatFilter,] #df for bigger points
dfs=df[!df$DG1 %in% treatFilter,] #df for smaller points

f_name=paste0('Fig3_Scores_12.pdf')
pdf(f_name,width=12,height=8)
plt=ggplot()+geom_point(data=dfs,aes_string(x='PC1',y='PC2'),size=7,alpha=0.7)+
  geom_point(data=dfb,aes_string(x='PC1',y='PC2'),size=12,alpha=0.7)+
  theme_bw()+xlab('PC1: 39% Variance Explained')+ylab('PC2: 20% Variance Explained')
print(plt)
dev.off()


f_name=paste0('Fig3_Scores_13.pdff')
pdf(f_name,width=12,height=8)
plt=ggplot()+geom_point(data=dfs,aes_string(x='PC1',y='PC3'),size=7,alpha=0.7)+
  geom_point(data=dfb,aes_string(x='PC1',y='PC3'),size=12,alpha=0.7)+
  theme_bw()+xlab('PC1: 39% Variance Explained')+ylab('PC2: 13% Variance Explained')
print(plt)
dev.off()


### Plots for the loadings ###


df=read_csv('PCA_Loadings.csv')
colnames(df)=c("X1","PC1","PC2","PC3","PC4","PC5","PC6")

#change the feature names to the wording used in the paper
featureLabels=c('Actin Long Angular Second Moment',
                'aSMA Long Angular Second Moment',
                'aSMA Long Correlation',
                'Actin Short Information Measure',
                'Actin Medium Information Measure',
                'aSMA Long Information Measure',
                'Actin Long Inverse Difference Moment',
                'aSMA Long Inverse Difference Moment',
                'aSMA Long Sum Average',
                'Actin Medium Sum Entropy',
                'Actin Integrated Intensity',
                'Collagen Integrated Intensity',
                'aSMA Integrated Intensity',
                'Collagen Lower Quartile Intensity',
                'Collagen Mean Intensity',
                'aSMA Mean Intensity',
                'Collagen Minimum Intensity',
                'Cell Area')
df$featLabels=featureLabels

f_name=paste0('Fig3S_Laodings_12.pdf')
pdf(f_name,width=12,height=8)
plt=ggplot()+geom_point(data=df,aes_string(x='PC1',y='PC2'),size=7,alpha=0.7)+
  geom_text_repel(data=df,aes_string(x='PC1',y='PC2',label='featLabels'),max.overlaps=10000,max.iter=10000,point.size = 6, size=4, box.padding = 0.3,point.padding=0.5)+
  theme_bw()+xlab('PC1: 39% Variance Explained')+ylab('PC2: 20% Variance Explained')
print(plt)
dev.off()

f_name=paste0('Fig3S_Laodings_13.pdf')
pdf(f_name,width=12,height=8)
plt=ggplot()+geom_point(data=df,aes_string(x='PC1',y='PC3'),size=7,alpha=0.7)+
  geom_text_repel(data=df,aes_string(x='PC1',y='PC3',label='featLabels'),max.overlaps=10000,max.iter=10000,point.size = 6, size=4, box.padding = 0.3,point.padding=0.5)+
  theme_bw()+xlab('PC1: 39% Variance Explained')+ylab('PC3: 13% Variance Explained')
print(plt)
dev.off()


featureFilter=c('AngularSecondMoment_Actin_10',
                'Intensity_IntegratedIntensity_Actin',
                'Intensity_IntegratedIntensity_Collagen',
                'Intensity_IntegratedIntensity_aSMA')


dfb=df[df$X1 %in% featureFilter,] #df for bigger points
dfs=df[!df$X1 %in% featureFilter,] #df for smaller points

f_name=paste0('Fig3_Laodings_12.pdf')
pdf(f_name,width=12,height=8)
plt=ggplot()+geom_point(data=dfs,aes_string(x='PC1',y='PC2'),size=7,alpha=0.7)+
  geom_point(data=dfb,aes_string(x='PC1',y='PC2'),size=12,alpha=0.7)+
  theme_bw()+xlab('PC1: 39% Variance Explained')+ylab('PC2: 20% Variance Explained')
print(plt)
dev.off()

f_name=paste0('Fig3_Laodings_13.pdf')
pdf(f_name,width=12,height=8)
plt=ggplot()+geom_point(data=dfs,aes_string(x='PC1',y='PC3'),size=7,alpha=0.7)+
  geom_point(data=dfb,aes_string(x='PC1',y='PC3'),size=12,alpha=0.7)+
  theme_bw()+xlab('PC1: 39% Variance Explained')+ylab('PC3: 13% Variance Explained')
print(plt)
dev.off()




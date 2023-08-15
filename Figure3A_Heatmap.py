# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 13:15:35 2022

@author: arn4va
"""
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt


df=pd.read_csv('ProcessedDrugScreen.csv')
from scipy.spatial import distance
from scipy.cluster import hierarchy
df=df.drop('AreaShape_Orientation',axis=1)
df2=df
df2.index=df2.DG1
df2=df2.drop('DG1',axis=1)
df2=df2.transform(stats.zscore,axis=0)

palette = ['#737373','#e31a1c','#3690c0','#88419d']
cytokines=[]
for c in df2.index:
    if 'control' in c:
        cytokines.append('#737373')
    if (('tgfb' in c.lower()) & ('il1' not in c.lower())):
        cytokines.append('#e31a1c')
    if (('tgfb' not in c.lower()) & ('il1' in c.lower())):
        cytokines.append('#3690c0')
    if (('tgfb' in c.lower()) & ('il1' in c.lower())):
        cytokines.append('#88419d')
        
        

channels=[]
for c in df2.columns:
    if 'Collagen' in c:
        channels.append('#41ab5d')
    elif 'Actin' in c:
        channels.append('#08306b')
    elif 'aSMA' in c:
        channels.append('#ec7014')
    else:
        channels.append('#9970ab')
        
        
features=df2.columns        
featNew=[]
for f in features:
    if 'AreaShape_' in f:
        new=f.split('_')[1]
    else:
        new=f
    new=new.replace('2','Short')
    new=new.replace('6','Medium') 
    new=new.replace('10','Long')
    featNew.append(new)
df4=df2
df4.columns=featNew
#rename the drug treatment groups
treatments=df2.index
newDrugs=[]
for d in treatments:
    if '_1_' in d:
        dose='Low'
    elif '_2_' in d:
        dose='Medium'
    elif '_3_' in d:
        dose='High'
    else:
        dose=''
    drug=d.split('_')[0]
    
    if (('TGFB' in d) & ('IL1' not in d)):
        cyto='TGFB'
    elif (('IL1' in d) & ('TGFB' not in d)):
        cyto='IL1'
    elif (('TGFB' in d) & ('IL1' in d)):
        cyto='TGFB+IL1'
    else:
        cyto=''
    
    if cyto=='':
        treat=drug+'_'+dose
    else:
        treat=drug+'_'+dose+'+'+cyto
    
    if 'Val_' in treat:
        treat.replace('Val_','ValBNP_')
    newDrugs.append(treat)
#row_colors=pd.DataFrame({'row_colors':row_colors.values}, index=df2.columns)

df4.index=newDrugs

fig,ax = plt.subplots(1,1,figsize=(35,35))
g=sns.clustermap(data=df4.T,cmap=sns.diverging_palette(220, 20, as_cmap=True), #need index reset to show row colors
                 figsize=(40,31),yticklabels=True,xticklabels=True,vmin=-2,vmax=2,col_colors=cytokines,row_colors=channels)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90,fontsize=18) 
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0,fontsize=18) 
#plt.setp(g.ax_heatmap.set_yticklabels(clusterNames),rotation=0,fontsize=14) 
plt.setp(g.ax_heatmap.set_ylabel(''),fontsize=16) 
plt.setp(g.ax_heatmap.set_xlabel(''),fontsize=14) 
plt.tight_layout()

plt.savefig('Fig3A.png')
#plt.savefig('Fig3A.svg')
plt.close(fig)


#Make reduced heatmap for the supplement

df=pd.read_csv('ProcessedDrugScreen.csv')
df=df.drop('AreaShape_Orientation',axis=1)
df2=df
df2.index=df2.DG1
df2=df2.drop('DG1',axis=1)
df2=df2.transform(stats.zscore,axis=0)

dfh=pd.read_csv('RetainedFeatures_Edited.csv')

df3=df2[df2.columns[df2.columns.isin(dfh[dfh.Keep=='Yes'].Feature)]]


features=df3.columns        
featNew=[]
for f in features:
    if 'AreaShape_' in f:
        new=f.split('_')[1]
    else:
        new=f
    new=new.replace('2','Short')
    new=new.replace('6','Medium') 
    new=new.replace('10','Long')
    featNew.append(new)
df4=df3
df4.columns=featNew

channels=[]
for c in df3.columns:
    if 'Collagen' in c:
        channels.append('#41ab5d')
    elif 'Actin' in c:
        channels.append('#08306b')
    elif 'aSMA' in c:
        channels.append('#ec7014')
    else:
        channels.append('#9970ab')
df4.index=newDrugs

fig,ax = plt.subplots(1,1,figsize=(48,48))
g=sns.clustermap(data=df4.T,cmap=sns.diverging_palette(220, 20, as_cmap=True), #need index reset to show row colors
                 figsize=(54,48),yticklabels=True,xticklabels=True,vmin=-2,vmax=2,col_colors=cytokines,row_colors=channels)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90,fontsize=21) 
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0,fontsize=28) 
#plt.setp(g.ax_heatmap.set_yticklabels(clusterNames),rotation=0,fontsize=14) 
plt.setp(g.ax_heatmap.set_ylabel(''),fontsize=16) 
plt.setp(g.ax_heatmap.set_xlabel(''),fontsize=14) 
plt.tight_layout()

plt.savefig('Figure_S2.png')
plt.savefig('Figure_S2.svg')
plt.close(fig)
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 15:17:37 2022

@author: arn4va
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import cluster, mixture, metrics


dt=pd.read_csv('ProcessedDrugScreen.csv')
dfh=pd.read_csv('RetainedFeatures_Edited.csv')
dt.index=dt.DG1
dt=dt.drop('DG1',axis=1)
dt=dt.drop(['AreaShape_Orientation'],axis=1)
#dt=dt.drop(['AreaShape_Orientation','plate_number','Metadata_Well'],axis=1)
df=dt
corr=df.corr()
x=corr.iloc[:,1:]
x=x.values



features=corr.columns        
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
corr.columns=featNew
corr.index=featNew



#plot the corr matrix with clustermap to get an idea of cluster number
fig=plt.figure(figsize=(30,30))
# my_colors=['#4292c6','#bdbdbd','#cb181d']


g=sns.clustermap(data=corr,cmap=sns.diverging_palette(220, 20, as_cmap=True),xticklabels=True,yticklabels=True,figsize=(30,30))
#set colors on yticks
tickLabels=g.ax_heatmap.axes.get_yticklabels()
row_colors=[]
for i in range(0,len(tickLabels)):
    tick_label=tickLabels[i]
    text=tick_label.get_text()
    filt=dfh[dfh.Feature==text]
    if filt.Keep.values == 'No':
        row_colors.append('#000000')
    else:
        row_colors.append('#cb181d')

for i in range(0,len(tickLabels)):
    tick_label=tickLabels[i]
    tick_text = tick_label.get_text()
    tick_label.set_color(row_colors[i])
    #set colors on xticks
tickLabels=g.ax_heatmap.axes.get_xticklabels()
for i in range(0,len(tickLabels)):
    tick_label=tickLabels[i]
    tick_text = tick_label.get_text()
    tick_label.set_color(row_colors[i])
    
# change the feature label names
## Xticks
tickLabels=g.ax_heatmap.axes.get_xticklabels()
for i in range(0,len(tickLabels)):
    tick_label=tickLabels[i]
    tick_text = tick_label.get_text()
    if 'AreaShape_' in tick_text:
        new=tick_text.split('_')[1]
    else:
        new=tick_text
    new=new.replace('2','Short')
    new=new.replace('6','Medium') 
    new=new.replace('10','Long')
    tick_label.set_text(new)
## Yticks
tickLabels=g.ax_heatmap.axes.get_yticklabels()
for i in range(0,len(tickLabels)):
    tick_label=tickLabels[i]
    tick_text = tick_label.get_text()
    if 'AreaShape_' in tick_text:
        new=tick_text.split('_')[1]
    else:
        new=tick_text
    new=new.replace('2','Short')
    new=new.replace('6','Medium') 
    new=new.replace('10','Long')
    tick_label.set_text(new)
    
    
# colorbar.set_ticks([-.99,0,.99])
# colorbar.set_ticklabels(['Decrease','No Change','Increase'])
plt.savefig('FigureSupp_CorrelationMatrix.png')
plt.close(fig)

#get Hierarchical av linkage clustering for 15 clusters, match with feature names and export
n_clusters=15

average_linkage = cluster.AgglomerativeClustering(
    linkage="average", affinity="cityblock",
    n_clusters=n_clusters) 


y_pred = average_linkage.fit_predict(x)   

clustDF=pd.DataFrame({'Feature':corr.columns,'Cluster':y_pred })
clustDF.to_csv('CorrClusteredFeatures.csv')

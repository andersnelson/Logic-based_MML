# -*- coding: utf-8 -*-
"""
Created on Wed May 17 15:50:14 2023

Run PCA on drug screen data and output scores and loadings

@author: arn4va
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#import statannot as sa
import matplotlib.ticker as ticker
import os
import statsmodels.api as sm
import statsmodels
from statsmodels.formula.api import ols
import re
from sklearn.metrics import r2_score
from scipy import stats
import statsmodels.formula.api as smf
from sklearn.linear_model import LinearRegression
from sklearn import metrics
from sklearn.neighbors import KernelDensity as get_kernel
from statsmodels.nonparametric.kde import KDEUnivariate
import math
import time
import scipy.stats as ss
import umap
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import random
from sklearn import cluster, mixture, metrics
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from functools import reduce
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
min_max_scaler=MinMaxScaler()
import statistics as stats
import csv
import os
import math
import scipy.stats as stats
from sklearn import preprocessing
from sklearn import linear_model
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statannotations.Annotator import Annotator
import functools
import scikit_posthocs
import statsmodels.api as sa
import statsmodels.formula.api as sfa
import scikit_posthocs as sp
import statsmodels.stats.multitest as multi
ttest=scikit_posthocs.posthoc_ttest


def unique(list1): 
        # intilize a null list 
        unique_list = [] 
        # traverse for all elements 
        for x in list1: 
            # check if exists in unique_list or not 
            if x not in unique_list: 
                unique_list.append(x) 
        # print list 
        return unique_list

''' Cluster feature correlation matrix'''
dt=pd.read_csv('ProcessedDrugScreen.csv')
dt.index=dt.DG1
dt=dt.drop('DG1',axis=1)
dt=dt.drop(['AreaShape_Orientation'],axis=1)
#dt=dt.drop(['AreaShape_Orientation','plate_number','Metadata_Well'],axis=1)
df=dt
corr=df.corr()
x=corr.iloc[:,1:]
x=x.values


#plot the corr matrix with clustermap to get an idea of cluster number
fig,ax = plt.subplots(1,1,figsize=(18,12))
# my_colors=['#4292c6','#bdbdbd','#cb181d']
sns.clustermap(data=corr,cmap=sns.diverging_palette(220, 20, as_cmap=True),xticklabels=True,yticklabels=True)

# colorbar.set_ticks([-.99,0,.99])
# colorbar.set_ticklabels(['Decrease','No Change','Increase'])
plt.show()


#import clustering algortihms
from sklearn import cluster, mixture

#get Hierarchical av linkage clustering for 15 clusters, match with feature names and export
n_clusters=15

average_linkage = cluster.AgglomerativeClustering(
    linkage="average", affinity="cityblock",
    n_clusters=n_clusters) 


y_pred = average_linkage.fit_predict(x)   

clustDF=pd.DataFrame({'Feature':corr.columns,'Cluster':y_pred })
clustDF.to_csv('CorrClusteredFeatures.csv')



''' Read in retained feature set and reduce using PCA '''

dfh=pd.read_csv('RetainedFeatures_Edited.csv')

keepFeatures=dfh[dfh.Keep=='Yes']
df2=df[[col for col in df.columns if col in keepFeatures.Feature.values]]

x=df2
x=x.values
x = StandardScaler().fit_transform(x)
pca = PCA()
pca.fit_transform(x)

pca = PCA(n_components=6)
principalComponents = pca.fit_transform(x)
pcf = pd.DataFrame(data = principalComponents)
PCA_frame=pcf #store PCA_reduced results
PCA_frame.insert(0,'DG1',df2.index)
pca.explained_variance_ratio_

PCA_frame.to_csv('PCA_Scores.csv')

loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2','PC3','PC4','PC5','PC6'], 
                        index=df2.columns)

loadings.to_csv('PCA_Loadings.csv')
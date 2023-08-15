# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 10:35:40 2021

@author: arn4va

Creates a heatmap of fold changes for drug screen experiment data
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
import statsmodels.stats.multitest as multi
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.graphics.gofplots import qqplot
import statsmodels.api as sa
import statsmodels.formula.api as sfa
import scikit_posthocs as sp

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
    


'''########################  Fold Change-Based heatmaps expanded with all concentrations #######################'''
#Read in drug screen data
df=pd.read_csv('RawScreenData.csv')
df1=df
#Drug names from screen
DrugNames=['Anakinra','Valsartan','BNP','Val_BNP',
           'CWHM12','Glutathione','Marimistat','Salbutamol',
           'Fasudil','Pirfenidone','SB203580','HY12289A','WH4023']
#Plate numbers corresponding to drug names above, plates where each drug was used
plateNumbers=[1,1,2,2,3,3,4,4,5,6,6,7,8]
def most_frequent(List):
    return max(set(List), key = List.count)
dfs=list()
for mets in [['Intensity_IntegratedIntensity_Collagen'],
             ['Intensity_IntegratedIntensity_Actin'],
             ['Intensity_UpperQuartileIntensity_aSMA']]:
    treats=list()
    changes=list()
    pvals=list()
    plate=list()
    metric=mets[0]
    controlSig=list()
    tgfbSig1=list()
    tgfbSig2=list()
    tgfbSig3=list()
    il1Sig1=list()
    il1Sig2=list()
    il1Sig3=list()
    comboSig=list()
    sig_level=0.1 #significance p-value cutoff
    for i in range(0,len(DrugNames)):
        drug=DrugNames[i]
        df2=df1[df1.plate_number == plateNumbers[i]]
        df3=pd.DataFrame(df2.groupby(['DG1','Metadata_Well'])[metric].median().groupby(['DG1']).mean())
        df3.reset_index(inplace=True)
        controlSig.append(df3[df3.DG1==drug+'_2_control'][mets[0]].values[0]/df3[df3.DG1=='control'][mets[0]].values[0])
        tgfbSig1.append(df3[df3.DG1==drug+'_1_TGFB'][mets[0]].values[0]/df3[df3.DG1=='TGFB'][mets[0]].values[0])
        tgfbSig2.append(df3[df3.DG1==drug+'_2_TGFB'][mets[0]].values[0]/df3[df3.DG1=='TGFB'][mets[0]].values[0])
        tgfbSig3.append(df3[df3.DG1==drug+'_3_TGFB'][mets[0]].values[0]/df3[df3.DG1=='TGFB'][mets[0]].values[0])
        il1Sig1.append(df3[df3.DG1==drug+'_1_IL1'][mets[0]].values[0]/df3[df3.DG1=='IL1'][mets[0]].values[0])
        il1Sig2.append(df3[df3.DG1==drug+'_2_IL1'][mets[0]].values[0]/df3[df3.DG1=='IL1'][mets[0]].values[0])
        il1Sig3.append(df3[df3.DG1==drug+'_3_IL1'][mets[0]].values[0]/df3[df3.DG1=='IL1'][mets[0]].values[0])
        comboSig.append(df3[df3.DG1==drug+'_2_TGFB_IL1'][mets[0]].values[0]/df3[df3.DG1=='TGFB_IL1'][mets[0]].values[0])

    Control_Exp=controlSig
    IL1_Exp1=il1Sig1
    IL1_Exp2=il1Sig2
    IL1_Exp3=il1Sig3
    TGFB_Exp1=tgfbSig1
    TGFB_Exp2=tgfbSig2
    TGFB_Exp3=tgfbSig3
    Combined_Exp=comboSig
    
    dfh=pd.DataFrame({'Control_Exp':Control_Exp,
                      'IL1_Exp1':IL1_Exp1,
                      'IL1_Exp2':IL1_Exp2,
                      'IL1_Exp3':IL1_Exp3,
                      'TGFB_Exp1':TGFB_Exp1,
                      'TGFB_Exp2':TGFB_Exp2,
                      'TGFB_Exp3':TGFB_Exp3,
                      'Combined_Exp':Combined_Exp})    
    

    dfh.index=DrugNames
    
    
    dfs.append(dfh)
dfz = {0: dfs[0], 1: dfs[1], 2: dfs[2]}
suffix = ('_Collagen', '_Factin', '_aSMA')
for i in dfz:
    dfz[i] = dfz[i].add_suffix(suffix[i])
import functools    
from functools import reduce
def agg_df(dfList):
    temp=reduce(lambda left, right: pd.merge(left, right, 
                                             left_index=True, right_index=True, 
                                             how='outer'), dfList)
    return temp





df = agg_df(dfz.values())

#plot the Exp data only


exp=df[[col for col in df.columns if 'Exp' in col]]

fig,ax = plt.subplots(1,1,figsize=(25,20))
# my_colors=['#4292c6','#bdbdbd','#cb181d']
#sns.heatmap(data=exp, ax=ax,cmap=sns.diverging_palette(220, 20, as_cmap=True),linewidth=0.1,vmin=0.5,vmax=1.5)
sns.clustermap(data=exp,col_cluster=False,cmap=sns.diverging_palette(220, 20, as_cmap=True),linewidth=0.1,vmin=0.5,vmax=1.5,
               figsize=(25,20))
#ax.set_title('Model-Experiment Comparisons for '+mets[1],fontsize=20)
plt.setp(ax.get_xticklabels(), rotation=90,fontsize=8)
plt.setp(ax.get_yticklabels(), fontsize=8)
#colorbar = ax.collections[0].colorbar
# colorbar.set_ticks([-.99,0,.99])
# colorbar.set_ticklabels(['Decrease','No Change','Increase'])
plt.savefig('Fig1A_Heatmap.png')  
plt.show()
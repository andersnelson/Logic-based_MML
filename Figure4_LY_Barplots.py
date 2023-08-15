# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 13:04:11 2022

Plot data and run stats for PI3K inhibitor experiment in human cardiac fibroblasts

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
import statsmodels.stats.multitest as multi
import scipy
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
    

#read in data file
df=pd.read_csv('ProcessedData_PI3K.csv')
#remove unwanted treatment groups
df=df[~df.TG.str.contains('Fak')]
df=df[~df.TG.str.contains('WH_TGFB')]
df=df[~df.TG.str.contains('SB_TGFB')]


for metric in ['Intensity_IntegratedIntensity_Collagen',
'Intensity_IntegratedIntensity_Actin',
'AngularSecondMoment_Actin_10',
]:
    df1=pd.DataFrame(df.groupby(['TG','Metadata_Well'])[metric].median())
    df1['TG1']=df1.index.get_level_values(0)
    # fig=plt.figure(figsize=(20,15))
    # g=sns.barplot(data=df1,x='TG1',y=metric,color='#7fcdbb',capsize=.2)
    # sns.swarmplot(data=df1,x='TG1',y=metric,color='black',size=10)
    # #set up the stat annottator for pvalues
    # g.set_title(metric, fontsize=30)
    # g.set_xlabel('Well', fontsize=20)
    # g.set_ylabel(metric, fontsize=20)
    # plt.setp(g.get_xticklabels(), fontsize=25,rotation=90)
    # plt.setp(g.get_yticklabels(), fontsize=30)
    
    df2=df1[(df1.TG1 =='control')|(df1.TG1 =='LY')]
    custom_dict = {'control': 0, 'LY': 1}
    df2=df2.sort_values(by=['TG1'], key=lambda x: x.map(custom_dict))
    fig=plt.figure(figsize=(20,15))
    g=sns.barplot(data=df2,x='TG1',y=metric,color='#7fcdbb',capsize=.2)
    #sns.swarmplot(data=df2,x='TG1',y=metric,color='black',size=10)
    #set up the stat annottator for pvalues
    g.set_title(metric, fontsize=30)
    g.set_xlabel('Well', fontsize=20)
    g.set_ylabel(metric, fontsize=20)
    plt.setp(g.get_xticklabels(), fontsize=25,rotation=90)
    plt.setp(g.get_yticklabels(), fontsize=30)
    plt.savefig('Figure4_LY_exp_'+metric+'.png')
    ttest=scipy.stats.ttest_ind(df2[df2.TG1=='control'][metric].values, df2[df2.TG1=='LY'][metric].values, axis=0)
    print(metric)
    print(ttest)



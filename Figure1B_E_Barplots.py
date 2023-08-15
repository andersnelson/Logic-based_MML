# -*- coding: utf-8 -*-
"""
Created on Tue May 31 14:09:04 2022

@author: arn4va
Generates barplots for used in figure 1 looking at model and experiment outputs
Experimental values for pirfenidone and fasudil treatment
Also preforms a TUKEY-HSD test across all comparisons
"""

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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
from sklearn import cluster, mixture
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import scipy.stats as stats
from sklearn import preprocessing
from sklearn import linear_model
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statannotations.Annotator import Annotator
import functools
import scikit_posthocs
import statsmodels.api as sa
import statsmodels.formula.api as sfa
ttest=scikit_posthocs.posthoc_ttest
import statsmodels.stats.multitest as multi
'''
### load in dataframes to be used ###
'''
#Read in the raw screen data
df=pd.read_csv('RawScreenData.csv')

#drugs of interest
drugs=['Fasudil','Pirfenidone','WH4023']
#plate numbers correspoinding to these drugs, used to get control conditions
plateNums=[5,6,8]

plt.rcParams.update({'font.size': 30})
# Loop to make plots for integrated intensities for Collagen, Actin, and aSMA
for i in range(0,len(drugs)):
    for metric in ['Intensity_IntegratedIntensity_Collagen',
               'Intensity_IntegratedIntensity_Actin',
               'Intensity_IntegratedIntensity_aSMA']:
        df1=df[df.plate_number==plateNums[i]]
        df1=df1[df1.DG1.isin(['control',drugs[i]+'_2_control','TGFB',drugs[i]+'_3_TGFB'])]
        custom_dict = {'control': 0, drugs[i]+'_2_control':1, 'TGFB': 2, drugs[i]+'_3_TGFB':3} 
        df1=pd.DataFrame(df1.groupby(['DG1','Metadata_Well'])[metric].median())
        df1.reset_index(inplace=True)
        tukey = pairwise_tukeyhsd(endog=df1[metric],groups=df1['DG1'],alpha=0.1)
        tdf = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
        formatted_pvalues = [f'p={pvalue:.3}' for pvalue in tdf['p-adj']] #allocate formatted p-values for test
        pairs=[]
        for k in range(0,len(tdf)):
            pairs.append((tdf.group1[k],tdf.group2[k]))
        df1=df1.sort_values(by=['DG1'], key=lambda x: x.map(custom_dict))
        err=df1.groupby(['DG1'],sort=False,)[metric].std()/math.sqrt(3) #SEM=stdev/sqrt(n), calculate S.E.M. for errorbars
        yVals=df1.groupby(['DG1'],sort=False,)[metric].mean()
        sns.set(font_scale = 3)
        sns.set_style("white")
        sns.set_style("ticks")
        fig=plt.figure(figsize=(20,15))
        g=sns.barplot(data=df1,x='DG1',y=metric,color='#bdbdbd',capsize=.2,ci=None)
        #add errorbars
        annotator = Annotator(g, pairs,data=df1,x='DG1',y=metric)
        annotator.set_pvalues(tdf['p-adj'].values) #star annotation default?
        annotator.annotate()
        g.errorbar([0,1,2,3], yVals, xerr=0, yerr=err.values,fmt='none',
                   color='k',capsize=10,elinewidth = 3)
        g.set_title('', fontsize=30)
        g.set_xlabel('', fontsize=20)
        g.set_ylabel('', fontsize=20)
       # plt.setp(g.get_xticklabels(), fontsize=30)
        g.set_xticklabels(['Control',drugs[i],'TGF'+r'$\beta$',drugs[i]+'+'+'\nTGF'+r'$\beta$'],size=40)
        plt.setp(g.get_yticklabels(), fontsize=40)
        plt.setp(g.get_yticklabels(), fontsize=30)
        # plt.show()
        plt.savefig('Figure1_Barplot_'+drugs[i]+'_'+metric+'_exp.png')
        plt.close(fig) 

'''
#### Section for WH plots for ASM
'''
#Read in the processed data
df=pd.read_csv('RawScreenData.csv')



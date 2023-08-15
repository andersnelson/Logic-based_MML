# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:13:21 2020


implement ridge regression for integrated F-actin intensity of and long actin angular second moment



@author: arn4va
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#import statannot as SA
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

#function to get unique items from list
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


''' Regression Forward predictions and show improved model '''
from scipy import stats
#load in the sensitivity analysis from Fib617
dfs=pd.read_csv('Fib617_AllInputSens.csv')
dfs=dfs.set_index('specID')
dfs=dfs.T
#load in the influence analysis from Fib617
dfi=pd.read_csv('Fib617_AllInputInfluence.csv')
dfi=dfi.set_index('specID')
dfi=dfi.T
#this output is used to train the regression model
dfm=pd.read_csv('SimulatedDrugScreen.csv')
#read in drug screen data from Experiment
df=pd.read_csv('ProcessedDrugScreen.csv')
df['Treatment']=df['DG1']
df=df[df.DG1.isin(dfm.Treatment.values)]
#cluster on combined influence and sensitivity analysis
for k in range(11,12):
    x1=dfs.values
    x1 = StandardScaler().fit_transform(x1).T
    x2=dfi.values
    x2 = StandardScaler().fit_transform(x2).T
    x=np.concatenate([x1,x2],axis=1)
    k_means=cluster.KMeans(n_clusters=k,random_state=1,n_init=20)
    y_pred = k_means.fit_predict(x)
sensClusters=y_pred

### Regression models for F-actin and angular second moment
mets=['Intensity_IntegratedIntensity_Actin',
  'AngularSecondMoment_Actin_10']
coefs=list()
for i in range(0,len(mets)): 
    #this output is used to train the regression model
    dfm=pd.read_csv('SimulatedDrugScreen.csv')
    #read in drug screen data from Experiment
    df=pd.read_csv('ProcessedDrugScreen.csv')
    df['Treatment']=df['DG1']
    dfm=dfm[~dfm.Treatment.str.contains('Gal')]
    if i == 0: #If statement to remove pirfenidone if Actin is treatment group
        dfm=dfm[~dfm.Treatment.str.contains('Pirf')]
    dfm=dfm.set_index('Treatment')    
    #run the regression
    col=mets[i]
    df1=df[['Treatment',col]]
    df1=df1.set_index('Treatment')
    df2=dfm.join(df1,on='Treatment')
    y_pred=sensClusters
    df4=df2
    Y=df4[col].values
    Y=min_max_scaler.fit_transform(Y.reshape(-1, 1))
    df3=dfm
    df3[df3.columns] = min_max_scaler.fit_transform(df3[df3.columns]) # set values across columns between 0-1
    df3=df3.T
    df3['Treatment']=df3.index
    df3['Module']=y_pred
    df3=df3.dropna()
    df3=df3.drop(['Treatment'],axis=1)
    df3.index=df3.Module.values
    df3=pd.DataFrame(df3.groupby(by='Module').mean())
    df3=df3.T
    X=df3
    X[X.columns]=min_max_scaler.fit_transform(X[X.columns])
    reg = linear_model.RidgeCV(alphas=(0.5),
                               cv=None,store_cv_values=True)
    reg.fit(X,Y)
    reg.score(X,Y)
    
    ## extension to Sklearn regresion to get paramas and p_Values for coefficients ###
    
    params = np.append(reg.intercept_,reg.coef_) #get the intercept and parameter values from
    predictions = reg.predict(X)
    
    newX = np.append(np.ones((len(X),1)), X, axis=1)
    MSE = (sum((Y-predictions)**2))/(len(newX)-len(newX[0]))
    
    var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
    sd_b = np.sqrt(var_b)
    ts_b = params/ sd_b
    
    p_values =[2*(1-stats.t.cdf(np.abs(i),(len(newX)-len(newX[0])))) for i in ts_b]
    
    sd_b = np.round(sd_b,3)
    ts_b = np.round(ts_b,3)
    p_values = np.round(p_values,3)
    params = np.round(params,4)
    
    myDF3 = pd.DataFrame()
    myDF3["Coefficients"],myDF3["Standard Errors"],myDF3["t values"],myDF3["Probabilities"] = [params,sd_b,ts_b,p_values]
    print(myDF3)

    vdf=myDF3.loc[1:]
    vdf['Modules']=['Smad3','Mech','Ras_Raf','Autocrine',
                   'PI3K','STAT','P38_Calcium','PKA',
                   'PDGF','MKK3','NP']
    
    def negLogConvert(row):
        val=row['Probabilities']
        val_new=-math.log(val)
        return val_new

    vdf['negLogP']=vdf.apply(lambda row: negLogConvert(row), axis=1)
     
    print('Regression LOOCV mean MSE =')
    print(reg.cv_values_.mean())
    coefs.append(reg.coef_[0])
    cf=pd.DataFrame({'Node':X.columns,
    mets[i]:reg.coef_[0]})
    
    # #run cross validation
    # from sklearn.model_selection import cross_val_score
    
    fig=plt.figure()
    fig.set_size_inches(20,18)
    cf=cf.sort_values(by='Node')
    g=sns.barplot(data=cf,y='Node',x=mets[i],
                  orient = 'h',color='#636363')
    g.set_title(mets[i],fontsize=24)
    g.set_ylabel('',fontsize=20)
    g.set_yticklabels(['Smad3','Mech','Ras_Raf','Autocrine',
                   'PI3K','STAT','P38_Calcium','PKA',
                   'PDGF','MKK3','NP'])
    plt.setp(g.get_xticklabels(), rotation=0,fontsize=40)
    plt.setp(g.get_yticklabels(), rotation=0,fontsize=40)
    plt.tight_layout() #prevents bottom of plot from being cutoff
    g.set_xlabel('',fontsize=20)
    
        # plt.show()
    plt.savefig('Fig2_4_RidgeCoefficients'+mets[i]+'.pdf')
    plt.savefig('Fig2_4_RidgeCoefficients'+mets[i]+'.png')
    plt.close(fig) 
    
       



''' F-ACTIN REGRESSION IMPROVES MODEL PREDICTIONS'''
### Train the model
dfm=pd.read_csv('SimulatedDrugScreen.csv')
#read in drug screen data from Experiment
df=pd.read_csv('ProcessedDrugScreen.csv')
df['Treatment']=df['DG1']
dfm=dfm[~dfm.Treatment.str.contains('Gal')]
dfm=dfm[~dfm.Treatment.str.contains('Pirf')]#This Line removes the pirfenidone data
dfm=dfm.set_index('Treatment')    

metric='Intensity_IntegratedIntensity_Actin'
#Model node corresponding to Exp data
mNode='Factin'
#Drug of interest
drugChoice='Pirfenidone'
stims=['Pirfenidone_3_TGFB','TGFB','control','Pirfenidone_2_control']
''' Training the model on single simulation of mean '''
col=metric
df=df[~df.DG1.str.contains('Pirf')] #This Line removes the pirfenidone data
df=df.reset_index()
df1=df[['Treatment',col]]
df1=df1.set_index('Treatment')
dfm1=dfm
df2=dfm1.join(df1,on='Treatment')
y_pred=sensClusters
df4=df2
Y=df4[col].values
Yscale=Y
Y=min_max_scaler.fit_transform(Y.reshape(-1, 1))
df3=dfm
df3[df3.columns] = min_max_scaler.fit_transform(df3[df3.columns]) # set values across columns between 0-1
df3=df3.T
df3['Treatment']=df3.index
df3['Module']=y_pred
df3=df3.dropna()
df3=df3.drop(['Treatment'],axis=1)
df3.index=df3.Module.values
df3=pd.DataFrame(df3.groupby(by='Module').mean())
df3=df3.T
X=df3
X[X.columns]=min_max_scaler.fit_transform(X[X.columns])
reg = linear_model.RidgeCV(alphas=0.5)
reg.fit(X,Y)

''' Test the model on Ensmeble simualtions '''

dfm=pd.read_csv('SimulatedDrugScreen_Ensembles.csv')
dfm=dfm[~dfm.Treatment.str.contains('Gal')]
dfm=dfm.set_index('Treatment')   
dfsNew=dfm
#### Get Prediciton Matrix for Metric ####
dfs0=dfsNew  
dfs0=dfs0.drop('Index',axis=1)

dfs0[dfs0.columns]=min_max_scaler.fit_transform(dfs0[dfs0.columns]) #normalize valuea cross full simulation set

dfs1=dfs0

drugs=dfs1.index.values

yPred=list()
dfs1['Ensemble']=np.repeat(np.linspace(1,100,num=100),56)
for k in np.linspace(1,100,num=100):
    dfs3=dfs1[dfs1.Ensemble==k]
    dfs3=dfs3.drop(columns=['Ensemble'])
    dfs2=dfs3.T
    dfs2['cluster']=y_pred
    xPred=dfs2.groupby('cluster').mean()
    xPred=xPred.T
    xPred[xPred.columns]=min_max_scaler.fit_transform(xPred[xPred.columns])
    xPred=xPred.T
    for node in xPred.columns:
        yPred.append(reg.predict(xPred[node].values.reshape(1, -1))[0][0])


#use ensmble for predictions

df_cell=pd.read_csv('DrugScreen_Well_Level.csv')

yf=pd.DataFrame({'Treatment':drugs,'Predicted Value':yPred})
yf1=yf[yf.Treatment.isin(stims)]
dfmSimsAll=dfm=pd.read_csv('SimulatedDrugScreen_Ensembles.csv')
dfmSimsAll=dfmSimsAll.set_index('Treatment')
dfm5=dfmSimsAll[dfmSimsAll.index.isin(stims)]
df5=df_cell[(df_cell.DG1.isin(stims)) & (df_cell.plate_number==6)]
yf1=yf1.sort_values(by='Treatment',ascending=False)
df5=df5.sort_values(by='DG1',ascending=False)

# def rescalePred(row):
#     val=row['Predicted Value']
#     val_new=val*(max(Yscale)-min(Yscale))+min(Yscale)
#     return val_new

def rescalePred(row):
    val=row['Predicted Value']
    val_new=val*(max(Yscale)-min(Yscale))+min(Yscale)
    return val_new

yf1['YpredScaled']=yf1.apply(lambda row: rescalePred(row), axis=1)

custom_dict = {'control': 0, 'TGFB': 1, drugChoice+'_2_control': 2,drugChoice+'_3_TGFB':3}
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20,12))
fig.suptitle('Model Comparisons for F-actin')
# ################################# Simulation ##########################
dfm5['DG1']=dfm5.index
dfm5=dfm5.sort_values(by=['DG1'], key=lambda x: x.map(custom_dict))
#get values for custom errorbars to add to plots
err=dfm5.groupby(['DG1'],sort=False,)[mNode].std()/math.sqrt(3) #SEM=stdev/sqrt(n)
yVals=dfm5.groupby(['DG1'],sort=False,)[mNode].mean()
#get tukey HSD values
tukey = pairwise_tukeyhsd(endog=dfm5[mNode],groups=dfm5['DG1'],alpha=0.5)
tdf = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
formatted_pvalues = [f'p={pvalue:.3}' for pvalue in tdf['p-adj']] #allocate formatted p-values for test
pairs=[]
for k in range(0,len(tdf)):
    pairs.append((tdf.group1[k],tdf.group2[k]))
sns.barplot(ax=ax1,data=dfm5, x='DG1', y=mNode, color='#c6dbef', ci=None)
annotator = Annotator(ax1, pairs, data=dfm5, x='DG1', y=mNode)
#annotator.set_custom_annotations(formatted_pvalues)
annotator.set_pvalues(tdf['p-adj'].values) #star annotation default?
annotator.annotate()
ax1.errorbar([0,1,2,3], yVals, xerr=0, yerr=err.values,fmt='none',
                  color='k',capsize=10,elinewidth = 3)
ax1.set_title('Fibroblast Model')

################################ Ridge Regression ###########################
#sns.barplot(ax=ax2,data=yf1, x='Treatment', y='Predicted Value')
yf1=yf1.sort_values(by=['Treatment'], key=lambda x: x.map(custom_dict))
err=yf1.groupby(['Treatment'],sort=False,)['Predicted Value'].std()/math.sqrt(3) #SEM=stdev/sqrt(n)
yVals=yf1.groupby(['Treatment'],sort=False,)['Predicted Value'].mean()
tukey = pairwise_tukeyhsd(endog=yf1['Predicted Value'],groups=yf1['Treatment'],alpha=0.5)
tdf = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
formatted_pvalues = [f'p={pvalue:.3}' for pvalue in tdf['p-adj']] #allocate formatted p-values for test
pairs=[]
for k in range(0,len(tdf)):
    pairs.append((tdf.group1[k],tdf.group2[k]))
sns.barplot(ax=ax3,data=yf1, x='Treatment', y='Predicted Value', color='#4292c6', ci=None)
annotator = Annotator(ax3, pairs, data=yf1, x='Treatment', y='Predicted Value')
#annotator.set_custom_annotations(formatted_pvalues)
annotator.set_pvalues(tdf['p-adj'].values) #star annotation default?
annotator.annotate()
ax3.errorbar([0,1,2,3], yVals, xerr=0, yerr=err.values,fmt='none',
                  color='k',capsize=10,elinewidth = 3)
ax3.set_title('Fibroblast Model + Regression')
#ax2.set_ylim([0.2, 0.4])
# ###############################  Experiment #################################
df5=df5.sort_values(by=['DG1'], key=lambda x: x.map(custom_dict))
#get values for custom errorbars to add to plots
err=df5.groupby(['DG1'],sort=False,)[metric].std()/math.sqrt(3) #SEM=stdev/sqrt(n)
yVals=df5.groupby(['DG1'],sort=False,)[metric].mean()
tukey = pairwise_tukeyhsd(endog=df5[metric],groups=df5['DG1'],alpha=0.5)
tdf = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
formatted_pvalues = [f'p={pvalue:.3}' for pvalue in tdf['p-adj']] #allocate formatted p-values for test
pairs=[]
for k in range(0,len(tdf)):
    pairs.append((tdf.group1[k],tdf.group2[k]))
sns.barplot(ax=ax2,data=df5, x='DG1', y=metric, color='#6baed6', ci=None)
annotator = Annotator(ax2, pairs, data=df5, x='DG1', y=metric)
annotator.set_pvalues(tdf['p-adj'].values) #star annotation default?
annotator.annotate()
ax2.errorbar([0,1,2,3], yVals, xerr=0, yerr=err.values,fmt='none',
                  color='k',capsize=10,elinewidth = 3)
ax2.set_title('Fibroblast Model + Regression')
ax2.set_title('Experiment')

plt.setp(ax1.get_xticklabels(), rotation=90,fontsize=18)
plt.setp(ax2.get_xticklabels(), rotation=90,fontsize=18)
plt.setp(ax3.get_xticklabels(), rotation=90,fontsize=18)

plt.tight_layout() #prevents bottom of plot from being cutoff
plt.savefig('Figure2_Actin_Regression_Bars_noPirfenidone.png')
plt.savefig('Figure2_Actin_Regression_Bars_noPirfenidone.pdf')
plt.close(fig)  
















''' KNOCKDOWN SNSITIVITY ANALYSES FOR RIDGE REGRESSIONS MODELS '''

import statistics as stats
mets=['Intensity_IntegratedIntensity_Actin',
      'AngularSecondMoment_Actin_10']
perturbs=['Pirfenidone_3_TGFB',
          'WH4023_3_TGFB']
#baselines=['TGFB','TGFB']
#baseline
drugLabs=['Pirfenidone','WH4023']
#metLabs=['F-actin', 'Actin Angular Second Moment (Long)']
metLabs=mets
for i in range(0,len(mets)): 
    df=pd.read_csv('ProcessedDrugScreen.csv')
    df['Treatment']=df.DG1
    dfm=pd.read_csv('SimulatedDrugScreen.csv')
    dfs=pd.read_csv('DrugSimulations_KDsens.csv')
    dfm=dfm[~dfm.Treatment.str.contains('Gal')]
    if i == 0: #If statement to remove pirfenidone if Actin is treatment group
        dfm=dfm[~dfm.Treatment.str.contains('Pirf')]
    dfm=dfm.set_index('Treatment')    
#### Determine which nodes change in reponse to drug perturbation #####
    metric=mets[i]
    drug=perturbs[i]
    #baseline=baselines[i]
    baseline='TGFB'
    
    #### Determine which nodes impact predicted Metric when KD
    ### fit Model
    
    #### Fit the model ####     
    col=metric
    df1=df[['Treatment',col]]
    df1=df1.set_index('Treatment')
    df2=dfm.join(df1,on='Treatment')
    y_pred=sensClusters
    df4=df2
    Y=df4[col].values
    Y=min_max_scaler.fit_transform(Y.reshape(-1, 1))
    df3=dfm
    df3[df3.columns] = min_max_scaler.fit_transform(df3[df3.columns]) # set values across columns between 0-1
    df3=df3.T
    df3['Treatment']=df3.index
    df3['Module']=y_pred
    df3=df3.dropna()
    df3=df3.drop(['Treatment'],axis=1)
    df3.index=df3.Module.values
    df3=pd.DataFrame(df3.groupby(by='Module').mean())
    df3=df3.T
    X=df3
    X[X.columns]=min_max_scaler.fit_transform(X[X.columns])
    reg = linear_model.Ridge(alpha=0.5)
    reg.fit(X,Y)
    #highCoefs=np.argsort(abs(reg.coef_))[0][8:11] #take the top 3 coefficients by absolute values
    sDevCutoff=1.645*(stats.stdev(abs(reg.coef_[0]))) #zcore 1.645 corresponds to p=0.95 in one tailed distribution
    highCoefs=np.where(abs(reg.coef_[0])>=sDevCutoff)[0] # take clusters with coefficients above abs value threshold
    #highCoefs=np.where(abs(reg.coef_[0])>=.2)[0]
    #stats.stdev(abs(reg.coef_[0]))
    #### Get Prediciton Matrix for Metric ####
    dfs0=dfs  
    dfs0=dfs0.drop('Drug',axis=1)
    dfs0=dfs0.drop('Node',axis=1)
    dfs0[dfs0.columns]=min_max_scaler.fit_transform(dfs0[dfs0.columns]) #normalize valuea cross full simulation set
    dfs0['Drug']=dfs.Drug.values
    dfs0['Node']=dfs.Node.values
    #dfs1=dfs0[dfs0.Node=='Control']
    dfs1=dfs0
    # drugs=['control','TGFB','IL1','TGFB_IL1','WH4023_3_TGFB',
    #        'Glutathione_3_TGFB','HY12289A_2_control','Glutathione_2_TGFB_IL1',
    #        'Fasudil_3_TGFB','Fasudil_2_control','HY12289A_3_TGFB','Salbutamol_2_TGFB_IL1']
    drugs=dfs1.Drug.values
    nodes=dfs1.Node.values
    dfs1=dfs1[dfs1.Drug.isin(drugs)]
    dfs2=dfs1
    dfs2=dfs2.drop('Drug',axis=1)
    dfs2=dfs2.drop('Node',axis=1)
    dfs2=dfs2.T
    dfs2['cluster']=y_pred
    xPred=dfs2.groupby('cluster').mean()
    xPred=xPred.T
    xPred[xPred.columns]=min_max_scaler.fit_transform(xPred[xPred.columns])
    yPred=list()
    xPred=xPred.T
    for node in xPred.columns:
        yPred.append(reg.predict(xPred[node].values.reshape(1, -1))[0][0])
        
    yf=pd.DataFrame({'Treatment':drugs,'Predicted Value':yPred,
                     'Node':nodes})    
    
    yf2=yf[yf.Treatment==drug]
    sc2=sensClusters
    sc2=np.insert(sc2,0,11)
    yf2['sensClusters']=sc2
    
    yf2=yf2[(yf2.sensClusters.isin(highCoefs)) | (yf2.sensClusters==11)] #filtering for nodes only in the top 3 clusters
   #list of nodes to remove
    rem_list=['inputMechanical','mechanical','B1int','MEKK1',
              'Factin','TNFaR','TNFa','inputTNFa','LOX','PI3K',
              'inputTGFB','TGFB','TGFB1R','CI','CIII',
              'periostin','contraction','CImRNA','CIIImRNA',
              'aSMA','xlinkFibers','PAI1','AngIImRNA','ET1','inputET1','ETAR',
             'CTGF','NFAT','calcineurin']
    yf2=yf2[~yf2.Node.isin(rem_list)]
    yf2=yf2.sort_values(by='Predicted Value')
    def Delta(row):
        ctrl=yf2[yf2.Node=='Control']['Predicted Value']
        val=row['Predicted Value']-ctrl
        return val
    yf2['Delta']=yf2.apply(lambda row: Delta(row), axis=1)
    yf2=yf2[yf2.Node!='Control']
    
    fig=plt.figure(figsize=(25,20))
    #fig.set_size_inches(15,25)
    g=sns.barplot(data=yf2,x='Node',y='Delta',hue='sensClusters',
                  dodge=False,palette=sns.color_palette(['#fb6a4a','#238443','#4292c6']) )
    plt.tight_layout() #prevents bottom of plot from being cutoff
    g.set_ylabel('',fontsize=24)
    g.set_xlabel('',fontsize=24)
    g.set_title(drugLabs[i]+'+'+mets[i],fontsize=24)
    plt.setp(g.get_xticklabels(), rotation=90,fontsize=50)
    plt.setp(g.get_yticklabels(), rotation=0,fontsize=45)
    g.legend_.remove()
    plt.tight_layout()
   # if i ==0:
        #lt.ylim(.4,.7)
    plt.subplots_adjust(left=0.28)
    plt.savefig('Figure_2_4_'+drug+'_'+metric+'_KD_Barplots.png')
    plt.savefig('Figure_2_4_'+drug+'_'+metric+'_KD_Barplots.pdf')
    plt.close(fig)    







# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:09:09 2020

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


#dfm = masked final filtered
#dfn = nuclei 

''' 
~~~~~~~~~~~~~~~~ READ IN AND PROCESS DATA FROM DRUG FILES ~~~~~~~~~~~
'''

   
meta_file='Y:\\Anders\\Experiments\\6-22-2022\\MetaData layout 6_27_2022.csv' #metadata plate file
#data_file='Y:\\Anders\\Experiments\\6_30_2022\\CellProfilerOutput_6_30_22_bothPlates\\Run_2_7_12_2022_test3_fullTexturesa_Cell.csv'
#data_file='Y:\\Anders\\Experiments\\6_30_2022\\CellProfilerOutput_6_30_22_bothPlates\\9_22_pipelineTestOutput_a_Cell.csv'
data_file='Y:\\Anders\\Experiments\\6_30_2022\\CellProfilerOutput_6_30_22_bothPlates\\9_23_pipelineTestOutput_unprocessedImagesa_Cell.csv' #data from cellProfiler
df=pd.read_csv(data_file)
meta=pd.read_csv(meta_file)

def get_TG (row):
    cols=row['Metadata_WellColumn']
    rows=row['Metadata_WellRow']
    treatment_group=meta.values[rows-2,cols-1]
    return treatment_group

df['TG']=df.apply (lambda row: get_TG(row),axis=1) #treatment group
df=df[~df['TG'].str.contains('Filler')] # remove the filler wells (PBS)


'''
~~~~~~~~~~~~ Get texture angle metrics for individual cells ~~~~~~~~~~~~~~~~`
'''

df2=df
# get list 'abb_names'
#col_select=[col for col in df.columns if 'Texture' in col] 
col_select=[col for col in df.columns if 'Angular' in col] 
# get list 'abb_names'
abb_names=list()
for k in range(0,len(col_select)):
    res = [i.start() for i in re.finditer('_', col_select[k])]
    abb_names.append(col_select[k][res[0]+1:res[2]])
#######
    
#max difference between numbers in a list
def maxDiff(dd):
    return abs(max(dd)-min(dd))

#big function to do everything    
#Subtracts smallest angle from largets angle value ot get final texture value for a given metric
def getTextureMetrics(row,abb_name):
    #tic=time.time()
    d2=row[[col for col in df.columns if abb_name+'_2' in col]].values #offset distance 2
    maxOffset2=maxDiff(d2)
    d6=row[[col for col in df.columns if abb_name+'_6' in col]].values #offset distance 6
    maxOffset6=maxDiff(d6)
    d10=row[[col for col in df.columns if abb_name+'_10' in col]].values #offset distance 10
    maxOffset10=maxDiff(d10)
    #print(tic-time.time())
    return maxOffset2, maxOffset6, maxOffset10
    
#Apply function across the matrix
for abb_name in unique(abb_names):
    tic=time.time()
    exec('df2["'+abb_name+'_2"]=df2.apply (lambda row: getTextureMetrics(row,abb_name)[0],axis=1)' )
    exec('df2["'+abb_name+'_6"]=df2.apply (lambda row: getTextureMetrics(row,abb_name)[1],axis=1)' )
    exec('df2["'+abb_name+'_10"]=df2.apply (lambda row: getTextureMetrics(row,abb_name)[2],axis=1)' )
    #exec('(df2["'+abb_name+'_2"] ,df2["'+abb_name+'_6"], df2["'+abb_name+'_10"])=df2.apply (lambda row: getTextureMetrics(row,abb_name),axis=1)' )
   
    print(time.time()-tic)

#df2.to_csv('Y:\Anders\Operreta Plate Analysis With Steve\\AngleTextureMatrix_7_13_2022.csv')
df2.to_csv('Y:\Anders\Operreta Plate Analysis With Steve\\AngleTextureMatrix_9_23_2022_test2.csv')


#rf=pd.read_csv('Y:\Anders\Operreta Plate Analysis With Steve\\AngleTextureMatrix_7_13_2022.csv')
rf=pd.read_csv('Y:\Anders\Operreta Plate Analysis With Steve\\AngleTextureMatrix_9_23_2022_test2.csv')
'''
###### reduced the number of dataframe columns to keep only the columns we need###
'''
colSelect=list()

for abb_name in unique(abb_names):
    colSelect.append(abb_name+'_2')
    colSelect.append(abb_name+'_6')
    colSelect.append(abb_name+'_10')

intense= [col for col in rf.columns if ('Intensity' in col) & ('Dapi' not in col) & ('Location' not in col) ]
areas= [col for col in rf.columns if ('AreaShape' in col) &  ('Center' not in col) & ('Euler' not in col)]
colSelect.extend(intense)
colSelect.extend(areas)
#colSelect.extend(['AreaShape_Area'])
colSelect.extend(['TG'])
colSelect.extend(['Metadata_Plate_Number'])

#remove orientation, its a nothing feature
rf1=rf[colSelect]
rf1=rf1.drop(['Metadata_Plate_Number'],axis=1)
#rf1.to_csv('Y:\Anders\Operreta Plate Analysis With Steve\\Pi3kExp_7_13_2022_texturesProcessed_cell_level.csv') #th
rf1.to_csv('Y:\Anders\Operreta Plate Analysis With Steve\\Pi3kExp_9_23_2022_test2_texturesProcessed_cell_level.csv')

rep_df=rf.groupby(['Metadata_Well','TG'])[colSelect].median()
rep_df=rep_df.drop(['Metadata_Plate_Number'],axis=1)
#rep_df.to_csv('Y:\Anders\Operreta Plate Analysis With Steve\\Pi3kExp_7_13_2022_texturesProcessed_well_level.csv')
rep_df.to_csv('Y:\Anders\Operreta Plate Analysis With Steve\\Pi3kExp_9_23_2022_test2_texturesProcessed_well_level.csv')

mf=rf.groupby(['Metadata_Plate_Number','Metadata_Well','TG'])[colSelect].median().groupby('TG').mean()
mf=mf.drop(['Metadata_Plate_Number'],axis=1)
#mf.to_csv('Y:\Anders\Operreta Plate Analysis With Steve\\Pi3kExp_7_13_2022_texturesProcessed_treatment_level.csv')
mf.to_csv('Y:\Anders\Operreta Plate Analysis With Steve\\Pi3kExp_9_23_2022_test2_texturesProcessed_treatment_level.csv')




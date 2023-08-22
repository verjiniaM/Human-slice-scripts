#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 15:11:58 2023

@author: rosie
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#%%load data

df = pd.read_excel('/Volumes/TOSHIBA EXT/Verjinia/Data/all_data_human 2.xlsx',
                   sheet_name='adult_repatch_all_QC')
col_drop = ['Unnamed: 0', 'area', 'patcher', 'filename', 'cell_ch', 'cell_ID',
            'repatch', 'comments', 'Unnamed: 0.1', 'Unnamed: 0.1.1', 'Unnamed: 25']

df.drop(columns=col_drop, inplace=True)

controls=df[df['treatment']=='Ctrl']
hiK= df.loc[df['treatment']=='high K']
#%%
controls['cell_ID_new'] = controls['cell_ID_new'].astype(str)
controls=df.loc[df['treatment']=='Ctrl']

#make list of ones for size - to set line thickness same for all
plt.figure(figsize=(2.5,5))
sns.lineplot(x='day', y='TH', data=controls, palette='viridis',hue='patient_age', size='cell_ID_new', sizes=list(np.ones(24)+1),legend=False, ci=None)
sns.despine()
#plt.gca().set_title('Threshold')
plt.gca().set_ylabel('TH, mV')
plt.gca().set_ylim([-50,-25])
plt.gca().set_xlabel('')

#%% Make Dfs with normalised values: Day 2 values normalised to Day1

columns=['Rin', 'resting_potential', 'max_spikes', 'Rheobase', 'AP_heigth', 'TH',
         'max_depol', 'max_repol', 'membra_time_constant_tau', 'capacitance']
#hiK
normd1HK=hiK.copy()
normd2HK=hiK.copy()

for col in columns:
    normd1HK[col] = hiK.loc[hiK['day']=='D1'][col]/hiK.loc[hiK['day']=='D1'][col]
    
for col in columns:
    normd2HK[col] = hiK.loc[hiK['day']=='D2'][col]/hiK.loc[hiK['day']=='D1'][col].values
    
#control
normd1=control.copy()
normd2=control.copy()

for col in columns:
    normd1[col] = controls.loc[controls['day']=='D1'][col]/controls.loc[controls['day']=='D1'][col]
    
for col in columns:
    normd2[col] = controls.loc[controls['day']=='D2'][col]/controls.loc[controls['day']=='D1'][col].values

#combined DF for day 2 normalised values
normd2_both = pd.merge(normd2,normd2HK, how='outer')

#%%
palette={'Ctrl':'w', 'high K': 'salmon'}

plt.figure(figsize=(1.5,5))
sns.violinplot(normd2_both['treatment'],normd2_both['resting_potential'], orient='vertical', hue=normd2_both['treatment'],  legend=False,palette=palette)
#plt.gca().set_title('Vm')
plt.legend([],[], frameon=False)
sns.despine(bottom=True)
plt.gca().tick_params(bottom=False)
plt.gca().set_ylabel('Normalised change, D2/D1')
plt.gca().set_xlabel('')
plt.gca().hlines(y=1, xmin=-0.5, xmax=1.5, linestyle='dotted')    
plt.gca().set_ylim([0,2])
plt.gca().set_yticks([0, 0.5,1,1.5,2])
#%%
def plot_paired_param(df, param1, param2, x2=2):
    if x2==2:
        x2 = np.repeat(x2, len(df))
    plt.rcParams.update({'font.size': 18, 'font.sans-serif': 'Arial'})
    fig, ax = plt.subplots(1,1,figsize=(5,5), dpi= 80)
    for p1, p2, p3, c in zip(df[param1], df[param2], x2, df['Color']):
        plt.plot([1,p3], [p1,p2],marker='o', markersize=8,c=c)
        ax.set_xticks([1,12,18,24])
        ax.set_xticklabels(['0', '12', '18', '24'])
        sns.despine()
    ylims=plt.gca().get_ylim()

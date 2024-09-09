import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import pyplot
import pyabf
import numpy as np
import os
import glob
plt.style.use('./style_plot_intrinsic.mplstyle')

def plot_trace(fn, sweeps_analyse, channel, scaling = 1.0, first_point = 0, last_point = 5500, unit = 'pA'):
    ''' 
    arguemnts : fn - filename, sweep - sweep number, channel - active channel from 1 to 8
    fast visualization of recordings and corresponding protocol
    accepts only 1 int for channel
    '''

    end_fn = fn.rfind('/') + 1
    abf_file = pyabf.ABF(fn)

    if len(abf_file.channelList) < 8:
        if '_Ipatch' in abf_file.adcNames:
            abf_file.adcNames[abf_file.adcNames.index('_Ipatch')] = 'Ch1'
        if 'IN0' in abf_file.adcNames:
            abf_file.adcNames[abf_file.adcNames.index('IN0')] = 'Ch1'
        if 'IN 0' in abf_file.adcNames:
            abf_file.adcNames[abf_file.adcNames.index('IN0')] = 'Ch1'
        channel_name = 'Ch' + str(channel)
        channel = abf_file.channelList[abf_file.adcNames.index(channel_name)]
    else:
        channel_name = 'Ch' + str(channel)
        channel = channel - 1
    
    abf_file.setSweep(0)
    swp_len = len(abf_file.sweepY)
    data_long = abf_file.data[channel] * scaling
    swp_num = int(len(data_long)/swp_len)

    reshape_data = data_long.reshape(swp_num, swp_len)
    
    fig, ax = plt.subplots(1,1, sharex = False, figsize = (14,8))
    if sweeps_analyse != 'all':
        sweeps_delete = list(set(list(range(swp_num))) - set([int(a) for a in sweeps_analyse[1:-1].split(',')]))
        reshape_data = np.delete(reshape_data, sweeps_delete, 0)
    else:
        reshape_data = np.delete(reshape_data, sweeps_delete, 0)

    reshape_data = np.delete(reshape_data, list(range(first_point, last_point)), 1)
    data_long = reshape_data.flatten()
    data_unit = unit if unit is not None else abf_file.adcUnits[channel]

    
    ax.plot(data_long)
    ax.set_ylabel(data_unit)

    #channel_name = 'Ch' + str(channel)
    ax.set_title('{0}, channel {1}'.format(fn[end_fn:],channel_name))
    
    return data_long, ax
    

meta_ = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_spontan/metadata_tables/2024-06-05spontan_meta_BACKUP_COMPLET.xlsx'
df_meta = pd.read_excel(meta_)
df_meta = df_meta.sort_values(['Name of recording'])
file_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_spontan/recordings/'
#path_indx = 51 # 43 no /Volumes

# pick random files from the results dir
data_files_all = sorted(os.listdir('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_spontan/output/cell_props/'))
files_check = [s for s in data_files_all if 'individual' in s]

medians, fns = [], []
for i in range(len(files_check[:10])):
    fn = files_check[i][:files_check[i].find('ch') - 1]
    chan_str = files_check[i][files_check[i].find('ch') + 2: files_check[i].find('ch') + 3]

    df_indv = pd.read_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_spontan/output/cell_props/' + fn + '_ch'+ chan_str +'_individual.csv')
    #avgs_ = pd.read_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_spontan/output/cell_props/' + fn + '_ch'+ chan_str +'_avgs.csv')

    file_path = file_dir + fn + '.abf'
    chan = int(chan_str) + 1

    swp_analysis = df_meta['swps_to_analyse'][(df_meta['Name of recording'] ==  fn + '.abf') \
    &  (df_meta['Channels to use'] == chan)].values
    if len(swp_analysis) == 0:
        print('problem with ' + fn + chan_str)
        continue
    else:
        swp_analysis = swp_analysis[0]

    %matplotlib qt  
    trace, ax = plot_trace(file_path, swp_analysis, chan)
    # fns.append(fn)
    # medians.append(np.median(trace))
    x = df_indv.iloc[0,:][1:].values
    y = [trace[int(a)] for a in x]
    #for col in df_indv.columns[1:]:
    ax.scatter(x,y, zorder = 10, c = 'k')
    #plt.show()
    #%matplotlib
    #input()
    
    
#%%
#%%

import os 
import pandas as pd
import numpy as np
import human_characterisation_functions as hchf
import connection_parameters as con_param
import fine_screen_diff_freq_functions as fine

work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/'
OP =  'OP220615' #input("OP to analyse") + '/'

file_list = sorted(os.listdir(work_dir+OP+'/'))

filenames = []
for file in range(len(file_list)):
    #pclamp files
    if file_list[file][-4:] == '.abf': 
        filenames.append(file_list[file])
    elif (file_list[file][-5:] == '.xlsx' and file_list[file][:2] != '~$'): 
        df_rec = pd.read_excel(work_dir + OP+'/' + file_list[file], header = 1)
        slice_indx = df_rec.index[df_rec['slice'].notnull()]
        slice_names = df_rec['slice'][slice_indx].tolist()

new_slice_names = []
for i in range(len(slice_names)):
    if i < len(slice_names)-1:
        new_slice_names.append([slice_names[i]]*(slice_indx[i+1]-slice_indx[i]))
    else :
        new_slice_names.append([slice_names[i]]*(15))
slice_names  = [x for xs in new_slice_names for x in xs]

con_data = pd.DataFrame(columns = ['fn', 'slice', 'chan_pre', 'chan_post', 'Vm pre', 'Vm post', 
    'Amp1',	'Amp2',	'Amp3',	'Amp4',	'Lat1',	'Lat2',	'Lat3',	'Lat4', 'num excluded swps', 'comments'])

#%%
file_num = 6

fine_con_screen = work_dir + OP + '/' + filenames[file_num-1] 

pre_cell = 2
post_cell = 4
slic = slice_names[file_num-1]

vm_pre = hchf.restingmembrane(fine_con_screen, pre_cell)

vm0, pre_signal = fine.get_pre_signal(fine_con_screen, pre_cell)
pa0, post_signal = fine.get_post_signal(fine_con_screen, post_cell)

hzs = ['200Hz', '100Hz', '50Hz','20Hz', '10Hz']
for d in range(5):
    hz = hzs[d]
    mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = fine.get_analysis_window(pre_signal, post_signal, hz)


post_peaks = con_param.find_postsynaptic_peaks(post_window, preAPs_shifted)
post_local_baseline = con_param.get_psp_baselines(post_window,preAPs_shifted)
onsets = con_param.get_onsets(preAPs_shifted, post_window, post_peaks, post_local_baseline)
latency = con_param.latencies(onsets, preAPs_shifted)


#%%
import neo
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import human_characterisation_functions as hchf

filename = fine_con_screen
r = neo.io.AxonIO(filename=filename)
block=r.read(signal_group_mode = 'split-all')[0]
sweeps = len(block.segments)
sweep_len = len(block.segments[0].analogsignals[0])
channels = len(block.segments[0].analogsignals) # no of channels in dep rec
channel_dict = {}
for ch in range(0, channels):
    name = block.segments[0].analogsignals[ch].name
    if name == '_Ipatch':
        name = 'Ch1'
        block.segments[0].analogsignals[ch].name = 'Ch1'
    channel_dict['AnalogSig%d' % (ch)] = name

dir_plots = "plots"
end_filename = filename.rfind('/')+1
path1 = os.path.join(filename[:end_filename], dir_plots)
if os.path.isdir(path1) == False: os.mkdir(path1)
path = os.path.join(filename[:end_filename],dir_plots + '/traces')
if os.path.isdir(path) == False: os.mkdir(path)

plt.style.use(['fast'])

ax = plt.figure(figsize=(6,6))

i = 1
cell_chan  = 3
ax = plt.subplot(1,1, i)

middle_swp_num = int(sweeps/2)
ch1 = block.segments[middle_swp_num].analogsignals[cell_chan-1].view(np.recarray).reshape(sweep_len).tolist()
# plot_swp_num = int(ch1.shape[1]/2+1)
ch_name = block.segments[0].analogsignals[cell_chan-1].name
sampl_rate = block.segments[middle_swp_num].analogsignals[cell_chan-1].sampling_rate
units = block.segments[middle_swp_num].analogsignals[cell_chan-1].units
times = np.linspace(0,sweep_len,sweep_len)/sampl_rate
ax.plot(times,ch1, lw=0.5)
ax.set_xlabel('sec')
ax.set_ylabel(str(units)[-2:])
ax.title.set_text(ch_name)
del ch1

ax.patch.set_facecolor('white')
plt.show()
# %%

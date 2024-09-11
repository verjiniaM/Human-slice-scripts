# %%
import os
import pandas as pd
from funcs_human_characterisation import load_traces
import numpy as np
import matplotlib.pyplot as plt
from ephys_analysis.stimulation_windows_ms import stim_window_con_screen
import neo
import funcs_sorting as sort
import funcs_human_characterisation as hcf

#%%

human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
OP = 'OP240215'
patcher = 'Verji'

work_dir, filenames, indices_dict, slice_names, pre_chans, post_chans = sort.get_OP_metadata(human_dir, OP, patcher)

#%%

## Fill in some things
file_index = [13,49]
rheos = [550, 550]
active_channels = [[1,2,5,6,8], [1,2,3,4,8]]
y_labels_D2 = ['repatch \n Ch 1', 'repatch \n Ch 3', 
'repatch \n  Ch 5','repatch \n Ch6', 'repatch \n Ch 8']

inj=[-300,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400,450,500,550,600,700,800,900,1000,1100,1200,1300,1400]
clrs = ['#377eb8', '#ff7f00', '#4daf4a',
        '#f781bf', '#984ea3','#999999', 
        '#e41a1c', '#dede00']
#match clrs for channels in 2 to be the same as repatched in 1
clrs2 = ['#377eb8', '#ff7f00', '#4daf4a',
        '#f781bf', '#984ea3']

#plot D1
fn = work_dir + filenames[file_index[0]]
x = len(active_channels[0])
fig, ax = plt.subplots(x,1,sharex=True, sharey=False,figsize=(10,12))
charact_data = load_traces(fn)

for i, ch in enumerate(active_channels[0]):
    key = 'Ch' + str(ch)
    ch1 = charact_data[key][0]

    ax[i].plot(ch1[0:25000,inj.index(rheos[0])], lw = 4, c=clrs[i]) 
    for h in range(25):
        ax[i].plot(ch1[0:25000,h], lw = 2, c=clrs[i], alpha = 0.25)
    ax[i].set_ylabel("Ch " + str(ch), fontsize = 19)
    ax[i].set_xticks([])
    ax[i].set_yticks([])
    ax[i].spines['top'].set_visible(False)
    ax[i].spines['right'].set_visible(False)
    ax[i].spines['bottom'].set_visible(False)
    ax[i].spines['left'].set_visible(False)

# ax[i+1].plot(ch1[0:25000,inj.index(600)], lw=4, c=clrs[i+1]) 
# for k in range(25):
#     ax[i+1].plot(ch1[0:25000,k], lw = 2, c=clrs[i+1], alpha = 0.25)
# #ax[i+1].plot(ch1[0:25000,0], lw=3, c=clrs[i+1]) 
# ax[i+1].set_ylabel('Ch ' + str(active_channels[0][-1]), fontsize = 19)
# ax[i+1].set_xticks([])
# ax[i+1].set_yticks([])
# ax[i+1].spines['top'].set_visible(False)
# ax[i+1].spines['right'].set_visible(False)
# ax[i+1].spines['bottom'].set_visible(False)
# ax[i+1].spines['left'].set_visible(False)


fig.suptitle('Day ' + str(0+1) + ',' + str(rheos[0])+'pA',fontsize=20)
fig.patch.set_facecolor('white')
fig.tight_layout()
#plt.show()
plt.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/for_lab_seminar_13.03.24/traces/charact_D1.pdf', format = 'pdf')
#plt.close()

fn2 = work_dir + filenames[file_index[1]]
charact_data2 = load_traces(fn2)

fig2, ax = plt.subplots(x,1,sharex=True, sharey=False,figsize=(10,12))
for j, ch2 in enumerate(active_channels[1]):
    key = 'Ch' + str(ch2)
    ch1 = charact_data2[key][0]
    
    ax[j].plot(ch1[0:25000,inj.index(rheos[1])], lw = 4, c=clrs2[j]) 
    for k in range(25):
        ax[j].plot(ch1[0:25000,k], lw = 2, c=clrs2[j], alpha = 0.25)
    #ax[j].plot(ch1[0:25000,0], lw = 3, c=clrs2[j]) 
    ax[j].set_ylabel(y_labels_D2[j], fontsize = 19)
    ax[j].set_xticks([])
    ax[j].set_yticks([])
    ax[j].spines['top'].set_visible(False)
    ax[j].spines['right'].set_visible(False)
    ax[j].spines['bottom'].set_visible(False)
    ax[j].spines['left'].set_visible(False)

# ax[j+1].plot(ch1[0:25000,inj.index(600)], lw=4, c=clrs2[j+1]) 
# for k in range(25):
#     ax[j+1].plot(ch1[0:25000,k], lw = 2, c=clrs2[j+1], alpha = 0.25)
# #ax[j+1].plot(ch1[0:25000,0], lw=3, c=clrs2[j+1]) 
# ax[j+1].set_ylabel(y_labels_D2[j+1], fontsize = 19)
# ax[j+1].set_xticks([])
# ax[j+1].set_yticks([])
# ax[j+1].spines['top'].set_visible(False)
# ax[j+1].spines['right'].set_visible(False)
# ax[j+1].spines['bottom'].set_visible(False)
# ax[j+1].spines['left'].set_visible(False)

fig2.suptitle('Day ' + str(1+1) + ',' + str(rheos[1])+'pA',fontsize=20)
fig2.patch.set_facecolor('white')
fig2.tight_layout()
#plt.show()
plt.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/for_lab_seminar_13.03.24/traces/charact_D2.pdf', format = 'pdf')



#%%
# plots connectivity
#Set the colors so that the repatched cells are in the same colors
fn = work_dir + filenames[15]
fn2 = work_dir + filenames[51]

z1 = 1.5
z2 = 40.5
stim_window = stim_window_con_screen #what are the stim windows???
con_screen_data = hcf.load_traces(fn)

x = len(active_channels[0]) 
stim_window = stim_window_con_screen

fig, ax = plt.subplots(x, x, sharex = True, sharey = False, figsize = (6,6))
for i, ch1 in enumerate(active_channels[0]):
    ch_name_i = 'Ch' + str(ch1)
    ch_data_i = con_screen_data[ch_name_i][0]
    avg = np.mean(ch_data_i, axis = 1)

    for j, ch2 in enumerate(active_channels[0]):
        # if i == j:
        #     ax[i,j].plot()
        ch_name_j = 'Ch' + str(ch2)
        plotwin = avg[stim_window[ch_name_j][0]:stim_window[ch_name_j][1]] #from the average takes the signal for this stim_window
        ax[i,j].plot(plotwin, clrs[i], lw=1.25)
        ax[i,0].set_ylabel(str(ch_name_i))
        ax[i,j].yaxis.label.set_color(clrs[i])
        ax[i,j].set_xticks([])
        ax[i,j].set_yticks([])
        ax[i,j].spines['top'].set_visible(False)
        ax[i,j].spines['right'].set_visible(False)
        ax[i,j].spines['bottom'].set_visible(False)
        ax[i,j].spines['left'].set_visible(False)
        if plotwin.max()-plotwin.min() < 10: 
            ax[i,j].set_ylim([plotwin.min() - z1, plotwin.max() + z1])
            v1 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z1, lw=0.7, color='k') 
        else:
            ax[i,j].set_ylim([plotwin.min() - z1, plotwin.max() + z2])
            v2 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z2, lw=0.7, color='k')

fig.suptitle('connections in ' + fn[-10:-4],fontsize=15)
fig.patch.set_facecolor('white')
fig.tight_layout()
# plt.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/for_lab_seminar_13.03.24/traces/con_screen_D1.pdf', format = 'pdf')


con_screen_data2 = hcf.load_traces(fn2)

x = len(active_channels[1]) 

fig, ax = plt.subplots(x, x, sharex = True, sharey = False, figsize = (6,6))
for i, ch1 in enumerate(active_channels[1]):
    ch_name_i = 'Ch' + str(ch1)
    ch_data_i = con_screen_data2[ch_name_i][0]
    avg = np.mean(ch_data_i, axis = 1)

    for j, ch2 in enumerate(active_channels[1]):
        # if i == j:
        #     ax[i,j].plot()
        ch_name_j = 'Ch' + str(ch2)
        plotwin = avg[stim_window[ch_name_j][0]:stim_window[ch_name_j][1]] #from the average takes the signal for this stim_window
        ax[i,j].plot(plotwin, clrs[i], lw=1.25)
        ax[i,0].set_ylabel(str(ch_name_i))
        ax[i,j].yaxis.label.set_color(clrs[i])
        ax[i,j].set_xticks([])
        ax[i,j].set_yticks([])
        ax[i,j].spines['top'].set_visible(False)
        ax[i,j].spines['right'].set_visible(False)
        ax[i,j].spines['bottom'].set_visible(False)
        ax[i,j].spines['left'].set_visible(False)
        if plotwin.max()-plotwin.min() < 10: 
            ax[i,j].set_ylim([plotwin.min() - z1, plotwin.max() + z1])
            v1 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z1, lw=0.7, color='k') 
        else:
            ax[i,j].set_ylim([plotwin.min() - z1, plotwin.max() + z2])
            v2 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z2, lw=0.7, color='k')

fig.suptitle('connections in ' + fn2[-10:-4],fontsize=15)
fig.patch.set_facecolor('white')
fig.tight_layout()
# plt.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/for_lab_seminar_13.03.24/traces/con_screen_D2.pdf', format = 'pdf')


# # results = work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.xlsx'
# # df = pd.read_excel(results)
# # df.to_csv(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.csv')

# # %%

# %%


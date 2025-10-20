import math
import pyabf
import numpy as np
import matplotlib.pyplot as plt
import ephys_analysis.funcs_sorting as sort
import ephys_analysis.funcs_human_characterisation as hcf

fn = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP250213/25213004.abf'
channel = 1

inj = hcf.get_inj_current_steps(fn)
end_fn = fn.rfind('/') + 1
trace = pyabf.ABF(fn)
mask = trace.sweepC != 0
inj_trace = np.zeros([len(inj[:-1]), len(mask)])
inj_trace[:, mask] = np.array(inj[:-1])[:, np.newaxis]

if len(trace.channelList) < 8:
    if '_Ipatch' in trace.adcNames:
        trace.adcNames[trace.adcNames.index('_Ipatch')] = 'Ch1'
    if 'IN0' in trace.adcNames:
        trace.adcNames[trace.adcNames.index('IN0')] = 'Ch1'
    channel_name = 'Ch' + str(channel)
    channel = trace.channelList[trace.adcNames.index(channel_name)]
else:
    channel = channel-1
    channel_name = 'Ch' + str(channel+1)

# all conseqitive
plt.style.use(['fast'])

fig, ax = plt.subplots(2,1, sharex = True, figsize = (14,8))
for swp in range(trace.sweepCount):
    trace.setSweep(sweepNumber = swp, channel = channel)
    ax[0].plot(trace.sweepX, trace.sweepY, c = '#0000FF', linewidth = 3, alpha = 0.7)
    ax[1].plot(trace.sweepX, inj_trace[swp,:], linewidth = 3)

ax[0].set_ylabel(trace.sweepLabelY)
ax[1].set_xlabel(trace.sweepLabelX)
ax[1].set_ylabel(trace.sweepLabelC)

sweep = 'all'
fig.suptitle('{0}, sweep num {1} , channel {2}'.format(fn[end_fn:], str(sweep),channel_name))

save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/monograph/sup_figs/char_overlay_'
if save_dir:
    plt.savefig(f'{save_dir}trace_{fn[end_fn:]}_{channel_name}_swp_{sweep}.svg',
    format = 'svg', bbox_inches = 'tight', dpi = 300)

# some overlay
fig, ax = plt.subplots(2,1, sharex = True, figsize = (14,8))
for swp in [0,1,2,3,4,5,6]:
    trace.setSweep(sweepNumber = swp, channel = channel)
    ax[0].plot(trace.sweepX, trace.sweepY, c = '#0000FF', linewidth = 3, alpha = 0.7)
    ax[1].plot(trace.sweepX, inj_trace[swp,:], linewidth = 3)

ax[0].set_ylabel(trace.sweepLabelY)
ax[1].set_xlabel(trace.sweepLabelX)
ax[1].set_ylabel(trace.sweepLabelC)

sweep = 'sup_3'
fig.suptitle('{0}, sweep num {1} , channel {2}'.format(fn[end_fn:], str(sweep),channel_name))

save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/monograph/sup_figs/char_min_max_'
if save_dir:
    plt.savefig(f'{save_dir}trace_{fn[end_fn:]}_{channel_name}_swp_{sweep}.svg',
    format = 'svg', bbox_inches = 'tight', dpi = 300)


# All overlay
fig, ax = plt.subplots(2,1, sharex = True, figsize = (24,8))

ax[0].plot(trace.data[0], c = '#0000FF', linewidth = 3, alpha = 0.7)
ax[1].plot(inj_trace.flatten(), linewidth = 3)
ax[0].set_ylabel(trace.sweepLabelY)
ax[1].set_xlabel(trace.sweepLabelX)
ax[1].set_ylabel(trace.sweepLabelC)
save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/monograph/sup_figs/char_flat_'
if save_dir:
    plt.savefig(f'{save_dir}trace_{fn[end_fn:]}_{channel_name}_swp_{sweep}.svg',
    format = 'svg', bbox_inches = 'tight', dpi = 300)



#%%
# FOR SUP FIG 3

fn = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP250213/25213004.abf'
channel = 1


inj = hcf.get_inj_current_steps(fn)
end_fn = fn.rfind('/') + 1
trace = pyabf.ABF(fn)
mask = trace.sweepC != 0
inj_trace = np.zeros([len(inj[:-1]), len(mask)])
inj_trace[:, mask] = np.array(inj[:-1])[:, np.newaxis]

if len(trace.channelList) < 8:
    if '_Ipatch' in trace.adcNames:
        trace.adcNames[trace.adcNames.index('_Ipatch')] = 'Ch1'
    if 'IN0' in trace.adcNames:
        trace.adcNames[trace.adcNames.index('IN0')] = 'Ch1'
    channel_name = 'Ch' + str(channel)
    channel = trace.channelList[trace.adcNames.index(channel_name)]
else:
    channel = channel-1
    channel_name = 'Ch' + str(channel+1)

plt.style.use('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/Human-slice-scripts/style_plot_paper.mplstyle')

w_cm = 18
h_cm = 11
fig, ax = plt.subplots(2, 1, figsize = (w_cm/2.54, h_cm / 2.54),
                       gridspec_kw={'height_ratios': [2, 0.5]})
for swp in [0,1,2,3,4,5,7]:
    trace.setSweep(sweepNumber = swp, channel = channel)
    ax[0].plot(trace.sweepX, trace.sweepY, c = "#2B2B45", linewidth = 1.6, alpha = 0.8)
    ax[1].plot(trace.sweepX, inj_trace[swp,:], c = "#2B2B45", linewidth = 1.9)

ax[0].set_ylabel(trace.sweepLabelY)
ax[0].set_xlabel(trace.sweepLabelX)
ax[1].set_ylabel(trace.sweepLabelC)

sweep = 'sup_3'
plt.subplots_adjust(hspace=0)  # Increase or decrease as needed
fig.suptitle('{0}, sweep num {1} , channel {2}'.format(fn[end_fn:], str(sweep),channel_name))

save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/monograph/sup_figs/char_min_max_'
if save_dir:
    plt.savefig(f'{save_dir}trace_{fn[end_fn:]}_{channel_name}_swp_{sweep}.svg',
    format = 'svg', bbox_inches = 'tight', dpi = 300)


# plot single AP
channels = [1]
filename = fn
w_cm2 = 7
h_cm2 = 7

end_fn = filename.rfind('/') + 1
dir_ap_props = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/monograph/sup_figs/'
    
charact_data = hcf.load_traces(filename)
inj = hcf.read_inj(inj)
max_spikes = hcf.get_max_spikes(charact_data, channels)
first_spikes, peaks_all, spike_counts_all, first_spiking_sweeps_all = hcf.get_ap_param_for_plotting(charact_data, channels, inj, max_spikes)
Rheobase_all, AP_all, THloc_all, TH_all, APheight_all, max_depol_all, max_repol_all = hcf.get_ap_param(charact_data, channels, inj, max_spikes)

for n, ch in enumerate(channels):
    key = 'Ch' + str(ch)
    ch_data = charact_data[key][0]

    if math.isnan(first_spikes[n]):
        fig, ax = plt.subplots(1,1,figsize= (w_cm/2.54, h_cm / 2.54))
        fig.suptitle(key, fontsize=15)
        plt.savefig(f'{dir_ap_props}{filename[end_fn:-4]}_AP#{ap+1}_{key}.png')
        plt.close()
        continue

    if np.max(spike_counts_all[n][:,1]) == 1:
        ap = 0
    else:
        ap = 1
    
    ap_trace = AP_all[n][120:]
    x_time = np.linspace(0, len(ap_trace)/20, len(ap_trace))
    
    fig, ax = plt.subplots(1,1,figsize= (w_cm/2.54, h_cm / 2.54))
    ax.plot(x_time, ap_trace, linewidth = 2.3, c = "#2B2B45")

    for f in [0.1, 0.5, 0.9]:
        y_param = TH_all[n] + f*APheight_all[n]
        # closest_index = np.argmin(np.abs(ap_trace[(THloc_all[n]-120):80] - y_param))

        x_depol = x_time[:80][(THloc_all[n]-120) +
                              (np.abs(ap_trace[(THloc_all[n]-120):80] - y_param)).argmin()]
        x_repol = x_time[80:][(np.abs(ap_trace[80:] - y_param)).argmin()]
        ax.scatter([x_depol, x_repol], [y_param, y_param], c = "#8748BE", s=80, marker = '_')

    plt.scatter((THloc_all[n]-120)/20, TH_all[n], c = "#730B4A", s = 80)
    plt.scatter(80/20, peaks_all[n][first_spiking_sweeps_all[n], ap, 2], c = "#0870CC", s = 45)
    
    fig.suptitle(f'Ch: {key}, AP# {ap+1}, \
    TH = {round(TH_all[n],2)}, \
    amp = {round(APheight_all[n],2)}')

    plt.savefig(f'{dir_ap_props}{filename[end_fn:-4]}_AP#{ap+1}_{key}.png')
    plt.close()




# plot Hyperpol
filename = fn
ch = 1
inj = hcf.get_inj_current_steps(fn)
onset = 2624
offset = 22624
w_cm2 = 18
h_cm2 = 8

end_fn = filename.rfind('/') + 1
dir_onset = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/monograph/sup_figs/'

charact_data = hcf.load_traces(filename)
tau_all, capacitance_all, mcs, V65s, RMPs_char = hcf.get_hyperpolar_param(charact_data, channels, inj)

key = 'Ch' + str(ch)
ch_data = charact_data[key][0]

fig, ax = plt.subplots(1,1,figsize= (w_cm2/2.54, h_cm2 /2.54))
for i in range(0,5):
    swp = list(ch_data[:,i])
    bl = np.median(ch_data[0:onset-20,i])

    res = list(filter(lambda ii: ii < V65s[n][i], swp))[0] #takes the first value in swp < V65
    tau65 = swp.index(res) #index of res
    tc = tau65 - onset
    ax.plot(np.linspace(0, len(ch_data[:,i])/20,len(ch_data[:,i])),ch_data[:,i], c = "#2B2B45", alpha = 0.75)
    ax.scatter((onset + tc)/20, V65s[n][i], c =  "#FC8739", marker = 'o', s = 60, zorder = 11)
    ax.annotate('V65  ', ((onset + tc)/20 - 20, V65s[n][i]), horizontalalignment='right')
    ax.scatter(onset/20, bl, c =  "#8748BE", marker = '_', s = 70, zorder = 10)
    ax.set_xlabel('ms')
    ax.set_ylabel('mV')
plt.annotate('  Baseline', (onset/20, bl))
fig.patch.set_facecolor('white')
plt.title(key)
plt.savefig(f'{dir_onset}/Char_onset_plot_{filename[end_fn:-4]}_{key}.png')


# zoom
end_ = 300*20
w_cm2 = 6
h_cm2 = 8

fig, ax = plt.subplots(1,1,figsize= (w_cm2/2.54, h_cm2/2.54))
for i in range(0,5):
    y_data = ch_data[:,i][:end_]
    swp = list(y_data)
    bl = np.median(ch_data[0:onset-20,i])

    res = list(filter(lambda ii: ii < V65s[n][i], swp))[0] #takes the first value in swp < V65
    tau65 = swp.index(res) #index of res
    tc = tau65 - onset
    ax.plot(np.linspace(0, len(y_data)/20,len(y_data)), y_data, c = "#2B2B45", alpha = 0.75)
    ax.scatter((onset + tc)/20, V65s[n][i], c =  "#FC8739", marker = 'o', s = 60, zorder = 11)
    ax.annotate('V65  ', ((onset + tc)/20 - 20, V65s[n][i]), horizontalalignment='right')
    ax.scatter(onset/20, bl, c =  "#8748BE", marker = '_', s = 70, zorder = 10)
    ax.set_xlabel('ms')
    ax.set_ylabel('mV')
plt.annotate('  Baseline', (onset/20, bl))
fig.patch.set_facecolor('white')
plt.title(key)
plt.savefig(f'{dir_onset}/Char_onset_ZOOM_plot_{filename[end_fn:-4]}_{key}.png')

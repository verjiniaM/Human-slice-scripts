from tkinter import font
import numpy as np
import matplotlib.pyplot as plt
import math
import ephys_analysis.stimulation_windows_ms as stim_win
import ephys_analysis.funcs_con_screen as con_param
import ephys_analysis.funcs_sorting as sort
import ephys_analysis.funcs_human_characterisation as hcf
import pyabf
from ipywidgets import interact

plt.style.use('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/Human-slice-scripts/style_plot_intrinsic.mplstyle')


# plots the middle sweep for each channel
# check the traces to see which channels were active or if the protocol names are enetered correctly
def plot_middle_sweep (filename): 
    end_fn = filename.rfind('/') + 1
    dir_plots = sort.make_dir_if_not_existing(filename[:end_fn], 'plots')
    dir_traces = sort.make_dir_if_not_existing(dir_plots, 'traces')

    data_dict = hcf.load_traces(filename)
    all_chans = list(data_dict.keys())
    any_chan = int((all_chans[0])[2])
    sweep_len = np.shape(data_dict[all_chans[0]][0])[0]
    sweep_count = np.shape(data_dict[all_chans[0]][0])[1]
    sampl_rate, units, times = hcf.get_abf_info(filename, any_chan, sweep_count, sweep_len)
    middle_swp_num = int(sweep_count/2)

    plt.style.use(['fast'])
    x = math.ceil(np.sqrt(len(all_chans)))
    fig = plt.figure(figsize=(6,6))
    plt.subplots_adjust(hspace=0.5)

    for i in range(1, len(all_chans)+1):
        key = all_chans[i-1]
        ch_data = data_dict[key][0]

        ax = plt.subplot(x, x, i)

        signal = ch_data[:,middle_swp_num]
        ax.plot(times,signal, lw=0.5)
        ax.set_xlabel('sec')
        ax.set_ylabel(str(units)[-2:])
        ax.title.set_text(key)

    fig.tight_layout() 
    fig.patch.set_facecolor('white')
    plt.savefig(dir_traces + '/trace_plot_' + filename[end_fn:-4] + '.png')
    plt.close(fig)
    return 'Trace plots saved in' + dir_traces

def plot_vc_holding (filename, channels):
    plt.style.use(['fast']) 
    end_fn = filename.rfind('/') + 1
    dir_vc_plots = sort.make_dir_if_not_existing(filename[:end_fn] + 'plots/', "vc_plots")

    vc_data = hcf.load_traces(filename)
    all_chans = list(vc_data.keys())
    sweep_count = np.shape(vc_data[all_chans[0]][0])[1]
    Res = np.ndarray([sweep_count,3])  

    for ch in channels:
        key = 'Ch' + str(ch)
        ch_data = vc_data[key][0]

        for n, trace in enumerate(ch_data.T):
            HC = np.median(trace[4000:]) #holding current
            avgBL = np.median(trace[10:990])
            minP = np.min(trace[990:1300])
            ss = np.median(trace[1400:1900])
            RaI = avgBL - minP
            RiI = avgBL - ss

            if RaI == 0 or RiI == 0: 
                THholding = Res[0,0]*0.1
                Res[n,0] = HC
                Res[n,1]= math.nan
                Res[n,2]= math.nan
                fig,ax = plt.subplots(1,1,sharex=True)
                ax.scatter(range(sweep_count),Res[:,0], color='b')
                ax.set_title('Holding current')
                ax.fill_between(range(sweep_count),Res[0,0]-THholding, Res[0,0]+THholding, color='b', alpha=0.1)
                ax.set_xlabel('Sweep num')
            else:
                ResRa = (0.004/(RaI*1e-12))/1000000
                ResRi = (0.004/(RiI*1e-12))/1000000
                Res[n,0] = HC
                Res[n,1] = ResRa
                Res[n,2] = ResRi
        
                THholding = Res[0,0] * 0.1
                THRs = Res[0,1] * 0.1 #10 %
                THRi = Res[0,2] * 0.1
        
                fig, ax = plt.subplots(3,1,sharex=True)
                hold = ax[0].scatter(range(sweep_count), Res[:,0], color='b', label = 'Holding current (median)')
                #ax[0].set_title('Holding current')
                hold_win = ax[0].fill_between(range(sweep_count), Res[0,0] - THholding, Res[0,0] + THholding, color='b', alpha=0.1, label = '20% window')
                ax[0].set_ylim([np.min(Res[0,0] - 2 * THholding), np.max(Res[0,0] + 2 * THholding)])
                ax[0].set_yticks([np.min(Res[0,0] - 2 * THholding), Res[0,0], np.max(Res[0,0] + 2 * THholding)])
                #ax[0].legend(loc = 'upper right')

                Rs = ax[1].scatter(range(sweep_count), Res[:,1], color='r', label = 'Series resistance (median)')
                #ax[1].set_title('Series resistance')
                Rs_win = ax[1].fill_between(range(sweep_count), Res[0,1] - THRs, Res[0,1] + THRs, color='r', alpha=0.1, label = '20% window')
                ax[1].set_ylim([np.min(Res[0,1] - 2 * THRs), np.max(Res[0,1] + 2 * THRs)])
                ax[1].set_yticks([np.min(Res[0,1] - 2 * THRs), Res[0,1], np.max(Res[0,1] + 2 * THRs)])
                #ax[1].legend(loc = 'upper right')

                Rin = ax[2].scatter(range(sweep_count), Res[:,2], color='g', label = 'Input resistance (median)')
                #ax[2].set_title('Input resistance')
                Rin_win = ax[2].fill_between(range(sweep_count), Res[0,2] - THRi, Res[0,2] + THRi, color='g', alpha=0.1, label = '20% window')
                ax[2].set_ylim([np.min(Res[0,2] - 2 * THRi), np.max(Res[0,2] + 2 * THRi)])
                ax[2].set_yticks([np.min(Res[0,2] - 2 * THRi), Res[0,2], np.max(Res[0,2] + 2 * THRi)])
                ax[2].set_xlabel('Sweep num')
            
                lgd = plt.legend(
                    [hold, hold_win, Rs, Rs_win, Rin, Rin_win],
                    ['Holding current (median)','20% window', 'Series resistance (median)', '20% window', 'Input resistance (median)', '20% window'],
                    ## by default, legend anchor to axis, but can
                    ## also be anchored to arbitary position
                    ## positions within [1,1] would be within the figure
                    ## all numbers are ratio by default

                    bbox_to_anchor=(1.05, 1),

                    ## loc indicates the position within the figure
                    ## it is defined consistent to the same Matlab function 
                    loc='center left',

                    ncol = 1
                    #mode="expand",
                    #borderaxespad=0.
                    )
                
            fig.patch.set_facecolor('white') 
            fig.suptitle(key, fontsize = 15)
            #fig.tight_layout()

            plt.savefig(dir_vc_plots + '/' + filename[end_fn:-4] + '_'+ key + '_VC_plot.png', \
                bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close(fig)


def plot_hyperpolar(filename, channels, inj, onset = 2624, offset = 22624,
clrs = ["b", "g", "r", "c", "m", "y", "#FF4500", "#800080"]):
    end_fn = filename.rfind('/') + 1
    dir_onset = sort.make_dir_if_not_existing(filename[:end_fn] + 'plots/', 'Onset')

    charact_data = hcf.load_traces(filename)
    inj = hcf.read_inj(inj)
    tau_all, capacitance_all, mcs, V65s, RMPs_char = hcf.get_hyperpolar_param(charact_data, channels, inj)

    for n, ch in enumerate(channels):
        key = 'Ch' + str(ch)
        ch_data = charact_data[key][0]

        fig = plt.figure()
        for i in range(0,5):
            swp = list(ch_data[:,i])
            bl = np.median(ch_data[0:onset-20,i])
            if list(filter(lambda ii: ii < V65s[n][i], swp)) == []:
                print("No hyperpolarization fig for " + filename[end_fn:-4] + key)
            else:
                res = list(filter(lambda ii: ii < V65s[n][i], swp))[0] #takes the first value in swp < V65
                tau65 = swp.index(res) #index of res    
                tc = tau65 - onset
                plt.plot(ch_data[:,i], c = clrs[i])
                plt.scatter(onset + tc, V65s[n][i], c=clrs[i])
                plt.annotate('V65  ', (onset + tc, V65s[n][i]), horizontalalignment='right')
                plt.scatter(onset, bl, c='r')
                plt.ylabel('mV')
        plt.annotate('  Baseline', (onset, bl))
        fig.patch.set_facecolor('white')    
        plt.title(key)
        plt.savefig(dir_onset + '/Char_onset_plot_' + filename[end_fn:-4]+'_'+ key + '.png')
        plt.close()
        
def plot_spikes (filename, channels, inj):
    '''
    plots alll sweeps with spikes
    used to for user to evaluate if the max_spikes function has worked   
    '''
    end_fn = filename.rfind('/') + 1
    dir_spikes = sort.make_dir_if_not_existing(filename[:end_fn] + 'plots/',  'Max_Spikes')

    charact_data = hcf.load_traces(filename)
    inj = hcf.read_inj(inj)
    max_spikes = hcf.get_max_spikes(charact_data, channels)
    first_spikes, peaks_all, spike_counts_all, fs = hcf.get_ap_param_for_plotting(charact_data, channels, inj, max_spikes)

    for n, ch in enumerate(channels):
        key = 'Ch' + str(ch)
        ch_data = charact_data[key][0]

        win = len(inj) - first_spikes[n]
        if math.isnan(first_spikes[n]):
            x = 1
            fig, ax = plt.subplots(x,x,sharex=True, sharey=False,figsize=(10,4))
            fig.suptitle(key, fontsize=15)
            plt.savefig(dir_spikes + '/char_spikes_plot_' + filename[end_fn:-4]+'_'+ key + '.png')
            plt.close(fig)
            continue
        x = math.ceil(np.sqrt(win))
        fig = plt.figure(figsize=(22,9))
        for i in range(win):
            ax = fig.add_subplot(x,x, i+1)
            if first_spikes[n] + i-1 < np.shape(ch_data)[1]: #plotting for last sweep with no spikes and all following sweeps
                ax.plot(ch_data[:,first_spikes[n] + i-1], lw=0.5, c='grey')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                #ax.spines['left'].set_visible(False)
                ax.scatter(peaks_all[n][first_spikes[n] + i-1, :, 1], 
                        peaks_all[n][first_spikes[n] + i-1, :, 2], marker='+', c='r')
                ax.annotate('+'+str(inj[first_spikes[n]+ i - 1])+' pA', (25500,0), (25500,0), color='b', rotation=90)
        
        fig.patch.set_facecolor('white')
        fig.suptitle(key, fontsize=15)
        #fig.tight_layout()
        plt.savefig(dir_spikes + '/char_spikes_plot_' + filename[end_fn:-4]+'_'+ key + '.png')
        plt.close(fig)

def plot_iv_curve (filename, channels, inj):
    end_fn = filename.rfind('/') + 1
    dir_iv_curve = sort.make_dir_if_not_existing(filename[:end_fn] + 'plots/',  'IV_curve')
     
    charact_data = hcf.load_traces(filename)
    inj = hcf.read_inj(inj)
    max_spikes = hcf.get_max_spikes(charact_data, channels)
    first_spikes, peaks_all, spike_counts_all, fs = hcf.get_ap_param_for_plotting(charact_data, channels, inj, max_spikes)

    for n, ch in enumerate(channels):
        key = 'Ch' + str(ch)
        ch_data = charact_data[key][0]

        if math.isnan(first_spikes[n]):
            fig, ax = plt.subplots(1,1,sharex=True, sharey=False,figsize=(10,4))
            fig.suptitle(key, fontsize=15)
            plt.savefig(dir_iv_curve + '/char_IV_curve_' + filename[end_fn:-4]+'_'+ key + '.png')
            plt.close()
            continue

        IOfit = np.polyfit(spike_counts_all[n][first_spikes[n]:,0], spike_counts_all[n][first_spikes[n]:,1],1)
        IO_slope = IOfit[0]
        fig = plt.figure()
        plt.plot(spike_counts_all[n][:,0], spike_counts_all[n][:,1])
        plt.gca().set_xlabel('Injection current, pA')
        plt.gca().set_ylabel('Number of spikes')
        Rheobase = inj[first_spikes[n]]   
        
        fig.suptitle(key, fontsize=15)
        fig.patch.set_facecolor('white')
        plt.savefig(dir_iv_curve + '/char_IV_curve_' + filename[end_fn:-4]+'_'+ key + '.png')
        plt.close()

def plot_ap_props(filename, channels, inj):
    end_fn = filename.rfind('/') + 1
    dir_ap_props = sort.make_dir_if_not_existing(filename[:end_fn] + 'plots/',  'AP_props')
     
    charact_data = hcf.load_traces(filename)
    inj = hcf.read_inj(inj)
    max_spikes = hcf.get_max_spikes(charact_data, channels)
    first_spikes, peaks_all, spike_counts_all, first_spiking_sweeps_all = hcf.get_ap_param_for_plotting(charact_data, channels, inj, max_spikes)
    Rheobase_all, AP_all, THloc_all, TH_all, APheight_all, max_depol_all, max_repol_all = hcf.get_ap_param(charact_data, channels, inj, max_spikes)

    for n, ch in enumerate(channels):
        key = 'Ch' + str(ch)
        ch_data = charact_data[key][0]

        if math.isnan(first_spikes[n]):
            x = 1
            fig, ax = plt.subplots(x,x,sharex=True, sharey=False,figsize=(10,4))
            fig.suptitle(key, fontsize=15)
            plt.savefig(dir_ap_props + '/' + filename[end_fn:-4] + '_AP#' + str(0) + '_' + key + '.png')
            plt.close()
            continue

        if np.max(spike_counts_all[n][:,1]) == 1:
            ap = 0
        else:
            ap = 1

        fig = plt.figure()
        plt.plot(AP_all[n])
        plt.scatter(THloc_all[n] ,TH_all[n], color = 'red')
        plt.scatter(200,peaks_all[n][first_spiking_sweeps_all[n], ap, 2], color = 'green')
        fig.suptitle('Ch: ' + key + ', AP#' + str(ap+1) + ', TH = ' + str(round(TH_all[n],2)) + ', amp = ' + str(round(APheight_all[n],2)))
        fig.patch.set_facecolor('white') 
        plt.savefig(dir_ap_props + '/' + filename[end_fn:-4] + '_AP#' + str(ap+1) + '_' + key + '.png')
        plt.close()

def plots_for_charact_file(filename, channels, inj):

    plot_hyperpolar(filename, channels, inj)
    plot_spikes(filename, channels, inj)
    plot_iv_curve(filename, channels, inj)
    plot_ap_props(filename, channels, inj)


def plot_mini_sweeps (filename, cell_chan, sweep):
    '''
    for a given sweep plots min and max values
    '''
    sweep = sweep - 1
    data_dict = hcf.load_traces(filename)
    #no test pulse
    ch_data = data_dict['Ch' + str(cell_chan)][0]
    signal_no_test_pulse = ch_data[5500:,sweep]
    
    min_val = np.amin(signal_no_test_pulse)
    max_val = np.amax(signal_no_test_pulse)
    loc_min = np.where(np.isclose(signal_no_test_pulse, min_val))
    loc_max = np.where(np.isclose(signal_no_test_pulse, max_val))
    num_points = np.shape(signal_no_test_pulse)[0]
    times = np.linspace(0,num_points,num_points)/20000
    plt.plot(times,signal_no_test_pulse, lw=0.5)
    plt.plot(int(loc_min[0])/20000, np.amin(signal_no_test_pulse),'ro')
    plt.plot(int(loc_max[0])/20000, np.amax(signal_no_test_pulse),'bo')
    plt.xlabel('sec')
    plt.ylabel('pA')
    plt.show()


## Connectivity plotting functions
def plot_con_screen_all(fn, active_chans):

    con_screen_data = hcf.load_traces(fn)

    fig,ax = plt.subplots(len(active_chans),1, sharex = True, figsize=(10,20))
    fig.patch.set_facecolor('white')

    for i, signal in enumerate(active_chans):
        ch_name = 'Ch' + str(signal)
        ch_data = con_screen_data[ch_name][0]
        avg_signal = np.mean(ch_data, axis = 1)

        ax[i].plot(avg_signal)
        ax[i].set_ylabel(ch_name)
    
    fig.suptitle(fn[fn.rfind('/') + 1:-4])
    #plt.savefig(save_folder + fn[fn.rfind('/') + 1:-4] )
    plt.show(block=True)


def plot_connect(fn, active_channels, z1=0.5, z2=35.5,
 clrs = ["b", "g", "k", "c", "k", "y", "#FF4500", "#800080"]):
    end_fn = fn.rfind('/') + 1
    dir_connect = sort.make_dir_if_not_existing(fn[:end_fn] + 'plots/',  'connectivity_plots')
    
    con_screen_data = hcf.load_traces(fn)
    x = len(active_channels) 
    stim_window = stim_win.stim_window_con_screen

    fig, ax = plt.subplots(x, x, sharex = True, sharey = False, figsize = (6,6))
    for i, ch1 in enumerate(active_channels):
        ch_name_i = 'Ch' + str(ch1)
        ch_data_i = con_screen_data[ch_name_i][0]
        avg = np.mean(ch_data_i, axis = 1)

        for j, ch2 in enumerate(active_channels):
            if i == j:
                ax[i,j].plot()
            ch_name_j = 'Ch' + str(ch2)
            plotwin = avg[stim_window[ch_name_j][0]:stim_window[ch_name_j][1]] #from the average takes the signal for this stim_window
            ax[i,j].plot(plotwin, clrs[i], lw=0.25)
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
                v1 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z1, lw=0.2, color='k') 
            else:
                ax[i,j].set_ylim([plotwin.min() - z1, plotwin.max() + z2])
                v2 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z2, lw=0.2, color='k')
    
    fig.suptitle('connections in ' + fn[end_fn:],fontsize=15)
    fig.patch.set_facecolor('white')
    fig.tight_layout()

    plt.savefig(dir_connect + '/' + fn[end_fn:-4] + 'con_screen_plot.png')
    plt.close(fig)

def  plot_connection_window(con_screen_file, preC, postC, pre_window, post_window, preAPs_shifted, postsig,\
                           onsets, preAPs, PSPs, bl):

    sampl_rate, units, times = hcf.get_abf_info(con_screen_file, preC, np.shape(postsig)[1], np.shape(postsig)[0])
    end_fn = con_screen_file.rfind('/') + 1
    dir_connect = sort.make_dir_if_not_existing(con_screen_file[:end_fn] + 'plots/',  'connectivity_plots')
    
    fig,axarr = plt.subplots(2,1, sharex=True, figsize=(12,12))
    fig.patch.set_facecolor('white')

    my_labels = {'l1' : 'peak pre_AP', 'l2' : 'baseline', 'l3': 'onset', 'l4': 'post synnaptic peak'}
    axarr[0].plot(times[:len(pre_window)], pre_window, color='k')   #plost pre APs 0 shifted
    
    for i in range(0,len(preAPs_shifted[0])):
        axarr[0].scatter(preAPs_shifted[0][i]/sampl_rate,\
             pre_window[preAPs_shifted[0][i]], marker='o', color='r', label = my_labels['l1'])
        my_labels['l1'] = "_nolegend_"

    for i in range(np.shape(postsig)[1]):
        y = postsig[:,i][preAPs[0][0]-750:preAPs[0][len(preAPs[0])-1]+750]
        axarr[1].plot(times[0:len(y)], y, lw=0.2, color='grey', alpha=0.4, zorder=0 )
        #axarr[1].plot(y, lw=0.2, color='grey', alpha=0.4)
    
    y = post_window
    axarr[1].plot(times[0:len(y)], post_window, color='k', lw=0.5, zorder=10) #post_window is the averaged 0 shifter post signal 

    for i in range(0,len(preAPs_shifted[0])):
        sample_rate = sampl_rate.item()
        axarr[1].hlines(bl[i,0], (preAPs_shifted[0][i]-170)/sample_rate, (preAPs_shifted[0][i]-30)/sample_rate, 
        linestyles = 'solid', lw = 3, color='b', label = my_labels['l2'],zorder=10, alpha = 0.6)
        my_labels['l2'] = "_nolegend_"

    axarr[1].scatter(onsets/sampl_rate, bl, marker='^', color='r', label = my_labels['l3'],  zorder=10)

    for i in range(0,len(PSPs)):
        axarr[1].scatter(PSPs[i][1]/sampl_rate, PSPs[i][0], marker='+',\
             color='g', s=30, linewidth=10, label = my_labels['l4'],zorder=10)
        my_labels['l4'] = "_nolegend_"

    axarr[0].set_title('Pre APs')
    axarr[0].set_ylabel(str(units)[-2:])
    axarr[1].set_title('Post cell responsees')
    axarr[1].set_xlabel('sec')
    axarr[1].set_ylabel(str(units)[-2:])
    plt.figlegend(loc = 'center right',  bbox_to_anchor=(1.05, 1))
    fig.patch.set_facecolor('white')

    plt.savefig(dir_connect + '/pre_post_events_' + con_screen_file[end_fn:-4] + '_Ch' + str(preC) + '#Ch' + str(postC) + '.png')
    plt.close()

    return fig, axarr

def plot_post_cell(con_screen_file, pre_cell_chan, post_cell_chan):
    end_fn = con_screen_file.rfind('/') + 1
    dir_connect = sort.make_dir_if_not_existing(con_screen_file[:end_fn] + 'plots/',  'connectivity_plots')

    stim_window = stim_win.stim_window_con_screen

    pre_sig, es, vm0_pre = con_param.presynaptic_screen(con_screen_file, pre_cell_chan)
    post_sig, vm0_post = con_param.postsynaptic_screen(con_screen_file, post_cell_chan, es)
    mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = con_param.get_analysis_window(pre_sig, post_sig)
    sampl_rate = 200_000

    fig = plt.figure(figsize=(40,24))
    plt.subplots(1,1,sharex=True)

    my_labels = {'l1' : 'peak pre_AP'}
    xvals = preAPs[0][0]-100
    for i in range(len(post_sig[0])):   
        plt.plot(post_sig[xvals:xvals+4000] - i * 4, color ='grey', lw = 0.5)
        plt.vlines(preAPs[0]-xvals,-50, -210, lw=0.25, color = 'r', label = my_labels['l1'])
        my_labels['l1'] = "_nolegend_"

    plt.figlegend(loc = 'upper right')
    fig.patch.set_facecolor('white')
    plt.savefig(dir_connect + '/post_swps_all_' + con_screen_file[end_fn:-4] + '_' + 
    'Ch' + str(pre_cell_chan) + '_to_' + 'Ch' + str(post_cell_chan) + '.png')
    plt.close()

    fig = plt.figure(figsize=(6,6))
    y = mean_post[xvals:xvals+4000]+10
    times = np.linspace(0,len(y), len(y))/(sampl_rate*1000) #foe ms
    plt.plot(times, y, lw=1.5, color='k')
    plt.xlabel('ms')
    plt.ylabel('mV')
    fig.patch.set_facecolor('white')
    plt.savefig(dir_connect + '/post_swps_mean_' + con_screen_file[end_fn:-4] + '_' + 
    'Ch' + str(pre_cell_chan) + '_to_' + 'Ch' + str(post_cell_chan) + '.png')
    plt.close()


def plot_connect_old_win(fn, active_channels, z1=0.5, z2=40.5,
 clrs = ["b", "g", "k", "c", "k", "y", "#FF4500", "#800080"]):
    end_fn = fn.rfind('/') + 1
    dir_connect = sort.make_dir_if_not_existing(fn[:end_fn] + 'plots/',  'connectivity_plots')
    
    con_screen_data = hcf.load_traces(fn)
    x = len(active_channels) 
    stim_window = stim_win.stim_window_con_screen_old

    fig, ax = plt.subplots(x, x, sharex = True, sharey = False, figsize = (6,6))
    for i, ch1 in enumerate(active_channels):
        ch_name_i = 'Ch' + str(ch1)
        ch_data_i = con_screen_data[ch_name_i][0]
        avg = np.mean(ch_data_i, axis = 1)

        for j, ch2 in enumerate(active_channels):
            # if i == j:
            #     ax[i,j].plot()
            ch_name_j = 'Ch' + str(ch2)
            plotwin = avg[stim_window[ch_name_j][0]:stim_window[ch_name_j][1]] #from the average takes the signal for this stim_window
            ax[i,j].plot(plotwin, clrs[i], lw=0.25)
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
                v1 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z1, lw=0.2, color='k') 
            else:
                ax[i,j].set_ylim([plotwin.min() - z1, plotwin.max() + z2])
                v2 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z2, lw=0.2, color='k')
    
    fig.suptitle('connections in ' + fn[end_fn:],fontsize=15)
    fig.patch.set_facecolor('white')
    fig.tight_layout()

    plt.savefig(dir_connect + '/' + fn[end_fn:-4] + 'con_screen_plot.png')
    plt.close(fig)

def plot_post_cell_old_win(con_screen_file, pre_cell_chan, post_cell_chan):
    end_fn = con_screen_file.rfind('/') + 1
    dir_connect = sort.make_dir_if_not_existing(con_screen_file[:end_fn] + 'plots/',  'connectivity_plots')

    stim_window = stim_win.stim_window_con_screen_old

    pre_sig, es, vm0_pre = con_param.presynaptic_screen(con_screen_file, pre_cell_chan)
    post_sig, vm0_post = con_param.postsynaptic_screen(con_screen_file, post_cell_chan, es)
    mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = con_param.get_analysis_window(pre_sig, post_sig)
    sampl_rate = 200_000

    fig = plt.figure(figsize=(20,12))
    plt.subplots(1,1,sharex=True)

    my_labels = {'l1' : 'peak pre_AP'}
    xvals = preAPs[0][0]-100
    for i in range(len(post_sig[0])):   
        plt.plot(post_sig[xvals:xvals+4000] - i * 4, color ='grey', lw = 0.5)
        plt.vlines(preAPs[0]-xvals,-50, -210, lw=0.25, color = 'r', label = my_labels['l1'])
        my_labels['l1'] = "_nolegend_"

    fig.patch.set_facecolor('white')
    plt.figlegend(loc = 'upper right')
    plt.savefig(dir_connect + '/post_swps_all_' + con_screen_file[end_fn:-4] + '_' + 
    'Ch' + str(pre_cell_chan) + '_to_' + 'Ch' + str(post_cell_chan) + '.png')
    plt.close()

    fig = plt.figure(figsize=(6,6))
    fig.patch.set_facecolor('white')
    y = mean_post[xvals:xvals+4000]+10
    times = np.linspace(0,len(y), len(y))/(sampl_rate*1000) #foe ms
    plt.plot(times, y, lw=1.5, color='k')
    plt.xlabel('ms')
    plt.ylabel('mV')
    plt.savefig(dir_connect + '/post_swps_mean_' + con_screen_file[end_fn:-4] + '_' + 
    'Ch' + str(pre_cell_chan) + '_to_' + 'Ch' + str(post_cell_chan) + '.png')
    plt.close()



#plots pre_cell in IC, post in VC 
def plot_connection_window_VC (con_screen_file, preC, postC, pre_window, post_window, preAPs_shifted, postsig,\
                           onsets, preAPs, PSPs, bl):

    sampl_rate, units, times = hcf.get_abf_info(con_screen_file, preC, np.shape(postsig)[1], np.shape(postsig)[0])
    end_fn = con_screen_file.rfind('/') + 1
    dir_connect = sort.make_dir_if_not_existing(con_screen_file[:end_fn] + 'plots/',  'connectivity_plots')
    
    fig,axarr = plt.subplots(2,1, sharex=True, figsize=(12,12))
    fig.patch.set_facecolor('white')

    my_labels = {'l1' : 'peak pre_AP', 'l2' : 'baseline', 'l3': 'onset', 'l4': 'post synnaptic peak'}
    axarr[0].plot(times[:len(pre_window)], pre_window, color='k')   #plost pre APs 0 shifted
    
    for i in range(0,len(preAPs_shifted[0])):
        axarr[0].scatter(preAPs_shifted[0][i]/sampl_rate,\
             pre_window[preAPs_shifted[0][i]], marker='o', color='r', label = my_labels['l1'])
        my_labels['l1'] = "_nolegend_"

    for i in range(np.shape(postsig)[1]):
        y = postsig[:,i][preAPs[0][0]-750:preAPs[0][len(preAPs[0])-1]+750]
        axarr[1].plot(times[0:len(y)], y, lw=0.2, color='grey', alpha=0.4, zorder=0 )
        #axarr[1].plot(y, lw=0.2, color='grey', alpha=0.4)
    
    y = post_window
    axarr[1].plot(times[0:len(y)], post_window, color='k', lw=0.5, zorder=10) #post_window is the averaged 0 shifter post signal 

    for i in range(0,len(preAPs_shifted[0])):
        sample_rate = sampl_rate.item()
        axarr[1].hlines(bl[i,0], (preAPs_shifted[0][i]-170)/sample_rate, (preAPs_shifted[0][i]-30)/sample_rate, 
        linestyles = 'solid', lw = 3, color='b', label = my_labels['l2'],zorder=10, alpha = 0.6)
        my_labels['l2'] = "_nolegend_"

    axarr[1].scatter(onsets/sampl_rate, bl, marker='^', color='r', label = my_labels['l3'],zorder=10)

    for i in range(0,len(PSPs)):
        axarr[1].scatter(PSPs[i][1]/sampl_rate, PSPs[i][0], marker='+',\
             color='g', s=30, linewidth=10, label = my_labels['l4'],zorder=10)
        my_labels['l4'] = "_nolegend_"

    axarr[0].set_title('Pre APs')
    axarr[0].set_ylabel(str(units)[-2:])
    axarr[1].set_title('Post cell responsees')
    axarr[1].set_xlabel('sec')
    axarr[1].set_ylabel(str(units)[-2:])
    plt.figlegend(loc = 'center right',  bbox_to_anchor=(0.92, 0.5))
    fig.patch.set_facecolor('white')

    plt.savefig(dir_connect + '/VC_pre_post_events_' + con_screen_file[end_fn:-4] + '_Ch' + str(preC) + '#Ch' + str(postC) + '.png')
    plt.close()


def plot_post_cell_VC(con_screen_file, pre_cell_chan, post_cell_chan):
    end_fn = con_screen_file.rfind('/') + 1
    dir_connect = sort.make_dir_if_not_existing(con_screen_file[:end_fn] + 'plots/',  'connectivity_plots')

    stim_window = stim_win.stim_window_con_screen

    pre_sig, es, vm0_pre = con_param.presynaptic_screen_IC(con_screen_file, pre_cell_chan)
    post_sig, vm0_post = con_param.postsynaptic_screen_VC(con_screen_file, post_cell_chan, es)
    mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = con_param.get_analysis_window_VC(pre_sig, post_sig)
    sampl_rate = 200_000
    
    xvals = preAPs[0][0]-100
    fig = plt.figure(figsize=(6,6))
    fig.patch.set_facecolor('white')
    y = mean_post[xvals:xvals+4000]+10
    times = np.linspace(0,len(y), len(y))/(sampl_rate*1000) #foe ms
    plt.plot(times, y, lw=1.5, color='k')
    plt.xlabel('ms')
    plt.ylabel('mV')
    plt.savefig(dir_connect + '/VC_post_swps_mean_' + con_screen_file[end_fn:-4] + '_' + 
    'Ch' + str(pre_cell_chan) + '_to_' + 'Ch' + str(post_cell_chan) + '.png')
    plt.close()


def plot_rheobase_trace(fn, chans):
    '''
    plots detected rheobase
    '''
    end_fn = fn.rfind('/') + 1
    dir_plots = sort.make_dir_if_not_existing(fn[:end_fn], 'plots')
    dir_rheobase = sort.make_dir_if_not_existing(dir_plots, 'rheobase')
    
    rheos, THs, THs_in_trace, swps = hcf.get_rheobase_from_ramp(fn, chans)
    
    ramp_abf = pyabf.ABF(fn)
    ramp_dict = hcf.load_traces(fn)

    for i, ch in enumerate(chans):
        fig, ax = plt.subplots(2,1)
        swp = swps[i]
        if swp is math.nan:
            plt.close(fig)
            continue
        end = THs_in_trace[i] + 2000
        ramp = ramp_dict['Ch' + str(ch)][0][:, swp]
        ramp_abf.setSweep(sweepNumber = swp, channel = 0)

        ax[0].plot(ramp_abf.sweepX[:end], ramp[:end]) 
        ax[0].set_ylabel(ramp_abf.sweepLabelY, fontsize=10)

        ax[1].plot(ramp_abf.sweepX[:end], ramp_abf.sweepC[:end])
        ax[1].set_xlabel(ramp_abf.sweepLabelX, fontsize=10)
        ax[1].set_ylabel(ramp_abf.sweepLabelC, fontsize=10)

        ax[0].scatter(ramp_abf.sweepX[THs_in_trace[i]], THs[i], c = '#ff7f00')
        ax[1].scatter(ramp_abf.sweepX[THs_in_trace[i]], rheos[i], c = '#ff7f00')
    
        fig.suptitle(str(round(rheos[i],2)) + ' pA' + ' rheobase for ' + fn[end_fn:] + 'ch '+ str(ch) +', swp #' + str(swp), fontsize=10)
        fig.patch.set_facecolor('white')
        fig.tight_layout()

        plt.savefig(dir_rheobase + '/' + fn[end_fn:-4] + '_Ch' + str(ch) +'_rheo.pdf')
        plt.close(fig)


def plot_trace(fn, sweep, channel, save_dir = None):
    ''' 
    arguemnts : fn - filename, sweep - sweep number, channel - active channel from 1 to 8
    fast visualization of recordings and corresponding protocol
    accepts only 1 int for channel
    '''

    end_fn = fn.rfind('/') + 1
    trace = pyabf.ABF(fn)

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

    plt.style.use(['fast'])
    fig, ax = plt.subplots(2,1, sharex = False, figsize = (14,8))
    if sweep == 'all':
        data_long = trace.data[channel]
        ax[0].plot(data_long)
    else:
        trace.setSweep(sweepNumber = sweep, channel = channel)
        ax[0].plot(trace.sweepX, trace.sweepY, c = '#0000FF', linewidth = 6, alpha = 0.7)
        x = trace.sweepX
        y = trace.sweepY
    
    ax[0].set_ylabel(trace.sweepLabelY)

    ax[1].plot(trace.sweepX, trace.sweepC, linewidth = 6)
    ax[1].set_xlabel(trace.sweepLabelX)
    ax[1].set_ylabel(trace.sweepLabelC)

    fig.suptitle('{0}, sweep num {1} , channel {2}'.format(fn[end_fn:], str(sweep),channel_name))

    if save_dir:
        plt.savefig('{0}trace_{1}_{2}_swp_{3}.png'.format(save_dir,fn[end_fn:], channel_name, sweep))

    plt.show()
    #return x,y

def plot_average_all_swps(fn, chans):
    '''
    plots the average of all sweeps
    for looking at inputs
    '''
    trace = pyabf.ABF(fn)
    if '_Ipatch' in trace.adcNames:
        trace.adcNames[trace.adcNames.index('_Ipatch')] = 'Ch1'
    if 'IN0' in trace.adcNames:
        trace.adcNames[trace.adcNames.index('IN0')] = 'Ch1'

    for channel in chans:
        if len(trace.channelList) < 8:
            channel_name = 'Ch' + str(channel)
            channel = trace.channelList[trace.adcNames.index(channel_name)]
        else:
            channel = channel-1

        data_long = trace.data[channel]
        trace.setSweep(sweepNumber = 0, channel = channel)
        data_wide = data_long.reshape(int(len(data_long)/ len(trace.sweepY)), len(trace.sweepY))
        data_wide_mean = np.mean(data_wide, axis = 0)
        fig, ax = plt.subplots(1,1, figsize = (7, 4))
        x = np.array(range(len(data_wide_mean)))/trace.sampleRate
        ax.plot(x, data_wide_mean)

        ax.set_xlabel(trace.sweepLabelX) 
        ax.set_ylabel(trace.sweepLabelY)

        fig.suptitle('Average across all sweeps')
        #plt.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/mouse/ca1-sub-channorhodopsin_mice/Calb-CreSubiculum/inputs/'+ filenames[i][:-4] + 'Ch.' +channel_str + '.jpg')
        plt.show()

def plt_trace_select_swps_cut_parts(fn, chan, swps_keep = 'all', start_point_cut = False, end_point_cut = False):
    ''' 
    arguemnts : fn - filename, chan - channel number from 1 to 8
    visualization of recordings with selection of channels and parts to cut
    accepts int for channel
    start_point_cut, end_point_cut - int, what aprt to cut from the trace
    '''
    trace = pyabf.ABF(fn)
    # fixing the naming
    if '_Ipatch' in trace.adcNames:
            trace.adcNames[trace.adcNames.index('_Ipatch')] = 'Ch1'
    if 'IN0' in trace.adcNames:
        trace.adcNames[trace.adcNames.index('IN0')] = 'Ch1'
    
    # channels are always consequtive nums in pyabf
    chan_name = 'Ch' + str(chan)
    if len(trace.channelList) < 8:
        chan = trace.channelList[trace.adcNames.index(chan_name)]
    else:
        chan = chan - 1

    data_long = trace.data[chan]
    swp_len = int(trace.sweepLengthSec * trace.sampleRate)
    data_wide = data_long.reshape(trace.sweepCount, swp_len) # each sweep is a row
    
    if len(swps_keep) == trace.sweepCount:
        data_plot = hcf.reshape_data(data_wide, 'all', start_point_cut, end_point_cut)
    else:
        data_plot = hcf.reshape_data(data_wide, swps_keep, start_point_cut, end_point_cut)
    
    fig, ax = plt.subplots(1,1, figsize = (7, 4))
    x = np.linspace(0, len(data_plot), len(data_plot))/trace.sampleRate
    ax.plot(x, data_plot)
    for i in range(1, len(swps_keep)):
        vline = i * trace.sweepLengthSec
        if i == 1:
            ax.axvline(x = vline, color='red', linestyle='--', linewidth=1, label = 'sweep end')
        else:
            ax.axvline(x = vline, color='red', linestyle='--', linewidth=1)

    ax.set_xlabel(trace.sweepLabelX) 
    ax.set_ylabel(trace.sweepLabelY)
    fig.suptitle(chan_name + ' from ' + fn[fn.rfind('/') + 1:] + ' nice sweeps')
    fig.legend(loc='lower left', bbox_to_anchor=(-0.1, -0.1), fontsize = 13)
    plt.show()
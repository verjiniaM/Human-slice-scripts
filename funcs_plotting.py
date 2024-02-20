
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import stimulation_windows_ms as stim_win
import funcs_con_screen as con_param
import funcs_sorting as sort
import funcs_human_characterisation as hcf
import pyabf

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
        
                #fig,ax = plt.subplots(3,2,sharex=True, figsize = (12,6))
                fig, ax = plt.subplots(3,2,sharex=True)
                hold = ax[0,0].scatter(range(sweep_count), Res[:,0], color='b', label = 'Holding current (median)')
                #ax[0].set_title('Holding current')
                hold_win = ax[0,0].fill_between(range(sweep_count), Res[0,0] - THholding, Res[0,0] + THholding, color='b', alpha=0.1, label = '20% window')
                ax[0,0].set_ylim([np.min(Res[0,0] - 1.5 * THholding), np.max(Res[0,0] + 1.5 * THholding)])
                ax[0,0].set_yticks([np.min(Res[0,0] - 1.5 * THholding), Res[0,0], np.max(Res[0,0] + 1.5 * THholding)])
                #ax[0].legend(loc = 'upper right')

                Rs = ax[1,0].scatter(range(sweep_count), Res[:,1], color='r', label = 'Series resistance (median)')
                #ax[1].set_title('Series resistance')
                Rs_win = ax[1,0].fill_between(range(sweep_count), Res[0,1] - THRs, Res[0,1] + THRs, color='r', alpha=0.1, label = '20% window')
                ax[1,0].set_ylim([np.min(Res[0,1] - 1.5 * THRs), np.max(Res[0,1] + 1.5 * THRs)])
                ax[1,0].set_yticks([np.min(Res[0,1] - 1.3 * THRs), Res[0,1], np.max(Res[0,1] + 1.3 * THRs)])
                #ax[1].legend(loc = 'upper right')

                Rin = ax[2,0].scatter(range(sweep_count), Res[:,2], color='g', label = 'Input resistance (median)')
                #ax[2].set_title('Input resistance')
                Rin_win = ax[2,0].fill_between(range(sweep_count), Res[0,2] - THRi, Res[0,2] + THRi, color='g', alpha=0.1, label = '20% window')
                ax[2,0].set_ylim([np.min(Res[0,2] - 1.5 * THRi), np.max(Res[0,2] + 1.5 * THRi)])
                ax[2,0].set_yticks([np.min(Res[0,2] - 1.3 * THRi), Res[0,2], np.max(Res[0,2] + 1.3 * THRi)])
                ax[2,0].set_xlabel('Sweep num')
                #ax[2].legend(loc = 'upper right')
            
            #leg = plt.figlegend(loc = 'center right', bbox_to_anchor=(1.4, 0.5))
                lgd = plt.legend(
                    [hold, hold_win, Rs, Rs_win, Rin, Rin_win],
                    ['Holding current (median)','20% window', 'Series resistance (median)', '20% window', 'Input resistance (median)', '20% window'],
                    ## by default, legend anchor to axis, but can
                    ## also be anchored to arbitary position
                    ## positions within [1,1] would be within the figure
                    ## all numbers are ratio by default

                    bbox_to_anchor=(1.4, 0.7),

                    ## loc indicates the position within the figure
                    ## it is defined consistent to the same Matlab function 
                    loc='center right',

                    ncol = 1
                    #mode="expand",
                    #borderaxespad=0.
                    )
                
            fig.patch.set_facecolor('white') 
            fig.suptitle(key, fontsize = 15)
            fig.tight_layout()

            plt.savefig(dir_vc_plots + '/' + filename[end_fn:-4] + '_'+ key + '_VC_plot.png')
            plt.close(fig)


def plot_hyperpolar (filename, channels, inj, onset = 2624, offset = 22624, 
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
        fig.tight_layout()
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

#plots the sweeps with the min value and the 
def plot_mini_sweeps (filename, cell_chan, sweep):
    sweep = sweep - 1
    ch_data, sweep_len, block = hchf.load_traces(filename, cell_chan) #ch_data: each sweep is a column
    #no test pulse

    signal_no_test_pulse = ch_data[4250:,sweep]
    
    min_val = np.amin(signal_no_test_pulse)
    max_val = np.amax(signal_no_test_pulse)
    #loc_max = np.where(signal_no_test_pulse == max_val)
    loc_min = np.where(signal_no_test_pulse == min_val)
    loc_max = np.where(signal_no_test_pulse == max_val)
    num_points = np.shape(signal_no_test_pulse)[0]
    times = np.linspace(0,num_points,num_points)/20000
    plt.plot(times,signal_no_test_pulse, lw=0.5)
    plt.plot(int(loc_min[0])/20000, np.amin(signal_no_test_pulse),'ro')
    plt.plot(int(loc_max[0])/20000, np.amax(signal_no_test_pulse),'k*')
    plt.xlabel('sec')
    plt.ylabel('pA')
    plt.show()

    #plot only min_interval
    min_interval = signal_no_test_pulse[int(loc_min[0][0])-2000:int(loc_min[0][0])+2000]
    min_val_int= np.amin(min_interval)
    loc_min_int = np.where(min_interval == min_val_int)

    num_points_int = np.shape(min_interval)[0]
    times_int = np.linspace(0,num_points_int,num_points_int)/20000
    plt.plot(times_int,min_interval)
    plt.plot(int(loc_min_int[0][0])/20000, min_val_int,'ro')
    plt.xlabel('sec')
    plt.ylabel('pA')
    plt.show()

    max_interval = signal_no_test_pulse[int(loc_max[0][0])-2000:int(loc_max[0][0])+2000]
    max_val_int= np.amin(max_interval)
    loc_max_int = np.where(max_interval == max_val_int)

    num_points_int = np.shape(max_interval)[0]
    plt.plot(times_int,max_interval)
    plt.plot(int(loc_max_int[0][0])/20000, max_val_int,'ro')
    plt.xlabel('sec')
    plt.ylabel('pA')
    plt.show()


## Connectivity plotting functions 

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

def plot_connection_window(con_screen_file, preC, postC, pre_window, post_window, preAPs_shifted, postsig,\
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
        y = postsig[:,i][preAPs[0][0]-750:preAPs[0][len(preAPs)-1]+750]
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
    y = mean_post[xvals:xvals+4000]+10
    times = np.linspace(0,len(y), len(y))/(sampl_rate*1000) #foe ms
    plt.plot(times, y, lw=1.5, color='k')
    plt.xlabel('ms')
    plt.ylabel('mV')
    plt.savefig(dir_connect + '/post_swps_mean_' + con_screen_file[end_fn:-4] + '_' + 
    'Ch' + str(pre_cell_chan) + '_to_' + 'Ch' + str(post_cell_chan) + '.png')
    fig.patch.set_facecolor('white')
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
        y = postsig[:,i][preAPs[0][0]-750:preAPs[0][len(preAPs)-1]+750]
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


def plot_trace(fn, sweep, channel):
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

    fig, ax = plt.subplots(2,1, sharex = False, figsize = (14,8))
    if sweep == 'all':
        data_long = trace.data[channel]
        ax[0].plot(data_long)
    else:
        trace.setSweep(sweepNumber = sweep, channel = channel)
        ax[0].plot(trace.sweepX, trace.sweepY)
        x = trace.sweepX
        y = trace.sweepY
    
    ax[0].set_ylabel(trace.sweepLabelY)

    ax[1].plot(trace.sweepX, trace.sweepC)
    ax[1].set_xlabel(trace.sweepLabelX)
    ax[1].set_ylabel(trace.sweepLabelC)

    channel_name = 'Ch' + str(channel)
    fig.suptitle('{0}, sweep num {1} , channel {2}'.format(fn[end_fn:], str(sweep),channel_name))

    plt.show()
    #return x,y

def plot_average_all_swps(fn, chans):
    '''
    plots the average of all sweeps
    for looking at inputs
    '''
    trace = pyabf.ABF(fn)

    for channel in chans:
        channel_str = str(channel)
        if len(trace.channelList) < 8:
            if '_Ipatch' in trace.adcNames:
                trace.adcNames[trace.adcNames.index('_Ipatch')] = 'Ch1'
            if 'IN0' in trace.adcNames:
                trace.adcNames[trace.adcNames.index('IN0')] = 'Ch1'
            channel_name = 'Ch' + str(channel)
            channel = trace.channelList[trace.adcNames.index(channel_name)]
        else:
            channel = channel-1

        data_long = trace.data[channel]
        trace.setSweep(sweepNumber = 0, channel = channel)
        data_wide = data_long.reshape(int(len(data_long)/ len(trace.sweepY)), len(trace.sweepY))
        data_wide_mean = np.mean(data_wide, axis = 0)
        fig2, ax = plt.subplots(1,1, figsize = (7, 4))
        x = np.array(range(len(data_wide_mean)))/20_000
        ax.plot(x, data_wide_mean)

        ax.set_xlabel(trace.sweepLabelX, fontsize = 14) 
        ax.set_ylabel(trace.sweepLabelY, fontsize = 14)

        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        fig.suptitle('Average across all sweeps', fontsize = 15)
        #plt.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/mouse/ca1-sub-channorhodopsin_mice/Calb-CreSubiculum/inputs/'+ filenames[i][:-4] + 'Ch.' +channel_str + '.jpg')
        plt.show()

####################
#for EPSP manually analyzed


def QC_filter_for_plotting(results_df, holding, temp, K_concentr = 8, min_age = 0):
    '''
    holding [str] - 'yes' or 'no'
    temp [str] - 'yes' or 'no'
    '''
    results_df_plot = results_df[(results_df['holding_minus_70_y_o_n'] == holding) &\
        (results_df['K_concentration'] == K_concentr) &\
            (results_df['recording_in'] != 'puff high K') &\
                #(results_df['Average amplitude (pA)'] > min_avg_amp) &\ 
                    (results_df['temperature'] == temp) &
                    (results_df['patient_age'] > min_age)]

    amp_non_negative = [-i for i in results_df_plot['Average amplitude (pA)']]
    results_df_plot.insert(20, 'Average amp (positive)', amp_non_negative)
    return results_df_plot

def QC_RMP_Ctrl(df, max_allowed_RMP_Ctrl):
    '''
    checks that the RMP in Ctrl condition is not above max_allowed_RMP_Ctrl
    exludes cells from both conditions, where it is
    '''

    #
    not_repeated_cells = []
    for cell in df['cell_ID'].unique():
        #print(len(df[df['cell_ID'] == cell]))
        if len(df[df['cell_ID'] == cell]) == 1: # or len(df[df['cell_ID'] == cell]) == 3:
            not_repeated_cells.append(cell)
    for cell in not_repeated_cells:
        df = df.drop(df.index[df['cell_ID'] == cell])     
    df.reset_index(inplace = True, drop = True)

    cell_ID_keep = df['cell_ID'][(df['recording_in'] == 'Ctrl') & (df['resting_potential'] < max_allowed_RMP_Ctrl)].to_list()
    cells_to_delete = list(set(df['cell_ID'].unique().tolist()) - set(cell_ID_keep))

    for cell in cells_to_delete:
        df = df.drop(df.index[df['cell_ID'] == cell])     
    df.reset_index(inplace = True, drop = True)
    return df

#plotting funcs
#mpl.rcParams - for all parameter settings
@mpl.rc_context({'axes.labelsize': 17, \
    'axes.spines.right': False, \
        'axes.spines.top': False, \
            'axes.titlesize': 15, \
                'xtick.labelsize': 15, \
                     'ytick.labelsize': 15, \
                        'figure.titlesize': 20})
def plot_ (df, title, params = ['Average amp (positive)', 'Average interevent interval (ms)']):
    OP_colors = ['#dede00', '#ff7f00', '#4daf4a', '#984ea3', 'violet']
    op_color_dict = {'OP230914':'#dede00', 'OP231005':'#ff7f00', 'OP231109':'#4daf4a', 'OP231123':'#984ea3', 'OP231130':'violet'}
    df = df.sort_values(by = ['recording_in', 'cell_ID']).reset_index(drop= True)
    fig, ax = plt.subplots(2,1, sharex = False, figsize=(8,10))
    for p, param in enumerate(params):
        x_vals = []
        for i, rec_solution in enumerate(sorted(df.recording_in.unique())):
            df_plot = df[df['recording_in'] == rec_solution].reset_index(drop= True) 
            x = np.linspace(2*i, 1 + 2*i, len(df_plot))
            ax[p].scatter(x, df_plot[param], alpha = 0.7, label = 'Ctrl')
            ax[p].plot([0.4 + 2*i, 0.6 + 2*i], [np.nanmean(df_plot[param]), np.nanmean(df_plot[param])], c = 'k')
            ax[p].text(0.4 + 2*i, np.nanmean(df_plot[param]) + p +0.03,  str(round(np.nanmean(df_plot[param]),2)), size = 15, c = 'k', zorder = 10)
            ax[p].set_title(param)

            x_vals.append(x)

            for j, OP in enumerate(sorted(df_plot.OP.unique())):
                indx = df_plot[df_plot['OP'] == OP].index
                x_op = x[indx]
                y_op = df_plot[param][indx]
                ax[p].scatter(x_op, y_op, c = op_color_dict[OP], s = 60, zorder = 5, label = OP)
                
        cell_IDs = df['cell_ID'][df['recording_in'] == 'Ctrl'].values
        for c, cell in enumerate(cell_IDs):
            #indx = index_s[c]
            #x_K = results_df_plot['x'].loc[results_df_plot['cell_ID'] == cell].tolist()[0]
            x = [x_vals[0][c], x_vals[1][c]]
            y = [df[param][df['recording_in'] == 'Ctrl'].tolist()[c], df[param][df['recording_in'] == 'high K'].tolist()[c]]
            #op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
            ax[p].plot(x, y, '-', color = 'darkgrey', alpha = 0.5, linewidth = 2, zorder = 1)
    
    ax[0].set_ylabel('Amplitude (pA)')
    ax[1].set_ylabel('IEI (ms)')

    ax[1].set_xticks(ticks = [0.5, 2.5], labels = ['Ctrl', 'high K'],)

    fig.suptitle(title)
    fig.tight_layout()
    fig.patch.set_facecolor('white')
    #fig.legend()

    plt.show()

@mpl.rc_context({'axes.labelsize': 17, \
    'axes.spines.right': False, \
        'axes.spines.top': False, \
            'axes.titlesize': 15, \
                'xtick.labelsize': 15, \
                     'ytick.labelsize': 15, \
                        'figure.titlesize': 20})
def sns_plot_MEA_data(df, title):
    colors = ['darkblue','#4daf4a', 'violet']
    customPalette = sns.set_palette(sns.color_palette(colors))

    fig, ax1 = plt.subplots(1,1, sharex = False, figsize=(8,5))
    sns.lineplot(
        data = df, x = 'Condition', y = "value",
        hue="OP", palette = customPalette , ax = ax1)
    sns.scatterplot(
        data = df, x = 'Condition', y = "value", 
        hue="OP", palette = customPalette , ax = ax1)

    ax1.set_xticks(ax1.get_xticks(), df['Condition'].unique(), rotation=30)
    
    ax1.set_xlabel('')
    ax1.set_ylabel('Network Activity \n(spikes\electrode\second')

    fig.suptitle(title)
    fig.tight_layout()
    fig.patch.set_facecolor('white')
    #fig.legend()
    fig.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.show()


def plot_from_full_results_table(results_df_plot, title):
    median_amp_no_hold_no_temp_Ctrl = results_df_plot['Average amp (positive)'][results_df_plot['recording_in'] == 'Ctrl'].values
    median_amp_no_hold_no_temp_highK = results_df_plot['Average amp (positive)'][results_df_plot['recording_in'] == 'high K'].values
    freq_no_hold_no_temp_Ctrl = results_df_plot['Average interevent interval (ms)'][results_df_plot['recording_in'] == 'Ctrl'].values
    freq_no_hold_no_temp_highK = results_df_plot['Average interevent interval (ms)'][results_df_plot['recording_in'] == 'high K'].values

    fig, ax = plt.subplots(2,1, sharex = True, figsize=(12,10))
    x1 = np.linspace(0, 1,len(median_amp_no_hold_no_temp_Ctrl))
    x2 = np.linspace(2, 3,len(median_amp_no_hold_no_temp_highK))
    ax[0].scatter(x1, median_amp_no_hold_no_temp_Ctrl, alpha = 0.7, label = 'Ctrl')
    ax[0].scatter(x2, median_amp_no_hold_no_temp_highK, alpha = 0.7, label = 'high K')
    ax[0].plot([0.4, 0.6], [np.nanmean(median_amp_no_hold_no_temp_Ctrl), np.nanmean(median_amp_no_hold_no_temp_Ctrl)], c = 'k')
    ax[0].plot([2.4, 2.6], [np.nanmean(median_amp_no_hold_no_temp_highK), np.nanmean(median_amp_no_hold_no_temp_highK)], c = 'k')

    ax[0].text(0.5, np.nanmean(median_amp_no_hold_no_temp_Ctrl) + 0.07,  str(round(np.nanmean(median_amp_no_hold_no_temp_Ctrl),2)), size = 15, c = 'k')
    ax[0].text(2.5, np.nanmean(median_amp_no_hold_no_temp_highK) + 0.07,  str(round(np.nanmean(median_amp_no_hold_no_temp_highK),2)), size = 15, c = 'k')
    ax[0].set_title('Average EPSP amplitude for cell per condition, no APs')

    cell_IDs = results_df_plot['cell_ID'][results_df_plot['recording_in'] == 'Ctrl'].values
    for c, cell in enumerate(cell_IDs):
        #indx = index_s[c]
        #x_K = results_df_plot['x'].loc[results_df_plot['cell_ID'] == cell].tolist()[0]
        x = [x1[c],x2[c]]
        y = [median_amp_no_hold_no_temp_Ctrl[c], median_amp_no_hold_no_temp_highK[c]]
        #op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
        ax[0].plot(x, y, '-', color = 'darkgrey', alpha = 0.5, linewidth = 2, zorder = 1)

    x3 = np.linspace(0, 1,len(freq_no_hold_no_temp_Ctrl))
    x4 = np.linspace(2, 3,len(freq_no_hold_no_temp_highK))
    ax[1].scatter(x3, freq_no_hold_no_temp_Ctrl, alpha = 0.7, label = 'Ctrl')
    ax[1].scatter(x4, freq_no_hold_no_temp_highK, alpha = 0.7, label = 'high K')
    ax[1].plot([0.4, 0.6], [np.nanmean(freq_no_hold_no_temp_Ctrl), np.nanmean(freq_no_hold_no_temp_Ctrl)], c = 'k')
    ax[1].plot([2.4, 2.6], [np.nanmean(freq_no_hold_no_temp_highK), np.nanmean(freq_no_hold_no_temp_highK)], c = 'k')

    ax[1].text(0.5, np.nanmean(freq_no_hold_no_temp_Ctrl) + 0.1,  str(round(np.nanmean(freq_no_hold_no_temp_Ctrl),2)), size = 15, c = 'k')
    ax[1].text(2.5, np.nanmean(freq_no_hold_no_temp_highK) + 0.1,  str(round(np.nanmean(freq_no_hold_no_temp_highK),2)), size = 15, c = 'k')
    ax[1].set_title('Average interevent interval (ms), no APs')

    cell_IDs = results_df_plot['cell_ID'][results_df_plot['recording_in'] == 'Ctrl'].values
    for c, cell in enumerate(cell_IDs):
        #indx = index_s[c]
        #x_K = results_df_plot['x'].loc[results_df_plot['cell_ID'] == cell].tolist()[0]
        x = [x3[c],x4[c]]
        y = [freq_no_hold_no_temp_Ctrl[c], freq_no_hold_no_temp_highK[c]]
        #op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
        ax[1].plot(x, y, '-', color = 'darkgrey', alpha = 0.5, linewidth = 2, zorder = 1)

    ax[0].set_ylabel('Amplitude (pA)')
    ax[1].set_ylabel('IEI (ms)')

    ax[1].set_xticks(ticks = [0.5, 2.5], labels = ['Ctrl', 'high K'])

    fig.suptitle(title)
    fig.tight_layout()
    fig.patch.set_facecolor('white')
    fig.legend()

    plt.show()

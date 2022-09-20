
import neo
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import human_characterisation_functions as hchf
import sorting_functions as sort


#%%
# =============================================================================
# some fast plotting functions
# =============================================================================

# plots the middle sweep for each channel
# check the traces to see which channels were active or if the protocol names are enetered correctly
def plot_traces(filename): 
    end_filename = filename.rfind('/') + 1
    dir_plots = sort.make_dir_if_not_existing(filename[:end_filename], 'plots')
    dir_traces = sort.make_dir_if_not_existing(dir_plots, 'traces')

    sweep_count, sweep_len, channels, channel_dict = hchf.get_analog_signals(filename)

    plt.style.use(['fast'])

    x = math.ceil(np.sqrt(channels))
    fig = plt.figure(figsize=(6,6))
    plt.subplots_adjust(hspace=0.5)

    for i in range(1,channels+1):
        ax = plt.subplot(x,x, i)
        cell_chan = i
        #ch1 = np.ndarray([1, sweeps])  #empty matrix
        # for key, val in channel_dict.items():
        #     if val == 'Ch'+ str(cell_chan): #-1
        #         cell_chan = int(key[-1])                  
        # for i in range(0,lens(block.segments)):
        middle_swp_num = int(sweep_count/2)
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

    fig.tight_layout() 
    fig.patch.set_facecolor('white')
    plt.savefig(dir_traces + '/trace_plot_' + filename[end_filename:-4] + '.png')
    plt.close(fig)
    return 'Trace plots saved in' + dir_traces

def plot_vc_holding (filename, cell_chan):
    end_filename = filename.rfind('/') + 1

    sweep_count, sweep_len, channels, channel_dict, vcdata = hchf.get_analog_signals(filename)
    dir_vc_plots = sor.make_dir_if_not_existing(filename[:end_filename]+'plots/', "vc_plots")

    plt.style.use(['fast'])   
    
    channel = 'Ch'+ str(cell_chan)
    swps=len(vcdata[channel][0])
    Res= np.ndarray([swps,3])
    for n, trace in enumerate(vcdata[channel][0]):
        # HC=np.median(trace[5000:].view(np.recarray)) 
        # avgBL=np.median(trace[2500:3000].view(np.recarray))
        # minP=np.min(trace[3000:3500].view(np.recarray))
        # ss=np.median(trace[3700:4000].view(np.recarray))
        HC = np.median(trace[4000:].view(np.recarray)) 
        avgBL = np.median(trace[10:990].view(np.recarray))
        minP = np.min(trace[990:1300].view(np.recarray))
        ss = np.median(trace[1400:1900].view(np.recarray))
        RaI = avgBL-minP
        RiI = avgBL-ss
        if RaI==0 or RiI == 0: 
            THholding=Res[0,0]*0.1
            Res[n,0]=HC
            Res[n,1]= math.nan
            Res[n,2]= math.nan
            fig,ax=plt.subplots(1,1,sharex=True)
            ax.scatter(range(swps),Res[:,0], color='b')
            ax.set_title('Holding current')
            ax.fill_between(range(swps),Res[0,0]-THholding, Res[0,0]+THholding, color='b', alpha=0.1)
            ax.set_xlabel('Sweep num')
            fig.suptitle(channel, fontsize=15)
            fig.tight_layout()

            plt.savefig(dir_vc_plots + '/' + filename[end_filename:-4]+'_'+channel + '_VC_plot.png')
            plt.close(fig)
        else:
            ResRa = (0.004/(RaI*1e-12))/1000000
            ResRi = (0.004/(RiI*1e-12))/1000000
            Res[n,0]=HC
            Res[n,1]=ResRa
            Res[n,2]=ResRi
    
            THholding=Res[0,0]*0.1
            THRs=Res[0,1]*0.1 #10 %
            THRi=Res[0,2]*0.1
    
            fig,ax=plt.subplots(3,1,sharex=True)
            ax[0].scatter(range(swps),Res[:,0], color='b')
            ax[1].scatter(range(swps),Res[:,1], color='r')
            ax[2].scatter(range(swps),Res[:,2], color='g')
            ax[0].set_title('Holding current')
            ax[1].set_title('Series resistance')
            ax[2].set_title('Input resistance')
            ax[0].fill_between(range(swps),Res[0,0]-THholding, Res[0,0]+THholding, color='b', alpha=0.1)
            ax[1].fill_between(range(swps),Res[0,1]-THRs, Res[0,1]+THRs, color='r', alpha=0.1)
            ax[2].fill_between(range(swps),Res[0,2]-THRi, Res[0,2]+THRi, color='g', alpha=0.1)
            ax[2].set_xlabel('Sweep num')

            fig.patch.set_facecolor('white') 
            fig.suptitle(channel, fontsize=15)
            fig.tight_layout()

            plt.savefig(dir_vc_plots + '/' + filename[end_filename:-4]+'_'+channel + '_VC_plot.png')
            plt.close(fig)


def plot_hyperpolar (filename, cell_chan, ch1, dir_plots, V65, mc, onset = 2624, offset = 22624, 
clrs = ["b", "g", "r", "c", "m", "y", "#FF4500", "#800080"]):
    end_fn = filename.rfind('/') + 1
    dir_onset = sort.make_dir_if_not_existing(dir_plots, 'Onset')
    fig = plt.figure()
    for i in range(0,5):
        swp = list(ch1[:,i])
        if list(filter(lambda ii: ii < V65, swp)) == []:
            print("No hyperpolarization fig for " + filename[end_fn:-4] + 'ch: ' + str(cell_chan))
        else:
            plt.plot(ch1[:,i], c = clrs[i])
            plt.scatter(onset + tc, V65, c=clrs[i])
            plt.annotate('V65  ', (onset + tc, V65), horizontalalignment='right')
            plt.scatter(onset, bl, c='r')
    fig.patch.set_facecolor('white')    
    plt.annotate('  Baseline', (onset, bl))
    plt.title('Ch ' + str(cell_chan))
    plt.savefig(dir_onset + '/Char_onset_plot_' + filename[end_fn:-4]+'_'+str(cell_chan) + '.png')
    plt.close()
     
def plot_spikes (filename, inj, ch1, dir_plots, first_spike, peaks):
    dir_spikes = sort.make_dir_if_not_existing(dir_plots, 'Max_Spikes')
    win = len(inj) - first_spike
    x = math.ceil(np.sqrt(win))
    fig, ax = plt.subplots(x,x,sharex=True, sharey=False,figsize=(6,6))
    for i in range(win):
        ax = fig.add_subplot(x,x, i+1)
        if first_spike+i-1 < np.shape(ch1)[1]:
            ax.plot(ch1[:,first_spike+i-1], lw=0.5, c='grey')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.scatter(peaks[first_spike+i-1, :, 1], 
                    peaks[first_spike+i-1, :, 2], marker='+', c='r')
            ax.annotate('+'+str(inj[first_spike+i-1])+' pA', (25500,0), (25500,0), color='b', rotation=90)
    
    fig.patch.set_facecolor('white')
    fig.suptitle('Ch ' +str(cell_chan), fontsize=15)
    fig.tight_layout()
    plt.savefig(dir_spikes + '/char_spikes_plot_' + filename[end_fn:-4]+'_'+str(cell_chan) + '.png')
    plt.close(fig)

def plot_iv_curve (filename, cell_chan, inj, dir_plots, first_spike, spike_counts):
    dir_iv_curve = sort.make_dir_if_not_existing(dir_plots, "IV_curve")
    IOfit = np.polyfit(spike_counts[first_spike:,0], spike_counts[first_spike:,1],1)
    IO_slope = IOfit[0]
    fig = plt.figure()
    plt.plot(spike_counts[:,0], spike_counts[:,1])
    plt.gca().set_xlabel('Injection current, pA')
    plt.gca().set_ylabel('Number of spikes')
    Rheobase = inj[first_spike]   
    
    fig.suptitle('Ch ' + str(cell_chan), fontsize=15)
    fig.patch.set_facecolor('white')
    plt.savefig(dir_iv_curve+ '/char_IV_curve_' + filename[end_fn:-4]+'_'+str(cell_chan) + '.png')
    plt.close()

def plot_ap_props(filename, cell_chan, dir_plots, AP, THloc, TH, APheight, first_spiking_sweep):
    dir_plots = sort.make_dir_if_not_existing(dir_plots, "AP_props")

    #if more than 1 AP fired, characterize the 2nd
    if np.max(spike_counts[:,1]) == 1:
        ap = 0
    else:
        ap = 1

    fig = plt.figure()
    plt.plot(AP)
    plt.scatter(THloc,TH, color = 'red')
    plt.scatter(200,peaks[first_spiking_sweep, ap, 2], color = 'green')
    fig.suptitle('Ch: ' + str(cell_chan) + ', AP#' + str(ap+1) + ', TH = ' + str(round(TH,2)) + ', amp = ' + str(round(APheight,2)))
    fig.patch.set_facecolor('white') 
    plt.savefig(dir_plots + '/' + filename[end_fn:-4] + '_AP#' + str(ap+1) + '_' + str(cell_chan) + '.png')
    plt.close()

#plots the sweeps with the min value and the 
def plot_mini_sweeps (filename, cell_chan, sweep):
    sweep = sweep - 1
    ch1, sweep_len, block = hchf.load_traces(filename, cell_chan) #ch1: each sweep is a column
    #no test pulse

    signal_no_test_pulse = ch1[4250:,sweep]
    
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

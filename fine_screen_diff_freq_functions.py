
import interpFL
import matplotlib.pyplot as plt
import neo
import numpy as np
from scipy.signal import find_peaks
import stimulation_windows_ms as stim_win
import os


def get_pre_signal(filename, cellchannel):
    data=neo.io.AxonIO(filename)    
    b1=data.read(signal_group_mode = 'split-all')[0]

    chan_name = 'Ch' + str(cellchannel)

    #finding the correct signal for ch 
    for i in range(len(b1.segments[1].analogsignals)):
        if b1.segments[0].analogsignals[i].name == chan_name:
            sig=i
        if b1.segments[0].analogsignals[i].name == '_Ipatch' and chan_name == 'Ch1':
            sig = i
        if b1.segments[0].analogsignals[i].name == 'IN0' and chan_name == 'Ch1':
            sig = i

    signal=[]
    for s in range(len(b1.segments)):
        signal.append(b1.segments[s].analogsignals[sig])

    vmO = np.mean(signal[0][0:950]) #original Vm from first sweep
    return vmO, signal 


def get_post_signal(filename, cellchannel):
    data=neo.io.AxonIO(filename)
    b1=data.read(signal_group_mode = 'split-all')[0]

    chan_name = 'Ch' + str(cellchannel)
    for i in range(len(b1.segments[1].analogsignals)):
        if b1.segments[0].analogsignals[i].name == chan_name:
            sig=i
        if b1.segments[0].analogsignals[i].name == '_Ipatch' and chan_name == 'Ch1':
            sig = i

    signal=[]
    for s in range(len(b1.segments)):
        signal.append(b1.segments[s].analogsignals[sig])
    pA0 = np.mean(signal[0][0:950]) #original Vm from first sweep
    return pA0, signal

def get_analysis_window(presig, postsig, hz):
    size_pre = np.shape(presig)[1]
    mean_pre = (np.mean(presig,axis=0)).reshape(size_pre)
    l=len(presig[0])

    stims = stim_win.stim_window_diff_freq
    win_start = stims[hz][0] #stimulation end
    win_end = stims[hz][1] #stimulation start

    mean_pre = mean_pre[win_start:win_end]
    preAPs=find_peaks(mean_pre[win_start:win_end], height=0)
    j = 0
    while len(preAPs[0]) == 0:
        j = j-5
        preAPs=find_peaks(mean_pre, height=j)
    print("Pre APs found with heigth " + str(j))
    num_aps = len(preAPs[0])

    pre_window = mean_pre[preAPs[0][0]-30:preAPs[0][num_aps-1]+300]
    
    size_post = np.shape(postsig)[1]
    mean_post = (np.mean(postsig, axis=0)).reshape(size_post)
  
    post_window = mean_post[preAPs[0][0]-30:preAPs[0][num_aps-1]+3000] 
    preAPs_shifted = find_peaks(pre_window, height=j)

    return mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs




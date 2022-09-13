#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 14:58:20 2021

@author: rosie
"""
import interpFL
import matplotlib.pyplot as plt
import neo
import numpy as np
from scipy.signal import find_peaks
import stimulation_windows_ms
import os

def plot_connect (fn, active_channels):
    clrs = ["b", "g", "k", "c", "k", "y", "#FF4500", "#800080"]
    #plot D1
    x = len(active_channels)

    z1=0.5
    z2=40.5
    stim_window = stimulation_windows_ms.stim_window #what are the stim windows???
    conscreen=neo.io.AxonIO(fn) #loading the .abf = axon binary file 
    b1=conscreen.read(signal_group_mode = 'split-all')[0] #reading it; the file consists of 50 sweeps (segments), each 10 s long
    x = len(b1.groups) 
    fig, ax = plt.subplots(x,x,sharex=True, sharey=False,figsize=(6,6))

    for i in range(len(b1.groups)): #for each active channel
        avg = np.mean(b1.groups[i].analogsignals, axis=0)
        ch_name_i = b1.groups[i].analogsignals[0].name
        if ch_name_i == '_Ipatch': ch_name_i = 'Ch1'
        for j in range(0, ax.shape[1]):
            if i == j:
                ax[i,j].plot()
            ch_name_j = b1.groups[j].analogsignals[0].name
            if ch_name_j == '_Ipatch': ch_name_j = 'Ch1'
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
                ax[i,j].set_ylim([plotwin.min()-z1, plotwin.min()+z1])
                v1 = ax[i,j].vlines(0,plotwin.min()+0.5, plotwin.min()+1, lw=0.2, color='k') 
            else:
                ax[i,j].set_ylim([plotwin.min()-0.5, plotwin.min()+z2])
                v2 = ax[i,j].vlines(0,plotwin.min()+0.5, plotwin.min()+40.5, lw=0.2, color='k')
    
    filename = fn
    end_filename = filename.rfind('/')+1
    
    fig.suptitle('connections in ' + filename[end_filename:],fontsize=15)
    fig.patch.set_facecolor('white')
    fig.tight_layout()

    dir_plots = "plots"
    path1 = os.path.join(filename[:end_filename], dir_plots)
    if os.path.isdir(path1) == False: os.mkdir(path1)
    path = os.path.join(filename[:end_filename],dir_plots + '/connectivity_plots')
    if os.path.isdir(path) == False: os.mkdir(path)

    plt.savefig(path + '/' + filename[end_filename:-4] + '.png')
    plt.close(fig)
    
    #plt.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/Human_tissue_meetings/' 
    # + fn[-12:-4]+'_connectivity.png')
    # plt.close()

#loads presynaptic cell and quality checks the sweeps
#threshold for acceptable Rm can be set within
def presynaptic_screen(filename,cellchannel):
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
    swps = len(b1.segments)
    vmO = np.mean(signal[0][0:4000]) #original Vm from first sweep
    es = []
    for i in range(0,swps): #intrasweep control for excessively depolarised
    #cells (above -50mV) or a drift in Vm of more than 10% from start to end of sweep
        vm1 = np.mean(signal[i][0:4000])
        vm2 = np.mean(signal[i][197999:199999])
        max_val = np.max(signal[i])
        if (vm1 > -50) or vm1-(vm1*-0.1)>vm2 or vm2>vm1+(vm1*-0.1):
            es.append(i)
            print("Excluding swp # " + str(i) + '; drift more than 0.1*RMP start to end or RMP > -50') 
        elif vmO-(vmO*-0.1)>vm1 or vm1>vmO+(vmO*-0.1):
            es.append(i)
            print("Excluding swp # " + str(i) + '; drift more than 0.1* RMP') 
        #this statement accounts for APs that reverse below 0
        elif max_val<0:
            es.append(i)
            print("Excluding swp # " + str(i) + '; APs < 0') 
    es.reverse()
    exclude = 'no'
    if len(es) == swps:
        print('All pre sweeps need to be excluded')
        es = []
        print('changing the es = [], to continue analysis')
        exclude = 'yes'
    for failedswp in es:
        del(signal[failedswp])
    
    return signal, es, exclude

#loads postsynaptic cell and quality checks sweeps. Automatically removes
#sweeps excluded from presynaptic screen (presynaptic function should
#therefore be run first)
def postsynaptic_screen(filename,cellchannel,es):
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
    for failedswp in es:
        del(signal[failedswp])
    swps = len(signal)
    vmO = np.mean(signal[0][0:4000]) #original Vm from first sweep
    es2 = []
    for i in range(0,swps): #intrasweep control for excessively depolarised
    #cells (above -50mV) or a drift in Vm of more than 10% from start to end of sweep
        vm1 = np.mean(signal[i][0:4000])
        vm2 = np.mean(signal[i][197999:199999])
        if (vm1 > -50) or vm1-(vm1*-0.1)>vm2 or vm2>vm1+(vm1*-0.1):
            es2.append(i)
        elif vmO-(vmO*-0.1)>vm1 or vm1>vmO+(vmO*-0.1):
            es2.append(i)
    es2.reverse()
    
    for failedswp in es2:
        del(signal[failedswp])
    
    return signal

#cuts out stimulation window for analysis
def get_analysis_window(presig, postsig):
    
    pre_window = mean_pre[preAPs[0][0]-750:preAPs[0][num_aps-1]+750]
    mean_post = (np.mean(postsig, axis=0)).reshape(200000)
    post_window = mean_post[preAPs[0][0]-750:preAPs[0][num_aps-1]+750] 
    preAPs_shifted = find_peaks(pre_window, height=j, distance=800) #shifts so that stim_sindow starts at 0

    return mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs

#removes sweeps where  POSTsynaptic cell fires APs during stimulation window
def remove_sweeps_with_POSTsynAPs(presig, postsig, preAPs):
    remove_sweeps=[]
    num_aps = len(preAPs[0])
    for sw in range(0,len(postsig)):
        if np.max(postsig[sw][preAPs[0][0]-750:preAPs[0][num_aps-1]+750])>0:
            remove_sweeps.append(sw)
            remove_sweeps.reverse()
    for val in remove_sweeps:
        del(postsig[val])
    for val in remove_sweeps:
        del(presig[val])
    
    return presig, postsig


#takes differential and uses this to find psp peaks
def find_postsynaptic_peaks(post_window, preAPs_shifted):
    post_d1 = np.diff(post_window,1)
    post_d1_peaks = find_peaks(post_d1, distance=800) 
    PSPs=np.ndarray([4,2])
    num_aps = len(preAPs_shifted)
    for p in range(0,num_aps):
        PSPs[p,0] = np.max(post_window[preAPs_shifted[0][p]+20:\
                           preAPs_shifted[0][p]+150])
        PSPs[p,1] = np.argmax(post_window[preAPs_shifted[0][p]+20:\
                              preAPs_shifted[0][p]+150])+preAPs_shifted[0][p]+20

    return PSPs


#gets local baseline before each psp
def get_psp_baselines(post_window,preAPs_shifted):
    num_aps = len(preAPs_shifted[0])
    bl = np.ndarray([num_aps,1])
    for h in range(num_aps):
        bl[h,0] = np.mean(post_window[preAPs_shifted[0][h]-70:\
                        preAPs_shifted[0][h]-30],axis=0)
    return bl


#calculates onset according to intersect of 20-80 slope with local baseline
def get_onsets(preAPs_shifted, post_window, PSPs, bl):
    
    num_aps = len(preAPs_shifted[0])
    onsets = np.ndarray([num_aps,1])
    for i in range(num_aps): #for each event
        PSP_window = post_window[preAPs_shifted[0][i]-5:\
                                 preAPs_shifted[0][i]+150] #more precise adjustment; smaller window
        #calc diff bewteen bl and peak
        delta = PSPs[i][0] - bl[i] #event size 
        x_20 = (delta*0.2 + bl[i][0]).item()
        x_80 = (delta*0.8 + bl[i][0]).item()
        fit_start = interpFL.findLevels(PSP_window, x_20, mode='rising')
        fit_end = interpFL.findLevels(PSP_window ,x_80, mode='rising')
        x = range(155)
        x1 = x[fit_start[0][0]:fit_end[0][0]]
        fit1 = np.polyfit(x1, PSP_window[x1], 1)
        foot = (bl[i][0]-fit1[1])/fit1[0] #where the fit crosses the local baseline
        onsets[i,0] = int(foot) + preAPs_shifted[0][i]-5
        
    return onsets

# x_start = np.where(PSP_window > x_20)[0][0]
# x_end = np.where(PSP_window > x_80)[0][0]
# plt.plot(PSP_window, color='k')
# plt.scatter(x_start, x_20, marker = 'o', color = 'b')
# plt.scatter(x_end, x_80, marker = 'o', color = 'g')
# plt.scatter(onsets[0], bl[0], marker='^', color='r')

# trendpoly = np.poly1d(fit1)  
# plt.plot(x1,trendpoly(x1))

#uses onset to calculate latency defined as time between presynaptic AP peak
#and foot of PSP (calculated as explained above in get_onsets)
def latencies(onsets, preAPs_shifted):
    latency = np.ndarray([4,1])
    num_aps = len(preAPs_shifted[0])
    for i in range(num_aps):
        latency[i,0] = int(onsets[i].item()) - preAPs_shifted[0][i]
    latency = latency/20
    
    return latency

#%%

#plots connections and presynaptic APs, marks AP peaks, PSP peaks,
#onset and baseline
def plot_connection_window(con_screen, preC, postC, pre_window, post_window, preAPs_shifted, postsig,\
                           onsets, preAPs, PSPs, bl):
    
    fig,axarr = plt.subplots(2,1, sharex=True, figsize=(12,12))
    fig.patch.set_facecolor('white')

    #plost pre APs 0 shifted
    axarr[0].plot(pre_window, color='k')
    for i in range(0,4):
        axarr[0].scatter(preAPs_shifted[0][i],\
             pre_window[preAPs_shifted[0][i]], marker='o', color='r')

    #plot all post signals 0 shifted 
    for i in range(len(postsig)):
        axarr[1].plot(postsig[i][preAPs[0][0]-750:\
                      preAPs[0][3]+750].view(np.recarray),\
    lw=0.2, color='grey', alpha=0.5)
    
    axarr[1].plot(post_window, color='k', lw=0.5) #post_window is the averaged 0 shifter post signal 
    #plot the baseline meausre 
    bl1 = axarr[1].hlines(bl[0,0],preAPs_shifted[0][0]-70,\
               preAPs_shifted[0][0]-30, lw=5, color='b')
    bl2 = axarr[1].hlines(bl[1,0],preAPs_shifted[0][1]-70,\
               preAPs_shifted[0][1]-30, lw=5, color='b')
    bl3 = axarr[1].hlines(bl[2,0],preAPs_shifted[0][2]-70,\
               preAPs_shifted[0][2]-30, lw=5, color='b')
    bl4 = axarr[1].hlines(bl[3,0],preAPs_shifted[0][3]-70,\
               preAPs_shifted[0][3]-30, lw=5, color='b')
    
    axarr[1].scatter(onsets, bl, marker='^', color='r')
    
    #plot peaks
    for i in range(0,4):
        axarr[1].scatter(PSPs[i][1], PSPs[i][0], marker='+',\
             color='g', s=30,linewidth=5)
    plt.close()

    # fig.subplots_adjust(left=0.2,
    #                 bottom=0.2, 
    #                 right=0.4, 
    #                 top=0.4, 
    #                 wspace=0.09, 
    #                 hspace=0.05)
    #plt.Axes(fig, [0., 0., 1., 1.])

    end_filename = con_screen.rfind('/')+1

    plt.savefig(con_screen[:end_filename] + 'plots/connectivity_plots/pre_post_events' + con_screen[end_filename:-4] + '_Ch' + str(preC) + '#Ch' + str(postC) + '.png')
    return fig, axarr
# %%

def plot_post_cell(screenfile, pre_cell, post_cell):
    preC = 'Ch' + str(pre_cell)
    postC = 'Ch' + str(post_cell)

    stim_window = stimulation_windows_ms.stim_window
    conscreen=neo.io.AxonIO(screenfile)
    b1=conscreen.read(signal_group_mode = 'split-all')[0]
    vm_screen={}
    for ch in range(len(b1.segments[0].analogsignals)):
        signal=[]
        for s in range(len(b1.segments)):
            signal.append(b1.segments[s].analogsignals[ch])
            name = b1.segments[0].analogsignals[ch].name
        vm_screen[name]=[signal,b1.segments[s].analogsignals[ch].annotations.values()]  #channel is zero-indexed so headstage1 == channel 0 etc etc
    
        if '_Ipatch' in vm_screen:
            vm_screen['Ch1'] = vm_screen.pop('_Ipatch')
    del ch, s, signal, name

    #%% plot sweeps with y offset
    mean_pre = (np.mean(vm_screen[preC][0],axis=0)).reshape(200000)
    mean_post = (np.mean(vm_screen[postC][0],axis=0)).reshape(200000)
    preAPs=find_peaks(mean_pre, height=0, distance=800)

    j = 0 #starting with min AP heigth 0
    while len(preAPs[0]) == 0:
        j = j-5
        preAPs=find_peaks(mean_pre, height=j, distance=800)
    print("Pre APs found with heigth " + str(j))

    filename = screenfile
    dir_plots = "plots"
    end_filename = filename.rfind('/')+1
    path1 = os.path.join(filename[:end_filename], dir_plots)
    if os.path.isdir(path1) == False: os.mkdir(path1)
    path = os.path.join(filename[:end_filename],dir_plots + '/connectivity_plots')
    if os.path.isdir(path) == False: os.mkdir(path)

    fig = plt.figure(figsize=(12,12))
    fig.patch.set_facecolor('white')
    plt.subplots(1,1,sharex=True)

    xvals = preAPs[0][0]-100
    for i in range(len(vm_screen[postC][0])):
        plt.plot(vm_screen[postC][0][i][xvals:xvals+4000].view(np.recarray)-i*3, color='k', lw = 0.5)
    #    plt.vlines(preAPs[0][0]-xvals,-50, -210, lw=0.25)
    #    plt.vlines(preAPs[0][1]-xvals,-50, -210, lw=0.25)
        plt.vlines(preAPs[0]-xvals,-50, -180, lw=0.25)

    plt.close()
    plt.savefig(path + '/post_swps_' + filename[end_filename:-4] + '_' + preC + '#' + postC + '.png')
    

    fig = plt.figure(figsize=(6,6))
    fig.patch.set_facecolor('white')
    plt.plot(mean_post[xvals:xvals+4000]+10, lw=1.5, color='k')
    plt.savefig(path + '/avg_post_sweeps' + filename[end_filename:-4] + '_' + preC + '#' + postC + '.png')
    plt.close()


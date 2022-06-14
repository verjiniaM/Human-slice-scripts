#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 15:08:24 2021

@author: rosie
"""

import neo
import numpy as np
import matplotlib.pyplot as plt
from detect_peaks import detect_peaks
import os
import math
#%%
# =============================================================================
# quality check for VC traces for mini & spon recordings   
# =============================================================================


# loads VC protocol file, outputs plot of holding current, series and input
# resistance across timeframe of recording - each point refers to single sweep
def get_holding_measures(filename, cell_chan):
    data=neo.io.AxonIO(filename)
    b1=data.read(signal_group_mode = 'split-all')[0]
    # create a dictionoary with keys ch names and value ch signal 
    # (name, annotation, sampling rate, time)
    vcdata={}
    for ch in range(len(b1.segments[0].analogsignals)):
        signal=[]
        for s in range(len(b1.segments)):
            signal.append(b1.segments[s].analogsignals[ch]) #channel signal
            name = b1.segments[0].analogsignals[ch].name #channel name
        vcdata[name]=[signal,b1.segments[s].analogsignals[ch].annotations.values()]  #channel is zero-indexed so headstage1 == channel 0 etc etc     
        if '_Ipatch' in vcdata:
            vcdata['Ch1'] = vcdata.pop('_Ipatch')
    
    dir_vc_plots = "vc_plots"
    end_filename = filename.rfind('/')+1
    path = os.path.join(filename[:end_filename]+'plots/', dir_vc_plots)
    if os.path.isdir(path) == False: os.mkdir(path)

    plt.style.use(['fast'])   
    
    channel = 'Ch'+ str(cell_chan)
    swps=len(vcdata[channel][0])
    Res= np.ndarray([swps,3])
    for n, trace in enumerate(vcdata[channel][0]):
        # HC=np.median(trace[5000:].view(np.recarray)) 
        # avgBL=np.median(trace[2500:3000].view(np.recarray))
        # minP=np.min(trace[3000:3500].view(np.recarray))
        # ss=np.median(trace[3700:4000].view(np.recarray))
        HC=np.median(trace[4000:].view(np.recarray)) 
        avgBL=np.median(trace[10:990].view(np.recarray))
        minP=np.min(trace[990:1300].view(np.recarray))
        ss=np.median(trace[1400:1900].view(np.recarray))
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

            plt.savefig(path + '/' + filename[end_filename:-4]+'_'+channel + '_VC_plot.png')
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

            plt.savefig(path + '/' + filename[end_filename:-4]+'_'+channel + '_VC_plot.png')
            plt.close(fig)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 14:33:18 2021

@author: rosie
"""
#%%
from cmath import sin
from tokenize import single_quoted
import neo
import numpy as np
import matplotlib.pyplot as plt
from detect_peaks import detect_peaks
#import grid_strategy as grds
import math
from matplotlib import gridspec
import os


#%%
# =============================================================================
# Basic functions for characterisation of human neurons
# =============================================================================

# loading traces
def load_traces(filename,cell_chan): 
    r = neo.io.AxonIO(filename=filename)
    block=r.read(signal_group_mode = 'split-all')[0]
    sweeps = len(block.segments) #number of sweeps
    sweep_len = len(block.segments[0].analogsignals[0])
    channels = len(block.segments[0].analogsignals) # no of channels in recording

    #create a dictionary with channel names 
    channel_dict = {}
    for ch in range(0, channels):
        name = block.segments[0].analogsignals[ch].name
        if name == '_Ipatch' or name == 'IN0':
            name = 'Ch1'
            block.segments[0].analogsignals[ch].name = 'Ch1'
        channel_dict['AnalogSig%d' % (ch)] = name

    #fixing the indexing of channels 
    for key, val in channel_dict.items():
        if val == 'Ch'+ str(cell_chan): #-1
            cell = int(key[-1])  

    #create a matrix where each sweep is a column
    ch1 = np.ndarray([sweep_len, sweeps])                
    for i in range(0,len(block.segments)):
        ch1[:,i] = block.segments[i].analogsignals[cell].view(np.recarray).reshape(sweep_len)
  

    return ch1, sweep_len, block


#series and input resistance from vc file
def access_resistance(vctpfile, channel):
    data=neo.io.AxonIO(vctpfile)
    b1o=data.read(signal_group_mode = 'split-all')[0]

    #creates a dictionary with channel name and signal
    chan={}
    for ch in range(len(b1o.segments[0].analogsignals)):
        signal=[]
        for s in range(len(b1o.segments)):
            signal.append(b1o.segments[s].analogsignals[ch])
            numb=ch+1
            chan['ch%d'%numb]=[signal,b1o.segments[s].analogsignals[ch].annotations.values()]    #channel is zero-indexed so headstage1 == channel 0 etc etc
    del ch,numb,s,signal

    key = 'ch'+str(channel)
    mean = np.mean(chan[key][0], axis=0)
    #changed
    mean = mean.reshape(np.shape(mean)[0])
    detect_peaks(mean,mpd = 999, edge='both', valley=True) #mpd - peaks separated by min peak distance
    minpeak = np.min(mean)
    avgBL = np.mean(mean[0:1000]) #baseline average
    #SS = steady state
    avgSS = np.mean(mean[np.argmax(mean)-500:np.argmax(mean)-50]) #argmax - the index of the max value
    RaI = avgBL-minpeak
    RiI = avgBL-avgSS
    ResRa = (0.004/(RaI*1e-12))/1000000 #0.004 - size of the step, Resistance in megaOhm
    ResRi = (0.004/(RiI*1e-12))/1000000
    # print(ResRa, ResRi)
    
    return ResRa, ResRi


# reads in characterisation file and outputs intrinsic and spiking properties
#input filename, cell channel to be analyzed and the characterization injection used
#if inj = "full", the default injection used in frequenz_analyse is used
#otherwise, a specified injection can be added (e.g. inj = [-300, -200, -150,-100,-50])
def APproperties(filename, cell_chan, inj):
    ch1, sweep_len, block = load_traces(filename, cell_chan) #in ch1 each sweep is a column
    if inj == "full":
        inj=[-300,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400,450,500,550,
         600,700,800,900,1000,1100,1200,1300,1400]
    else:
        inj = inj

    #start and end of the pulses
    onset=2624
    offset=22624
    mc = np.ndarray([5,3])

    #creating folder to save the plots
    dir_plots = "Onset"
    end_filename = filename.rfind('/')+1
    path = os.path.join(filename[:end_filename]+'plots/', dir_plots)
    if os.path.isdir(path) == False: os.mkdir(path)

    # plotting hyperpolarization
    clrs = ["b", "g", "r", "c", "m", "y", "#FF4500", "#800080"]
    fig = plt.figure()
    for i in range(0,5):
        I = inj[i]*1e-12 #check step size
        bl = np.median(ch1[0:onset-20,i])
        ss = np.median(ch1[offset-2000:offset-1,i]) #steady state, during the current step
        swp = list(ch1[:,i])
        Vdiff = bl-ss #voltage deflection size
        v65 = Vdiff*0.63
        V65 = bl-v65
        if list(filter(lambda ii: ii < V65, swp)) == []:
            print("No hyperpolarization fig for " + filename[end_filename:-4] + 'ch: ' + str(cell_chan))
        else:
            res = list(filter(lambda ii: ii < V65, swp))[0] #takes the first value in swp < V65
            tau65 = swp.index(res) #index of res
            R = (Vdiff/1000)/-I     
            tc = tau65 - onset
            mc[i,0] = tc*0.05 #membranec capacitance; tc - time constant
            mc[i,1] = R*1e-6 #resistance
            mc[i,2] = tc*5e-5/R #capacitance
            plt.plot(ch1[:,i], c=clrs[i])
            plt.scatter(onset+tc, V65, c=clrs[i])
            plt.annotate('V65  ', (onset+tc, V65), horizontalalignment='right')
            plt.scatter(onset,bl,c='r')

    fig.patch.set_facecolor('white')    
    plt.annotate('  Baseline', (onset,bl))
    plt.title('Ch ' + str(cell_chan))
    plt.savefig(path + '/Char_onset_plot_' + filename[end_filename:-4]+'_'+str(cell_chan) + '.png')
    plt.close()
    mc[:,2] = mc[:,2]/1e-12  
    capacitance=mc[1,2]
    tau=mc[1,0]    

    #calculating resting membrane potential (RMP)
    Vm = np.median(ch1[:,5]) #when the inj = 0mV

    #maximum number of spikes
    max_spikes = 0
    for i, j in enumerate(range(0, len(ch1[0]))): #loop through all swps
        pks = detect_peaks(ch1[:,j], mph=0,mpd=50) # detects the peaks for each timepoint? 
        if len(pks)> max_spikes:
            max_spikes = len(pks) #find the max number of spikes (peaks)
    if max_spikes == 0: 
        print("No spikes found for" + filename[end_filename:-4] + ' Ch: ' + str(cell_chan))
        TH = math.nan
        max_depol = math.nan
        max_repol = math.nan
        APheight = math.nan
        Rheobase = math.nan
        return Vm, max_spikes, Rheobase,APheight, max_depol, max_repol, TH, capacitance, tau
    else:
        spike_counts = np.ndarray([len(inj),2])
        peaks = np.empty([len(inj), max_spikes, 3])
        peaks.fill(np.nan)
        for i, j in enumerate(range(0, len(ch1[0]))): #forr all swps
            pks = detect_peaks(ch1[:,j], mph=0,mpd=50) 
            peaks[i,0:len(pks), 0] = inj[i] #injected current
            peaks[i,0:len(pks), 1] = pks #
            peaks[i,0:len(pks), 2] = ch1[pks,j] #sweep number
        # number of spikes in each step
        for i, j in enumerate(inj):
            spike_counts[i,0] = j
            spike_counts[i,1] = (np.sum(np.isfinite(peaks[i,:,:])))/3
        
        spikes = np.where(np.isfinite(peaks)) 
        
        first_spike = spikes[0][0]
        if np.max(spike_counts[:,1]) == 1:
            print('MAX number of AP = 1 for ' + filename[end_filename:-4] + ' Ch ' + str(cell_chan))
            first_spiking_sweep = np.where(spike_counts[:,1]==1)[0][0]
        else:
            first_spiking_sweep=np.where(spike_counts[:,1]>1)[0][0] #where there is more than 1 AP

        dir_plots = "Max_Spikes"
        path = os.path.join(filename[:end_filename]+'plots/', dir_plots)
        if os.path.isdir(path) == False: os.mkdir(path)

        #plots all of the sweeps with APs, each of the detected spikes are marked with a cross
        win = len(inj) - first_spike
        x = math.ceil(np.sqrt(win))
        fig, ax = plt.subplots(x,x,sharex=True, sharey=False,figsize=(6,6))
        for i in range(0,win):
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
        plt.savefig(path + '/char_spikes_plot' + filename[end_filename:-4]+'_'+str(cell_chan) + '.png')
        plt.close(fig)

        #IV curve 
        dir_plots = "IV_curve"
        path = os.path.join(filename[:end_filename]+'plots/', dir_plots)
        if os.path.isdir(path) == False: os.mkdir(path)

        #IV cuve; input output function of the cell
        IOfit = np.polyfit(spike_counts[first_spike:,0], spike_counts[first_spike:,1],1)
        IO_slope = IOfit[0]
        plt.figure()
        plt.plot(spike_counts[:,0], spike_counts[:,1])
        plt.gca().set_xlabel('Injection current, pA')
        plt.gca().set_ylabel('Number of spikes')
        Rheobase = inj[first_spike]   
        
        fig.suptitle('Ch ' +str(cell_chan), fontsize=15)
        fig.patch.set_facecolor('white')
        plt.savefig(path + '/char_IV_curve_' + filename[end_filename:-4]+'_'+str(cell_chan) + '.png')
        plt.close()
        
        dir_plots = "AP_props"
        path = os.path.join(filename[:end_filename]+'plots/', dir_plots)
        if os.path.isdir(path) == False: os.mkdir(path)
        #calc TH from second AP in first spiking sweep. TH at vm when dV/dt > 10 mV/ms   
        # TH = threshold;
        
        if np.max(spike_counts[:,1]) == 1:
            peak_loc = np.where(ch1[:,first_spiking_sweep] == peaks[first_spiking_sweep,0,2])[0][0]
            AP1 = ch1[:, first_spiking_sweep][peak_loc - 200:peak_loc+200]
            d1_AP1 = np.diff(AP1)*20
            THloc_all = np.where(d1_AP1[:195] > 10)
            if THloc_all[0].size == 0:
                print("Only 1 SLOW AP found for " + filename[end_filename:-4] + ' Ch : ' + str(cell_chan))
                TH = math.nan
                max_depol = math.nan
                max_repol = math.nan
                APheight = math.nan
                Rheobase = math.nan
                return Vm, max_spikes, Rheobase,APheight, max_depol, max_repol, TH, capacitance, tau
            
            THloc = THloc_all[0][0]
            TH=AP1[THloc-1]
            APheight=peaks[first_spiking_sweep,0,2]-TH #amplitude of AP, the 2nd AP
            max_depol=np.max(d1_AP1[:200]) #how quickly is the voltage change occuring
            max_repol=np.min(d1_AP1[-200:])

            fig = plt.figure()
            plt.plot(AP1)
            plt.scatter(THloc,TH, color = 'red')
            plt.scatter(200,peaks[first_spiking_sweep,0,2], color = 'green')
            fig.suptitle('Ch: ' + str(cell_chan) + ', AP#1' +', TH = ' + str(round(TH,2)) + ', amp = ' + str(round(APheight,2)))
            fig.patch.set_facecolor('white') 
            plt.savefig(path + '/' + filename[end_filename:-4]+'_AP1_'+str(cell_chan) + '.png')
            plt.close()
        else:
            peak_loc = np.where(ch1[:,first_spiking_sweep] == peaks[first_spiking_sweep,1,2])[0][0]
            AP2 = ch1[:, first_spiking_sweep][peak_loc - 200:peak_loc+200]
            d1_AP2 = np.diff(AP2)*20
            THloc_all = np.where(d1_AP2[:195] > 10)
            if THloc_all[0].size == 0:
                print("No TH for SLOW AP in " + filename[end_filename:-4] + ' Ch : ' + str(cell_chan))
                TH = math.nan
                max_depol = math.nan
                max_repol = math.nan
                APheight = math.nan
                Rheobase = math.nan
                return Vm, max_spikes, Rheobase,APheight, max_depol, max_repol, TH, capacitance, tau
            THloc = THloc_all[0][0]
            TH=AP2[THloc-1]
            APheight=peaks[first_spiking_sweep,1,2]-TH #amplitude of AP, the 2nd AP
            max_depol=np.max(d1_AP2[:200]) #how quickly is the voltage change occuring
            max_repol=np.min(d1_AP2[-200:])

            fig = plt.figure()
            plt.plot(AP2)
            plt.scatter(THloc,TH, color = 'red')
            plt.scatter(200,peaks[first_spiking_sweep,1,2], color = 'green')
            fig.suptitle('Ch: ' + str(cell_chan) + ', AP#2' + ', TH = ' + str(round(TH,2)) + ', amp = ' + str(round(APheight,2)))
            fig.patch.set_facecolor('white') 
            plt.savefig(path + '/AP2_' + filename[end_filename:-4]+'_'+str(cell_chan) + '.png')
            plt.close()


        return Vm, max_spikes, Rheobase,APheight, max_depol, max_repol, TH, capacitance, tau


# calculates vm from 5th sweep of characterisation protocol (ie. 0pA step)
#lookes at membrane potential during characterization, no 
def vm(ch1):
    Vm = np.median(ch1[:,5])
    return Vm

# alternative to above for calculating Vm - now often a vm sweep is taken
# before holding current is added to hold cells at -60mV
# sweep where there is no holding
def restingmembrane(vmfile,cell_channel):
    ch1, sweep_len, block = load_traces(vmfile, cell_channel)
    restingmem=np.median(ch1[:,0])
    return restingmem


#quality control step, deciding whether spontaneous. minis recording is going to be analyzed
#checking on which sweeps the dynamic range (max - min signal, disregarding minis) is smaller than a value
def rec_stability (filename, cell_chan, max_range):
    ch1, sweep_len, block = load_traces(filename, cell_chan) #ch1: each sweep is a column

    min_vals = []
    max_vals = []
    good_swps = []
    bad_swps = []
    #finding the max/min value for each sweep (along each column axis = 0)
    for swp in range(len(ch1[0])):
        signal_no_test_pulse = ch1[4250:,swp]

        max_val = np.amax(signal_no_test_pulse)
        loc_max = np.where(signal_no_test_pulse == max_val)
        #if more than 1 location have the same max, find the biggest average interval (200ms around max_val)
        avg_max = []

        #print(len(loc_max))

        for i in range(len(loc_max)):
            max_interval = signal_no_test_pulse[int(loc_max[0][i])-2000:int(loc_max[0][i])+2000]
            avg_max_interval = np.mean(max_interval)
            avg_max.append(avg_max_interval)

        max_val = np.amax(np.array(avg_max))
        max_vals.append(max_val)
        
        min_val = np.amin(signal_no_test_pulse)
        loc_min = np.where(signal_no_test_pulse == min_val)
        #do the same for min val
        avg_min = []
        for j in range(len(loc_min)):
            min_interval = signal_no_test_pulse[int(loc_min[0][i])-2000:int(loc_min[0][i])+2000]
            avg_min_int = np.mean(min_interval)
            avg_min.append(avg_min_int)

        min_val = np.amin(np.array(avg_min))
        min_vals.append(min_val)

        dyn_range = max_val - min_val
        #sweep nums that can be analysed
        if dyn_range < max_range:
            good_swps.append(swp+1)
        else:
            bad_swps.append(swp+1)

    return min_vals, max_vals, good_swps, bad_swps, dyn_range


# def main():
#     filename = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP220426/22427042.abf'
#     cell_chan = 3
#     inj = 'full'
#     APproperties(filename, cell_chan, inj)

# main()
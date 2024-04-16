import interpFL
import matplotlib.pyplot as plt
import neo
import numpy as np
from scipy.signal import find_peaks
import stimulation_windows_ms as stim_win
import os
import funcs_human_characterisation as hcf
import math


# con_screen_file = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP220615/22616005.abf'
# pre_cell_chan = 4
# post_cell_chan = 7


#loads presynaptic cell and quality checks the sweeps
#threshold for acceptable Rm can be set within
def presynaptic_screen(con_screen_file, pre_cell_chan):
    con_screen_data = hcf.load_traces(con_screen_file)
    chan_name = 'Ch' + str(pre_cell_chan)

    pre_sig = con_screen_data[chan_name][0]
    sweep_count = np.shape(pre_sig)[1]
    vmO = np.mean(pre_sig[:,0][0:4000]) #baseline 10 ms 

    es = []
    for i in range(sweep_count): #intrasweep control for excessively depolarised
    #cells (above -50mV) or a drift in Vm of more than 10% from start to end of sweep
        vm1 = np.mean(pre_sig[:,i][0:4000])
        vm2 = np.mean(pre_sig[:,i][197999:199999])
        max_val = np.max(pre_sig[:,i])
        if (vm1 > -50) or vm1-(vm1 * -0.1) > vm2 or vm2 >vm1 + (vm1 * -0.1):
            es.append(i)
            print("Excluding swp # " + str(i) + '; drift more than 0.1*RMP start to end or RMP > -50') 
        elif vmO - (vmO * -0.1) > vm1 or vm1 > vmO + (vmO * -0.1):
            es.append(i)
            print("Excluding swp # " + str(i) + '; drift more than 0.1* RMP') 
        #this statement accounts for APs that reverse below 0
        elif max_val < 0:
            es.append(i)
            print("Excluding swp # " + str(i) + '; APs < 0') 
    exclude = 'no'
    if len(es) == sweep_count:
        print('Stop analsis')
        # es = []
        # print('changing the es = [], to continue analysis')
        # exclude = 'yes'

    pre_sig = np.delete(pre_sig, es, axis=1)
    return pre_sig, es, vmO

#loads postsynaptic cell and quality checks sweeps. Automatically removes
#sweeps excluded from presynaptic screen (presynaptic function should
#therefore be run first)
def postsynaptic_screen(con_screen_file, post_cell_chan, es):
    con_screen_data = hcf.load_traces(con_screen_file)
    chan_name = 'Ch' + str(post_cell_chan)

    post_sig = con_screen_data[chan_name][0]
    post_sig = np.delete(post_sig, es, axis=1)
    sweep_count = np.shape(post_sig)[1]

    if post_sig[0].size == 0:
        return post_sig, []
    vmO = np.mean(post_sig[:,0][0:4000]) #from column 0, 0:4000

    es2 = []
    for i in range(0,sweep_count): #intrasweep control for excessively depolarised
    #cells (above -50mV) or a drift in Vm of more than 10% from start to end of sweep
        vm1 = np.mean(post_sig[:,i][0:4000]) #baseline 
        vm2 = np.mean(post_sig[:,i][197999:199999])
        post_amp = np.max(post_sig[:,i]) - np.min(post_sig[:,i])  
        if (vm1 > -50) or vm1 - (vm1 * -0.1) > vm2 or vm2 > vm1 + (vm1 * -0.1):
            es2.append(i)
        elif vm1 < -80:
            es2.append(i)
        elif vmO-(vmO*-0.1)>vm1 or vm1>vmO+(vmO*-0.1):
            es2.append(i)

    post_sig = np.delete(post_sig, es2, axis=1)
    return post_sig, vmO

#cuts out stimulation window for analysis
def get_analysis_window(pre_sig, post_sig):

    mean_pre = np.mean(pre_sig, axis = 1)

    preAPs = find_peaks(mean_pre, height=0, distance=800)
    j = 0 #starting with min AP heigth 0
    while len(preAPs[0]) == 0:
        j = j-5
        preAPs=find_peaks(mean_pre, height=j, distance=800)
    print("Pre APs found with heigth " + str(j))

    num_APs = len(preAPs[0])
    pre_window = mean_pre[preAPs[0][0]-750:preAPs[0][num_APs-1]+750]

    mean_post = np.mean(post_sig, axis = 1)
    post_window = mean_post[preAPs[0][0]-750:preAPs[0][num_APs-1]+750]    
    preAPs_shifted = find_peaks(pre_window, height=j, distance=800) #shifts so that stim_sindow starts at 0
    
    return mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs

#removes sweeps where  POSTsynaptic cell fires APs during stimulation window
def remove_sweeps_with_POSTsynAPs(pre_sig, post_sig, preAPs):
    remove_sweeps=[]
    num_aps = len(preAPs[0])

    for sw in range(0,np.shape(post_sig)[1]):
        if np.max(post_sig[:,sw][preAPs[0][0] - 750:preAPs[0][num_aps-1] + 750]) > 0:
            remove_sweeps.append(sw)
            remove_sweeps.reverse()

    pre_sig = np.delete(pre_sig, remove_sweeps, axis=1)
    post_sig = np.delete(post_sig, remove_sweeps, axis=1)
    
    return pre_sig, post_sig

#takes differential and uses this to find psp peaks
def find_postsynaptic_peaks(post_window, preAPs_shifted):
    post_d1 = np.diff(post_window,1)
    post_d1_peaks = find_peaks(post_d1 ,distance=800) 

    num_aps = len(preAPs_shifted[0])
    PSPs = np.ndarray([num_aps,2])
    for p in range(num_aps):
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
        x_20 = (delta * 0.2 + bl[i][0]).item()
        x_80 = (delta * 0.8 + bl[i][0]).item()
        fit_start = interpFL.findLevels(PSP_window, x_20, mode='rising')
        fit_end = interpFL.findLevels(PSP_window, x_80, mode='rising')
        if fit_start[0].size == 0 or fit_end[0].size == 0:
            onsets[i,0] = math.nan
            continue
        if fit_start[0][0] == fit_end[0][0]:
            onsets[i,0] = math.nan
            continue
        x = range(155)
        if fit_end[0][0] < fit_start[0][0]:
            onsets[i,0] = math.nan
            continue
            # end_corr = fit_end[0][fit_end[0] > fit_start[0][0]][0]
            # x1 = x[fit_start[0][0] : end_corr]
        else:
            x1 = x[fit_start[0][0] : fit_end[0][0]]
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
        if math.isnan(onsets[i].item()):
            latency[i,0] = math.nan
            continue
        latency[i,0] = int(onsets[i].item()) - preAPs_shifted[0][i]
    latency = latency/20
    return latency

def get_amps(PSPs, bl):
    amps = []
    for u in range(len(PSPs)):
        amps.append(PSPs[u][0] - bl[u][0])
    while len(amps) - 4 < 0:
                amps.append(math.nan)
    return amps


#for stim at diff freqs


def get_pre_aps_diff_freqs (con_screen_file, pre_cell_chan, post_cell_chan, hz):
    con_screen_data = hcf.load_traces(con_screen_file)
    pre_ch, post_ch = 'Ch' + str(pre_cell_chan), 'Ch' + str(post_cell_chan)

    sweep_len = np.shape(con_screen_data[pre_ch][0])[0]

    stims = stim_win.stim_window_diff_freq
    win_start, win_end = stims[hz][0], stims[hz][1] 

    mean_pre = np.mean(con_screen_data[pre_ch][0], axis = 1)[win_start:win_end]
    vm0 = np.mean(con_screen_data[pre_ch][0][:,0][0:950]) #from column 0, 0:950
    preAPs = find_peaks(mean_pre, height=0)
    j = 0
    while len(preAPs[0]) == 0:
        j = j-5
        preAPs=find_peaks(mean_pre, height=j)
    print("Pre APs found with heigth " + str(j))
    num_aps = len(preAPs[0])

    pre_window = mean_pre[preAPs[0][0]-30:preAPs[0][num_aps-1]+300]

    mean_post = np.mean(con_screen_data[post_ch][0], axis = 1)
    post_window = mean_post[preAPs[0][0]-30:preAPs[0][num_aps-1]+3000] 
    preAPs_shifted = find_peaks(pre_window, height=j)
    pA0 = np.mean(con_screen_data[post_ch][0][:,0][0:950]) #from column 0, 0:950

    return vm0, pA0, mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs

    # def find_postsynaptic_peaks(post_window, preAPs_shifted):
    #     post_d1 = np.diff(post_window,1)
    #     post_d1_peaks = find_peaks(post_d1, distance=800) 

    #     num_aps = len(preAPs_shifted[0])
    #     PSPs = np.ndarray([num_aps,2])
    #     for p in range(num_aps):
    #         PSPs[p,0] = np.max(post_window[preAPs_shifted[0][p]+20:\
    #                         preAPs_shifted[0][p]+150])
    #         PSPs[p,1] = np.argmax(post_window[preAPs_shifted[0][p]+20:\
    #                             preAPs_shifted[0][p]+150])+preAPs_shifted[0][p]+20
    #     return PSPs


#functions for pre_cahn in IC; post chan in VC

def presynaptic_screen_IC(con_screen_file, pre_cell_chan):
    con_screen_data = hcf.load_traces(con_screen_file)
    chan_name = 'Ch' + str(pre_cell_chan)

    pre_sig = con_screen_data[chan_name][0]
    sweep_count = np.shape(pre_sig)[1]
    vmO = np.mean(pre_sig[:,0][10_000:14_000]) #from column 0, 0:4000

    es = []
    for i in range(sweep_count): #intrasweep control for excessively depolarised
    #cells (above -50mV) or a drift in Vm of more than 10% from start to end of sweep
        vm1 = np.mean(pre_sig[:,i][7000:11000])
        vm2 = np.mean(pre_sig[:,i][40_000:44_000])
        max_val = np.max(pre_sig[:,i])
        if (vm1 > -50) or vm1-(vm1 * -0.1) > vm2 or vm2 >vm1 + (vm1 * -0.1):
            es.append(i)
            print("Excluding swp # " + str(i) + '; drift more than 0.1*RMP start to end or RMP > -50') 
        elif vmO - (vmO * -0.1) > vm1 or vm1 > vmO + (vmO * -0.1):
            es.append(i)
            print("Excluding swp # " + str(i) + '; drift more than 0.1* RMP') 
        #this statement accounts for APs that reverse below 0
        elif max_val < 0:
            es.append(i)
            print("Excluding swp # " + str(i) + '; APs < 0') 
    exclude = 'no'
    if len(es) == sweep_count:
        print('Stop analsis')
        # es = []
        # print('changing the es = [], to continue analysis')
        # exclude = 'yes'

    pre_sig = np.delete(pre_sig, es, axis=1)
    return pre_sig, es, vmO


def postsynaptic_screen_VC (con_screen_file, post_cell_chan, es):
    con_screen_data = hcf.load_traces(con_screen_file)
    chan_name = 'Ch' + str(post_cell_chan)

    post_sig = con_screen_data[chan_name][0]
    post_sig = np.delete(post_sig, es, axis=1)
    sweep_count = np.shape(post_sig)[1]

    if post_sig[0].size == 0:
        return post_sig, [] 
    vmO = np.mean(post_sig[:,0][30_000:34_000]) #from column 0, 0:4000

    # es2 = []
    # for i in range(0,sweep_count): #intrasweep control for excessively depolarised
    # #cells (above -50mV) or a drift in Vm of more than 10% from start to end of sweep
    #     vm1 = np.mean(post_sig[:,i][7_000:11_000])
    #     vm2 = np.mean(post_sig[:,i][40_000:44_000])
    #     if (vm1 < - 800) or vm1 - (vm1 * -0.1) > vm2 or vm2 > vm1 + (vm1 * -0.1):
    #         es2.append(i)
    #     elif vmO-(vmO*-0.1) > vm1 or vm1 > vmO+(vmO*-0.1):
    #         es2.append(i)
    # post_sig = np.delete(post_sig, es2, axis=1)
    return post_sig, vmO

def get_analysis_window_VC (pre_sig, post_sig):

    mean_pre = np.mean(pre_sig, axis = 1)

    preAPs = find_peaks(mean_pre, height=0, distance=400)
    j = 0 #starting with min AP heigth 0
    while len(preAPs[0]) == 0:
        j = j-5
        preAPs=find_peaks(mean_pre, height=j, distance=400)
    print("Pre APs found with heigth " + str(j))

    num_APs = len(preAPs[0])
    pre_window = mean_pre[preAPs[0][0]-750:preAPs[0][num_APs-1]+750]

    mean_post = np.mean(post_sig, axis = 1)
    post_window = mean_post[preAPs[0][0]-750:preAPs[0][num_APs-1]+750]    
    preAPs_shifted = find_peaks(pre_window, height=j, distance=400) #shifts so that stim_sindow starts at 0
    
    return mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs

def find_postsynaptic_peaks_VC(post_window, preAPs_shifted):
    post_window_negative = post_window * -1 #so that the find_peaks_works

    post_d1 = np.diff(post_window_negative,1)
    post_d1_peaks = find_peaks(post_d1 ,distance=800) 

    num_aps = len(preAPs_shifted[0])
    PSPs = np.ndarray([num_aps,2])
    for p in range(num_aps):
        PSPs[p,0] = np.min(post_window[preAPs_shifted[0][p]+20:\
                           preAPs_shifted[0][p]+150])
        PSPs[p,1] = np.argmin(post_window[preAPs_shifted[0][p]+20:\
                              preAPs_shifted[0][p]+150])+preAPs_shifted[0][p]+20
    return PSPs

def get_onsets_VC(preAPs_shifted, post_window, PSPs, bl):
    
    num_aps = len(preAPs_shifted[0])
    onsets = np.ndarray([num_aps,1])
    for i in range(num_aps): #for each event
        PSP_window = post_window[preAPs_shifted[0][i]-5:\
                                 preAPs_shifted[0][i]+150] #more precise adjustment; smaller window
        #calc diff bewteen bl and peak
        delta = PSPs[i][0] - bl[i] #event size
        x_20 = (delta * 0.2 + bl[i][0]).item()
        x_80 = (delta * 0.8 + bl[i][0]).item()
        fit_start = interpFL.findLevels(PSP_window, x_20, mode='falling')
        fit_end = interpFL.findLevels(PSP_window, x_80, mode='falling')
        if fit_start[0].size == 0 or fit_end[0].size == 0:
            onsets[i,0] = math.nan
            continue
        if fit_start[0][0] == fit_end[0][0]:
            onsets[i,0] = math.nan
            continue
        x = range(155)
        if fit_end[0][0] < fit_start[0][0]:
            onsets[i,0] = math.nan
            continue
            # end_corr = fit_end[0][fit_end[0] > fit_start[0][0]][0]
            # x1 = x[fit_start[0][0] : end_corr]
        else:
            x1 = x[fit_start[0][0] : fit_end[0][0]]
        fit1 = np.polyfit(x1, PSP_window[x1], 1)
        foot = (bl[i][0]-fit1[1])/fit1[0] #where the fit crosses the local baseline
        onsets[i,0] = int(foot) + preAPs_shifted[0][i]-5
        
    return onsets

def get_amps_VC (PSPs, bl):
    amps = []
    for u in range(len(bl)):
        amps.append(bl[u][0] - PSPs[u][0])
    while len(amps) - 4 < 0:
                amps.append(math.nan)
    return amps

#%%
import neo
import numpy as np
from detect_peaks import detect_peaks
import math
import pyabf
import matplotlib.pyplot as plt
import funcs_sorting as sort


#%%
# =============================================================================
# Basic functions for characterisation of human neurons
# =============================================================================
def read_abf(filename):
    '''
    returns a neo block object, where analog signals data is stored
    '''
    r = neo.io.AxonIO(filename=filename)
    block = r.read(signal_group_mode = 'split-all')[0] 
    channels = len(block.segments[0].analogsignals)

    for ch in range(0, channels):
        name = block.segments[0].analogsignals[ch].name
        if name == '_Ipatch' or name == 'IN0':
            block.segments[0].analogsignals[ch].name = 'Ch1'
    return block

def get_recording_time(filename):
    '''
    gives the times of the recording of this filename
    '''
    block = read_abf(filename)
    rec_time = block.rec_datetime
    return rec_time

def get_cell_IDs (filename_char, slic, active_channels):
    end_fn = filename_char.rfind('/') + 1
    cell_IDs = []
    for ch in active_channels:
        cellID = filename_char[end_fn:-7] + slic + 'c' + str(ch)
        cell_IDs.append(cellID)
    return cell_IDs

def get_new_cell_IDs(filename_char, slic, active_channels, patcher):
    patcher_dict = {'Verji':'vm', 'Rosie': 'rs'}
    
    end_fn = filename_char.rfind('/') + 1
    cell_IDs = []
    for ch in active_channels:
        cellID = patcher_dict[patcher] + filename_char[end_fn:-7] + slic + 'c' + str(ch)
        cell_IDs.append(cellID)
    return cell_IDs

def get_new_cell_IDs_fn(fn, slic, active_channels, patcher):
    patcher_dict = {'Verji':'vm', 'Rosie': 'rs'}
    
    cell_IDs = []
    for ch in active_channels:
        cellID = patcher_dict[patcher] + fn[:-7] + slic + 'c' + str(ch)
        cell_IDs.append(cellID)
    return cell_IDs


def get_connection_ID (filename_con_screen, slic, pre_chan, post_chan):
    end_fn = filename_con_screen.rfind('/') + 1
    con_ID = filename_con_screen[end_fn:-7] + slic + 'c' + str(pre_chan) + '#' + str(post_chan)
    return con_ID

def read_inj(inj):
    if inj == "full":
        inj=[-300,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400,450,500,550,
        600,700,800,900,1000,1100,1200,1300,1400]
    else:
        inj = inj
    return inj

def find_charact_onset_offset(char_fn):
    ''' 
    returrs offset and onset of the charact steps
    '''
    trace = pyabf.ABF(char_fn)
    trace.setSweep(sweepNumber = 0, channel = 0)

    onset = np.where(trace.sweepC != 0)[0][0] - 1
    offset = np.where(trace.sweepC != 0)[0][-1]
     
    return onset, offset

def get_inj_current_steps(fn):
    '''
    loads the characterization fn
    returns the steps currents used
    assuming the same steps were given to all channels
    '''
    char = pyabf.ABF(fn)

    if char.userList is None:
        inj = []
        for swp in range(char.sweepCount):
            char.setSweep(swp)
            if len(np.unique(char.sweepC)) == 1:
                inj.append(np.unique(char.sweepC)[0])
                continue
            inj.append(np.unique(char.sweepC[char.sweepC != 0])[0])
            #plotting_funcs.plot_trace(fn, swp, 2)
    else:
        inj = char.userList
    return inj

def load_traces(filename):
    '''
    returns a dictionary; key: channel name (e.g. 'Ch1')
    value: np.array(sweep_len, sweep_count), each column is a sweep, len column - sweep len
    '''
    block = read_abf(filename)
    sweep_count = len(block.segments) #number of sweeps
    sweep_len = len(block.segments[0].analogsignals[0])
    channels = len(block.segments[0].analogsignals)

    data_dict = {}
    for ch in range(0, channels):
        name = block.segments[0].analogsignals[ch].name
        ch_data = np.ndarray([sweep_len, sweep_count])                
        for i in range(0,sweep_count):
            ch_data[:,i] = block.segments[i].analogsignals[ch].view(np.recarray).reshape(sweep_len)

        data_dict[name]=[ch_data]  
    return data_dict
    #len(data_dict) # number of channels
    #np.shape(data_dict1['Ch1'][0])[0] - sweep len
    #np.shape(data_dict1['Ch1'][0])[1] - number of sweeps

def get_abf_info(filename, cell_chan, sweep_count, sweep_len):
    '''
    returns descriptive parameters about the block object
    '''
    block = read_abf(filename)
    middle_swp_num = int(sweep_count/2)
    if len(block.segments[middle_swp_num].analogsignals) < 8:
        all_chans_len = len(block.segments[middle_swp_num].analogsignals)
        active_channels = []
        for ch in range(all_chans_len):
            chan_name = int(block.segments[0].analogsignals[ch].name[-1])
            active_channels.append(chan_name)
        chan_index = active_channels.index(int(cell_chan))
        #have to find the index of this chan in all chans
        sampl_rate = block.segments[middle_swp_num].analogsignals[chan_index].sampling_rate
        units = block.segments[middle_swp_num].analogsignals[chan_index].units
        times = np.linspace(0,sweep_len,sweep_len)/sampl_rate
        return sampl_rate, units, times
    #signal = block.segments[middle_swp_num].analogsignals[cell_chan-1].view(np.recarray).reshape(sweep_len).tolist()
    sampl_rate = block.segments[middle_swp_num].analogsignals[int(cell_chan)-1].sampling_rate
    units = block.segments[middle_swp_num].analogsignals[int(cell_chan)-1].units
    times = np.linspace(0,sweep_len,sweep_len)/sampl_rate
    return sampl_rate, units, times

#series and input resistance from vc file
def get_access_resistance(vctpfile, channels):
    '''
    returns two lists with length channels
    ResRa - access resistance
    ResRi - input resistance 
    '''
    data_dict = load_traces(vctpfile)
    ResRa_all, ResRi_all = [], []
    for ch in channels:
        key = 'Ch' + str(ch)
        mean = np.mean(data_dict[key][0], axis = 1) 
        
        #detect_peaks(mean,mpd = 999, edge='both', valley=True) #mpd - peaks separated by min peak distance
        minpeak = np.min(mean)
        avgBL = np.mean(mean[0:1000]) #baseline average
        avgSS = np.mean(mean[np.argmax(mean)-500:np.argmax(mean)-50]) #argmax - the index of the max value
        if avgBL == avgSS: 
            ResRa = math.nan
            ResRi = math.nan
            print('trying to compute access resistance of a dead cell :(')
            ResRa_all.append(ResRa)
            ResRi_all.append(ResRi)
            continue
        RaI = avgBL - minpeak
        RiI = avgBL - avgSS
        ResRa = (0.004/(RaI*1e-12))/1000000 #0.004 - size of the step, Resistance in megaOhm
        ResRi = (0.004/(RiI*1e-12))/1000000
        
        ResRa_all.append(ResRa)
        ResRi_all.append(ResRi)
    return ResRa_all, ResRi_all

# Vm sweep before holding current is added to hold cells at -60mV
def get_RMP (vmfile, channels):
    vm_dict = load_traces(vmfile)
    resting_mems = []
    for ch in channels:
        key = 'Ch' + str(ch)
        ch1 = vm_dict[key][0]
        resting_mems.append(np.median(ch1[:,0]))
    return resting_mems

def find_charact_onset_offset(char_fn):
    ''' 
    returrs offset and onset of the charact steps
    '''
    trace = pyabf.ABF(char_fn)
    trace.setSweep(sweepNumber = 0, channel = 0)

    onset = np.where(trace.sweepC != 0)[0][0] - 1
    offset = np.where(trace.sweepC != 0)[0][-1]
     
    return onset, offset

def get_hyperpolar_param(charact_data, channels, inj, onset = 2624, offset = 22624,mc = np.ndarray([5,3])):
    '''
    returns 5 lists with length channels
    '''
    params = tau_all, capacitance_all, mc_all, V65_all, RMPs_char = [], [], [], [], []
    #RMPs_char = []
    for ch in channels:
        key = 'Ch' + str(ch)
        ch1 = charact_data[key][0]
        resting_char = np.median(ch1[:,5]) #when the inj = 0mV
        V65s = []
        mc = np.ndarray([21,3])
        for i in range(0,5):
            I = inj[i]*1e-12 #check step size
            bl = np.median(ch1[0:onset-20,i])
            ss = np.median(ch1[offset-2000:offset-1,i]) #steady state, during the current step
            swp = list(ch1[:,i])
            Vdiff = bl-ss #voltage deflection size
            v65 = Vdiff * 0.63
            V65 = bl - v65
            V65s.append(V65)
            if list(filter(lambda ii: ii < V65, swp)) == []:
                continue
            else:
                res = list(filter(lambda ii: ii < V65, swp))[0] #takes the first value in swp < V65
                tau65 = swp.index(res) #index of res
                R = (Vdiff/1000)/-I     
                tc = tau65 - onset
                mc[i,0] = tc * 0.05 #membranec capacitance; tc - time constant
                mc[i,1] = R * 1e-6 #resistance
                mc[i,2] = tc * 5e-5 / R #capacitance
        mc[:,2] = mc[:,2]/1e-12  
        tau = mc[1,0]
        capacitance = mc[1,2]
        ch_params = [tau, capacitance, mc, V65s, resting_char]

        for i, param in enumerate(params):
            param.append(ch_params[i])
    return tau_all, capacitance_all, mc_all, V65_all, RMPs_char

def get_max_spikes(charact_data, channels):
    '''
    returns a list of max number of spikes for each channel
    '''
    max_spikes_all = []
    for ch in channels:
        key = 'Ch' + str(ch)
        ch1 = charact_data[key][0]
        max_spikes = 0
        for j in range(len(ch1[0])): #loop through all swps
            pks = detect_peaks(ch1[:,j], mph=2, mpd=5) # detects the peaks for each timepoint? 
            if len(pks)> max_spikes:
                max_spikes = len(pks) #find the max number of spikes (peaks)
        max_spikes_all.append(max_spikes)
    return max_spikes_all 

def get_ap_param(charact_data, channels, inj, max_spikes):
    '''
    returns a nd.array with inj, peak loc, sweep num
    spike_counts for each inj
    array with all peak locations [inj, peak, sweep num]
    '''
    params = [Rheobase_all, AP_all, THloc_all, TH_all, APheight_all, max_depol_all, max_repol_all] = [[], [], [], [], [], [], []]
    for a, ch in enumerate(channels):
        key = 'Ch' + str(ch)
        ch1 = charact_data[key][0]
        
        if max_spikes[a] == 0:
            AP = []
            Rheobase = math.nan
            TH, THloc = math.nan, math.nan
            max_depol = math.nan
            max_repol = math.nan
            APheight = math.nan
            ch_params = [Rheobase, AP, THloc, TH, APheight, max_depol, max_repol]

            for i, param in enumerate(params):
                param.append(ch_params[i])
            continue
        peaks = np.empty([len(inj), max_spikes[a], 3])
        peaks.fill(np.nan)
        for i in range(len(ch1[0])): #for all swps
            pks = detect_peaks(ch1[:,i], mph=20,mpd=50) 
            peaks[i,0:len(pks), 0] = inj[i] #injected current
            peaks[i,0:len(pks), 1] = pks #
            peaks[i,0:len(pks), 2] = ch1[pks,i] #sweep number
        # number of spikes in each step
        spike_counts = np.ndarray([len(inj),2])
        for i, j in enumerate(inj):
            spike_counts[i,0] = j
            spike_counts[i,1] = (np.sum(np.isfinite(peaks[i,:,:])))/3
        
        spikes = np.where(np.isfinite(peaks)) 
        if spikes[0].size == 0:
            AP = []
            Rheobase = math.nan
            TH, THloc = math.nan, math.nan
            max_depol = math.nan
            max_repol = math.nan
            APheight = math.nan
            ch_params = [Rheobase, AP, THloc, TH, APheight, max_depol, max_repol]

            for i, param in enumerate(params):
                param.append(ch_params[i])
            continue

        first_spike = spikes[0][0]
        Rheobase = inj[first_spike]  
        
        if np.max(spike_counts[:,1]) == 1:
            print('MAX number of AP = 1 for ' + key)
            first_spiking_sweep = np.where(spike_counts[:,1]==1)[0][0]
        else:
            first_spiking_sweep = np.where(spike_counts[:,1]>1)[0][0] #where there is more than 1 AP
        
        if np.max(spike_counts[:,1]) == 1:
            ap = 0
        else:
            ap = 1

        peak_loc = np.where(ch1[:,first_spiking_sweep] == peaks[first_spiking_sweep, ap ,2])[0][0]
        AP = ch1[:, first_spiking_sweep][peak_loc - 200:peak_loc+200]
        d1_AP = np.diff(AP)*20
        THlocs = np.where(d1_AP[:195] > 10)
        if THlocs[0].size != 0:
            THloc = THlocs[0][0]
            TH = AP[THloc-1]
            APheight = peaks[first_spiking_sweep, ap ,2]-TH #amplitude of AP, the 2nd AP
            max_depol = np.max(d1_AP[:200]) #how quickly is the voltage change occuring
            max_repol = np.min(d1_AP[-200:])
            ch_params = [Rheobase, AP, THloc, TH, APheight, max_depol, max_repol]

            for i, param in enumerate(params):
                param.append(ch_params[i])
            continue
        else:
            print("Only 1 SLOW AP found for " + key)
            TH, THloc = math.nan, math.nan
            max_depol = math.nan
            max_repol = math.nan
            APheight = math.nan
            ch_params = [Rheobase, AP, THloc, TH, APheight, max_depol, max_repol]
            for i, param in enumerate(params):
                param.append(ch_params[i])

    return Rheobase_all, AP_all, THloc_all, TH_all, APheight_all, max_depol_all, max_repol_all 

def get_AP_HW(channels, AP_all, APheight_all, TH_all):
    '''
    returns AP halfwidt in ms at 1/2 of the AP amplitude,
    measured from AP threshold to AP peak
    '''
    AP_HWs = []
    for i, ch in enumerate(channels):
        if len(AP_all[i]) == 0:
            AP_HWs.append(math.nan)
            continue

        half_amp = TH_all[i] + APheight_all[i]/2

        t30 = TH_all[i] + APheight_all[i]*0.3
        t70 = TH_all[i] + APheight_all[i]*0.7

        rise_start = np.where(AP_all[i] > t30)[0][0]
        rise_end = np.where(AP_all[i] > t70)[0][0]

        x = range(rise_start, rise_end)
        y = AP_all[i][rise_start:rise_end]
         
        model_rise = np.polyfit(x, y, 1)
        hw1 = (half_amp - model_rise[1])/model_rise[0] 

        fall_start = np.where(AP_all[i] > t30)[0][-1]
        fall_end = np.where(AP_all[i] > t70)[0][-1]

        x_fall = range(fall_end, fall_start)
        y_fall = AP_all[i][fall_end:fall_start]

        model_fall = np.polyfit(x_fall, y_fall, 1)
        hw2 = (half_amp - model_fall[1])/model_fall[0]

        AP_HWs.append((hw2-hw1)/20) #conversion to ms

        # predict_rise = np.poly1d(model_rise)
        # rise_prediction = predict_rise(np.linspace(rise_start, rise_end, 30))
        # difference_array = np.absolute(np.linspace(rise_start, rise_end, 30)-hw1)
        # index_rise = difference_array.argmin()
        # hw1_y = rise_prediction[index_rise]

        # predict_fall = np.poly1d(model_fall)
        # fall_prediction = predict_fall(np.linspace(fall_end, fall_start, 30))
        # difference_array = np.absolute(np.linspace(fall_end, fall_start, 30)-hw2)
        # index_fall = difference_array.argmin()
        # hw2_y = fall_prediction[index_fall]

    return AP_HWs #rise_start, rise_end, rise_prediction, fall_end, fall_start, fall_prediction, hw1, hw2, hw1_y, hw2_y

def all_chracterization_params(filename, channels, inj):
    charact_dict = load_traces(filename)
    onset, offset = find_charact_onset_offset(filename)

    inj = get_inj_current_steps(filename)
    inj2 = read_inj(inj)
    if inj != inj2:
        print('fail to read inj from fn')
        inj = inj2

    tau_all, capacitance_all, mc_all, V65_all, RMPs_char = get_hyperpolar_param(charact_dict, channels, inj, onset, offset)
    max_spikes_all = get_max_spikes(charact_dict, channels)
    Rheobase_all, AP_all, THloc_all, TH_all, APheight_all, max_depol_all, max_repol_all = get_ap_param(charact_dict, channels, inj, max_spikes_all)
    AP_HWs = get_AP_HW(channels, AP_all, APheight_all, TH_all)

    keys = ['max_spikes', 'Rheobase', 'AP_heigth', 'TH', 'max_depol', 'max_repol', 'membra_time_constant_tau', 'capacitance', 'AP_halfwidth']
    vals = [max_spikes_all, Rheobase_all, APheight_all, TH_all, max_depol_all, max_repol_all, tau_all, capacitance_all, AP_HWs]
    params_dict = {}
    for i, param in enumerate(vals):
        key = keys[i]
        params_dict[key] = param
    
    return params_dict

def get_ap_param_for_plotting (charact_data, channels, inj, max_spikes):
    '''
    '''
    params = [first_spike_all, peaks_all, spike_counts_all, first_spiking_sweep_all] = [[], [], [], []]
    for a, ch in enumerate(channels):
        key = 'Ch' + str(ch)
        ch1 = charact_data[key][0]
        
        if max_spikes[a] == 0:
            peaks_all.append([])
            first_spike_all.append(float('nan'))
            spike_counts_all.append(0)
            first_spiking_sweep_all.append([])
            # peaks = []
            # first_spike = float('nan')
            # ch_params = [first_spike, peaks]
            # for i, param in enumerate(params):
            #     param.append(ch_params[i])
            continue
        peaks = np.empty([len(inj), max_spikes[a], 3])
        peaks.fill(np.nan)
        for i in range(len(ch1[0])): #for all swps
            pks = detect_peaks(ch1[:,i], mph=2,mpd=5) 
            peaks[i,0:len(pks), 0] = inj[i] #injected current
            peaks[i,0:len(pks), 1] = pks #
            peaks[i,0:len(pks), 2] = ch1[pks,i] #sweep number
        # number of spikes in each step
        spike_counts = np.ndarray([len(inj),2])
        for i, j in enumerate(inj):
            spike_counts[i,0] = j
            spike_counts[i,1] = (np.sum(np.isfinite(peaks[i,:,:])))/3
        
        spikes = np.where(np.isfinite(peaks)) 
        first_spike = spikes[0][0]

        if np.max(spike_counts[:,1]) == 1:
            print('MAX number of AP = 1 for ' + key)
            first_spiking_sweep = np.where(spike_counts[:,1]==1)[0][0]
        else:
            first_spiking_sweep=np.where(spike_counts[:,1]>1)[0][0] 

        ch_params = [first_spike, peaks, spike_counts, first_spiking_sweep]
        for i, param in enumerate(params):
            param.append(ch_params[i])
    return first_spike_all, peaks_all, spike_counts_all, first_spiking_sweep_all


def get_rheobase_from_ramp(fn, chans):
    ''' 
    finds where dv/dt > 20 mv/ms
    finds the first 
    '''
    ramp_dict = load_traces(fn)
    num_swps = np.shape(ramp_dict['Ch' + str(chans[0])][0])[1]
    swp = 0
    rheos, THs, THs_in_trace, swps = [], [], [], []
    for ch in chans:
        key = 'Ch' + str(ch)
        ramp = ramp_dict[key][0][:, swp]        
        AP1_peak = detect_peaks(ramp, mph = 20, mpd = 30)
        if len(AP1_peak) == 0:
            #runs through the rest of the swps to find peaks
            swp, AP1_peak = try_other_swps(num_swps, ramp_dict, key)
        if swp == num_swps-1 and len(AP1_peak) == 0:
            print('not possible to find rheobase, no APs fired')
            rheos.append(math.nan)
            THs.append(math.nan)
            THs_in_trace.append(math.nan)
            swps.append(math.nan)
            continue
            
        search_frame = np.arange(AP1_peak[0] - 300,AP1_peak[0])
        AP1 = ramp[search_frame]
        d1_AP1 = np.diff(AP1)*20
        THloc = np.where(d1_AP1[:295] > 20)
        if len(THloc[0]) < 1:
            print('trying 10 mV/ms. APs seems slow')
            THloc = np.where(d1_AP1[:295] > 10)
            if len(THloc[0]) < 1 :
                print('No Rheobase found')
                rheos.append(math.nan)
                THs.append(math.nan)
                THs_in_trace.append(math.nan)
                swps.append(math.nan)
                continue

        TH = AP1[THloc[0][0]-1]
        THloc_in_trace = np.where(ramp == TH)[0][0]

        ramp_abf = pyabf.ABF(fn)
        ramp_abf.setSweep(sweepNumber = swp, channel = 0) #because the stimulus is the same for all channels
        rheobase = ramp_abf.sweepC[THloc_in_trace-1]

        rheos.append(rheobase)
        THs.append(TH)
        THs_in_trace.append(THloc_in_trace)
        swps.append(swp)
    return rheos, THs, THs_in_trace, swps

def try_other_swps(num_swps, ramp_dict, key):
    '''
    to be used when no firing on the first sweep
    or the first AP is slow/ not picked up
    '''
    for swp in range(1,num_swps):
        ramp = ramp_dict[key][0][:, swp] 
        AP1_peak = detect_peaks(ramp, mph = 20, mpd = 30)
        if len(AP1_peak) > 0:
            return swp, AP1_peak
    return swp, AP1_peak


#quality control step, deciding whether spontaneous or mini recording is going to be analyzed
#checking on which sweeps the dynamic range (max - min signal, disregarding minis) is smaller than a value
def rec_stability (mini_filename, active_chans, max_drift, min_holding = -800):
    mini_dict = load_traces(mini_filename)

    rec_stability = {}
    for ch in active_chans:
        key = 'Ch' + str(ch)
        ch1 = mini_dict[key][0]

        all_swps = np.arange(0, len(ch1[0]), 1).tolist()
        vm_1st_swp = np.mean(ch1[4250:,0][0:2000])
        vm_last_swp = np.mean(ch1[4250:, len(ch1[0])-1][0:2000])

        min_vals, max_vals, bad_swps, dyn_range  = [], [], [], []
 
        for swp in range(len(ch1[0])):
            signal_no_test_pulse = ch1[4250:,swp]
            vm_start = np.mean(signal_no_test_pulse[0:2000])
            if vm_start < min_holding :
                bad_swps.append(swp)
                continue

        #if more than 1 location have the same max, find the biggest average interval (200ms around max_val)
            max_val = np.amax(signal_no_test_pulse)   
            loc_max = np.where(signal_no_test_pulse == max_val)[0][-1] #take the last loc of max
            if loc_max <= 1000:
                avg_max_interval = np.mean(signal_no_test_pulse[loc_max: loc_max + 2000])
            else:
                avg_max_interval = np.mean(signal_no_test_pulse[loc_max - 1000 : loc_max + 1000])
            max_vals.append(avg_max_interval)
            
            min_val = np.amin(signal_no_test_pulse)
            loc_min = np.where(signal_no_test_pulse == min_val)[0][-1]
            if loc_min <= 1000:
                avg_min_int = np.mean(signal_no_test_pulse[loc_min:loc_min + 2000])
            else:
                avg_min_int = np.mean(signal_no_test_pulse[loc_min - 1000 : loc_min + 1000])
            min_vals.append(avg_min_int)

            dyn_range = abs(avg_max_interval - avg_min_int)
            if dyn_range > max_drift :
                bad_swps.append(swp)
            # else:
            #    good_swps.append(swp)
        good_swps = [x for x in all_swps if x not in bad_swps]
            
        rec_stability[str(ch)] = {'holding_sweep1': vm_1st_swp, 'holding_last_sweep' : vm_last_swp, 
        'swps_to_analyse':good_swps, 'swps_to_discard':bad_swps, 'min_vals_each_swp':min_vals, 
        'max_vals_each_swp':max_vals, 'drift':dyn_range}
    return rec_stability


def get_RMP_over_time(filename, channels):
    data_dict = load_traces(filename)

    RMPs = []
    for ch in channels:
        key = 'Ch' + str(ch)
        ch1 = data_dict[key][0]
        num_swps = np.shape(ch1)[1]
        
        RMP_all = []
        if num_swps > 1:
            for swp in range(0, num_swps, 10):
                RMP = np.mean(ch1[:,swp][0:500])
                RMP_all.append(RMP)
        else: 
            data_points = int(np.shape(ch1)[0] / 100_000)
            for j in range(0,data_points):
                start = 99_500 + j*100_000
                RMP = np.mean(ch1[start : start + 500])
                RMP_all.append(RMP)
        RMPs.append(RMP_all)
    return RMPs


def get_initial_firing_rate(fn, channels, inj = 'full'):
    '''
    function creats a matrix
    'channels' --> rows in matrices, 'inj' --> columns [-300, -300, -200, -200, etc.]
    value1 = IFF, value2 = num_aps
    each injection appears twice 
    '''

    charact_data = load_traces(fn)
    inj = read_inj(inj)
    ch_0 = 'Ch' + str(channels[0])
    sampl_rate, units, times = get_abf_info(fn, channels[0], np.shape(charact_data[ch_0][0])[1], np.shape(charact_data[ch_0][0])[0])

    params = np.zeros((len(channels),len(read_inj(inj))*2))
    for i, ch in enumerate(channels):
        key = 'Ch' + str(ch)
        ch1 = charact_data[key][0]
        for j in range(0, len(ch1[0])): #loop through all swps
            k = 2*j
            pks = detect_peaks(ch1[:,j], mph=20, mpd=50) # detects peaks
            if len(pks) > 1 :
                params[i,k] = (1/(times[pks[1]] - times[pks[0]])) #inital firing freq = 1/(time 2nd AP - time 1st AP)
                params[i,k+1] = (len(pks)) #num_aps
            else :
                params[i,k] = 0 #inital firing freq = 1/(time 2nd AP - time 1st AP)
                params[i,k+1] = len(pks) #num_aps

    return params
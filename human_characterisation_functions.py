
#%%
import neo
import numpy as np
from detect_peaks import detect_peaks
import math
import sorting_functions as sort
import intrinsic_props_plotting_funcs as in_props_plot

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
    return block

def get_cell_IDs (filename_char, slic, active_channels):
    end_fn = filename_char.rfind('/') + 1
    cell_IDs = []
    for ch in active_channels:
        cellID = filename_char[end_fn:-7] + slic + 'c' + str(ch)
        cell_IDs.append(cellID)
    return cell_IDs

def read_inj(inj):
    if inj == "full":
        inj=[-300,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400,450,500,550,
        600,700,800,900,1000,1100,1200,1300,1400]
    else:
        inj = inj
    return inj

def load_traces (filename):
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
        if name == '_Ipatch' or name == 'IN0':
            name = 'Ch1'
            block.segments[0].analogsignals[ch].name = 'Ch1'
        
        ch_data = np.ndarray([sweep_len, sweep_count])                
        for i in range(0,sweep_count):
            ch_data[:,i] = block.segments[i].analogsignals[ch].view(np.recarray).reshape(sweep_len)

        data_dict[name]=[ch_data]  
    return data_dict
    #len(data_dict) # number of channels
    #np.shape(data_dict1['Ch1'][0])[0] - sweep len
    #np.shape(data_dict1['Ch1'][0])[1] - number of sweeps

def get_abf_info (filename, cell_chan, sweep_count, sweep_len):
    '''
    returns descriptive parameters about the block object
    '''
    block = read_abf(filename)
    middle_swp_num = int(sweep_count/2)
    #signal = block.segments[middle_swp_num].analogsignals[cell_chan-1].view(np.recarray).reshape(sweep_len).tolist()
    sampl_rate = block.segments[middle_swp_num].analogsignals[cell_chan-1].sampling_rate
    units = block.segments[middle_swp_num].analogsignals[cell_chan-1].units
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
        
        detect_peaks(mean,mpd = 999, edge='both', valley=True) #mpd - peaks separated by min peak distance
        minpeak = np.min(mean)
        avgBL = np.mean(mean[0:1000]) #baseline average
        avgSS = np.mean(mean[np.argmax(mean)-500:np.argmax(mean)-50]) #argmax - the index of the max value
        RaI = avgBL-minpeak
        RiI = avgBL-avgSS
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

def get_hyperpolar_param(charact_data, channels, inj, onset = 2624, offset = 22624, mc = np.ndarray([5,3])):
    '''
    returns 4 lists with length channels
    '''
    params = tau_all, capacitance_all, mc_all, V65_all = [], [], [], []
    for ch in channels:
        key = 'Ch' + str(ch)
        ch1 = charact_data[key][0]
        V65s = []
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
        ch_params = [tau, capacitance, mc, V65s]

        for i, param in enumerate(params):
            param.append(ch_params[i])
    return tau_all, capacitance_all, mc_all, V65_all

def get_max_spikes(charact_data, channels):
    '''
    returns a list of max number of spikes for each channel
    '''
    max_spikes_all = []
    for ch in channels:
        key = 'Ch' + str(ch)
        ch1 = charact_data[key][0]
        max_spikes = 0
        for i, j in enumerate(range(0, len(ch1[0]))): #loop through all swps
            pks = detect_peaks(ch1[:,j], mph=0,mpd=50) # detects the peaks for each timepoint? 
            if len(pks)> max_spikes:
                max_spikes = len(pks) #find the max number of spikes (peaks)
        max_spikes_all.append(max_spikes)
    return max_spikes_all 

def get_ap_param (charact_data, channels, inj, max_spikes):
    '''
    returns a nd.array with inj, peak loc, sweep num
    spike_counts for each inj
    array with all peak locations [inj, peak, sweep num]
    '''
    params = [Rheobase_all, AP_all, THloc_all, TH_all, APheight_all, max_depol_all, max_repol_all] = [[], [], [], [], [], [], []]
    for i, ch in enumerate(channels):
        key = 'Ch' + str(ch)
        ch1 = charact_data[key][0]
        
        if max_spikes[i] == 0:
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
        peaks = np.empty([len(inj), max_spikes[i], 3])
        peaks.fill(np.nan)
        for i, j in enumerate(range(0, len(ch1[0]))): #for all swps
            pks = detect_peaks(ch1[:,j], mph=0,mpd=50) 
            peaks[i,0:len(pks), 0] = inj[i] #injected current
            peaks[i,0:len(pks), 1] = pks #
            peaks[i,0:len(pks), 2] = ch1[pks,j] #sweep number
        # number of spikes in each step
        spike_counts = np.ndarray([len(inj),2])
        for i, j in enumerate(inj):
            spike_counts[i,0] = j
            spike_counts[i,1] = (np.sum(np.isfinite(peaks[i,:,:])))/3
        
        spikes = np.where(np.isfinite(peaks)) 
        first_spike = spikes[0][0]
        Rheobase = inj[first_spike]  
        
        if np.max(spike_counts[:,1]) == 1:
            print('MAX number of AP = 1 for ' + filename[end_fn:-4] + key)
            first_spiking_sweep = np.where(spike_counts[:,1]==1)[0][0]
        else:
            first_spiking_sweep=np.where(spike_counts[:,1]>1)[0][0] #where there is more than 1 AP
        
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
            print("Only 1 SLOW AP found for " + filename[end_fn:-4] + key)
            TH, THloc = math.nan, math.nan
            max_depol = math.nan
            max_repol = math.nan
            APheight = math.nan
            ch_params = [Rheobase, AP, THloc, TH, APheight, max_depol, max_repol]
            for i, param in enumerate(params):
                param.append(ch_params[i])

    return Rheobase_all, AP_all, THloc_all, TH_all, APheight_all, max_depol_all, max_repol_all 


def all_chracterization_params (filename, channels, inj, onset = 2624, offset = 22624):
    charact_dict = load_traces(filename)
    inj = read_inj(inj)

    tau_all, capacitance_all, mc_all, V65_all = get_hyperpolar_param(charact_dict, channels, inj)
    max_spikes_all = get_max_spikes(charact_dict, channels)
    Rheobase_all, AP_all, THloc_all, TH_all, APheight_all, max_depol_all, max_repol_all = get_ap_param (charact_dict, channels, inj, max_spikes_all)

    keys = ['max_spikes', 'Rheobase', 'AP_heigth', 'TH', 'max_depol', 'max_repol', 'membra_time_constant_tau', 'capacitance']
    vals = [max_spikes_all, Rheobase_all, APheight_all, TH_all, max_depol_all, max_repol_all, tau_all, capacitance_all]
    params_dict = {}
    for i, param in enumerate(vals):
        key = keys[i]
        params_dict[key] = param
    
    return params_dict

def get_ap_param_for_plotting (charact_data, channels, inj, max_spikes):
    '''
    '''
    params = [first_spike_all, peaks_all, spike_counts_all, first_spiking_sweep_all] = [[], [], [], []]
    for i, ch in enumerate(channels):
        key = 'Ch' + str(ch)
        ch1 = charact_data[key][0]
        
        if max_spikes[i] == 0:
            peaks = []
            first_spike = float('nan')
            ch_params = [first_spike, peaks]
            for i, param in enumerate(params):
                param.append(ch_params[i])
            continue
        peaks = np.empty([len(inj), max_spikes[i], 3])
        peaks.fill(np.nan)
        for i, j in enumerate(range(0, len(ch1[0]))): #for all swps
            pks = detect_peaks(ch1[:,j], mph=0,mpd=50) 
            peaks[i,0:len(pks), 0] = inj[i] #injected current
            peaks[i,0:len(pks), 1] = pks #
            peaks[i,0:len(pks), 2] = ch1[pks,j] #sweep number
        # number of spikes in each step
        spike_counts = np.ndarray([len(inj),2])
        for i, j in enumerate(inj):
            spike_counts[i,0] = j
            spike_counts[i,1] = (np.sum(np.isfinite(peaks[i,:,:])))/3
        
        spikes = np.where(np.isfinite(peaks)) 
        first_spike = spikes[0][0]

        if np.max(spike_counts[:,1]) == 1:
            print('MAX number of AP = 1 for ' + filename[end_fn:-4] + key)
            first_spiking_sweep = np.where(spike_counts[:,1]==1)[0][0]
        else:
            first_spiking_sweep=np.where(spike_counts[:,1]>1)[0][0] 

        ch_params = [first_spike, peaks, spike_counts, first_spiking_sweep]
        for i, param in enumerate(params):
            param.append(ch_params[i])
    return first_spike_all, peaks_all, spike_counts_all, first_spiking_sweep_all


#quality control step, deciding whether spontaneous or mini recording is going to be analyzed
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


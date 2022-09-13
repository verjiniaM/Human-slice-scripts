import human_characterisation_functions as hchf
import numpy as np
import pandas as pd
import os


#quality control step, deciding whether spontaneous or mini recording is going to be analyzed
#checking on which sweeps the dynamic range, aka drift, (max - min signal, disregarding minis) is smaller than a value - max_range
def rec_stability (filename, cell_chan, max_range):
    ch1, sweep_len, block = hchf.load_traces(filename, cell_chan) #ch1: each sweep is a column

    min_vals = []
    max_vals = []
    good_swps = []
    bad_swps = []
    drifts = []
    #finding the max/min value for each sweep (along each column axis = 0)
    for swp in range(len(ch1[0])):
        signal_no_test_pulse = ch1[4250:,swp]

        max_val = np.amax(signal_no_test_pulse)
        loc_max = np.where(signal_no_test_pulse == max_val)
        if len(loc_max) > 1:
        #if more than 1 location have the same max, find the biggest average interval (200ms around max_val)
            avg_max = []
            for i in range(len(loc_max)):
                max_interval = signal_no_test_pulse[int(loc_max[0][i])-2000:int(loc_max[0][i])+2000]
                avg_max_interval = np.mean(max_interval)
                avg_max.append(avg_max_interval)
            max_val = np.amax(np.array(avg_max))
        max_vals.append(max_val)
            
        min_val = np.amin(signal_no_test_pulse)
        loc_min = np.where(signal_no_test_pulse == min_val)
        #do the same for min val
        if len(loc_min) > 1:
            avg_min = []
            for j in range(len(loc_min)):
                min_interval = signal_no_test_pulse[int(loc_min[0][i])-2000:int(loc_min[0][i])+2000]
                avg_min_int = np.mean(min_interval)
                avg_min.append(avg_min_int)

            min_val = np.amin(np.array(avg_min))
        min_vals.append(min_val)

        dyn_range = max_val - min_val
        drifts.append(dyn_range)
        #sweep nums that can be analysed
        if dyn_range < max_range:
            good_swps.append(swp+1)
        else:
            bad_swps.append(swp+1)

    return min_vals, max_vals, good_swps, bad_swps, drifts

#calculate the average drifts
#decide upon a value for the drift
#


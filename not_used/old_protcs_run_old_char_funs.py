#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 10:18:49 2018

@author: rosie
"""

import old_protcs_human_characterisation_functions_old_protocols as ncf
import time
import numpy as np
import pandas as pd
import xlrd
import matplotlib.pyplot as plt
#%% get Cell ID and filenames|

CellID = input('Enter the cell ID:')

book = xlrd.open_workbook('/Volumes/TOSHIBA EXT/PaS Electrophysiology/Recording_indexing2 2.xlsx')
#book = xlrd.open_workbook('/Volumes/AG-Schmitz-2/Rosie/PaS Connectivity/Recording_indexing2 2.xlsx')

sheet = book.sheet_by_name('Sheet1')

#create list of cell names to search for CellID and get row of cellID (=entry)
cell_names = []
for r in range(sheet.nrows-1):
    cell_names.append(str(sheet.cell_value(r+1,0))+"_"+ str(int(sheet.cell_value(r+1,1))))
entry = cell_names.index(CellID) + 1

# check how many dep and hyp traces were collected, as default load first - can add a checking/changing option later

hyp = sheet.cell(entry,5).value
hyps = hyp.split(',')
hyp_file = str(sheet.cell_value(entry,4))+str(hyps[0])+'.abf'

dep = sheet.cell(entry,6).value
deps = dep.split(',')
dep_file = str(sheet.cell_value(entry, 4))+str(deps[0])+'.abf'

vctp = sheet.cell(entry,9).value
vctps = vctp.split(',')
vctpfile = str(sheet.cell_value(entry,4)+str(vctps[0])+'.abf')
channel = int(sheet.cell(entry, 11).value)
#%%
start = time.time()

ch1, hyp_sweeps, dep_sweeps, sweep_len, bl, message1, message2 = ncf.load_traces(hyp_file, dep_file, channel)

Ra, Ri = ncf.access_resistance(vctpfile, channel, vctp)

baseline, onset = ncf.get_baseline(ch1, bl, dep_sweeps, hyp_sweeps)

Vm = ncf.vm(ch1, bl, dep_sweeps, hyp_sweeps)

sag_ratio, peak_deflection = ncf.sag(ch1, hyp_sweeps, bl, dep_sweeps)

inj, spike_counts, max_spikes, peaks, first_spike, IO_slope, ch1, dep_sweeps = ncf.apcount(ch1, dep_sweeps, hyp_sweeps, dep_file, channel, deps)    

ax_spec = ncf.APpanelplot(first_spike, inj)

Rheobase = ncf.rheobase(inj, first_spike)

ivplot, subTHiv_fit, IR = ncf.IV_curve(hyp_sweeps, first_spike, peak_deflection, ch1, Vm)

TH, th_fig, d1 = ncf.threshold(ch1, hyp_sweeps, inj, spike_counts, peaks, dep_sweeps, first_spike, max_spikes)

Firing_freq = ncf.rel_firing_freq(first_spike, spike_counts)

max_firing_freq = ncf.max_firing_freq(spike_counts)

AP1_max_deriv, AP1_min_deriv, slope_ratio = ncf.slopes(TH, first_spike, d1, peaks)

mAHP = ncf.mahp(dep_sweeps, ch1, hyp_sweeps, baseline)

Latency = ncf.latency(TH, inj, onset)

AP_height = ncf.ap_height(peaks, first_spike, TH)

AP12int, AP910int, adaptation, ISI = ncf.intervals(TH, spike_counts, max_spikes)

cut_AP = ncf.cut_ap(TH, ch1, hyp_sweeps, spike_counts, first_spike,peaks,sweep=None, spike=1)

HW = ncf.halfwidth(first_spike, TH, ch1, hyp_sweeps, spike_counts,peaks,sweep=None, spike=1)

AHP = ncf.ahp(TH, ch1, hyp_sweeps, spike_counts, first_spike,peaks,sweep=None, spike=1)


print("---%s seconds ---" % (time.time() - start))


#%% add new entries and resave dict

if type(AP12int) == np.ndarray:
    AP12 = AP12int[0][0].item()
elif np.isnan(AP12int):
    AP12 = np.nan
else:
    AP12 = AP12int


if type(AP910int) == np.ndarray:
    AP910 = AP910int[0][0].item()
elif np.isnan(AP910int):
    AP910 = np.nan
else:
    AP910 = AP910int


if type(adaptation) == np.ndarray:
    adaptn = adaptation[0][0].item()
elif np.isnan(adaptation):
    adaptn = np.nan
else:
    adaptn = adaptation

if type(ISI) == np.ndarray:
    if np.nanmin(ISI[:,:,1]) <= 3:
        ISI_res = np.nanmin(np.ma.masked_less(ISI,3)[:,:,1])
    else: 
        ISI_res = np.nanmin(ISI[:,:,1])
elif np.isnan(ISI):
    ISI_res = np.nan
else:
    ISI_res = ISI

notes = ''

params = [CellID, Vm, sag_ratio[-1], IR, mAHP[-1][2], Latency[np.where(~np.isnan(Latency[:,1]))[0][0],1],
          AP12, AP910, adaptn, AP_height, AP1_max_deriv, AP1_min_deriv, slope_ratio, TH[first_spike, 0,3],
          AHP, HW, Rheobase, Firing_freq, ISI_res, subTHiv_fit[0], notes]

count = np.max(PaS_Cells.keys())
new_entry = count+1
PaS_Cells[new_entry] = {'CellID': params[0], 'Vm':params[1], 'Sag Ratio':params[2], 'Input Resistance': IR, 'mAHP':params[4], 'Latency': params[5], 'AP1-AP2int':params[6], 'AP9-AP10int': params[7], 'Adaptation': params[8], 'AP height': params[9], 'Max dV/dt': params[10], 'Min dV/dt': params[11], 'dV/dt ratio': params[12], 'TH': params[13], 'AHP': params[14], 'Halfwidth': params[15], 'Rheobase': params[16], 'Firing freq': params[17], 'minISI': params[18], 'IV slope': params[19], 'Notes': notes}

np.save('/Volumes/TOSHIBA EXT/Ephys_analysis/RS_PaS_Cells.npy', PaS_Cells)

plt.close('all')

#%% Load existing dict

PaS_Cells = np.load('/Volumes/TOSHIBA EXT/Ephys_analysis/RS_PaS_Cells.npy').item()
#%%

#parameters = ['CellID', 'Vm', 'Sag ratio', 'mAHP', 'Latency', 'AP1-AP2', 'AP9-AP10', 
#              'Adaptation', 'AP height', 'Max_deriv', 'Min_deriv', 'dV/dt ratio',
#              'Threshold', 'AHP', 'Halfwidth', 'Rheobase', 'firing freq', 'ISI', 'IV Slope', 'IR','theta power', 'theta freq', 'Cluster', 'staining']

PaS_DFData = pd.DataFrame.from_dict(PaS_Cells)

#%%
#run fast for verji data check

hyp_file = work_dir + filenames[indices_dict['hyp'][i]]
dep_file = work_dir + filenames[indices_dict['dep'][i]]
dep_file2 = work_dir + filenames[indices_dict['dep2'][i]]
channel = 2     #[2, 3, 4, 5, 6]

data_dict_hyp = hcf.load_traces(hyp_file)
#data_dict_dep = hcf.load_traces(dep_file)
data_dict_dep = hcf.load_traces(dep_file2)

ch1, hyp_sweeps, dep_sweeps, sweep_len, bl, message1, message2 = ncf.load_traces(hyp_file, dep_file, channel)
baseline, onset, hyp_onset = ncf.get_baseline(ch1, bl, dep_sweeps, hyp_sweeps, 'cell')

inj, spike_counts, max_spikes, peaks, first_spike, IO_slope, ch1, dep_sweeps = ncf.apcount(ch1, dep_sweeps, hyp_sweeps, dep_file, channel, [2])    
#not working with spikes = 0
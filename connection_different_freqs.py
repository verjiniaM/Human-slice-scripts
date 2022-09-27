
import os 
import pandas as pd
import numpy as np
import human_characterisation_functions as hcf
import connection_parameters as con_param
import fine_screen_diff_freq_functions as fine
import sorting_functions as sort
import funcs_for_results_tables as results

human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
patcher = 'Verji'
OP =  'OP220615' 
diff_freq_file_indices = [6,7]
pre_chans, post_chans = [4,5], [7,6]


def get_post_cell_responses_diff_freqs(human_dir, OP, patcher, diff_freq_file_indices, pre_chans, post_chans):
    work_dir, filenames, indices_dict, slice_names = results.get_OP_metadata(human_dir, OP, patcher)

    sort.to_json (work_dir, OP, '_diff_freq_file_indices.json', diff_freq_file_indices, pre_chans, post_chans)
    diff_freqs_dict = sort.from_json(work_dir, OP, '_diff_freq_file_indices.json')

    con_data = pd.DataFrame(columns = ['OP', ' patcher', 'fn', 'slice', 'chan_pre', 'chan_post', 'hz', 'Vm pre', 'Vm post', 
    'Amp1',	'Amp2',	'Amp3',	'Amp4',	'Lat1',	'Lat2',	'Lat3',	'Lat4', 'num excluded swps', 'comments'])


#needs debuggging and some tests
#last time checked, pree_peaks off with around 500 data points
#stim wondow not correctly read
    for i, indx in enumerate(diff_freqs_dict['indices']):
        fn_freqs = work_dir + filenames[indx]
        pre_cell = diff_freqs_dict['pre_chans'][i]
        post_cell = diff_freqs_dict['post_chans'][i]
        slic = slice_names[indx]

        hzs = ['200Hz', '100Hz', '50Hz', '20Hz', '10Hz']
        for hz in hzs:
            vm_sweep0, pa_sweep0, mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = \
                con_param.get_pre_aps_diff_freqs (fn_freqs, pre_cell, post_cell, hz)
            post_peaks = con_param.find_postsynaptic_peaks(post_window, preAPs_shifted)
            post_local_baseline = con_param.get_psp_baselines(post_window,preAPs_shifted)
            onsets = con_param.get_onsets(preAPs_shifted, post_window, post_peaks, post_local_baseline)
            latency = con_param.latencies(onsets, preAPs_shifted)

            data_to_add = pd.DataFrame({'OP':OP[:-1], 'patcher': patcher, 'fn':filenames[indx], 'slice':slic, 
            'chan_pre' : str(pre_cell), 'chan_post': str(post_cell), 'hz': hz, 'Vm pre': vm_sweep0, 'pA post': pa_sweep0, 
            'Amp1':,	'Amp2',	'Amp3',	'Amp4',	'Lat1',	'Lat2',	'Lat3',	'Lat4'
         })
   


#%%
file_num = 6

fine_con_screen = work_dir + OP + '/' + filenames[file_num-1] 

pre_cell = 2
post_cell = 4




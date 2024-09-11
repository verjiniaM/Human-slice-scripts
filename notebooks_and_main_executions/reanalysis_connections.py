import pandas as pd
import glob
import intrinsic_props_and_connectivity.funcs_sorting as sort
import src.intrinsic_props_and_connectivity.funcs_human_characterisation as hcf
import numpy as np
import intrinsic_props_and_connectivity.funcs_plotting_raw_traces as funcs_plotting_raw_traces
import src.intrinsic_props_and_connectivity.funcs_con_screen as con_param


human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
results_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/'

exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0]) 

OP = 'OP240507' 
patcher = 'Verji'

work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher, 'old')
dir_plots = sort.make_dir_if_not_existing (work_dir, 'connection_analysis')

con_screen_json = sort.get_json_meta_connect_all(human_dir, OP, patcher)

for a, i in enumerate(con_screen_json[0]['con_screen_file']):
    %matplotlib qt
    funcs_plotting_raw_traces.plot_con_screen_all(work_dir + filenames[i], con_screen_json[0]['active_chans'][a])
    %matplotlib


#%%
# analyse connection properties

#sort.to_json(work_dir, OP, '_con_screen_only.json', con_screen_json)
con_screen_json = sort.get_json_meta_connect_all(human_dir, OP, patcher)

con_data = pd.DataFrame(columns = ['OP', 'fn', 'slice', 'day', 'connection_ID', 'repatch', 'treatment', 
'hrs_after_OP', 'hrs_incubation','chan_pre', 'chan_post', 'Vm pre', 'Vm post', 
'Amp 1','Amp 2','Amp 3', 'Amp 4',	'Lat1',	'Lat2',	'Lat3',	'Lat4', 'num excluded swps', 'comments'])

for i, indx in enumerate(con_screen_json[0]['con_screen_file']):
    con_screen_file = work_dir + filenames[indx]

    if len(con_screen_json[0]['pre_chans']) == 0:
        con_data.to_excel(work_dir + '/connection_analysis/' + OP + '_connections.xlsx', index=False) 
        print('no connections in this dataset')
        continue

    pre_cells = con_screen_json[0]['pre_chans'][i]
    post_cells = con_screen_json[0]['post_chans'][i]
    
    time_after_op = sort.get_time_after_OP(con_screen_file, exp_view.cortex_out[exp_view['OP'] == OP].tolist()[0])
    
    exclude = ''
    for j, pre_cell in enumerate(pre_cells):
        post_cell = post_cells[j]
        slic = slice_names[indx]
        if slic in con_screen_json[0]['slices']:
            treatment = con_screen_json[0]['treatment'][con_screen_json[0]['slices'].index(slic)]
        else:
            treatment  = ' '
        con_ID = hcf.get_connection_ID(con_screen_file, slic, pre_cell, post_cell)
        day = 'D1'
        if slic[-2:] == 'D2': 
            day = 'D2'

        pre_signal, es, vm_pre = con_param.presynaptic_screen(con_screen_file, pre_cell)
        post_signal, vm_post = con_param.postsynaptic_screen(con_screen_file, post_cell, es)
        if (np.array(vm_post)).size == 0:
            print('QC not passed!!')
            exclude = 'all'
            es2 = []
            post_signal, vm_post = con_param.postsynaptic_screen(con_screen_file, post_cell, es2)
            continue
        mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = \
            con_param.get_analysis_window(pre_signal, post_signal)
        pre_signal, post_signal = con_param.remove_sweeps_with_POSTsynAPs(pre_signal, post_signal, preAPs)
        post_peaks = con_param.find_postsynaptic_peaks(post_window, preAPs_shifted)
        post_local_baseline = con_param.get_psp_baselines(post_window,preAPs_shifted)
        onsets = con_param.get_onsets(preAPs_shifted, post_window, post_peaks, post_local_baseline)
        latency = con_param.latencies(onsets, preAPs_shifted)
        amps = con_param.get_amps(post_peaks, post_local_baseline)

        df_add = pd.DataFrame({'OP': OP, 'fn': filenames[indx], 'slice': slic, 'day':day, 'treatment': treatment, 
        'connection_ID' : con_ID, 'hrs_after_OP' : time_after_op,
        'chan_pre': pre_cell, 'chan_post': post_cell, 'Vm pre' :vm_pre, 'Vm post': vm_post,
        'Amp 1': amps[0], 'Amp 2': amps[1],	'Amp 3': amps[2],	'Amp 4': amps[3],	
            'Lat1': latency[0][0],	'Lat2' : latency[1][0],	'Lat3': latency[2][0],	'Lat4': latency[3][0], 
        'num excluded swps': len(es), 'comments': exclude}, index=[0])
        con_data = pd.concat([con_data.loc[:], df_add]).reset_index(drop=True)

        #plotting
        funcs_plotting_raw_traces.plot_connection_window(con_screen_file, pre_cell, post_cell, pre_window, \
            post_window, preAPs_shifted, post_signal, onsets, preAPs, post_peaks, post_local_baseline)
        funcs_plotting_raw_traces.plot_post_cell(con_screen_file, pre_cell, post_cell)

con_data.to_excel(work_dir + '/connection_analysis/' + OP + '_connections.xlsx', index=False) 

#%%
 
con_data_VC = pd.DataFrame(columns = ['OP', 'fn', 'slice', 'day', 'connection_ID', 'repatch', 'treatment', 
'hrs_after_OP', 'hrs_incubation','chan_pre', 'chan_post', 'Vm pre', 'Holding post', 
'Amp 1',	'Amp 2',	'Amp 3',	'Amp 4',	'Lat1',	'Lat2',	'Lat3',	'Lat4', 'num excluded swps', 'comments'])

for i, indx in enumerate(con_screen_json[1]['con_screen_IC_file_indices']):
    con_screen_file_IC = work_dir + filenames[indx]
    
    pre_cells = con_screen_json[1]['pre_chans_IC'][i]
    post_cells = con_screen_json[1]['post_chans_VC'][i]
    if pre_cells == [] or post_cells == []:
        continue
    time_after_op = sort.get_time_after_OP(con_screen_file_IC, exp_view.cortex_out[exp_view['OP'] == OP].tolist()[0])
    
    exclude = ''
    for j, pre_cell in enumerate(pre_cells):
        post_cell = post_cells[j]
        slic = slice_names[indx]
        if slic in con_screen_json[0]['slices']:
            treatment = con_screen_json[0]['treatment'][con_screen_json[0]['slices'].index(slic)]
        else:
            treatment  = ' '

        con_ID = hcf.get_connection_ID (con_screen_file_IC, slic, pre_cell, post_cell)

        day = 'D1'
        if slic[-2:] == 'D2': 
            day = 'D2'

        pre_sig, es, vm_pre = con_param.presynaptic_screen_IC(con_screen_file_IC, pre_cell)
        post_sig, holding_post = con_param.postsynaptic_screen_VC (con_screen_file_IC, post_cell, es)
        if (np.array(holding_post)).size == 0:
            print('QC not passed!!')
            exclude = 'all'
            es2 = []
            post_sig, holding_post = con_param.postsynaptic_screen_VC (con_screen_file_IC, post_cell, es)
            continue
        mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = con_param.get_analysis_window_VC(pre_sig, post_sig)
        post_peaks = con_param.find_postsynaptic_peaks_VC(post_window, preAPs_shifted)
        bl = con_param.get_psp_baselines(post_window,preAPs_shifted)
        onsets = con_param.get_onsets_VC(preAPs_shifted, post_window, post_peaks, bl)
        latency = con_param.latencies(onsets, preAPs_shifted)
        amps = con_param.get_amps_VC(post_peaks, bl)

        df_add = pd.DataFrame({'OP': OP, 'fn': filenames[indx], 'slice': slic, 'day':day, 'treatment': treatment, 
        'connection_ID' : con_ID, 'hrs_after_OP' : time_after_op,
        'chan_pre': pre_cell, 'chan_post': post_cell, 'Vm pre' :vm_pre, 'Holding post': holding_post,
        'Amp 1': amps[0], 'Amp 2': amps[1],	'Amp 3': amps[2],	'Amp 4': amps[3],	
            'Lat1': latency[0][0],	'Lat2' : latency[1][0],	'Lat3': latency[2][0],	'Lat4': latency[3][0], 
        'num excluded swps': len(es), 'comments': exclude}, index=[0])
        con_data_VC = pd.concat([con_data_VC.loc[:], df_add]).reset_index(drop=True)

        #plotting
        funcs_plotting_raw_traces.plot_connection_window_VC(con_screen_file_IC, pre_cell, post_cell, pre_window, \
                post_window, preAPs_shifted, post_sig, onsets, preAPs, post_peaks, bl)
        funcs_plotting_raw_traces.plot_post_cell_VC(con_screen_file_IC, pre_cell, post_cell)

con_data_VC.to_excel(work_dir + '/connection_analysis/' + OP + '_connected_cell_properties_post_in_VC.xlsx', index=False) 



# %%
import src.intrinsic_props_and_connectivity.funcs_plot_intrinsic_props as plotting_funcs
import pandas as pd 




df = pd.read_excel(work_dir + '/connection_analysis/' + OP + '_connections.xlsx')
df = df[df['repatch'] == 'yes']
op_color_dict = plotting_funcs.get_op_color_dict(df)
plotting_funcs.plot_connect_amplitude(df, 'repatch', op_color_dict,  results_ = 'amp',
destination_dir = work_dir + 'connection_analysis/')

# %%

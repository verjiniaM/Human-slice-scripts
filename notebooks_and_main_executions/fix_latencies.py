import pandas as pd
import ephys_analysis.funcs_sorting as sort
import ephys_analysis.funcs_human_characterisation as hcf
import numpy as np
import ephys_analysis.funcs_plotting_raw_traces as funcs_plotting_raw_traces
import ephys_analysis.funcs_con_screen as con_param
import PyQt5

# use this script to check all the traces that have an artifact or
# the EPSP parameters haven't been estimated correctly
# because small events or anything else


def add_con_ID_match(df):

    df = df.copy()
    patcher_dict = {'Verji': 'vm',
                    'Rosie': 'rs'}
    if 'con_ID_match' not in df.columns:
        con_IDs_match = []
        for i, OP in enumerate(df.OP):
            patcher_prefix = patcher_dict[df.patcher[i]]
            slic = df.slice[i]
            pre_cell = df.chan_pre[i]
            post_cell = df.chan_post[i]

            con_IDs_match.append(f'{patcher_prefix}{OP[2:]}{slic}c{pre_cell}#{post_cell}')
        df.insert(0, 'con_ID_match', con_IDs_match)
    else:
        return df
    return df
# useful when you want to add/ check connectivity for a given surgery
# plots interactively each entered connectiity file and lets user inspect\

#%% # LOAD the DATA
work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/'
save_plots_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/connectivity_paper/plots/figs_draft1/'
VM_data_dir = f'{work_dir}results/human/connectivity_paper/data/'
human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'

df_all_connections = pd.read_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+ \
                            '/connectivity_paper/data/connection_properties_cross_checked_with_intr.csv')
df_all_connections = df_all_connections[df_all_connections.treatment != 'fixed']
df_slice_newest = pd.read_csv(f'{VM_data_dir}slice_cons_with_failure_rates_Nov_2026.csv')
df_repatch_newest = pd.read_csv(f'{VM_data_dir}repatched_connections_looser_QC_Nov_2026.csv')

df_all_connections = add_con_ID_match(df_all_connections)
df_slice_newest = add_con_ID_match(df_slice_newest)
df_repatch_newest = add_con_ID_match(df_repatch_newest)

cons_all = set(df_all_connections.con_ID_match.unique())
cons_slice = set(df_slice_newest.con_ID_match.unique())
cons_repatch = set(df_repatch_newest.con_ID_match.unique())

cons_exclude = cons_slice | cons_repatch
df_all_cons_check = df_all_connections[~df_all_connections['con_ID_match'].isin(cons_exclude)].reset_index(drop=True)

# keep only where RS2 and treatment TTX
df_all_cons_check = df_all_cons_check[~((df_all_cons_check['day'] == 'D1') & \
                                        (df_all_cons_check['treatment'].isin(['Ctrl', 'high K'])))]
df_all_cons_check = df_all_cons_check.sort_values(['OP', 'day', 'slice', 'connection_ID'])
df_all_cons_check = df_all_cons_check.loc[:, ~df_all_cons_check.columns.str.contains('Unnamed', case=False)].reset_index(drop=True)

# df_all_cons_check.to_excel(f'{VM_data_dir}/manual_fixed_latencies/manual_check_cons_TTX_and_D2_orig.xlsx')

cells_exclude = ['vm211123S2_D2c5#2', 'vm211123S2_D2c5#7', 'vm211209S3c5#7',
                 'vm211209S4c3#4', 'vm211209S4_D2c8#2', 'vm211209S5_D2c7#6',
                 'vm220602S3c8#3', 'vm220602S3_D2c8#7', 'vm220811S2_D2c1#2',
                 'vm220811S3_D2c5#2', 'vm220914S2c3#4', 'vm221024S3_D2c4#7',
                 'vm221027S3_D2c3#1', 'vm221027S3_D2c3#2', 'vm230523S1c2#3',
                 'vm230523S1c4#5', 'vm230808S1c1#3', 'vm230808S1c1#4',
                 'vm230808S1c2#4', 'vm230808S2_D2c6#7', 
                ' vm230810S2c2#3', 'vm230810S2c2#4', 'vm230810S2c2#6',
                'vm230810S2c2#7', 'vm230810S2c3#2', 'vm230810S2c4#2', 
                'vm230810S2c5#3', 'vm230810S2c6#2', 'vm230810S2c7#2',
                'vm240117S1_D2c5#8', 'vm240117S1_D2c8#6', 'vm240321S1_D2c7#8',
                'vm240411S2_D2c3#4', 'vm240417S1_D2c7#6', 'vm240417S2_D2c8#2', 
                'vm240503S1_D2c5#2', 'vm240503S2_D2c5#2', 'vm240503S3_D2c6#2',
                'vm240503S3_D2c2#5', 'vm250207S2_D2c1#5', 'vm240926S2_D2c1#3']
df_all_cons_check = df_all_cons_check[~df_all_cons_check['con_ID_match'].isin(cells_exclude)].reset_index(drop=True)


### !!!!
# cons that might need an extra check vm220602S3c8#6 and vm220602S3_D2c7#3 - it's a repatch but not sure they are connections

#%%

df_check = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/connectivity_paper/data/fix_repatch_latencies_last_time.xlsx')

patcher_dict = {'Verji': 'vm',
                'Rosie': 'rs'}

# # for OP in OP_list:
# # # import importlib
# # # importlib.reload(funcs_plotting_raw_traces)
# # #     i = 0
# #     i_s = df_check[df_check.OP == OP].index



# 'vm230420S2_D2c7#6' manual plot and lat/ amp, 'vm230420S1_D2c7#8'

# 'vm240417S2_D2c7#8' check la

for i, OP in enumerate(df_check.OP):
    patcher = df_check.patcher[i]
    fn = df_check.fn[i]
    ch_pre = df_check.chan_pre[i]
    ch_post = df_check.chan_post[i]
    work_dir, _, _, _ = sort.get_OP_metadata(human_dir, OP, patcher, 'old')
    print(f'plotting {fn} {OP} {ch_pre} {ch_post}')
    %matplotlib qt
    funcs_plotting_raw_traces.plot_con_screen_click(f'{work_dir}{fn}', [ch_pre, ch_post])
    %matplotlib


#%% 
# AFTER manually checking which
## 
# PLOT AND Calcualte new values

df_manual_correction = pd.read_excel(f'{VM_data_dir}/manual_fixed_latencies/manual_check_cons_TTX_and_D2.xlsx')
df_manual_correction[['peak_1', 'peak_2', 'peak_3', 'peak_4']] = \
            df_manual_correction[['peak_1', 'peak_2', 'peak_3', 'peak_4']].astype('Int64')

con_data = pd.DataFrame({})
cells_exclude = []
for i, OP in enumerate(df_manual_correction.OP):
    if 'exclude' in df_manual_correction.comments[i]:
        print(f'excluding {df_manual_correction.con_ID_match[i]}')
        cells_exclude.append(df_manual_correction.con_ID_match[i])
        continue
    fn = df_manual_correction.fn[i]
    patcher = df_manual_correction.patcher[i]
    pre_cell = df_manual_correction.chan_pre[i]
    post_cell = df_manual_correction.chan_post[i]
    slic = df_manual_correction.slice[i]
    work_dir, _, _, _ = sort.get_OP_metadata(human_dir, OP, patcher, 'old')
    con_screen_file = work_dir + fn

    con_ID = hcf.get_connection_ID(con_screen_file, slic, pre_cell, post_cell)
    id_pre = f'{patcher_dict[patcher]}{OP[2:]}{slic}{pre_cell}'
    id_post = f'{patcher_dict[patcher]}{OP[2:]}{slic}{post_cell}'

    pre_signal, es, vm_pre = con_param.presynaptic_screen(con_screen_file, pre_cell)
    post_signal, vm_post = con_param.postsynaptic_screen(con_screen_file, post_cell, es)

    mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = \
        con_param.get_analysis_window(pre_signal, post_signal)
    if 'skip' in df_manual_correction.comments[i]:
        print(f'skipping because {df_manual_correction.loc[i, ['comments']]}')
        continue

    if df_manual_correction.loc[i, ['pre_peak1', 'pre_peak2', 'pre_peak3', 'pre_peak4']].notna().any():
        # mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = \
        #         con_param.get_analysis_window(pre_signal, post_signal, True)
        AP_peaks = df_manual_correction.loc[i, ['pre_peak1', 'pre_peak2', 'pre_peak3', 'pre_peak4']].astype(int).to_numpy()
        num_APs = len(AP_peaks)
        AP_peak_vals = [mean_pre[pk_indx] for pk_indx in AP_peaks]

        preAPs = (AP_peaks, {'peak_heights': np.array(AP_peak_vals)})
        pre_window = mean_pre[preAPs[0][0]-750:preAPs[0][num_APs-1]+750]
        post_window = mean_post[preAPs[0][0]-750:preAPs[0][num_APs-1]+750]
        APs_shifted = [(749 + AP_peaks[i] - AP_peaks[0]) for i in range(num_APs)]
        AP_shifted_peak_vals = [pre_window[shifted_peak] for shifted_peak in APs_shifted]
        preAPs_shifted = (APs_shifted, {'peak_heights': np.array(AP_shifted_peak_vals)})

    if df_manual_correction.loc[i, ['pre_peak1', 'pre_peak2', 'pre_peak3', 'pre_peak4']].notna().any():
        post_peaks_indx = df_manual_correction.loc[i, ['peak_1', 'peak_2', 'peak_3', 'peak_4']].dropna().astype(int).tolist()
        # post_peaks_indx = [int(v) for v in post_peaks_indx]
        if len(post_peaks_indx) < 4:
            preAPs_shifted = (preAPs_shifted[0][:len(post_peaks_indx)], 
            {'peak_heights': preAPs_shifted[1]['peak_heights'][:len(post_peaks_indx)]})

        post_peak_values =  mean_post[post_peaks_indx]
        post_peaks_indx = [int(post_peak - preAPs[0][0] + 750) for post_peak in post_peaks_indx]
        post_peaks = np.column_stack((post_peak_values, post_peaks_indx))
    else:
        post_peaks = con_param.find_postsynaptic_peaks(post_window, preAPs_shifted)
    
    post_local_baseline = con_param.get_psp_baselines(post_window, preAPs_shifted)
    onsets = con_param.get_onsets(preAPs_shifted, post_window, post_peaks, post_local_baseline)
    if 'lat' in df_manual_correction.comments[i]:
        latency = df_manual_correction.loc[3, ['Lat1', 'Lat2', 'Lat3', 'Lat4']].to_numpy().reshape(-1, 1)
    else:
        latency = con_param.latencies(onsets, preAPs_shifted)
    amps = con_param.get_amps(post_peaks, post_local_baseline)

    df_add = pd.DataFrame({'OP': OP, 'fn': fn, 'slice': slic, 'day':df_manual_correction.day[i],
                           'patcher': patcher, 'chan_pre' : pre_cell, 'chan_post': post_cell,  
                            'treatment': df_manual_correction.treatment[i], 'connection_ID' : con_ID, 
                            'Amp 1': amps[0], 'Amp 2': amps[1],	'Amp 3': amps[2],	'Amp 4': amps[3],
                            'Lat1': latency[0][0],	'Lat2' : latency[1][0],	'Lat3': latency[2][0],	'Lat4': latency[3][0],
                            'comments': 'manual_fix_latency', 
                            'cell_IDs_match_pre': id_pre, 'cell_IDs_match_post': id_post},index=[0])
    con_data = pd.concat([con_data.loc[:], df_add]).reset_index(drop=True)

    #plotting
    funcs_plotting_raw_traces.plot_connection_window(con_screen_file, pre_cell, post_cell, pre_window, \
        post_window, preAPs_shifted, post_signal, onsets, preAPs, post_peaks, post_local_baseline, \
            save_dir=f'{VM_data_dir}/second_analysis_EPSPs/')
    funcs_plotting_raw_traces.plot_post_cell(con_screen_file, pre_cell, post_cell, \
                                             save_dir=f'{VM_data_dir}/second_analysis_EPSPs/')

# con_data.to_excel(f'{VM_data_dir}/con_analysis_clean_all/connections_fixed_lats.xlsx', index=False)


#%%
# matching connection_IDs

# add_connection_ID_match
df_all_cons_check = df_all_cons_check[~df_all_cons_check['con_ID_match'].isin(cells_exclude)].reset_index(drop=True)

if 'con_ID_match' not in con_data.columns:
    con_IDs_match = []
    for i, OP in enumerate(con_data.OP):
        patcher_prefix = patcher_dict[con_data.patcher[i]]
        slic = con_data.slice[i]
        pre_cell = con_data.chan_pre[i]
        post_cell = con_data.chan_post[i]

        con_IDs_match.append(f'{patcher_prefix}{OP[2:]}{slic}c{pre_cell}#{post_cell}')
    con_data.insert(0, 'con_ID_match', con_IDs_match)



cells_exclude_tripple_check = ['vm220602S3c8#6', 'vm220602S3_D2c7#3', 'vm220602S3c8#7']

# exclude cells that manually fixed
print(len(df_all_cons_check))
df_all_filtered = df_all_cons_check[~df_all_cons_check['con_ID_match'].isin(con_data['con_ID_match'])]
print(len(df_all_filtered))
print(f'Adding {len(con_data)} cells')
df_combined = pd.concat([df_all_filtered, con_data], ignore_index=True)
df_combined = df_combined[~df_combined['con_ID_match'].isin(cells_exclude_tripple_check)]
print(f'combined df with len {len(df_combined)}')

# change the repatch status
print(df_combined[df_combined.repatch=='yes'])
not_a_repatch  = ['vm220602S3_D2c7#6', 'vm230523S1c5#7', 'vm240321S2_D2c8#1', 'vm240507S3_D2c8#6']

index_change_repatch = df_combined[df_combined.con_ID_match.isin(not_a_repatch)].index
df_combined.loc[index_change_repatch,'repatch'] = 'no'
print(df_combined[df_combined.repatch=='yes'])

####
# ADD DATATYPE COLUMNS
df_combined = df_combined.sort_values(['patcher', 'OP', 'fn']).reset_index(drop=True)
df_combined.insert(len(df_combined.columns), 'dataset', 'TTX; inc_only, unchecked')
index_repatch_yes = df_combined[df_combined.repatch =='yes'].index
df_combined.loc[index_repatch_yes,'dataset'] = 'repatch'

df_repatch_newest.insert(len(df_repatch_newest.columns), 'dataset', 'repatch')
df_slice_newest.insert(len(df_slice_newest.columns), 'dataset', 'slice')

df_to_save_all_cons = pd.concat([df_combined, df_repatch_newest, df_slice_newest])
# df_combined.to_csv(f'{VM_data_dir}connection_properties_cross_checked_with_intr_QC_checked.csv')
# df_to_save_all_cons.to_csv(f'{VM_data_dir}all_connections_QC_checked_complete.csv')





# df_crazy = df_combined[(df_combined.Lat1 < 0) |
#         (df_combined.Lat2 < 0) |
#         (df_combined.Lat3 < 0) |
#         (df_combined.Lat4 < 0)]

# df_crazy = df_crazy.reset_index(drop = True)
# for i, OP in enumerate(df_crazy.OP):
    
#     patcher = df_crazy.patcher[i]
#     fn = df_crazy.fn[i]
#     ch_pre = df_crazy.chan_pre[i]
#     ch_post = df_crazy.chan_post[i]
#     work_dir, _, _, _ = sort.get_OP_metadata(human_dir, OP, patcher, 'old')
#     print(f'plotting {fn} {OP} {ch_pre} {ch_post}')
#     %matplotlib qt
#     funcs_plotting_raw_traces.plot_con_screen_click(f'{work_dir}{fn}', [ch_pre, ch_post])
#     %matplotlib


df_crazy = pd.read_excel(f'{VM_data_dir}/con_analysis_clean_all/fix_latencies_2_copy.xlsx')
print(len(df_combined))
df_all_filtered2 = df_combined[~df_combined['con_ID_match'].isin(df_crazy['con_ID_match'])]
print(len(df_all_filtered2))
print(f'Adding {len(df_crazy)} cells')
df_combined_final = pd.concat([df_all_filtered2, df_crazy], ignore_index=True)
print(f'combined df with len {len(df_combined_final)}')


# get the repatch df ready 

df_repatch_newest_ =  add_con_ID_match(df_repatch_newest_)

# df_repatch_add = df_repatch_newest_[df_repatch_newest_['comment'].str.contains('appearing|disappearing', case=False, na=False)].reset_index(drop = True)
# df_repatched_fixed = pd.concat([df_combined_final[df_combined_final['repatch'] == 'yes'], df_repatch_add])
# df_repatch_removed = df_repatch_newest_[~df_repatch_newest_['con_ID_match'].isin(df_repatched_fixed['con_ID_match'])].reset_index(drop=True)

# repatched_fixed = pd.read_excel(f'{VM_data_dir}con_analysis_clean_all/fix_latencies_3_repatch.xlsx')
# df_repatched_save = pd.concat([repatched_fixed, df_repatched_fixed])
# df_repatched_save.to_excel(f'{VM_data_dir}/con_analysis_clean_all/repatched_connections_looser_QC_Nov_2026.xlsx')

# df_slice_final = df_combined_final[df_combined_final['repatch'] == 'no']
# df_slice_final.to_excel(f'{VM_data_dir}/con_analysis_clean_all/slice_cons_with_failure_rates_Nov_2026.xlsx')


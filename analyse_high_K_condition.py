#%%
import sorting_functions as sort
import plotting_funcs
import human_characterisation_functions as hcf
import pandas as pd


human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
OP = 'OP221116' # and also OP220615 and OP220413
patcher = 'Verji'

work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)
file_list = sort.get_sorted_file_list(work_dir)
jsons = sort.get_json_files(file_list)
if OP + '_indices_dict.json' in jsons:
    indices_dict = sort.from_json(work_dir, OP, '_indices_dict.json')
else: 
    sort.to_json(work_dir, OP, '_indices_dict.json', indices_dict)

#creating a dir to save plots and data_tables (if not existing)
dir_plots = sort.make_dir_if_not_existing (work_dir, 'plots')
sort.make_dir_if_not_existing (work_dir, 'data_tables')

#check if the traces dir is empty and only then plot the mimiddle sweep for each filename
traces_folder =  os.path.join(dir_plots, "traces/")
if os.path.isdir(traces_folder) == 0 :
    for rec in range(len(filenames)):
        filename = work_dir + filenames[rec]
        plotting_funcs.plot_middle_sweep(filename)
else:
        print("skipping plotting")

active_chans_meta = sort.get_json_meta_high_K(human_dir, OP, patcher, '_meta_active_chans_high_K.json')
cortex_out_time = sort.get_datetime_from_input(active_chans_meta[0]['OP_time'][0])

df_vc = pd.DataFrame(columns=['filename', 'slice', 'cell_ch', 'hrs_after_OP', 'recording_time',
'circumstance', 'Rs', 'Rin']) 
#circumstance: 'before', 'wash in high K', 'puff high K' 'wash out'
for i, vc in enumerate(indices_dict['vc']):

    slic = slice_names[vc]

    filename_vc = work_dir + filenames[vc]
    time_after_op = sort.get_time_after_OP(filename_vc, cortex_out_time)

    active_channels = active_chans_meta[0]['active_chans'][i]
    
    time_after_op = sort.get_time_after_OP(filename_vc, cortex_out_time)
    rec_time = str(hcf.get_recording_time(filename_vc))
    Rs, Rin = hcf.get_access_resistance(filename_vc, active_channels)

    add_df = pd.DataFrame({'filename': filenames[vc], 'slice' : slic, 'cell_ch': active_channels,
    'hrs_after_OP' : time_after_op, 'recording_time': rec_time, 'Rs' : Rs, 'Rin': Rin })

    #plotting_funcs.plot_vc_holding (filename_vc, active_channels)

    df_vc = pd.concat([df_vc.loc[:], add_df]).reset_index(drop=True)
#df_vc.to_excel(work_dir + 'data_tables/' + OP + '_resistance.xlsx', index=False) 

df_AP_props = pd.DataFrame(columns=['filename', 'slice', 'cell_ch', 'hrs_after_OP', 'circumstance',
'recording_time', 'max_spikes', 'Rheobase', 'AP_heigth', 'TH', 
'max_depol', 'max_repol', 'membra_time_constant_tau', 'capacitance'])

for j, char in enumerate(indices_dict['freq analyse']):
    slic = slice_names[char]

    filename_char = work_dir + filenames[char]
    time_after_op = sort.get_time_after_OP(filename_char, cortex_out_time)

    active_channels = active_chans_meta[1]['active_chans'][j]
    
    time_after_op = sort.get_time_after_OP(filename_char, cortex_out_time)
    rec_time = str(hcf.get_recording_time(filename_char))

    meta_df = pd.DataFrame({'filename': filenames[char], 'slice' : slic, 'cell_ch': active_channels,
    'hrs_after_OP' : time_after_op, 'recording_time': rec_time })

    charact_params  = hcf.all_chracterization_params(filename_char, active_channels, 'full')
    df_char = pd.DataFrame.from_dict(charact_params)

    df_to_add = pd.concat([meta_df, df_char], axis = 1)
    df_AP_props = pd.concat([df_AP_props.loc[:], df_to_add]).reset_index(drop=True)

    #plotting_funcs.plots_for_charact_file(filename_char, active_channels, 'full')
#df_AP_props.to_excel(work_dir + 'data_tables/' + OP + '_AP_props.xlsx', index=False) 

df_resting = pd.DataFrame(columns=['filename', 'slice', 'cell_ch', 'hrs_after_OP', 'circumstance',
'recording_time', 'resting_potential'])

all_RMP_files = indices_dict['resting'] + indices_dict['resting long']
all_RMP_active_chans = active_chans_meta[2]['active_chans'] + active_chans_meta[3]['active_chans']
for k, RMP in enumerate(all_RMP_files):

    slic = slice_names[RMP]
    filename_resting = work_dir + filenames[RMP]
    time_after_op = sort.get_time_after_OP(filename_resting, cortex_out_time)

    active_channels = all_RMP_active_chans[k]
    
    time_after_op = sort.get_time_after_OP(filename_resting, cortex_out_time)
    rec_time = str(hcf.get_recording_time(filename_resting))

    RMPs  = hcf.get_RMP_over_time(filename_resting, active_channels)

    df_to_add = pd.DataFrame({'filename': filenames[RMP], 'slice' : slic, 'cell_ch': active_channels,
    'hrs_after_OP' : time_after_op, 'recording_time': rec_time, 'resting_potential': RMPs })

    df_resting = pd.concat([df_resting.loc[:], df_to_add]).reset_index(drop=True)

df_resting.sort_values(by = 'recording_time', axis=0, ascending=True, inplace=True)
#df_resting.to_excel(work_dir + 'data_tables/' + OP + '_RMPs_over_time.xlsx', index=False) 



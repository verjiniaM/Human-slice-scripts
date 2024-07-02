import funcs_sorting as sort
import funcs_sorting as sort
import pandas as pd
import funcs_human_characterisation as hcf
import funcs_plotting as plotting_funcs
import funcs_ as plot_intrinsic_props


def main():
    human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'

    OP = 'OP231207'
    patcher = 'Verji'
    tissue_source = 'Bielefeld'
    patient_age = 23

    get_data_high_K(human_dir, OP, patcher, tissue_source, patient_age)
    input('Press enter when AP props and VC data_tables have been manually checked')

    work_dir = sort.get_work_dir(human_dir, OP, patcher)
    vc_path = work_dir + 'data_tables/VC_high_K_' + OP + '.xlsx'
    AP_props_path = work_dir + 'data_tables/AP_props_high_K_' + OP + '.xlsx'
    RMP_path = work_dir + 'data_tables/RMP_high_K_' + OP + '.xlsx'
    spontan_path = work_dir + 'data_tables/spontan_high_K_' + OP + '_QC_measures.xlsx'
    ramp_path = work_dir + 'data_tables/ramp_high_K_' + OP + '.xlsx'

    check_vc(vc_path)
    exclude_cell_ID_only_in_second_df(AP_props_path , vc_path)
    exclude_cell_ID_only_in_second_df(AP_props_path , RMP_path)
    exclude_cell_ID_only_in_second_df(AP_props_path , spontan_path)
    exclude_cell_ID_only_in_second_df(AP_props_path, ramp_path)
    

def get_data_high_K(human_dir, OP, patcher, tissue_source, patient_age):

    work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)

    #creating a dir to save plots and data_tables (if not existing)
    dir_plots = sort.make_dir_if_not_existing (work_dir, 'plots')
    sort.make_dir_if_not_existing (work_dir, 'data_tables')

    sort.plot_trace_if_not_done(work_dir, dir_plots, filenames)
    
    active_chans_meta = sort.get_json_meta(human_dir, OP, patcher, '_meta_active_chans.json')[0]  
    cortex_out_time = sort.get_datetime_from_input(active_chans_meta['OP_time'][0])

    df_vc = pd.DataFrame(columns=['tissue_source', 'OP', 'patcher', 'patient_age','filename', 
    'slice', 'cell_ch', 'cell_ID', 'incubation_solution','hrs_incubation', 'recording_time','recording_in',
    'hrs_after_OP', 'K concentration', 'Rs', 'Rin'])
    for vc in sorted(indices_dict['vc'] + indices_dict['vc_end']):
        filename_vc = work_dir + filenames[vc]
        slic = slice_names[vc]
        indx_slice = active_chans_meta['slices'].index(slic)
        active_chans = active_chans_meta['active_chans'][indx_slice] #getting the active chans for this slice

        cell_IDs = hcf.get_new_cell_IDs(filename_vc, slic, active_chans, patcher)
        time_after_op = sort.get_time_after_OP(filename_vc, cortex_out_time)
        rec_time = str(hcf.get_recording_time(filename_vc))

        Rs, Rin = hcf.get_access_resistance(filename_vc, active_chans)
        df_add = pd.DataFrame({'tissue_source': tissue_source, 'OP': OP, 'patcher': patcher, 'patient_age': patient_age,
        'filename':filenames[vc], 'slice':slic, 'cell_ch':active_chans, 'cell_ID':cell_IDs, 'recording_time' :rec_time,
        'hrs_after_OP':time_after_op, 'Rs':Rs, 'Rin':Rin})
        plotting_funcs.plot_vc_holding(filename_vc, active_chans)
        
        df_vc = pd.concat([df_vc.loc[:], df_add]).reset_index(drop = True)
    df_vc.to_excel(work_dir + 'data_tables/VC_high_K_' + OP + '.xlsx', index=False) 

    df_resting = pd.DataFrame(columns=['tissue_source', 'OP', 'patcher', 'patient_age','filename', 
    'slice', 'cell_ch', 'cell_ID', 'incubation_solution','hrs_incubation', 'recording_time','recording_in',
    'hrs_after_OP', 'holding_minus_70_y_o_n', 'temperature', 'K concentration', 'resting_potential'])
    for RMP in sorted(indices_dict['resting'] + indices_dict['resting_long']):
        slic = slice_names[RMP]
        filename_resting = work_dir + filenames[RMP]

        indx_slice = active_chans_meta['slices'].index(slic)
        active_channels = active_chans_meta['active_chans'][indx_slice] #getting the active chans for this slic

        cell_IDs = hcf.get_new_cell_IDs(filename_resting, slic, active_channels, patcher)
        time_after_op = sort.get_time_after_OP(filename_resting, cortex_out_time)
        rec_time = str(hcf.get_recording_time(filename_resting))

        RMPs  = hcf.get_RMP_over_time(filename_resting, active_channels)

        df_to_add = pd.DataFrame({'tissue_source': tissue_source, 'OP': OP, 'patcher': patcher, 'patient_age': patient_age,
        'filename': filenames[RMP], 'slice' : slic,'cell_ch': active_channels, 'cell_ID':cell_IDs, 'recording_time':rec_time,
        'hrs_after_OP' : time_after_op, 'recording_time': rec_time, 'resting_potential': RMPs })

        df_resting = pd.concat([df_resting.loc[:], df_to_add]).reset_index(drop=True)
        
        sort.make_dir_if_not_existing(work_dir, 'plots/RMP/')
        plot_intrinsic_props.plot_spiking_when_RMP_increases(filename_resting, work_dir + 'plots/RMP/')

    df_resting.sort_values(by = 'recording_time', axis=0, ascending=True, inplace=True)
    df_resting.to_excel(work_dir + 'data_tables/RMP_high_K_' + OP + '.xlsx', index=False) 

    df_AP_props = pd.DataFrame(columns=['tissue_source', 'OP', 'patcher', 'patient_age','filename', 'slice', 'cell_ch', 'cell_ID', 
    'incubation_solution','hrs_incubation','holding_minus_70_y_o_n', 'temperature','recording_time','recording_in','hrs_after_OP', 
    'K concentration', 'max_spikes', 'Rheobase', 'AP_heigth', 'TH', 'max_depol', 'max_repol', 'membra_time_constant_tau','capacitance', 'Rheobse_ramp'])

    for j, char in enumerate(indices_dict['freq analyse']):
        
        filename_char = work_dir + filenames[char]
        slic = slice_names[char]
        indx_slice = active_chans_meta['slices'].index(slic)
        active_channels = active_chans_meta['active_chans'][indx_slice]
        
        cell_IDs = hcf.get_new_cell_IDs(filename_char, slic, active_channels, patcher)
        time_after_op = sort.get_time_after_OP(filename_char, cortex_out_time)
        rec_time = str(hcf.get_recording_time(filename_char))

        meta_df = pd.DataFrame({'tissue_source': tissue_source, 'OP': OP, 'patcher': patcher, 'patient_age': patient_age,
        'filename': filenames[char], 'slice' : slic, 'cell_ch': active_channels, 'cell_ID':cell_IDs,
        'hrs_after_OP' : time_after_op, 'recording_time': rec_time })

        charact_params  = hcf.all_chracterization_params(filename_char, active_channels, 'full')
        df_char = pd.DataFrame.from_dict(charact_params)

        df_to_add = pd.concat([meta_df, df_char], axis = 1)
        df_AP_props = pd.concat([df_AP_props.loc[:], df_to_add]).reset_index(drop=True)

        plotting_funcs.plots_for_charact_file(filename_char, active_channels, 'full')
    df_AP_props.to_excel(work_dir + 'data_tables/AP_props_high_K_' + OP + '.xlsx', index=False) 

    df_ramp = pd.DataFrame(columns=['tissue_source', 'OP', 'patcher', 'patient_age','filename', 'slice', 'cell_ch', 'cell_ID', 
    'incubation_solution','hrs_incubation', 'holding_minus_70_y_o_n', 'temperature','recording_time','recording_in','hrs_after_OP', 'K concentration', 'Rheobse_ramp'])
    for j, ramp in enumerate(indices_dict['ramp']):
        
        filename_ramp = work_dir + filenames[ramp]
        slic = slice_names[ramp]
        indx_slice = active_chans_meta['slices'].index(slic)
        active_channels = active_chans_meta['active_chans'][indx_slice]
        # active_channels = active_chans_ramp[j]
        cell_IDs = hcf.get_new_cell_IDs(filename_ramp, slic, active_channels, patcher)
        time_after_op = sort.get_time_after_OP(filename_ramp, cortex_out_time)
        rec_time = str(hcf.get_recording_time(filename_ramp))

        rheos, THs, THs_in_trace, swps = hcf.get_rheobase_from_ramp(filename_ramp, active_channels)

        df_to_add = pd.DataFrame({'tissue_source': tissue_source, 'OP': OP, 'patcher': patcher, 'patient_age': patient_age,
        'filename': filenames[ramp], 'slice' : slic, 'cell_ch': active_channels, 'cell_ID':cell_IDs,
        'hrs_after_OP' : time_after_op, 'recording_time': rec_time, 'Rheobse_ramp':rheos})

        df_ramp = pd.concat([df_ramp.loc[:], df_to_add]).reset_index(drop=True)
    df_ramp.to_excel(work_dir + 'data_tables/ramp_high_K_' + OP + '.xlsx', index=False) 

    record_sorting = pd.DataFrame()
    for i, u in enumerate(indices_dict['spontan']):
        filename_spontan = work_dir + filenames[u]
        slic = slice_names[u]
        indx_slice = active_chans_meta['slices'].index(slic)
        active_channels = active_chans_meta['active_chans'][indx_slice]

        cell_IDs = hcf.get_new_cell_IDs(filename_spontan, slic, active_channels, patcher)

        spontan_QC = hcf.rec_stability(filename_spontan, active_channels , 60)
        df_QC = pd.DataFrame(spontan_QC).T
        df_QC['cell_ch'] = df_QC.index 
        df_QC.insert(0, 'cell_ch', df_QC.pop('cell_ch'))
        df_QC.insert(0, 'slice', slic)
        df_QC.insert(0, 'cell_ID', cell_IDs)
        df_QC.insert(0, 'filename', filenames[u])
        df_QC.insert(0, 'patcher', patcher)
        df_QC.insert(0, 'OP', OP)

        record_sorting = pd.concat([record_sorting.loc[:], df_QC]).reset_index(drop=True)
    if len(record_sorting) > 0:
        mask = []
        for i in record_sorting.swps_to_analyse:
            if len(i) <= 0:
                mask.append(False)
                continue
            mask.append(True)
        record_sorting = record_sorting.loc[mask,:]
    record_sorting.to_excel(work_dir + 'data_tables/spontan_high_K_' + OP + '_QC_measures.xlsx', index=False) 

def check_vc(df_vc_path):
    ''' 
    checks that the Rs is not > 30
    where Rs > 30 adds a comment
    '''
    vc_df = pd.read_excel(df_vc_path)

    comments = list(range(0, len(vc_df)))
    vc_df.insert(len(vc_df.columns), 'comments', comments)

    mask = vc_df['Rs'] > 30
    vc_df.loc[mask, 'comments'] == 'exclude, Rs > 30'
    vc_df.to_excel(df_vc_path)


def exclude_cell_ID_only_in_second_df(AP_props_path , df_to_QC_path):
    '''checks for cell_IDs present only in AP_props_df
    !! which have been manually QC controlled!!
    keeps only data about these cells in df_to_QC'''
    AP_props_df = pd.read_excel(AP_props_path)
    df_to_QC = pd.read_excel(df_to_QC_path)

    if len(df_to_QC) == 0:
        return
    
    exclude_list = list(set(df_to_QC['cell_ID'].tolist()) - set(AP_props_df['cell_ID'].tolist()))

    for j in exclude_list:
        indx_RMP = df_to_QC[df_to_QC['cell_ID'] == j].index
        indx_AP = AP_props_df[AP_props_df['cell_ID'] == j].index
        AP_props_df.drop(indx_AP, axis=0, inplace = True)
        df_to_QC.drop(indx_RMP, axis=0, inplace=True)

    AP_props_df.reset_index(inplace = True,  drop = True)
    AP_props_df.to_excel(AP_props_path)
    df_to_QC.reset_index(inplace = True, drop = True)
    df_to_QC.to_excel(df_to_QC_path)


if __name__ == '__main__':
    main()








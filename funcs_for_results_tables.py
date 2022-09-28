import os
import human_characterisation_functions as hcf
import plotting_funcs
import pandas as pd
import sorting_functions as sort
import connection_parameters as con_param

def get_OP_metadata (human_dir, OP, patcher):
    OP_folder = OP + '/'
    if patcher == 'Verji': 
        work_dir = human_dir + 'data_verji/'+ OP_folder
    else:
        work_dir = human_dir + 'data_rosie/'+ OP_folder 

    file_list = sort.get_sorted_file_list(work_dir)
    df_rec= sort.get_lab_book(work_dir)
    filenames = sort.get_abf_files(file_list)
    slice_indx, def_slice_names, indices_dict = sort.sort_protocol_names (file_list, df_rec)
    slice_names = sort.fix_slice_names (def_slice_names, slice_indx)
    return work_dir, filenames, indices_dict, slice_names


def get_json_meta (human_dir, OP, patcher, out_fn): #  out_fn = '_meta_active_chans.json'
    work_dir, filenames, indices_dict, slice_names = get_OP_metadata(human_dir, OP, patcher)

    file_list = sort.get_sorted_file_list(work_dir)
    jsons = sort.get_json_files(file_list)

    if OP + out_fn in jsons:
        json_meta = sort.from_json(work_dir, OP, out_fn)
        return json_meta

    slice_names_dict = []
    vc_indices = indices_dict['vc']
    active_chans_all =[]
    slice_names_dict = []
    for i in vc_indices:
        active_channels = [int(item) for item in input('Used channels in ' + OP + ' ' + slice_names[i]).split()]
        active_chans_all.append(active_channels)
        slice_names_dict.append(slice_names[i])

    con_screen_indices = indices_dict['con_screen']
    pre_chans_all, post_chans_all = [], []
    for indx in con_screen_indices:
        pre_chans = [int(item) for item in input('Pre channels in ' + filenames[indx]).split()]
        post_chans = [int(item) for item in input('Post channels in ' + filenames[indx]).split()]
        pre_chans_all.append(pre_chans)
        post_chans_all.append(post_chans)

    sort.to_json (work_dir, OP, out_fn, con_screen_indices, pre_chans_all, post_chans_all, slice_names_unique, vc_indices, active_chans_all)
    json_meta = sort.from_json(work_dir, OP, out_fn)
    return json_meta

def get_intrinsic_properties_df(human_dir, OP, tissue_source, patcher, age, inj):
    '''
    saves a pandas dataframe with intrinsic cell peoperties for OP
    '''
    work_dir, filenames, indices_dict, slice_names = get_OP_metadata(human_dir, OP, patcher)

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

    #Correct indices if needed
    [print(key,':',value) for key, value in indices_dict.items()]
    if len(indices_dict['freq analyse']) != len(indices_dict['vc']): 
        print('Fix protocol names. Unequal number of VC and freq analyse protocols')
        input('Press enter when indices have been fixed')
        work_dir, filenames, indices_dict, slice_names = get_OP_metadata(human_dir, OP, patcher)
        [print(key,':',value) for key, value in indices_dict.items()]
        # index_vc_in = [int(item) for item in input('Vc files corresponding to characterization files for ' + OP +' (input with spaces in between)').split()]
        # #saved the original indices
        # indices_dict['vc_orig'] = indices_dict['vc']
        # indices_dict['vc'] = index_vc_in
        # # index_char = [int(item) for item in input('corresponding characterization files for ' + OP +' (input with spaces in between)').split()]
        # # indices_dict['freq analyse'] = index_char
 
    active_chans_meta = get_json_meta(human_dir, OP, patcher, '_meta_active_chans.json')

    #creating the dataframe
    df_OP = pd.DataFrame(columns=['slice', 'cell_ID','Rs', 'Rin', 'resting_potential', 'max_spikes', 'Rheobase', 'AP_heigth', 'TH', 'max_depol', 
    'max_repol', 'membra_time_constant_tau', 'capacitance'])

    for i in range(len(indices_dict['vc'])):
        vc = indices_dict['vc'][i]
        vm = indices_dict['vm'][i]
        char = indices_dict['freq analyse'][i]
        slic = slice_names[vc]

        filename_vc = work_dir + filenames[vc]
        filename_vm = work_dir + filenames[vm]
        filename_char = work_dir + filenames[char]
        filename_con_screen = work_dir + filenames[indices_dict['con_screen'][i]]

        active_channels = active_chans_meta[0]['active_chans'][i]
        
        cell_IDs = hcf.get_cell_IDs(filename_char, slic, active_channels)
        Rs, Rin = hcf.get_access_resistance(filename_vc, active_channels) 
        RMPs = hcf.get_RMP(filename_vm, active_channels)
        params1_df = pd.DataFrame({'slice' : slic, 'cell_ID': cell_IDs, 'Rs' : Rs, 
        'Rin': Rin, 'resting_potential': RMPs })

        charact_params  = hcf.all_chracterization_params(filename_char, active_channels, inj)
        df_char = pd.DataFrame.from_dict(charact_params)

        df_to_add = pd.concat([params1_df, df_char], axis = 1)
        df_OP = pd.concat([df_OP.loc[:], df_to_add]).reset_index(drop=True)

        #plotting function
        plotting_funcs.plot_vc_holding (filename_vc, active_channels)
        plotting_funcs.plots_for_charact_file(filename_char, active_channels, inj)
        plotting_funcs.plot_connect(filename_con_screen, active_channels)
     
    tissue = pd.Series(tissue_source).repeat(len(df_OP))
    OPs = pd.Series(OP).repeat(len(df_OP))
    researcher = pd.Series(patcher).repeat(len(df_OP))
    patient_age = pd.Series(age).repeat(len(df_OP))
    series_df = pd.DataFrame({'tissue_source': tissue, 'OP': OPs, 'patcher': researcher, 'patient_age': patient_age}).reset_index(drop=True)

    df_intrinsic = pd.concat([series_df, df_OP], axis = 1)

    df_intrinsic.to_excel(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.xlsx') 
    df_intrinsic.to_csv(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.csv')

    print('Intrinsic properties DataFrame for  ' + OP + 'saved successfully. ' + '\n' + 'Exclude recordings if necessary.')

def get_QC_access_resistance_df (human_dir, OP, patcher):
    work_dir, filenames, indices_dict, slice_names = get_OP_metadata(human_dir, OP, patcher)

    [print(key,':',value) for key, value in indices_dict.items()]
    if len(indices_dict['vc']) != len(indices_dict['vc_end']): 
        print('Fix protocol names. Unequal number of VC and freq analyse protocols')
    #     index_vc_in = [int(item) for item in input('Vc files corresponding to vc_end files for ' + OP +' (input with spaces in between)').split()]
    #     index_vc_end_in = [int(item) for item in input('Vc_end files ' + OP +' (input with spaces in between)').split()]
    # #saved the original indices
    #     indices_dict['vc_orig'] = indices_dict['vc']
    #     indices_dict['vc_end_org'] = indices_dict['vc_end']
    #     indices_dict['vc'] = index_vc_in
    #     indices_dict['vc_end'] = index_vc_end_in
        input('Press enter when indices have been fixed')
        work_dir, filenames, indices_dict, slice_names = get_OP_metadata(human_dir, OP, patcher)
        [print(key,':',value) for key, value in indices_dict.items()]

    active_chans_meta = get_json_meta(human_dir, OP, patcher, '_meta_active_chans.json')
    #date_frame for quality control - change in Rin, Rs 
    df_qc = pd.DataFrame(columns=['OP','patcher', 'filename', 'slice', 'cell_ch', 'Rs_start', 'Rin_start', 'Rs_end', 'Rin_end', 
    'chage_rs', 'change_rin'])

    for i in range(len(indices_dict['vc'])):
        vc = indices_dict['vc'][i]
        vc_end = indices_dict['vc_end'][i]
        slic = slice_names[vc]

        filename_vc = work_dir + filenames[vc]
        filename_vc_end = work_dir + filenames[vc_end]
       
        active_chans_meta[0]['active_chans'][i]        
        cell_IDs = hcf.get_cell_IDs(filename_vc, slic, active_channels)
        Rs, Rin = hcf.get_access_resistance(filename_vc, active_channels)
        Rs_end, Rin_end = hcf.get_access_resistance(filename_vc_end, active_channels)

        Rs_diff = [x - y for x,y in zip(Rs, Rs_end)]
        Rin_diff = [x - y for x,y in zip(Rin, Rin_end)]

        data_to_add = pd.DataFrame({'OP':OP[:-1], 'patcher':patcher, 'filename':filenames[vc], 'slice':slic, 
            'cell_ch':active_channels, 'Rs_start': Rs, 'Rin_start': Rin, 'Rs_end': Rs_end, 'Rin_end': Rin_end, 
            'chage_rs': Rs_diff, 'change_rin':Rin_diff})
        df_qc = pd.concat([df_qc.loc[:], data_to_add]).reset_index(drop=True)

    df_qc.to_excel(work_dir + 'data_tables/' + OP[:-1] + '_QC_measures_rs.xlsx') 
    df_qc.to_csv(work_dir + 'data_tables/' + OP[:-1] + '_QC_measures_rs.csv')

def get_con_params_df (human_dir, OP, patcher):
    work_dir, filenames, indices_dict, slice_names = get_OP_metadata(human_dir, OP, patcher)

    con_sccreen_connected_chans = get_json_meta(human_dir, OP, patcher, '_meta_active_chans.json')[1]

    con_data = pd.DataFrame(columns = ['fn', 'slice', 'chan_pre', 'chan_post', 'Vm pre', 'Vm post', 
    'Amp1',	'Amp2',	'Amp3',	'Amp4',	'Lat1',	'Lat2',	'Lat3',	'Lat4', 'num excluded swps', 'comments'])

    for i, indx in enumerate(con_sccreen_connected_chans['con_screen_file_indices']):
        con_screen_file = work_dir + filenames[indx]
        
        pre_cells = con_sccreen_connected_chans['pre_chans'][i]
        post_cells = con_sccreen_connected_chans['post_chans'][i]
        
        exclude = ''
        for j, pre_cell in enumerate(pre_cells):
            post_cell = post_cells[j]
            slic = slice_names[j]

            pre_signal, es, vm_pre = con_param.presynaptic_screen(con_screen_file, pre_cell)
            post_signal, vm_post = con_param.postsynaptic_screen(con_screen_file, post_cell, es)
            if vm_post == []:
                print('QC not passed!!')
                exclude = 'all'
                es2 = []
                post_signal, vm_post = con_param.postsynaptic_screen(con_screen_file, post_cell, es2)
            mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = \
                con_param.get_analysis_window(pre_signal, post_signal)
            pre_signal, post_signal = con_param.remove_sweeps_with_POSTsynAPs(pre_signal, post_signal, preAPs)
            post_peaks = con_param.find_postsynaptic_peaks(post_window, preAPs_shifted)
            post_local_baseline = con_param.get_psp_baselines(post_window,preAPs_shifted)
            onsets = con_param.get_onsets(preAPs_shifted, post_window, post_peaks, post_local_baseline)
            latency = con_param.latencies(onsets, preAPs_shifted)

            amps = []
            for u in range(4):
                amps.append(post_peaks[u][0] - post_local_baseline[0][0])

            df_add = pd.DataFrame({'fn': filenames[indx], 'slice': slic, 'chan_pre': pre_cell,
            'chan_post': post_cell, 'Vm pre' :vm_pre, 'Vm post': vm_post,
            'Amp1': amps[0], 'Amp2': amps[1],	'Amp3': amps[2],	'Amp4': amps[3],	
             'Lat1': latency[0][0],	'Lat2' : latency[1][0],	'Lat3': latency[2][0],	'Lat4': latency[3][0], 
            'num excluded swps': len(es), 'comments': exclude}, index=[0])
            con_data = pd.concat([con_data.loc[:], df_add]).reset_index(drop=True)

            #plotting
            plotting_funcs.plot_connection_window(con_screen_file, pre_cell, post_cell, pre_window, \
                post_window, preAPs_shifted, post_signal, onsets, preAPs, post_peaks, post_local_baseline)
            plotting_funcs.plot_post_cell(con_screen_file, pre_cell, post_cell)

    con_data.to_excel(work_dir + '/data_tables/' + OP + '_connected_cell_properties.xlsx') 




import os
import pandas as pd
import sorting_functions as sort
import human_characterisation_functions as hcf
import plotting_funcs
import connection_parameters as con_param

def get_intrinsic_properties_df(human_dir, OP, tissue_source, patcher, age, inj):
    '''
    saves a pandas dataframe with intrinsic cell peoperties for OP
    '''
    work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)

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
        work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)
        [print(key,':',value) for key, value in indices_dict.items()]
        # index_vc_in = [int(item) for item in input('Vc files corresponding to characterization files for ' + OP +' (input with spaces in between)').split()]
        # #saved the original indices
        # indices_dict['vc_orig'] = indices_dict['vc']
        # indices_dict['vc'] = index_vc_in
        # # index_char = [int(item) for item in input('corresponding characterization files for ' + OP +' (input with spaces in between)').split()]
        # # indices_dict['freq analyse'] = index_char
 
    active_chans_meta = sort.get_json_meta(human_dir, OP, patcher, '_meta_active_chans.json')[0]  
    cortex_out_time = sort.get_datetime_from_input(active_chans_meta['OP_time'])

    #creating the dataframe
    df_OP = pd.DataFrame(columns=['filename', 'slice', 'cell_ch', 'cell_ID', 'day', 'treatment', 
    'hrs_incubation', 'repatch', 'hrs_after_OP', 'Rs', 'Rin', 'resting_potential', 'max_spikes', 
    'Rheobase', 'AP_heigth', 'TH', 'max_depol', 'max_repol', 'membra_time_constant_tau', 'capacitance'])

    for i in range(len(indices_dict['vc'])):
        vc = indices_dict['vc'][i]
        vm = indices_dict['resting'][i]
        char = indices_dict['freq analyse'][i]
        slic = slice_names[vc]
        treatment = active_chans_meta['treatment'][i][0]
        day = 'D1'
        if slic[-2:] == 'D2': 
            day = 'D2'

        filename_vc = work_dir + filenames[vc]
        filename_vm = work_dir + filenames[vm]
        filename_char = work_dir + filenames[char]
        filename_con_screen = work_dir + filenames[indices_dict['con_screen'][i]]

        active_channels = active_chans_meta['active_chans'][i]
        
        cell_IDs = hcf.get_cell_IDs(filename_char, slic, active_channels)
        time_after_op = sort.get_time_after_OP(filename_char, cortex_out_time)
        Rs, Rin = hcf.get_access_resistance(filename_vc, active_channels) 
        RMPs = hcf.get_RMP(filename_vm, active_channels)
        params1_df = pd.DataFrame({'filename': filenames[char], 'slice' : slic, 'cell_ch': active_channels,
        cell_ID:'cell_IDs', 'day' : day , 'treatment': treatment, 'Rs' : Rs, 'Rin': Rin, 'resting_potential': RMPs })

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

    df_intrinsic.to_excel(work_dir + 'data_tables/' + OP + '_Intrinsic_and_synaptic_properties.xlsx') 
    #df_intrinsic.to_csv(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.csv')

    print('Intrinsic properties DataFrame for  ' + OP + 'saved successfully. ' + '\n' + 'Exclude recordings if necessary.')

def get_QC_access_resistance_df (human_dir, OP, patcher):
    work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)

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
        work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)
        [print(key,':',value) for key, value in indices_dict.items()]

    active_chans_meta = sort.get_json_meta(human_dir, OP, patcher, '_meta_active_chans.json')
    #date_frame for quality control - change in Rin, Rs 
    df_qc = pd.DataFrame(columns=['OP','patcher', 'filename', 'slice', 'cell_ID', 'cell_ch', 'Rs_start', 'Rin_start', 'Rs_end', 'Rin_end', 
    'change_rs', 'change_rin'])

    for i in range(len(indices_dict['vc'])):
        vc = indices_dict['vc'][i]
        vc_end = indices_dict['vc_end'][i]
        slic = slice_names[vc]

        filename_vc = work_dir + filenames[vc]
        filename_vc_end = work_dir + filenames[vc_end]
       
        active_channels = active_chans_meta[0]['active_chans'][i]        
        cell_IDs = hcf.get_cell_IDs(filename_vc, slic, active_channels)
        Rs, Rin = hcf.get_access_resistance(filename_vc, active_channels)
        Rs_end, Rin_end = hcf.get_access_resistance(filename_vc_end, active_channels)

        Rs_diff = [x - y for x,y in zip(Rs, Rs_end)]
        Rin_diff = [x - y for x,y in zip(Rin, Rin_end)]

        data_to_add = pd.DataFrame({'OP':OP[:-1], 'patcher':patcher, 'filename':filenames[vc], 'slice':slic, 
            'cell_ID': cell_IDs, 'cell_ch':active_channels, 'Rs_start': Rs, 'Rin_start': Rin, 
            'Rs_end': Rs_end, 'Rin_end': Rin_end, 'chagne_rs': Rs_diff, 'change_rin':Rin_diff})
        df_qc = pd.concat([df_qc.loc[:], data_to_add]).reset_index(drop=True)

    df_qc.to_excel(work_dir + 'data_tables/' + OP + '_QC_measures_rs.xlsx') 
    #df_qc.to_csv(work_dir + 'data_tables/' + OP[:-1] + '_QC_measures_rs.csv')

def get_con_params_df (human_dir, OP, patcher):
    work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)

    con_sccreen_connected_chans = sort.get_json_meta(human_dir, OP, patcher, '_meta_active_chans.json')[1]
    cortex_out_time = sort.get_datetime_from_input(con_sccreen_connected_chans['OP_time'])

    con_data = pd.DataFrame(columns = ['OP', 'fn', 'slice', 'day', 'cell_ID', 'repatch', 'treatment', 
    'hrs_after_OP', 'hrs_incubation','chan_pre', 'chan_post', 'Vm pre', 'Vm post', 
    'Amp1',	'Amp2',	'Amp3',	'Amp4',	'Lat1',	'Lat2',	'Lat3',	'Lat4', 'num excluded swps', 'comments'])

    for i, indx in enumerate(con_sccreen_connected_chans['con_screen_file_indices']):
        con_screen_file = work_dir + filenames[indx]
        
        pre_cells = con_sccreen_connected_chans['pre_chans'][i]
        post_cells = con_sccreen_connected_chans['post_chans'][i]
        treatment = con_sccreen_connected_chans['treatment'][i]
        time_after_op = sort.get_time_after_OP(con_screen_file, cortex_out_time)
        
        exclude = ''
        for j, pre_cell in enumerate(pre_cells):
            post_cell = post_cells[j]
            slic = slice_names[indx]
            cell_ID = hcf.get_cell_IDs(con_screen_file, slic, [pre_cell])
            day = 'D1'
            if slic[-2:] == 'D2': 
                day = 'D2'

            pre_signal, es, vm_pre = con_param.presynaptic_screen(con_screen_file, pre_cell)
            post_signal, vm_post = con_param.postsynaptic_screen(con_screen_file, post_cell, es)
            if vm_post.size == 0:
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

            df_add = pd.DataFrame({'OP': OP, 'fn': filenames[indx], 'slice': slic, 'day':day, 'treatment': treatment, 
            'cell_ID' : cell_ID, 'hrs_after_OP' : time_after_op,
            'chan_pre': pre_cell, 'chan_post': post_cell, 'Vm pre' :vm_pre, 'Vm post': vm_post,
            'Amp1': amps[0], 'Amp2': amps[1],	'Amp3': amps[2],	'Amp4': amps[3],	
             'Lat1': latency[0][0],	'Lat2' : latency[1][0],	'Lat3': latency[2][0],	'Lat4': latency[3][0], 
            'num excluded swps': len(es), 'comments': exclude}, index=[0])
            con_data = pd.concat([con_data.loc[:], df_add]).reset_index(drop=True)

            #plotting
            plotting_funcs.plot_connection_window(con_screen_file, pre_cell, post_cell, pre_window, \
                post_window, preAPs_shifted, post_signal, onsets, preAPs, post_peaks, post_local_baseline)
            #plotting_funcs.plot_post_cell(con_screen_file, pre_cell, post_cell)
            plotting_funcs.plot_post_cell_old_win(con_screen_file, pre_cell, post_cell)

    con_data.to_excel(work_dir + '/data_tables/' + OP + '_connected_cell_properties.xlsx') 

def get_spontan_QC(human_dir, OP, patcher):
    work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)

    active_chans_meta = sort.get_json_meta(human_dir, OP, patcher, '_meta_active_chans.json')

    #record_sorting = pd.DataFrame(columns = ['OP', 'patcher', 'filename', 'cell_ch','swps_to_analyse', 'swps_to_discard', 
    #'min_vals_each_swp', 'max_vals_each_swp', 'drift'])
    record_sorting = pd.DataFrame()
    for u in indices_dict['spontan']:
        filename_spontan = work_dir + filenames[u]
        slic = slice_names[u]   

        index_chans = active_chans_meta[0]['slices'].index(slic)
        active_channels = active_chans_meta[0]['active_chans'][index_chans]
        cell_IDs = hcf.get_cell_IDs(filename_spontan, slic, active_channels)

        spontan_QC = hcf.rec_stability (filename_spontan, active_channels , 60)
        df_QC = pd.DataFrame(spontan_QC).T
        df_QC['cell_ch'] = df_QC.index 
        df_QC.insert(0, 'cell_ch', df_QC.pop('cell_ch'))
        df_QC.insert(0, 'slice', slic)
        df_QC.insert(0, 'cell_ID', cell_IDs)
        df_QC.insert(0, 'filename', filenames[u])
        df_QC.insert(0, 'patcher', patcher)
        df_QC.insert(0, 'OP', OP[:-1])

        record_sorting = pd.concat([record_sorting.loc[:], df_QC]).reset_index(drop=True)

    record_sorting.to_excel(work_dir + 'data_tables/' + OP + '_QC_measures_spontan.xlsx') 
    #record_sorting.to_csv(work_dir + 'data_tables/' + OP + '_QC_measures_spontan.csv')


def get_minis_QC(human_dir, OP, patcher):
    work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)

    active_chans_meta = sort.get_json_meta(human_dir, OP, patcher, '_meta_active_chans.json')
    cortex_out_time = sort.get_datetime_from_input(active_chans_meta['OP_time'])

    record_sorting = pd.DataFrame(columns = ['OP','patcher', 'filename', 'cell_ch','swps_to_analyse', 'swps_to_discard', 
    'min_vals_each_swp', 'max_vals_each_swp', 'drift', 'Rs_cahnge', 'Rin_change', 'Rs_start'])

    for i, indx in enumerate(indices_dict['minis']):
        mini_file = work_dir + filenames[indx]
        slic = slice_names[indx]

        index_chans = active_chans_meta[2]['mini_slices'].index(slic)
        active_channels = active_chans_meta[2]['mini_chans'][index_chans]
        cell_IDs = hcf.get_cell_IDs(mini_file, slic, active_channels)
        time_after_op = sort.get_time_after_OP(mini_file, cortex_out_time)
        mini_QC = hcf.rec_stability(mini_file, active_channels , 60)

        if len(indices_dict['vc_mini']) == len(indices_dict['minis']):
            vc_indx = indices_dict['vc_mini'][i]
            Rs, Rin = hcf.get_access_resistance(work_dir + filenames[vc_indx], active_channels)
        
        if len(indices_dict['vc_mini_end']) == len(indices_dict['minis']):
            vc_end_indx = indices_dict['vc_mini_end'][i]
            Rs_end, Rin_end = hcf.get_access_resistance(work_dir + filenames[vc_end_indx], active_channels)
        
        Rs_diff = [x - y for x,y in zip(Rs, Rs_end)]
        Rin_diff = [x - y for x,y in zip(Rin, Rin_end)]

        df_QC = pd.DataFrame(mini_QC).T
        df_QC['cell_ch'] = df_QC.index 
        df_QC.insert(0, 'cell_ch', df_QC.pop('cell_ch'))
        df_QC.insert(0, 'slice', slic)
        df_QC.insert(0, 'cell_ID', cell_IDs)
        df_QC.insert(0, 'filename', filenames[indx])
        df_QC.insert(0, 'patcher', patcher)
        df_QC.insert(0, 'OP', OP[:-1])
        df_QC['Rs_change'], df_QC['Rin_change'], df_QC['Rs_start'] = Rs_diff, Rin_diff, Rs

        record_sorting = pd.concat([record_sorting.loc[:], df_QC]).reset_index(drop=True)
    record_sorting.to_excel(work_dir + 'data_tables/' + OP + '_QC_measures_minis.xlsx') 

def get_intrinsic_properties_df_no_VM_file (human_dir, OP, tissue_source, patcher, age, inj):
    '''
    saves a pandas dataframe with intrinsic cell peoperties for OP
    '''
    work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)

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
        work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)
        [print(key,':',value) for key, value in indices_dict.items()]
        # index_vc_in = [int(item) for item in input('Vc files corresponding to characterization files for ' + OP +' (input with spaces in between)').split()]
        # #saved the original indices
        # indices_dict['vc_orig'] = indices_dict['vc']
        # indices_dict['vc'] = index_vc_in
        # # index_char = [int(item) for item in input('corresponding characterization files for ' + OP +' (input with spaces in between)').split()]
        # # indices_dict['freq analyse'] = index_char
 
    active_chans_meta = sort.get_json_meta(human_dir, OP, patcher, '_meta_active_chans.json')[0]
    cortex_out_time = sort.get_datetime_from_input(active_chans_meta['OP_time'])

    #creating the dataframe
    df_OP = pd.DataFrame(columns=['filename', 'slice', 'cell_ch', 'cell_ID', 'day', 'treatment', 
    'hrs_incubation', 'repatch', 'hrs_after_OP', 'Rs', 'Rin', 'resting_potential', 'max_spikes', 'Rheobase', 
    'AP_heigth', 'TH', 'max_depol', 'max_repol', 'membra_time_constant_tau', 'capacitance'])

    for i in range(len(indices_dict['vc'])):
        vc = indices_dict['vc'][i]
        char = indices_dict['freq analyse'][i]
        slic = slice_names[vc]
        day = 'D1'
        if slic[-2:] == 'D2': 
            day = 'D2'
        treatment = active_chans_meta['treatment'][i][0]

        filename_vc = work_dir + filenames[vc]
        filename_char = work_dir + filenames[char]
        filename_con_screen = work_dir + filenames[indices_dict['con_screen'][i]]
        time_after_op = sort.get_time_after_OP(filename_char, cortex_out_time)

        active_channels = active_chans_meta['active_chans'][i]
        
        cell_IDs = hcf.get_cell_IDs(filename_char, slic, active_channels)
        Rs, Rin = hcf.get_access_resistance(filename_vc, active_channels) 
    
        charact_dict = hcf.load_traces(filename_char)
        inj = hcf.read_inj('full')
        tau_all, capacitance_all, mc_all, V65_all, RMPs = hcf.get_hyperpolar_param(charact_dict, active_channels, inj)
        
        params1_df = pd.DataFrame({'filename': filenames[char], 'slice' : slic, 'cell_ch': active_channels,
        'hrs_after_OP' : time_after_op,
        'cell_ID': cell_IDs, 'day' : day , 'treatment': treatment, 'Rs' : Rs, 'Rin': Rin, 'resting_potential': RMPs })        

        charact_params  = hcf.all_chracterization_params(filename_char, active_channels, inj)
        df_char = pd.DataFrame.from_dict(charact_params)

        df_to_add = pd.concat([params1_df, df_char], axis = 1)
        df_OP = pd.concat([df_OP.loc[:], df_to_add]).reset_index(drop=True)

        #plotting function
        plotting_funcs.plot_vc_holding (filename_vc, active_channels)
        plotting_funcs.plots_for_charact_file(filename_char, active_channels, inj)
        #plotting_funcs.plot_connect(filename_con_screen, active_channels)
        plotting_funcs.plot_connect_old_win(filename_con_screen, active_channels)
     
    tissue = pd.Series(tissue_source).repeat(len(df_OP))
    OPs = pd.Series(OP).repeat(len(df_OP))
    researcher = pd.Series(patcher).repeat(len(df_OP))
    patient_age = pd.Series(age).repeat(len(df_OP))
    series_df = pd.DataFrame({'tissue_source': tissue, 'OP': OPs, 'patcher': researcher, 'patient_age': patient_age}).reset_index(drop=True)

    df_intrinsic = pd.concat([series_df, df_OP], axis = 1)

    df_intrinsic.to_excel(work_dir + 'data_tables/' + OP + '_Intrinsic_and_synaptic_properties.xlsx') 
    #df_intrinsic.to_csv(work_dir + 'data_tables/' + OP + '_Intrinsic_and_synaptic_properties.csv')

    print('Intrinsic properties DataFrame for  ' + OP + 'saved successfully. ' + '\n' + 'Exclude recordings if necessary.')


def remove_bad_data (OP, patcher, human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)

    intr_df = pd.read_excel(work_dir + 'data_tables/' + OP + '_Intrinsic_and_synaptic_properties.xlsx')
    df_vc_qc =  pd.read_excel(work_dir + 'data_tables/' + OP + '_QC_measures_rs.xlsx')

    comments = list(range(0, len(intr_df)))
    big_change_in_Rs = df_vc_qc['cell_ID'][abs(df_vc_qc['change_rs']) > 15].tolist()
    for i in range(len(big_change_in_Rs)):
        cell_indx = intr_df.index[intr_df['cell_ID'] == big_change_in_Rs[i]][0]
        comments[cell_indx] = 'exclude; big change in Rs from start to end'

    #remove cells where Rs > 30
    Rs_QC_fail = df_vc_qc['cell_ID'][df_vc_qc['Rs_end'] > 30].tolist()
    for j in range(len(Rs_QC_fail)):
        cell_indx = intr_df.index[intr_df['cell_ID'] == Rs_QC_fail[j]][0]
        comments[cell_indx] = str(comments[cell_indx]) + '; exclude; Rs_end > 30'

    intr_df.insert(len(intr_df.columns), 'comments', comments)

    intr_df.to_excel(work_dir + 'data_tables/' + 'QC_passed_' + OP + '_Intrinsic_and_synaptic_properties.xlsx') 






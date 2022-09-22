import os
import human_characterisation_functions as hcf
import intrinsic_props_plotting_funcs as in_props_plot
import pandas as pd
import sorting_functions as sort

def get_intrinsic_properties_df(OP, tissue_source, patcher, age, inj):
    '''
    saves a pandas dataframe with intrinsic cell peoperties for OP
    '''
    OP_folder = OP + '/'
    
    if patcher == 'Verji': 
        work_dir = human_dir + 'data_verji/'+ OP_folder
    else:
        work_dir = human_dir + 'data_rosie/'+ OP_folder 
    
    file_list = sort.get_sorted_file_list(work_dir)
    df_rec= sort.get_lab_book(work_dir)
    filenames = sort.get_abf_files(file_list)
    slice_indx, def_slice_names, indices_dict = sort.sort_protocol_names (file_list, df_rec)

    #adjusting the slice names 
    slice_names = sort.fix_slice_names (def_slice_names, slice_indx)

    #creating a dir to save plots and data_tables (if not existing)
    dir_plots = sort.make_dir_if_not_existing (work_dir, 'plots')
    sort.make_dir_if_not_existing (work_dir, 'data_tables')

    #check if the traces dir is empty and only then plot the mimiddle sweep for each filename
    traces_folder =  os.path.join(dir_plots, "traces/")
    if os.path.isdir(traces_folder) == 0 :
        for rec in range(len(filenames)):
            filename = work_dir + filenames[rec]
            in_props_plot.plot_middle_sweep(filename)
    else:
         print("skipping plotting")

    #QC indices
    [print(key,':',value) for key, value in indices_dict.items()]
    if len(indices_dict['freq analyse']) != len(indices_dict['vc']): 
        print('Fix protocol names. Unequal number of VC and freq analyse protocols')
        index_vc_in = [int(item) for item in input('Vc files corresponding to characterization files for ' + OP +' (input with spaces in between)').split()]
        #saved the original indices
        indices_dict['vc_orig'] = indices_dict['vc']
        indices_dict['vc'] = index_vc_in
        # index_char = [int(item) for item in input('corresponding characterization files for ' + OP +' (input with spaces in between)').split()]
        # indices_dict['freq analyse'] = index_char
 
    proceed_y_n = input("do all traces correspond to specified filenames in lab book for " + OP +  "(y/n)?")
    if proceed_y_n == 'n': 
        print('correct the lab book entries. Continuing to next OP')
        pass

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

        active_channels = [int(item) for item in input('Channels used in ' + filenames[vc] +'(input with spaces in between)').split()]
        
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
        in_props_plot.plot_vc_holding (filename_vc, active_channels)
        in_props_plot.plots_for_charact_file(filename_char, active_channels, inj)
     
    tissue = pd.Series(tissue_source).repeat(len(df_OP))
    OPs = pd.Series(OP).repeat(len(df_OP))
    researcher = pd.Series(patcher).repeat(len(df_OP))
    patient_age = pd.Series(age).repeat(len(df_OP))
    series_df = pd.DataFrame({'tissue_source': tissue, 'OP': OPs, 'patcher': researcher, 'patient_age': patient_age}).reset_index(drop=True)

    df_intrinsic = pd.concat([series_df, df_OP], axis = 1)

    df_intrinsic.to_excel(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.xlsx') 
    df_intrinsic.to_csv(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.csv')

    print('Intrinsic properties DataFrame for  ' + OP + 'saved successfully. ' + '\n' + 'Exclude recordings if necessary.')

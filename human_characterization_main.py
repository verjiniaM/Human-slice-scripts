
#%%
import os
import human_characterisation_functions as hcf
import intrinsic_props_plotting_funcs as in_props_plot
import pandas as pd
import glob
import sorting_functions as sort
import itertools as it
import funcs_for_results_tables as get_results

#loading the updated experiments_overview and the old summary_data_table
human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
results_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'

exp_view = pd.read_excel(human_dir + 'experiemnts_overview.xlsx') 
#loading the latest data table, should be the last one in the folder
latest_sum_data_table = sorted(glob.glob(results_dir + 'summary_data_tables/' + '*.xlsx'))[-1]
end_path = latest_sum_data_table.rfind('/')
summary_data_table = pd.read_excel(latest_sum_data_table)
print('Using summary data table from ' + latest_sum_data_table[end_path + 1 : end_path + 11])

#%%
# updating full OP list from latest summmary based on all OPs in exp_overview
last_update_op_list = summary_data_table.OP.unique().tolist()
newest_op_list = exp_view.OP.unique().tolist()
all_OPs = last_update_op_list + newest_op_list
op_to_analyse = [i for i in all_OPs if all_OPs.count(i)==1]


#%%
#Decide which OPs need to be analyzed 


tissue_source = "Bielefeld"
patcher = "Verji"
age = "A"
OP = op_to_analyse[-1]
    

# for i in op_to_analyse:
#     OP = op_to_analyse[i]
#     y_or_n = input('Do you want to start analysis of ' + OP + '(y/n)?') 
#     if y_or_n == "n":
#         continue

#     tissue_source = input('Tissue source for ' + OP + '(Bielefeld, Mitte, Virchow): ')
#     patcher = input('Patcher ' + OP + '(Verji or Rosie): ')
#     age = input('Patient age ' + OP + '(if not known: A for adult, J for juvenile): ')
#     inj = "full"
    
#     print('starting analysis for '+ OP)
#     get_results.get_intrinsic_properties_df(OP, tissue_source, patcher, age, inj)

    for file in range(len(file_list)):
        #pclamp files
        if  (file_list[file][-5:] == '.xlsx' and file_list[file][:2] != '~$'): 
            df_rec = pd.read_excel(work_dir + file_list[file], header = 1)
            index_vc = df_rec.index[df_rec['protocol'] == 'vc'].tolist()
            index_vc_end = df_rec.index[df_rec['protocol'] == 'vc_end'].tolist()
            index_spontan = df_rec.index[df_rec['protocol'] == 'vm'].tolist()
            index_min = df_rec.index[df_rec['protocol'] == 'minis'].tolist()

    print('VC files:     ' + ' '.join(str(h) for h in index_vc))
    print('VC_end files: ' + ' '.join(str(h) for h in index_vc_end))
    print('Spontaneous files: ' + str(index_spontan))
    print('Mini files: ' + str(index_mini))

    proceed = input('Proceeding with QC steps? All mvc, vc_end, mini, and spontaneous files detected? (y/n)')

    if proceed == 'n': 
        
        input('waiting to continue. Do not forget to close the excel file, if you change something')
        print('updating the excel file')

        print('Fix protocol names. Unequal number of VC and freq analyse protocols')
        index_vc = [int(item) for item in input('Vc files with appropriate vc_end files for ' + OP +' (input with spaces in between)').split()]
        index_vc_end = [int(item) for item in input('Vc_end files corresponding to characterization files for ' + OP +' (input with spaces in between)').split()]

        print('VC files:     ' + ' '.join(str(h) for h in index_vc))
        print('VC_end files: ' + ' '.join(str(h) for h in index_vc_end))
        print('Spontaneous files: ' + str(index_spontan))
        print('Mini files: ' + str(index_mini))

        if len(index_vc) != len(index_vc_end): 
            print('Fix protocol names. Unequal number of VC and vc_end files')
            index_vc = [int(item) for item in input('Vc files corresponding to vc_end files for ' + OP +' (input with spaces in between)').split()]
            index_vc_end = [int(item) for item in input('Vc_end files ' + OP +' (input with spaces in between)').split()]

    #date_frame for quality control - change in Rin, Rs 
    df_qc = pd.DataFrame(columns=['OP','patcher', 'filename', 'slice', 'cell_ch', 'Rs_start', 'Rin_start', 'Rs_end', 'Rin_end', 
    'chage_rs', 'change_rin'])

    for k in range(len(index_vc_end)):
        vc = index_vc[k]
        vc_end = index_vc_end[k]
        slice = slice_names[vc]

        filename_vc = work_dir+filenames[vc]
        filename_vc_end = work_dir+filenames[vc_end]

        active_channels = [int(item) for item in input('Channels used in ' + filenames[vc] +'(input with spaces in between)').split()]

        for l in range(len(active_channels)):
            ch = active_channels[l]
            Rs, Rin = hcf.get_access_resistance(filename_vc, ch) 
            Rs_end, Rin_end = hcf.access_resistance(filename_vc_end, ch) 
            cellID = filenames[vc][:-7]+slice+'c'+str(ch)

            data_to_add = pd.DataFrame({'OP':OP[:-1], 'patcher':patcher, 'filename':filenames[vc],'slice':slice, 
            'cell_ch':ch, 'Rs_start': Rs, 'Rin_start': Rin, 'Rs_end': Rs_end, 'Rin_end': Rin_end, 
            'chage_rs': Rs-Rs_end, 'change_rin':Rin-Rin_end}, index=[0])
            df_qc = pd.concat([df_qc.loc[:], data_to_add]).reset_index(drop=True)

    df_qc.to_excel(work_dir + 'data_tables/' + OP[:-1] + '_QC_measures_rs.xlsx') 
    df_qc.to_csv(work_dir + 'data_tables/' + OP[:-1] + '_QC_measures_rs.csv')

    #quality check for minis
    #calcualted the signal drift and sorts suitable/unsuitable recordings

    record_sorting = pd.DataFrame(columns = ['OP','patcher', 'filename', 'mini/spontan?', 'cell_ch','swps_to_analyse', 'swps_to_discard', 
    'min_vals_each_swp', 'max_vals_each_swp', 'drift'])

    for u in range(len(index_spontan)):
        spontan = index_spontan[u]
        slice = slice_names[spontan]

        filename_spontan = work_dir+filenames[spontan]
        
        active_channels = [int(item) for item in input('Channels used in ' + filenames[spontan] +'(input with spaces in between)').split()]

        for q in range(len(active_channels)):
            ch = active_channels[q]
            min_vals, max_vals, good_swps, bad_swps, drift = hcf.rec_stability (filename_spontan, ch , 50)

            data_to_add = pd. DataFrame({'OP':OP[:-1], 'patcher':patcher, 'filename':filenames[spontan], 'mini/spontan?':'spontan',
            'slice':slice, 'cell_ch':ch, 'swps_to_analyse':[good_swps], 'swps_to_discard':[bad_swps], 'min_vals_each_swp':[min_vals], 
            'max_vals_each_swp':[max_vals], 'drift':[drift]})
            record_sorting = pd.concat([record_sorting.loc[:], data_to_add]).reset_index(drop=True)



    for h in range(len(index_mini)):
        mini = index_mini[h]
        slice = slice_names[mini]

        mini_file = work_dir+filenames[mini]

        active_channels = [int(item) for item in input('Channels used in ' + filenames[mini] +'(input with spaces in between)').split()]

        for b in range(len(active_channels)):
            ch = active_channels[b]
            min_vals, max_vals, good_swps, bad_swps, drift = hcf.rec_stability (mini_file, ch , 50)

            data_to_add = pd. DataFrame({'OP':OP[:-1], 'patcher':patcher, 'filename':filenames[mini], 'mini/spontan?':'minis',
            'slice':slice, 'cell_ch':ch, 'swps_to_analyse':[good_swps], 'swps_to_discard':[bad_swps], 'min_vals_each_swp':[min_vals], 
            'max_vals_each_swp':[max_vals], 'drift':[drift]})
            record_sorting = pd.concat([record_sorting.loc[:], data_to_add]).reset_index(drop=True)

    record_sorting.to_excel(work_dir + 'data_tables/' + OP + '_QC_measures_events.xlsx') 
    record_sorting.to_csv(work_dir + 'data_tables/' + OP + '_QC_measures_events.csv')

    print('The data for ' + OP + 'is fully analyzed and saved in ' + work_dir + 'data_tables/')

# %%


#%%
import os
import human_characterisation_functions as hcf
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
#     get_results.get_intrinsic_properties_df(human_dir, OP, tissue_source, patcher, age, inj)
#     get_results.get_QC_access_resistance_df (human_dir, OP, patcher)
#     get_results.get_con_params_df(human_dir, OP, patcher)

    #quality check for minis
    #calcualted the signal drift and sorts suitable/unsuitable recordings
def get_minis_QC():
    work_dir, filenames, indices_dict, slice_names = get_OP_metadata(human_dir, OP, patcher)
    

    record_sorting = pd.DataFrame(columns = ['OP','patcher', 'filename', 'mini/spontan?', 'cell_ch','swps_to_analyse', 'swps_to_discard', 
    'min_vals_each_swp', 'max_vals_each_swp', 'drift'])


    for u in range(len(index_spontan)):
        filename_spontan = index_spontan[u]
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

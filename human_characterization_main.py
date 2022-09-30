
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

# updating full OP list from latest summmary based on all OPs in exp_overview
last_update_op_list = summary_data_table.OP.unique().tolist()
newest_op_list = exp_view.OP.unique().tolist()
all_OPs = last_update_op_list + newest_op_list
op_to_analyse = [i for i in all_OPs if all_OPs.count(i)==1]

for i in op_to_analyse:
    OP = op_to_analyse[i]
    y_or_n = input('Do you want to start analysis of ' + OP + '(y/n)?') 
    if y_or_n == "n":
        continue

    tissue_source = input('Tissue source for ' + OP + '(Bielefeld, Mitte, Virchow): ')
    patcher = input('Patcher ' + OP + '(Verji or Rosie): ')
    age = input('Patient age ' + OP + '(if not known: A for adult, J for juvenile): ')
    inj = "full"
    
    print('starting analysis for '+ OP)
    get_results.get_intrinsic_properties_df(human_dir, OP, tissue_source, patcher, age, inj)
    get_results.get_QC_access_resistance_df (human_dir, OP, patcher)
    get_results.get_con_params_df(human_dir, OP, patcher)
    get_results.get_spontan_QC(human_dir, OP, patcher)
    get_results.get_mini_QC(human_dir, OP, patcher)

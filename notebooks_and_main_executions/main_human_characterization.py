import os
import pandas as pd
import glob
import intrinsic_props_and_connectivity.funcs_sorting as sort
import src.intrinsic_props_and_connectivity.funcs_for_results_tables as get_results


# #%%
OP = 'OP231130'
patcher = 'Verji'
tissue_source = 'Bielefeld'
inj = 'full'
age = 39
#%%
#loading the updated experiments_overview and the old summary_data_table
human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
results_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/'

exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0]) 
exp_view_new = sort.update_op_list(human_dir, exp_view)
sort.add_cortex_out_time(human_dir, exp_view_new)

#%%
#loading the latest data table, should be the last one in the folder
latest_sum_data_table = sorted(glob.glob(results_dir + 'summary_data_tables/intrinsic_properties/' + '*.xlsx'))[-1]
end_path = latest_sum_data_table.rfind('/')
summary_data_table = pd.read_excel(latest_sum_data_table)
print('Using summary data table from ' + latest_sum_data_table[end_path + 1 : end_path + 11])

# updating full OP list from latest summmary based on all OPs in exp_overview
last_update_op_list = summary_data_table.OP.unique().tolist()
newest_op_list = exp_view.OP.unique().tolist()
all_OPs = last_update_op_list + newest_op_list
op_to_analyse = [i for i in all_OPs if all_OPs.count(i)==1]

for OP in op_to_analyse:
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
    input() # to give the user time to check and correct the results
    get_results.get_con_params_df(human_dir, OP, patcher)
    get_results.get_con_screen_VC(human_dir, OP, patcher)
    get_results.get_spontan_QC(human_dir, OP, patcher)
    get_results.get_minis_QC(human_dir, OP, patcher)
    get_results.check_cell_IDs(human_dir, OP, patcher) #this could be run multiple times
    input('fix cell IDs in spontan df if needed')
    get_results.check_cell_IDs(human_dir, OP, patcher)
    get_results.create_IFF_data_table(OP, patcher, '_meta_active_chans.json', inj, human_dir)


#%%
#run after QC checked df has been manually compared to the Onsets and the AP_prop plots
#checks if somewhhre the capacitance was not calculated properly 
#or if some reaptched cells are wrongly labeled

work_dir = sort.get_work_dir(human_dir, OP, patcher)
results_df = pd.read_excel(glob.glob(work_dir + '/data_tables/' + '*QC_passed' + '*.xlsx')[0])

mask = (results_df['capacitance'] > 900) | (results_df['capacitance'] < 10) | (results_df['capacitance'].isnull())
cap_mistakes = results_df.loc[mask, :].reset_index()
if len(cap_mistakes) > 0:
    for i in range(len(cap_mistakes)):
        print('capacitance manual for '+ cap_mistakes.filename[i] + ' Ch' +str(cap_mistakes.cell_ch[i]))

repatch_df = results_df[results_df['repatch'] == 'yes']
repatch_df.reset_index(inplace = True, drop = True)

not_repatched_cells = []
for cell in repatch_df['cell_ID'].unique():
    if len(repatch_df[repatch_df['cell_ID'] == cell]) < 2:
        not_repatched_cells.append(cell)
if len(not_repatched_cells) > 0:
    print('Not repatched cells ' + str(not_repatched_cells)[1:-1])

#%%

# get_results.prapare_for_event_analysis(human_dir)
# connections_IC, connections_VC = get_results.collect_connections_df(human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/')

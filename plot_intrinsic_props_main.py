import plot_intrinsic_props as pl_intr
import funcs_for_results_tables as get_results
import pandas as pd

df_intr_props = get_results.collect_intrinsic_df()

adult_df = pl_intr.get_adult_data(df_intr_props) #only patient age >18 and other QC criteria 
adult_df = adult_df[adult_df['area'] == 'temporal']
adult_df = pl_intr.create_new_cell_IDs(adult_df)
adult_df_8mM = adult_df[(adult_df['high K concentration'] != '15 mM')]
pl_intr.plot_param_for_days(adult_df_8mM, 'all', '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/intrinsic_properties/all_QC_passed/') 

repatch_df = pl_intr.get_repatch_df(adult_df) #adds also the new cell_IDs
repatch_df = repatch_df.sort_values(['cell_ID_new', 'treatment'])
pl_intr.plot_param_for_days(repatch_df, '8mM_repatch') #for all params repatch!

#%%
#for connectivity summary plot

all_cons_IC, all_cons_VC =  get_results.collect_connections_df()

exclude_not_adult = pl_intr.in_list1_not_in_list2(all_cons_IC['OP'].unique(), adult_df['OP'].unique())
    
for OP in exclude_not_adult:
    all_cons_IC = all_cons_IC.drop(all_cons_IC.index[all_cons_IC['OP'] == OP]) 
    all_cons_VC = all_cons_IC.drop(all_cons_IC.index[all_cons_IC['OP'] == OP])         
all_cons_IC.reset_index(inplace = True, drop = True)
all_cons_VC.reset_index(inplace = True, drop = True)

connect_QC_passed = pl_intr.get_QC_connectivity_df(all_cons_IC)
pl_intr.plot_connect_amplitude(connect_QC_passed, 'all_IC_adult')
repatch_connect = pl_intr.get_repatch_connectivity_df(connect_QC_passed)
pl_intr.plot_connect_amplitude(repatch_connect, 'repatch_IC_adult')

connect_QC_passed_VC = pl_intr.get_QC_connectivity_VC(all_cons_VC)
pl_intr.plot_connect_amplitude(connect_QC_passed_VC, 'all_VC_adult')
# repatch_connect_VC = pl_intr.get_repatch_connectivity_df(connect_QC_passed_VC)
# pl_intr.plot_connect_amplitude(repatch_connect_VC, 'repatch_VC')

#%%
# Plot Initial firing frequency and number of APs against injected current

IFF_collected = get_results.collect_IFF_dfs(human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/')
#adult_df = adult_df[(adult_df['high K concentration'] != '15 mM')]

IFF_adult, IFF_repatch = get_results.get_QC_for_IFF_df(adult_df, IFF_collected) #takes only cell_IDs present in 
IFF_repatch_firing = get_results.remove_non_firing_cells_D1(IFF_repatch)

#pl_intr.plot_IFF_distribution(IFF_repatch, 'repatch', 'num_aps') #data_type - all, repatch, repatch firing cells; DV - num_aps or IFF
#pl_intr.plot_IFF_distribution(IFF_repatch, 'repatch', 'IFF')

pl_intr.plot_IFF_avg_against_current_inj(IFF_repatch, 'repatch', 'num_aps')
pl_intr.plot_IFF_avg_against_current_inj(IFF_repatch, 'repatch', 'IFF')


IFF_collected = get_results.collect_IFF_dfs(human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/')

adult_df_8mM = adult_df[(adult_df['high K concentration'] != '15 mM')]
IFF_adult, IFF_repatch_8mM = get_results.get_QC_for_IFF_df(adult_df_8mM, IFF_collected) #takes only cell_IDs present in 

pl_intr.plot_IFF_avg_against_current_inj(IFF_repatch_8mM, 'repatch_8mM', 'num_aps')
pl_intr.plot_IFF_avg_against_current_inj(IFF_repatch_8mM, 'repatch_8mM', 'IFF')

# pl_intr.plot_IFF_distribution(IFF_repatch_firing, 'repatch_firing_D1_8mM', 'num_aps') #data_type - all, repatch, repatch firing cells; DV - num_aps or IFF
# pl_intr.plot_IFF_distribution(IFF_repatch_firing, 'repatch_firing_D1_8mM', 'IFF')

# pl_intr.plot_IFF_avg_against_current_inj(IFF_repatch_firing, 'repatch_firing_D1_8mM', 'num_aps')
# pl_intr.plot_IFF_avg_against_current_inj(IFF_repatch_firing, 'repatch_firing_D1_8mM', 'IFF')

#%%
#analysis high K 

df_resting = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP221116/data_tables/OP221116_RMPs_over_time.xlsx')
save_dir_RMP = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/high K evaluation/'

df_resting = pl_intr.format_RMP_column(df_resting)
df_resting = pl_intr.get_hh_mm_from_rec_time(df_resting)
pl_intr,plot_RMP_time(df_resting, save_dir_RMP)


df_ap_props = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP221116/data_tables/OP221116_changes_in_AP_props.xlsx')
df_ap_props = pl_intr.get_hh_mm_from_rec_time(df_ap_props)

save_dir_ap_props = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/high K evaluation/'
ap_props_dict = pl_intr.dict_for_plotting()
pl_intr.plot_ap_props(df_ap_props, ap_props_dict, save_dir_ap_props) 


#plotting whole RMP changes of a Ch over recordings 
file_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP221116/'
plots_destination = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP221116/plots/'
files = ['22n16037.abf', '22n16039.abf', '22n16040.abf']
channel = 6

pl_intr.plot_full_RMP_trace(file_folder, files, plots_destination, channel) #files = ['fn1', 'fn2']

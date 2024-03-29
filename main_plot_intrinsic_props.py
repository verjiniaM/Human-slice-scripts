import funcs_for_results_tables as get_results
import pandas as pd
import glob
import funcs_plot_intrinsic_props as pl_intr

human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0]) 

#collect all QC intrinsic datatables
df_intr_props = get_results.collect_intrinsic_df()
#add patient info
# df_op_info = pd.DataFrame(columns = ['OP', 'patient_age', 'resected area'])
# for OP in df_intr_props['OP'].unique():
#     OP = df_intr_props['OP'].loc[df_intr_props['OP'] == OP].to_list()[0]
#     area = exp_overview['region'].loc[exp_overview['OP'] == OP].to_list()[0]
#     patient_age = df_intr_props['patient_age'].loc[df_intr_props['OP'] == OP].tolist()[0]
#     df_op_add = pd.DataFrame({'OP':OP, 'patient_age':[patient_age], 'resected area':area})
#     df_op_info = pd.concat([df_op_info.loc[:], df_op_add]).reset_index(drop=True)
  
#filter on age and hrs incubation
adult_df = pl_intr.filter_adult_hrs_incubation_data(df_intr_props, min_age = 18, hrs_inc = 16, max_age = 151) # QC criteria 
adult_df = adult_df[adult_df['area'] == 'temporal']
op_color_dict = pl_intr.get_op_color_dict(adult_df)
adult_df = pl_intr.create_new_cell_IDs(adult_df)
#adult_df = pl_intr.get_precise_treatment(adult_df)

#create slice comparison and repatch dataframes
adult_df_slice_comparison = adult_df.loc[adult_df['repatch'] == 'no']
adult_df_repatch = pl_intr.get_repatch_df(adult_df)
adult_df_repatch = adult_df_repatch.sort_values(['cell_ID_new', 'treatment'])
adult_repatch_norm = pl_intr.get_normalized_df(adult_df_repatch)

#plot all parameters 
pl_intr.plot_param_for_days_slice(adult_df_slice_comparison, op_color_dict)
pl_intr.plot_param_for_days_slice_TTX_incl(adult_df_slice_comparison, op_color_dict)
pl_intr.plot_param_for_days_repatch(adult_df_repatch, op_color_dict)
pl_intr.plot_param_for_days_repatch_plus_TTX(adult_df_repatch, op_color_dict)
pl_intr.plot_param_for_days_repatch_all_params(adult_df_repatch, op_color_dict)
#pl_intr.plot_param_for_days_repatch_norm(adult_repatch_norm, op_color_dict)

###########################################################################################
#%%
###### Plot 3 h incubation data
results_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/intrinsic_properties/'

adult_df_short = pl_intr.filter_adult_hrs_incubation_data(df_intr_props, min_age = 18, hrs_inc = 0, max_age = 151, max_hrs_incubation = 5) # QC criteria 
#adult_df_short = adult_df_short[adult_df_short['area'] == 'temporal']
op_color_dict = pl_intr.get_op_color_dict(adult_df)
adult_df_short = pl_intr.create_new_cell_IDs(adult_df_short)
#adult_df = pl_intr.get_precise_treatment(adult_df)

#create slice comparison and repatch dataframes
adult_df_slice_comparison_short = adult_df_short.loc[adult_df_short['repatch'] == 'no']
adult_df_repatch_short = pl_intr.get_repatch_df(adult_df_short)
adult_df_repatch_short = adult_df_repatch_short.sort_values(['cell_ID_new', 'treatment'])
adult_repatch_norm_short = pl_intr.get_normalized_df(adult_df_repatch_short)

#plot all parameters 
pl_intr.plot_param_for_days_slice(adult_df_slice_comparison_short, op_color_dict, results_dir + 'adult_slice/3h_inc/')
pl_intr.plot_param_for_days_repatch(adult_df_repatch_short, op_color_dict,results_dir + 'adult_repatch/3h_inc/')
pl_intr.plot_param_for_days_repatch_all_params(adult_df_repatch_short, op_color_dict, results_dir + 'adult_repatch/3h_inc/')
pl_intr.plot_param_for_days_repatch_norm(adult_repatch_norm_short, op_color_dict,results_dir + 'adult_repatch/3h_inc/')

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

IFF_QC, IFF_repatch = pl_intr.get_QC_for_IFF_df(adult_df, IFF_collected) 

IFF_QC = pl_intr.remove_non_firing_cells_D1(IFF_QC)
IFF_QC = IFF_QC.loc[IFF_QC['treatment'] != 'TTX']
IFF_QC_slice = IFF_QC.loc[IFF_QC['repatch'] == 'no']
IFF_QC_slice_firing = pl_intr.remove_non_firing_cells_D1(IFF_QC_slice)
IFF_QC_slice_firing = get_precise_treatment(IFF_QC_slice_firing)

pl_intr.plot_IFF_avg_against_current_inj(IFF_QC_slice_firing, 'slice_firing', 'IFF')


IFF_repatch_firing = pl_intr.remove_non_firing_cells_D1(IFF_repatch)
IFF_repatch_firing = IFF_repatch_firing.loc[IFF_repatch_firing['treatment'] != 'TTX']
IFF_repatch_firing = get_precise_treatment(IFF_repatch_firing)

pl_intr.plot_IFF_avg_against_current_inj(IFF_repatch_firing, 'repatch_firing', 'IFF')

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

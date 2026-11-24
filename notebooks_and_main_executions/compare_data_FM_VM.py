import numpy as np
import pandas as pd
import ephys_analysis.funcs_plot_conn_paper_figures as con_plots
import ephys_analysis.funcs_con_paper_analysis_as_FM as funcs_FM
import connectivity_estimation as con_estm

work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/'
save_plots_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/connectivity_paper/plots/V2_'

#%%
# ANALYZE DATA AS FM
VM_data_dir = f'{work_dir}results/human/connectivity_paper/data/'

# GET THE DATA FRAME IF LOST
# intr_all = pd.read_csv(f'{VM_data_dir}con_analysis_clean_all/all_intr_props_all_treatments.csv')
# intr_all_VM = funcs_FM.get_data_analyzed_as_FM(intr_all)

# read df with fixed tau
intr_all_VM = pd.read_csv(f'{VM_data_dir}all_intr_as_FM.csv')

df_con_VM = pd.read_csv(f'{VM_data_dir}all_connect_cross_with_intr_with_etypes.csv')

df_FM = pd.read_csv(f'{work_dir}data/human/humandata/data_tables/tcell_all_psdn.csv')
df_FM = df_FM[(df_FM.removal_boolean == 1) & \
            (df_FM.is_pc == 1) & \
                (df_FM.piadistance < 1200) & \
                    (df_FM.type == 1)]
# len(df_FM) # should be 901


#%%
# dictionary for comparison
dict_intrinsic_VM_FM = funcs_FM.get_dict_VM_FM()
df_intr_VM_f = funcs_FM.apply_filter_properties_mask(intr_all_VM)

# df_intr_VM_f = funcs_FM.remove_nans_from_dict(df_intr_VM_f)

# # t-test and visualize params across dataframes
# con_plots.intr_params_VM_vs_FM(df_FM, df_intr_VM_f, dict_intrinsic_VM_FM, save_plots_dir)
# con_plots.t_test_and_plot_params(df_FM, df_intr_VM_f, dict_intrinsic_VM_FM, save_plots_dir)

# # create KNN and plot results only with complete data
# paper_pca, my_pca, df_VM_e_types = funcs_FM.knn_analysis(dict_intrinsic_VM_FM, df_intr_VM_f, df_FM)
# df_paper_VM_complete_0 = funcs_FM.create_df_plot(paper_pca, df_FM, my_pca, df_intr_VM_f)
# con_plots.plot_e_types_across_datasets(df_paper_VM_complete_0)



# USE IMPUTATION COMPLETE DATASET

df_paper_VM_complete, df_VM_etypes = funcs_FM.impute_on_own_data_assign_e_types(
    df_intr_VM_f, df_FM, dict_intrinsic_VM_FM)

con_plots.plot_e_types_across_datasets(df_paper_VM_complete)

# check percentages
etype_distribution = funcs_FM.compare_etype_distributions(
    df_FM, df_paper_VM_complete[df_paper_VM_complete.dataset == 'my data'])

etype_distribution

#%%
# ADD ETYPE TO CONNECTIVITY DATA

df_con_VM = df_con_VM[~(df_con_VM.treatment == 'fixed')].copy()

cell_ID_lookup = {(row['cell_IDs_match']): row['e_type']
     for _, row in df_VM_etypes.iterrows()}

pre_etypes = [cell_ID_lookup.get((cell_id), np.nan)
    for cell_id in df_con_VM['cell_IDs_match_pre']]

post_etypes = [cell_ID_lookup.get((cell_id), np.nan)
    for cell_id in df_con_VM['cell_IDs_match_post']]

df_con_VM['pre_etype'] = pre_etypes
df_con_VM['post_etype'] = post_etypes

# remove rows with missing e types
df_con_VM_etypes = df_con_VM.dropna(subset=['pre_etype', 'post_etype']).copy()

#%%

# plot
df_slice_connections = funcs_FM.get_slice_connection_summary(df_VM_etypes, df_con_VM_etypes)

df_slice_connections = df_slice_connections[df_slice_connections.treatment != 'fixed']
con_estm.quick_dirty_plot_connectivity(df_slice_connections, 'none', save_dir = False)


con_percent_saved = pd.read_csv(f'{VM_data_dir}/con_analysis_clean_all/con_percentages_from_cross_checked_with_intr.csv')
con_percent_saved = con_percent_saved[con_percent_saved.treatment != 'fixed']
con_estm.quick_dirty_plot_connectivity(con_percent_saved, 'none', save_dir = False)


#%%

# ESTIMATION OF CONNECTION TYPES


import importlib
importlib.reload(funcs_FM)
df_all_possible_cons = funcs_FM.get_all_possible_connection_from_intr(df_VM_etypes, df_con_VM_etypes)

df_D1 = df_all_possible_cons[df_all_possible_cons['day'] == 'D1'].copy()
df_CTR_D2 = df_all_possible_cons[(df_all_possible_cons['day'] == 'D2') &
                             (df_all_possible_cons['treatment'] == 'Ctrl')].copy()
df_HiK_D2 = df_all_possible_cons[(df_all_possible_cons['day'] == 'D2') &
                             (df_all_possible_cons['treatment'] == 'high K')].copy()
df_TTX_D2 = df_all_possible_cons[(df_all_possible_cons['day'] == 'D2') &
                             (df_all_possible_cons['treatment'] == 'TTX')].copy()

con_matr_d1 = funcs_FM.get_connection_matrix(df_D1, df_FM)
con_matr_CTR = funcs_FM.get_connection_matrix(df_CTR_D2, df_FM)
con_matr_HiK = funcs_FM.get_connection_matrix(df_HiK_D2, df_FM)
con_matr_TTX = funcs_FM.get_connection_matrix(df_TTX_D2, df_FM)

print(f'All RS1 no incubation connection types \n {con_matr_d1} \n')
print(f'All CTR incubation    connection types \n {con_matr_CTR}\n')
print(f'All HiK incubation    connection types \n {con_matr_HiK}\n')
print(f'All TTX incubation    connection types \n {con_matr_TTX}\n')

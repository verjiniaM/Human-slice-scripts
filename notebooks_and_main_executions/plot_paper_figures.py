import importlib
import pandas as pd
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
import ephys_analysis.funcs_plot_paper_figures as paper_plot

#%%
#start from fig 2 --> decide size of text, itcks,labels, ratios
data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
     'paper_figs_collected_checked/data/'
# destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
#         'paper_figs_collected_checked/figures/trying_out/'

# # fig 1
# df_mea = pd.read_excel(data_dir + 'summary_pMEA_r.xlsx')
# mea_fig1(df_mea, destination_dir)

# fig 2
# df_slice = pd.read_excel(data_dir + 'slice_data_temporal.xlsx')

# paper_plot.firing_props_fig2_slice(df_slice, 'firing_cells_slice', 'num_aps', 'CI', destination_dir)
# paper_plot.firing_props_fig2_slice(df_slice, 'firing_cells_slice', 'IFF', 'CI', destination_dir)

# paper_plot.firing_props_fig3(df_slice, 'CI', destination_dir)


# # cell_IDs_slice = {'pre': '24201S1c8',
# #             'post CTR': '24117S2_D2c1',
# #             'post HiK': '23420S2_D2c7'}
# # plot_example_firing('slice', df_slice, destination_dir, cell_IDs_slice)

#paper_plot.plot_op_slice_cell_histogram(df_slice, 'slice', destination_dir)

# paper_plot.intr_params_fig2(df_slice, 'AP', destination_dir)

# get nnew version of no-firing cells
# slice_non_fire = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/trying_out/slice_non_firing_cells_count_V2.xlsx')
# paper_plot.plot_non_firing_across_inj_count(slice_non_fire,  'prop', 'slice', destination_dir )

# fig 3
# df_repatch = pd.read_excel(data_dir + 'repatch_data_temporal.xlsx')

# cell_IDs_repatch = {'CTR': 'vm24321S2c8',
#             'HiK': 'vm24926S4c1'}
# plot_example_firing_fig3('repatch', df_repatch, destination_dir, cell_IDs_repatch)
# paper_plot.firing_props_fig3(adult_df_repatch_TTX, 'SE', destination_dir+ 'TTX')

#paper_plot.firing_props_fig3(df_repatch, 'SE', destination_dir)
# paper_plot.firing_props_fig3(df_repatch, 'CI', destination_dir)

# Usage with your data

# paper_plot.plot_op_slice_cell_histogram(df_repatch, 'repatch', destination_dir)

# paper_plot.intr_params_fig3(df_repatch, destination_dir)

# repatch_non_fire = pl_intr.count_non_firing_cells_d2(df_repatch, destination_dir, 'repatch')
# paper_plot.plot_non_firing_across_inj_count(repatch_non_fire, 'prop', 'cell_ID_new', destination_dir )

# slice CTR all vs hrs_after_OP
# intr_slice_CTR(destination_dir, er_type = 'SE')
# intr_slice_CTR(destination_dir, er_type = 'CI') # doesn't complete the figure

# df_slice_all = pd.read_excel(data_dir + 'slice_all.xlsx')

# paper_plot.plot_example_firing_corr('slice', df_slice_all, destination_dir,
#                                     {'pre': '24215S1c4',
#                                      'Ctrl': '25213S5c4',
#                                      'CtrlD2': '2311S3_D2c8'})

# paper_plot.slice_all_CTR_hist(df_slice_all, destination_dir) # removes high K D2


# fig 5
# df_AIS = pd.read_excel(data_dir + 'AIS_all_data.xlsx')
# AIS_fig5(df_AIS, destination_dir)

# # incubation only
df_slice_short_inc = pd.read_excel(data_dir + '10.07.25_slice_incubation_only_no_move.xlsx')
df_slice_short_inc = df_slice_short_inc[df_slice_short_inc['hrs_after_OP'] < 35]
df_slice_short_inc = df_slice_short_inc[df_slice_short_inc['hrs_incubation'] > 15.9]
df_slice_short_inc.reset_index(drop = True, inplace = True)
df_slice_short_inc = df_slice_short_inc.loc[[i for i, sl in enumerate(df_slice_short_inc.slice) if len(sl) <= 2], :]

# paper_plot.plot_op_slice_cell_histogram(df_slice_short_inc, 'incubation', destination_dir)
# paper_plot.intr_params_inc_only(df_slice_short_inc, destination_dir + 'inc_ap_2')

# inc_non_fire = pl_intr.count_non_firing_cells_d2(df_slice_short_inc, destination_dir, 'inc_only')
# inc_non_fire = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/trying_out/.xlsx')
# paper_plot.plot_non_firing_across_inj_count(inc_non_fire,  'prop', 'inc_only', destination_dir )

import importlib
importlib.reload(paper_plot)
paper_plot.firing_props_inc_only_V2(df_slice_short_inc, 'CI', destination_dir + 'CI_')
# paper_plot.firing_props_inc_only_V2(df_slice_short_inc, 'SE', destination_dir + 'new_')

# cell_IDs_inc = {'Ctrl': '25220S3c3',
#             'high K': '25220S2c3'}
# plot_example_firing('inc_only', df_slice_short_inc, destination_dir, cell_IDs_inc)

# FIG - correlations
# data_dir =  '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
# save_dir_figs = data_dir[0:89] +  '/figures/draft5_figs/parts/'
# x = 'hrs_after_OP'

# df_ctrl_0_inc = pd.read_csv(data_dir + 'df_ctrl_0_inc.csv')
# df_ctrl_short_inc = pd.read_csv(data_dir + 'df_ctrl_short_inc.csv')
# df_ctrl_long_inc = pd.read_csv(data_dir + 'df_ctrl_long_inc.csv')
# df_list = {'no_inc': [df_ctrl_0_inc, ],
#            'short_inc': [df_ctrl_short_inc, ],
#             'long_inc': [df_ctrl_long_inc, ]}

# same as combining the 3 dfs together
# df_all_inc_times = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/'+\
#                                  'human/paper_figs_collected_checked/stats/correlations/'+\
#                                     'intrinsic_time_all/used_for_analysis_all_intr.xlsx')
# save_dir_corr = data_dir[0:89] + 'stats/correlations/intrinsic_time_all/'

# c_all = paper_plot.figure_colors()['slice_all_CTR']['CTR']
# df_list =  {'all_intr' : [df_all_inc_times, c_all]}

# paper_plot.plot_corrs('hrs_after_OP', df_list, save_dir_corr, save_dir_figs)

# # plot non-firing cells
# repatch_nf = pd.read_excel(data_dir+ 'repatch_non_firing_cells_count.xlsx')
# slice_nf = pd.read_excel(data_dir+ 'slice_non_firing_cells_count.xlsx')
# inc_nf = pd.read_excel(data_dir+ 'inc_only_non_firing_cells_count.xlsx')

# df_dict = {'cell_ID_new': repatch_nf,
#            'slice': slice_nf,
#            'inc_only': inc_nf}
# inj_vals = [200, 400, 600]
# save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/parts/'

# for data_type, df in df_dict.items():
#     plot_non_firing_across_inj_proportion(df, data_type, save_dir)


# # plot correlatioons
# data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
# dir_read_corrs = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/sup_tables/stats/'

# df_corrs = pd.read_csv(data_dir + 'adult_complete_0_hrs_incubation.csv')
# df_corrs = df_corrs[~((df_corrs.day == 'D2') & (df_corrs.treatment == 'high K'))]
# df_corrs = df_corrs[~((df_corrs.day == 'D2') & (df_corrs.repatch == 'yes'))]

# # # FIRING
# params_list = [ 'X..200pA....num_aps..',
#                 'X..200pA....IFF..',
#                 'X..400pA....num_aps..',
#                 'X..400pA....IFF..',
#                 'X..600pA....num_aps..',
#                 'X..600pA....IFF..']

# fir_cors = {'CTR_firing_cells_correlations': [df_corrs, 'black']}

# paper_plot.plot_corrs_firing(x = 'hrs_after_OP',
#                   df_list = fir_cors,
#                   dir_read_corrs = dir_read_corrs,
#                   save_dir_figs = destination_dir + 'fig_freq_',
#                   params_list = params_list)


# INTRINSIC

# params_mem = ['hrs_after_OP',
#                  'Rin',
#                  'resting_potential',
#                  'membra_time_constant_tau',
#                  'sag']
# df_list = {'CTR_correlations_intr': [df_corrs, 'black']}
# paper_plot.plot_corrs('hrs_after_OP', df_list, dir_read_corrs, destination_dir + 'mem', params_mem)

# params_AP = ['hrs_after_OP',
#                  'AP_heigth',
#                  'AP_halfwidth',
#                  'max_depol',
#                  'max_repol']
# df_list = {'CTR_correlations_intr': [df_corrs, 'darkgrey']}
# paper_plot.plot_corrs('hrs_after_OP', df_list, dir_read_corrs, destination_dir + 'ap', params_AP)


#%%
### PART DOS - SYNAPTIC STUFF

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/connectivity/'
# save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/parts/connectivity/'
save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/trying_out/'
df_connections_QC = pd.read_excel(data_dir + 'all_QCed_connections.xlsx')
df_connections_QC = df_connections_QC[df_connections_QC.repatch == 'no'].reset_index(drop = True)
df_con_percentage_all = pd.read_excel(data_dir + 'connectivity_percent_estimation.xlsx')

df_con_percentage_all.groupby(['treatment', 'day'])['found_connections'].agg(['sum'])
df_con_percentage_all.groupby(['treatment', 'day'])['num_possible_connections'].agg(['sum'])

# paper_plot.plot_connectivity_percentage(df_con_percentage_all, 'slice', save_dir + 'NEW_')

# OVER TIME
# plot Ctrl connectivity over time
df_CTR_con = df_con_percentage_all[(df_con_percentage_all.treatment == 'Ctrl') | \
                               ((df_con_percentage_all.treatment == 'high K') & \
                                (df_con_percentage_all.day == 'D1'))]
df_CTR_con = df_CTR_con.dropna(subset=['hrs_after_OP', 'con_percentage_fixed'])
df_CTR_con = df_CTR_con[df_CTR_con.hrs_after_OP < 41.05] # THE OTHER CONNECTIONS ARE REPATCHED

paper_plot.plot_corr_single('hrs_after_OP', 'con_percentage_fixed', df_CTR_con, 'all', 'Ctrl', save_dir + 'TEMP')

# df_con_percentage_no_0 = df_con_percentage_all[df_con_percentage_all.con_percentage_fixed != 0]
# df_con_percentage_no_0_CTR = df_con_percentage_no_0[(df_con_percentage_no_0.treatment == 'Ctrl') | \
#                                ((df_con_percentage_no_0.treatment == 'high K') & \
#                                 (df_con_percentage_no_0.day == 'D1'))]
# paper_plot.plot_corr_single('hrs_after_OP', 'con_percentage_fixed', df_con_percentage_no_0_CTR, 'all', 'Ctrl', save_dir + 'TEMP')

# Amp 1 over time in CTR
df_CTR = df_connections_QC[(df_connections_QC.treatment == 'Ctrl') | ((df_connections_QC.treatment == 'high K') & (df_connections_QC.day == 'D1'))]
df_CTR = df_CTR.dropna(subset=['hrs_after_OP', 'Lat1'])

paper_plot.plot_corr_single('hrs_after_OP', 'Amp 1', df_CTR, 'all', 'Ctrl', save_dir)
paper_plot.plot_corr_single('hrs_after_OP', 'Lat1', df_CTR, 'all', 'Ctrl', save_dir + 'lat')

## FOR PAPER
paper_plot.plot_distribution_1st_amp_lat_only(df_connections_QC , 'slice', save_dir)
paper_plot.plot_distribution_1st_amp_lat_only(df_connections_QC_lat , 'slice', save_dir + 'Lat_')



paper_plot.plot_norm_con_params(df_connections_QC, 'amp', 'slice', save_dir+'_')
paper_plot.plot_norm_con_params(df_connections_QC, 'lat', 'slice', save_dir + '_')
# paper_plot.plot_amp_lat_cumulative(df_connections_QC, 'slice', save_dir)

# AMP and LAT
paper_plot.plot_connect_params(df_connections_QC, 'slice', 'amp', save_dir + 's_')

# Remove rows with missing values or negative values in latency columns
latency_cols = ['Lat1', 'Lat2', 'Lat3', 'Lat4']
df_connections_QC_lat = df_connections_QC.dropna(subset=latency_cols)
df_connections_QC_lat = df_connections_QC_lat[(df_connections_QC_lat[latency_cols] >= 0).all(axis=1)].reset_index(drop=True)
paper_plot.plot_connect_params(df_connections_QC_lat, 'slice', 'lat', save_dir+'l')

import importlib
importlib.reload(paper_plot)

paper_plot.plot_connection_ratio_bar(df_con_percentage_all, 'slice', save_dir + 'for_defense')


# connection ratios
paper_plot.plot_connection_ratio(df_con_percentage_all, 'slice', save_dir + 'NEW')
# paper_plot.plot_connection_ratio2(df_con_percentage_no_0, 'slice', save_dir + 'no 0 cons')

# create supplementary figures
# sup_save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/sup_figures/sup_figs_parts/connectivity/'
# paper_plot.create_hist_sup_figures(df_connections_QC, 'Amp 1', 'Amplitude 1st EPSP (mV)', sup_save_dir)
# paper_plot.create_hist_sup_figures(df_connections_QC_lat, 'Lat1', 'Latency 1st EPSP (ms)', sup_save_dir)

# better function for sup figures
# paper_plot.plot_distribution_1st_amp_lat_only(df_connections_QC , 'slice', destination_dir)
# paper_plot.plot_distribution_1st_amp_lat_only(df_connections_QC_lat , 'slice', destination_dir + 'Lat_')



## REPATCHED CONNECTIONS

#### df_repatched_cons = pd.read_excel(data_dir + 'repatch_QCed_connections.xlsx')

#samefor theloosely QC repatched
df_con_repatch_col_QC = pd.read_excel(f'{data_dir}/repatched_connections_looser_QC.xlsx')
df_repatch_con_both_days = df_con_repatch_col_QC[~df_con_repatch_col_QC['comment'].str.contains('appearing|disappearing', case=False, na=False)].reset_index(drop = True)

# latency_cols = ['Lat1']
# df_repatch_con_both_days_l = df_repatch_con_both_days.dropna(subset=latency_cols)
# df_repatch_con_both_days_l = df_repatch_con_both_days_l[(df_repatch_con_both_days_l[latency_cols] >= 0).all(axis=1)].reset_index(drop=True)
# df_repatch_con_both_days_l = df_repatch_con_both_days_l.groupby('connection_ID').filter(lambda x: len(x) > 1)

# paper_plot.plot_connect_params(df_repatch_con_both_days_l, 'cell_ID_new', 'lat', save_dir + '_')

# df_con_repatch_col_QC.groupby(['treatment', 'day'])['Amp 1'].agg(['mean', 'std', 'count'])
# df_repatch_con_both_days.groupby(['treatment', 'day'])['Amp 1'].agg(['mean', 'std', 'count'])
# df_repatch_con_both_days_l.groupby(['treatment', 'day'])['Lat1'].agg(['mean', 'std', 'count'])


# # create donut plots for appeaing, disappearing
# paper_plot.create_donut_plots(size_C = 15, size_K = 28, save_dir = save_dir + 'repatch/')

# paper_plot.plot_norm_con_params(df_repatch_con_both_days, 'amp', 'cell_ID_new', save_dir + 'repatch/connected_on_both_days')
# paper_plot.plot_norm_con_params(df_repatch_con_both_days, 'lat', 'cell_ID_new', save_dir + 'repatch/connected_on_both_days')
# paper_plot.plot_amp_lat_cumulative(df_repatch_con_both_days,  'cell_ID_new', save_dir + 'repatch/connected_on_both_days')

# paper_plot.plot_connect_1st_amp_lat_only(df_repatch_con_both_days,  'cell_ID_new', save_dir + 'repatch/connected_on_both_days')
# paper_plot.plot_connect_1st_amp_lat_only(df_repatch_con_both_days,  'cell_ID_new', destination_dir)


# repatch_connect_percentage = paper_plot.get_repatch_connect_estimation(df_con_repatch_col_QC, df_con_percentage_all)
# paper_plot.plot_connectivity_percentage(repatch_connect_percentage, 'cell_ID_new', save_dir + 'repatch/')
# paper_plot.plot_connection_ratio(repatch_connect_percentage, 'cell_ID_new', save_dir + 'repatch/')
# paper_plot.plot_connection_ratio2(repatch_connect_percentage, 'cell_ID_new', save_dir + 'repatch/')


# PLOT APPEARING / DISAPPEARING
save_dir_special = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/trying_out'
# for sup fig summary
df_dis_appearing = df_con_repatch_col_QC[df_con_repatch_col_QC['comment'].str.contains('appearing|disappearing', case=False, na=False)].reset_index(drop = True)
df_dis_appearing.groupby(['treatment', 'day'])['Amp 1'].agg(['mean', 'std', 'count'])

indx_appearing = [4, 5, 9, 13]
appear_only = df_dis_appearing.iloc[indx_appearing, :]
disappear_only = df_dis_appearing.drop(indx_appearing).reset_index(drop = True)
paper_plot.plot_appear_disappear(disappear_only, appear_only, save_dir_special)


# plots the connectivity of OP240215 to be used for repatch figure
paper_plot.plot_connection_window_fig()

# # plotting appearing connection
# OP = 'OP240215'
# patcher = 'Verji'
# save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/parts/connectivity/'
# file_index_char = [13,49]
# active_channels = [[1,2,5,8],
#                 [1,2,3,8]]

# paper_plot.plot_repatch_char_D1_D2(OP, patcher, file_index_char, active_channels, save_dir_special)
# paper_plot.plot_repatch_connect_D1_D2(OP, patcher, file_index_connect, active_channels, save_dir_special)

# # plot disappearing. appearing cons
# patcher = 'Verji'
# dict_plot_cons = {
#     # 'OP220217': [[20,20], [7,7,7,7], [4,3,3,2]] # only evidence in vc
#     # 'OP240417': [[19, 53], [1,1], [5,5], ['pre', 'post HiK']],
#     # 'OP240503_S1': [[5, 44, 5, 44], [3, 5, 3, 5], [2, 2, 6, 6],\
#     #                 ['pre', 'gone post HiK', 'pre', 'remains post HiK']],
#     'OP230420': [[17, 91], [4,4], [7,8], ['pre', 'post HiK']],
#     # 'OP240503_S1': [[5, 44], [3, 5], [2, 2],\
#     #                 ['pre', 'post HiK']],
#     # 'OP240503_S2': [[18, 57], [5, 5], [2, 2], ['pre', 'post CTR']],
#     'OP240503_S3': [[33, 69], [2, 2], [5, 5], ['pre', 'post HiK']],
#     'OP240503_S3_1': [[33, 69], [6, 6], [2, 2], ['pre', 'post HiK']],
#     # 'OP220426': [[20, 74], [2, 2], [4, 4], ['pre', 'post HiK']],
#     # 'OP240215': [[15, 51], [8, 8], [1, 1], ['pre', 'post HiK']],
#     'OP240926': [[22, 70], [3, 3], [4, 4], ['pre', 'post CTR']]}

# for (OP, item) in dict_plot_cons.items():
#     if OP != 'OP240503_S3_1':
#         continue
#     fn_indx = item[0]
#     pre_chans = item[1]
#     post_chans = item[2]
#     labels = item[3]
#     if len(OP) > 8:
#         paper_plot.plot_connection_window(OP[:8], patcher, labels, fn_indx,\
#                                       pre_chans, post_chans, save_dir_special + 'all' + OP[8:])
#     else:
#         print(f'plotting {OP}')
#         paper_plot.plot_connection_window(OP, patcher, labels, fn_indx,\
#                                         pre_chans, post_chans, save_dir_special + 'QC')
#     if len(fn_indx) > 3:
#         for i in range(2):
#             if i == 1:
#                 i +=1
#             fn_indx_1 = [fn_indx[i], fn_indx[i + 1]]
#             active_channels = [[pre_chans[i], pre_chans[i + 1]],
#                                [post_chans[i], post_chans[i + 1]]]
#             paper_plot.plot_repatch_connect_D1_D2(OP, patcher, fn_indx_1, active_channels, save_dir_special)

#     else:
#         active_channels = [[pre_chans[0], post_chans[0]],
#                             [pre_chans[1], post_chans[1]]]
#         paper_plot.plot_repatch_connect_D1_D2(OP, patcher, fn_indx, active_channels, save_dir_special)


# %%

df_repatch_cons_plus_intr = pd.read_excel(f'{data_dir}all_repatched_and_intr_manually_checked.xlsx')

paper_plot.QC_plot_repatched_connections(df_repatch_cons_plus_intr, save_dir + 'repatch/')

paper_plot.QC_plot_single_cell_param_repatched_cons(df_repatch_cons_plus_intr, save_dir, ['Rin'])
paper_plot.QC_plot_single_cell_param_repatched_cons(df_repatch_cons_plus_intr, save_dir, ['Rs_end_pre'])
paper_plot.QC_plot_single_cell_param_repatched_cons(df_repatch_cons_plus_intr, save_dir, ['Rs_end_post'])

importlib.reload(paper_plot)
paper_plot.QC_plot_Rs_vs_EPSP_amp(df_repatch_cons_plus_intr, save_dir_special+ '/')

## PLOT FAILURE RATES SLICE
data_dir_ = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/connectivity/'
df_slice_newest = pd.read_excel(f'{data_dir_}/all_slice_cons_with_failure_rates.xlsx')
paper_plot.plot_failure_rates_slice(df_slice_newest, save_dir + '/ALL_slice_')


#### PLOT FAILURE RATES repatch
df_repatch_newest_ = pd.read_excel(f'{data_dir}repatched_connections_looser_QC_MARCH_2026.xlsx')
df_repatch_newest = df_repatch_newest_[~df_repatch_newest_['comment'].str.contains('appearing|disappearing', case=False, na=False)].reset_index(drop = True)

paper_plot.plot_failure_rates_repatch(df_repatch_newest_, save_dir_special + '/ALL')



# %%

# some extra mini/spontan plots
df_minis_all = pd.read_excel('/Users/verjim/Desktop/minid_df_complete.xlsx')
save_dir_mini_corr = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/sup_figures/sup_figs_parts/'

mini_CTR = df_minis_all[df_minis_all.treatment == 'Ctrl']
paper_plot.plot_corr_single('hrs_incubation', 'avg_amp_adj', mini_CTR, 'minis', 'Ctrl', save_dir_mini_corr + 'minis_CTR_only_')
paper_plot.plot_corr_single('hrs_incubation', 'avg_IEI_adj', mini_CTR, 'minis', 'Ctrl', save_dir_mini_corr + 'minis_CTR_only_')

mini_HiK = df_minis_all[df_minis_all.treatment == 'high K']
paper_plot.plot_corr_single('hrs_incubation', 'avg_amp_adj', mini_HiK , 'minis', 'high K', save_dir_mini_corr + 'minis_HiK_only_')
paper_plot.plot_corr_single('hrs_incubation', 'avg_IEI_adj', mini_HiK , 'minis', 'high K', save_dir_mini_corr + 'minis_HiK_only_')


# plot trace spontan

import ephys_analysis.funcs_plotting_raw_traces as raw_plot
op_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP240926/'
raw_plot.plot_trace(op_dir + '24926037.abf', 2, 3, save_dir)

# plot pH // osmol
paper_plot.plot_ph_osmol(save_dir)


####
# PLOT HIGH K AP fired trace
from ephys_analysis.funcs_sorting import get_OP_metadata

importlib.reload(paper_plot)
OP = 'OP221116'
human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
work_dir, filenames, indices_dict, slice_names = get_OP_metadata(human_dir, OP, 'Verji', 'old')


save_trace_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
        'paper_figs_collected_checked/figures/trying_out/'
paper_plot.plt_washin_trace(work_dir + filenames[31], 5, [6], 0, 400_000, 
                            save_dir = save_trace_dir, ticks = True)

paper_plot.plt_washin_trace(work_dir + filenames[31], 5, [6], 377_450, 377_750,
                            save_dir = save_trace_dir + 'AP_')

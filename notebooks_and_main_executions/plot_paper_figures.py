import pandas as pd
import ephys_analysis.funcs_plot_paper_figures as paper_plot

#%%
#start from fig 2 --> decide size of text, itcks,labels, ratios
data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
    'paper_figs_collected_checked/data/'
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
        'paper_figs_collected_checked/figures/draft5_figs/parts/'

# # fig 1
# df_mea = pd.read_excel(data_dir + 'summary_pMEA_r.xlsx')
# mea_fig1(df_mea, destination_dir)

# # fig 2
# df_slice = pd.read_excel(data_dir + 'slice_data_temporal.xlsx')
# firing_props_fig2_slice(df_slice, 'firing_cells_slice', 'num_aps', 'SE', destination_dir)
# # firing_props_fig2(df_slice, 'firing_cells_slice', 'num_aps', 'CI', destination_dir)

# cell_IDs_slice = {'pre': '24201S1c8',
#             'post CTR': '24117S2_D2c1',
#             'post HiK': '23420S2_D2c7'}
# plot_example_firing('slice', df_slice, destination_dir, cell_IDs_slice)
# intr_params_fig2(df_slice, destination_dir)

# # fig 3
# df_repatch = pd.read_excel(data_dir + 'repatch_data_temporal.xlsx')
# intr_params_fig3(df_repatch, destination_dir)
# cell_IDs_repatch = {'CTR': 'vm24321S2c8',
#             'HiK': 'vm24926S4c1'}
# plot_example_firing_fig3('repatch', df_repatch, destination_dir, cell_IDs_repatch)
# firing_props_fig3(df_repatch, 'SE', destination_dir)
# firing_props_fig3(df_repatch, 'CI', destination_dir)

# # slice CTR all vs hrs_after_OP
# intr_slice_CTR(destination_dir, er_type = 'SE')
# intr_slice_CTR(destination_dir, er_type = 'CI') # doesn't complete teh figure

# df_slice_all = pd.read_excel(data_dir + 'slice_all.xlsx')
# slice_all_CTR_hist(df_slice_all, destination_dir)
   
# fig 5
# df_AIS = pd.read_excel(data_dir + 'AIS_all_data.xlsx')
# AIS_fig5(df_AIS, destination_dir)

# # incubation only
df_slice_short_inc = pd.read_excel(data_dir + '10.07.25_slice_incubation_only_no_move.xlsx')
df_slice_short_inc = df_slice_short_inc[df_slice_short_inc['hrs_after_OP'] < 35]
df_slice_short_inc = df_slice_short_inc[df_slice_short_inc['hrs_incubation'] > 15.9]
df_slice_short_inc.reset_index(drop = True, inplace = True)
df_slice_short_inc = df_slice_short_inc.loc[[i for i, sl in enumerate(df_slice_short_inc.slice) if len(sl) <= 2], :]
# intr_params_inc_only(df_slice_short_inc, destination_dir)
# firing_props_inc_only_V2(df_slice_short_inc, 'SE', destination_dir)

# cell_IDs_inc = {'Ctrl': '25220S3c3',
#             'high K': '25220S2c3'}
# plot_example_firing('inc_only', df_slice_short_inc, destination_dir, cell_IDs_inc)

# FIG - correlations
data_dir =  '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
save_dir_figs = data_dir[0:89] +  '/figures/draft5_figs/parts/'
x = 'hrs_after_OP'

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

# c_all = figure_colors()['slice_all_CTR']['CTR']
# df_list =  {'all_intr' : [df_all_inc_times, c_all]}
# plot_corrs(x, df_list, save_dir_corr, save_dir_figs)

# 0 hrs_incubation, 5 hrs after OP
# 17 hrs_incubation, 24 hrs after OP
# 24 hrs_incubation, 48 hrs after OP
# plot_example_firing_corr('slice', adult_df_slice_all, dest_dir,
#                     {'pre': '23n09S1c6',
#                     'Ctrl': '25207S3c5',
#                      'CtrlD2': '22427S4c4'})


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

# plot correlatioons
data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
df_cors_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/sup_tables/stats/'
save_plot_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/parts/'

df_corrs = pd.read_csv(data_dir + 'adult_complete_0_hrs_incubation.csv')
df_corrs = df_corrs[~((df_corrs.day == 'D2') & (df_corrs.treatment == 'high K'))]
df_corrs = df_corrs[~((df_corrs.day == 'D2') & (df_corrs.repatch == 'yes'))]

# FIRING
params_list = [ 'X..200pA....num_aps..',
                'X..200pA....IFF..',
                'X..400pA....num_aps..',
                'X..400pA....IFF..',
                'X..600pA....num_aps..',
                'X..600pA....IFF..']

fir_cors = {'CTR_firing_cells_correlations': [df_corrs, 'pink']}

# plot_corrs_firing(x = 'hrs_after_OP',
#                   df_list = fir_cors,
#                   dir_read_corrs = df_cors_dir,
#                   save_dir_figs = save_plot_dir + 'fig_freq_',
#                   params_list = params_list)


# INTRINSIC
dir_read_corrs = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/sup_tables/stats/'

params = ['hrs_after_OP',
                 'Rin',
                 'resting_potential',
                 'rheo_ramp_c',
                 'TH', 
                 'AP_halfwidth',
                 'max_depol',
                 'max_repol',
                 'sag',
                 'AP_heigth',
                 'membra_time_constant_tau']
df_list = {'CTR_correlations_intr': [df_corrs, 'blue']}

# plot_corrs('hrs_after_OP', df_list, dir_read_corrs, save_plot_dir , params)







#%%
### PART DOS SYNAPTIC STUFF
data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/connectivity/'
save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/parts/connectivity/'

df_connections_QC = pd.read_excel(data_dir + 'all_QCed_connections.xlsx')
df_connections_QC = df_connections_QC[df_connections_QC.repatch == 'no']
df_con_percentage_all = pd.read_excel(data_dir + 'connectivity_percent_estimation.xlsx')

paper_plot.plot_connectivity_percentage(df_con_percentage_all, 'slice', save_dir)

# OVER TIME
# plot Ctrl connectivity over time
df_CTR_con = df_con_percentage_all[(df_con_percentage_all.treatment == 'Ctrl') | \
                               ((df_con_percentage_all.treatment == 'high K') & \
                                (df_con_percentage_all.day == 'D1'))]
df_CTR_con = df_CTR_con.dropna(subset=['hrs_after_OP', 'con_percentage_fixed'])
paper_plot.plot_corr_single('hrs_after_OP', 'con_percentage_fixed', df_CTR_con, 'all', 'Ctrl', save_dir)

# Amp 1 over time in CTR
df_CTR = df_connections_QC[(df_connections_QC.treatment == 'Ctrl') | ((df_connections_QC.treatment == 'high K') & (df_connections_QC.day == 'D1'))]
df_CTR = df_CTR.dropna(subset=['hrs_after_OP', 'Amp 1'])
paper_plot.plot_corr_single('hrs_after_OP', 'Amp 1', df_CTR, 'all', 'Ctrl', save_dir)

paper_plot.plot_norm_con_params(df_connections_QC, 'amp', 'slice', save_dir)
paper_plot.plot_norm_con_params(df_connections_QC, 'lat', 'slice', save_dir)
paper_plot.plot_amp_lat_cumulative(df_connections_QC, 'slice', save_dir)

# AMP and LAT
paper_plot.plot_connect_params(df_connections_QC, 'slice', 'amp', save_dir)
# Remove rows with missing values or negative values in latency columns
latency_cols = ['Lat1', 'Lat2', 'Lat3', 'Lat4']
df_connections_QC_lat = df_connections_QC.dropna(subset=latency_cols)
df_connections_QC_lat = df_connections_QC_lat[(df_connections_QC_lat[latency_cols] >= 0).all(axis=1)].reset_index(drop=True)
paper_plot.plot_connect_params(df_connections_QC_lat, 'slice', 'lat', save_dir)

# connection ratios
paper_plot.plot_connection_ratio(df_con_percentage_all, 'slice', save_dir)
paper_plot.plot_connection_ratio2(df_con_percentage_all, 'slice', save_dir)




## REPATCHED CONNECTIONS

# df_repatched_cons = pd.read_excel(data_dir + 'repatch_QCed_connections.xlsx')
# plot_connect_params(df_repatched_cons, 'cell_ID_new', 'amp', save_dir)


#samefor theloosely QC repatched
df_con_repatch_col_QC = pd.read_excel(f'{data_dir}/repatched_connections_looser_QC.xlsx')
con_IDs_appearing_disapearing_cons = ['22427S3c2#4',
                                      '24503S1c3#2',
                                      '24503S2c5#2']

df_con_repatch_col_QC = df_con_repatch_col_QC[~df_con_repatch_col_QC.connection_ID\
                                              .isin(con_IDs_appearing_disapearing_cons)].reset_index(drop = True)

df_con_repatch_col_QC.groupby(['treatment', 'day'])['Amp 1'].agg(['mean', 'std', 'count'])

paper_plot.plot_norm_con_params(df_con_repatch_col_QC, 'amp', 'cell_ID_new', save_dir + 'repatch/')
paper_plot.plot_norm_con_params(df_con_repatch_col_QC, 'lat', 'cell_ID_new', save_dir + 'repatch/')
paper_plot.plot_amp_lat_cumulative(df_con_repatch_col_QC, 'cell_ID_new', save_dir  + 'repatch/')

paper_plot.plot_connect_params(df_con_repatch_col_QC, 'cell_ID_new', 'amp', save_dir + 'repatch/')
paper_plot.plot_connect_params(df_con_repatch_col_QC, 'cell_ID_new', 'lat', save_dir + 'repatch/')

paper_plot.plot_connect_1st_amp_lat_only(df_con_repatch_col_QC,  'cell_ID_new', save_dir = save_dir + 'repatch/')










# PLOT APPEARING / DISAPPEARING

save_dir_special = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/disappearing_connections/'
# plotting appearing connection
OP = 'OP240215'
patcher = 'Verji'
save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/monograph/parts/connectivity/'
file_index_char = [13,49]
active_channels = [[1,2,5,8],
                [1,2,3,8]]

paper_plot.plot_repatch_char_D1_D2(OP, patcher, file_index_char, active_channels, save_dir_special)
paper_plot.plot_repatch_connect_D1_D2(OP, patcher, file_index_connect, active_channels, save_dir_special)


# plot disappearing. appearing cons
patcher = 'Verji'
dict_plot_cons = {
    # 'OP220217': [[20,20], [7,7,7,7], [4,3,3,2]] # only evidence in vc
    'OP240417': [[19, 53], [1,1], [5,5], ['pre', 'post HiK']],
    'OP240503_S1': [[5, 44, 5, 44], [3, 5, 3, 5], [2, 2, 6, 6],\
                    ['pre', 'gone post HiK', 'pre', 'remains post HiK']],
    'OP240503_S2': [[18, 57], [5, 5], [2, 2], ['pre', 'post CTR']],
    'OP220426': [[20, 74], [2, 2], [4, 4], ['pre', 'post HiK']],
  #  'OP240215': [[15, 51], [1,2,5,8], [1,2,5,8]]
}

import importlib
importlib.reload(paper_plot)

for (OP, item) in dict_plot_cons.items():
    fn_indx = item[0]
    pre_chans = item[1]
    post_chans = item[2]
    labels = item[3]
    if len(OP) > 8:
        OP = OP[:8]
    print(f'plotting {OP}')
    paper_plot.plot_connection_window(OP, patcher, labels, fn_indx,\
                                      pre_chans, post_chans, save_dir_special)
    if len(fn_indx) > 3:
        for i in range(2):
            if i == 1:
                i +=1
            fn_indx_1 = [fn_indx[i], fn_indx[i + 1]]
            active_channels = [[pre_chans[i], pre_chans[i + 1]],
                               [post_chans[i], post_chans[i + 1]]]
            paper_plot.plot_repatch_connect_D1_D2(OP, patcher, fn_indx_1, active_channels, save_dir_special)

    else:
        active_channels = [[pre_chans[0], post_chans[0]],
                            [pre_chans[1], post_chans[1]]]
        paper_plot.plot_repatch_connect_D1_D2(OP, patcher, fn_indx, active_channels, save_dir_special)


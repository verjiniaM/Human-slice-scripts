import importlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ephys_analysis.funcs_plot_conn_paper_figures as con_plots
import ephys_analysis.funcs_plot_paper_figures as paper_plot
import ephys_analysis.funcs_con_paper_analysis_as_FM as funcs_FM

#importlib.reload(con_plots)

work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/'
save_plots_dir = f'{work_dir}results/human/connectivity_paper/plots/figs_draft2/'
VM_data_dir = f'{work_dir}results/human/connectivity_paper/data/'
stats_dir = f'{work_dir}results/human/connectivity_paper/stats/for_plotting/'

plt.style.use('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/Human-slice-scripts/style_plot_paper.mplstyle')

#%%

### CONNECTIVITY
con_percent_saved = pd.read_csv(f'{VM_data_dir}con_percentages_from_cross_checked_with_intr.csv')
# df_con_percentage_all = con_percent_saved[(con_percent_saved.treatment == 'high K') |
#                                            (con_percent_saved.treatment == 'Ctrl')].reset_index(drop = True)
df_con_percentage_all = con_percent_saved[(con_percent_saved.treatment != 'fixed')].reset_index(drop = True)

con_plots.plot_connection_ratio_bar(df_con_percentage_all, 'slice', save_plots_dir +'TTX')
con_plots.plot_connectivity_percentages(df_con_percentage_all, save_plots_dir +'TTX')


### EVENT PARAMETERS
all_event_props = pd.read_csv(f'{VM_data_dir}all_connections_QC_checked_complete.csv')
all_event_props = all_event_props[(all_event_props.day =='D1') |\
                                  ((all_event_props.day =='D2') & (all_event_props.hrs_incubation > 16))]
df_slice = all_event_props[all_event_props.dataset != 'repatch']

df_repatch_newest = pd.read_csv(f'{VM_data_dir}repatched_connections_looser_QC_Nov_2026.csv')
df_repatch_no_app_disap = df_repatch_newest[~df_repatch_newest['comment'].str.contains('appearing|disappearing', case=False, na=False)].reset_index(drop = True)


# PLOT BOXPLOTS AMP / LAT
# for repatch- using all possible connections that are repatched (incl. disap/appearing)

dfs_dict = {
    'slice': df_slice,
    'cell_ID_new': df_repatch_newest,
    'lat_specific': df_repatch_no_app_disap
}
# plot boxplots, descriptive
counts = {}
for data_type, df in dfs_dict.items():

    if data_type == 'lat_specific':
        con_plots.plot_connect_params(df, 'cell_ID_new', 'lat', save_plots_dir +'_for_fig')
        continue
    elif data_type == 'cell_ID_new':
        con_plots.plot_connect_params(df, data_type, 'amp', save_plots_dir +'_for_fig')
    else:
        con_plots.plot_connect_params(df, data_type, 'amp', save_plots_dir +'_for_fig')
        con_plots.plot_connect_params(df, data_type, 'lat', save_plots_dir +'_for_fig')

    # sanity check for connections per OP
    df = df.copy()
    df['condition'] = np.where(df['day'] == 'D1', 'RS0', df['treatment'])
    counts[data_type] = (df.groupby(['OP', 'condition'])['connection_ID']
        .nunique().unstack(fill_value=0))


# PLOT connection properties
data_type_dict = {
    'slice' : 'slice_ID',
    'cell_ID_new': 'connection_ID'
}
for data_name, data_type in data_type_dict.items():
    TTX = 'TTX' if data_type == 'slice_ID' else ''

    # DESCRIPTIVE vs LLM ESTIMATES
    for var in ('amp', 'lat'):
        descriptive = pd.read_csv(f'{stats_dir}{TTX}EPSP_paired_log_95CI_bootstrap_{data_type}{var}.csv')
        model = pd.read_csv(f'{stats_dir}{TTX}EPSP_paired_log_95CI_from_model_{data_type}{var}.csv')

        # main figure
        con_plots.plot_norm_con_params_log_ci(descriptive, var, data_name, save_plots_dir)

        # sup. figure
        fig, axes = con_plots.plot_ppr_comparison(descriptive, model, data_type, DV=var,  save_dir=save_plots_dir)
        fig.savefig(f'{save_plots_dir}sup_{TTX}PPR_comparison_{data_type}{var}.svg', bbox_inches='tight')
        plt.show()

    # DONOR CONTRIBUTION HEATMAP
    df_long = pd.read_csv(f'{stats_dir}{TTX}EPSP_paired_log_data_{data_type}amp.csv')
    con_plots.plot_ppr_donor_heatmap(df_long, data_type, save_dir=save_plots_dir)
    
    # sanity check
    print(counts[data_name])

    # FACILITATING/ DEPRESSING
    df_fac = pd.read_csv(f'{stats_dir}{TTX}EPSP_facilitation_CI_binomial_observed{data_type}.csv')
    con_plots.plot_con_facil(df_fac, data_name, save_plots_dir)
    descriptive_amp = pd.read_csv(f'{stats_dir}{TTX}EPSP_paired_log_95CI_bootstrap_{data_type}amp.csv')
    con_plots.plot_con_ppr_scatter_ci(dfs_dict[data_name], descriptive_amp, data_name, save_plots_dir)



### FAILURE RATES
# slice
con_plots.plot_failure_rates_slice(df_slice, save_plots_dir)
# repatch
con_plots.plot_failure_rates_repatch(df_repatch_no_app_disap, save_plots_dir)


dfs_dict = {'slice': df_slice,
            'cell_ID_new': df_repatch_no_app_disap}


df_fit_repatch = pd.read_excel(f'{work_dir}results/human/connectivity_paper/stats/connection_ID_amplitude_failure.xlsx', sheet_name='Amp_fitted')
df_fit_slice = pd.read_excel(f'{work_dir}results/human/connectivity_paper/stats/slice_ID_amplitude_failure.xlsx', sheet_name='Amp_fitted')

df_slice['failure_rate'] = df_slice.failure_rate / 100
con_plots.plot_amp1_vs_failure_rate(df_repatch_no_app_disap, df_fit_repatch, 'cell_ID_new', save_dir=save_plots_dir)
con_plots.plot_amp1_vs_failure_rate(df_slice, df_fit_slice, 'slice', save_dir=save_plots_dir)


#%%
### E TYPES
df_FM = pd.read_csv(f'{work_dir}data/human/humandata/data_tables/tcell_all_psdn.csv')
df_FM = df_FM[(df_FM.removal_boolean == 1) & \
            (df_FM.is_pc == 1) & \
                (df_FM.piadistance < 1200) & \
                    (df_FM.type == 1)]
df_VM_etypes = pd.read_csv(f'{VM_data_dir}VM_intr_with_etype_QCed.csv')

# MERGE ETYPES FOR CONNECTED CELLS
df_con_VM_etypes = all_event_props.copy()
for side in ['pre', 'post']:
    df_con_VM_etypes = df_con_VM_etypes.merge(
        df_VM_etypes[['cell_IDs_match', 'e_type']].drop_duplicates('cell_IDs_match'),
        left_on=f'cell_IDs_match_{side}', right_on='cell_IDs_match',
        how='left', validate='m:1'
    ).rename(columns={'e_type': f'e_type_{side}'}).drop(columns='cell_IDs_match')

# keep only where post cell
df_con_post_etypes = df_con_VM_etypes.dropna(subset=['e_type_post']).copy()

# keep where both
df_con_pre_post_etypes = df_con_post_etypes.dropna(subset=['e_type_pre']).copy()


# POSSIBLE SUPPLEMANTARY ANALySIS POST ETYPE FACIL?DEPRESSING CONNECTION
# # plot PPR per etypes post
# fig, axs = con_plots.plot_ppr_etypes(df_con_post_etypes, DV='amp')
# con_plots.plot_con_ppr_scatter_etypes(df_con_post_etypes, save_dir=False)

df_all_possible_cons = funcs_FM.get_all_possible_connection_from_intr(df_VM_etypes, df_con_VM_etypes)

# PLOT CONNECTION MATRIX

fig, axes = plt.subplots(1, 5, figsize=(10, 2.5)) # full connection matr
fig_small, axs_small = plt.subplots(1, 5, figsize=(16/2.54, 5/2.54)) # small connect matr

conditions = ['D1', 'Ctrl', 'high K', 'TTX']
for i, condition in enumerate(conditions):
    if condition == 'D1':
        df_plot = df_all_possible_cons[df_all_possible_cons['day'] == 'D1'].copy()
        con_matr = funcs_FM.get_connection_matrix(df_plot, df_FM)
        bl_matrix = con_matr.copy()
    else:
        df_plot = df_all_possible_cons[(df_all_possible_cons['day'] == 'D2') &
                             (df_all_possible_cons['treatment'] == condition)].copy()

    con_matr = funcs_FM.get_connection_matrix(df_plot, df_FM)
    
    con_plots.plot_connection_matrix(con_matr, title=condition, ax=axes[i], vmax=25, colorbar=False)

    con_plots.plot_connectivity_change_matrix(connection_matrix = con_matr.loc[['SlowAP', 'LowRin'], ['SlowAP', 'LowRin']],\
        baseline_matrix = bl_matrix.loc[['SlowAP', 'LowRin'], ['SlowAP', 'LowRin']], \
                                          title=condition, ax=axs_small[i], colorbar=False)

con_plots.plot_connection_matrix(con_matr, title='CTR', ax=axes[4], vmax=25, colorbar=True)
con_plots.plot_connectivity_change_matrix(connection_matrix = con_matr.loc[['SlowAP', 'LowRin'], ['SlowAP', 'LowRin']],\
    baseline_matrix = bl_matrix.loc[['SlowAP', 'LowRin'], ['SlowAP', 'LowRin']], \
                                        title=condition, ax=axs_small[4], colorbar=False)

fig.tight_layout()
fig_small.tight_layout()

fig.savefig(f'{save_plots_dir}connection_etype_NEW.svg', bbox_inches='tight')
fig_small.savefig(f'{save_plots_dir}connection_etype_changes.svg', bbox_inches='tight')
plt.show()


# %%
# Plot changes in synapse size per condition and e_type

etype_names = (df_FM[['etype_idx', 'etype']].drop_duplicates().
                set_index('etype_idx')['etype'].to_dict())


# selecting only slow AP and Rin etypes
df_con_slow_Rin = df_con_VM_etypes[(df_con_VM_etypes.e_type_pre == 2) |
                                   (df_con_VM_etypes.e_type_pre == 5)]
df_con_slow_Rin = df_con_slow_Rin[(df_con_slow_Rin.e_type_post == 2) |
                                  (df_con_slow_Rin.e_type_post == 5)]
df_con_slow_Rin = df_con_slow_Rin[df_con_slow_Rin['Amp 1'] > 0]


ETYPES = [2, 5]


# PLOT THE AMPLITUDES PER CONDITION
for param in ['Amp 1', 'Lat1']:
    fig_param, axes = plt.subplots(1, 5, figsize=(8, 2))

    for i, condition in enumerate(conditions):
        if condition == 'D1':
            amp_df, n_amp_df = funcs_FM.get_first_param_matrix(df_con_slow_Rin[df_con_slow_Rin.day == 'D1'],
                                                                df_FM, param, ETYPES)
            bl_matr = amp_df
        else:
            amp_df, n_amp_df = funcs_FM.get_first_param_matrix(df_con_slow_Rin[(df_con_slow_Rin.day == 'D2') &
                                                                                (df_con_slow_Rin.treatment == condition)],
                                                                df_FM, param, ETYPES)

        con_plots.plot_param_matrix(amp_df, n_amp_df, param, title = condition, ax = axes[i],\
                                colorbar = False, baseline_matrix=bl_matr)


    con_plots.plot_param_matrix(amp_df, n_amp_df, param, title = 'RS1', ax = axes[4],\
                            colorbar = True, baseline_matrix=bl_matr)

    fig_param.tight_layout()
    plt.savefig(f'{save_plots_dir}{param}s.svg',
        format='svg', bbox_inches='tight', dpi=300)
    plt.show()




# check proportions of. etypes in the different groups
df_VM_etypes = df_VM_etypes[df_VM_etypes.treatment != 'fixed']
etype_prop = (
    df_VM_etypes.groupby(['day', 'treatment'])['e_type']
    .value_counts(normalize=True)
    .mul(100)
    .unstack(fill_value=0)
    .reindex(columns=[2, 3, 4, 5], fill_value=0)
)

print(etype_prop)
# plot them
con_plots.plot_etype_proportions(df_VM_etypes, df_FM, save_plots_dir)



#%%
# TRACES 
OP = 'OP240926' # 'OP220602'
patcher = 'Verji' #'Verji'
file_index = [46, 90] # [26, 69]
active_channels = [[1,2,3,4,5,6,7], [1,2,3,4,5,6,7]] # [[1,2,4,5,6,7], [3,4,6,7, 8]]

importlib.reload(con_plots)

con_plots.plot_connect_D1_D2(OP, patcher, file_index, active_channels, save_plots_dir)


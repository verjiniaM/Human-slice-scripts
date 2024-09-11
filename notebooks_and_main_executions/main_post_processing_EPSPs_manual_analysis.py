import pandas as pd
import intrinsic_props_and_connectivity.funcs_sorting as sort
import intrinsic_props_and_connectivity.funcs_plotting_raw_traces as plotting_funcs

#ANALYSIS MAIN

meta_keep_path = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/EPSPs/meta_dfs/meta_keep.xlsx'
results_QC_path = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/EPSPs/results_df/QC_ed/'
verji_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/'

results_df = sort.create_full_results_df(meta_keep_path, results_QC_path, verji_dir)
results_df  = sort.get_time_of_recording_sweeps(results_df)

results_df_plot = plotting_funcs.QC_filter_for_plotting(results_df, holding = 'no', temp = 'no', min_age = 10)
results_df_plot = plotting_funcs.QC_RMP_Ctrl(results_df_plot, max_allowed_RMP_Ctrl = -50)

plotting_funcs.plot_ (results_df_plot, 'age > 18, no holding, no temp, RMP Ctrl < -50 mV', params = ['Average amp (positive)', 'Average interevent interval (ms)', 'resting_potential'])

#%%
'''
Try to see if there's a difference if Ctrl first or later on
'''
Ctrl_first = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/EPSPs/results_df/QC_ed/Ctrl_then_high_K.xlsx')
Ctrl_plot = plotting_funcs.QC_filter_for_plotting(Ctrl_first, holding = 'no', temp = 'no')
Ctrl_plot = plotting_funcs.QC_RMP_Ctrl(Ctrl_plot, max_allowed_RMP_Ctrl = -50)

plotting_funcs.plot_ (Ctrl_plot, 'Ctrl aCSF --> High K, age > 18, no temp, RMP Ctrl < -50 mV', params = ['Average amp (positive)', 'Average interevent interval (ms)', 'resting_potential'])


high_K_first = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/EPSPs/results_df/QC_ed/high_k_then_ctrl_aCSF.xlsx')
high_K_plot = plotting_funcs.QC_filter_for_plotting(high_K_first, holding = 'no', temp = 'no')
high_K_plot = plotting_funcs.QC_RMP_Ctrl(high_K_plot, max_allowed_RMP_Ctrl = -50)


plotting_funcs.plot_(high_K_plot, 'high K --> Ctrl aCSF, age > 18, no temp, RMP Ctrl < -50 mV', params = ['Average amp (positive)', 'Average interevent interval (ms)', 'resting_potential'])


#%%
#plottig of the MEA data

df_MEA = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/EPSPs_highK/MEA_Overview_MIchael_sns.xlsx')
df_MEA_Ctrl_storage = df_MEA[df_MEA['Storage'] == 'CTRL']

plotting_funcs.sns_plot_MEA_data(df_MEA, 'Absolute change in detected spikes')





# median_amp_no_hold_no_temp_highK = results_df['amplitude median'][(results_df['holding_minus_70_y_o_n'] == 'no') &\
#     (results_df['K_concentration'] == 8) &  (results_df['recording_in'] == 'high K')].values
# median_amp_no_hold_no_temp_Ctrl = results_df['amplitude median'][(results_df['holding_minus_70_y_o_n'] == 'no') &\
#     (results_df['K_concentration'] == 8) &  (results_df['recording_in'] == 'Ctrl')].values
# #remove APs
# median_amp_no_hold_no_temp_Ctrl = median_amp_no_hold_no_temp_Ctrl[median_amp_no_hold_no_temp_Ctrl < 10] 

# freq_no_hold_no_temp_highK = results_df['frequency'][(results_df['holding_minus_70_y_o_n'] == 'no') &\
#     (results_df['K_concentration'] == 8) &  (results_df['recording_in'] == 'high K')].values
# freq_no_hold_no_temp_Ctrl = results_df['frequency'][(results_df['holding_minus_70_y_o_n'] == 'no') &\
#     (results_df['K_concentration'] == 8) &  (results_df['recording_in'] == 'Ctrl')].values

median_amp_no_hold_no_temp_Ctrl = results_df_plot['amplitude median'][results_df['recording_in'] == 'Ctrl'].values
median_amp_no_hold_no_temp_highK = results_df_plot['amplitude median'][results_df['recording_in'] == 'high K'].values
freq_no_hold_no_temp_Ctrl = results_df_plot['frequency'][results_df['recording_in'] == 'Ctrl'].values
freq_no_hold_no_temp_highK = results_df_plot['frequency'][results_df['recording_in'] == 'high K'].values

fig, ax = plt.subplots(2,1, sharex = False, figsize=(12,10))
x1 = np.linspace(0, 1,len(median_amp_no_hold_no_temp_Ctrl))
x2 = np.linspace(2, 3,len(median_amp_no_hold_no_temp_highK))
ax[0].scatter(x1, median_amp_no_hold_no_temp_Ctrl, alpha = 0.7, label = 'Ctrl')
ax[0].scatter(x2, median_amp_no_hold_no_temp_highK, alpha = 0.7, label = 'high K')
ax[0].plot([0.4, 0.6], [np.nanmean(median_amp_no_hold_no_temp_Ctrl), np.nanmean(median_amp_no_hold_no_temp_Ctrl)], c = 'k')
ax[0].plot([2.4, 2.6], [np.nanmean(median_amp_no_hold_no_temp_highK), np.nanmean(median_amp_no_hold_no_temp_highK)], c = 'k')

ax[0].text(0.5, np.nanmean(median_amp_no_hold_no_temp_Ctrl) + 0.07,  str(round(np.nanmean(median_amp_no_hold_no_temp_Ctrl),2)), size = 15, c = 'k')
ax[0].text(2.5, np.nanmean(median_amp_no_hold_no_temp_highK) + 0.07,  str(round(np.nanmean(median_amp_no_hold_no_temp_highK),2)), size = 15, c = 'k')
ax[0].set_title('Median EPSP amplitude for cell per condition, no APs')

cell_IDs = results_df_plot['cell_ID'][results_df['recording_in'] == 'Ctrl'].values
for c, cell in enumerate(cell_IDs):
    #indx = index_s[c]
    #x_K = results_df_plot['x'].loc[results_df_plot['cell_ID'] == cell].tolist()[0]
    x = [x1[c],x2[c]]
    y = [median_amp_no_hold_no_temp_Ctrl[c], median_amp_no_hold_no_temp_highK[c]]
    #op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
    ax[0].plot(x, y, '-', color = 'darkgrey', alpha = 0.5, linewidth = 2, zorder = 1)

x3 = np.linspace(0, 1,len(freq_no_hold_no_temp_Ctrl))
x4 = np.linspace(2, 3,len(freq_no_hold_no_temp_highK))
ax[1].scatter(x3, freq_no_hold_no_temp_Ctrl, alpha = 0.7, label = 'Ctrl')
ax[1].scatter(x4, freq_no_hold_no_temp_highK, alpha = 0.7, label = 'high K')
ax[1].plot([0.4, 0.6], [np.nanmean(freq_no_hold_no_temp_Ctrl), np.nanmean(freq_no_hold_no_temp_Ctrl)], c = 'k')
ax[1].plot([2.4, 2.6], [np.nanmean(freq_no_hold_no_temp_highK), np.nanmean(freq_no_hold_no_temp_highK)], c = 'k')

ax[1].text(0.5, np.nanmean(freq_no_hold_no_temp_Ctrl) + 0.1,  str(round(np.nanmean(freq_no_hold_no_temp_Ctrl),2)), size = 15, c = 'k')
ax[1].text(2.5, np.nanmean(freq_no_hold_no_temp_highK) + 0.1,  str(round(np.nanmean(freq_no_hold_no_temp_highK),2)), size = 15, c = 'k')
ax[1].set_title('Median EPSP frequency for cell per condition, no APs')

cell_IDs = results_df_plot['cell_ID'][results_df['recording_in'] == 'Ctrl'].values
for c, cell in enumerate(cell_IDs):
    #indx = index_s[c]
    #x_K = results_df_plot['x'].loc[results_df_plot['cell_ID'] == cell].tolist()[0]
    x = [x3[c],x4[c]]
    y = [freq_no_hold_no_temp_Ctrl[c], freq_no_hold_no_temp_highK[c]]
    #op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
    ax[1].plot(x, y, '-', color = 'darkgrey', alpha = 0.5, linewidth = 2, zorder = 1)

ax[0].set_ylabel('Amplitude (mV)', fontsize = 15)
ax[1].set_ylabel('Frequency (Hz)', fontsize = 15)

ax[0].tick_params(axis='y', labelsize=15)
ax[0].tick_params(axis='x', labelsize=15)
ax[1].tick_params(axis='y', labelsize=15)
ax[1].tick_params(axis='x', labelsize=15)


fig.tight_layout()
fig.patch.set_facecolor('white')
fig.legend()

plt.show()


# #seabornb plots to see trend

# sns.lineplot(x='recording_in', y='amplitude median', data=results_df_plot, palette='viridis', size='cell_ID', sizes=list(np.ones(24)+1),legend=False, errorbar=None)
# plt.show()
# sns.lineplot(x='recording_in', y='frequency', data=results_df_plot, palette='viridis', size='cell_ID', sizes=list(np.ones(24)+1),legend=False, errorbar=None)
# plt.show()


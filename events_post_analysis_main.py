import pandas as pd
import events_post_analysis_funcs as event_funcs
import plot_intrinsic_props as plot_intr
import datetime

#%%
#Pre-processing steps

# event_funcs.copy_event_analysis_data_to_analysis_folder('spontan')
# spontan_new = event_funcs.post_events_analysis_add_metadata('spontan')

#event_funcs.copy_event_analysis_data_to_analysis_folder('minis')
# mini_new = event_funcs.post_events_analysis_add_metadata('minis')

#%%
results_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/results/'

spontan_old = pd.read_excel(results_dir + 'spontan_all_01.06.2023.xlsx')
#spontan_df_orig = pd.concat([spontan_old.loc[:], spontan_new]).reset_index(drop=True)
#date = str(datetime.date.today())
#spontan_df_orig.to_excel(results_dir + 'spontan_all_' +  date + '.xlsx')

visual_QC_spontan = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/notes/analysis/events_detection/spontan_rough_evaluation_algorithm_all_results_01.02.2023.xlsx'
QC = pd.read_excel(visual_QC_spontan)

mask = QC['comment_elaboration'] == 'shifted events peaks'
QC['comment'].loc[mask] = None

spontan_df = event_funcs.post_event_analysis_main(QC, spontan_df_orig, 3) #min_event_size = 3, min_hrs = 20h
#spontan_8mM = event_funcs.remove_high_K_15mM(spontan_df)

n_nums_all = event_funcs.get_n_nums_per_day(spontan_df)
event_funcs.plot_event_by_day_spontan(spontan_df, n_nums_all, 'all')

# n_nums_8mM_high_K = event_funcs.get_n_nums_per_day(spontan_8mM)
# event_funcs.plot_event_by_day_spontan(spontan_8mM, n_nums_8mM_high_K , '8mM_high_K')

repatch_spontan = event_funcs.get_repatched_cells(spontan_df)
n_nums_repatch = event_funcs.get_n_nums_per_day(repatch_spontan)
event_funcs.plot_event_by_day_spontan(repatch_spontan, n_nums_repatch, 'repatch')

# repatch_spontan_8mM = event_funcs.get_repatched_cells(spontan_8mM )
# n_nums_repatch_8 = event_funcs.get_n_nums_per_day(repatch_spontan_8mM)
# event_funcs.plot_event_by_day_spontan(repatch_spontan_8mM, n_nums_repatch_8, 'repatch_8mM')


#%%

results_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/results/'

date = str(datetime.date.today())

#mini_df_orig = pd.concat([mini_old.loc[:], mini_new]).reset_index(drop=True)
#mini_df_orig.to_excel(results_dir + 'mini_all_' + date + '.xlsx')

visual_QC_minis = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/notes/analysis/events_detection/minis_rough_evaluation_algorithm_all_results_13.02.2023.xlsx'
QC = pd.read_excel(visual_QC_minis)
mask = QC['comment_elaboration'] == 'shifted events peaks'
QC['comment'].loc[mask] = None

mini_df = event_funcs.post_event_analysis_main(QC, mini_df_orig, min_event_size = 1, min_hrs = 20) #
#mini_8mM = event_funcs.remove_high_K_15mM(mini_df)

n_nums_minis = event_funcs.get_n_nums_per_day(mini_df)
event_funcs.plot_event_by_day_mini(mini_df, n_nums_minis)


#%%

fig,axarr = plt.subplots(1,1, figsize=(8,8))
fig.patch.set_facecolor('white')
for tr in results_df['treatment'].unique().tolist():
    plot_data = results_df[results_df['treatment'] == tr]

    count, bins_count = np.histogram(plot_data['Average amplitude (pA)'], bins = 15)
    pdf = count / sum(count) #probabily distribution function
    cdf = np.cumsum(pdf) #cumulative distribution function
    plt.plot(bins_count[1:], cdf, label = tr)
plt.figlegend()

labels = results_df['treatment'].unique().tolist()
fig,axarr = plt.subplots(1,1, figsize=(8,8))
fig.patch.set_facecolor('white')
for tr in results_df['treatment'].unique().tolist():
    plot_data = results_df[results_df['treatment'] == tr]
    sorted_plot_data = plot_data.sort_values('Average amplitude (pA)')

    sorted_plot_data['cum_sum'] = sorted_plot_data['Average amplitude (pA)'].cumsum()
    #calculate cumulative percentage of column (rounded to 2 decimal places)
    sorted_plot_data['cum_percent'] = round(100 * sorted_plot_data.cum_sum / sorted_plot_data['Average amplitude (pA)'].sum(),2)
    plt.plot(sorted_plot_data['Average amplitude (pA)'], sorted_plot_data['cum_percent'])

    #plt.plot(base[:-1], len(data)-cumulative, c='green') #survival function

#%%

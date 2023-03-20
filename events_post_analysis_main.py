import pandas as pd
import events_post_analysis_funcs as event_funcs


#%%
results_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/results/'

spontan_1 = pd.read_excel(results_dir + 'spontan_2023-02-07_summary_results_pt1.xlsx', 'summary')
spontan_2 = pd.read_excel(results_dir + 'spontan_2023-02-07_summary_results_pt2.xlsx', 'summary')
spontan_df_orig = pd.concat([spontan_1.loc[:], spontan_2]).reset_index(drop=True)

visual_QC_spontan = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/notes/analysis/events_detection/spontan_rough_evaluation_algorithm_all_results_01.02.2023.xlsx'
QC = pd.read_excel(visual_QC_spontan)

spontan_df = event_funcs.post_event_analysis_main(QC, spontan_df_orig) #min_event_size = 3, min_hrs = 24

repatch_df = event_funcs.get_repatched_cells(spontan_df)
n_nums_repatch = event_funcs.get_n_nums_per_day(repatch_df)
event_funcs.plot_event_by_day_spontan(repatch_df, n_nums_repatch, 'spontan')

#%%

results_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/results/'

mini_1 = pd.read_excel( results_dir + 'minis_2023-02-14_summary_results_pt1.xlsx', 'summary')
mini_2 = pd.read_excel( results_dir + 'minis_2023-02-14_summary_results_pt2.xlsx', 'summary')
mini_df_orig = pd.concat([mini_1.loc[:], mini_2]).reset_index(drop=True)

visual_QC_minis = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/notes/analysis/events_detection/minis_rough_evaluation_algorithm_all_results_13.02.2023.xlsx'
QC = pd.read_excel(visual_QC_minis)

mini_df = event_funcs.post_event_analysis_main(QC, mini_df_orig, min_event_size = 3, min_hrs = 22)

n_nums_minis = event_funcs.get_n_nums_per_day(mini_df)
event_funcs.plot_event_by_day_mini(mini_df, n_nums_minis, 'minis')



#%%


#plot

# #remove outliers
# for i in reversed(range(len(results_df))):       
#     if results_df['Average interevent interval (ms)'][i] > 1000:
#         results_df = results_df.drop([results_df.index[i]])

inc_time_D1 = results_df[results_df['Incubation time'] == 0]
inc_time_D2 = results_df[results_df['Incubation time'] >= 22]

boxplot = inc_time_D1.boxplot(column=['Average amplitude (pA)'], by=['treatment', ])
boxplot = inc_time_D2.boxplot(column=['Average amplitude (pA)'], by=['treatment'])


#remove outliers
for i in reversed(range(len(results_df))):       
    if results_df['Average interevent interval (ms)'][i] > 1000:
        results_df = results_df.drop([results_df.index[i]])

inc_time_D1 = results_df[results_df['Incubation time'] == 0]
inc_time_D2 = results_df[results_df['Incubation time'] >= 22]

boxplot = inc_time_D1.boxplot(column=['Average interevent interval (ms)'], by=['treatment'])
boxplot = inc_time_D2.boxplot(column=['Average interevent interval (ms)'], by=['treatment'])


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
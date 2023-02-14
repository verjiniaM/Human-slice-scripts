import pandas as pd
import events_post_analysis_funcs as event_funcs


#%%
results_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/results/'

spontan_1 = pd.read_excel(results_dir + 'spontan_2023-02-07_summary_results_pt1.xlsx', 'summary')
spontan_2 = pd.read_excel(results_dir + 'spontan_2023-02-07_summary_results_pt2.xlsx', 'summary')
spontan_df_orig = pd.concat([spontan_1.loc[:], spontan_2]).reset_index(drop=True)

visual_QC_spontan = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/notes/analysis/events_detection/rough_evaluation_algorithm_all_results_01.02.2023.xlsx'
QC = pd.read_excel(visual_QC_spontan)

QC = event_funcs.exclude_fn_only_in_QC(spontan_df_orig, QC)
spontan_df = event_funcs.get_non_excluded_traces(spontan_df_orig, QC)
spontan_df = event_funcs.add_day_of_recording_column(spontan_df)
#spontan_df = event_funcs.filter_min_hrs_incubation(spontan_df, 24)
spontan_df = event_funcs.remove_noisy_traces(spontan_df)
#spontan_df = event_funcs.remove_small_events(df = spontan_df, min_event_size = 3)

repatch_df = event_funcs.get_repatched_cells(spontan_df)
n_nums_repatch = event_funcs.get_n_nums_per_day(repatch_df)
event_funcs.    plot_event_by_day(repatch_df, n_nums_repatch, 'spontan')

#%%

results_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/results/'

mini_1 = pd.read_excel( results_dir + 'minis_2022-11-04_summary_results.xlsx', 'summary')
mini_2 = pd.read_excel( results_dir + 'minis_2022-11-08_summary_results.xlsx', 'summary')
mini_df = pd.concat([mini_1.loc[:], mini_2]).reset_index(drop=True)

for i in reversed(range(len(mini_df))):
    if mini_df['hrs_incubation'][i]  == 0:
        continue
    if float(mini_df['hrs_incubation'][i]) < 20 and float(mini_df['hrs_incubation'][i]) > 30:
        mini_df = mini_df.drop([mini_df.index[i]])
mini_df.reset_index(inplace = True, drop = True)

for i in reversed(range(len(mini_df))):
    if float(mini_df['Average interevent interval (ms)'][i]) > 1000 :
        mini_df = mini_df.drop([mini_df.index[i]])
        continue
mini_df.reset_index(inplace = True, drop = True)

for i in reversed(range(len(mini_df))):
    if float(mini_df['Average amplitude (pA)'][i]) > 30 :
        mini_df = mini_df.drop([mini_df.index[i]])
        continue
mini_df.reset_index(inplace = True, drop = True)

mini_df.boxplot(column=['Average interevent interval (ms)'], by=['treatment'])
mini_df.boxplot(column=['Average amplitude (pA)'], by=['treatment'])

treatments = ['Ctrl', 'high K', 'TTX']
nums_dict = {}
for tr in treatments:
    nums_dict[tr] = len(mini_df[(mini_df['treatment'] == tr)])


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

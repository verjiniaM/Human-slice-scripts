import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import os
import glob
import shutil

def copy_event_analysis_data_to_analysis_folder(event_type):
    '''
    moves current metadata.xlsx and the results.xlsx from spontaneous analysis folders to .../human/meta_events/
    loads results.xlsx - summary of event analysis
    '''
    results = '/Users/verjim/spontaneous-postsynaptic-currents-detection/results/results.xlsx'
    metadata = '/Users/verjim/spontaneous-postsynaptic-currents-detection/metadata/metadata.xlsx'
    date = str(datetime.date.today())

    human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
    shutil.copy(results, human_dir + '/meta_events/results/output_algirithm/' + event_type + '_' + date + '_results.xlsx')
    shutil.copy(metadata, human_dir + '/meta_events/analyzed_metadata/' + event_type + '_' + date + '_metadata.xlsx')

def post_events_analysis_add_metadata(event_type, human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/',
    results = '/Users/verjim/spontaneous-postsynaptic-currents-detection/results/results.xlsx',
    metadata = '/Users/verjim/spontaneous-postsynaptic-currents-detection/metadata/metadata.xlsx'):
    date = str(datetime.date.today())

    results_df = pd.read_excel(results, 'Summary results')
    meta_events = pd.read_excel(metadata)

    treatments, patchers, slices, OPs, incubation_times, cell_IDs = [], [], [], [], [], []
    for i in range(len(results_df['Recording filename'])):
        fn = results_df['Recording filename'][i]
        chan = results_df['Channel'][i]

        tr = meta_events['treatment'][meta_events['Name of recording'] == fn]
        patcher = meta_events['patcher'][meta_events['Name of recording'] == fn]
        slic = meta_events['slice'][meta_events['Name of recording'] == fn]
        hrs = meta_events['hrs_incubation'][meta_events['Name of recording'] == fn]
        OP = meta_events['OP'][meta_events['Name of recording'] == fn]

        cell_ID = meta_events['cell_ID'][(meta_events['Name of recording'] == fn) & \
                                        (meta_events['Channels to use'] == chan)].tolist()[0]

        treatments.append(tr.tolist()[0])
        patchers.append(patcher.tolist()[0])
        slices.append(slic.tolist()[0])
        incubation_times.append(hrs.tolist()[0])
        OPs.append(OP.tolist()[0])
        cell_IDs.append(cell_ID)

    results_df.insert(2, 'treatment', treatments)
    results_df.insert(3, ' patcher', patchers)
    results_df.insert(4, ' slice', slices)
    results_df.insert(5, 'hrs_incubation', incubation_times)
    results_df.insert(6, 'OP', OPs)
    results_df.insert(7, 'cell_ID', cell_IDs)

    #remove files where less than 2 min analysis
    for i in reversed(range(len(results_df))):
        if results_df['# analyzd sweeps'][i] < 12:
            results_df = results_df.drop([results_df.index[i]])
    results_df.reset_index(inplace = True, drop = True)


    with pd.ExcelWriter(human_dir + '/meta_events/results/' + event_type + '_' + date + '_summary_results.xlsx') as writer:
        results_df.to_excel(writer, sheet_name = "summary", index=False)

    return results_df

 
def get_events_numbers(results_df_long):
    # human_dir + '/meta_events/results/output_algirithm/' + date + '_results.xlsx'
    num_events = {}
    for tr in results_df['treatment'].unique().tolist():
        fns = results_df['Recording filename'][results_df['treatment'] == tr].tolist()
        chans = results_df['Channel'][results_df['treatment'] == tr].tolist()
        number_of_events = []
        for i, fn in enumerate(fns):
            chan = chans[i]
            sheet_name = fn[:-4] + '_' + str(chan)
            fn_df = pd.read_excel(results_path, sheet_name)
            number_of_events.append(len(fn_df))
        num_events[tr] = number_of_events
    return num_events

def remove_small_events(df, min_event_size):
    mask = (df['Average amplitude (pA)'] > min_event_size)
    df = df.loc[mask, :]
    df.reset_index(inplace = True, drop = True)
    return df

def add_day_of_recording_column(df):
    days = []
    for i in range(len(df)):
        if len(df[' slice'][i]) > 3:
            days.append('D2')
        else:
            days.append('D1')
    df['day'] = days
    return df

def filter_min_hrs_incubation(df, min_hrs):
    for i in reversed(range(len(df))):
        if df['hrs_incubation'][i]  == 0:
            continue
        if float(df['hrs_incubation'][i]) < min_hrs:
            df = df.drop([df.index[i]])
    df.reset_index(inplace = True, drop = True)
    return df

def exclude_fn_only_in_QC(spontan_df, QC):
    #exclude fn only in QC
    exclude_list = list(set(QC['recording'].tolist()) - set(spontan_df['Recording filename'].tolist()))

    for j in exclude_list:
        indx = QC[QC['recording'] == j].index
        QC.drop(indx, axis=0, inplace=True)
    QC.reset_index(inplace = True, drop = True)
    return QC

def get_non_excluded_traces(spontan_df, QC):
    for i in range(len(QC)):
        if QC.comment[i] == 'exclude':
            fn = QC.recording[i]
            chan = QC.channel[i]
            indx = spontan_df[(spontan_df['Recording filename'] == fn) & \
                             (spontan_df['Channel'] == chan)].index
            spontan_df.drop(indx, axis=0, inplace=True)
    spontan_df.reset_index(inplace = True, drop = True)
    return spontan_df

def remove_noisy_traces(df):
    '''
    removes traces with more noise than 1.5*SDs of the average noise
    plots a histogram for visualization
    '''
    avg_noise = np.mean(df['Stdev of the baseline signal (pA)'])
    sd_noise = np.std(df['Stdev of the baseline signal (pA)'])

    plt.hist(df['Stdev of the baseline signal (pA)'])
    plt.axvline(avg_noise + 1.5*sd_noise, color = 'r', label = '1.5 * SD noise')
    plt.legend()
    plt.show()

    mask = (df['Stdev of the baseline signal (pA)'] < (avg_noise + 1.5*sd_noise))
    df = df.loc[mask, :]
    df.reset_index(inplace = True, drop = True)
    return df

def get_n_nums_per_day(spontan_df):
    treatments = spontan_df['treatment'].unique()
    days = spontan_df['day'].unique() 
    
    n_nums_dict = {}
    for tr in treatments:
        for day in days:
            n_nums_dict[day + ' ' + tr] = len(spontan_df[(spontan_df['day'] == day) & (spontan_df['treatment'] == tr)])
    return n_nums_dict

def get_repatched_cells(spontan_df):
    '''
    Counts if the cell_ID is repeated
    adds counts column in the dataframe 
    '''
    count_column = []
    for i in range(len(spontan_df)):
        cell = spontan_df['cell_ID'][i]
        count_cell = spontan_df['cell_ID'].tolist().count(cell)
        count_column.append(count_cell)

    spontan_df['cell_count'] = count_column
    mask = (spontan_df['cell_count'] > 1)
    spontan_df = spontan_df.loc[mask, :]
    spontan_df.reset_index(inplace = True, drop = True)

    return spontan_df

def plot_event_by_day(spontan_df, n_nums, event_type, 
save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/event_analysis'):
    '''
    Takes the n_nums dictionary for each condition [day x treatemtn]
    scatter plot of the datapoints and median
    '''
    colors = ['moccasin', 'red', 'moccasin', 'cadetblue', 'moccasin', 'mediumpurple']
    fig = plt.figure(figsize=(10,5))
    ax = plt.subplot(1,1,1)
    day_label = []

    for i, comb in enumerate(n_nums.keys()):
        k = 0 + i

        day = comb[:2]
        treatment = comb[3:]
        plot_df = spontan_df[(spontan_df['treatment'] == treatment) & (spontan_df['day'] == day)]
        x = np.linspace(0.65+k, 1.35+k, len(plot_df))
        y = plot_df['Average amplitude (pA)']
        median = np.median(y)
        ax.scatter(x, y, alpha = 0.8, c = colors[k], s = 40)
        ax.plot([0.8+k, 1.2+k], [median, median], c = 'k', linestyle = 'solid', linewidth = 2)
        day_label.append(comb)
        ax.text(k+0.85, median + 0.5, str(round(median, 2)), size = 10)
        ax.text(k+0.85, int(np.max(spontan_df['Average amplitude (pA)'])+1), 'n = ' + str(n_nums[comb]), size = 10)


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks(ticks = list(range(1,7)), labels = day_label, size = 12)
    ax.set_yticks(range(int(np.min(spontan_df['Average amplitude (pA)'])),
                        int(np.max(spontan_df['Average amplitude (pA)'])+4),5), size = 12) 
    plt.title('Average amplitude (pA)', fontsize = 19)
    ax.set_xlabel('Day', fontsize = 15)
    ax.set_ylabel('pA', fontsize = 15)

    fig.patch.set_facecolor('white')
    plt.savefig(save_dir  + event_type + 'avg_amplitude_not_excluded.png')
    plt.show(fig)

def post_event_analysis_main(event_type):
    copy_event_analysis_data_to_analysis_folder(event_type)
    results_df = post_events_analysis_add_metadata(event_type)

#%%

#plot

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

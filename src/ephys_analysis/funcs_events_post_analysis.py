import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import shutil
import ephys_analysis.funcs_human_characterisation as hcf
import ephys_analysis.funcs_plot_intrinsic_props  as pl_intr

def copy_event_analysis_data_to_analysis_folder(event_type):
    '''
    moves current metadata.xlsx and the results.xlsx from spontaneous analysis folders to .../human/meta_events/
    loads results.xlsx - summary of event analysis
    '''
    results = '/Users/verjim/spontaneous-postsynaptic-currents-detection/results/results.xlsx'
    metadata = '/Users/verjim/spontaneous-postsynaptic-currents-detection/metadata/metadata.xlsx'
    date = str(datetime.date.today())

    human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
    shutil.copy(results, human_dir + 'meta_events/results/output_algirithm/' + event_type + '_' + date + '_results.xlsx')
    shutil.copy(metadata, human_dir + 'meta_events/analyzed_metadata/' + event_type + '_' + date + '_metadata.xlsx')

def post_events_analysis_add_metadata(event_type, human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/',
    results = '/Users/verjim/spontaneous-postsynaptic-currents-detection/results/results.xlsx',
    metadata = '/Users/verjim/spontaneous-postsynaptic-currents-detection/metadata/metadata.xlsx'):
    
    ''' 
    Takes the overview of the results and the metadataframe from the last analysis in the events direcotry
    Adds treatment, patcher, slice, incubation time and OP info
    removes files with less than 2 min of analysis
    saves and returns a results dataframe
    '''

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

        #cell_ID = meta_events['cell_ID'][(meta_events['Name of recording'] == fn) & \
        #                               (meta_events['Channels to use'] == chan)].tolist()[0]

        treatments.append(tr.tolist()[0])
        patchers.append(patcher.tolist()[0])
        slices.append(slic.tolist()[0])
        incubation_times.append(hrs.tolist()[0])
        OPs.append(OP.tolist()[0])
        #cell_IDs.append(cell_ID)

    results_df.insert(2, 'treatment', treatments)
    results_df.insert(3, ' patcher', patchers)
    results_df.insert(4, ' slice', slices)
    results_df.insert(5, 'hrs_incubation', incubation_times)
    results_df.insert(6, 'OP', OPs)
    #results_df.insert(7, 'cell_ID', cell_IDs)

    #remove files where less than 2 min analysis
    for i in reversed(range(len(results_df))):
        if results_df['# analyzed sweeps'][i] < 12:
            results_df = results_df.drop([results_df.index[i]])
    results_df.reset_index(inplace = True, drop = True)

    results_df = add_num_events_col(results_df, results)

    with pd.ExcelWriter(human_dir + '/meta_events/results/summaries/' + event_type + '_' + date + '_summary_results.xlsx') as writer:
        results_df.to_excel(writer, sheet_name = "summary", index=False)

    return results_df

 
def add_num_events_col(results_df, results_long_path):
    ''' 
    number of events = len of the sheet in results_long_path
    results_df - summary results file in human directory
    results_long_path - full results table in spontan_analysis folder
    '''

    number_of_events = []
  
    for i in range(len(results_df)):
        dat = results_df.loc[i,:]
        fn = dat['Recording filename']
        chan = dat['Channel']
        sheet_name = fn[:-4] + '_' + str(chan)
        fn_df = pd.read_excel(results_long_path, sheet_name)
        number_of_events.append(len(fn_df))

    results_df.insert(16, 'num_events', number_of_events)
    return results_df


## 
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

def exclude_fn_only_in_QC(df, QC):
    #exclude fn only in QC
    exclude_list = list(set(QC['recording'].tolist()) - set(df['Recording filename'].tolist()))

    for j in exclude_list:
        indx = QC[QC['recording'] == j].index
        QC.drop(indx, axis=0, inplace=True)
    QC.reset_index(inplace = True, drop = True)
    return QC

def get_non_excluded_traces(df, QC):

    for i in range(len(QC)):
        if QC.comment[i] == 'exclude':
            fn = QC.recording[i]
            chan = QC.channel[i]
            indx = df[(df['Recording filename'] == fn) & \
                             (df['Channel'] == chan)].index
            df.drop(indx, axis=0, inplace=True)
    df.reset_index(inplace = True, drop = True)
    return df

def remove_noisy_traces(df):
    '''
    removes traces with more noise than 1.5*SDs of the average noise
    plots a histogram for visualization
    '''
    avg_noise = np.mean(df['Stdev of the baseline signal (pA)'])
    sd_noise = np.std(df['Stdev of the baseline signal (pA)'])

    plt.hist(df['Stdev of the baseline signal (pA)'])
    plt.axvline(avg_noise + 1.5*sd_noise, color = 'r', label = '1.5 * SD noise')
    plt.xlabel('SD of baseline signal')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()

    mask = (df['Stdev of the baseline signal (pA)'] < (avg_noise + 1.5*sd_noise))
    df = df.loc[mask, :]
    df.reset_index(inplace = True, drop = True)
    return df

def get_n_nums_per_day(df, treatments = ['Ctrl', 'high K']):
    
    days = df['day'].unique() 
    
    n_nums_dict = {}
    for tr in treatments:
        for day in days:
            n_nums_dict[day + ' ' + tr] = len(df[(df['day'] == day) & (df['treatment'] == tr)])
    return n_nums_dict

def get_repatched_cells(df):
    '''
    Counts if the cell_ID is repeated
    adds counts column in the dataframe 
    '''
    count_column = []
    for i in range(len(df)):
        cell = df['cell_ID'][i]
        count_cell = df['cell_ID'].tolist().count(cell)
        count_column.append(count_cell)

    df['cell_count'] = count_column
    mask = (df['cell_count'] > 1)
    df = df.loc[mask, :]
    df.reset_index(inplace = True, drop = True)

    return df

def plot_event_by_day_spontan(df, n_nums, data_type, 
save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/event_analysis/'):
    '''
    Takes the n_nums dictionary for each condition [day x treatemtn]
    scatter plot of the datapoints and median
    '''

    colors = ['#dede00', '#ff7f00', '#dede00', '#DA67BA', '#dede00', '#7FEE8F'] #Ctrl, TTX, 8mM Kcl, 15 mM Kcl

    title1 = 'Spontaneous EPSCs'
    cmap = plt.cm.get_cmap('tab20')
    op_color_dict = {}
    for h, op in enumerate(df['OP'].unique()):
        op_color_dict[op] = cmap((h+1)/10)
    
    fig = plt.figure(figsize=(10,5))
    ax = plt.subplot(1,1,1)
    day_label = []

    min_h_incubation = str(min(df['hrs_incubation']))
    min_amplitude = str(min(df['Average amplitude (pA)']))[:4]
    x_plot = []
    for i, comb in enumerate(n_nums.keys()):
        k = 0 + i
        day = comb[:2]
        treatment = comb[3:]
        plot_df = df[(df['treatment'] == treatment) & (df['day'] == day)]
        x = np.linspace(0.65+k, 1.35+k, len(plot_df))
        x_plot.append(x)
        y = plot_df['Average amplitude (pA)']
        median = np.median(y)
        ax.scatter(x, y, alpha = 0.7, c = colors[k], s = 40)
        ax.plot([0.8+k, 1.2+k], [median, median], c = 'k', linestyle = 'solid', linewidth = 2)
        day_label.append(comb)
        ax.text(k+0.85, median + 0.5, str(round(median, 2)), size = 15)
        ax.text(k+0.85, int(np.max(df['Average amplitude (pA)'])+1), 'n = ' + str(n_nums[comb]), size = 12)
        if k in [1,3,5] and data_type == 'repatch':
            for c, cell in enumerate(plot_df['cell_ID']):
                x1 = [x_plot[0][c], x[c]] 
                y = df['Average amplitude (pA)'][df['cell_ID'] == cell]
                op = plot_df['OP'][plot_df['cell_ID'] == cell].tolist()[0]
                plt.plot(x1, y, '-', color = op_color_dict[op], alpha = 0.6, label = op, zorder = 0)
            title1 = 'Spontaneous EPSCs (repatch)'
            x_plot = []

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks(ticks = list(range(1,5)), labels = day_label, size = 15)
    ax.set_yticks(range(int(np.min(df['Average amplitude (pA)'])),
                        int(np.max(df['Average amplitude (pA)'])+4),5), size = 15) 
    plt.title(title1, fontsize = 19, x = 0.5, y = 1)
    ax.set_xlabel('Condition', fontsize = 15)
    ax.set_ylabel('Average amplitude (pA)', fontsize = 15)
    fig.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    #plt.figlegend(loc = 'upper right',  bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.patch.set_facecolor('white')
    #plt.show(fig)
    plt.savefig(save_dir  + data_type + '_spontan_min_amplitude_' + min_amplitude + '_min_hrs_incubation_' + min_h_incubation + '_avg_amplitude_not_excluded.pdf')
    plt.close(fig)

def plot_event_by_day_mini(df, n_nums, save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/event_analysis/'):
    '''
    Takes the n_nums dictionary for each condition 
    scatter plot of the datapoints and median
    '''
    #colors = ['red', 'cadetblue', 'mediumpurple']
    colors = ['#ff7f00',  '#DA67BA', '#7FEE8F', '#186F25'] #Ctrl, TTX, 8mM Kcl, 15 mM Kcl

    fig = plt.figure(figsize=(10,5))
    ax = plt.subplot(1,1,1)
    tr_label = []

    min_h_incubation = str(min(df['hrs_incubation']))
    min_amplitude = str(min(df['Average amplitude (pA)']))[:4]

    for i, comb in enumerate(n_nums.keys()):
        treatment = comb[3:]
        tr_label.append(treatment)
        plot_df = df[(df['treatment'] == treatment)]
        x = np.linspace(0.65+i, 1.35+i, len(plot_df))
        y = plot_df['Average amplitude (pA)']
        median = np.median(y)
        ax.scatter(x, y, alpha = 0.8, c = colors[i], s = 40)
        ax.plot([0.8+i, 1.2+i], [median, median], c = 'k', linestyle = 'solid', linewidth = 2)
        ax.text(i+0.85, median + 0.5, str(round(median, 2)), size = 15)
        ax.text(i+0.85, int(np.max(df['Average amplitude (pA)'])+1), 'n = ' + str(n_nums[comb]), size = 15)


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks(ticks = list(range(1,4)), labels = tr_label, size = 15)
    ax.set_yticks(range(int(np.min(df['Average amplitude (pA)'])),
                        int(np.max(df['Average amplitude (pA)'])+4),5), size = 15) 
    plt.title('mEPSCs in 1 µm TTX, 1 µm Gabazine and 50 µm APV', fontsize = 19, x = 0.5, y = 1.2)
    ax.set_xlabel('Condition', fontsize = 17)
    ax.set_ylabel('Average amplitude (pA)', fontsize = 17)

    fig.patch.set_facecolor('white')
    fig.tight_layout()
    plt.savefig(save_dir + 'minis_min_amplitude_' + min_amplitude + '_min_hrs_incubation_' + min_h_incubation + '_avg_amplitude_not_excluded.png')
    plt.show(fig)


def post_event_analysis_main(QC, df_orig, min_event_size, min_hrs = 20):
    ''' 
    
    '''
    QC_analyzed = exclude_fn_only_in_QC(df_orig, QC)
    df_analysis = get_non_excluded_traces(df_orig, QC_analyzed)
    df_analysis = add_day_of_recording_column(df_analysis)
    df_analysis = filter_min_hrs_incubation(df_analysis, min_hrs)
    df_analysis = remove_noisy_traces(df_analysis)
    df_analysis = remove_small_events(df = df_analysis, min_event_size = min_event_size)
    return df_analysis

def get_meta_data(results_df, meta_df):
    '''
    gets in
    '''

    event_types,cell_IDs, treatments, keep, exclusion_reasons = [],[],[], [], []
    for i, fn in enumerate(results_df['Recording filename']):
        chan = results_df['Channel'][i]
        dat_to_add = meta_df[(meta_df['Name of recording'] == fn) & 
                    (meta_df['Channels to use'] == chan)]

        event_types.append(dat_to_add['event_type'].item())
        cell_IDs.append(dat_to_add['cell_ID'].item())
        #treatments.append(dat_to_add['treatment'].item())
        keep.append(dat_to_add['comment'].item())
        exclusion_reasons.append(dat_to_add['exclusion_reason'].item())

    list_of_lists = [event_types, cell_IDs, keep, exclusion_reasons]
    list_col_names = ['event_type', 'cell_ID', 'comment', 'exclusion_reason']

    #columns_to_include = list(set(list(RMP_high_K_all.columns)) - set(list(results_df.columns)))
    for i, col in enumerate(list_of_lists):
        results_df.insert(i, list_col_names[i], col)
    
    return results_df

def QC_post_event_analysis (df_analysis, meta_df, min_event_size, min_hrs = 16):
    ''' 
    QC for manually analyzed events
    '''
    df_analysis = get_meta_data(df_analysis, meta_df)
    df_analysis = add_day_of_recording_column(df_analysis)
    df_analysis = filter_min_hrs_incubation(df_analysis, min_hrs)
    #df_analysis = remove_noisy_traces(df_analysis) # removes traces with more noise than 1.5*SDs
    df_analysis = remove_small_events(df = df_analysis, min_event_size = min_event_size)

    return df_analysis

def add_high_K_concentration(spontan_df):
    exp_overview = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/2023-05-08_experiments_overview.xlsx')

    dict_higk_K_concentr = {}
    for OP in spontan_df_orig['OP'].unique():
        K_concentr = exp_overview['K concentration'][exp_overview['OP'] == OP].tolist()[0]
        dict_higk_K_concentr[OP] = K_concentr

    K_concentr_all = []
    OP_list = spontan_df['OP'].tolist()
    for i in range(len(spontan_df)):
        OP = OP_list[i]
        K_concentr_all.append(dict_higk_K_concentr[OP])

    spontan_df['high K concentration'] = K_concentr_all
    spontan_df['K_concentr_bool'] = K_concentr_all

    spontan_df = spontan_df.replace({'K_concentr_bool': '8 mM'}, 0)
    spontan_df = spontan_df.replace({'K_concentr_bool': '15 mM'}, 1)

    return spontan_df



def QC_spontan_with_intrinsic(intrinsic_df, spontan_df):
    cell_IDs, repatches, not_in_intrinsic = [], [], []
    for i in range(len(spontan_df)):
        OP = spontan_df['OP'][i]
        patcher = spontan_df[' patcher'][i]
        fn = spontan_df['Recording filename'][i]
        slic = spontan_df[' slice'][i]
        day = spontan_df['day'][i]
        cell_ch = spontan_df['Channel'][i]

        cell_ID = intrinsic_df['cell_ID_new'][(intrinsic_df['OP'] == OP) & (intrinsic_df['patcher'] == patcher) & \
        (intrinsic_df['slice'] == slic) & (intrinsic_df['day'] == day) & \
        (intrinsic_df['cell_ch'] == cell_ch)]

        if len(cell_ID) == 0:
            cell_ID = hcf.get_new_cell_IDs_fn(fn, slic, [cell_ch], patcher)[0]
            
            not_in_intrinsic.append(cell_ID)
            cell_IDs.append(cell_ID)
            repatches.append('no')
            continue

        cell_IDs.append(cell_ID.tolist()[0])
        repatches.append(intrinsic_df['repatch'][intrinsic_df['cell_ID_new'] == cell_ID.tolist()[0]].tolist()[0])


    spontan_df.insert(23, 'cell_ID_new', cell_IDs)
    spontan_df.insert(24, 'repatch', repatches)

    #remove cells which are not in the intrinsic_df (not QC passed)
    not_QC_passed_cells = pl_intr.in_list1_not_in_list2(spontan_df['cell_ID_new'].tolist(), intrinsic_df['cell_ID_new'].tolist())
    for cell in not_QC_passed_cells :
        spontan_df = spontan_df.drop(spontan_df.index[spontan_df['cell_ID_new'] == cell])
    spontan_df.reset_index(inplace = True, drop = True)

    spontan_repatch = pl_intr.get_repatch_df(spontan_df)
    return spontan_df, spontan_repatch

# def remove_high_K_15mM (df):
#     indx1 = sorted(df['OP'].unique()).index('OP220602')
#     exclude_list = sorted(df['OP'].unique())[indx1:]
#     for j in exclude_list:
#         indx = df[df['OP'] == j].index
#         df.drop(indx, axis=0, inplace=True)
#     df.reset_index(inplace = True, drop = True)
#     return df

#%%
#%%

def plot_event_by_day_mini_2_IEI(df, n_nums,
save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/event_analysis/lab_seminar_september_2023/'):
    colors = ['#ff7f00', '#DA67BA','#186F25','#7FEE8F'] #Ctrl, TTX, 8mM Kcl,15 mM Kcl

    title1 = 'Mini EPSCs in slice'

    cmap = plt.cm.get_cmap('tab20')
    op_color_dict = {}
    for h, op in enumerate(df['OP'].unique()):
        op_color_dict[op] = cmap((h+1)/10)

    fig = plt.figure(figsize=(7,7))
    ax = plt.subplot(1,1,1)
    tr_label = []

    min_h_incubation = str(min(df['hrs_incubation']))
    min_amplitude = str(min(df['Average amplitude (pA)']))[:4]
    day_label = []
    x_plot = []
    for i, comb in enumerate(n_nums.keys()):
        k = 0 + i

        treatment = comb[3:]
        plot_df = df[df['treatment'] == treatment]
        x = np.linspace(0.8+k, 1.35+k, len(plot_df))
        x_plot.append(x)
        y = plot_df['Average interevent interval (ms)']
        median = np.median(y)
        if treatment == 'high K':   
            y_8mm = plot_df['Average interevent interval (ms)'].loc[plot_df['high K concentration'] == '8 mM'] 
            x_8 = np.linspace(0.8+k, 1.05+k, len(y_8mm))
            median_8 = np.median(y_8mm)
            y_15mm = plot_df['Average interevent interval (ms)'].loc[plot_df['high K concentration'] == '15 mM'] 
            x_15 = np.linspace(1.10+k, 1.35+k, len(y_15mm))
            median_15 = np.median(y_15mm)
            ax.scatter(x_8, y_8mm, alpha = 0.9, c = colors[2], s = 80)
            ax.scatter(x_15, y_15mm, alpha = 0.9, c = colors[3], s = 80)
            ax.boxplot(y_8mm, positions = [k+0.65], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[2], alpha = 0.75),
            medianprops = dict(linewidth=2.3, color = 'k'))
            ax.boxplot(y_15mm, positions = [k+1.55], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[3], alpha = 0.75),
            medianprops = dict(linewidth=2.3, color = 'k'))
            # ax.text(k+0.65, int(np.max(df['Average interevent interval (ms)'])+1), 'n = ' + str(len(y_8mm)), size = 15, color = colors[1])
            # ax.text(k+1.15, int(np.max(df['Average interevent interval (ms)'])+1), 'n = ' + str(len(y_15mm)), size = 15, color = colors[2])

        else:
            ax.boxplot(y, positions = [k + 0.7], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[k], alpha = 0.75),
            medianprops = dict(linewidth=2.3, color = 'k'))
            y = plot_df['Average interevent interval (ms)']
            ax.scatter(x, y, alpha = 0.9, c = colors[k], s = 80)
            #ax.plot([0.9+k, 1.1+k], [median, median], c = 'k', linestyle = 'solid', linewidth = 1)
            #ax.text(k+0.85, median + 0.5, str(round(median, 2)), size = 15)
            day_label.append(comb)
            #ax.text(k+0.85, int(np.max(df['Average interevent interval (ms)'])+1), 'n = ' + str(n_nums[comb]), size = 15,color = colors[0])

    #ax.boxplot([[1,4],[2,5],[3,6]], positions=[2,4,5.5])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks(ticks = [0.7, 1.7, 2.7, 3.55], labels = ['Ctrl aCSF \n incub.', 
    'aCSF + \n 1 um TTX', 'aCSF + \n 8 mM KCL', 'aCSF + \n 15 mM KCL'], size = 15)
    plt.title(title1, fontsize = 19, x = 0.5, y = 1)
    #ax.set_xlabel('Condition', fontsize = 15)
    ax.set_ylabel('Average interevent interval (ms)', fontsize = 24)
    #ax.set_yticks([5,10, 15,20,25,30])
    ax.tick_params(axis='y', labelsize=22)

    plt.figlegend(loc = 'upper right',  bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.patch.set_facecolor('white')
    plt.savefig(save_dir + 'minis_IEI_plot.pdf')




def plot_event_by_day_mini_2_amp(df, n_nums,
save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/event_analysis/lab_seminar_september_2023/'):
    colors = ['#ff7f00', '#DA67BA','#186F25','#7FEE8F'] #Ctrl, TTX, 8mM Kcl,15 mM Kcl

    title1 = 'Mini EPSCs in slice'

    cmap = plt.cm.get_cmap('tab20')
    op_color_dict = {}
    for h, op in enumerate(df['OP'].unique()):
        op_color_dict[op] = cmap((h+1)/10)

    fig = plt.figure(figsize=(7,7))
    ax = plt.subplot(1,1,1)
    tr_label = []

    min_h_incubation = str(min(df['hrs_incubation']))
    min_amplitude = str(min(df['Average amplitude (pA)']))[:4]
    x_plot = []
    for i, comb in enumerate(n_nums.keys()):
        k = 0 + i

        treatment = comb[3:]
        plot_df = df[df['treatment'] == treatment]
        x = np.linspace(0.8+k, 1.35+k, len(plot_df))
        x_plot.append(x)
        y = plot_df['Average amplitude (pA)']
        median = np.median(y)
        if treatment == 'high K':   
            y_8mm = plot_df['Average amplitude (pA)'].loc[plot_df['high K concentration'] == '8 mM'] 
            x_8 = np.linspace(0.8+k, 1.05+k, len(y_8mm))
            median_8 = np.median(y_8mm)
            y_15mm = plot_df['Average amplitude (pA)'].loc[plot_df['high K concentration'] == '15 mM'] 
            x_15 = np.linspace(1.10+k, 1.35+k, len(y_15mm))
            median_15 = np.median(y_15mm)
            ax.scatter(x_8, y_8mm, alpha = 0.9, c = colors[2], s = 80)
            ax.scatter(x_15, y_15mm, alpha = 0.9, c = colors[3], s = 80)
            ax.boxplot(y_8mm, positions = [k+0.65], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[2], alpha = 0.75),
            medianprops = dict(linewidth=2.3, color = 'k'))
            ax.boxplot(y_15mm, positions = [k+1.55], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[3], alpha = 0.75),
            medianprops = dict(linewidth=2.3, color = 'k'))
            ax.text(k+0.65, int(np.max(df['Average amplitude (pA)'])+1), 'n = ' + str(len(y_8mm)), size = 15, color = colors[2])
            ax.text(k+1.15, int(np.max(df['Average amplitude (pA)'])+1), 'n = ' + str(len(y_15mm)), size = 15, color = colors[3])

        else:
            ax.boxplot(y, positions = [k + 0.7], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[k], alpha = 0.75),
            medianprops = dict(linewidth=2.3, color = 'k'))
            y = plot_df['Average amplitude (pA)']
            ax.scatter(x, y, alpha = 0.9, c = colors[k], s = 80)
            #ax.plot([0.9+k, 1.1+k], [median, median], c = 'k', linestyle = 'solid', linewidth = 1)
            #ax.text(k+0.85, median + 0.5, str(round(median, 2)), size = 15)
            day_label.append(comb)
            ax.text(k+0.85, int(np.max(df['Average amplitude (pA)'])+1), 'n = ' + str(n_nums[comb]), size = 15,color = colors[k])

    #ax.boxplot([[1,4],[2,5],[3,6]], positions=[2,4,5.5])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks(ticks = [0.7, 1.7, 2.7, 3.55], labels = ['Ctrl aCSF \n incub.', 
    'aCSF + \n 1 um TTX', 'aCSF + \n 8 mM KCL', 'aCSF + \n 15 mM KCL'], size = 15)
    #plt.title(title1, fontsize = 19, x = 0.5, y = 1)
    #ax.set_xlabel('Condition', fontsize = 15)
    ax.set_ylabel('Average amplitude (pA)', fontsize = 24)
    ax.set_yticks([5,10, 15,20,25,30])
    ax.tick_params(axis='y', labelsize=22)

    plt.figlegend(loc = 'upper right',  bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.patch.set_facecolor('white')
    plt.savefig(save_dir + '/mini_slice_avg_amp.pdf')


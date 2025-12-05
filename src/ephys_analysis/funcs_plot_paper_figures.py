import os
import pandas as pd
import numpy as np
import seaborn as sns
import datetime
import pyabf
import matplotlib.pyplot as plt
import ephys_analysis.funcs_con_screen as con_param
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
import ephys_analysis.funcs_human_characterisation as hcf
import ephys_analysis.funcs_sorting as sort
import ephys_analysis.stimulation_windows_ms as stim_win
from scipy import stats

plt.style.use('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/Human-slice-scripts/style_plot_paper.mplstyle')

#%%
def figure_colors():
    color_dict =  {'inc_only':{
                    'Ctrl': "#570b77",
                    'high K':  "#1bd7e0"},
                    'slice_all_CTR': {'CTR': "#0e62ff"},
                   'slice':
                    {'CTRD1':'darkgrey', #'#dede00', 'MediumPurple2'
                    'CTRD2': "#d383cf", # pink
                    'HiKD1': 'darkgrey', # '#dede00',
                    'HiKD2': "#3a9837",
                    'CtrlD1':'darkgrey', 
                    'CtrlD2': "#d383cf", 
                    'high KD1': 'darkgrey', 
                    'high KD2': "#3a9837",
                    'Ctrl': "#d383cf",
                    'HighK': "#3a9837",
                    'high K': "#3a9837",
                    'NI':'darkgrey', 
                    'pre': 'darkgrey',
                    'post HiK': "#3a9837",
                    'post CTR': "#d383cf"},
                    'cell_ID_new':
                    {'CTRD1': 'darkgrey',
                    'CTRD2': '#ff7f00', 
                    'HiKD1': 'darkgrey', 
                    'HiKD2': "#4A74DD",
                    'CtrlD1': 'darkgrey', 
                    'CtrlD2': '#ff7f00', 
                    'high KD1': 'darkgrey', # '#dede00',
                    'high KD2': "#4A74DD",
                    'pre': 'darkgrey',
                    'post HiK': "#4A74DD",
                    'post CTR': '#ff7f00'},
                    'all':
                    {'Ctrl': '#496B7F'}}
    return color_dict

def trace_colors_trace():
    color_dict =  {'slice':
                    {'D1.3':"#A9A3A3",
                     'D1.2': "#8E8D8D",
                     'D1.1': "#636262", 
                    'CTRD2.3': "#e4a1e1",
                    'CTRD2.2': "#d383cf",
                    'CTRD2.1': "#91518e",
                    'HiKD2.3': "#65db60",
                    'HiKD2.2': "#3a9837",
                    'HiKD2.1': "#1a7717"}, 
                    'cell_ID_new':
                    {'D1.3':"#A9A3A3",
                     'D1.2': "#908686",
                     'D1.1': "#2C2B2B",
                    'CTRD2.3': "#f65f08", 
                    'CTRD2.2': '#ff7f00',
                    'CTRD2.1': "#f1b62be6",
                    'HiKD2.1': "#3655A1",
                    'HiKD2.2': "#4A74DD",
                    'HiKD2.3': "#5C77EECA"}}
    return color_dict

def mea_fig1(df, dest_dir, w_cm = 6, h_cm = 5.7):

    dict_label = {'Ctrl':'CTR', 'high K':'HiK'}
    colors = figure_colors()['slice']

    fig, ax = plt.subplots(1, 1, figsize=(w_cm/2.54, h_cm/2.54))
    labels_ = []
    for i, treatment in enumerate(df.treatment.unique()):
        df_tr = df[df['treatment'] == treatment]
        x = np.linspace(0.7 + i, 1.3 + i, len(df_tr))
        ax.scatter(x, df_tr.value, c = colors[treatment], zorder = 5, alpha = 0.3)
        ax.scatter(x, df_tr.value, edgecolor = colors[treatment ],
                   facecolors = 'none', alpha = 0.75)
     
        ax.boxplot(df_tr.value, positions = [i + 1], zorder = 10)
        
        labels_.append(dict_label[treatment])

    ax.set_xticks(ticks = [1,2], labels = labels_)
    ax.set_ylabel('Normalized \n (Incub. Sol - BL)/BL')
    ax.set_yticks(ticks = np.linspace(0,20,5), labels = [0, '', 10, '', 20])

    date = str(datetime.date.today())
    plt.savefig(dest_dir  + date + '_plot_MEA.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.close(fig)

def fig_order():

    """
    Returns the order of figures for intr
    """
    order_dict = {'inc_only':{
                    'Rin': (0,1),
                    'resting_potential': (0,2),
                    'rheo_ramp_c' : (1,0),
                    'TH': (1,1),
                    'sag': (1,2)},
                  'other':{
                    'resting_potential': (1,2),
                    'rheo_ramp_c' : (2,0),
                    'Rin': (0,2),
                    'sag': (2,2),
                    'TH': (2,1)},
                  'slice_all_CTR':{
                      'Rin': (0,1),
                      'resting_potential': (0,2),
                      'TH': (1,0),
                      'rheo_ramp_c' : (1,1),
                      'sag': (1,2)}}

    return order_dict

def intr_params_fig2(df, dest_dir, w_cm = 20, h_cm = 15.6):
    """
    Plots parameters for different treatments and days 
    """
    # for single plots - w_cm = 6.5, h_cm = 5.7
    mask_exclude = (df['treatment'] != 'high K') & (df['treatment'] != 'Ctrl')
    df = df[~mask_exclude]
    params = ['Rin', 'resting_potential', 'rheo_ramp_c', 'TH', 'sag']
    ax_elements = fig_order()['other']
    df = df.sort_values(['treatment', 'day'])
    titles_all = pl_intr.dict_for_plotting()
    titles_dict = {key: titles_all[key] for key in params if key in titles_all}

    dict_label = {'Ctrl':'CTR', 'high K':'HiK'}
    colors = figure_colors()['slice']
    fig2, axs = plt.subplots(3,3, figsize = (w_cm/2.54, h_cm / 2.54))

    for param, item in titles_dict.items():
        ax_el = ax_elements[param]
        label_, data_boxplot = [], []
        num_cels = {}
        for i, tr in enumerate(df.treatment.unique()):
            for j, day in enumerate(df['day'].unique()):
                k = j + 2.5*i
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                df_plot = df_plot[~df_plot[param].isna()]  # remove NaN values
                x = np.linspace(0.65+k, 1.35+k, len(df_plot))

                axs[ax_el].scatter(x, df_plot[param], alpha = 0.3,
                            color = colors[tr + day], zorder = 5) # dots inside
                axs[ax_el].scatter(x, df_plot[param], alpha = 0.75,
                            edgecolor = colors[tr + day],
                            facecolors = 'none')

                axs[ax_el].boxplot(df_plot[param], positions = [k + 1], zorder = 10, widths = 0.5)

                label_.append(day + '\n' + dict_label[tr])
                num_cels[tr + ' ' + day] = len(df_plot)
                data_boxplot.append(df_plot[param])

        axs[ax_el].set_xticks(ticks = [1,2,3.5,4.5], labels = ['pre', 'post', 'pre', 'post'])

        axs[ax_el].set_ylabel(item[0] + ' (' + item[1] + ')')
        axs[ax_el].set_yticks(ticks = item[2], labels = item[3])
        # axs[ax_el].set_title('CTR                            HiK')

    date = str(datetime.date.today())
    os.makedirs(dest_dir, exist_ok = True)
    plt.savefig(dest_dir  + date + 'V2_plot_all' + '.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 1000)

def firing_props_fig2(iff_df, data_type, dv, error_type, dest_dir,
                 background_dot_size = 3, w_cm = 13, h_cm = 6):
    """
    Plots the initial firing frequency and number of APs
    """
    df_summary = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results'+ \
                           '/human/paper_figs_collected_checked/stats/archive/'+ \
                            'sum_data_CIs_SEs_for_plotting/firing_plot_CIs_temporal.xlsx')
    iff_df = iff_df.sort_values(['treatment', 'day'])
    mask_exclude = (iff_df['treatment'] != 'high K') & (iff_df['treatment'] != 'Ctrl')
    iff_df = iff_df[~mask_exclude]
    num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(iff_df)

    colors = figure_colors()['slice']
    dt_sum = 'slice'

    dv_dict = {'IFF' : [iff_indx, 'Initial firing frequency (AP#1 to AP#2)',
                'Initial firing \n frequency (Hz)'],
                'num_aps' : [num_aps_indx, 'Number of fired action potentials', 'AP count']}

    # fig, ax = plt.subplots(3, len(iff_df.treatment.unique()), figsize=(w_cm / 2.54, h_cm/2.54),\
    #                         sharex = True, gridspec_kw = {'height_ratios': [5, 1, 1]})
    fig, ax = plt.subplots(1, len(iff_df.treatment.unique()), figsize=(w_cm / 2.54, h_cm/2.54), sharex = True)

    counts_firing = {}
    for k, treatment in enumerate(iff_df.treatment.unique()):
        tr_df = iff_df[iff_df['treatment'] == treatment]

        for d, day in enumerate(tr_df['day'].unique()):
            day_df = tr_df[tr_df['day'] == day]
            print(treatment, day)
            avgs, inj = [], []

            counts = []
            for i, col in enumerate(dv_dict[dv][0][5:5+13]):
                data = day_df.iloc[:,col]
                inj_val = int(day_df.columns[col][2:day_df.columns[col].rfind('pA')])

                # removing non-firing cells
                data.replace(0, np.nan, inplace = True)
                data.dropna(inplace = True)
                avg = np.mean(data)

                if error_type == 'SE':
                    x = np.linspace(0.75 + i , 1.25 + i, len(data))
                    # ax[k].scatter(x, data, alpha = 0.15, s = background_dot_size,
                    #             c = colors[treatment + day])
                    sum_data = df_summary[(df_summary['data_type'] == dt_sum) & \
                                    (df_summary['var_firing'] == dv) & \
                                    (df_summary['treatment'] == treatment) & \
                                    (df_summary['day'] == d) & \
                                    (df_summary['inj_current'] == inj_val)]
                    # if there is enough data for SE
                    if len(sum_data['SE']) > 0:
                        sem = sum_data['SE'].values[0]
                        ax[k].errorbar(i + 1, avg, yerr = sem, color = colors[treatment + day])

                counts.append(len(data))
                avgs.append(avg)
                inj.append(inj_val)

            counts_firing[treatment + day] = counts

            sum_data_tr = df_summary[(df_summary['data_type'] == dt_sum) & \
                                    (df_summary['var_firing'] == dv) & \
                                    (df_summary['treatment'] == treatment) & \
                                    (df_summary['day'] == d)]

            ax[k].plot(range(1, len(inj)+1), avgs, label = day, color = colors[treatment + day])

            if error_type == 'CI':
                start_summary = len(inj) - len(sum_data_tr) + 1
                ax[k].scatter(range(1, len(inj)+1), avgs, label = day, color = colors[treatment + day])
                ax[k].fill_between(range(start_summary, len(inj)+1), sum_data_tr['CI_L'], 
                                    sum_data_tr['CI_U'], color = colors[treatment + day],alpha = 0.3)
                if dv == 'num_aps':
                    ax[k].set_yticks(ticks = [0, 10, 20, 30])
                else:
                    ax[k].set_yticks(ticks = [0, 20, 40, 60, 80, 100])
            else:
                if dv == 'num_aps':
                    ax[k].set_yticks(ticks = [0, 10, 20, 30, 40, 50])
                else:
                    ax[k].set_yticks(ticks = [0, 40, 80, 120, 160, 200])

            ax[0].set_ylabel(dv_dict[dv][2])

    # bar graph after counting firing cells
    # x = np.arange(1, 14, 1)
    # ax[1,0].bar(x, counts_firing['CtrlD1'], color = colors['CtrlD1'])
    # ax[2,0].bar(x, counts_firing['CtrlD2'], color = colors['CtrlD2'])
    # ax[1,1].bar(x, counts_firing['high KD1'], color = colors['high KD1'])
    # ax[2,1].bar(x, counts_firing['high KD2'], color = colors['high KD2'])

    ticks_ = np.arange(1, 14, 2)
    labels_ = [inj[a] for a in np.arange(0,  13, 2)]
    ax[0].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
    ax[1].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
    ax[0].set_xlabel('Current injection (pA)')
    # ax[1,0].set_ylabel('pre')
    # ax[2,0].set_ylabel('Firing cells count \n post')
    ax[0].set_ylabel('Firing cells count')
    ax[0].set_title('CTR', color = colors['CtrlD2'])
    ax[1].set_title('HiK', color = colors['high KD2'])

    # # sharing only between the first row
    # for row in ax:
    #     for ax in row[0:]:
    #         ax.sharey(row[0])

    plt.subplots_adjust(hspace = 0.12, wspace = 0.15)
    
    plt.savefig(dest_dir + error_type + data_type + dv + '_vs_curernt_inj_slice.svg', \
                format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.close(fig)

def firing_props_fig2_slice(iff_df, data_type, dv, error_type, dest_dir,
                 background_dot_size = 3, w_cm = 13, h_cm = 6):
    """
    Plots the initial firing frequency and number of APs
    """

    iff_df = iff_df.sort_values(['treatment', 'day'])
    mask_exclude = (iff_df['treatment'] != 'high K') & (iff_df['treatment'] != 'Ctrl')
    iff_df = iff_df[~mask_exclude]
    num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(iff_df)

    colors = figure_colors()['slice']
    dt_sum = 'slice'

    dv_dict = {'IFF' : [iff_indx, 'Initial firing \n frequency (Hz)'],
                'num_aps' : [num_aps_indx, 'Firing frequency (Hz)']}

    fig, ax = plt.subplots(1, len(iff_df.treatment.unique()), figsize=(w_cm / 2.54, h_cm/2.54),\
                            sharex = True)

    counts_firing = {}
    for k, treatment in enumerate(iff_df.treatment.unique()):
        tr_df = iff_df[iff_df['treatment'] == treatment]

        for d, day in enumerate(tr_df['day'].unique()):
            day_df = tr_df[tr_df['day'] == day]
            print(treatment, day)
            avgs, inj, sems = [], [], []

            counts = []
            for i, col in enumerate(dv_dict[dv][0][5:5+13]):
                data = day_df.iloc[:,col]
                inj_val = int(day_df.columns[col][2:day_df.columns[col].rfind('pA')])

                # removing non-firing cells
                data.replace(0, np.nan, inplace = True)
                data.dropna(inplace = True)
                avg = np.mean(data)

                if error_type == 'SE':
                    x = np.linspace(0.75 + i , 1.25 + i, len(data))
                    # ax[k].scatter(x, data, alpha = 0.15, s = background_dot_size,
                    #             c = colors[treatment + day])
                    sem = np.std(data.values) / np.sqrt(len(data))
                    ax[k].errorbar(i + 1, avg, yerr = sem, color = colors[treatment + day])

                counts.append(len(data))
                avgs.append(avg)
                inj.append(inj_val)
                sems.append(sem)
                ## TO DO which lines does this have to be on
                ax[k].fill_between(range(1, len(inj)+1), np.array(avgs) + np.array(sems), \
                                    np.array(avgs) - np.array(sems), color = colors[treatment + day], alpha = 0.3)
            ax[k].plot(range(1, len(inj)+1), avgs, label = day, color = colors[treatment + day])
            counts_firing[treatment + day] = counts

            if dv == 'num_aps':
                ax[k].set_yticks(ticks = [0, 7.5, 15, 22.5, 30], labels = [0, '', 15, '', 30])
            else:
                ax[k].set_yticks(ticks = np.linspace(0, 100, 5), labels = [0, '', 50, '', 100])

    ticks_ = np.arange(1, 14, 2)
    labels_ = [inj[a] for a in np.arange(0,  13, 2)]
    ax[0].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
    ax[1].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
    ax[1].set_xlabel('Current injection (pA)')
    ax[0].set_ylabel(dv_dict[dv][1])
    ax[0].set_title('CTR', color = colors['CtrlD2'])
    ax[1].set_title('HiK', color = colors['high KD2'])

    plt.subplots_adjust(hspace = 0.12, wspace = 0.15)
    plt.savefig(dest_dir + error_type + data_type + dv + '_vs_curernt_inj_slice_V2.svg', \
                format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.close(fig)

def plot_example_firing(data_type, df, dest_dir, cell_IDs):
    '''
    cell_IDs - dictionary {'treatment': 'cell_ID}
    '''
    w_cm = 13
    h_cm = 9
    work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/'

    swps = [9, 13, 17]
    # colors = trace_colors_trace()[data_type]

    colors = figure_colors()[data_type]

    fig, ax = plt.subplots(3, 3, sharex = True, sharey = 'row',
                           figsize = (w_cm/2.54, h_cm/2.54), \
                        gridspec_kw = {'height_ratios': [6, 6, 1]})
    fig.subplots_adjust(hspace = 0.01,  wspace = 0.01)
    # Plot examples traces
    
    for j, item in enumerate(cell_IDs.items()):
        df_plot = df[df['cell_ID'] == item[1]]
        fn = work_dir + df_plot.OP.values[0] + '/' + df_plot.filename.values[0]
        channel = df_plot.cell_ch.values[0]
        inj = hcf.get_inj_current_steps(fn)
        trace = pyabf.ABF(fn)

        if len(trace.channelList) < 8:
            if '_Ipatch' in trace.adcNames:
                trace.adcNames[trace.adcNames.index('_Ipatch')] = 'Ch1'
            if 'IN0' in trace.adcNames:
                trace.adcNames[trace.adcNames.index('IN0')] = 'Ch1'
            channel_name = 'Ch' + str(channel)
            channel = trace.channelList[trace.adcNames.index(channel_name)]
        else:
            channel = channel-1
            channel_name = 'Ch' + str(channel+1)

        for i, sweep in enumerate(swps):
            trace.setSweep(sweepNumber = sweep, channel = channel)
            # ax[0, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
            #                 color = colors[item[0]+ '.' + str(i+1)], alpha = 0.7)
            ax[0+j, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
                            color = colors[item[0]], alpha = 0.7)
            if i == 0:
                ax[0+j, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
                            color = colors[item[0]], alpha = 0.7, label = item[0])
            # x = trace.sweepX
            # y = trace.sweepY

            ax[0+j,i].set_ylabel(trace.sweepLabelY)
            # ax[0,i].set_title(item[0])
            ax[0+j,i].axis('off')

            stim_inj = np.concatenate((np.zeros(2625), np.repeat(inj[sweep], 20_000), np.zeros(9375)))
            ax[2, i].plot(trace.sweepX [:-8000], stim_inj, c = 'black')
            # ax[1, i].set_xlabel(trace.sweepLabelX)
            # ax[1, i].set_ylabel(trace.sweepLabelC)
            ax[1+j, i].axis('off')
            ax[0, i].set_title(str(int(inj[sweep])) + ' pA')
    
    plt.subplots_adjust(hspace = 0.1, wspace = 0.01)
    # fig.legend(loc = 'upper center', ncol = 3, bbox_to_anchor=(0.5, 1.3))
    plt.show()
    plt.savefig(dest_dir + data_type + 'I_O_cureve.svg', \
                format = 'svg', bbox_inches = 'tight', dpi = 1000)

def intr_params_fig3(df, dest_dir, w_cm = 20, h_cm = 15.6):
    '''
    plots the repatch data
    no coloring based on age
    cold be added but plot become too much
    '''
    params = ['Rin', 'resting_potential', 'rheo_ramp_c', 'TH', 'sag']
    ax_elements = fig_order()['other']
    df = df.sort_values(['treatment', 'day'])
    titles_all = pl_intr.dict_for_plotting()
    titles_dict = {key: titles_all[key] for key in params if key in titles_all}
    colors = figure_colors()['cell_ID_new']
    dict_label = {'Ctrl':'CTR', 'high K':'HiK'}

    fig2, axs = plt.subplots(3,3, figsize = (w_cm/2.54, h_cm / 2.54))

    for param, item in titles_dict.items():
        ax_el = ax_elements[param]
        label_ = []
        for i, tr in enumerate(df.treatment.unique()):
            x_plot =[]
            for j, day in enumerate(df['day'].unique()):
                k = j + 2.5*i
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                df_plot.reset_index(drop=True, inplace=True)
                df_plot = df_plot[~df_plot[param].isna()]
                mean = np.nanmedian(df_plot[param])

                x = np.linspace(0.7+k, 1.3+k, len(df_plot))
                x_plot.append(x)

                axs[ax_el].scatter(x, df_plot[param], c=colors[tr + day], zorder = 5, alpha = 0.3)
                axs[ax_el].scatter(x, df_plot[param], edgecolor = colors[tr + day], alpha = 0.75,
                                   facecolors = 'none')
                axs[ax_el].scatter(1 + k, mean, color = 'k', marker = '_', s = 300, zorder = 10)

                if param == 'rheo_ramp_c':
                    axs[ax_el].text(0.75 + k, (1.05*mean), str(int(round(mean,0))), 
                                    size = 11, zorder = 11)
                elif param == 'sag':
                    axs[ax_el].text(0.75 + k, (1.05*mean), str(round(mean,2)), 
                                    size = 11, zorder = 11)
                else:
                    axs[ax_el].text(0.75 + k, (1.05*mean), str(round(mean,1)), 
                                    size = 11, zorder = 11)

                x_plot.append(x)
                if k in [1,3.5,5]:
                    for c, cell in enumerate(df_plot['cell_ID_new']):
                        x1 = [x_plot[0][c], x[c]]
                        y = df[param][df['cell_ID_new'] == cell]
                        axs[ax_el].plot(x1, y, '-', color = colors[tr + day],
                                alpha = 0.4, zorder = 1)
                if day =='D1':
                    x_ax = 'pre'
                elif day =='D2':
                    x_ax = 'post'
                label_.append(x_ax)

                # ax.text(k + 0.65, int((1 + 0.1)*np.max(df[param])), 'n = ' + str(len(df_plot)), size = 12, c = colors[tr + day])

        axs[ax_el].set_xticks(ticks = [1,2,3.5,4.5], labels = label_)
        axs[ax_el].set_ylabel(item[0] + ' (' + item[1] + ')')
        axs[ax_el].set_yticks(ticks = item[2], labels = item[3])
        # axs[ax_el].set_title('CTR                            HiK')

    date = str(datetime.date.today())
    plt.savefig(dest_dir + date + '_plot_all_repatch.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.close(fig2)

def firing_props_fig3(iff_df, error_type, dest_dir,
                 background_dot_size = 3, w_cm = 14.7, h_cm = 12):
    """
    Plots the initial firing frequency and number of APs
    """
    df_summary = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results'+ \
                           '/human/paper_figs_collected_checked/stats/'+ \
                            'sum_data_CIs_SEs_for_plotting/firing_plot_CIs_temporal.xlsx')
    iff_df = iff_df.sort_values(['treatment', 'day'])
    num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(iff_df)
    dt_sum = 'cell_ID_new'
    colors = figure_colors()[dt_sum]

    dv_dict = {'IFF' : [iff_indx, 'Initial firing frequency (AP#1 to AP#2)',
                'Initial firing \n frequency (Hz)'],
                'num_aps' : [num_aps_indx, 'Number of fired action potentials', 'AP count']}

    fig, ax = plt.subplots(4, 2, figsize=(w_cm / 2.54, h_cm/2.54),\
                            sharex = True, gridspec_kw = {'height_ratios': [5, 5, 1, 1]})

    for v, dv in enumerate(['num_aps', 'IFF']):
        counts_firing = {}
        for k, treatment in enumerate(iff_df.treatment.unique()):
            tr_df = iff_df[iff_df['treatment'] == treatment]

            for d, day in enumerate(tr_df['day'].unique()):
                day_df = tr_df[tr_df['day'] == day]
                print(treatment, day)

                avgs, inj, counts = [], [], []
                for i, col in enumerate(dv_dict[dv][0][5:5+13]):
                    data = day_df.iloc[:,col]
                    inj_val = int(day_df.columns[col][2:day_df.columns[col].rfind('pA')])

                    # removing non-firing cells
                    data.replace(0, np.nan, inplace = True)
                    data.dropna(inplace = True)
                    avg = np.mean(data)

                    if error_type == 'SE':
                        x = np.linspace(0.75 + i , 1.25 + i, len(data))
                        ax[v,k].scatter(x, data, alpha = 0.15, s = background_dot_size,
                                    c = colors[treatment + day])
                        sum_data = df_summary[(df_summary['data_type'] == dt_sum) & \
                                        (df_summary['var_firing'] == dv) & \
                                        (df_summary['treatment'] == treatment) & \
                                        (df_summary['day'] == d) & \
                                        (df_summary['inj_current'] == inj_val)]
                        # if there is enough data for SE
                        if len(sum_data['SE']) > 0:
                            sem = sum_data['SE'].values[0]
                            ax[v,k].errorbar(i + 1, avg, yerr = sem, color = colors[treatment + day])

                    counts.append(len(data))
                    avgs.append(avg)
                    inj.append(inj_val)

                counts_firing[treatment + day] = counts

                sum_data_tr = df_summary[(df_summary['data_type'] == dt_sum) & \
                                        (df_summary['var_firing'] == dv) & \
                                        (df_summary['treatment'] == treatment) & \
                                        (df_summary['day'] == d)]

                ax[v,k].plot(range(1, len(inj)+1), avgs, label = day, color = colors[treatment + day])

                if error_type == 'CI':
                    start_summary = len(inj) - len(sum_data_tr) + 1
                    ax[v,k].scatter(range(1, len(inj)+1), avgs, label = day, color = colors[treatment + day])
                    ax[v,k].fill_between(range(start_summary, len(inj)+1), sum_data_tr['CI_L'], 
                                        sum_data_tr['CI_U'], color = colors[treatment + day],alpha = 0.3)
                    if dv == 'num_aps':
                        ax[v,k].set_yticks(ticks = [0, 10, 20, 30])
                    else:
                        ax[v,k].set_yticks(ticks = [0, 20, 40, 60, 80, 100])
                else:
                    if dv == 'num_aps':
                        ax[v,k].set_yticks(ticks = [0, 10, 20, 30, 40, 50])
                    else:
                        ax[v,k].set_yticks(ticks = [0, 40, 80, 120, 160, 200])

                ax[v,0].set_ylabel(dv_dict[dv][2])

        if dv == 'IFF':
            # bar graph after counting firing cells
            x = np.arange(1, 14, 1)
            ax[2,0].bar(x, counts_firing['CtrlD1'], color = colors['CtrlD1'])
            ax[3,0].bar(x, counts_firing['CtrlD2'], color = colors['CtrlD2'])
            ax[2,1].bar(x, counts_firing['high KD1'], color = colors['high KD1'])
            ax[3,1].bar(x, counts_firing['high KD2'], color = colors['high KD2'])

        ticks_ = np.arange(1, 14, 2)
        labels_ = [inj[a] for a in np.arange(0,  13, 2)]

        # top plot
        # ax[0,0].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
        # ax[0,1].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
        
        ax[0,0].set_xticks(ticks_)
        ax[0,0].set_xticklabels(labels_, rotation=30)
        ax[0,1].set_xticks(ticks_)
        ax[0,1].set_xticklabels(labels_, rotation=30)

        # bottom
        ax[3,0].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
        ax[3,1].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
        ax[3,0].set_xlabel('Current injection (pA)')
        ax[3,0].set_ylabel('Firing cells count')
        ax[0,0].set_title('CTR', color = colors['CtrlD2'])
        ax[0,1].set_title('HiK', color = colors['high KD2'])

    plt.subplots_adjust(hspace = 0.12, wspace = 0.15)
    plt.savefig(dest_dir + error_type + dv + '_vs_curernt_inj_repatch.svg', \
                format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.close(fig)

def AIS_fig5(df, dest_dir, w_cm = 6.5, h_cm = 5.7):
    '''
    plots the AIS length data
    '''
    df = df.sort_values(['treatment'])
    dict_label = {'NI':'after \n slicing', 'Ctrl':'18h+ \n CTR', 'HighK':'18h+ \n HiK'}
    colors = figure_colors()['slice']

    fig, ax = plt.subplots(1, 1, figsize=(w_cm/2.54, h_cm / 2.54))
    labels_ = []
    for i, treatment in enumerate(['NI', 'Ctrl', 'HighK']):

        df_tr = df[df['treatment'] == treatment]
        x = np.linspace(0.7 + i, 1.3 + i, len(df_tr))
        ax.scatter(x, df_tr.Length, c = colors[treatment], zorder = 5, alpha = 0.3)
        ax.scatter(x, df_tr.Length, edgecolor = colors[treatment ],
                   facecolors = 'none', alpha = 0.75)

        ax.boxplot(df_tr.Length, positions = [i + 1], zorder = 10, widths = 0.45)
        labels_.append(dict_label[treatment])

    ax.set_xticks(ticks = [1,2,3], labels = labels_)
    ax.set_ylabel('AIS length (µm)')
    ax.set_yticks(ticks = [10, 20, 30, 40, 50, 60])

    date = str(datetime.date.today())
    plt.savefig(dest_dir  + date + '_plot_AIS.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.close(fig)

def intr_ext_single(dest_dir, w_cm = 17, h_cm = 9.5, er_type = 'SE'):
    '''
    plots the mean and CI/SE of intrinsic properties against age and time after OP
    '''
    stats_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
                  'paper_figs_collected_checked/stats/sum_data_CIs_SEs_for_plotting/'
    df = pd.read_excel(stats_dir + 'sum_data_age_hrs_temporal.xlsx')

    color_dict = figure_colors()
    shape_dict = {'D1': 'o', 'D2': 'v'}

    y_plot_dict = pl_intr.dict_for_plotting_reduced(er_type)
    x_plot_dict = {'hrs_after_OP':
                 {'ticks': np.linspace(10, 50, 5),
                  'labels': ['5-15', '16-25', '26-35', '36-45', '46-55'],
                  'axis_label': 'Time after tissue extraction (hrs)'},
                 'patient_age':
                 {'ticks': [15, 25, 35, 45, 60],
                  'labels': ['11-20', '21-30', '31-40', '41-50', '51+'],
                  'axis_label': 'Age (years)'}}

    for var in df['var'].unique():
        fig, ax = plt.subplots(2, 2,figsize = (w_cm/2.54, h_cm/2.54))
        for i, data_type in enumerate(df['data_type'].unique()):
            for j, x_param in enumerate(df.param.unique()):
                for d, day in enumerate(df.day.unique()):
                    for t, tr in enumerate(df.treatment_word.unique()):
                        k = d + 1.5*t
                        df_plot = df[(df['var'] == var) & \
                                        (df['data_type'] == data_type) & \
                                        (df['param'] == x_param) & \
                                        (df['day']== day) & \
                                        (df['treatment_word'] == tr)]
                        ax[j,i].scatter(df_plot['group'] + k, df_plot['mean'],
                                        marker = shape_dict[day],
                                        color = color_dict[data_type][tr+day])
                        
                        # !!! check for sure that the CIs are calcualteed correctly!
                        # for the python formula
                        if er_type == 'SE':
                            ebars = df_plot['SE']
                            ax[j,i].fill_between(df_plot['group'] + k, df_plot['mean'] - df_plot['SE'], 
                                        df_plot['mean'] + df_plot['SE'], 
                                        color = color_dict[data_type][tr+day], alpha = 0.3)
                        elif er_type == 'CI':
                            upper = abs(abs(df_plot['CI.U']) - abs(df_plot['mean'])).values
                            lower = abs(abs(df_plot['mean']) - abs(df_plot['CI.L'])).values
                            ebars = np.concat([upper, lower]).reshape(2,(len(upper)))
                            ax[j,i].fill_between(df_plot['group'] + k, df_plot['CI.L'], 
                                        df_plot['CI.U'], color = color_dict[data_type][tr+day],
                                        alpha = 0.3)

                        ax[j,i].errorbar(df_plot['group'] + k, df_plot['mean'],
                                            yerr = ebars, alpha = 0.8,
                                            color = color_dict[data_type][tr+day])
                        ax[j,i].set_xticks(ticks = x_plot_dict[x_param]['ticks'],
                                            labels = x_plot_dict[x_param]['labels'], rotation = 25)
                        ax[j,i].set_yticks(ticks = y_plot_dict[var][2])
                        ax[j,i].set_xlabel(x_plot_dict[x_param]['axis_label'])

        ax[0,0].set_ylabel(y_plot_dict[var][0] + ' (' + y_plot_dict[var][1] + ')')
        ax[1,0].set_ylabel(y_plot_dict[var][0] + ' (' + y_plot_dict[var][1] + ')')

        date = str(datetime.date.today())
        plt.savefig(dest_dir  + er_type + '_' + date + '_' +var + '.svg', \
            format = 'svg', bbox_inches = 'tight', dpi = 1000)
        plt.close(fig)

def OLD_intr_ext_fig4(dest_dir, w_cm = 22, h_cm = 21, er_type = 'SE'):
    '''
    plots the mean and CI/SE of intrinsic properties against age and time after OP
    '''
    stats_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
                  'paper_figs_collected_checked/stats/sum_data_CIs_SEs_for_plotting/'
    df = pd.read_excel(stats_dir + 'sum_data_age_hrs_temporal.xlsx')

    color_dict = figure_colors()
    shape_dict = {'D1': 's', 'D2': 'v'}

    y_plot_dict = pl_intr.dict_for_plotting_reduced(er_type)
    x_plot_dict = {'hrs_after_OP':
                 {'ticks': np.linspace(10, 50, 5),
                  'labels': ['5-15', '16-25', '26-35', '36-45', '46-55'],
                  'axis_label': 'Time after tissue extraction (hrs)'},
                 'patient_age':
                 {'ticks': [15, 25, 35, 45, 60],
                  'labels': ['11-20', '21-30', '31-40', '41-50', '51+'],
                  'axis_label': 'Age (years)'}}

    fig, ax = plt.subplots(3, 4,figsize = (w_cm/2.54, h_cm/2.54), sharey='row',
                           gridspec_kw={'hspace': 0.4, 'wspace': 0.21})
    for v, var in enumerate(['TH', 'rheo_ramp_c', 'sag']):
        for i, data_type in enumerate(['slice', 'cell_ID_new']):
            for j, x_param in enumerate(df.param.unique()):
                a =  i + 2*j
                for d, day in enumerate(df.day.unique()):
                    for t, tr in enumerate(df.treatment_word.unique()):
                        k = d + 1.5*t
                        df_plot = df[(df['var'] == var) & \
                                        (df['data_type'] == data_type) & \
                                        (df['param'] == x_param) & \
                                        (df['day']== day) & \
                                        (df['treatment_word'] == tr)]
                        ax[v,a].scatter(df_plot['group'] + k, df_plot['mean'],
                                        marker = shape_dict[day],
                                        color = color_dict[data_type][tr+day])
                        
                        # !!! check for sure that the CIs are calcualteed correctly!
                        # for the python formula
                        if er_type == 'SE':
                            df_plot['SE'] = df_plot['SE'].fillna(0)
                            ebars = df_plot['SE']
                            ax[v,a].fill_between(df_plot['group'] + k, df_plot['mean'] - df_plot['SE'],
                                        df_plot['mean'] + df_plot['SE'],
                                        color = color_dict[data_type][tr+day], alpha = 0.3)
                        elif er_type == 'CI':
                            upper = abs(abs(df_plot['CI.U']) - abs(df_plot['mean'])).values
                            lower = abs(abs(df_plot['mean']) - abs(df_plot['CI.L'])).values
                            ebars = np.concat([upper, lower]).reshape(2,(len(upper)))
                            ax[v,a].fill_between(df_plot['group'] + k, df_plot['CI.L'],
                                        df_plot['CI.U'], color = color_dict[data_type][tr+day],
                                        alpha = 0.3)

                        ax[v,a].errorbar(df_plot['group'] + k, df_plot['mean'],
                                            yerr = ebars, alpha = 0.8,
                                            color = color_dict[data_type][tr+day])

                        ax[v,a].set_xticks(ticks = x_plot_dict[x_param]['ticks'],
                                             labels = x_plot_dict[x_param]['labels'], rotation=35)
            
                ax[2,j+i].set_xlabel(x_plot_dict[x_param]['axis_label'])
        
        ax[v,0].set_yticks(ticks = y_plot_dict[var][2], labels = y_plot_dict[var][3])
        ax[v,0].set_ylabel(y_plot_dict[var][0] + '\n'  + '(' + y_plot_dict[var][1] + ')')

        date = str(datetime.date.today())
        plt.savefig(dest_dir  + er_type + '_' + date + '_complete.svg', \
            format = 'svg', bbox_inches = 'tight', dpi = 1000)
        # plt.close(fig)

def intr_params_inc_only(df, dest_dir, w_cm = 20, h_cm = 11):
    """
    Plots parameters for different treatments and days 
    """

    mask_exclude = (df['treatment'] != 'high K') & (df['treatment'] != 'Ctrl')
    df = df[~mask_exclude]
    params = ['Rin', 'resting_potential', 'rheo_ramp_c', 'TH', 'sag']
    # params = ['TH']
    ax_elements = fig_order()['inc_only']
    df = df.sort_values(['treatment', 'day'])
    titles_all = pl_intr.dict_for_plotting()
    titles_dict = {key: titles_all[key] for key in params if key in titles_all}

    dict_label = {'Ctrl':'CTR', 'high K':'HiK'}
    colors = figure_colors()['inc_only']

    fig2, axs = plt.subplots(2,3, figsize = (w_cm/2.54, h_cm / 2.54))
    for param, item in titles_dict.items():
        ax_el = ax_elements[param]
        label_, data_boxplot = [], []
        num_cels = {}
        for i, tr in enumerate(df.treatment.unique()):
            df_plot = df[(df['treatment'] == tr) & (df['day'] == 'D2')]
            k = i + 1
            df_plot = df_plot[~df_plot[param].isna()]  # remove NaN values
            x = np.linspace(0.45+k, 1.15+k, len(df_plot))

            axs[ax_el].scatter(x, df_plot[param], alpha = 0.3,
                        color = colors[tr], zorder = 5) # dots inside
            axs[ax_el].scatter(x, df_plot[param], alpha = 0.75,
                        edgecolor = colors[tr],
                        facecolors = 'none')

            axs[ax_el].boxplot(df_plot[param], positions = [k + 0.8], zorder = 10, widths = 0.5)

            label_.append( dict_label[tr])
            num_cels[tr] = len(df_plot)
            data_boxplot.append(df_plot[param])

        axs[ax_el].set_xticks(ticks = [1.8, 2.8], labels = ['post', 'post'])

        axs[ax_el].set_ylabel(item[0] + ' (' + item[1] + ')')
        axs[ax_el].set_yticks(ticks = item[2], labels = item[3])

    date = str(datetime.date.today())
    os.makedirs(dest_dir, exist_ok = True)
    plt.savefig(dest_dir  + date + 'inc_only' + '.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 1000)

def firing_props_inc_only(iff_df, error_type, dest_dir,
                 background_dot_size = 3, w_cm = 8, h_cm = 6.7):
    """
    Plots the initial firing frequency and number of APs for only incubated slices
    """
    df_summary = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/' + \
                              'paper_figs_collected_checked/stats/archive/sum_data_CIs_SEs_for_plotting' + \
                                '/firing_plot_CIs_inc_only.xlsx')
    iff_df = iff_df.sort_values(['treatment', 'day'])
    num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(iff_df)
    dt_sum = 'slice'
    colors = figure_colors()[dt_sum]

    dv_dict = {'IFF' : [iff_indx, 'Initial firing frequency (AP#1 to AP#2)',
                'Initial firing \n frequency (Hz)'],
                'num_aps' : [num_aps_indx, 'Number of fired action potentials', 'AP frequency (Hz)']}


    for v, dv in enumerate(['num_aps', 'IFF']):
        fig, ax = plt.subplots(1, 1, figsize=(w_cm / 2.54, h_cm/2.54))
        counts_firing = {}
        for k, treatment in enumerate(iff_df.treatment.unique()):
            

            for d, day in enumerate(iff_df['day'].unique()):
                if day == 'D1':
                    day_df = iff_df[iff_df['day'] == day]
                else:
                    day_df = iff_df[(iff_df['treatment'] == treatment) & \
                                (iff_df['day'] == day)]
                print(treatment, day)

                avgs, inj, counts = [], [], []
                for i, col in enumerate(dv_dict[dv][0][5:5+13]):
                    
                    data = day_df.iloc[:,col]
                    inj_val = int(day_df.columns[col][2:day_df.columns[col].rfind('pA')])

                    # removing non-firing cells
                    data.replace(0, np.nan, inplace = True)
                    data.dropna(inplace = True)
                    avg = np.mean(data)

                    if error_type == 'SE':
                        x = np.linspace(0.75 + i , 1.25 + i, len(data))
                        # ax.scatter(x, data, alpha = 0.15, s = background_dot_size,
                        #             c = colors[treatment + day])
                        sum_data = df_summary[(df_summary['var_firing'] == dv) & \
                                        (df_summary['treatment'] == treatment) & \
                                        (df_summary['day'] == d) & \
                                        (df_summary['inj_current'] == inj_val)]
                        # if there is enough data for SE
                        if len(sum_data['SE']) > 0:
                            sem = sum_data['SE'].values[0]
                            ax.errorbar(i + 1, avg, yerr = sem, color = colors[treatment + day])

                    inj.append(inj_val)
                    counts.append(len(data))
                    avgs.append(avg)

                counts_firing[treatment + day] = counts

                sum_data_tr = df_summary[(df_summary['var_firing'] == dv) & \
                                        (df_summary['treatment'] == treatment) & \
                                        (df_summary['day'] == d)]

                ax.plot(range(1, len(inj)+1), avgs, label = day, color = colors[treatment + day])

                if error_type == 'CI':
                    start_summary = len(inj) - len(sum_data_tr) + 1
                    # ax.scatter(range(1, len(inj)+1), avgs, label = day, color = colors[treatment + day])
                    ax.fill_between(range(start_summary, len(inj)+1), sum_data_tr['CI_L'], 
                                        sum_data_tr['CI_U'], color = colors[treatment + day],alpha = 0.3)
                #     if dv == 'num_aps':
                #         ax.set_yticks(ticks = [0, 10, 20, 30])
                #     else:
                #         ax.set_yticks(ticks = [0, 20, 40, 60, 80, 100])
                # else:
                #     if dv == 'num_aps':
                #         ax.set_yticks(ticks = [0, 10, 20, 30, 40, 50])
                #     else:
                #         ax.set_yticks(ticks = [0, 40, 80, 120, 160, 200])

                ax.set_ylabel(dv_dict[dv][2])

        if dv == 'num_aps':
            ax.set_yticks(ticks = [0, 5, 10, 15, 20, 25])
        else:
            ax.set_yticks(ticks = [0, 20, 40,60, 80, 100])

        ticks_ = np.arange(1, 14, 2)
        labels_ = [inj[a] for a in np.arange(0,  13, 2)]
        
        ax.set_xticks(ticks_)
        ax.set_xticklabels(labels_, rotation=30)

        ax.set_xlabel('Current injection (pA)')

        plt.subplots_adjust(hspace = 0.12, wspace = 0.15)
        plt.savefig(dest_dir + error_type + dv + '_vs_curernt_inj_inc_only.svg', \
                    format = 'svg', bbox_inches = 'tight', dpi = 1000)
        plt.close(fig)

def firing_props_inc_only_V2(iff_df, error_type, dest_dir,
                 background_dot_size = 3, w_cm = 6, h_cm = 6):
    """
    Plots the initial firing frequency and number of APs for only incubated slices
    """

    iff_df = iff_df.sort_values(['treatment', 'day'])
    num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(iff_df)
    colors = figure_colors()['inc_only']

    dv_dict = {'IFF' : [iff_indx, 'Initial firing frequency (AP#1 to AP#2)',
                'Initial firing \n frequency (Hz)'],
                'num_aps' : [num_aps_indx, 'Number of fired action potentials', 'AP frequency (Hz)']}

    for v, dv in enumerate(['num_aps', 'IFF']):
        fig, ax = plt.subplots(1, 1, figsize=(w_cm / 2.54, h_cm/2.54))
        counts_firing = {}
        for k, treatment in enumerate(iff_df.treatment.unique()):
                tr_df = iff_df[iff_df['treatment'] == treatment]
                avgs, inj, counts, sems = [], [], [], []
                for i, col in enumerate(dv_dict[dv][0][5:5+13]):
                    
                    data = tr_df.iloc[:,col]
                    inj_val = int(tr_df.columns[col][2:tr_df.columns[col].rfind('pA')])

                    # removing non-firing cells
                    data.replace(0, np.nan, inplace = True)
                    data.dropna(inplace = True)
                    avg = np.mean(data)

                    if error_type == 'SE':
                        x = np.linspace(0.75 + i , 1.25 + i, len(data))
                        sem = np.std(data) / np.sqrt(len(data))

                        ax.errorbar(i + 1, avg, yerr = sem, color = colors[treatment])

                    sems.append(sem)
                    inj.append(inj_val)
                    counts.append(len(data))
                    avgs.append(avg)

                counts_firing[treatment] = counts
                ax.plot(range(1, len(inj)+1), avgs, label = treatment, color = colors[treatment])
                ax.set_ylabel(dv_dict[dv][2])
            
                ax.fill_between(range(1, len(inj)+1), np.array(avgs) + np.array(sems), \
                        np.array(avgs) - np.array(sems), color = colors[treatment], alpha = 0.3)
        if dv == 'num_aps':
            ax.set_yticks(ticks = [0, 5, 10, 15, 20, 25])
        else:
            ax.set_yticks(ticks = [0, 20, 40,60, 80, 100])

        ticks_ = np.arange(1, 14, 2)
        labels_ = [inj[a] for a in np.arange(0,  13, 2)]
        
        ax.set_xticks(ticks_)
        ax.set_xticklabels(labels_, rotation=30)
        ax.set_xlabel('Current injection (pA)')

        plt.subplots_adjust(hspace = 0.12, wspace = 0.15)
        plt.savefig(dest_dir + error_type + dv + '_vs_curernt_inj_inc_only_V2.svg', \
                    format = 'svg', bbox_inches = 'tight', dpi = 1000)
        plt.close(fig)

def intr_slice_CTR(dest_dir, w_cm = 21, h_cm = 12, er_type = 'SE',
                   params = ['Rin', 'resting_potential', 'TH', 'rheo_ramp_c', 'sag']):
    '''
    plots the mean and CI/SE of intrinsic properties against age and time after OP
    '''
    stats_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
                  'paper_figs_collected_checked/stats/for_plot/'
    
    df = pd.read_excel(stats_dir + 'sum_data_slice_CTR_hrs.xlsx')

    color_dict = figure_colors()['slice_all_CTR']
    ax_elements = fig_order()['slice_all_CTR']
    y_plot_dict = pl_intr.dict_for_plotting_reduced(er_type)
    x_plot_dict = {'ticks': np.linspace(10, 50, 5),
                  'labels': ['5-15', '16-25', '26-35', '36-45', '46-55'],
                  'axis_label': 'Time after tissue extraction (hrs)'}

    fig, ax = plt.subplots(2, 3,figsize = (w_cm/2.54, h_cm/2.54),
                           gridspec_kw={'hspace': 0.7, 'wspace': 0.65})
    for v, var in enumerate(params):
        pos = ax_elements[var]

        df_plot = df[df['var'] == var]
        ax[pos].scatter(df_plot['group'], df_plot['mean'],
                        marker = 's',
                        color = color_dict['CTR'])
        
        # for the python formula
        if er_type == 'SE':
            df_plot['SE'] = df_plot['SE'].fillna(0)
            ebars = df_plot['SE']
            ax[pos].fill_between(df_plot['group'], df_plot['mean'] - df_plot['SE'],
                        df_plot['mean'] + df_plot['SE'],
                        color = color_dict['CTR'], alpha = 0.3)
        elif er_type == 'CI':
            upper = abs(abs(df_plot['CI.U']) - abs(df_plot['mean'])).values
            lower = abs(abs(df_plot['mean']) - abs(df_plot['CI.L'])).values
            ebars = np.concat([upper, lower]).reshape(2,(len(upper)))
            ax[pos].fill_between(df_plot['group'], df_plot['CI.L'],
                        df_plot['CI.U'], color = color_dict['CTR'],
                        alpha = 0.3)

        ax[pos].errorbar(df_plot['group'], df_plot['mean'],
                            yerr = ebars, alpha = 0.8,
                            color = color_dict['CTR'])
        ax[pos].set_xticks(ticks = x_plot_dict['ticks'],
                                labels = x_plot_dict['labels'], rotation=35)
        ax[pos].set_yticks(ticks = y_plot_dict[var][2], labels = y_plot_dict[var][3])
        ax[pos].set_ylabel(y_plot_dict[var][0] + '\n'  + '(' + y_plot_dict[var][1] + ')')
        ax[0,1].set_xlabel(x_plot_dict['axis_label'])
        ax[1,1].set_xlabel(x_plot_dict['axis_label'])

        date = str(datetime.date.today())
        plt.savefig(dest_dir  + er_type + '_' + date + '_slice_all_CTR.svg', \
            format = 'svg', bbox_inches = 'tight', dpi = 1000)
        # plt.close(fig)

def slice_all_CTR_hist(df, dest_dir, w_cm = 6.5, h_cm = 5.5):
    '''
    plots a histogram of the time after OP for all CTR cells in slice condition
    '''
    df = df[ ((df['treatment'] == 'Ctrl') & (df['day'] == 'D2')) |
                                    (df['day'] == 'D1')]
    df.reset_index(drop = True, inplace = True)
    print('Number of cells in slice CTR:', len(df))

    fig, ax = plt.subplots(1,1, figsize = (w_cm/2.54, h_cm/2.54))
    ax.hist(df.hrs_after_OP, bins = np.arange(5, 60, 5), alpha = 0.8,
            color = figure_colors()['slice_all_CTR']['CTR'], edgecolor = figure_colors()['slice_all_CTR']['CTR'])

    ax.set_xticks(ticks = np.arange(0,60,10), labels = ['', 10, '', 30, '', 50])
    ax.set_yticks(ticks = np.arange(0,100,10), labels = ['', 10, '', 30, '', 50, '', 70, '',90])
    ax.set_ylabel('Number of cells')

    date = str(datetime.date.today())
    plt.savefig(dest_dir + 'hist' + date + '_slice_all_CTR.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 1000)

def plot_firing_example_trace_messy_grid():
        
    # plot firing in a grid
    cell_IDs = {'pre': '24201S1c8',
                'post CTR': '24117S2_D2c1',
                'post HiK': '23420S2_D2c7'}
    data_type = 'slice'
    df = df_slice
    dest_dir = destination_dir

    w_cm = 15
    h_cm = 4
    work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/'

    swps = [9, 13, 17]
    # colors = trace_colors_trace()[data_type]
    colors = figure_colors()[data_type]

    fig, ax = plt.subplots(4, 3, sharex = True, sharey = 'row',
                            figsize = (w_cm/2.54, h_cm/2.54), \
                        gridspec_kw = {'height_ratios': [6, 6, 6, 2]})
    fig.subplots_adjust(hspace = 0.01,  wspace = 0.01)
    # Plot examples traces

    for j, item in enumerate(cell_IDs.items()):
        df_plot = df[df['cell_ID'] == item[1]]
        fn = work_dir + df_plot.OP.values[0] + '/' + df_plot.filename.values[0]
        channel = df_plot.cell_ch.values[0]
        inj = hcf.get_inj_current_steps(fn)
        trace = pyabf.ABF(fn)

        if len(trace.channelList) < 8:
            if '_Ipatch' in trace.adcNames:
                trace.adcNames[trace.adcNames.index('_Ipatch')] = 'Ch1'
            if 'IN0' in trace.adcNames:
                trace.adcNames[trace.adcNames.index('IN0')] = 'Ch1'
            channel_name = 'Ch' + str(channel)
            channel = trace.channelList[trace.adcNames.index(channel_name)]
        else:
            channel = channel-1
            channel_name = 'Ch' + str(channel+1)

        for i, sweep in enumerate(swps):
            trace.setSweep(sweepNumber = sweep, channel = channel)
            # ax[0, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
            #                 color = colors[item[0]+ '.' + str(i+1)], alpha = 0.7)
            ax[j, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
                            color = colors[item[0]], alpha = 0.7)
            if i == 0:
                ax[j, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
                            color = colors[item[0]], alpha = 0.7, label = item[0])
            # x = trace.sweepX
            # y = trace.sweepY

            ax[j,i].set_ylabel(trace.sweepLabelY)
            # ax[0,i].set_title(item[0])
            ax[j,i].axis('off')

            stim_inj = np.concatenate((np.zeros(2625), np.repeat(inj[sweep], 20_000), np.zeros(9375)))
            ax[3, i].plot(trace.sweepX [:-8000], stim_inj, c = 'black')
            # ax[1, i].set_xlabel(trace.sweepLabelX)
            # ax[1, i].set_ylabel(trace.sweepLabelC)
            ax[3, i].axis('off')
            ax[0,i].set_title(str(int(inj[sweep])) + ' pA')
        # fig.legend(loc = 'upper center', ncol = 3, bbox_to_anchor=(0.5, 1.3))
    # plt.savefig(dest_dir + data_type + 'I_O_curve_grid.svg', \
    #            format = 'svg', bbox_inches = 'tight', dpi = 1000)

def plot_corrs(x, df_list, dir_read_corrs, save_dir_figs, params_list = []):
    '''
    x - x_axis parameter
    df_list - dictionary. 'df_name' = [df, 'color']
    save_dir_corr - where the correlations dataframe is saved
    '''
    w_cm =  6 * len(params_list)
    h_cm = 6
    
    titles_dict = pl_intr.dict_for_plotting()
    for df_name, df_info in df_list.items():

        df = df_info[0]
        col = df_info[1]

        df_cor = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'cor_estimate')
        df_p_val = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'p_vals')
        df_l_CI = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'lower_CI')
        df_u_CI = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'upper_CI')

        df_cor.index = [df_cor.columns]
        df_p_val.index = [df_p_val.columns]
        df_l_CI.index = [df_l_CI.columns]
        df_u_CI.index = [df_u_CI.columns]

        # cor.loc['hrs_after_OP', 'hrs_after_OP']
        params_list = df_u_CI.columns.to_list()
        params_list.remove(x)

        fig, ax = plt.subplots(1, len(params_list), figsize = (w_cm/2.54, h_cm/2.54))
        for i, param in enumerate(params_list):

            cor = df_cor.loc[x, param].iloc[0]
            p_val = df_p_val.loc[x, param].iloc[0]
            sns.regplot(x=df[x], y=df[param], ax = ax[i], color = col, scatter = False)
            ax[i].scatter(df[x], df[param], c = 'darkgrey', alpha = 0.65, zorder = 0, s = 3.5)
        
            # Add correlation info to plot
            ax[i].text(0.95, 0.95, f'r = {cor:.3f}', transform=ax[i].transAxes,
                    va = 'top', ha = 'right')
            ax[i].text(0.95, 0.85, f'p = {p_val:.3f}', transform=ax[i].transAxes,
                    va = 'top', ha = 'right')
            # ax[i].text(0.05, 0.75, f'CI: [{l_CI:.3f}, {u_CI:.3f}]',
            #            transform=ax[i].transAxes, verticalalignment='top', fontsize=8)
            
            # Set labels and title
            ax[i].set_xlabel(f'{titles_dict[x][0]} ({titles_dict[x][1]})')
            ax[i].set_xticks(ticks = titles_dict[x][2], labels = [int(x) for x in titles_dict[x][2]])
            # ax[i].set_ylabel(f'{titles_dict[col_name][0]} \n ({titles_dict[col_name][1]})')
            # ax[i].set_yticks(ticks = titles_dict[param][2], labels = titles_dict[param][3])
            # ax[i].set_title(f'{param}')
        
        plt.savefig(f'{save_dir_figs}{df_name}corr_plots.svg',
                    format='svg', bbox_inches='tight', dpi=300)

def plot_example_firing_corr(data_type, df, dest_dir, cell_IDs):
    '''
    hi
    '''
    w_cm = 13
    h_cm = 9
    work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/'

    swps = [9, 13, 17]
    # colors = trace_colors_trace()[data_type]

    colors = plot_figs.figure_colors()[data_type]
    num_rows = len(cell_IDs)
    fig, ax = plt.subplots(num_rows+1, 3, sharex = True, sharey = 'row',
                           figsize = (w_cm/2.54, h_cm/2.54), \
                        gridspec_kw = {'height_ratios': [6] * num_rows + [1]})
    fig.subplots_adjust(hspace = 0.01,  wspace = 0.01)
    # Plot examples traces
    
    for j, item in enumerate(cell_IDs.items()):
        df_plot = df[df['cell_ID'] == item[1]]
        fn = work_dir + df_plot.OP.values[0] + '/' + df_plot.filename.values[0]
        channel = df_plot.cell_ch.values[0]
        inj = hcf.get_inj_current_steps(fn)
        trace = pyabf.ABF(fn)

        if len(trace.channelList) < 8:
            if '_Ipatch' in trace.adcNames:
                trace.adcNames[trace.adcNames.index('_Ipatch')] = 'Ch1'
            if 'IN0' in trace.adcNames:
                trace.adcNames[trace.adcNames.index('IN0')] = 'Ch1'
            channel_name = 'Ch' + str(channel)
            channel = trace.channelList[trace.adcNames.index(channel_name)]
        else:
            channel = channel-1
            channel_name = 'Ch' + str(channel+1)

        for i, sweep in enumerate(swps):
            trace.setSweep(sweepNumber = sweep, channel = channel)
            # ax[0, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
            #                 color = colors[item[0]+ '.' + str(i+1)], alpha = 0.7)
            ax[0+j, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
                            color = colors[item[0]], alpha = 0.7)
            if i == 0:
                ax[0+j, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
                            color = colors[item[0]], alpha = 0.7, label = item[0])
            # x = trace.sweepX
            # y = trace.sweepY

            ax[0+j,i].set_ylabel(trace.sweepLabelY)
            # ax[0,i].set_title(item[0])
            ax[0+j,i].axis('off')

            stim_inj = np.concatenate((np.zeros(2625), np.repeat(inj[sweep], 20_000), np.zeros(9375)))
            ax[num_rows, i].plot(trace.sweepX [:-8000], stim_inj, c = 'black')
            # ax[1, i].set_xlabel(trace.sweepLabelX)
            # ax[1, i].set_ylabel(trace.sweepLabelC)
            ax[1+j, i].axis('off')
            ax[0, i].set_title(str(int(inj[sweep])) + ' pA')
    
    plt.subplots_adjust(hspace = 0.1, wspace = 0.01)
    # fig.legend(loc = 'upper center', ncol = 3, bbox_to_anchor=(0.5, 1.3))
    plt.savefig(dest_dir + data_type + 'I_O_curve.svg', \
                format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.show()

def plot_non_firing_cells_bar(df, data_type, save_dir, inj_vals):
    '''
    plots the counts of non-firing cells at given injection value
    '''
    w_cm = 5
    h_cm = 5
    colors = figure_colors()[data_type]
    # titles_dict = pl_intr.dict_for_plotting()
    offset_dict = {'Ctrl': -0.2, 'high K': 0.2}


    fig, ax = plt.subplots(1, 1, figsize = (w_cm/2.54, h_cm/2.54))
    for tr in sorted(df.treatment.unique()):
        for i, inj in enumerate(inj_vals):
            x = 1 + i + offset_dict[tr]
            if data_type != 'inc_only':
                y_D1 = df.count_non_firing_cells[(df.treatment == tr) &
                                                    (df.day == 'D1') &
                                                    (df.injection == inj)].values[0]
                y_D2 = df.count_non_firing_cells[(df.treatment == tr) &
                                                    (df.day == 'D2') &
                                                    (df.injection == inj)].values[0]
            
                ax.bar(x, y_D2/y_D1, width = 0.15, color = colors[tr + 'D2'])
                ax.set_ylabel('proportion non-firing cells \n (D2/D1)')
                ax.set_yticks(ticks = np.linspace(0, 2, 5), labels = [0, '', 1, '', 2])
            else:
                y = df.count_non_firing_cells[(df.treatment == tr) &
                                                    (df.injection == inj)].values[0]
                ax.bar(x, y, width = 0.15, color = colors[tr])
                ax.set_ylabel('Non-firing cells (count)')
                ax.set_yticks(ticks = np.linspace(0, 12, 5), labels = [0, '', 6, '', 12])

    ax.set_xlabel('Current injection (pA)')
    ax.set_xticks(ticks  = np.arange(1,len(inj_vals)+1), labels = inj_vals)

    plt.savefig(f'{save_dir}{data_type}_non_firing_cells_bar.svg',
            format='svg', bbox_inches='tight', dpi=300)

def plot_non_firing_across_inj_proportion(df, data_type, save_dir):
    '''
    '''
    w_cm = 6
    h_cm = 6
    colors = figure_colors()[data_type]

    fig, ax = plt.subplots(1, 1, figsize = (w_cm/2.54, h_cm/2.54))
    x = df.injection.unique()[5:18]
    for tr in sorted(df.treatment.unique()):
        if data_type != 'inc_only':
            y_D1 = df.count_non_firing_cells[(df.treatment == tr) &
                                                (df.day == 'D1')].values[5:18]
            y_D2 = df.count_non_firing_cells[(df.treatment == tr) &
                                                (df.day == 'D2')].values[5:18]
            # replace 0 with small number
            y_D1 = y_D1 + 1
            y_D2 = y_D2 + 1
            ax.plot(x, y_D2/y_D1, color = colors[tr + 'D2'], alpha = 0.75)
            ax.scatter(x, y_D2/y_D1, color = colors[tr + 'D2'], alpha = 0.85)
            ax.set_ylabel('proportion non-firing cells \n (D2/D1)')
            #ax.set_yticks(ticks = np.linspace(0, 2, 5), labels = [0, '', 1, '', 2])
        else:
            y = df.count_non_firing_cells[(df.treatment == tr)].values[5:18]
            ax.plot(x, y, color = colors[tr], alpha = 0.75)
            ax.scatter(x, y, color = colors[tr], alpha = 0.85)
            ax.set_ylabel('Non-firing cells \n (count)')
            #ax.set_yticks(ticks = np.linspace(0, 12, 5), labels = [0, '', 6, '', 12])

    ax.set_xlabel('Current injection (pA)')
    labels_ = []
    for l, val in enumerate(x):
        if l % 2 == 0:
            labels_.append(val)
        else:
            labels_.append('')
    ax.set_xticks(ticks  = x, labels = labels_, rotation = 45)

    plt.savefig(f'{save_dir}{data_type}_non_firing_cells.svg',
            format='svg', bbox_inches='tight', dpi=300)

def plot_non_firing_across_inj_count(df, data_type, save_dir):
    '''
    '''
    w_cm = 12
    h_cm = 6
    colors = figure_colors()[data_type]

    fig, ax = plt.subplots(1, 2, figsize = (w_cm/2.54, h_cm/2.54))
    x = df.injection.unique()[5:18]
    for i, tr in enumerate(sorted(df.treatment.unique())):
        for day in ['D1', 'D2']:
            if data_type != 'inc_only':
                y = df.count_non_firing_cells[(df.treatment == tr) &
                                                    (df.day == day)].values[5:18]

                ax[i].plot(x, y, color = colors[tr + day], alpha = 0.75)
                ax[i].scatter(x, y, color = colors[tr + day], alpha = 0.75)
                ax[i].set_ylabel('Non-firing cells \n (count)')
                ax[i].set_yticks(ticks = np.linspace(0, 36, 5), labels = [0, '', 18, '', 36])
            else:
                y = df.count_non_firing_cells[(df.treatment == tr)].values[5:18]
                ax[0].plot(x, y, color = colors[tr], alpha = 0.75)
                ax[0].scatter(x, y, color = colors[tr], alpha = 0.85)
                ax[0].set_ylabel('Non-firing cells \n (count)')
                ax[0].set_yticks(ticks = np.linspace(0, 36, 5), labels = [0, '', 18, '', 36])

    ax[0].set_xlabel('Current injection (pA)')
    labels_ = []
    for l, val in enumerate(x):
        if l % 2 == 0:
            labels_.append(val)
        else:
            labels_.append('')
    ax[0].set_xticks(ticks  = x, labels = labels_, rotation = 45)

    plt.savefig(f'{save_dir}{data_type}_non_firing_cells.svg',
            format='svg', bbox_inches='tight', dpi=300)

def plot_corrs_firing(x, df_list, dir_read_corrs, save_dir_figs, params_list = []):
    '''
    x - x_axis parameter
    df_list - dictionary. 'df_name' = [df, 'color']
    save_dir_corr - where the correlations dataframe is saved
    '''
    w_cm = 16
    h_cm = 10
    
    titles_dict = pl_intr.dict_for_plotting()
    for df_name, df_info in df_list.items():

        df = df_info[0]
        col = df_info[1]

        df_cor = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'cor_estimate')
        df_p_val = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'p_vals')
        df_l_CI = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'lower_CI')
        df_u_CI = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'upper_CI')

        # cor.loc['hrs_after_OP', 'hrs_after_OP']
        params_list = df_u_CI.columns.to_list()
        params_list.remove(x)

        fig, ax = plt.subplots(2, int(len(params_list)/2), figsize = (w_cm/2.54, h_cm/2.54))
        for p, param in enumerate(params_list):
            cor = df_cor[param].iloc[0]
            p_val = df_p_val[param].iloc[0]
            if 'num_aps' in param:
                col_name = f"('{param[3:6]}pA', 'num_aps')"
                i = 0
                j = int(p /2)
            else:
                col_name = f"('{param[3:6]}pA', 'IFF')"
                i = 1
                j = int(p/2 - 0.5)
            
            if 'ALL' not in df_name:
                # take only firing cells
                df_plot = df[df[col_name].notna() & (df[col_name] != 0)]
            
            sns.regplot(x=df_plot[x], y=df_plot[col_name], ax = ax[i, j], color = col, scatter = False)
            ax[i, j].scatter(df_plot[x], df_plot[col_name], c = 'darkgrey', alpha = 0.65, zorder = 0, s = 3.5)
  
            # Add correlation info to plot
            ax[i, j].text(0.95, 0.95, f'r = {cor:.3f}', transform=ax[i, j].transAxes,
                    va = 'top', ha = 'right')
            ax[i, j].text(0.95, 0.85, f'p = {p_val:.3f}', transform=ax[i, j].transAxes,
                    va = 'top', ha = 'right')
            # ax[i].text(0.05, 0.75, f'CI: [{l_CI:.3f}, {u_CI:.3f}]',
            #            transform=ax[i].transAxes, verticalalignment='top', fontsize=8)
            
            # Set labels and title
            ax[i, j].set_xlabel(f'{titles_dict[x][0]} ({titles_dict[x][1]})')
            ax[i, j].set_xticks(ticks = titles_dict[x][2], labels = [int(x) for x in titles_dict[x][2]])
            # ax[i].set_ylabel(f'{titles_dict[col_name][0]} \n ({titles_dict[col_name][1]})')
            # ax[i].set_yticks(ticks = titles_dict[param][2], labels = titles_dict[param][3])
            # ax[i].set_title(f'{param}')
        
        plt.savefig(f'{save_dir_figs}{df_name}corr_plots.svg',
                    format='svg', bbox_inches='tight', dpi=300)


# FUNCTIONS PLOT CONNECT

def plot_connectivity_percentage(df_connections, data_type, save_dir = False):

    w_cm = 6
    h_cm = 6
    fig, ax = plt.subplots(1, 1,figsize=(w_cm/2.54,h_cm/2.54))
    colors = figure_colors()[data_type]

    for j, day in enumerate(sorted(df_connections.day.unique())):
        for i, tr in enumerate(sorted(df_connections.treatment.unique())):

            k_ = j + 2.5*i
            df_plot = df_connections.loc[(df_connections.treatment == tr)\
                                         & (df_connections.day == day)]
            
            x = np.linspace(0.65+k_, 1.35+k_, len(df_plot))
            ax.scatter(x, df_plot.con_percentage_fixed, c = colors[tr+day])

            avg = np.mean(df_plot.con_percentage_fixed)
            ax.scatter(k_+1, avg, marker = '_', s = 300, c = 'k')
            
            ax.text(k_+0.85, 1.2*avg, str(round(avg,0)), size = 10)
    
    ax.set_xticks(ticks = [1,2,3.5,4.5], labels = ['pre', 'post', 'pre', 'post'])
    ax.set_yticks(ticks = np.linspace(0, 16, 5), labels = [0, '', 8, '', 16])
    
    if save_dir:
        plt.savefig(f'{save_dir}{data_type}connectivity_percentage.svg',
            format='svg', bbox_inches='tight', dpi=300)
    plt.show()

def plot_connect_trace_repatched_cell(save_dir = []):
    """
    Plot connection screens for D1 and D2 side by side in a single figure.
    """
    # Configuration
    line_width = 3
    human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
    OP = 'OP240926'
    patcher = 'Verji'
    active_channels = [3, 6]
    colors = ['#f781bf', '#984ea3']
    file_indices = [46, 90]  # D1, D2
    w_cm = 16
    h_cm = 12
    
    stim_windows = {
        'Ch1': [5000, 14600], 'Ch2': [25000, 34600], 'Ch3': [45000, 54600],
        'Ch4': [65000, 74600], 'Ch5': [85000, 94600], 'Ch6': [105000, 114600],
        'Ch7': [125000, 134600], 'Ch8': [142300, 154800]
    }
    
    # Get metadata
    work_dir, filenames, _, _, _, _ = sort.get_OP_metadata(human_dir, OP, patcher)
    filenames_to_plot = [work_dir + filenames[idx] for idx in file_indices]    
    # Create single figure with both plots side by side
    fig, ax = plt.subplots(len(active_channels), len(active_channels),
                            figsize=(w_cm/2.54, h_cm/2.54))

    for day_indx, filename in enumerate(filenames_to_plot):
        # Load data
        con_screen_data = hcf.load_traces(filename)
        
        row_indx_pre = 0
        row_indx_post = 1

        pre_ch = active_channels[row_indx_pre]
        pre_data = con_screen_data[f'Ch{pre_ch}'][0]
        pre_avg = np.mean(pre_data, axis=1)
        start, end = stim_windows[f'Ch{pre_ch}']

        post_ch = active_channels[row_indx_post]
        post_data = con_screen_data[f'Ch{post_ch}'][0]
        post_avg = np.mean(post_data, axis=1)
                
        ax[row_indx_pre, day_indx].plot(pre_avg[start:end], colors[row_indx_pre], lw=line_width)
        ax[row_indx_post, day_indx].plot(post_avg[start:end], colors[row_indx_post], lw=line_width)
        ax[0, day_indx].set_title(f'Day {day_indx+1}')
    # Set consistent y-limits for each row
    padding = [0.1, 0.3]
    for j, signal_ in enumerate([pre_avg[start:end], post_avg[start:end]]):

        y_min, y_max = signal_.min(), signal_.max()
        y_range = (y_max - y_min) * padding[j]

        for col_indx in [0,1]:
            ax[j, col_indx].set_ylim([y_min - y_range, y_max + y_range])

    plt.savefig(f'{save_dir}/con_screen_{OP}Ch{pre_ch}_Ch{post_ch}.svg',
            format='svg', bbox_inches='tight', dpi=300)

def plot_corr_single(x, y, df, data_type = 'all', tr = 'Ctrl', save_dir = False):
    '''
    
    '''
    colors = figure_colors()[data_type][tr]
    w_cm = 6
    h_cm = 6
    
    titles_dict = pl_intr.dict_for_plotting_conn()

    cor, p_val = stats.pearsonr(df[x], df[y])

    fig, ax = plt.subplots(1, 1, figsize = (w_cm/2.54, h_cm/2.54))
    sns.regplot(x = df[x], y = df[y], ax = ax, color = colors, scatter = False)
    ax.scatter(df[x], df[y], c = 'darkgrey', alpha = 0.65, zorder = 0, s = 3.5)

    # Add correlation info to plot
    ax.text(0.95, 0.95, f'r = {cor:.3f}', transform=ax.transAxes,
            va = 'top', ha = 'right')
    ax.text(0.95, 0.85, f'p = {p_val:.3f}', transform=ax.transAxes,
            va = 'top', ha = 'right')

    # Set labels and title
    ax.set_xlabel(f'{titles_dict[x][0]} ({titles_dict[x][1]})')
    ax.set_xticks(ticks = titles_dict[x][2], labels = [int(x) for x in titles_dict[x][2]])
    ax.set_ylabel(f'{titles_dict[y][0]} \n ({titles_dict[y][1]})')
    ax.set_yticks(ticks = titles_dict[y][2], labels = titles_dict[y][3])

    if save_dir:
        plt.savefig(f'{save_dir}{data_type}{x}{y}corr_plots.svg',
                    format='svg', bbox_inches='tight', dpi=300)

def plot_norm_con_params(df, DV, data_type, save_dir = False):
    '''
    DV - 'amp' or 'lat'
    data_type = 'slice' for all; cell_ID_new for repatch (colors)
    '''
    titles_dict = pl_intr.dict_for_plotting_conn()
    w_cm, h_cm = 12, 8
    colors = figure_colors()[data_type]

    if DV == 'amp':
        cols = ['Amp 1', 'Amp 2', 'Amp 3', 'Amp 4']
        param_1 = 'Amp 1'
    else:
        cols = ['Lat1', 'Lat2', 'Lat3', 'Lat4']
        param_1 = 'Lat1'

    fig, ax = plt.subplots(1, 1, figsize=(w_cm/2.54, h_cm/2.54))
    x_positions = np.arange(1, len(cols) + 1)  # 1, 2, 3, 4 for Amp1-4

    for day in sorted(df.day.unique()):
        for tr in sorted(df.treatment.unique()):
            if day == 'D1':
                df_plot = df[df.day == day]
            else:
                df_plot = df[(df.day == day) & (df.treatment == tr)]
            
            if len(df_plot) == 0:
                continue
            df_plot = df_plot.dropna(subset=cols)
            df_plot = df_plot[(df_plot[cols] > 0).all(axis=1)]
                
            # Calculate normalized amplitudes (divide by Amp 1)
            normalized_means = []
            normalized_errors = []
            
            for amp_col in cols:
                
                # Normalize to Amp 1
                normalized_values = df_plot[amp_col] / df_plot[param_1]
                
                # Calculate mean and standard error
                mean_norm = np.mean(normalized_values)
                se_norm = stats.sem(normalized_values)  # Standard error
                
                normalized_means.append(mean_norm)
                normalized_errors.append(se_norm)
            
            # Plot line with error bars
            color_key = f'{tr}{day}' if f'{tr}{day}' in colors else tr
            
            ax.errorbar(x_positions, normalized_means, yerr=normalized_errors,
                        color=colors[color_key], marker='o', linewidth=2,
                        markersize=6, capsize=4, capthick=2,
                        label=f'{tr} {day}')
            
            # Print sample sizes
            valid_count = len(df_plot[(df_plot[param_1].notna()) & 
                                    (df_plot[param_1] > 0)])
            print(f"{tr} {day}: n={valid_count}")


    ax.set_xticks(ticks = x_positions, labels = [f'EPSP {x}' for x in x_positions])
    ax.set_ylabel(f'{titles_dict[DV][0]} \n ({titles_dict[DV][1]})')
    ax.set_yticks(ticks = titles_dict[DV][2], labels = titles_dict[DV][3])

    # Add horizontal line at y=1 for reference
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.7)

    # Save if directory provided
    if save_dir:
        plt.savefig(f'{save_dir}/normalized_{DV}_{data_type}.svg',
                    format='svg', bbox_inches='tight', dpi=300)

    plt.show()

def plot_connect_params(df, data_type, DV = 'amp',save_dir = False):

    h_cm = 6
    w_cm = 20
    titles_dict = pl_intr.dict_for_plotting_conn()

    if DV == 'amp':
        cols = ['Amp 1', 'Amp 2', 'Amp 3', 'Amp 4']
        param_1 = 'Amp 1'
    elif DV == 'lat':
        cols = ['Lat1', 'Lat2', 'Lat3', 'Lat4']
        param_1 = 'Lat1'

    colors = figure_colors()[data_type]

    fig, ax = plt.subplots(1, len(cols), figsize = (w_cm /2.54,h_cm/2.54), sharey = True)
    ax = ax.flatten()
    for a, param in enumerate(cols):
        
        label_ = []
        num_cels = {}
        
        for i, tr in enumerate(sorted(df.treatment.unique())):
            x_plot =[]
            for j, day in enumerate(sorted(df.day.unique())):
                k = j + 2.5*i
                df_plot = df[(df['treatment'] == tr) &\
                              (df['day'] == day)]
                # df_plot = df_plot.dropna(subset=cols)
                df_plot = df_plot[(df_plot[cols] > 0).all(axis=1)]
                
                x = np.linspace(0.65+k, 1.35+k, len(df_plot))
                x_plot.append(x)
                ax[a].scatter(x, df_plot[param], alpha = 0.3,
                            color = colors[tr + day], zorder = 5)
                ax[a].scatter(x, df_plot[param], alpha = 0.75,
                            edgecolor = colors[tr + day],
                            facecolors = 'none')
                if 'cell' not in data_type:
                    ax[a].boxplot(df_plot[param], positions = [k + 1], zorder = 10, widths = 0.5)

                    label_.append(day + tr)
                    num_cels[tr + ' ' + day] = len(df_plot)
                    #data_boxplot.append(df_plot[amp)

            if k in [1,3.5,5] and 'cell' in data_type:
                for c, cell in enumerate(df_plot['connection_ID']):
                    if k == 1 and len(df_plot['connection_ID'])-1 == c:
                        print(str(c+1) + 'Gencho')
                    x1 = [x_plot[0][c], x[c]]
                    y = df[param][df['connection_ID'] == cell]

                    ax[a].plot(x1, y, '-', color = colors[tr+'D2'], alpha = 0.5,
                               linewidth = 2, zorder = 0)

        
        ax[a].set_xticks(ticks = [1,2,3.5,4.5], labels = ['pre', 'post', 'pre', 'post'])

        ax[0].set_ylabel(titles_dict[cols[0]][0] + ' (' + titles_dict[cols[0]][1] + ')')
        ax[0].set_yticks(ticks = titles_dict[cols[0]][2], labels = titles_dict[cols[0]][3])
            
    if save_dir:
        plt.savefig(f'{save_dir}/con_all_{DV}_{data_type}.svg',
                    format='svg', bbox_inches='tight', dpi=300)
    plt.show(fig)


def plot_connect_1st_amp_lat_only(df, data_type, save_dir = False):

    h_cm = 6
    w_cm = 12
    titles_dict = pl_intr.dict_for_plotting_conn()

    cols = ['Amp 1', 'Lat1']
    colors = figure_colors()[data_type]

    fig, ax = plt.subplots(1, 2, figsize = (w_cm /2.54,h_cm/2.54))
    ax = ax.flatten()
    for a, param in enumerate(cols):
        if 'cell' in data_type:
        # df_plot = df_plot.dropna(subset = cols)
            # Get connection_IDs from rows that would be dropped due to NaN in specific columns
            nan_mask = df[param].isnull()# .any(axis=1)
            con_ids_drop = df[nan_mask]['connection_ID'].unique()
            # Remove all rows that have those connection_IDs
            df = df[~df['connection_ID'].isin(con_ids_drop)]
        
        label_ = []
        num_cels = {}
        
        for i, tr in enumerate(sorted(df.treatment.unique())):
            x_plot =[]
            for j, day in enumerate(sorted(df.day.unique())):
                k = j + 2.5*i
                df_plot = df[(df['treatment'] == tr) &\
                              (df['day'] == day)]
                # df_plot = df_plot.dropna(subset=cols)
                # df_plot = df_plot[(df_plot[cols] > 0).all(axis=1)]
                
                x = np.linspace(0.65+k, 1.35+k, len(df_plot))
                x_plot.append(x)
                ax[a].scatter(x, df_plot[param], alpha = 0.3,
                            color = colors[tr + day], zorder = 5)
                ax[a].scatter(x, df_plot[param], alpha = 0.75,
                            edgecolor = colors[tr + day],
                            facecolors = 'none')
                
                if 'cell' not in data_type:
                    ax[a].boxplot(df_plot[param], positions = [k + 1], zorder = 10, widths = 0.5)

                    label_.append(day + tr)
                    num_cels[tr + ' ' + day] = len(df_plot)
                    #data_boxplot.append(df_plot[amp)
                else:
                    ax[a].plot(k+1, df_plot[param].median(), zorder = 10, marker = '_',
                               markersize = 15, c = 'k')

            if k in [1,3.5,5] and 'cell' in data_type:
                for c, cell in enumerate(df_plot['connection_ID']):
                    if k == 1 and len(df_plot['connection_ID'])-1 == c:
                        print(str(c+1) + 'Gencho')
                    x1 = [x_plot[0][c], x[c]]
                    y = df[param][df['connection_ID'] == cell]

                    ax[a].plot(x1, y, '-', color = colors[tr+'D2'], alpha = 0.5,
                               linewidth = 2, zorder = 0)

        ax[a].set_xticks(ticks = [1,2,3.5,4.5], labels = ['pre', 'post', 'pre', 'post'])

        ax[a].set_ylabel(titles_dict[param][0] + ' (' + titles_dict[param][1] + ')')
        ax[a].set_yticks(ticks = titles_dict[param][2], labels = titles_dict[param][3])
            
    if save_dir:
        plt.savefig(f'{save_dir}/con_params_{data_type}_amp1_lat1.svg',
                    format='svg', bbox_inches='tight', dpi=300)
    plt.show(fig)


def plot_amp_lat_cumulative(df, data_type, save_dir = False):
    w_cm = 12
    h_cm = 6

    titles_dict = pl_intr.dict_for_plotting_conn()
    colors_ = figure_colors()

    params = ['Amp 1', 'Lat1']
    fig, ax = plt.subplots(1, 2, figsize = (w_cm/2.54, h_cm/2.54))
    for u, param in enumerate(params):
        for day in ['D1', 'D2']:
            for tr in ['Ctrl', 'high K']:

                col = colors_[data_type][tr + day]
                if day == 'D1':
                    data_plot = df[(df.day == day)]
                else:
                    data_plot = df[(df.day == day) &
                            (df.treatment == tr)]
                data_plot = data_plot.dropna(subset=param)

                if len(data_plot[param]) == 0:
                    continue

                limit_x = np.max(data_plot[param])

                # y = np.arange(1, len(data_plot) +1 ) / len(data_plot)
                # ax[u].plot(data_plot[param], y, c = col)
                ax[u].ecdf(data_plot[param], c = col, linewidth = 2.5, label = tr)
                # ax[u].text(50, 0.1 * u, f'{tr+day} n = {len(data_plot)}', c = col)
                ax[u].set_title(data_type)
                if param == 'Amp 1':
                    ax[u].set_xticks(ticks = titles_dict[param][2])
                else:
                    ax[u].set_xticks(ticks = np.linspace(0, 5, 5))
                ax[u].set_xlabel(f'{titles_dict[param][0]} ({titles_dict[param][1]})')

    # plt.figlegend(loc = 'upper right',  bbox_to_anchor=(1, 1))

    if save_dir:
        plt.savefig(f'{save_dir}{data_type}_cum_amps_lat.svg', \
                    format = 'svg', bbox_inches = 'tight', dpi = 300)
        
def plot_connection_ratio(df_con_percentage_all, data_type='slice', save_dir=False):
    """
    Plot the ratio of found connections to all possible connections for each treatment-day combination
    """
    w_cm = 8
    h_cm = 6
    fig, ax = plt.subplots(1, 1, figsize=(w_cm/2.54, h_cm/2.54))
    colors = figure_colors()[data_type]
    
    for j, day in enumerate(sorted(df_con_percentage_all.day.unique())):
        for i, tr in enumerate(sorted(df_con_percentage_all.treatment.unique())):
            k_ = j + 2.5*i
            df_plot = df_con_percentage_all[(df_con_percentage_all.treatment == tr) & 
                                          (df_con_percentage_all.day == day)]
            
            if len(df_plot) == 0:
                continue
            
            # Calculate overall ratio for this treatment-day combination
            total_found = df_plot.found_connections.sum()
            total_possible = df_plot.num_possible_connections.sum()
            overall_ratio = total_found / total_possible
                
            # Plot the overall ratio as a single point
            ax.scatter(k_+1, overall_ratio, c=colors[tr+day], s=100, alpha=0.8)
            
            # Add text with ratio and counts
            ax.text(k_+0.85, 1.1*overall_ratio, f'{overall_ratio:.3f}\n({total_found}/{total_possible})', 
                   size=9, ha='center')
    
    ax.set_xticks(ticks=[1, 2, 3.5, 4.5], labels=['pre', 'post', 'pre', 'post'])
    ax.set_ylabel('Connection ratio\n(Found / All possible)')
    ax.set_ylim(0, ax.get_ylim()[1] * 1.2)  # Add space for text
    
    # Add treatment labels
    ax.text(1.5, ax.get_ylim()[1] * 0.9, 'CTR', ha='center', fontsize=12, 
           color=colors['CtrlD2'])
    ax.text(4, ax.get_ylim()[1] * 0.9, 'HiK', ha='center', fontsize=12, 
           color=colors['high KD2'])
    
    if save_dir:
        plt.savefig(f'{save_dir}{data_type}_connection_ratio.svg',
                   format='svg', bbox_inches='tight', dpi=300)
    plt.show()

def plot_connection_ratio2(df_con_percentage_all, data_type='slice', save_dir=False):
    """
    Plot the ratio of found connections to all possible connections for each treatment-day combination
    """
    w_cm = 8
    h_cm = 6
    fig, ax = plt.subplots(1, 1, figsize=(w_cm/2.54, h_cm/2.54))
    colors = figure_colors()[data_type]
    
    # Group by treatment and day, then sum the counts
    grouped = df_con_percentage_all.groupby(['treatment', 'day']).agg({
        'found_connections': 'sum',
        'num_possible_connections': 'sum'
    }).reset_index()
    
    # Calculate ratio for each group
    grouped['ratio'] = grouped['found_connections'] / grouped['num_possible_connections']
    
    for j, day in enumerate(sorted(grouped.day.unique())):
        for i, tr in enumerate(sorted(grouped.treatment.unique())):
            k_ = j + 2.5*i
            row = grouped[(grouped.treatment == tr) & (grouped.day == day)]
            
            if len(row) == 0:
                continue
            
            ratio = row['ratio'].iloc[0]
            found = row['found_connections'].iloc[0]
            possible = row['num_possible_connections'].iloc[0]
                
            # Plot the ratio as a single point
            ax.scatter(k_+1, ratio, c=colors[tr+day], s=100, alpha=0.8)
            
            # Add text with ratio and counts
            ax.text(k_+0.85, 1.1*ratio, f'{ratio:.3f}\n({found}/{possible})', 
                   size=9, ha='center')
    
    ax.set_xticks(ticks=[1, 2, 3.5, 4.5], labels=['pre', 'post', 'pre', 'post'])
    ax.set_ylabel('Connection ratio\n(Found / All possible)')
    ax.set_ylim(0, ax.get_ylim()[1] * 1.2)  # Add space for text
    
    # Add treatment labels
    ax.text(1.5, ax.get_ylim()[1] * 0.9, 'CTR', ha='center', fontsize=12, 
           color=colors['CtrlD2'])
    ax.text(4, ax.get_ylim()[1] * 0.9, 'HiK', ha='center', fontsize=12, 
           color=colors['high KD2'])
    
    if save_dir:
        plt.savefig(f'{save_dir}{data_type}_connection_ratio.svg',
                   format='svg', bbox_inches='tight', dpi=300)
    plt.show()

def plot_repatch_char_D1_D2(OP, patcher, file_index, active_channels, save_dir = False,
        human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    '''
    plots the characterization on D1 and D2
    plots inj 0, 5 (baseline), and the second sweep with more than 1 AP
    '''
    w_cm = h_cm = len(active_channels[0]) * 3
    clrs = [['#377eb8', '#ff7f00', '#4daf4a',
            '#984ea3', '#f781bf', '#999999', 
            '#e41a1c', '#dede00'],
             ['#377eb8', '#ff7f00', '#4daf4a',
               '#984ea3',  '#f781bf']]
    y_labels_D2 = [f're Ch{ch}' for ch in active_channels[0]]

    work_dir, filenames, _, _, _, _ = sort.get_OP_metadata(human_dir, OP, patcher)
    fn = work_dir + filenames[file_index[0]]
    fn2 = work_dir + filenames[file_index[1]]
    charact_data = [hcf.load_traces(fn), hcf.load_traces(fn2)]

    fig, ax = plt.subplots(len(active_channels[0]), 2, sharex=True, sharey=True,\
                            figsize = (w_cm/2.54, h_cm / 2.54))
    for f, data in enumerate(charact_data):
        for i, ch in enumerate(active_channels[f]):
            ch1 = data[f'Ch{ch}'][0]

            first_col = np.where((ch1 > 0).any(axis=0))[0][0] + 1
            # last_indx = inj.index(rheos[0])
            # plt_swps = np.concatenate([[0], np.arange(5, last_indx, 3)])
            for h in [0, 5, first_col]:
                ax[i, f].plot(ch1[0:25000,h], c = clrs[f][i], alpha = 0.75) # lw = 2
            
            if f == 0:
                ax[i, f].set_ylabel(f'Ch{ch}')
            else:
                ax[i, f].set_ylabel(y_labels_D2[i])
            ax[i, f].set_xticks([])
            ax[i, f].set_yticks([])
            ax[i, f].spines['bottom'].set_visible(False)
            ax[i, f].spines['left'].set_visible(False)
            # x[i, f].set_title(f'Day {i+1}, {rheos[0]} pA')
    if save_dir:
        plt.savefig(f'{save_dir}char_trace_{OP}.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 500)
    plt.show()

def plot_repatch_connect_D1_D2(OP, patcher, file_index, active_channels, save_dir = False,
        human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    w_cm = len(active_channels[0]) * 4.5
    h_cm = len(active_channels[0]) * 3
    z1 = 1.5
    z2 = 40.5
    stim_window = stim_win.stim_window_con_screen
    clrs = [['#377eb8', '#ff7f00', '#4daf4a',
            '#984ea3', '#f781bf', '#999999', 
            '#e41a1c', '#dede00'],
             ['#377eb8', '#ff7f00', '#4daf4a',
               '#984ea3',  '#f781bf']]

    y_labels_D2 = [f're Ch{ch}' for ch in active_channels[0]]

    work_dir, filenames, _, _, _, _ = sort.get_OP_metadata(human_dir, OP, patcher)
    fn = work_dir + filenames[file_index[0]]
    fn2 = work_dir + filenames[file_index[1]]
    con_screen_data = [hcf.load_traces(fn), hcf.load_traces(fn2)]

    fig, ax = plt.subplots(len(active_channels[0]), 2 * len(active_channels[0]), sharex = True, sharey = False,\
                            figsize = (w_cm/2.54, h_cm / 2.54))
    for f, data in enumerate(con_screen_data):
        for i, ch1 in enumerate(active_channels[f]):
            ch1_data = data[f'Ch{ch1}'][0]
            avg = np.mean(ch1_data, axis = 1)

            for j, ch2 in enumerate(active_channels[0]):
                if f==1:
                    j = j +len(active_channels[0])
                plotwin = avg[stim_window[f'Ch{ch2}'][0]:stim_window[f'Ch{ch2}'][1]] #from the average takes the signal for this stim_window
                ax[i,j].plot(plotwin, c = clrs[f][i])
                
                if j == 0:
                    ax[i, j].set_ylabel(f'Ch{ch1}')
                if j == len(active_channels[f]):
                    ax[i, j].set_ylabel(y_labels_D2[i])
                ax[i, j].yaxis.label.set_color(clrs[f][i])
                ax[i, j].set_xticks([])
                ax[i, j].set_yticks([])
                ax[i, j].spines['bottom'].set_visible(False)
                ax[i, j].spines['left'].set_visible(False)

                if plotwin.max()-plotwin.min() < 10: 
                    ax[i,j].set_ylim([plotwin.min() - z1, plotwin.max() + z1])
                    v1 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z1, color='k') # lw=0.7
                else:
                    ax[i,j].set_ylim([plotwin.min() - z1, plotwin.max() + z2])
                    v2 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z2, color='k')

    if save_dir:
        plt.savefig(f'{save_dir}connect_trace_{OP}.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 500)
    plt.show()

def  plot_connection_window(OP, patcher, labels, fn_indx, pre_chans, post_chans, save_dir = False,
                            human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    '''
    
    '''
    h_cm = 3.5 * len(pre_chans)
    w_cm = 14

    work_dir, filenames, _, _, _, _ = sort.get_OP_metadata(human_dir, OP, patcher)

    fig, ax = plt.subplots(2, len(fn_indx), sharex = True, # sharey = 'row',
                           figsize = ((w_cm/2.54, h_cm / 2.54)))

    for f, indx in enumerate(fn_indx):
        con_screen_file = f'{work_dir}{filenames[indx]}'
        pre_ch = pre_chans[f]
        post_ch = post_chans[f]
        # getting needed data
        pre_signal = con_param.presynaptic_screen_ALLOW_ALL(con_screen_file, pre_ch)
        post_signal, _ = con_param.postsynaptic_screen_ALLOW_ALL(con_screen_file, post_ch, [])
        _, _, pre_window, post_window, preAPs_shifted, preAPs = con_param.get_analysis_window(pre_signal, post_signal)
        post_peaks = con_param.find_postsynaptic_peaks(post_window, preAPs_shifted)
        bl = con_param.get_psp_baselines(post_window,preAPs_shifted)
        # onsets = con_param.get_onsets(preAPs_shifted, post_window, post_peaks, bl)
        sampl_rate, units, times = hcf.get_abf_info(con_screen_file, pre_ch, np.shape(post_signal)[1], np.shape(post_signal)[0])

        my_labels = {'l1' : 'peak pre_AP', 'l2' : 'baseline', 'l3': 'onset', 'l4': 'post synnaptic peak'}
        ax[0, f].plot(times[:len(pre_window)], pre_window, color='k')   #plost pre APs 0 shifted
        
        for i in range(0,len(preAPs_shifted[0])):
            ax[0, f].scatter(preAPs_shifted[0][i]/sampl_rate,\
                pre_window[preAPs_shifted[0][i]], marker='o', color='r', label = my_labels['l1'])
            my_labels['l1'] = "_nolegend_"

        for i in range(np.shape(post_signal)[1]):
            y = post_signal[:,i][preAPs[0][0]-750:preAPs[0][len(preAPs[0])-1]+750]
            ax[1, f].plot(times[0:len(y)], y, lw=0.2, color='grey', alpha=0.4, zorder=0 )
            #ax[1].plot(y, lw=0.2, color='grey', alpha=0.4)
        
        y = post_window
        ax[1, f].plot(times[0:len(y)], post_window, color='k', lw=0.5, zorder=10) #post_window is the averaged 0 shifter post signal 
        ax[1, f].set_xlabel('Time (sec)')
        ax[0, f].set_yticks(ticks = np.linspace(-75, 55, 5), labels = [-75, '', -10, '', 55])
        ax[0, f].set_ylabel(str(units)[-2:])
        ax[0, f].set_title(f'{labels[f]} Pre ch{pre_ch}')
        ax[1, f].set_title(f'Post ch{post_ch}')

        # add legend markers
        # for i in range(0,len(preAPs_shifted[0])):
        #     sample_rate = sampl_rate.item()
        #     ax[f, 1].hlines(bl[i,0], (preAPs_shifted[0][i]-170)/sample_rate, (preAPs_shifted[0][i]-30)/sample_rate, 
        #     linestyles = 'solid', lw = 3, color='b', label = my_labels['l2'],zorder=10, alpha = 0.6)
        #     my_labels['l2'] = "_nolegend_"

        # ax[1, f].scatter(onsets/sampl_rate, bl, marker='^', color='r', label = my_labels['l3'],  zorder=10)

        # for i in range(0,len(post_peaks)):
        #     ax[1, f].scatter(post_peaks[i][1]/sampl_rate, post_peaks[i][0], marker='+',\
        #         color='g', s=30, linewidth=10, label = my_labels['l4'],zorder=10)
        #     my_labels['l4'] = "_nolegend_"

    # ax[1, 0].set_title('Post cell responsees')
    ax[1, 1].set_xticks(ticks = np.linspace(0, 0.24, 5), labels = [0, '', 1.2, '', 2.4])

    # plt.figlegend(loc = 'center right',  bbox_to_anchor=(0.95, 1))

    if save_dir:
        plt.savefig(f'{save_dir}{OP}-repatch_con_trace.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 500)
    plt.show()

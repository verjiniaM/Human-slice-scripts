import os
import datetime
import math
from scipy import stats
import pyabf
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ephys_analysis.funcs_con_screen as con_param
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
import ephys_analysis.funcs_human_characterisation as hcf
import ephys_analysis.funcs_sorting as sort
import ephys_analysis.stimulation_windows_ms as stim_win

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
                    'post CTR': "#d383cf",
                    'TTXD1': 'darkgrey',
                    'TTXD2': "#0e62ff"},
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
                    {'Ctrl': '#496B7F'},
                    'minis':{
                    'Ctrl': "#570b77",
                    'high K':  "#1bd7e0",
                    'correlation': {'CTR': "#0e62ff"}}}
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
                    'Rin': (0,2),
                    'resting_potential': (0,3),
                    'rheo_ramp_c' : (1,1),
                    'TH': (1,0),
                    'membra_time_constant_tau': (0,4),
                    'sag': (0,5)},
                    'inc_only_2':{
                    'Rin': (0,0),
                    'TH': (0,1),
                    'AP_heigth': (0,2),
                    'sag': (0,3)},
                    'inc_only_su':{
                    'resting_potential': (0,2),
                    'membra_time_constant_tau': (0,3),
                    'rheo_ramp_c' : (1,0),
                    'max_depol': (1,1),
                    'AP_halfwidth': (1,2),
                    'max_repol': (1,3)},
                'membr':{'Rin': (0,0),
                    'resting_potential': (0,1),
                    'membra_time_constant_tau': (0,2),
                    'sag': (0,3)},
                'AP':{'AP_heigth': (0,0),
                    'AP_halfwidth': (0,1),
                    'max_depol': (0,2),
                    'max_repol': (0,3),
                    'TH': (1,1),
                    'rheo_ramp_c' : (1,0),
                    'Rs': (1,2)},
                'other':{
                    'resting_potential': (0,2),
                    'rheo_ramp_c' : (1,0),
                    'Rin': (0,0),
                    'sag': (1,2),
                    'TH': (1,1),
                    'Rs':(0,1)},
                  'slice_all_CTR':{
                      'Rin': (0,1),
                      'resting_potential': (0,2),
                      'TH': (1,0),
                      'rheo_ramp_c' : (1,1),
                      'sag': (1,2)}}

    return order_dict

def intr_params_fig2(df, fig_type, dest_dir, w_cm = 18, h_cm = 8):
    """
    Plots parameters for different treatments and days 
    fig_type: 'other', 'AP', 'membr'
    """
    # for single plots - w_cm = 6.5, h_cm = 5.7
    mask_exclude = (df['treatment'] != 'high K') & (df['treatment'] != 'Ctrl')
    df = df[~mask_exclude]
    # params = ['Rin', 'resting_potential', 'rheo_ramp_c', 'TH', 'sag']
    params = ['sag', 'AP_heigth', 'max_depol', 'max_repol',
              'membra_time_constant_tau', 'AP_halfwidth',
              'Rin', 'resting_potential', 'rheo_ramp_c', 'TH']
    params = ['Rs']
    ax_elements = fig_order()[fig_type]
    df = df.sort_values(['treatment', 'day'])
    titles_all = pl_intr.dict_for_plotting()
    titles_dict = {key: titles_all[key] for key in  ax_elements.keys() if key in titles_all}

    dict_label = {'Ctrl':'CTR', 'high K':'HiK'}
    colors = figure_colors()['slice']
    fig2, axs = plt.subplots(2,4, figsize = (w_cm/2.54, h_cm / 2.54))

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
                if param == 'Rs':
                    print(f'Avg. for {tr} {day} {np.round(np.mean(df_plot[param]), 2)}')
                label_.append(day + '\n' + dict_label[tr])
                num_cels[tr + ' ' + day] = len(df_plot)
                data_boxplot.append(df_plot[param])

        axs[ax_el].set_xticks(ticks = [1,2,3.5,4.5], labels = ['pre', 'post', 'pre', 'post'])

        axs[ax_el].set_ylabel(item[0] +'\n(' + item[1] + ')')
        axs[ax_el].set_yticks(ticks = item[2], labels = item[3])
        # axs[ax_el].set_title('CTR                            HiK')

    date = str(datetime.date.today())
    os.makedirs(dest_dir, exist_ok = True)
    plt.savefig(f'{dest_dir}{date}{fig_type}_V2_plot.svg', \
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

    df_summary = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results'+ \
                        '/human/paper_figs_collected_checked/OLD___stats/archive/'+ \
                        'sum_data_CIs_SEs_for_plotting/firing_plot_CIs_temporal.xlsx')

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
    # fig, ax = plt.subplots(1, 2, figsize=(w_cm / 2.54, h_cm/2.54),\
    #                     sharex = True)

    counts_firing = {}
    for k, treatment in enumerate(iff_df.treatment.unique()):
    # for k, treatment in enumerate(['Ctrl', 'Ctrl']):
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
                    sems.append(sem)
                
                counts.append(len(data))
                avgs.append(avg)
                inj.append(inj_val)

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

                # ax[0].set_ylabel(dv_dict[dv][1])

                ## TO DO which lines does this have to be on
                # ax[k].fill_between(range(1, len(inj)+1), np.array(avgs) + np.array(sems), \
                #                     np.array(avgs) - np.array(sems), color = colors[treatment + day], alpha = 0.3)
           
            # ax[k].plot(range(1, len(inj)+1), avgs, label = day, color = colors[treatment + day])
            # counts_firing[treatment + day] = counts

            # if dv == 'num_aps':
            #     ax[k].set_yticks(ticks = [0, 7.5, 15, 22.5, 30], labels = [0, '', 15, '', 30])
            # else:
            #     ax[k].set_yticks(ticks = np.linspace(0, 100, 5), labels = [0, '', 50, '', 100])

    ticks_ = np.arange(1, 14, 2)
    labels_ = [inj[a] for a in np.arange(0,  13, 2)]
    ax[0].set_xticks(ticks = ticks_, labels = labels_)
    ax[1].set_xticks(ticks = ticks_, labels = labels_)
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

def intr_params_fig3(df, dest_dir, w_cm = 18, h_cm = 8):
    '''
    plots the repatch data
    plots means with line and value
    no coloring based on age
    could be added but plot become too much
    '''
    # params = ['Rin', 'resting_potential', 'rheo_ramp_c', 'TH', 'sag']
    params = ['max_depol', 'max_repol',
              'AP_halfwidth', 'AP_heigth']
    params = ['Rs']
    ax_elements = fig_order()['AP']
    df = df.sort_values(['treatment', 'day'])
    titles_all = pl_intr.dict_for_plotting()
    titles_dict = {key: titles_all[key] for key in params if key in titles_all}
    colors = figure_colors()['cell_ID_new']
    dict_label = {'Ctrl':'CTR', 'high K':'HiK'}

    fig2, axs = plt.subplots(2,4, figsize = (w_cm/2.54, h_cm / 2.54))

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
                mean = np.nanmean(df_plot[param])

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
    plt.savefig(dest_dir + date + '_plot_all_repatch_MEAN_sup.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.show(fig2)

def firing_props_fig3(iff_df, error_type, dest_dir,
                 background_dot_size = 3, w_cm = 18, h_cm = 15):
    """
    Plots the initial firing frequency and number of APs
    """
    df_summary = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results'+ \
                           '/human/paper_figs_collected_checked/OLD___stats/archive/'+ \
                            'sum_data_CIs_SEs_for_plotting/firing_plot_CIs_temporal.xlsx')
    
    iff_df = iff_df.sort_values(['treatment', 'day'])
    num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(iff_df)
    dt_sum = 'slice'
    colors = figure_colors()[dt_sum]

    dv_dict = {'IFF' : [iff_indx, 'Initial firing frequency',
                'Initial firing \n frequency (Hz)'],
                'num_aps' : [num_aps_indx, 'bluff', 'Firing frequency \n (Hz)']}

    fig, ax = plt.subplots(4, 4, figsize=(w_cm / 2.54, h_cm/2.54),\
                            sharex = False) # gridspec_kw = {'height_ratios': [5, 5, 1, 1]})

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
                        ax[v,k].set_yticks(ticks = np.linspace(0, 30, 5), labels = [0, '', 15, '', 30])
                    else:
                        ax[v,k].set_yticks(ticks = np.linspace(0, 100, 5), labels = [0, '', 50, '', 100])
                else:
                    if dv == 'num_aps':
                        ax[v,k].set_yticks(ticks = np.linspace(0, 50, 5), labels = [0, '', 25, '', 50])
                    else:
                        ax[v,k].set_yticks(ticks = [0, 40, 80, 120, 160, 200])

                ax[v,0].set_ylabel(dv_dict[dv][2])

        if dv == 'IFF':
            # bar graph after counting firing cells
            x = np.arange(1, 14, 1)
            ax[2,0].bar(x, counts_firing['CtrlD1'], color = colors['CtrlD1'])
            ax[3,0].bar(x, counts_firing['CtrlD2'], color = colors['CtrlD2'])
            # ax[2,1].bar(x, counts_firing['high KD1'], color = colors['high KD1'])
            # ax[3,1].bar(x, counts_firing['high KD2'], color = colors['high KD2'])

        ticks_ = np.arange(1, 14, 2)
        labels_ = [inj[a] for a in np.arange(0,  13, 2)]

        # top plot
        # ax[0,0].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
        # ax[0,1].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
        
        ax[0,0].set_xticks(ticks_)
        ax[0,0].set_xticklabels(labels_, rotation=30)
        ax[0,1].set_xticks(ticks_)
        ax[0,1].set_xticklabels(labels_, rotation=30)
        
        ax[v,k].set_xticks(ticks_)
        ax[v,k].set_xticklabels(labels_, rotation=30)
        ax[v,k].set_xticks(ticks_)
        ax[v,k].set_xticklabels(labels_, rotation=30)


        # bottom
        ax[3,0].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
        ax[3,1].set_xticks(ticks = ticks_, labels = labels_, rotation = 30)
        ax[3,0].set_xlabel('Current injection (pA)')
        ax[3,0].set_ylabel('Firing cells count')
        # ax[0,0].set_title('CTR', color = colors['CtrlD2'])
        # ax[0,1].set_title('HiK', color = colors['high KD2'])

    plt.subplots_adjust(hspace = 0.12, wspace = 0.15)
    plt.savefig(dest_dir + error_type + dv + f'_vs_curernt_inj_{dt_sum}.svg', \
                format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.show(fig)

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

def intr_params_inc_only(df, dest_dir, w_cm = 22, h_cm = 8):
    """
    Plots parameters for different treatments and days 
    """

    mask_exclude = (df['treatment'] != 'high K') & (df['treatment'] != 'Ctrl')
    df = df[~mask_exclude]
    params = ['rheo_ramp_c', 'TH', 'Rin', 'resting_potential', 'membra_time_constant_tau', 'sag']
    params = ['AP_heigth','AP_halfwidth', 'max_depol', 'max_repol', 'rheo_ramp_c', 'TH']
    ax_elements = fig_order()['AP'] #'AP'
    df = df.sort_values(['treatment', 'day'])
    titles_all = pl_intr.dict_for_plotting()
    titles_dict = {key: titles_all[key] for key in params if key in titles_all}

    dict_label = {'Ctrl':'CTR', 'high K':'HiK'}
    colors = figure_colors()['inc_only']

    fig2, axs = plt.subplots(2,6, figsize = (w_cm/2.54, h_cm / 2.54))
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
    plt.savefig(dest_dir  + date + 'inc_only_sub_' + '.svg', \
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
                 background_dot_size = 3, w_cm = 6, h_cm = 4):
    """
    Plots the initial firing frequency and number of APs for only incubated slices
    """

    df_summary = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results'+ \
                        '/human/paper_figs_collected_checked/OLD___stats/archive/'+ \
                        'sum_data_CIs_SEs_for_plotting/22.firing_plot_CIs_temporal_inc_only.xlsx')

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

                    counts.append(len(data))
                    avgs.append(avg)
                    inj.append(inj_val)

                    counts_firing[treatment] = counts


                    if error_type == 'CI':
                        sum_data_tr = df_summary[(df_summary['var_firing'] == dv) & \
                                        (df_summary['treatment'] == treatment) & \
                                        (df_summary['inj_current'] == inj_val)]

                        start_summary = len(inj) - len(sum_data_tr) + 1
        
                        ax.scatter(range(1, len(inj)+1), avgs, color = colors[treatment])
                        ax.fill_between(range(start_summary, len(inj)+1), sum_data_tr['CI_L'],
                                            sum_data_tr['CI_U'], color = colors[treatment],alpha = 0.3)

                ax.plot(range(1, len(inj)+1), avgs, color = colors[treatment])
                ax.scatter(range(1, len(inj)+1), avgs, color = colors[treatment])

                if error_type == 'CI':
                    sum_data_tr = df_summary[(df_summary['var_firing'] == dv) & \
                                    (df_summary['treatment'] == treatment)]

                    start_summary = len(inj) - len(sum_data_tr) + 1
      
                    ax.scatter(range(1, len(inj)+1), avgs, color = colors[treatment])
                    ax.fill_between(range(start_summary, len(inj)+1), sum_data_tr['CI_L'],
                                        sum_data_tr['CI_U'], color = colors[treatment],alpha = 0.3)
                    
                    if dv == 'num_aps':
                        ax.set_yticks(ticks = np.linspace(0, 30, 5))
                    else:
                        ax.set_yticks(ticks = np.linspace(0, 60, 5))
                else:
                    if dv == 'num_aps':
                        ax.set_yticks(ticks = np.linspace(0, 30, 5))
                    else:
                        ax.set_yticks(ticks = np.linspace(0, 600, 5))
                    ax.fill_between(range(1, len(inj)+1), np.array(avgs) + np.array(sems), \
                        np.array(avgs) - np.array(sems), color = colors[treatment], alpha = 0.3)
                    ax.plot(range(1, len(inj)+1), avgs, label = treatment, color = colors[treatment])

                    counts_firing[treatment] = counts

        ticks_ = np.arange(1, 14, 2)
        labels_ = [inj[a] for a in np.arange(0,  13, 2)]
       
        ax.set_ylabel(dv_dict[dv][2])
        ax.set_xticks(ticks_[::2])
        ax.set_xticklabels(labels_[::2])
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
    removes the high K D2 rows.
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
    w_cm =  18 # 6 * len(params_list)
    h_cm = 4
    
    titles_dict = pl_intr.dict_for_plotting()
    for df_name, df_info in df_list.items():

        df = df_info[0]
        col = df_info[1]

        df_cor = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'cor_estimate')
        df_p_val = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'p_vals')
        df_p_val = df_p_val.iloc[0:11]
        df_l_CI = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'lower_CI')
        df_u_CI = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'upper_CI')

        df_cor.index = [df_cor.columns[1:]]
        df_p_val.index = [df_p_val.columns[1:]]
        df_l_CI.index = [df_l_CI.columns[1:]]
        df_u_CI.index = [df_u_CI.columns[1:]]

        # cor.loc['hrs_after_OP', 'hrs_after_OP']
        # params_list = df_u_CI.columns.to_list()
        params_list.remove(x)
        # params_list.remove('Unnamed: 0')

        # fig, ax = plt.subplots(1, len(params_list), figsize = (w_cm/2.54, h_cm/2.54))
        fig, ax = plt.subplots(1, 4, figsize = (w_cm/2.54, h_cm/2.54))
        ax = ax.flatten()
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
            print(f'For {param}: r = {cor:.3f} \n p = {p_val:.3f}')

            # Set labels and title
            ax[i].set_xlabel(f'{titles_dict[x][0]} ({titles_dict[x][1]})')
            ax[i].set_xticks(ticks = titles_dict[x][2], labels = [int(x) for x in titles_dict[x][2]])
            if param == 'AP_halfwidth':
                ax[i].set_ylabel(f'{titles_dict[param][0]} \n ({titles_dict[param][1]})')
                ax[i].set_yticks(ticks = [0, 0.65, 1.3, 1.95, 2.6], labels = [0, '', 1.3, '', 2.6])
            else:
                ax[i].set_ylabel(f'{titles_dict[param][0]} \n ({titles_dict[param][1]})')
                ax[i].set_yticks(ticks = titles_dict[param][2], labels = titles_dict[param][3])
            #ax[i].set_title(f'{param}')
        
        plt.savefig(f'{save_dir_figs}{df_name}2corr_plots.svg',
                    format='svg', bbox_inches='tight', dpi=300)
        plt.show()

def plot_example_firing_corr(data_type, df, dest_dir, cell_IDs):
    '''
    hi
    '''
    w_cm = 20
    h_cm = 9
    work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/'

    swps = [9, 13, 17]
    # colors = trace_colors_trace()[data_type]

    colors = figure_colors()[data_type]
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

def plot_non_firing_across_inj_count(df, prop_or_count, data_type, save_dir):
    '''
    prop_or_count: refers to the proportion of non_firing to all
    '''
    w_cm = 9
    h_cm = 4
    colors = figure_colors()[data_type]

    if prop_or_count == 'prop':
        param = 'prop_non_firing'
        y_axis = 'Non-firing cells \n (proportion)'
        y_ticks = np.linspace(0, 1, 5)
        y_labels = [0, '', 0.5, '', 1]
        col_ = (df['count_non_firing_cells'] / (df['count_non_firing_cells'] + df['count_firing_cells'])).values
        df[param] = col_
    else:
        param = 'count_non_firing_cells'
        y_axis = 'Non-firing cells \n (count)'
        y_ticks = np.linspace(0, 70, 5) # (0, 36, 5)
        y_labels = [0, '', 35, '', 70] #  [0, '', 18, '', 36]

    fig, ax = plt.subplots(1, 2, figsize = (w_cm/2.54, h_cm/2.54))
    x = df.injection.unique()[5:18]
    for i, tr in enumerate(sorted(df.treatment.unique())):
        for day in ['D1', 'D2']:
            if data_type != 'inc_only':
                y = df[param][(df.treatment == tr) &
                                                    (df.day == day)].values[5:18]

                ax[i].plot(x, y, color = colors[tr + day], alpha = 0.75)
                ax[i].scatter(x, y, color = colors[tr + day], alpha = 0.75)
                ax[i].set_ylabel(y_axis)
                ax[i].set_yticks(ticks = y_ticks, labels = y_labels)
            else:
                y = df[param][(df.treatment == tr)].values[5:18]
                ax[0].plot(x, y, color = colors[tr], alpha = 0.75)
                ax[0].scatter(x, y, color = colors[tr], alpha = 0.85)
                ax[0].set_ylabel(y_axis)
                ax[0].set_yticks(ticks = y_ticks, labels = y_labels)

    ax[0].set_xlabel('Current injection (pA)')
    # labels_ = []
    # for l, val in enumerate(x):
    #     if l % 2 == 0:
    #         labels_.append(val)
    #     else:
    #         labels_.append('')
    ax[0].set_xticks(ticks  = np.linspace(0, 600, 13), 
    labels = [0, '', '', '', 200, '', '', '', 400, '', '', '',600])
    # ax[1].set_xticks(ticks  = x[::4], labels = labels_[::4])

    plt.savefig(f'{save_dir}{data_type}{prop_or_count}_non_firing_cells.svg',
            format='svg', bbox_inches='tight', dpi=300)

def plot_corrs_firing(x, df_list, dir_read_corrs, save_dir_figs, params_list = []):
    '''
    x - x_axis parameter
    df_list - dictionary. 'df_name' = [df, 'color']
    save_dir_corr - where the correlations dataframe is saved
    '''
    w_cm = 20
    h_cm = 11
    
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
                ticks_y = [0, 12.5, 25, 37.5, 50]
                labels_y = [0, '', 25, '', 50]
            else:
                col_name = f"('{param[3:6]}pA', 'IFF')"
                i = 1
                j = int(p/2 - 0.5)
                ticks_y = [0, 50, 100, 150, 200]
                labels_y = [0, '', 100, '', 200]
            
            if 'ALL' not in df_name:
                # take only firing cells
                df_plot = df[df[col_name].notna() & (df[col_name] != 0)]
            
            sns.regplot(data=df_plot, x=df_plot[x], y=df_plot[col_name], ax = ax[i, j], color = col, scatter = False)
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
            ax[i, j].set_ylabel(f'{param} \n ()')
            ax[i, j].set_yticks(ticks = ticks_y, labels = labels_y)
            ax[i, j].set_title(f'{param[3:6]} pA')
        
        plt.savefig(f'{save_dir_figs}{df_name}corr_plots.svg',
                    format='svg', bbox_inches='tight', dpi=300)


# FUNCTIONS PLOT CONNECT

def plot_connectivity_percentage(df_connections, data_type, save_dir=False):

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

def plot_connect_trace_repatched_cell(save_dir=[]):
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

def plot_corr_single(x, y, df, data_type='all', tr ='Ctrl', save_dir=False):
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
    if 'mini' in data_type:
        ax.set_xlabel(f'Incubation time {titles_dict['hrs_after_OP'][1]})')
        ax.set_xticks(ticks = [10, 20, 30, 40, 50])
        ax.set_ylabel(y)
        ax.set_yticks(ticks = np.linspace(0, np.max(df[y]), 5))
    else:
        ax.set_xlabel(f'{titles_dict[x][0]} ({titles_dict[x][1]})')
        ax.set_xticks(ticks = titles_dict[x][2], labels = [int(x) for x in titles_dict[x][2]])
        ax.set_ylabel(f'{titles_dict[y][0]} \n ({titles_dict[y][1]})')
        ax.set_yticks(ticks = titles_dict[y][2], labels = titles_dict[y][3])

    if save_dir:
        plt.savefig(f'{save_dir}{data_type}{x}{y}corr_plots.svg',
                    format='svg', bbox_inches='tight', dpi=300)

def plot_norm_con_params(df, DV, data_type, save_dir=False):
    '''
    DV - 'amp' or 'lat'
    data_type = 'slice' for all; cell_ID_new for repatch (colors)
    '''
    titles_dict = pl_intr.dict_for_plotting_conn()
    w_cm, h_cm = 6, 4.5
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
            
            for f, amp_col in enumerate(cols):
                
                # Normalize to Amp 1
                normalized_values = df_plot[amp_col] / df_plot[param_1]
                
                # Calculate mean and standard error
                mean_norm = np.mean(normalized_values)
                se_norm = stats.sem(normalized_values)  # Standard error
                
                normalized_means.append(mean_norm)
                normalized_errors.append(se_norm)

                x_temp = np.linspace(1+f-0.3, 1+f+0.3, len(normalized_values))
                # ax.scatter(x_temp, normalized_values, alpha = 0.3)

            
            # Plot line with error bars
            color_key = f'{tr}{day}' if f'{tr}{day}' in colors else tr
            
            ax.errorbar(x_positions, normalized_means, yerr=normalized_errors,
                        color=colors[color_key], marker='o', linewidth=1,
                        markersize=4, capsize=2.2, capthick=1,
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
        plt.savefig(f'{save_dir}_normalized_{DV}_{data_type}.svg',
                    format='svg', bbox_inches='tight', dpi=300)

    plt.show()

def plot_connect_params(df, data_type, DV='amp',save_dir=False):

    h_cm = 4.5
    w_cm = 20
    titles_dict = pl_intr.dict_for_plotting_conn()

    if DV == 'amp':
        cols = ['Amp 1', 'Amp 2', 'Amp 3', 'Amp 4']
        param_1 = 'Amp 1'
    elif DV == 'lat':
        cols = ['Lat1']# , 'Lat2'] # , 'Lat3', 'Lat4']
        param_1 = 'Lat1'

    if 'TTX' in df.treatment.unique():
        treatments = ['Ctrl', 'high K', 'TTX']
    else:
        treatments = ['Ctrl', 'high K']

    colors = figure_colors()[data_type]

    fig, ax = plt.subplots(1, 4, figsize = (w_cm /2.54,h_cm/2.54), sharey = True)
    ax = ax.flatten()
    for a, param in enumerate(cols):
        ticks_ = []
        label_ = []
        num_cels = {}
        
        for i, tr in enumerate(treatments):
            x_plot =[]
            for j, day in enumerate(sorted(df.day.unique())):
                k = j + 2.5*i
                df_plot = df[(df['treatment'] == tr) &\
                              (df['day'] == day)]
                # df_plot = df_plot.dropna(subset=cols)
                df_plot = df_plot[(df_plot[cols] > 0).all(axis=1)]
                
                x = np.linspace(0.65+k, 1.35+k, len(df_plot))
                x_plot.append(x)
                ax[a].scatter(x, df_plot[param], alpha = 0.6,
                            color = colors[tr + day], zorder = 5)
                ax[a].scatter(x, df_plot[param], alpha = 0.7,
                            edgecolor = colors[tr + day],
                            facecolors = 'none', linewidths=1.2)
                ticks_.append(1+k)
                label_.append(day[-1] + tr)
                if 'cell' not in data_type:
                    ax[a].boxplot(df_plot[param], positions = [k + 1], zorder = 10, widths = 0.5)

                    num_cels[tr + ' ' + day] = len(df_plot)
                    #data_boxplot.append(df_plot[amp)
                else:
                    avg = np.mean(df_plot[param])
                    ax[a].text(0.75 + k, (1.05*avg), str(round(avg,2)),
                                    size = 11, zorder = 11)
                    ax[a].scatter(1 + k, avg, color = '#430704',
                                  marker = '_', s = 300, zorder = 10)

            if k in [1,3.5,5] and 'cell' in data_type:
                for c, cell in enumerate(df_plot['connection_ID']):
                    if k == 1 and len(df_plot['connection_ID'])-1 == c:
                        print(str(c+1) + 'Gencho')
                    x1 = [x_plot[0][c], x[c]]
                    y = df[param][df['connection_ID'] == cell]

                    ax[a].plot(x1, y, '-', color = colors[tr+'D2'], alpha = 0.5,
                               linewidth = 2, zorder = 0)

        
        ax[a].set_xticks(ticks = ticks_, labels = label_)
        
        ax[0].set_ylabel(titles_dict[cols[0]][0] + ' (' + titles_dict[cols[0]][1] + ')')
        ax[0].set_yticks(ticks = titles_dict[cols[0]][2], labels = titles_dict[cols[0]][3])
        # ax[0].set_yticks(ticks = np.linspace(0, 5.4, 5), labels = [0, '', 2.7, '', 5.4])
            
    if save_dir:
        plt.savefig(f'{save_dir}con_all_{DV}_{data_type}.svg',
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
                
                num_cels[tr + ' ' + day] = len(df_plot)
                if 'cell' not in data_type:
                    ax[a].boxplot(df_plot[param], positions = [k + 1], zorder = 10, widths = 0.5)

                    label_.append(day + tr)
                    #data_boxplot.append(df_plot[amp)
                else:
                    ax[a].plot(k+1, df_plot[param].mean(), zorder = 10, marker = '_',
                               markersize = 15, c = 'k')
                    ax[a].text(0.75 + k, (1.05 * df_plot[param].mean()), str(round(df_plot[param].mean(),2)),
                                    size = 11, zorder = 11)

            if k in [1,3.5,5] and 'cell' in data_type:
                for c, cell in enumerate(df_plot['connection_ID']):
                    if k == 1 and len(df_plot['connection_ID'])-1 == c:
                        print(str(c+1) + 'Gencho')
                    x1 = [x_plot[0][c], x[c]]
                    y = df[param][df['connection_ID'] == cell]

                    ax[a].plot(x1, y, '-', color = colors[tr+'D2'], alpha = 0.5,
                               linewidth = 2, zorder = 1)
        print(num_cels)
        ax[a].set_xticks(ticks = [1,2,3.5,4.5], labels = ['pre', 'post', 'pre', 'post'])

        ax[a].set_ylabel(titles_dict[param][0] + '\n (' + titles_dict[param][1] + ')')
        if param == 'Lat1':
            ax[a].set_yticks(ticks = np.linspace(0, 4, 5), labels = [0, '', 2, '', 4 ])
        else:
            ax[a].set_yticks(ticks = titles_dict[param][2], labels = titles_dict[param][3])
    if save_dir:
        plt.savefig(f'{save_dir}_con_params_{data_type}_amp1_lat1.svg',
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
    w_cm = 6
    h_cm = 5
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
            ax.scatter(k_+1, ratio, c=colors[tr+day], s = 300, marker = '_', alpha=0.8)
            
            # Add text with ratio and counts
            ax.text(k_+0.85, 1.1*ratio, f'{ratio:.3f}\n({found}/{possible})',
                   size=9, ha='center')
    
    ax.set_xticks(ticks=[0, 1, 2, 3.5, 4.5], labels=['', 'pre', 'post', 'pre', 'post'])
    ax.set_ylabel('Connection ratio\n(found / all possible)')
    ax.set_yticks(np.linspace(0, 0.2, 5), [0, '', 0.1, '', 0.2])
    
    if save_dir:
        plt.savefig(f'{save_dir}{data_type}_connection_ratio2.svg',
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

def plot_connection_window(OP, patcher, labels, fn_indx, pre_chans, post_chans, save_dir = False,
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
        # pre_signal, _, _ = con_param.presynaptic_screen(con_screen_file, pre_ch)
        post_signal, _ = con_param.postsynaptic_screen_ALLOW_ALL(con_screen_file, post_ch, [])
        # post_signal, _ = con_param.postsynaptic_screen(con_screen_file, post_ch, [])
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

        # post_temp = []
        for i in range(np.shape(post_signal)[1]):
            y = post_signal[:,i][preAPs[0][0]-750:preAPs[0][len(preAPs[0])-1]+750]
            #in case want only specific sweeps
                # if np.mean(y) > -63:
                #     continue
                # if len(post_temp) == 0:
                #     post_temp = y
                # else:
                #     upd_post_temp = []
                #     for j, tt in enumerate(y):
                #         upd_post_temp.append((post_temp[j] + tt) /2)
                #     post_temp = upd_post_temp
            ax[1, f].plot(times[0:len(y)], y, lw=0.2, color='grey', alpha=0.4, zorder=0 )
            #ax[1].plot(y, lw=0.2, color='grey', alpha=0.4)
        
        y = post_window
        ax[1, f].plot(times[0:len(y)], post_window, color='k', lw=0.5, zorder=10) #post_window is the averaged 0 shifter post signal
        # ax[1, f].plot(times[0:len(y)], post_temp, color='k', lw=0.5, zorder=10)
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

def plot_connection_window_fig(OP = 'OP240215', patcher = 'Verji', fn_indx = [15, 51], \
                                PRE_chans = [[2, 2, 5, 8, 8], [2, 2, 3, 8, 8]], \
                                POST_chans = [[1, 5, 1, 1, 5], [1, 3, 1, 1, 3]], \
                                save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/figures/trying_out/',\
                                human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    '''
    
    '''
    h_cm = 18
    w_cm = 20

    work_dir, filenames, _, _, _, _ = sort.get_OP_metadata(human_dir, OP, patcher)

    for p in [0,1]:
        pre_chans = PRE_chans[p]
        post_chans = POST_chans[p]
        fig, ax = plt.subplots(3, 5, sharex = 'col', sharey = 'row',
                            figsize = ((w_cm/2.54, h_cm / 2.54)))

        for f, pre_ch in enumerate(pre_chans):
            con_screen_file = f'{work_dir}{filenames[fn_indx[p]]}'
            post_ch = post_chans[f]
            
            # getting needed data
            pre_signal = con_param.presynaptic_screen_ALLOW_ALL(con_screen_file, pre_ch)
            # pre_signal, _, _ = con_param.presynaptic_screen(con_screen_file, pre_ch)
            post_signal, _ = con_param.postsynaptic_screen_ALLOW_ALL(con_screen_file, post_ch, [])
            # post_signal, _ = con_param.postsynaptic_screen(con_screen_file, post_ch, [])
            _, _, pre_window, post_window, preAPs_shifted, preAPs = con_param.get_analysis_window(pre_signal, post_signal)
            _, units, times = hcf.get_abf_info(con_screen_file, pre_ch, np.shape(post_signal)[1], np.shape(post_signal)[0])

            # plotting
            ax[0, f].plot(times[:len(pre_window)], pre_window, color='k')   # plots pre APs 0 shifted

            if post_ch not in [3,5]:
                ax[1, f].plot(times[0:len(post_window)], post_window, color='k', lw=0.5, zorder=10)
                ax[1, f].set_yticks(ticks = np.linspace(-58, -60, 5))
            elif post_ch == 5:
                ax[2, f].plot(times[0:len(post_window)], post_window, color='k', lw=0.5, zorder=10)
                ax[2, f].set_yticks(ticks = np.linspace(-64, -62, 5))
                ax[2, 0].set_yticks(ticks = np.linspace(-64, -62, 5), labels = [-64, '', -61, '', -62])
            elif post_ch == 3:
                ax[2, f].plot(times[0:len(post_window)], post_window, color='k', lw=0.5, zorder=10)
                ax[2, f].set_yticks(ticks = np.linspace(-67.5, -65.5, 5))
                ax[2, 0].set_yticks(ticks = np.linspace(-67.5, -65.5, 5), labels = [-67.5, '', -66.5, '', -65.5])
            
            ax[1, f].set_title(f'Post ch{post_ch}')
            ax[0, f].set_title(f'Pre ch{pre_ch}')
            ax[0, f].set_yticks(ticks = np.linspace(-80, 50, 5))
            
        ax[0, 0].set_yticks(ticks = np.linspace(-80, 50, 5), labels = [-80, '', -15, '', 50])
        ax[0, 0].set_ylabel(str(units)[-2:])

        ax[1, 0].set_yticks(ticks = np.linspace(-58, -60, 5), labels = [-58, '', -59, '', -60])
        
        ax[2, 2].set_xlabel('Time (sec)')
        ax[2, f].set_xticks(ticks = np.linspace(0, 0.24, 5), labels = [0, '', 0.12, '', 0.24])

        if save_dir:
            plt.savefig(f'{save_dir}{OP}-{p}_repatch_con_trace.svg', \
            format = 'svg', bbox_inches = 'tight', dpi = 500)
        plt.show()

def get_repatch_connect_estimation(df_con_repatch_col_QC, df_con_percentage_all):

    '''
    gets a df with connectivity percentage on D1 and D2
    for repatched cells those are obviously ihghly similar
    '''
    # manually constructed dictionary
    # [0] = number of discovered conections,
    # [1] - overall number of repatched cells
    num_repatched = {'OP201027_S1': [1, 2],
                    'OP211209_S5': [5, 5],
                    'OP220426_S3': [1, 2],
                    'OP221027_S1': [2, 4],
                    'OP221027_S3': [1, 4],
                    'OP230208_S2': [1, 3],
                    'OP230420_S2': [1, 5],
                    'OP231109_S2': [1, 2],
                    'OP240201_S2': [1, 4],
                    'OP240215_S2': [5, 5],
                    'OP240411_S2': [1, 4],
                    'OP240417_S1': [1, 3],
                    'OP240503_S1': [3, 4],
                    'OP240503_S2': [7, 4],
                    'OP240503_S3': [5, 5],
                    'OP240507_S2': [1, 2],
                    'OP240926_S2': [2, 4],
                    'OP240926_S3': [1, 4],
                    'OP240926_S4': [3, 6],
                    'OP250207_S1': [1, 3],
                    'OP250207_S2': [1, 4]}

    repatch_connect_percentage = pd.DataFrame(columns = list(df_con_percentage_all.columns.values[1:-1]))
    for op in sorted(df_con_repatch_col_QC.OP.unique()):
        df_op = df_con_repatch_col_QC[df_con_repatch_col_QC.OP == op]
        for slic in df_op.slice.unique():
            op_slic_key = f'{op}_{slic[:2]}'
            df_slic = df_op[df_op.slice == slic]
            if num_repatched[op_slic_key][0] != len(df_slic):
                print('P A N I C')
                break
            if df_slic['comment'].str.contains('appearing', na=False).any():
                # for D1, only remove appearing
                if len(slic) == 2:
                    num_disappearing = len(df_slic[df_slic['comment'].str.contains('disappearing', na=False)])
                    num_appearing = len(df_slic[df_slic['comment'].str.contains('appearing', na=False)])
                    if num_appearing == num_disappearing:
                        num_found = len(df_slic)
                    else:
                        num_found = len(df_slic) - num_appearing + num_disappearing
                else: # D2 only remove disappearing
                    num_disappearing = len(df_slic[df_slic['comment'].str.contains('disappearing', na=False)])
                    num_found = len(df_slic) - num_disappearing
            else:
                num_found = len(df_slic)

            num_possible = num_repatched[op_slic_key][1] * (num_repatched[op_slic_key][1] - 1)
            to_add = pd.DataFrame({'OP': op,
                                'slices': slic,
                                'day': df_slic.day.values[0],
                                'treatment': df_slic.treatment.values[0],
                                'hrs_after_OP': df_slic.hrs_after_OP.values[0],
                                'num_channels': num_repatched[op_slic_key][1],
                                'num_possible_connections': num_possible,
                                'found_connections':  num_found,
                                'con_percentage': num_possible / len(df_slic)}, index=[0])
            repatch_connect_percentage = pd.concat([repatch_connect_percentage[:], to_add]).reset_index(drop=True)
    # fix connectivity percentage
    con_percent_fixed = []
    for con_percent in repatch_connect_percentage.con_percentage:
        if math.isinf(con_percent):
            con_percent_fixed.append(0)
        else:
            con_percent_fixed.append(con_percent)
    repatch_connect_percentage.insert(len(repatch_connect_percentage.columns), 'con_percentage_fixed', con_percent_fixed)
    return repatch_connect_percentage

# funcs to collect intr data for connected cells

def collect_intrinsic_dfs(repatch_con_df_OPs, QC, ops = [],
                          human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    '''
    collects *QC_passed data tables from all OP folders
    adds area and high K concentratioon columns from exp_view
    '''
    exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0]) 

    area_dict, high_k_dict = {}, {}
    for OP in exp_view['OP'].unique():
        area_dict[OP] = exp_view['region'][exp_view['OP'] == OP].tolist()[0]
        high_k_dict[OP] = exp_view['K concentration'][exp_view['OP'] == OP].tolist()[0]

    if QC:
        intr_props_dirs = glob.glob(f'{human_dir}data_*/OP*/data_tables/*QC_passed*.xlsx')
    else:
        intr_props_dirs = glob.glob(f'{human_dir}data_*/OP*/data_tables/*Intrinsic_and_synaptic_properties*.xlsx')

    for i, intr_path in enumerate(intr_props_dirs):
        # if 'QC_passed' in intr_path:
        #     continue
        if '~$' in intr_path:
            continue
        df = pd.read_excel(intr_path)
        if df.empty or df.isna().all().all():
            continue
        if i == 0:
            all_intr_props = df
            continue
        all_intr_props = pd.concat([all_intr_props.loc[:], df]).reset_index(drop=True)
    
    if len(ops) > 0:
        all_intr_props = all_intr_props[all_intr_props.OP.isin(ops)]
    all_intr_props.insert(0, 'area', '')
    for OP in all_intr_props['OP'].unique():
        mask = all_intr_props['OP'] == OP
        all_intr_props.loc[mask, 'area'] = area_dict[OP]
        all_intr_props.loc[mask, 'high K concentration'] = high_k_dict[OP]

    all_intr_props.loc[all_intr_props['treatment'] == 'high k'] = 'high K'

    # remove OPs that are not in repatch connectivity
    repatch_ops = repatch_con_df_OPs.OP.unique().tolist()
    intr_proprs_repatch = all_intr_props[all_intr_props['OP'].isin(repatch_ops)].reset_index(drop = True)
    
    # create match cell IDs
    intr_proprs_repatch = create_cell_ID_connect_match(intr_proprs_repatch)

    print(f'length dataframe: {len(intr_proprs_repatch)}')
    print(f'length cell_ID unique: {len(intr_proprs_repatch.cell_ID_match.unique())}')

    return intr_proprs_repatch

def create_cell_ID_connect_match(df):
    '''
    creates unique cell_IDs to mactch the cells in the connectivity_df
    '''
    cell_IDs = []
    for i, ch_cell in enumerate(df['cell_ch']):
        op = df.OP[i][2:]
        slic = df.slice[i]
        patcher = df.patcher[i]
        cell_IDs.append(f'{patcher}_{op}_{slic}_{int(ch_cell)}')
    df.insert(len(df.columns), 'cell_ID_match', cell_IDs)
    return df

def create_missing_con_IDs(try_1):
    '''
    try_1 - merged df of repatched connections and collected intrinsic dfs QCed
    returns a dataframe of the collected missing cells. This df needs manual fixing!
    '''
    # DEALING WITH MISSING CELLS. MANuaL Check
    connection_counts = try_1['connection_ID'].value_counts()
    # Find connection_IDs with count < 2
    connections_less_than_2 = connection_counts[connection_counts < 2]
    # Get the actual connection_ID values (as a list)
    connection_ids_less_than_2 = connections_less_than_2.index.tolist()
    # get misssing connections from original df
    missing_connections_QC = list(set(df_con_repatch_col_QC.connection_ID) - set(try_1.connection_ID))

    all_missing_cons = list(set(set(connection_ids_less_than_2) | set(missing_connections_QC)))

    ops_add = df_con_repatch_col_QC['OP'][df_con_repatch_col_QC.connection_ID.isin(all_missing_cons)].unique()

    intr_proprs_repatch_no_QC = collect_intrinsic_dfs(df_con_repatch_col_QC, QC = False, ops = ops_add)

    all_missing_cells = list(set(df_con_repatch_col_QC['pre_cell_ID'][df_con_repatch_col_QC.connection_ID.isin(all_missing_cons)].unique().tolist() +\
    df_con_repatch_col_QC['post_cell_ID'][df_con_repatch_col_QC.connection_ID.isin(all_missing_cons)].unique().tolist()))

    missing_cells_df = intr_proprs_repatch_no_QC[intr_proprs_repatch_no_QC.cell_ID_match.isin(all_missing_cells)]

    return missing_cells_df

def get_match_cell_id_connectivity(df):
    '''
    cell_ID = patcher_initial + OP_date + slice + channel
    '''
    if 'cell_IDs_match_pre' in df.columns:
        df = df.drop(columns = ['cell_IDs_match_pre'])

    if 'cell_IDs_match_post' in df.columns:
        df = df.drop(columns = ['cell_IDs_match_post'])

    patcher_dict = {'Verji':'vm', 'Rosie': 'rs'}

    cell_IDs_match_pre, cell_IDs_match_post = [], []
    for i in range(len(df)):

        # print(df.OP[i])
        chan_pre = df.chan_pre[i]
        chan_post = df.chan_post[i]

        slic = df.slice[i]
        
        if isinstance(chan_pre, (float, int, np.int64)):
            chan_pre = str(int(chan_pre))
            chan_post = str(int(chan_post))
        else:
            chan_pre = str(chan_pre)[-1] # ch repatched form D1
            chan_post = str(chan_post)[-1]

        cell_IDs_match_pre.append(patcher_dict[df.patcher[i]] +
                              df.OP[i][2:] +
                              slic +
                              'c' + chan_pre)
        cell_IDs_match_post.append(patcher_dict[df.patcher[i]] +
                              df.OP[i][2:] +
                              slic +
                              'c' + chan_post)

    df.insert(len(df.columns), 'cell_IDs_match_pre', cell_IDs_match_pre)
    df.insert(len(df.columns), 'cell_IDs_match_post', cell_IDs_match_post)

    return df

def delete_repeating_cols(df_repatch_cons_plus_intr):
    '''
    runs through columns of df_repatch_cons_plus_intr and removes the duplicated ones
    that emerged due to merging _pre and _post
    returns cols_delete - list and 
            indx_not_identical - dictionary with cols with non-identical data
    '''
    # get pre_ and post_ columns
    pre_cols = df_repatch_cons_plus_intr.columns[df_repatch_cons_plus_intr.columns.str.contains('pre')]
    post_cols = df_repatch_cons_plus_intr.columns[df_repatch_cons_plus_intr.columns.str.contains('post')]
    other_cols = df_repatch_cons_plus_intr.columns[~((df_repatch_cons_plus_intr.columns.str.contains('pre') | \
                                                    (df_repatch_cons_plus_intr.columns.str.contains('post'))))]
    # identify potential duplicates
    df_repatch_cons_plus_intr = df_repatch_cons_plus_intr.reset_index(drop = True)
    rep_pre_cols, rep_post_cols = [],[]
    for col in other_cols:
        if f'{col}_pre' in pre_cols:
            rep_pre_cols.append(f'{col}_pre')
        if f'{col}_post' in post_cols:
            rep_post_cols.append(f'{col}_post')

    # check if data the same in potentially duplicated columns
    cols_delete = []
    indx_not_identical = {}
    for rep_col in rep_pre_cols:
        data_orig_col = df_repatch_cons_plus_intr[rep_col[:-4]].to_list()
        data_pre_col = df_repatch_cons_plus_intr[rep_col[:-4]].to_list()

        not_ident_indx = []
        for i, dat in enumerate(data_orig_col):
            if dat != data_pre_col[i]:
                not_ident_indx.append(i)
        if len(not_ident_indx) > 0:
            indx_not_identical[rep_col] = not_ident_indx
        else:
            cols_delete.append(rep_col)
    # same for post
    for rep_col in rep_post_cols:
        data_orig_col = df_repatch_cons_plus_intr[rep_col[:-5]].to_list()
        data_post_col = df_repatch_cons_plus_intr[rep_col[:-5]].to_list()

        not_ident_indx = []
        for i, dat in enumerate(data_orig_col):
            if dat != data_post_col[i]:
                not_ident_indx.append(i)
        if len(not_ident_indx) > 0:
            indx_not_identical[rep_col] = not_ident_indx
        else:
            cols_delete.append(rep_col)
    return cols_delete, indx_not_identical

def get_Rs_end(df_repatch_cons_plus_intr, human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    '''
    loads the Rs start end check for
    all of the OPs in df_repatch_cons_plus_intr
    returns df with added columns Rs start and end for pre_  and post_ cells
    '''
    ops = df_repatch_cons_plus_intr.OP.unique().tolist()
    rs_end_dirs = glob.glob(f'{human_dir}data_*/OP*/data_tables/*_QC_measures_rs.xlsx')

    df_merge = pd.DataFrame(['pre_cell_ID', 'post_cell_ID', 'Rs_pre', 'Rs_post'])
    for fn_dir in rs_end_dirs:
        for op in ops:
            if op in fn_dir:
                df_rs = pd.read_excel(fn_dir)
                df_rs = create_cell_ID_connect_match(df_rs)
                if len(df_rs.cell_ID_match.unique()) != len(df_rs):
                    print(f'fix unique IDs for {op}')
                    continue
                if len(df_rs.OP[0]) < 8:
                    print(f'fix OP length for {op}')
                    continue
                pre_cells = df_repatch_cons_plus_intr.pre_cell_ID[df_repatch_cons_plus_intr.OP == op]
                post_cells = df_repatch_cons_plus_intr.post_cell_ID[df_repatch_cons_plus_intr.OP == op]
                if pre_cells.index.equals(post_cells.index):
                    pre_rs_start, post_rs_start, pre_rs_end, post_rs_end = [], [], [], []
                    for i, pre_cell in enumerate(pre_cells.values):
                        post_cell = post_cells.values[i]
                        pre_rs_start.append(df_rs.Rs_start[df_rs.cell_ID_match == pre_cell].values[0])
                        post_rs_start.append(df_rs.Rs_start[df_rs.cell_ID_match == post_cell].values[0])
                        pre_rs_end.append(df_rs.Rs_end[df_rs.cell_ID_match == pre_cell].values[0])
                        post_rs_end.append(df_rs.Rs_end[df_rs.cell_ID_match == post_cell].values[0])
           
                    df_add = pd.DataFrame( {'pre_cell_ID': pre_cells.values,
                                            'post_cell_ID': post_cells.values,
                                            'Rs_start_pre': pre_rs_start,
                                            'Rs_start_post':post_rs_start,
                                            'Rs_end_pre': pre_rs_end,
                                            'Rs_end_post':post_rs_end})
                    df_merge = pd.concat([df_merge, df_add]).reset_index(drop=True)
                else:
                    print(f'weird mistake with {op}')
    
    out_df_complete = pd.merge(df_repatch_cons_plus_intr, df_merge, on = ['pre_cell_ID', 'post_cell_ID'],
                               suffixes = (None, '_added')).reset_index(drop=True)
    
    #out_df_complete = out_df_complete.loc[:, ~out_df_complete.columns.str.contains('Unnamed')]
    out_df_complete = pd.concat([out_df_complete, \
                                 df_repatch_cons_plus_intr[df_repatch_cons_plus_intr.OP == 'OP211209']]).reset_index(drop = True)
    return out_df_complete

def rename_pre_cols(df_repatch_cons_plus_intr):
    '''
    because of previous merge pre cols,don't have'_pre' suffinx
    renaming them
    '''
    post_cols = df_repatch_cons_plus_intr.columns[df_repatch_cons_plus_intr.columns.str.contains('post')]

    for post_col in post_cols:
        root_col = post_col[:-5]
        if root_col in df_repatch_cons_plus_intr.columns:
            print(f'renaming {root_col} to {root_col}_pre')
            df_repatch_cons_plus_intr = df_repatch_cons_plus_intr.rename(columns={root_col: f'{root_col}_pre'})
    
    pre_cols = df_repatch_cons_plus_intr.columns[df_repatch_cons_plus_intr.columns.str.contains('pre')]
    if len(pre_cols) == len(post_cols):
        return df_repatch_cons_plus_intr

def just_code_to_get_intr_with_repeatch_cons(data_dir):
    '''
    THIS CAN'T be USED AS A FUNCTION
    '''
    df_con_repatch_col_QC = pd.read_excel(f'{data_dir}/repatched_connections_looser_QC.xlsx')
    intr_proprs_repatch = collect_intrinsic_dfs(df_con_repatch_col_QC, QC = True)

    try_1 = pd.merge(df_con_repatch_col_QC, intr_proprs_repatch, left_on = 'pre_cell_ID', \
                            right_on = 'cell_ID_match', suffixes = (None, '_pre'), \
                                indicator = True).reset_index(drop=True)

    print(f'merging on pre_cells. Length {len(try_1)}')

    # merge on post synaptic
    try_1 = pd.merge(try_1, intr_proprs_repatch, left_on = 'post_cell_ID', \
                            right_on = 'cell_ID_match',suffixes = (None, '_post')).reset_index(drop=True)
    print(f'merging on post_cells. Length {len(try_1)}')
    try_1 = get_match_cell_id_connectivity(try_1)


    data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/connectivity/'
    # loading the manually checked df. alternative use create_missing_con_IDs(try_1)
    missing_cells_df = pd.read_excel(f'{data_dir}missing_cells_manual_check.xlsx')
    missing_cells_df['cell_ID_match'].value_counts()

    # MERGIING WITH MISSING CELLS
    try_2 = pd.merge(df_con_repatch_col_QC, missing_cells_df, left_on = 'pre_cell_ID', \
                            right_on = 'cell_ID_match', suffixes = (None, '_pre'), \
                                indicator = True).reset_index(drop=True)

    print(f'merging on pre_cells. Length {len(try_2)}')

    # merge on post synaptic
    try_2 = pd.merge(try_2, missing_cells_df, left_on = 'post_cell_ID', \
                            right_on = 'cell_ID_match',suffixes = (None, '_post')).reset_index(drop=True)
    print(f'merging on post_cells. Length {len(try_2)}')
    try_2 = get_match_cell_id_connectivity(try_2)


    # loading repatch dataframe creted by matching adult intr QCed
    df_connections_QC = pd.read_excel(data_dir + 'all_QCed_connections.xlsx')
    df_repatch_QCed = df_connections_QC[df_connections_QC.repatch == 'yes']

    print(f'most strict QC connectivity df length: {len(df_repatch_QCed)}\n\
        number of unique connection in strict QCed: {len(df_repatch_QCed.connection_ID.unique())} \n\
        colleced QC connectivity df length: {len(try_1)}\n\
        number of unique connection in collected QCed: {len(try_1.connection_ID.unique())} \n\
        colleced non-QC connectivity df length: {len(try_2)}\n\
        number of unique connection in collected QCed: {len(try_2.connection_ID.unique())} \n')


    # used the following ones already and a cleaner dataframe was saved
    # df_repatch_cons_plus_intr = pd.concat([try_1, try_2]).sort_values('connection_ID')
    # df_repatch_cons_plus_intr = df_repatch_cons_plus_intr.loc[:, ~df_repatch_cons_plus_intr.columns.str.contains('Unnamed')]
    # df_repatch_cons_plus_intr.connection_ID.value_counts() > 2
    # cols_delete, indx_not_identical = delete_repeating_cols(df_repatch_cons_plus_intr)
    # df_repatch_cons_plus_intr = df_repatch_cons_plus_intr.loc[:, ~df_repatch_cons_plus_intr.columns.isin(cols_delete)]
    # df_repatch_cons_plus_intr = get_Rs_end(df_repatch_cons_plus_intr)
    ### df_repatch_cons_plus_intr.to_excel(f'{data_dir}all_repatched_and_intr_manually_checked.xlsx')

def QC_plot_repatched_connections(df, dest_dir, w_cm = 20, h_cm = 15):
    '''
    plots intrinsic parameters for pre and post synaptic cell
    '''

    params = ['Rs_start', 'Rs_end', 'Vm', 'Rheobse_ramp', 'TH']

    df = df.sort_values(['treatment', 'day'])
    titles_all = pl_intr.dict_for_plotting_conn()
    titles_dict = {key: titles_all[key] for key in params if key in titles_all}
    colors = figure_colors()['cell_ID_new']

    fig2, axs = plt.subplots(2,3, figsize = (w_cm/2.54, h_cm / 2.54))
    axs = axs.flatten()
    for a, (param, item) in enumerate(titles_dict.items()):
        nan_connection_ids = []
        label_, ticks_ = [], []
        for t, tr in enumerate(df.treatment.unique()):
            x_plot =[]
            for j, day in enumerate(df['day'].unique()):
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                df_plot.reset_index(drop=True, inplace=True)
                
                for c, cell in enumerate(['_pre', '_post']):
                    k = 2.5*c + j
                    col = f'{param}{cell}'

                    # Check for NaNs and collect connection_IDs to exclude
                    nan_mask = df_plot[col].isna()
                    if nan_mask.any():
                        # Get connection_IDs of rows with NaN values
                        nan_connection_ids.append(df_plot.loc[nan_mask, 'connection_ID'].tolist())

                    df_plot = df_plot[~nan_mask]
                    mean = np.nanmean(df_plot[col])

                    x = np.linspace(0.7+k, 1.3+k, len(df_plot))
                    x_plot.append(x)
                    axs[a].scatter(x, df_plot[col], c = colors[tr + day], zorder = 5, alpha = 0.3)
                    print(a)
                    axs[a].scatter(x, df_plot[col], edgecolor = colors[tr + day], alpha = 0.75,facecolors = 'none')
                    axs[a].scatter(1 + k, mean, c = colors[tr + day], edgecolor = 'k', marker = '_', s = 300, zorder = 10)

                    axs[a].text(0.75 + k, (1.05*mean), str(round(mean,1)),
                                    size = 11, c = colors[tr + 'D2'], zorder = 11)

                    # x_plot.append(x)
                    # if k in [1,3.5,5]:
                    #     # df = df[~df['connection_ID'].isin(nan_connection_ids)]
                    #     for s, con in enumerate(df_plot[f'connection_ID']):
                    #         x1 = [x_plot[0][s], x[s]]
                    #         y = df[col][df['connection_ID'] == con]
                    #         axs[a].plot(x1, y, '-', color = colors[tr + day],
                    #                 alpha = 0.4, zorder = 1)
                    if t == 0:
                        ticks_.append(k)
                        label_.append(f'{cell[1:]}synaptic cell \n {day}')

                # ax.text(k + 0.65, int((1 + 0.1)*np.max(df[param])), 'n = ' + str(len(df_plot)), size = 12, c = colors[tr + day])

        axs[a].set_title(param)
        axs[a].set_xticks(ticks = ticks_, labels = label_, rotation = 30)
        axs[a].set_ylabel(item[0] + ' (' + item[1] + ')')
        axs[a].set_yticks(ticks = item[2], labels = item[3])

    plt.savefig(dest_dir + '_repatched_connections_QC.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 300)
    plt.close(fig2)

def QC_plot_single_cell_param_repatched_cons(df, dest_dir, params = ['Rin']):
    '''
    plots the rin for the presynaotic cells
    '''
    df = df.sort_values(['treatment', 'day'])
    titles_all = pl_intr.dict_for_plotting()
    if 'Rin' in params[0]:
        titles_dict = {key: titles_all[key] for key in params if key in titles_all}
    else:
        param = 'Rs'
        titles_dict = {param: titles_all[param]}
        print(titles_dict)
    colors = figure_colors()['cell_ID_new']

    fig2, axs = plt.subplots(1,2, figsize = (14/2.54, 8/2.54))
    axs = axs.flatten()
    for g, (param, item) in enumerate(titles_dict.items()):
        nan_connection_ids = []
        label_, ticks_ = [], []
        for t, tr in enumerate(df.treatment.unique()):
            x_plot =[]
            for j, day in enumerate(df['day'].unique()):
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                df_plot.reset_index(drop=True, inplace=True)
                
                for a, cell in enumerate(['_end_pre', '_end_post']):
                    k = j
                    col = f'{param}{cell}'

                    # Check for NaNs and collect connection_IDs to exclude
                    nan_mask = df_plot[col].isna()
                    if nan_mask.any():
                        # Get connection_IDs of rows with NaN values
                        nan_connection_ids.append(df_plot.loc[nan_mask, 'connection_ID'].tolist())

                    df_plot = df_plot[~nan_mask]
                    mean = np.nanmean(df_plot[col])

                    x = np.linspace(0.7+k, 1.3+k, len(df_plot))
                    x_plot.append(x)
                    axs[a].scatter(x, df_plot[col], c = colors[tr + day], zorder = 5, alpha = 0.3)
                    print(a)
                    axs[a].scatter(x, df_plot[col], edgecolor = colors[tr + day], alpha = 0.75,facecolors = 'none')
                    axs[a].scatter(1 + k, mean, c = colors[tr + day], edgecolor = 'k', marker = '_', s = 300, zorder = 10)

                    axs[a].text(0.75 + k, (1.05*mean), str(round(mean,1)),
                                    size = 11, c = colors[tr + 'D2'], zorder = 11)

                    if t == 0:
                        ticks_.append(k)
                        label_.append(f'presynaptic cell \n {day}')

                # ax.text(k + 0.65, int((1 + 0.1)*np.max(df[param])), 'n = ' + str(len(df_plot)), size = 12, c = colors[tr + day])

        axs[a].set_title(param)
        axs[a].set_xticks(ticks = ticks_, labels = label_, rotation = 30)
        axs[a].set_ylabel(item[0] + ' (' + item[1] + ')')
        axs[a].set_yticks(ticks = item[2], labels = item[3])

    plt.savefig(f'{dest_dir}{params[0]}_repatched_connections_QC.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 300)
    plt.close(fig2)


def QC_plot_Rs_vs_EPSP_amp(df, dest_dir):
    '''
    plots the rin for the presynaotic cells
    '''
    df = df.sort_values(['treatment', 'day'])
    titles_con = pl_intr.dict_for_plotting_conn()

    titles_dict = {'Rs': titles_con['Rs'], 'Amp 1': titles_con['Amp 1']}
    colors = figure_colors()['cell_ID_new']

    fig2, axs = plt.subplots(1,2, figsize = (8/2.54, 4/2.54))
    axs = axs.flatten()
  
    for t, tr in enumerate(df.treatment.unique()):
    
        df_plot_D1 = df[(df['treatment'] == tr) & (df['day'] == 'D1')]
        df_plot_D2 = df[(df['treatment'] == tr) & (df['day'] == 'D2')]

        df_plot_D1 = df_plot_D1.sort_values('connection_ID')
        df_plot_D2 = df_plot_D2.sort_values('connection_ID')
        
        if df_plot_D1.connection_ID.tolist() != df_plot_D2.connection_ID.tolist():
            print('some data missing connections D1 != D2')
            return
        
        for a, cell in enumerate(['_end_pre', '_end_post']):
            col = f'Rs{cell}'

            x = df_plot_D1['Amp 1'].values - df_plot_D2['Amp 1'].values
            y = df_plot_D1[col].values - df_plot_D2[col].values
            axs[a].scatter(x, y, c = colors[tr + 'D2'], zorder = 5 , alpha = 0.3)
            axs[a].scatter(x, y, edgecolor = colors[tr + 'D2'], alpha = 0.75, facecolors = 'none')

            axs[a].set_xticks(ticks = np.linspace(-1, 1, 5))
            axs[a].set_xlabel(rf'$\Delta$ 1st EPSP amplitude \n (D1 - D2)')

            axs[a].set_yticks(ticks = np.linspace(-35, 35, 5), labels = [-35, '', 0, '', 35])
            axs[a].set_ylabel(rf'$\Delta$ Series resistance \n (D1 - D2)')

        axs[0].set_title('Presynaptic cell')
        axs[1].set_title('Postsynaptic cell')

    plt.savefig(f'{dest_dir}corr_delta_Rs_EPSP_amp.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 300)
    plt.close(fig2)

def get_tidy_data_for_repatch_con_plus_intd():
    '''
    MAYBE TPOTALY USELESS AND WRONG
    '''
    # define columns
    pre_cols = [col for col in df_repatch_cons_plus_intr.columns if col.endswith('_pre')]
    post_cols = [col for col in df_repatch_cons_plus_intr.columns if col.endswith('_post')]
    id_cols = [col for col in df_repatch_cons_plus_intr.columns if not (col.endswith('_pre') or col.endswith('_post'))]

    # Remove repeating cell_IDs
    df_pre = df_repatch_cons_plus_intr.drop_duplicates(subset='cell_IDs_match_pre').reset_index(drop=True)
    pre_cells = df_pre.cell_IDs_match_pre.unique().tolist()
    # Melt pre columns
    df_pre = df_pre[id_cols + pre_cols].melt(
        id_vars=id_cols,
        value_vars=pre_cols,
        var_name='parameter',
        value_name='value'
    )
    df_pre['cell_type'] = 'pre'
    df_pre['parameter'] = df_pre['parameter'].str.replace('_pre', '')

    # Remove repeating cell_IDs. Also ones from df_pre
    df_post = df_repatch_cons_plus_intr.drop_duplicates(subset='cell_IDs_match_post').reset_index(drop=True)
    indx_in_pre = [i for i, cell in enumerate(df_post.cell_IDs_match_post) if cell in pre_cells]
    df_post = df_post.drop(index=indx_in_pre).reset_index(drop=True)

    # Melt post columns
    df_post = df_post[id_cols + post_cols].melt(
        id_vars=id_cols,
        value_vars=post_cols,
        var_name='parameter',
        value_name='value'
    )


    df_post['cell_type'] = 'post'
    df_post['parameter'] = df_post['parameter'].str.replace('_post', '')

    # Combine
    df_tidy = pd.concat([df_pre, df_post], ignore_index=True)
    df_tidy = df_tidy.loc[:, df_tidy.columns[~df_tidy.columns.str.contains('Unnamed')]]

    return df_tidy

def get_long_data_for_repatch_con_plus_intr(df_repatch_cons_plus_intr):
    '''
    gets all unique cells that are in the repatched connections
    gets kind of a tidy data for that
    '''
    # Remove repeating cell_IDs
    df_pre = df_repatch_cons_plus_intr.drop_duplicates(subset='cell_IDs_match_pre').reset_index(drop=True)
    df_pre = df_pre.loc[:, df_pre.columns[df_pre.columns.str.contains('pre')]]
    pre_cells = df_pre.cell_IDs_match_pre.unique().tolist()

    # Remove repeating cell_IDs. Also ones from df_pre
    df_post = df_repatch_cons_plus_intr.drop_duplicates(subset='cell_IDs_match_post').reset_index(drop=True)
    post_cells = df_post.cell_IDs_match_pre.unique().tolist()

    indx_in_pre = [i for i, cell in enumerate(df_post.cell_IDs_match_post) if cell in pre_cells]
    df_post = df_post.drop(index=indx_in_pre).reset_index(drop=True)

    # Combine
    df_tidy = pd.concat([df_pre, df_post], ignore_index=True)
    df_tidy = df_tidy.loc[:, df_tidy.columns[~df_tidy.columns.str.contains('Unnamed')]]

    # # get cell_ID
    # cell_types = []
    # for i, cell in enumerate(df_tidy.cell_ID_match):
    #     if cell in pre_cells and cell in post_cells:
    #         cell_types.append('both')
    #     elif cell in pre_cells:
    #         cell_types.append('pre')
    #     elif cell in post_cells:
    #         cell_types.append('post')
    # df_tidy.insert(len(df_tidy.columns), 'cell_type', cell_types)

    return df_tidy

def create_hist_sup_figures(df, param, x_axis_label, save_dir):
    '''
    creates a simple histogram with mean +/- SD in title
    '''
    # Create histogram

    plt.figure(figsize=(6/2.54, 6/2.54))
    plt.hist(df[param], bins=15, color='grey', alpha=0.7, edgecolor='black')
    plt.xlabel(x_axis_label)
    plt.ylabel('Frequency')
    # Calculate mean and standard deviation
    mean_val = np.mean(df[param])
    std_val = np.std(df[param])

    # Add text in top right corner
    plt.text(1, 1, f'Mean: {mean_val:.2f}\nSD: {std_val:.2f}', 
            transform=plt.gca().transAxes,
            verticalalignment='top', horizontalalignment='right')

    plt.savefig(f'{save_dir}{x_axis_label}_hist.svg', 
                format='svg', bbox_inches='tight', dpi=300)

def create_donut_plots(size_C = 15, size_K = 28, save_dir = False):
    '''
    creates donut plots for repatched connections/ appearing and disappearing
    '''
    # Create figure with 2 subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12/2.54, 6/2.54))
    
    # FIRST DONUT PLOT (left)
    # Create data for 15 equal parts
    sizes1 = [1] * size_C
    colors1 = ['#ff7f00'] * size_C
    alphas1 = [0.75] * size_C

    # Set 2 adjacent parts to have alpha 1.0 (parts 7 and 8, which are next to each other)
    highlight_indices1 = [1, 2]  # 0-based indexing for parts 7 and 8
    for idx in highlight_indices1:
        alphas1[idx] = 0.5

    # Create the first pie chart without labels
    wedges1, texts1 = ax1.pie(sizes1, colors=colors1, startangle=90, counterclock=False)

    # Apply alpha values to each wedge
    for i, wedge in enumerate(wedges1):
        wedge.set_alpha(alphas1[i])

    # Add text to non-highlighted parts
    for i, wedge in enumerate(wedges1):
        if i not in highlight_indices1:
            # Get the center angle of the wedge
            angle = (wedge.theta1 + wedge.theta2) / 2
            # Convert to radians
            angle_rad = np.radians(angle)
            # Position text at 75% of the radius (between center and edge)
            radius = 0.75
            x = radius * np.cos(angle_rad)
            y = radius * np.sin(angle_rad)
            
            # Add text
            # ax1.text(x, y, 'repatched\nconnections\ncontrol', 
            #         ha='center', va='center', fontsize=8, 
            #         rotation=angle-90 if angle > 90 and angle < 270 else angle+90,
            #         color='white', weight='bold')

    # Create the donut hole for first plot
    centre_circle1 = plt.Circle((0, 0), 0.5, fc='white')
    ax1.add_artist(centre_circle1)
    ax1.axis('equal')
    
    # SECOND DONUT PLOT (right)
    # Create data for 15 equal parts
    sizes2 = [1] * size_K  # 15 equal parts
    colors2 = ["#4A74DD"] * size_K
    alphas2 = [0.75] * size_K
    
    # Set 6 adjacent parts: 4 fainter (alpha 0.3) and 2 black
    highlight_indices2 = list(range(6))  # parts 0-5 (6 adjacent parts)
    faint_indices = highlight_indices2[:4]  # first 4 parts fainter
    black_indices = highlight_indices2[4:]  # last 2 parts black
    
    # Set alpha values and colors
    for idx in faint_indices:
        alphas2[idx] = 0.5  # fainter
    
    for idx in black_indices:
        # colors2[idx] = 'black'  # black color
        alphas2[idx] = 1.0

    # Create the second pie chart without labels
    wedges2, texts2 = ax2.pie(sizes2, colors=colors2, startangle=90, counterclock=False)

    # Apply alpha values to each wedge
    for i, wedge in enumerate(wedges2):
        wedge.set_alpha(alphas2[i])

    # Add text to non-highlighted parts (parts 6-14)
    for i, wedge in enumerate(wedges2):
        if i not in highlight_indices2:
            # Get the center angle of the wedge
            angle = (wedge.theta1 + wedge.theta2) / 2
            # Convert to radians
            angle_rad = np.radians(angle)
            # Position text at 75% of the radius (between center and edge)
            radius = 0.75
            x = radius * np.cos(angle_rad)
            y = radius * np.sin(angle_rad)
            
            # Add text
            # ax2.text(x, y, 'appearing\ndisappearing\nconnections',
            #         ha='center', va='center', fontsize=8,
            #         rotation=angle-90 if angle > 90 and angle < 270 else angle+90,
            #         color='white', weight='bold')

    # Create the donut hole for second plot
    centre_circle2 = plt.Circle((0, 0), 0.5, fc='white')
    ax2.add_artist(centre_circle2)
    ax2.axis('equal')

    # Save if directory provided
    if save_dir:
        plt.savefig(f'{save_dir}_donut_num_connections.svg',
                    format='svg', bbox_inches='tight', dpi=300)

    plt.show()

    plt.show()

def scale_bar(abf=None, hideTicks=True, hideFrame=True, fontSize=8, scaleXsize=None, \
              scaleYsize=None, scaleXunits="", scaleYunits="", lineWidth=2):
    """ 
    Add an L-shaped scalebar to the current figure. 
    This removes current axis labels, ticks, and the figure frame. 
    """ 

    # if an ABF objet is given, use its sweep units
    if abf:
        scaleXunits = abf.sweepUnitsX
        scaleYunits = abf.sweepUnitsY 

    # calculate the current data area
    x1, x2, y1, y2 = plt.axis()  # bounds
    xc, yc = (x1+x2)/2, (y1+y2)/2  # center point
    xs, ys = abs(x2-x1), abs(y2-y1)  # span

    # determine how big we want the scalebar to be
    if not scaleXsize:
        scaleXsize = abs(plt.xticks()[0][1]-plt.xticks()[0][0])/2
    if not scaleYsize: 
        scaleYsize = abs(plt.yticks()[0][1]-plt.yticks()[0][0])/2

    # create the scale bar labels
    lblX = str(scaleXsize)
    lblY = str(scaleYsize)

    # prevent units unecessarially ending in ".0"
    if lblX.endswith(".0"): 
        lblX = lblX[:-2] 
    if lblY.endswith(".0"): 
        lblY = lblY[:-2] 

    if scaleXunits == "sec" and "." in lblX: 
        lblX = str(int(float(lblX)*1000)) 
        scaleXunits = "ms"

    # add units to the labels
    lblX = lblX+" "+scaleXunits
    lblY = lblY+" "+scaleYunits
    lblX = lblX.strip()
    lblY = lblY.strip()

    # determine the dimensions of the scalebar
    scaleBarPadX = 0.10
    scaleBarPadY = 0.05
    scaleBarX = x2-scaleBarPadX*xs
    scaleBarX2 = scaleBarX-scaleXsize
    scaleBarY = y1+scaleBarPadY*ys
    scaleBarY2 = scaleBarY+scaleYsize

    # determine the center of the scalebar (where text will go)
    scaleBarXc = (scaleBarX+scaleBarX2)/2
    scaleBarYc = (scaleBarY+scaleBarY2)/2

    # create a scalebar point array suitable for plotting as a line
    scaleBarXs = [scaleBarX2, scaleBarX, scaleBarX]
    scaleBarYs = [scaleBarY, scaleBarY, scaleBarY2]

    # the text shouldn't touch the scalebar, so calculate how much to pad it
    lblPadMult = .005
    lblPadMult += .002*lineWidth
    lblPadX = xs*lblPadMult
    lblPadY = ys*lblPadMult

    # hide the old tick marks
    if hideTicks:
        plt.gca().get_yaxis().set_visible(False)
        plt.gca().get_xaxis().set_visible(False)

    # hide the square around the image
    if hideFrame:
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)

    # now do the plotting
    plt.plot(scaleBarXs, scaleBarYs, 'k-', lw=lineWidth)
    plt.text(scaleBarXc, scaleBarY-lblPadY, lblX,
            ha='center', va='top', fontsize=fontSize)
    plt.text(scaleBarX+lblPadX, scaleBarYc, lblY,
            ha='left', va='center', fontsize=fontSize)

def plot_appear_disappear_OLD(disappear_only, appear_only, save_dir = False):
    '''
    
    '''
    colors = figure_colors()['cell_ID_new']
    titles_dict = pl_intr.dict_for_plotting_conn()
    size_ = 33

    fig, ax = plt.subplots(1,1, figsize = (6/2.54, 6/ 2.54))
    df = disappear_only
    x_d1 , x_d2, y_d1, y_d2 = [], [], [], []
    for a, day in enumerate(['D1', 'D2']):
        for tr in ['Ctrl', 'high K']:
            df_tr = df['Amp 1'][(df.day == day) & (df.treatment == tr)]
            x = np.linspace(a - 0.2, a + 0.2, len(df_tr))
            if day == 'D1':
                ax.scatter(x, df_tr, c = colors[tr + day], zorder = 5, alpha = 0.88,  s = size_)
                x_d1.extend(x.tolist())
                y_d1.extend(df_tr.to_list())
            else:
                ax.scatter(x, np.repeat(0, len(df_tr)), c = colors[tr + day], zorder = 5,
                           alpha = 0.3,  s =size_)
                ax.scatter(x, np.repeat(0, len(df_tr)), edgecolor = colors[tr + day],
                   facecolors = 'none', alpha = 0.5,  s = size_)
                x_d2.extend(x.tolist())
                y_d2.extend(np.repeat(0, len(df_tr)).tolist())
    for j, x in enumerate(x_d1):
        ax.plot([x, x_d2[j]], [y_d1[j], y_d2[j]], colors[tr+'D2'], alpha = 0.5,
                linewidth = 2, zorder = 0)
            
    df = appear_only
    x_d1 , x_d2, y_d1, y_d2 = [], [], [], []
    for a, day in enumerate(['D1', 'D2']):
        df_tr = df['Amp 1'][(df.day == day)]
        x = np.linspace(a - 0.2, a + 0.2, len(df_tr))
        if day == 'D2':
            ax.scatter(x, df_tr, c = colors[tr + day], zorder = 5,alpha = 0.88, s = size_)
            x_d2.extend(x.tolist())
            y_d2.extend(df_tr.to_list())
        else:
            ax.scatter(x, np.repeat(0, len(df_tr)), c = colors[tr + day], zorder = 5,
                       alpha = 0.3, s = size_)
            ax.scatter(x, np.repeat(0, len(df_tr)), edgecolor = colors[tr + day],
                facecolors = 'none', alpha = 0.5,  s = size_)
            x_d1.extend(x.tolist())
            y_d1.extend(np.repeat(0, len(df_tr)).tolist())
    for j, x in enumerate(x_d1):
        ax.plot([x, x_d2[j]], [y_d1[j], y_d2[j]], colors[tr+'D2'], alpha = 0.5, linewidth = 2, zorder = 0)
            
    ax.set_xticks(ticks = [0, 1], labels = ['pre', 'post'])

    ax.set_ylabel('EPSP amplitude (mV)')
    ax.set_yticks(ticks = np.linspace(0, 0.6, 5), labels = [0, '', 0.3, '', 0.6])
    
    if save_dir:
        plt.savefig(f'{save_dir}/disappear_appear_amplitudes.svg',
                    format='svg', bbox_inches='tight', dpi=300)
    plt.show(fig)

def plot_appear_disappear(disappear_only, appear_only, save_dir = False):
    '''
    
    '''
    colors = figure_colors()['cell_ID_new']
    titles_dict = pl_intr.dict_for_plotting_conn()
    size_ = 33

    fig, ax = plt.subplots(1,1, figsize = (4.5/2.54, 6.2/2.54))
    df = disappear_only
    a = 0 
    for tr in ['Ctrl', 'high K']:
        df_tr = df['Amp 1'][(df.day == 'D1') & (df.treatment == tr)]
        x = np.linspace(a - 0.1, a + 0.1, len(df_tr))
        ax.scatter(x, df_tr, c = colors[tr + 'D1'], zorder = 5, alpha = 0.88,  s = size_)
        ax.scatter(x, df_tr, edgecolor = colors[tr + 'D2'], facecolors = 'none',
                    zorder = 5, alpha = 0.88,  s = size_)

    df = appear_only
    a = 1
    df_tr = df['Amp 1'][(df.day == 'D2')]
    x = np.linspace(a - 0.1, a + 0.1, len(df_tr))
    ax.scatter(x, df_tr, c = colors[tr + 'D2'], zorder = 5,alpha = 0.88, s = size_)
       
    ax.set_xticks(ticks = [-0.4 ,0, 1, 1.4], labels = ['', 'disappearing', 'appearing', ''],
                  rotation = 20)

    ax.set_ylabel('EPSP amplitude (mV)')
    ax.set_yticks(ticks = np.linspace(0, 0.6, 5), labels = [0, '', 0.3, '', 0.6])
    
    if save_dir:
        plt.savefig(f'{save_dir}/disappear_appear_amplitudes.svg',
                    format='svg', bbox_inches='tight', dpi=300)
    plt.show(fig)


def plot_distribution_1st_amp_lat_only(df, data_type, dest_dir = False):
    '''
    Plots side-by-side histogram distributions for Amp 1 and Lat1.
    If 'cell' in data_type, overlays disappearing connections (dotted lines) 
    and appearing connections (solid lines).
    '''
    if 'cell' in data_type:
        df_repatch_all = df
        df = df_repatch_all[~df_repatch_all['comment'].str.contains('appearing|disappearing', case=False, na=False)].reset_index(drop = True)
        df_special = df_repatch_all[df_repatch_all['comment'].str.contains('appearing|disappearing', case=False, na=False)].reset_index(drop = True)
        df_special_dis = df_special[(df_special.day == 'D1') & (df_special['Amp 1'] > 0.08)].reset_index(drop = True)
        df_special_app = df_special[(df_special.day == 'D2') & (df_special['Amp 1'] > 0.15)].reset_index(drop = True)

    h_cm = 6
    w_cm = 12
    titles_dict = pl_intr.dict_for_plotting_conn()
    
    cols = ['Amp 1', 'Lat1']
    colors = figure_colors()[data_type]
    
    if 'cell' in data_type:
        # Get connection_IDs from rows that would be dropped due to NaN in specific columns
        nan_mask = df[cols].isnull().any(axis=1)
        con_ids_drop = df[nan_mask]['connection_ID'].unique()
        # Remove all rows that have those connection_IDs
        df = df[~df['connection_ID'].isin(con_ids_drop)]
    
    fig, ax = plt.subplots(1, 2, figsize=(w_cm/2.54, h_cm/2.54))
    ax = ax.flatten()
    
    for a, param in enumerate(cols):
        # Define bins
        bins = np.linspace(df[param].min(), df[param].max(), 15)
        bin_width = bins[1] - bins[0]
        bar_width = bin_width / 3.5
        
        for i, tr in enumerate(sorted(df.treatment.unique())):
            for j, day in enumerate(sorted(df.day.unique())):
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                df_plot = df_plot.dropna(subset=[param])
                
                # Calculate histogram
                counts, _ = np.histogram(df_plot[param], bins=bins)
                
                # Calculate x positions for bars
                bin_centers = (bins[:-1] + bins[1:]) / 2
                k = j + 2.5*i
                x_pos = bin_centers + (k - 1.5) * bar_width
                
                # Plot bars
                ax[a].bar(x_pos, counts, width=bar_width, alpha=0.8, color=colors[tr + day], label=f'{tr} {day}')
                
        # Plot disappearing connections (dotted lines) if 'cell' in data_type
        if 'cell' in data_type and param == 'Amp 1':
            for idx, row in df_special_dis.iterrows():
                treatment = row['treatment']
                amp_value = row[param]
                if pd.notna(amp_value):
                    ax[a].axvline(x=amp_value, linestyle=':', color= 'k',
                                 linewidth=0.9, alpha=0.8, zorder = 0)
            
            # Plot appearing connections (solid lines)
            for idx, row in df_special_app.iterrows():
                treatment = row['treatment']
                amp_value = row[param]
                if pd.notna(amp_value):
                    ax[a].axvline(x=amp_value, linestyle='-', color='k', 
                                 linewidth=0.9, alpha=0.8, zorder = 0)
        
        ax[a].set_xlabel(titles_dict[param][0] + '\n(' + titles_dict[param][1] + ')')
        ax[a].set_ylabel('Count')
        
        if param == 'Lat1':
            ax[a].set_xticks(np.linspace(0, 4, 5))
        else:
            ax[a].set_xticks(titles_dict[param][2])
    
    plt.tight_layout()
    if dest_dir:
        plt.savefig(f'{dest_dir}_histogram_{data_type}_amp1_lat1.svg', dpi=300)
    plt.show()

def plot_op_slice_cell_histogram(df, data_type, save_dir=None):

    """
    Creates a histogram with bars for each OP showing:
    - Bar height = total number of cells per OP
    - Different colors for each slice within each OP
    - Numbers on bars = cell counts per slice
    """
    
    h = 6
    w = 3.5
    if data_type == 'repatch':
        df = df[df['day'] == 'D1']
        col_name = 'cell_ID_new'
    else: 
        col_name = 'cell_ID'
    for day in sorted(df.day.unique()):

        fig, ax = plt.subplots(figsize=(w/2.54, h/2.54))

        df_g = df[df['day'] == day]
        grouped = df_g.groupby(['OP', 'slice'])[col_name].nunique().reset_index()
        grouped.columns = ['OP', 'slice', 'cell_count']

        # Get unique OPs and create color map for slices
        ops = sorted(grouped['OP'].unique())
        all_slices = sorted(grouped['slice'].unique())
        slice_colors = plt.cm.Accent(np.linspace(0, 1, len(all_slices)))
        
        # Create color mapping for slices
        slice_color_map = {slice_name: slice_colors[i] for i, slice_name in enumerate(all_slices)}
    
        x_pos = np.arange(len(ops))
        bottom = np.zeros(len(ops))
        
        # Create stacked bars for each slice within each OP
        for i, op in enumerate(ops):
            op_data = grouped[grouped['OP'] == op].sort_values('slice')
            if i == 0:
                max_y = op_data.cell_count.sum()
            if op_data.cell_count.sum() > max_y:
                max_y = op_data.cell_count.sum()

            for j, (_, row) in enumerate(op_data.iterrows()):
                slice_name = row['slice']
                cell_count = row['cell_count']
                
                # Plot bar for this slice - height = cell_count (not 1)
                bar = ax.bar(x_pos[i], cell_count, bottom=bottom[i],
                            color=slice_color_map[slice_name],
                            alpha=0.7)
                
                # Add cell count text on the bar
                # ax.text(x_pos[i], bottom[i] + cell_count/2, str(cell_count),
                #     ha='center', va='center', fontweight='bold')
                
                bottom[i] += cell_count  # Stack next slice on top
        
        # Customize plot
        ax.set_xlabel('Patients')
        ax.set_ylabel('Cells (count)')  # Y-axis shows cell count
        ax.set_xticks(x_pos, ['' for s in range(len(x_pos))])
        if max_y%2 == 1:
            max_y+=1
        ax.set_yticks(np.linspace(0, max_y, 5), \
                      [int(a) if i%2 == 0 else '' for i, a in enumerate(np.linspace(0, max_y, 5))])
        if save_dir:
            plt.savefig(f'{save_dir}hist_op_cell_count_{data_type}{day}.svg', dpi=300,
                        format='svg', bbox_inches='tight')
        
        plt.show()

## PLOT PH - Osmol

def plot_ph_osmol(save_dir = False):
    '''
    loads the data and plots pH and Osmol plots
    '''
    # setting up
    col = figure_colors()['slice']
    param_dict = {'osmol': ['Osmolarity (mOsm/L)', np.linspace(290, 320, 5), [290, '', 305 , '', 320]],
                    'pH':['pH', np.linspace(6.9, 7.5, 5), [6.9, '', 7.2, '', 7.5]]}

    data = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/solutions_QC.xlsx')
    filtered_data = data[(data['temp'] == 'RT') & (data['lid'] == 'closed')]
    filtered_data = filtered_data[(filtered_data.time == 0) | (filtered_data.time > 16)]
    indx_remove = filtered_data[(filtered_data['OP'] == 'OP240417') & ((filtered_data['time'] == 17))].index
    filtered_data = filtered_data.drop([9, 15]).reset_index(drop = True)
    filtered_data['time_binary'] = filtered_data['time'].apply(lambda x: 0 if x == 0 else 1)

    # Separating the data into 'Ctrl' and 'high K'
    ctrl_data = filtered_data[filtered_data['solution'] == 'Ctrl']
    high_k_data = filtered_data[filtered_data['solution'] == 'high K (+ 8mM)']

    # Grouping the filtered data by 'OP' column
    ctrl_groups = ctrl_data.groupby('OP')
    high_k_groups = high_k_data.groupby('OP')

    fig, ax = plt.subplots(1, 2, figsize=(9/2.54, 4/2.54))
    # plotting
    for a, param in enumerate(['osmol', 'pH']):

        # Plotting 'Ctrl' data
        for name, group in ctrl_groups:
            ax[a].scatter(group.time, group[param], color = col['Ctrl'], zorder = 5, alpha = 0.5)
            ax[a].scatter(group.time, group[param], edgecolor = col['Ctrl'], facecolor = 'none', zorder = 5, alpha = 0.9)
            ax[a].plot(group.time, group[param], linestyle='-', alpha = 0.5, color = col['Ctrl'])

        # Plotting 'high K' data
        for name, group in high_k_groups:
            ax[a].scatter(group.time + 43, group[param], color = col['high K'], zorder = 5, alpha = 0.5, marker='x')
            ax[a].scatter(group.time + 43, group[param], edgecolor = col['high K'], facecolor = 'none', zorder = 5, alpha = 0.9)
            ax[a].plot(group.time + 43, group[param], linestyle='-', alpha = 0.5, color = col['high K'])

        # Calculating the average pH for each time point within each group
        ctrl_avg = ctrl_data.groupby('time_binary')[param].mean().reset_index()
        high_k_avg = high_k_data.groupby('time_binary')[param].mean().reset_index()

        # # Plotting average data
        ax[a].scatter([0, 23], ctrl_avg[param], marker = '_', s = 75, c = 'k', alpha = 0.7, zorder = 10)
        ax[a].scatter([43, 43 + 23], high_k_avg[param], marker='_', s = 75, c = 'k', alpha = 0.7, zorder = 10)

        ax[a].set_xticks(ticks = [-7, 0, 23, 43, 66, 73], labels = ['', 'pre', 'post CTR', 'pre', 'post HiK', ''])
        ax[a].set_ylabel(f'{param_dict[param][0]}')
        ax[a].set_yticks(ticks = param_dict[param][1], labels = param_dict[param][2])
        # plt.title(param + ' over time (RT and closed lid)')

    if save_dir:
        plt.savefig(f'{save_dir}/QC_sum.svg',
                    format = 'svg', bbox_inches = 'tight', dpi = 300)
    plt.show()

def plot_failure_rates_repatch(df, save_dir, w_cm = 16, h_cm = 5):
    '''
    plots the failure rates accross the repatch dataset
    '''
    df = df.copy()
    if 'num_swps' not in df.columns:
        return print('reduqred columns not in dataframe')

    failure_rates = 100*(df.num_swps - df.num_swps_events) / df.num_swps
    df.insert(0, 'failure_rates', failure_rates)

    cols = ['failure_rates', 'num_swps', 'num_swps_events']

    data_type = 'cell_ID_new'

    ax_elements = {'failure_rates' : (0,0),
                   'num_swps': (0,1),
                   'num_swps_events': (0,2)}
    df = df.sort_values(['treatment', 'day'])

    titles_dict = pl_intr.dict_for_plotting_conn()
    colors = figure_colors()['cell_ID_new']

    fig, ax = plt.subplots(1, 3, figsize = (w_cm /2.54, h_cm/2.54))
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
                                
                x = np.linspace(0.65+k, 1.35+k, len(df_plot))
                x_plot.append(x)
                ax[a].scatter(x, df_plot[param], alpha = 0.3,
                            color = colors[tr + day], zorder = 5,
                            s = 20)
                ax[a].scatter(x, df_plot[param], alpha = 0.75,
                            edgecolor = colors[tr + day],
                            facecolors = 'none', s = 20)
                if 'cell' not in data_type:
                    ax[a].boxplot(df_plot[param], positions = [k + 1], zorder = 10, widths = 0.5)

                    label_.append(day + tr)
                    num_cels[tr + ' ' + day] = len(df_plot)
                else:
                    avg = np.mean(df_plot[param])
                    ax[a].text(0.75 + k, (1.05*avg), str(round(avg,2)),
                                    size = 11, zorder = 11)
                    ax[a].scatter(1 + k, avg, color = '#430704',
                                  marker = '_', s = 300, zorder = 10)

                    print(f'{param} avg for {tr} {day} {np.mean(df_plot[param])}')
                    print(f' n = {len(df_plot[param])}')
            if k in [1,3.5,5] and 'cell' in data_type:
                for c, cell in enumerate(df_plot['connection_ID']):
                    if k == 1 and len(df_plot['connection_ID'])-1 == c:
                        print(str(c+1) + 'Gencho')
                    x1 = [x_plot[0][c], x[c]]
                    y = df[param][df['connection_ID'] == cell]

                    ax[a].plot(x1, y, '-', color = colors[tr+'D2'], alpha = 0.5,
                               linewidth = 1.5, zorder = 0)
        
            ax[a].set_xticks(ticks = [1,2,3.5,4.5], labels = ['pre', 'post', 'pre', 'post'])

            ax[a].set_ylabel(titles_dict[param][0] + ' (' + titles_dict[param][1] + ')')
            ax[a].set_yticks(ticks = titles_dict[param][2], labels = titles_dict[param][3])
                
    if save_dir:
        plt.savefig(f'{save_dir}_failure_rates_{data_type}.svg',
                    format='svg', bbox_inches='tight', dpi=300)
    plt.show(fig)

def plot_failure_rates_slice(df, dest_dir, w_cm = 16, h_cm = 5):
    '''
    plots the failure rates accross the repatch dataset
    '''
    if 'num_swps' not in df.columns:
        return print('required columns not in dataframe')
    df = df.copy()
    failure_rates = 100*(df.num_swps - df.num_swps_events) / df.num_swps
    df.insert(0, 'failure_rates', failure_rates)

    cols = ['failure_rates', 'num_swps', 'num_swps_events']

    data_type = 'slice'

    ax_elements = {'failure_rates' : (0,0),
                   'num_swps': (0,1),
                   'num_swps_events': (0,2)}
    df = df.sort_values(['treatment', 'day'])
    if 'TTX' in sorted(df.treatment.unique()):
        treatments = ['Ctrl', 'high K', 'TTX']
    else:
        treatments = ['Ctrl', 'high K']

    titles_dict = pl_intr.dict_for_plotting_conn()
    colors = figure_colors()[data_type]

    fig, ax = plt.subplots(1, 3, figsize = (w_cm /2.54, h_cm/2.54))
    ax = ax.flatten()
    for a, param in enumerate(cols):
        
        label_ = []
        num_cels = {}
        ks, labels_ = [], []
        for i, tr in enumerate(treatments):
            x_plot =[]
            for j, day in enumerate(sorted(df.day.unique())):
                k = j + 2.5*i
                df_plot = df[(df['treatment'] == tr) &\
                              (df['day'] == day)]
                                
                x = np.linspace(0.65+k, 1.35+k, len(df_plot))
                x_plot.append(x)
                ax[a].scatter(x, df_plot[param], alpha = 0.3,
                            color = colors[tr + day], zorder = 5,
                            s = 20)
                ax[a].scatter(x, df_plot[param], alpha = 0.75,
                            edgecolor = colors[tr + day],
                            facecolors = 'none', s = 20)
                if 'cell' not in data_type:
                    ax[a].boxplot(df_plot[param], positions = [k + 1], zorder = 10, widths = 0.5)

                    label_.append(day + tr)
                    num_cels[tr + ' ' + day] = len(df_plot)
                else:
                    avg = np.mean(df_plot[param])
                    ax[a].text(0.75 + k, (1.05*avg), str(round(avg,2)),
                                    size = 11, zorder = 11)
                    ax[a].scatter(1 + k, avg, color = 'k', marker = '_', s = 300, zorder = 10)

                print(f'{param} avg for {tr} {day} {np.mean(df_plot[param])}')
                print(f' n = {len(df_plot[param])}')
                ks.append(k+1)
                labels_.append(f'{day} {tr}')
            ax[a].set_xticks(ticks = ks, labels = labels_)

            ax[a].set_ylabel(f'{titles_dict[param][0]}{titles_dict[param][1]}')
            ax[a].set_yticks(ticks = titles_dict[param][2], labels = titles_dict[param][3])

    date = str(datetime.date.today())
    if dest_dir:
        plt.savefig(dest_dir + date + '_plot_all_SLICE_MEAN_sup_V2.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.show(fig)

def plt_washin_trace(fn, chan, swps_keep = 'all', start_point_cut = False, \
                     end_point_cut = False, save_dir = False, ticks = False):
    ''' 
    arguemnts : fn - filename, chan - channel number from 1 to 8
    visualization of recordings with selection of channels and parts to cut
    accepts int for channel
    start_point_cut, end_point_cut - int, what aprt to cut from the trace
    '''
    trace = pyabf.ABF(fn)
    # fixing the naming
    if '_Ipatch' in trace.adcNames:
            trace.adcNames[trace.adcNames.index('_Ipatch')] = 'Ch1'
    if 'IN0' in trace.adcNames:
        trace.adcNames[trace.adcNames.index('IN0')] = 'Ch1'
    
    # channels are always consequtive nums in pyabf
    chan_name = 'Ch' + str(chan)
    if len(trace.channelList) < 8:
        chan = trace.channelList[trace.adcNames.index(chan_name)]
    else:
        chan = chan - 1

    data_long = trace.data[chan]
    swp_len = int(trace.sweepLengthSec * trace.sampleRate)
    data_wide = data_long.reshape(trace.sweepCount, swp_len) # each sweep is a row
    
    if isinstance(swps_keep, str):
        if swps_keep == 'all':
            data_plot = hcf.reshape_data(data_wide, 'all', start_point_cut, end_point_cut) 
    else:
        data_plot = hcf.reshape_data(data_wide, swps_keep, start_point_cut, end_point_cut)

    fig, ax = plt.subplots(1,1, figsize = (6/2.54, 6/2.54))
    x = np.linspace(0, len(data_plot), len(data_plot))/trace.sampleRate
    ax.plot(x, data_plot)
    for i in range(1, len(swps_keep)):
        vline = i * trace.sweepLengthSec
        if i == 1:
            ax.axvline(x = vline, color='red', linestyle='--', linewidth=1, label = 'sweep end')
        else:
            ax.axvline(x = vline, color='red', linestyle='--', linewidth=1)

    if ticks:
        ax.set_yticks(np.linspace(-60, 50, 5), labels = [-60, '', -5, '', 50])
        ax.set_xlabel('Time (sec)')
        ax.set_ylabel('Membrane potential \n (mV)')

    if save_dir:
        plt.savefig(f'{save_dir}_wash_in_ap_trace.svg',
                    format='svg', bbox_inches='tight', dpi=300)
    plt.show(fig)

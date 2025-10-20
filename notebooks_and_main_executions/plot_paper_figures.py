import os
import pandas as pd
import numpy as np
import datetime
import pyabf
import matplotlib.pyplot as plt
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
import ephys_analysis.funcs_human_characterisation as hcf
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
                    'post CTR': '#ff7f00'}}# 'MediumPurple2'}}
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
                 background_dot_size = 3, w_cm = 14.5, h_cm = 7):
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
                 background_dot_size = 3, w_cm = 14, h_cm = 5.2):
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
                 background_dot_size = 3, w_cm = 6.7, h_cm = 6):
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

# %%

#start from fig 2 --> decide size of text, itcks,labels, ratios
data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
    'paper_figs_collected_checked/data/'
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
        'paper_figs_collected_checked/figures/draft4_figs/fig_parts/'

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
# df_slice_short_inc = pd.read_excel(data_dir + '10.07.25_slice_incubation_only_no_move.xlsx')
# df_slice_short_inc = df_slice_short_inc[df_slice_short_inc['hrs_after_OP'] < 35]
# df_slice_short_inc = df_slice_short_inc[df_slice_short_inc['hrs_incubation'] > 15.9]
# df_slice_short_inc.reset_index(drop = True, inplace = True)
# df_slice_short_inc = df_slice_short_inc.loc[[i for i, sl in enumerate(df_slice_short_inc.slice) if len(sl) <= 2], :]
# intr_params_inc_only(df_slice_short_inc, destination_dir)
# firing_props_inc_only_V2(df_slice_short_inc, 'SE', destination_dir)

# cell_IDs_inc = {'Ctrl': '25220S3c3',
#             'high K': '25220S2c3'}
# plot_example_firing('inc_only', df_slice_short_inc, destination_dir, cell_IDs_inc)



# VERY UNREADALE CODE PROBABLY 
# TO do - check carefully the proper tables that are the sum,m,ary and don't use them!!!
# do not use the presaves summary tables

# # for sag
# work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/'

# df = df_slice
# w_cm = 8
# h_cm = 8

# # SAG
# # example_cell = '24201S2c5'
# fn = work_dir + 'OP240201' + '/'+ '24201018.abf'
# channel = 5
# w_cm = 8
# h_cm = 8

# onset = 2624
# offset = 22624
# charact_data = hcf.load_traces(fn)
# inj = hcf.get_inj_current_steps(fn)
# tau_all, capacitance_all, mcs, V65s, RMPs_char = hcf.get_hyperpolar_param(charact_data, [channel], inj)

# key = 'Ch' + str(channel)
# ch_data = charact_data[key][0]

# fig = plt.figure(figsize = (w_cm/2.54, h_cm/2.54))
# i = 0
# swp = list(ch_data[:,i])
# bl = np.median(ch_data[0:onset-20,i])

# ss = np.median(ch_data[15_000:20_000, i])
# plt.plot(ch_data[:,i], c = 'darkgrey')

# plt.xticks(ticks = [0, 10_000, 20_000], labels = [0, 10_000/20, 20_000/20])
# plt.axhline(y = ss, c = 'darkgrey', linestyle = '--')
# plt.axhline(y = min(ch_data[:,i]), c = 'darkgrey', linestyle = '--')
# plt.axhline(y = bl,c = 'darkgrey', linestyle = '--')
# plt.axvline(x = onset + 2000, c = 'darkgrey', linestyle = '--')
# # plt.savefig(destination_dir + 'sag.svg', \
# #             format = 'svg', bbox_inches = 'tight', dpi = 1000)


# # AP TH
# w_cm = 8
# h_cm = 8
# example_cell = '24201S1c8'
# df_plot = df[df['cell_ID'] == example_cell]
# channel = df_plot.cell_ch.values[0]
# fn = work_dir + df_plot.OP.values[0] + '/' + df_plot.filename.values[0]

# inj = hcf.get_inj_current_steps(fn)
# max_spikes = hcf.get_max_spikes(charact_data, [channel])
# first_spikes, peaks_all, spike_counts_all, first_spiking_sweeps_all = hcf.get_ap_param_for_plotting(charact_data, [channel], inj, max_spikes)
# AP_all, THloc_all, TH_all, = hcf.get_ap_param(charact_data, [channel], inj, max_spikes)[1:4]

# key = 'Ch' + str(channel)
# ch_data = charact_data[key][0]

# fig = plt.figure(figsize = (w_cm/2.54, h_cm/2.54))
# plt.plot(AP_all[0])
# plt.scatter(THloc_all[0] ,TH_all[0], color = 'red')
# plt.scatter(200,peaks_all[0][first_spiking_sweeps_all[0], 1, 2], color = 'green')
# # plt.savefig(destination_dir + 'TH.svg', \
# #             format = 'svg', bbox_inches = 'tight', dpi = 1000)



# # REHOBASE PLOT
# w_cm = 8
# h_cm = 8

# fn = work_dir + df_plot.OP.values[0] + '/'+ '24201002.abf'
# channel = 8

# rheos, THs, THs_in_trace, swps = hcf.get_rheobase_from_ramp(fn, [channel])
# ramp_abf = pyabf.ABF(fn)
# ramp_dict = hcf.load_traces(fn)

# fig, ax = plt.subplots(2,1,figsize = (w_cm/2.54, h_cm/2.54))
# swp = swps[0]
# start = 12_500
# end = THs_in_trace[0] + 500
# ramp = ramp_dict['Ch' + str(channel)][0][:, swp]
# ramp_abf.setSweep(sweepNumber = swp, channel = 0)
# ax[0].plot(ramp_abf.sweepX[start:end], ramp[start:end]) 
# ax[1].plot(ramp_abf.sweepX[start:end], ramp_abf.sweepC[start:end])
# ax[0].scatter(ramp_abf.sweepX[THs_in_trace[0]], THs[0], c = '#ff7f00')
# ax[1].scatter(ramp_abf.sweepX[THs_in_trace[0]], rheos[0], c = '#ff7f00')
# # for a in ax:
# #     a.axis('off')
# # plt.savefig(destination_dir + 'rheo.svg', \
# #             format = 'svg', bbox_inches = 'tight', dpi = 1000)


# #%% 

# # plot firing in a grid
# cell_IDs = {'pre': '24201S1c8',
#             'post CTR': '24117S2_D2c1',
#             'post HiK': '23420S2_D2c7'}
# data_type = 'slice'
# df = df_slice
# dest_dir = destination_dir

# w_cm = 15
# h_cm = 4
# work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/'

# swps = [9, 13, 17]
# # colors = trace_colors_trace()[data_type]
# colors = figure_colors()[data_type]

# fig, ax = plt.subplots(4, 3, sharex = True, sharey = 'row', 
#                         figsize = (w_cm/2.54, h_cm/2.54), \
#                     gridspec_kw = {'height_ratios': [6, 6, 6, 2]})
# fig.subplots_adjust(hspace = 0.01,  wspace = 0.01)
# # Plot examples traces

# for j, item in enumerate(cell_IDs.items()):
#     df_plot = df[df['cell_ID'] == item[1]]
#     fn = work_dir + df_plot.OP.values[0] + '/' + df_plot.filename.values[0]
#     channel = df_plot.cell_ch.values[0]
#     inj = hcf.get_inj_current_steps(fn)
#     trace = pyabf.ABF(fn)

#     if len(trace.channelList) < 8:
#         if '_Ipatch' in trace.adcNames:
#             trace.adcNames[trace.adcNames.index('_Ipatch')] = 'Ch1'
#         if 'IN0' in trace.adcNames:
#             trace.adcNames[trace.adcNames.index('IN0')] = 'Ch1'
#         channel_name = 'Ch' + str(channel)
#         channel = trace.channelList[trace.adcNames.index(channel_name)]
#     else:
#         channel = channel-1
#         channel_name = 'Ch' + str(channel+1)

#     for i, sweep in enumerate(swps):
#         trace.setSweep(sweepNumber = sweep, channel = channel)
#         # ax[0, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
#         #                 color = colors[item[0]+ '.' + str(i+1)], alpha = 0.7)
#         ax[j, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
#                         color = colors[item[0]], alpha = 0.7)
#         if i == 0:
#             ax[j, i].plot(trace.sweepX[:-8000], trace.sweepY[:-8000], \
#                         color = colors[item[0]], alpha = 0.7, label = item[0])
#         # x = trace.sweepX
#         # y = trace.sweepY

#         ax[j,i].set_ylabel(trace.sweepLabelY)
#         # ax[0,i].set_title(item[0])
#         ax[j,i].axis('off')

#         stim_inj = np.concatenate((np.zeros(2625), np.repeat(inj[sweep], 20_000), np.zeros(9375)))
#         ax[3, i].plot(trace.sweepX [:-8000], stim_inj, c = 'black')
#         # ax[1, i].set_xlabel(trace.sweepLabelX)
#         # ax[1, i].set_ylabel(trace.sweepLabelC)
#         ax[3, i].axis('off')
#         ax[0,i].set_title(str(int(inj[sweep])) + ' pA')
#     # fig.legend(loc = 'upper center', ncol = 3, bbox_to_anchor=(0.5, 1.3))
# # plt.savefig(dest_dir + data_type + 'I_O_curve_grid.svg', \
# #            format = 'svg', bbox_inches = 'tight', dpi = 1000)




# #%%
# # count types of firing cells

# def count_firing_cells(df, data_type):
#     '''
#     clasification of the firing cells based on type
#     no AP, single AP, multiple APs
#     '''

#     df_slice = pd.read_excel(data_dir + 'slice_data_temporal.xlsx')
#     df_repatch = pd.read_excel(data_dir + 'repatch_data_temporal.xlsx')

#     df = df_slice
#     data_type = 'slice'


#     firing_type = []
#     for i, spikes in enumerate(df.max_spikes):
#         if spikes == 0 or np.isnan(spikes):
#             firing_type.append('no_APs')
#         elif spikes == 1:
#             firing_type.append('single_AP')
#         elif spikes > 1:
#             firing_type.append('multiple_APs')
#     df.insert(len(df.columns), 'firing_type', firing_type)


#     w_cm = 8
#     h_cm = 8
#     colors = figure_colors()['slice']


#     count_firing_cells = {}
#     labels_, ticks_ = [], []
#     fig, ax = plt.subplots(1,1, figsize = (w_cm/2.54, h_cm/2.54))
#     for j, day in enumerate(df.day.unique()):
#         for i, treatment in enumerate(['Ctrl', 'high K']):
#             condition_df = df[(df['day'] == day) & (df['treatment'] == treatment)]
#             for f, f_type in enumerate(['no_APs', 'single_AP', 'multiple_APs']):
#                 df_plot = df[(df['day'] == day) &\
#                                 (df['treatment'] == treatment) &\
#                                 (df['firing_type'] == f_type)]
#                 x = f + 0.2 * (i + 2*j)
#                 print(x)
#                 ax.bar(x, len(df_plot) / len(condition_df), 
#                     width = 0.1, color = figure_colors()[data_type][treatment + day])
#                 labels_.append(treatment + day + f_type)
#                 ticks_.append(x)
#                 count_firing_cells[treatment + day + f_type] = len(df_plot)

#     ax.set_xticks(ticks = ticks_, labels = labels_, rotation = 45)
#     plt.show()
#     plt.close()










# #%% might be useful

# # filtered_df = adult_df_slice_all[~adult_df_slice_all['cell_IDs_match'].isin(adult_df_slice_no_repeats['cell_IDs_match'])]


# # df_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
# # human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
# # patchers_dict = {'Rosie': 'data_rosie/', 'Verji':'data_verji/'}

# # unused_df = pd.read_excel(df_folder + 'slice_incubation_only_no_move.xlsx')
# # slice_all = pd.read_excel(df_folder + 'slice_all.xlsx')

# # df_dict = {'slice_all.xlsx' : slice_all, 
# #            'slice_incubation_only_no_move.xlsx': unused_df}

# # for key, df in df_dict.items():
# #     df = df.drop(columns=['sag'])
# #     sag = []
# #     for i in range(len(df)):
# #         patcher = df.patcher[i]
# #         fn = human_dir + patchers_dict[patcher] + df.OP[i] + '/' + df.filename[i]
# #         chan = [int(df.cell_ch[i])]
# #         sag.append(hcf.sag(fn, chan)[0][0])

# #     df.insert(len(df.columns), 'sag', sag)
# #     df.to_excel(df_folder + key)
# #     df.to_csv(df_folder + key[:-5] + '.csv')



# ## PLAYGROUND
# df_sllice_short_inc = pd.read_excel(data_dir + 'slice_all.xlsx')
# df_sllice_short_inc = df_sllice_short_inc[df_sllice_short_inc['hrs_after_OP'] < 31]
# iff_df = df_sllice_short_inc
# background_dot_size = 3
# error_type = 'CI'
# dest_dir = destination_dir
# def firing_props_inc_only(iff_df, error_type, dest_dir,
#                  background_dot_size = 3, w_cm = 14.7, h_cm = 12):
#     """
#     Plots the initial firing frequency and number of APs for only incubated slices
#     """
#     df_summary = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results'+ \
#                             '/human/paper_figs_collected_checked/stats/'+ \
#                             'sum_data_CIs_SEs_for_plotting/firing_plot_CIs_inc_only.xlsx')
#     iff_df = iff_df.sort_values(['treatment', 'day'])
#     num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(iff_df)
#     dt_sum = 'slice'
#     colors = figure_colors()[dt_sum]

#     dv_dict = {'IFF' : [iff_indx, 'Initial firing frequency (AP#1 to AP#2)',
#                 'Initial firing \n frequency (Hz)'],
#                 'num_aps' : [num_aps_indx, 'Number of fired action potentials', 'AP frequency (Hz)']}


#     for v, dv in enumerate(['num_aps', 'IFF']):
#         fig, ax = plt.subplots(1, 1, figsize=(w_cm / 2.54, h_cm/2.54))
#         counts_firing = {}
#         for k, treatment in enumerate(iff_df.treatment.unique()):
            
#             for d, day in enumerate(iff_df['day'].unique()):
#                 if day == 'D1':
#                     day_df = iff_df[iff_df['day'] == day]
#                 else:
#                     day_df = iff_df[(iff_df['treatment'] == treatment) & \
#                                 (iff_df['day'] == day)]
#                 print(treatment, day)

#                 avgs, inj, counts = [], [], []
#                 for i, col in enumerate(dv_dict[dv][0][5:5+13]):
                    
#                     data = day_df.iloc[:,col]
#                     inj_val = int(day_df.columns[col][2:day_df.columns[col].rfind('pA')])

#                     # removing non-firing cells
#                     data.replace(0, np.nan, inplace = True)
#                     data.dropna(inplace = True)
#                     avg = np.mean(data)

#                     if error_type == 'SE':
#                         x = np.linspace(0.75 + i , 1.25 + i, len(data))
#                         # ax.scatter(x, data, alpha = 0.15, s = background_dot_size,
#                         #             c = colors[treatment + day])
#                         sum_data = df_summary[(df_summary['var_firing'] == dv) & \
#                                         (df_summary['treatment'] == treatment) & \
#                                         (df_summary['day'] == d) & \
#                                         (df_summary['inj_current'] == inj_val)]
#                         # if there is enough data for SE
#                         if len(sum_data['SE']) > 0:
#                             sem = sum_data['SE'].values[0]
#                             ax.errorbar(i + 1, avg, yerr = sem, color = colors[treatment + day])

#                     inj.append(inj_val)
#                     counts.append(len(data))
#                     avgs.append(avg)

#                 counts_firing[treatment + day] = counts

#                 sum_data_tr = df_summary[(df_summary['var_firing'] == dv) & \
#                                         (df_summary['treatment'] == treatment) & \
#                                         (df_summary['day'] == d)]

#                 ax.plot(range(1, len(inj)+1), avgs, label = day, color = colors[treatment + day])

#                 if error_type == 'CI':
#                     start_summary = len(inj) - len(sum_data_tr) + 1
#                     # ax.scatter(range(1, len(inj)+1), avgs, label = day, color = colors[treatment + day])
#                     ax.fill_between(range(start_summary, len(inj)+1), sum_data_tr['CI_L'], 
#                                         sum_data_tr['CI_U'], color = colors[treatment + day],alpha = 0.3)
#                 #     if dv == 'num_aps':
#                 #         ax.set_yticks(ticks = [0, 10, 20, 30])
#                 #     else:
#                 #         ax.set_yticks(ticks = [0, 20, 40, 60, 80, 100])
#                 # else:
#                 #     if dv == 'num_aps':
#                 #         ax.set_yticks(ticks = [0, 10, 20, 30, 40, 50])
#                 #     else:
#                 #         ax.set_yticks(ticks = [0, 40, 80, 120, 160, 200])

#                 ax.set_ylabel(dv_dict[dv][2])

#         if dv == 'num_aps':
#             ax.set_yticks(ticks = [0, 5, 10, 15, 20, 25])
#         else:
#             ax.set_yticks(ticks = [0, 20, 40,60, 80, 100])

#         ticks_ = np.arange(1, 14, 2)
#         labels_ = [inj[a] for a in np.arange(0,  13, 2)]
        
#         ax.set_xticks(ticks_)
#         ax.set_xticklabels(labels_, rotation=30)

#         ax.set_xlabel('Current injection (pA)')


#         plt.subplots_adjust(hspace = 0.12, wspace = 0.15)
#         # plt.savefig(dest_dir + error_type + dv + '_vs_curernt_inj_inc_only.svg', \
#         #            format = 'svg', bbox_inches = 'tight', dpi = 1000)
#         plt.close(fig)





# #%% 


# # # trying to figure out firing error
# # iff_df = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/inc_only_for_python_plot_delete_soon2.xlsx')

# # treat_iff = ['CTR' if a == 0 else 'HiK' for a in iff_df.treatment_r]
# # iff_df.drop('treatment', axis = 1, inplace = True)
# # iff_df.insert(0, 'treatment', treat_iff)

# # iff_df = iff_df.sort_values(['treatment', 'day'])
# # num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(iff_df)
# # dt_sum = 'slice'
# # colors = figure_colors()[dt_sum]

# # dv_dict = {'IFF' : [iff_indx, 'Initial firing frequency (AP#1 to AP#2)',
# #             'Initial firing \n frequency (Hz)'],
# #             'num_aps' : [num_aps_indx, 'Number of fired action potentials', 'AP frequency (Hz)']}

# # df_summary = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/stats/sum_data_CIs_SEs_for_plotting/firing_plot_temporary_inc_only_temporal.xlsx')

# # treat = ['CTR' if a == 0 else 'HiK' for a in df_summary.treatment_r]
# # df_summary.insert(0, 'treatment', treat)

# # treat_iff = ['CTR' if a == 'Ctrl' else 'HiK' for a in iff_df.treatment]
# # iff_df.drop('treatment', axis = 1, inplace = True)
# # iff_df.insert(0, 'treatment', treat_iff)

# # for v, dv in enumerate(['num_aps', 'IFF']):
# #     fig, ax = plt.subplots(1, 1, figsize=(w_cm / 2.54, h_cm/2.54))
# #     counts_firing = {}
# #     for k, treatment in enumerate(iff_df.treatment.unique()):
        
# #         day_df = iff_df[(iff_df['treatment'] == treatment)]
      
# #         avgs, inj, counts = [], [], []
# #         for i, col in enumerate(dv_dict[dv][0][5:5+13]):
            
# #             data = day_df.iloc[:,col]
# #             # inj_val = int(day_df.columns[col][2:day_df.columns[col].rfind('pA')])
# #             # when dataframes comes from R:
# #             inj_val = int(day_df.columns[col][3:day_df.columns[col].rfind('pA')])

# #             # removing non-firing cells
# #             data.replace(0, np.nan, inplace = True)
# #             data.dropna(inplace = True)
# #             avg = np.mean(data)

# #             if error_type == 'SE':
# #                 x = np.linspace(0.75 + i , 1.25 + i, len(data))

# #                 sum_data = df_summary[(df_summary['var_firing'] == dv) & \
# #                                 (df_summary['treatment'] == treatment) & \
# #                                 (df_summary['inj_current'] == inj_val)]
# #                 # if there is enough data for SE
# #                 if len(sum_data['SE']) > 0:
# #                     sem = sum_data['SE'].values[0]
# #                     ax.errorbar(i + 1, avg, yerr = sem, color = colors[treatment + 'D2'])

# #             inj.append(inj_val)
# #             counts.append(len(data))
# #             avgs.append(avg)

# #         counts_firing[treatment + 'D2'] = counts

# #         sum_data_tr = df_summary[(df_summary['var_firing'] == dv) & \
# #                                 (df_summary['treatment'] == treatment)]

# #         ax.plot(range(1, len(inj)+1), avgs, label = 'D2', color = colors[treatment + 'D2'])

# #         if error_type == 'CI':
# #             start_summary = len(inj) - len(sum_data_tr) + 1
# #             ax.fill_between(range(start_summary, len(inj)+1), sum_data_tr['CI_L'], 
# #                                 sum_data_tr['CI_U'], color = colors[treatment + 'D2'],alpha = 0.3)

# #         ax.set_ylabel(dv_dict[dv][2])

# #     if dv == 'num_aps':
# #         ax.set_yticks(ticks = [0, 5, 10, 15, 20, 25])
# #     else:
# #         ax.set_yticks(ticks = [0, 20, 40,60, 80, 100])

# #     ticks_ = np.arange(1, 14, 2)
# #     labels_ = [inj[a] for a in np.arange(0,  13, 2)]
    
# #     ax.set_xticks(ticks_)
# #     ax.set_xticklabels(labels_, rotation=30)

# #     ax.set_xlabel('Current injection (pA)')




# # data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/sanity_check/'

# # melt_df = pd.read_excel(data_dir + 'slice_melt_df_temporal.xlsx')
# # iff_df = pd.read_csv(data_dir + 'slice_data_temporal.csv')


# # fig, ax = plt.subplots(1, 1)
# # for treatment in melt_df.treatment_r.unique():
# #     avgs = []
# #     for i, inj in enumerate(melt_df.inj_current.unique()):
# #         df_plot = melt_df[(melt_df['treatment_r'] == treatment) & \
# #                           (melt_df['inj_current'] == inj)]
# #         avg = df_plot['VAL'].mean()
# #         sem = df_plot['VAL'].std() / np.sqrt(len(df_plot))

# #         x = np.linspace(0.75 + inj , 1.25 + inj, len(df_plot))
# #         ax.scatter(x, df_plot['VAL'], alpha = 0.15)
# #         ax.errorbar(inj + 1, avg, yerr = sem, label = treatment)
# #         avgs.append(df_plot['VAL'].mean())

# #     ax.plot(melt_df.inj_current.unique(), avgs, label = treatment)

# # # data pre-processing for plotting
# # iff_df = iff_df.sort_values(['treatment'])
# # num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(iff_df)

# # dv_dict = {'IFF' : [iff_indx, 'Initial firing frequency (AP#1 to AP#2)',
# #             'Initial firing \n frequency (Hz)'],
# #             'num_aps' : [num_aps_indx, 'Number of fired action potentials', 'AP frequency (Hz)']}

# # fig, ax = plt.subplots(1, 1)
# # dv = 'num_aps'
# # for k, treatment in enumerate(iff_df.treatment.unique()):
    
# #     avgs, inj= [], []
# #     for i, col in enumerate(dv_dict[dv][0][5:5+13]):
# #         data = iff_df[(iff_df['treatment'] == treatment)].iloc[:,col]
# #         # if error 3 when from R, 2 when from python. try both
# #         inj_val = int(iff_df.columns[col][2:iff_df.columns[col].rfind('pA')])

# #         # removing non-firing cells
# #         data.replace(0, np.nan, inplace = True)
# #         data.dropna(inplace = True)
        
# #         avg = np.mean(data)
# #         sem = np.std(data.values) / np.sqrt(len(data))

# #         x = np.linspace(0.75 + inj_val , 1.25 + inj_val, len(data))
# #         ax.scatter(x, data, alpha = 0.15)
# #         ax.errorbar(inj_val + 1, avg, yerr = sem, label = treatment)
# #         avgs.append(avg)
# #         inj.append(inj_val)

# #     ax.plot(inj, avgs, label = treatment)


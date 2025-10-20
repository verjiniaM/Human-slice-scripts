#%%
import numpy as np
import pandas as pd
import importlib
import matplotlib.pyplot as plt
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
importlib.reload(pl_intr)
plt.style.use('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/Human-slice-scripts/style_plot_paper.mplstyle')


# FUNCTTIONS
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

def plot_slice_params(df, params, dest_dir, w_cm = 7.3, h_cm = 30):
    """
    Plots parameters for different treatments and days 
    """
    # for single plots - w_cm = 6.5, h_cm = 5.7
    mask_exclude = (df['treatment'] != 'high K') & (df['treatment'] != 'Ctrl')
    df = df[~mask_exclude]
    df = df.sort_values(['treatment', 'day'])
    titles_all = pl_intr.dict_for_plotting()
    titles_dict = {key: titles_all[key] for key in params if key in titles_all}

    dict_label = {'Ctrl':'CTR', 'high K':'HiK'}
    colors = figure_colors()['slice']
    fig2, axs = plt.subplots(len(params), 1, figsize = (w_cm/2.54, h_cm / 2.54),
                             sharex = True)

    for u, (param, item) in enumerate(titles_dict.items()):
        label_, data_boxplot = [], []
        num_cels = {}
        for i, tr in enumerate(df.treatment.unique()):
            for j, day in enumerate(df['day'].unique()):
                k = j + 2.5*i
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                df_plot = df_plot[~df_plot[param].isna()]  # remove NaN values
                x = np.linspace(0.65+k, 1.35+k, len(df_plot))

                axs[u].scatter(x, df_plot[param], alpha = 0.3,
                            color = colors[tr + day], zorder = 5) # dots inside
                axs[u].scatter(x, df_plot[param], alpha = 0.75,
                            edgecolor = colors[tr + day],
                            facecolors = 'none')

                axs[u].boxplot(df_plot[param], positions = [k + 1], zorder = 10, widths = 0.5)

                label_.append(day + '\n' + dict_label[tr])
                num_cels[tr + ' ' + day] = len(df_plot)
                data_boxplot.append(df_plot[param])
                axs[u].set_ylabel(item[0] + '\n(' + item[1] + ')')
                axs[u].set_yticks(ticks = item[2], labels = item[3])

    axs[u].set_xticks(ticks = [1,2,3.5,4.5], labels = ['pre', 'post', 'pre', 'post'])
    axs[0].set_title('CTR                        HiK')

    plt.savefig(dest_dir  + 'slice_intr_plot' + '.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 300)

def plot_repatch_params(df, params, dest_dir, w_cm = 7.3, h_cm = 30):
    '''
    plots the repatch data
    no coloring based on age
    cold be added but plot become too much
    '''
    df = df.sort_values(['treatment', 'day'])
    titles_all = pl_intr.dict_for_plotting()
    titles_dict = {key: titles_all[key] for key in params if key in titles_all}
    colors = figure_colors()['cell_ID_new']

    fig2, axs = plt.subplots(len(params),1, figsize = (w_cm/2.54, h_cm / 2.54),
                             sharex = True)

    for u, (param, item) in enumerate(titles_dict.items()):
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

                axs[u].scatter(x, df_plot[param], c=colors[tr + day], zorder = 5, alpha = 0.3)
                axs[u].scatter(x, df_plot[param], edgecolor = colors[tr + day], alpha = 0.75,
                                   facecolors = 'none')
                axs[u].scatter(1 + k, mean, color = 'k', marker = '_', s = 300, zorder = 10)

                if param == 'rheo_ramp_c':
                    axs[u].text(0.75 + k, (1.05*mean), str(int(round(mean,0))),
                                    size = 11, zorder = 11)
                elif param == 'sag':
                    axs[u].text(0.75 + k, (1.05*mean), str(round(mean,2)),
                                    size = 11, zorder = 11)
                else:
                    axs[u].text(0.75 + k, (1.05*mean), str(round(mean,1)),
                                    size = 11, zorder = 11)

                x_plot.append(x)
                if k in [1,3.5,5]:
                    for c, cell in enumerate(df_plot['cell_ID_new']):
                        x1 = [x_plot[0][c], x[c]]
                        y = df[param][df['cell_ID_new'] == cell]
                        axs[u].plot(x1, y, '-', color = colors[tr + day],
                                alpha = 0.4, zorder = 1)
                if day =='D1':
                    x_ax = 'pre'
                elif day =='D2':
                    x_ax = 'post'
                label_.append(x_ax)

                # ax.text(k + 0.65, int((1 + 0.1)*np.max(df[param])), 'n = ' + str(len(df_plot)), size = 12, c = colors[tr + day])

        axs[u].set_xticks(ticks = [1,2,3.5,4.5], labels = label_)
        axs[u].set_ylabel(item[0] + '\n (' + item[1] + ')')
        axs[u].set_yticks(ticks = item[2], labels = item[3])
        # axs[u].set_title('CTR                            HiK')

    plt.savefig(dest_dir  + 'repatch_intr_plot' + '.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 1000)
    
def intr_params_inc_only(df, params, dest_dir, w_cm = 6, h_cm = 30):
    """
    Plots parameters for different treatments and days 
    """

    mask_exclude = (df['treatment'] != 'high K') & (df['treatment'] != 'Ctrl')
    df = df[~mask_exclude]
    # params = ['TH']
    df = df.sort_values(['treatment', 'day'])
    titles_all = pl_intr.dict_for_plotting()
    titles_dict = {key: titles_all[key] for key in params if key in titles_all}

    dict_label = {'Ctrl':'CTR', 'high K':'HiK'}
    colors = figure_colors()['inc_only']

    fig2, axs = plt.subplots(len(params), 1, figsize = (w_cm/2.54, h_cm / 2.54),
                             sharex = True)
    for u, (param, item) in enumerate(titles_dict.items()):
        label_, data_boxplot = [], []
        num_cels = {}
        for i, tr in enumerate(df.treatment.unique()):
            df_plot = df[(df['treatment'] == tr) & (df['day'] == 'D2')]
            k = i + 1
            df_plot = df_plot[~df_plot[param].isna()]  # remove NaN values
            x = np.linspace(0.45+k, 1.15+k, len(df_plot))

            axs[u].scatter(x, df_plot[param], alpha = 0.3,
                        color = colors[tr], zorder = 5) # dots inside
            axs[u].scatter(x, df_plot[param], alpha = 0.75,
                        edgecolor = colors[tr],
                        facecolors = 'none')

            axs[u].boxplot(df_plot[param], positions = [k + 0.8], zorder = 10, widths = 0.5)

            label_.append( dict_label[tr])
            num_cels[tr] = len(df_plot)
            data_boxplot.append(df_plot[param])

        axs[u].set_xticks(ticks = [1.8, 2.8], labels = ['post', 'post'])

        axs[u].set_ylabel(item[0] + '\n (' + item[1] + ')')
        axs[u].set_yticks(ticks = item[2], labels = item[3])

        plt.savefig(dest_dir  + 'incub_intr_plot' + '.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 1000)

def plot_hrs_after_op_vs_fir_params(df_type, treatments, df, dest_dir):

    num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(df)
    dv_dict = {
        'IFF':[iff_indx,'Initial firing frequency (Hz)'],
        'num_aps': [num_aps_indx,  'AP frequency (Hz)']}

    colors = figure_colors()[df_type]
    w_cm = 7
    h_cm = 11
    fig2, axs = plt.subplots(2, 1, figsize = (w_cm/2.54, h_cm / 2.54),
                                sharex = True)
    for i, key in enumerate(dv_dict.keys()):
        for day in df.day.unique():
            for tr in treatments:
                if len(df.day.unique()) == 1:
                    day = ''
                    df_plot = df[df.treatment == tr]
                else:
                    df_plot = df[(df.day == day) & (df.treatment == tr)]
                df_plot = df_plot.sort_values(['hrs_after_OP'])
                max_param = df_plot.iloc[:,dv_dict[key][0]].max(axis=1)

                if key == 'num_aps':
                    indx = [a for a, p in enumerate(max_param) if p < 50]
                    max_param = max_param.iloc[indx]
                    axs[i].scatter(df_plot.hrs_after_OP.iloc[indx], max_param, color = colors[tr+day])
                else:
                    axs[i].scatter(df_plot.hrs_after_OP, max_param, color = colors[tr+day])

                axs[i].set_xlabel('Time tissue extraction (hrs)')
                axs[i].set_xticks(ticks = np.linspace(0, 50 , 5), labels = [0, '', 25, '', 50])
                axs[i].set_ylabel(dv_dict[key][1])
                if key == 'num_aps':
                    axs[i].set_yticks(ticks = np.linspace(0, 50 ,5), labels = [0, '', 25, '', 50])
                else:
                    axs[i].set_yticks(ticks = np.linspace(0, 250, 5), labels = [0, '', 125, '', 250])
        
    plt.savefig(dest_dir  + df_type + 'time_vs_max_firing_plot.svg', \
    format = 'svg', bbox_inches = 'tight', dpi = 1000)

def plot_hrs_after_op_vs_fir_params_600(df_type, treatments, df, dest_dir):

    dv_dict = {'IFF':'Initial firing frequency (Hz)',
               'num_aps': 'AP frequency (Hz)'}

    colors = figure_colors()[df_type]
    w_cm = 7
    h_cm = 11
    fig2, axs = plt.subplots(2, 1, figsize = (w_cm/2.54, h_cm / 2.54),
                                sharex = True)
    for i, key in enumerate(dv_dict.keys()):
        for day in df.day.unique():
            for tr in treatments:
                if len(df.day.unique()) == 1:
                    day = ''
                    df_plot = df[df.treatment == tr]
                else:
                    df_plot = df[(df.day == day) & (df.treatment == tr)]

                indx_600 = df_plot.columns.str.contains("('600pA', " + "'" + key + "')").nonzero()[0]

                axs[i].scatter(df_plot.hrs_after_OP, df_plot.iloc[:,indx_600], color = colors[tr+day])

                axs[i].set_xlabel('Time tissue extraction (hrs)')
                axs[i].set_xticks(ticks = np.linspace(0, 50 , 5), labels = [0, '', 25, '', 50])
                axs[i].set_ylabel(dv_dict[key])
                if key == 'num_aps':
                    axs[i].set_yticks(ticks = np.linspace(0, 50 ,5), labels = [0, '', 25, '', 50])
                else:
                    axs[i].set_yticks(ticks = np.linspace(0, 250, 5), labels = [0, '', 125, '', 250])
        
    plt.savefig(dest_dir  + df_type + 'time_vs_firing_at_600pA_plot.svg', \
    format = 'svg', bbox_inches = 'tight', dpi = 1000)



data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
    'paper_figs_collected_checked/data/'
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/'+\
    'human/paper_figs_collected_checked/figures/all_params_reference/add/'
df_fns = ['10.07.25_slice_incubation_only_no_move.xlsx',
          'slice_data_temporal.xlsx',
          'repatch_data_temporal.xlsx']
df_types = ['inc_only', 'slice', 'cell_ID_new']

params = ['membra_time_constant_tau', 'resting_potential', 'cap_adj']

for df_name in df_fns:
    df = pd.read_excel(data_dir + df_name)
    if  'incubation' in df_name:
        df = df[df['hrs_after_OP'] < 35]
        df = df[df['hrs_incubation'] > 15.9]
        df.reset_index(drop = True, inplace = True)
        df = df.loc[[i for i, sl in enumerate(df.slice) if len(sl) <= 2], :]
        intr_params_inc_only(df, params, destination_dir)
        # plot_hrs_after_op_vs_fir_params(df_types[0], ['Ctrl'], df, destination_dir)
        # plot_hrs_after_op_vs_fir_params_600(df_types[0], ['Ctrl'], df, destination_dir)
    if 'slice_data' in df_name:
        mask = df.loc[(df.AP_heigth > 50) & (df.AP_halfwidth < 1.65), :]
        df.reset_index(drop=True, inplace=True)
        plot_slice_params(df, params, destination_dir)
        # plot_hrs_after_op_vs_fir_params(df_types[1], ['Ctrl'], df, destination_dir)
        # plot_hrs_after_op_vs_fir_params_600(df_types[1], ['Ctrl'], df, destination_dir)
    if 'repatch' in df_name:
        # df = pd.read_excel(data_dir + df_fns[2])
        df = df.loc[df.AP_halfwidth < 1.9, :]
        df.reset_index(drop=True, inplace=True)
        plot_repatch_params(df, params, destination_dir)
        # plot_hrs_after_op_vs_fir_params(df_types[2], ['Ctrl'], df, destination_dir)
        # plot_hrs_after_op_vs_fir_params_600(df_types[2], ['Ctrl'], df, destination_dir)
    

#%%
df = pd.read_excel(data_dir + df_fns[0])
df = df[df['hrs_after_OP'] < 35]
df = df[df['hrs_incubation'] > 15.9]
df.reset_index(drop = True, inplace = True)
df = df.loc[[i for i, sl in enumerate(df.slice) if len(sl) <= 2], :]
df_inc = df
df = pd.read_excel(data_dir + df_fns[1])

# dv_dict = {'IFF':'Initial firing frequency (Hz)',
#                'num_aps': 'AP frequency (Hz)'}

# colors = figure_colors()
# w_cm = 7
# h_cm = 11
# fig2, axs = plt.subplots(2, 1, figsize = (w_cm/2.54, h_cm / 2.54),
#                             sharex = True)
# for i, key in enumerate(dv_dict.keys()):
#     for day in ['D1', 'D2']:
#         for tr in ['Ctrl']:
#             if len(df.day.unique()) == 1:
#                 day = ''
#                 df_plot = df[df.treatment == tr]
#             else:
#                 df_plot = df[(df.day == day) & (df.treatment == tr)]
#             df_plot = df_plot.sort_values(['hrs_after_OP'])

#             indx_600 = df_plot.columns.str.contains("('600pA', " + "'" + key + "')").nonzero()[0]
#             axs[i].scatter(df_plot.hrs_after_OP, df_plot.iloc[:,indx_600], color = colors['slice'][tr+day])
            
#             df_plot_inc = df_inc[df_inc.treatment == tr]
#             indx_600_inc = df_plot_inc.columns.str.contains("('600pA', " + "'" + key + "')").nonzero()[0]
#             axs[i].scatter(df_plot_inc.hrs_after_OP, df_plot_inc.iloc[:,indx_600_inc], color = colors['inc_only'][tr])

#             axs[i].set_xlabel('Time tissue extraction (hrs)')
#             axs[i].set_xticks(ticks = np.linspace(0, 50 , 5), labels = [0, '', 25, '', 50])
#             axs[i].set_ylabel(dv_dict[key])
#             if key == 'num_aps':
#                 axs[i].set_yticks(ticks = np.linspace(0, 50 ,5), labels = [0, '', 25, '', 50])
#             else:
#                 axs[i].set_yticks(ticks = np.linspace(0, 250, 5), labels = [0, '', 125, '', 250])

# plt.savefig(destination_dir  + 'slice_all_time_vs_600_firing_plot.svg', \
# format = 'svg', bbox_inches = 'tight', dpi = 1000)
    


# PLOT hrs_after OP vs max_firing slice_ALL

# num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(df)
# num_aps_indx_inc, iff_indx_inc = pl_intr.get_num_aps_and_iff_data_culumns(df_inc)
# dv_dict = {
#     'IFF':[iff_indx,'Initial firing frequency (Hz)'],
#     'num_aps': [num_aps_indx,  'AP frequency (Hz)']}

# dv_dict_inc = {
#     'IFF':[iff_indx_inc,'Initial firing frequency (Hz)'],
#     'num_aps': [num_aps_indx_inc,  'AP frequency (Hz)']}

# colors = figure_colors()
# w_cm = 7
# h_cm = 11
# fig2, axs = plt.subplots(2, 1, figsize = (w_cm/2.54, h_cm / 2.54),
#                             sharex = True)
# for i, key in enumerate(dv_dict.keys()):
#     for day in ['D1', 'D2']:
#         for tr in ['Ctrl']:
#             if len(df.day.unique()) == 1:
#                 day = ''
#                 df_plot = df[df.treatment == tr]
#             else:
#                 df_plot = df[(df.day == day) & (df.treatment == tr)]
#             df_plot = df_plot.sort_values(['hrs_after_OP'])
#             max_param = df_plot.iloc[:,dv_dict[key][0]].max(axis=1)
#             axs[i].scatter(df_plot.hrs_after_OP, max_param, color = colors['slice'][tr+day])
            
#             df_plot_inc = df_inc[df_inc.treatment == tr]
#             max_param_inc = df_plot_inc.iloc[:,dv_dict_inc[key][0]].max(axis=1)
#             axs[i].scatter(df_plot_inc.hrs_after_OP, max_param_inc, color = colors['inc_only'][tr])

#             axs[i].set_xlabel('Time tissue extraction (hrs)')
#             axs[i].set_xticks(ticks = np.linspace(0, 50 , 5), labels = [0, '', 25, '', 50])
#             axs[i].set_ylabel(dv_dict[key][1])
#             if key == 'num_aps':
#                 axs[i].set_yticks(ticks = np.linspace(0, 50 ,5), labels = [0, '', 25, '', 50])
#             else:
#                 axs[i].set_yticks(ticks = np.linspace(0, 250, 5), labels = [0, '', 125, '', 250])

# # plt.savefig(destination_dir  + 'slice_all_time_vs_max_firing_plot.svg', \
# # format = 'svg', bbox_inches = 'tight', dpi = 1000)




# plot repatch
df = pd.read_excel(data_dir + df_fns[2])

num_aps_indx, iff_indx = pl_intr.get_num_aps_and_iff_data_culumns(df)
dv_dict = {
    'IFF':[iff_indx,'Initial firing frequency (Hz)'],
    'num_aps': [num_aps_indx,  'AP frequency (Hz)']}

colors = figure_colors()['cell_ID_new']
w_cm = 7
h_cm = 11
fig2, axs = plt.subplots(2, 1, figsize = (w_cm/2.54, h_cm / 2.54),
                            sharex = True)
for i, key in enumerate(dv_dict.keys()):
    for day in df.day.unique():
        for tr in ['Ctrl']:
            if len(df.day.unique()) == 1:
                day = ''
                df_plot = df[df.treatment == tr]
            else:
                df_plot = df[(df.day == day) & (df.treatment == tr)]
            df_plot = df_plot.sort_values(['hrs_after_OP'])
            max_param = df_plot.iloc[:,dv_dict[key][0]].max(axis=1)
            if key == 'num_aps':
                indx = [a for a, p in enumerate(max_param) if p < 50]
                max_param = max_param.iloc[indx]
                axs[i].scatter(df_plot.hrs_after_OP.iloc[indx], max_param, color = colors[tr+day])
            else:
                axs[i].scatter(df_plot.hrs_after_OP, max_param, color = colors[tr+day])

            axs[i].set_xlabel('Time tissue extraction (hrs)')
            axs[i].set_xticks(ticks = np.linspace(0, 50 , 5), labels = [0, '', 25, '', 50])
            axs[i].set_ylabel(dv_dict[key][1])
            if key == 'num_aps':
                axs[i].set_yticks(ticks = np.linspace(0, 50 ,5), labels = [0, '', 25, '', 50])
            else:
                axs[i].set_yticks(ticks = np.linspace(0, 250, 5), labels = [0, '', 125, '', 250])
    
# plt.savefig(dest_dir  + df_type + 'time_vs_max_firing_plot.svg', \
# format = 'svg', bbox_inches = 'tight', dpi = 1000)
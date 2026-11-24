import os
import re
import datetime
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
import ephys_analysis.stimulation_windows_ms as stim_win
import ephys_analysis.funcs_human_characterisation as hcf
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
import ephys_analysis.funcs_sorting as sort
from scipy.signal import butter, sosfiltfilt


plt.style.use('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/Human-slice-scripts/style_plot_paper.mplstyle')

def figure_colors():
    '''
    outputs a color dictonary
    'inc_only', 'slice', 'all', 'cell_ID_new', 'correlation'
    '''
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
                    'TTXD2': '#0e62ff'},
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

def plot_corrs(x, df_list, dir_read_corrs, save_dir_figs, params_list = []):
    '''
    x - x_axis parameter
    df_list - dictionary. 'df_name' = [df, 'color']
    save_dir_corr - where the correlations dataframe is saved
    '''
    w_cm =  18 # 6 * len(params_list)
    h_cm = 10

    titles_dict = pl_intr.dict_for_plotting_FM()
    for df_name, df_info in df_list.items():

        df = df_info[0]
        col = df_info[1]

        df_cor = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'cor_estimate')
        df_p_val = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'p_vals')
        # df_p_val = df_p_val.iloc[0:11]
        df_l_CI = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'lower_CI')
        df_u_CI = pd.read_excel(dir_read_corrs + df_name + '.xlsx',  sheet_name = 'upper_CI')

        df_cor.index = [df_cor.columns]
        df_p_val.index = [df_p_val.columns]
        df_l_CI.index = [df_l_CI.columns]
        df_u_CI.index = [df_u_CI.columns]

        # cor.loc['hrs_after_OP', 'hrs_after_OP']
        # params_list = df_u_CI.columns.to_list()
        params_list.remove(x)
        # params_list.remove('Unnamed: 0')

        fig, ax = plt.subplots(2, 3, 
                               figsize = (w_cm/2.54, h_cm/2.54),
                               constrained_layout = True)

        fig.set_constrained_layout_pads(
            w_pad = 0.03,
            h_pad = 0.03,
            wspace = 0.15,
            hspace = 0.15)
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
                    format='svg')
        plt.show()

def intr_params_VM_vs_FM(df_FM, df_VM, dict_params, dest_dir = False, w_cm = 18, h_cm = 11):
    """
    Plots parameters from the data of FM and VM
    boxplots to compare the data
    """
    titles_dict = pl_intr.dict_for_plotting_FM()
    colors = ['green', 'purple']

    fig2, axs = plt.subplots(4,4, figsize = (w_cm/2.54, h_cm / 2.54))
    axs = axs.flatten()

    for ax_el, (param_VM, param_FM) in enumerate(dict_params.items()):
        if param_VM == "('700pA', 'num_aps')":
            plot_spec = titles_dict['max_spikes']
        elif param_VM == 'rheo_complete':
            plot_spec = titles_dict['rheo_ramp_c']
        else:
            plot_spec = titles_dict[param_VM]

        if len(param_FM) == 0:
            continue

        for i, df in enumerate([df_FM, df_VM]):
            param = [param_FM, param_VM][i]
            col = colors[i]

            i = i+1
            df_plot = df[~df[param].isna()]  # remove NaN values
            x = np.linspace(i-0.6, i+0.3, len(df_plot))

            axs[ax_el].scatter(x, df_plot[param], alpha = 0.3,
                        color = col, zorder = 5) # dots inside
            axs[ax_el].scatter(x, df_plot[param], alpha = 0.4,
                        edgecolor = col, facecolors = 'none')

            axs[ax_el].boxplot(df_plot[param], positions = [i],
                            zorder = 10, widths = 0.5)


        axs[ax_el].set_xticks(ticks = [1, 2], labels = ['FM', 'VM'])

        axs[ax_el].set_ylabel(plot_spec[0] +'\n(' + plot_spec[1] + ')')
        axs[ax_el].set_yticks(ticks = plot_spec[2], labels = plot_spec[3])

    date = str(datetime.date.today())
    os.makedirs(dest_dir, exist_ok = True)
    # plt.savefig(f'{dest_dir}{date}comparison_plot.svg', \
    #     format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.savefig(f'{dest_dir}{date}comparison_plot.png', \
            dpi = 300, format = 'png')

    plt.show()

def t_test_and_plot_params(df_FM, df_VM, dict_params, save_dir = False,
                                 w_cm = 18, h_cm = 11):
    '''
    normalizes the parameters fromthe two dataframes and 
    performs a t-test to see if they are different/identical
    '''
    colors_dict = {'FM': ['purple'],
                   'VM': ['green']}

    fig2, axs = plt.subplots(4,4, figsize = (w_cm/2.54, h_cm / 2.54))
    axs = axs.flatten()

    for ax_el, (param_VM, param_FM) in enumerate(dict_params.items()):
        if len(param_FM) == 0:
            continue

        params = {'FM': df_FM[param_FM],
                  'VM': df_VM[param_VM]}

        # perform t-test
        t_stat, p_value = stats.ttest_ind(df_VM[param_VM], df_FM[param_FM], equal_var=False)
        if p_value < 0.05:
            p_value = p_value * 12
        if p_value < 0.05:
            print(f'significant differences in {param_VM}')
        for i, (key, param) in enumerate(params.items()):
            col = colors_dict[key]

            i = i+1
            x = np.linspace(i-0.6, i+0.3, len(param))

            axs[ax_el].scatter(x, param, alpha = 0.3,
                        color = col, zorder = 5) # dots inside
            axs[ax_el].scatter(x, param, alpha = 0.4,
                        edgecolor = col, facecolors = 'none')

            axs[ax_el].boxplot(param, positions = [i],
                            zorder = 10, widths = 0.5)

        axs[ax_el].set_title(f't = {t_stat.round(2)}, p = {p_value.round(4)}')
        axs[ax_el].set_xticks(ticks = [1, 2], labels = list(params.keys()))
        axs[ax_el].set_ylabel(f'norm {param_VM}')

    if save_dir:
        date = str(datetime.date.today())
        os.makedirs(save_dir, exist_ok = True)
        plt.savefig(f'{save_dir}{date}_t_test_comparison_plot.png', \
            format = 'png', bbox_inches = 'tight', dpi = 300)

    plt.show()

def plot_e_types_across_datasets(df_paper_VM_complete, save_dir=False):
    df_paper_VM_complete['e_type'] = df_paper_VM_complete['e_type'].astype(int)
    colors = sns.color_palette('colorblind')
    palette = {
        2: colors[1], # yellow
        3: colors[2], # green
        4: colors[0], # blue
        5: colors[3]  # orange
    }
    fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)

    sns.scatterplot(
        data=df_paper_VM_complete[df_paper_VM_complete['dataset'] == 'paper'],
        x='PC1',
        y='PC2',
        hue='e_type',
        palette=palette,
        marker='o',
        alpha=0.75,
        s=80,
        ax=ax[0]
    )

    sns.scatterplot(
        data=df_paper_VM_complete[df_paper_VM_complete['dataset'] == 'my data'],
        x='PC1',
        y='PC2',
        hue='e_type',
        palette=palette,
        marker='X',
        alpha=0.75,
        s=80,
        ax=ax[1]
    )

    ax[0].set_title('Paper')
    ax[1].set_title('My data')
    plt.tight_layout()
    if save_dir:
        plt.savefig(f'{save_dir}_PCA_across_dfs_plot.png', \
        dpi = 300, format = 'png')

    plt.show()

def plot_connection_ratio_bar(df_con_percentage_all, data_type='slice', save_dir=False):
    """
    Plot the ratio of found connections to all possible connections for each treatment-day combination
    as a filled bar up to the measured ratio.
    the numbers show pairs with connection (below) and 
    pairs where no connection was found(in bar)
    """
    w_cm = 9
    h_cm = 6
    fig, ax = plt.subplots(1, 1, figsize=(w_cm/2.54, h_cm/2.54))
    colors = figure_colors()[data_type]

    bar_width = 0.67
    ticks_, labels_ = [], []
    for j, day in enumerate(sorted(df_con_percentage_all.day.unique())):
        if 'TTX' in df_con_percentage_all.treatment.unique():
            treatments = ['Ctrl', 'high K', 'TTX']
        else:
            treatments = ['Ctrl', 'high K']
        for i, tr in enumerate(treatments):
            k_ = j + 2.5 * i
            df_plot = df_con_percentage_all[(df_con_percentage_all.treatment == tr) &
                                            (df_con_percentage_all.day == day)]

            if len(df_plot) == 0:
                continue

            total_found = df_plot.found_connections.sum()
            total_possible = df_plot.num_possible_connections.sum()
            overall_ratio = total_found / total_possible

            x_pos = k_ + 1

            # full bar to 100%
            ax.bar(x_pos, 1.0, width=bar_width, color='none',
                   edgecolor='black', linewidth=1.0, zorder=2)

            # filled portion up to the observed ratio
            ax.bar(x_pos, overall_ratio, width=bar_width,
                   color=colors[tr + day], alpha=0.85, zorder=1)

            # write counts inside the bar
            text_y = 0.06 # max(overall_ratio * 0.5, 0.04)
            text_color = 'white' if overall_ratio > 0.2 else 'black'
            ax.text(x_pos, text_y, f'{total_found}',
                    size=9, ha='center', va='center', color=text_color, zorder=3)
            ax.text(x_pos, 0.55, f'{total_possible-total_found}',
                                size=9, ha='center', va='center', color='k', zorder=3)
            labels_.append(f'RS {day[-1]} \n {tr}')
            ticks_.append(x_pos)

    ax.set_xticks(ticks = ticks_, labels = labels_)
    ax.set_ylabel('Connection ratio\n(found / tested)')
    ax.set_ylim(0, 1.15)

    if save_dir:
        plt.savefig(f'{save_dir}{data_type}_connection_ratio_bar.svg',
                    format='svg', bbox_inches='tight', dpi=300)
    plt.show()

def plot_connectivity_percentages(df_connections, save_dir = False):

    w_cm = 10
    h_cm = 6
    fig, ax = plt.subplots(1, 1, figsize=(w_cm/2.54, h_cm/2.54))
    cols_num_cells = plt.get_cmap('YlGnBu')(np.linspace(0.25, 0.9, 7))

    ticks_, labels = [], []
    for j, day in enumerate(['D1', 'D2']):
        if 'TTX' in df_connections.treatment.unique():
            treatments = ['Ctrl', 'high K', 'TTX']
        else:
            treatments = ['Ctrl', 'high K']
        for i, tr in enumerate(treatments):
            k_ = j + 2.5 * i
            df_plot = df_connections[(df_connections.treatment == tr)\
                                         & (df_connections.day == day)].reset_index(drop = True)
            avg = np.nanmean(df_plot.con_percentage_fixed)
            #print(median)
            x = np.linspace(k_ - 0.3, k_ + 0.3, len(df_plot))
            for i, con_percent in enumerate(df_plot.con_percentage_fixed):
                num_cells = df_plot.loc[i, 'num_channels']
                ax.scatter(x[i], con_percent, c = cols_num_cells[num_cells-2], alpha = 0.85)
                ax.scatter(x[i], con_percent, facecolors='none',
                           edgecolors = cols_num_cells[num_cells-2], zorder = 8,
                           linewidths=2.6)
            # ax.boxplot(df_plot.con_percentage_fixed, positions = [k_], zorder = 10,
            #           widths = 0.4)
            ax.plot([k_-0.3, k_+0.3], [avg, avg], zorder = 10, color = '#430704',
                    linewidth = 1.7)

            labels.append(f'RS{day[-1]} \n{tr}')
            ticks_.append(k_)
            # ax.text(k_-0.25, 1.2*avg, str(round(avg,2)), size = 10)


    ax.set_xticks(ticks = ticks_, labels = labels)
    ax.set_yticks(ticks = np.linspace(0, 100, 5))
    ax.set_ylabel('Connectivity (%)')
    legend_handles = [Line2D([0], [0], color=cols_num_cells[i], lw=2, label=str(i))
                      for i in range(len(cols_num_cells))]
    plt.legend(handles=legend_handles, title='#cells \n per slice',
               loc='upper left', bbox_to_anchor=(1.02, 1))

    if save_dir:
        plt.savefig(f'{save_dir}connectivity_percentage.svg',
            format='svg', bbox_inches='tight', dpi=300)
    plt.show()

def plot_connection_matrix(connection_matrix, title=None, ax=None,
                           cmap='Blues', vmin=0, vmax=20,
                           colorbar=True):
    '''
    Plot connection probabilities as a connectivity matrix.

    Cell labels show:
    found/tested
    (connection probability)
    '''

    if ax is None:
        fig, ax = plt.subplots(figsize=(5/2.54, 5/2.54))

    probabilities = np.full(connection_matrix.shape, np.nan)
    labels = np.empty(connection_matrix.shape, dtype=object)

    for i in range(connection_matrix.shape[0]):
        for j in range(connection_matrix.shape[1]):

            value = str(connection_matrix.iloc[i, j])

            match = re.match(r'([\d.]+)%\s*\((\d+)/(\d+)\)', value)

            if match:
                probability = float(match.group(1))
                found = match.group(2)
                tested = match.group(3)

                probabilities[i, j] = probability
                labels[i, j] = (f'{found}/{tested}\n'f'({probability:g}%)')
            else:
                labels[i, j] = '0/0\n(NA)'

    if vmax is None:
        vmax = np.nanmax(probabilities)

    cmap_plot = plt.get_cmap(cmap).copy()
    cmap_plot.set_bad('white')

    im = ax.imshow(probabilities, cmap=cmap_plot, vmin=vmin, vmax=vmax)

    # cell borders
    ax.set_xticks(np.arange(-0.5, connection_matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, connection_matrix.shape[0], 1), minor=True)

    ax.grid(which='minor', linewidth=1.5)
    ax.tick_params(which='minor', bottom=False, left=False)

    # axis labels
    ax.set_xticks(range(connection_matrix.shape[1]))
    ax.set_yticks(range(connection_matrix.shape[0]))

    ax.set_xticklabels(connection_matrix.columns, rotation = 45)
    ax.set_yticklabels(connection_matrix.index)

    ax.set_xlabel('Postsynaptic e-type')
    ax.set_ylabel('Presynaptic e-type')

    if title is not None:
        ax.set_title(title)

    # annotations
    for i in range(connection_matrix.shape[0]):
        for j in range(connection_matrix.shape[1]):

            probability = probabilities[i, j]

            if np.isnan(probability):
                text_color = 'black'
            else:
                relative_value = (probability - vmin) / (vmax - vmin)
                text_color = 'white' if relative_value > 0.55 else 'black'

            ax.text(j, i, labels[i, j], ha='center', va='center', color=text_color)

    if colorbar:
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Connection probability (%)')

    return ax

def plot_param_matrix(param_matrix, n_matrix, parameter, title=None, ax=None,
                      cmap=None, vmin=-60, vmax=60, colorbar=True,
                      baseline_matrix=None):
    '''
    Plot mean EPSP amplitude (mV) or latency (ms) for each
    pre/post e-type combination, colored by percentage change from RS1.
    '''

    if ax is None:
        fig, ax = plt.subplots(figsize=(5/2.54, 5/2.54))

    param_values = param_matrix.to_numpy(dtype=float)

    if baseline_matrix is None:
        values = np.full_like(param_values, np.nan)  # RS1: white
    else:
        baseline = baseline_matrix.reindex(index=param_matrix.index, columns=param_matrix.columns).to_numpy(dtype=float)
        values = np.full_like(param_values, np.nan)
        np.divide(100 * (param_values - baseline), baseline, out=values,
                where=(baseline > 0) & np.isfinite(param_values) & (n_matrix.to_numpy() > 0))

    # 10% color intervals
    bounds = np.arange(vmin, vmax + 10, 10)

    if cmap is None:
        cmap_plot = LinearSegmentedColormap.from_list(
            'orange_white_green',
            ['#d95f02', '#f7f7f7', '#1b9e77'],
            N=len(bounds) - 1
        )
    else:
        cmap_plot = plt.get_cmap(cmap, len(bounds) - 1).copy()

    cmap_plot.set_bad('white')

    im = ax.imshow(values, cmap=cmap_plot,
                norm=BoundaryNorm(bounds, cmap_plot.N))


    # cell borders
    ax.set_xticks(np.arange(-0.5, param_matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, param_matrix.shape[0], 1), minor=True)
    ax.grid(which='minor', linewidth=1.5)
    ax.tick_params(which='minor', bottom=False, left=False)

    # axis labels
    ax.set_xticks(range(param_matrix.shape[1]))
    ax.set_yticks(range(param_matrix.shape[0]))
    ax.set_xticklabels(param_matrix.columns, rotation=45)
    ax.set_yticklabels(param_matrix.index)

    ax.set_xlabel('Postsynaptic e-type')
    ax.set_ylabel('Presynaptic e-type')

    if title is not None:
        ax.set_title(title)

    # annotations
    for i in range(param_matrix.shape[0]):
        for j in range(param_matrix.shape[1]):

            param_val = param_values[i, j]
            n = int(n_matrix.iloc[i, j])

            if n == 0 or not np.isfinite(param_val):
                label = f'NA\n(n={n})'
            elif baseline_matrix is None:
                label = f'{param_val:.2f}\n(n={n})'
            else:
                # change = values[i, j]
                label = f'{param_val:.2f}\n(n={n})'
                # label = (f'{amp:.2f}\n({change:+.1f}%, n={n})'
                #         if np.isfinite(change) else f'{amp:.2f}\n(change NA, n={n})')

            text_color = 'white' if np.isfinite(values[i, j]) and abs(values[i, j]) >= 70 else 'black'
            ax.text(j, i, label, ha='center', va='center', color=text_color)

    if colorbar and baseline_matrix is not None:
        cbar = plt.colorbar(im, ax=ax, ticks=np.arange(vmin, vmax + 1, 20))
        cbar.set_label(f'Change in {parameter} from RS1 (%)')

    return ax

def plot_etype_proportions(df, df_FM, save_dir=False):
    '''
    Plot e-type proportions for each day x treatment condition
    as 100% stacked bars.
    '''
    fig, ax = plt.subplots(1, 1, figsize=(9/2.54, 6/2.54))

    etype_names = df_FM[['etype_idx', 'etype']].drop_duplicates().set_index('etype_idx')['etype'].to_dict()
    etype_order = ['LowRin', 'SlowAP', 'FastAP', 'DoubSpk']
    ETYPES = [k for name in etype_order for k, v in etype_names.items() if v == name]

    proportions = (
        df.groupby(['day', 'treatment'])['e_type']
        .value_counts(normalize=True)
        .mul(100)
        .unstack(fill_value=0)
        .reindex(columns=ETYPES, fill_value=0)
    )

    day_order = ['D1', 'D2']
    if 'TTX' in df.treatment.unique():
        treatment_order = ['Ctrl', 'high K', 'TTX']
    else:
        treatment_order = ['Ctrl', 'high K']

    conditions = sorted(
        proportions.index.tolist(),
        key=lambda x: (treatment_order.index(x[1]), day_order.index(x[0]))
    )
    labels = [f'{day} {treatment}' for day, treatment in conditions]

    etype_colors = {
        2: '#4E79A7',
        3: '#F28E2B',
        4: '#59A14F',
        5: '#E15759'
    }

    x = np.arange(len(conditions))
    bottom = np.zeros(len(conditions))

    for etype in ETYPES:
        values = proportions[etype].values
        ax.bar(x, values, bottom=bottom, width=0.67, color=etype_colors[etype],
               label=etype_names[etype])
        bottom += values

    ax.set_xticks(x, labels=labels, rotation = 45)
    ax.set_ylabel('E-type proportion (%)')
    ax.set_ylim(0, 100)
    ax.legend(title='E-type', loc='upper left', bbox_to_anchor=(1.02, 1))

    if save_dir:
        plt.savefig(f'{save_dir}etype_proportions.svg', format='svg', bbox_inches='tight', dpi=300)

    plt.show()

    return proportions


def plot_connectivity_change_matrix(connection_matrix, baseline_matrix=None,
                                    title=None, ax=None, vmin=-15, vmax=15,
                                    step=5, colorbar=True):
    """Connectivity matrix colored by percentage-point change from baseline."""

    from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm

    def parse_matrix(matrix):
        parsed = matrix.stack(dropna=False).astype(str).str.extract(
            r'([\d.]+)%\s*\((\d+)/(\d+)\)')
        found = pd.to_numeric(parsed[1]).unstack().reindex_like(matrix)
        tested = pd.to_numeric(parsed[2]).unstack().reindex_like(matrix)
        prob = 100 * found / tested.where(tested > 0)
        return prob, found, tested

    if ax is None:
        fig, ax = plt.subplots(figsize=(5/2.54, 5/2.54))

    prob, found, tested = parse_matrix(connection_matrix)

    if baseline_matrix is None:
        changes = np.full(prob.shape, np.nan)
    else:
        baseline, _, _ = parse_matrix(baseline_matrix)
        baseline = baseline.reindex(index=prob.index, columns=prob.columns)
        changes = (prob - baseline).to_numpy(dtype=float)

    bounds = np.arange(vmin, vmax + step, step)
    cmap = LinearSegmentedColormap.from_list(
        'orange_white_green', ['#d95f02', '#f7f7f7', '#1b9e77'],
        N=len(bounds)-1
    )
    cmap.set_bad('white')
    norm = BoundaryNorm(bounds, cmap.N, clip=True)

    im = ax.imshow(changes, cmap=cmap, norm=norm)

    # Cell borders
    ax.set_xticks(np.arange(-0.5, prob.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, prob.shape[0], 1), minor=True)
    ax.grid(which='minor', color='0.7' if baseline_matrix is None else 'white',
            linewidth=1.5)
    ax.tick_params(which='minor', bottom=False, left=False)

    # Labels
    ax.set_xticks(range(prob.shape[1]), prob.columns, rotation=45)
    ax.set_yticks(range(prob.shape[0]), prob.index)
    ax.set_xlabel('Postsynaptic e-type')
    ax.set_ylabel('Presynaptic e-type')
    ax.set_title(title or '')

    for i in range(prob.shape[0]):
        for j in range(prob.shape[1]):
            n = tested.iloc[i, j]

            if pd.isna(n) or n == 0:
                label = 'NA\n(n=0)'
            else:
                label = f'{int(found.iloc[i, j])}/{int(n)}\n({prob.iloc[i, j]:.1f}%)'

            change = changes[i, j]
            text_color = 'black'
            if np.isfinite(change):
                r, g, b, _ = cmap(norm(change))
                if 0.2126*r + 0.7152*g + 0.0722*b < 0.48:
                    text_color = 'white'

            ax.text(j, i, label, ha='center', va='center', color=text_color)

    if colorbar and baseline_matrix is not None:
        cbar = plt.colorbar(im, ax=ax, ticks=bounds[::2])
        cbar.set_label('Change in connectivity (percentage points)')

    return ax



def plot_connect_params(df, data_type, DV='amp',save_dir=False):

    h_cm = 4.5
    w_cm = 20
    titles_dict = pl_intr.dict_for_plotting_conn()
    df = df.sort_values(['patcher', 'OP', 'day', 'connection_ID'])
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

def plot_norm_con_params_log_ci(df, DV, data_type, save_dir=False):
    '''
    df - summary dataframe exported from R (EPSP_log_95CI_amp/lat.csv)
    DV - 'amp' or 'lat'
    data_type - 'slice' or 'cell_ID_new'
    '''
    titles_dict = pl_intr.dict_for_plotting_conn()
    w_cm, h_cm = 7, 4.5
    colors = figure_colors()[data_type]

    fig, ax = plt.subplots(1, 1, figsize=(w_cm/2.54, h_cm/2.54))
    x_positions = np.arange(1, 5)

    # Map the R conditions to your original colour keys
    if 'TTX' in df.condition.unique():
        conditions = {'RS0': ('Ctrl', 'D1'),
                        'Ctrl': ('Ctrl', 'D2'),
                        'high K': ('high K', 'D2'),
                        'TTX': ('TTX', 'D2')}
    else:
        conditions = {'RS0': ('Ctrl', 'D1'),
                    'Ctrl': ('Ctrl', 'D2'),
                    'high K': ('high K', 'D2')}

    for i, (condition, (tr, day)) in enumerate(conditions.items()):

        df_plot = df[df.condition == condition].sort_values('EPSP_number')

        if len(df_plot) == 0:
            continue

        # Mean log-transformed values and 95% CI from R
        #x = df_plot.EPSP_number.to_numpy()
        x = df_plot.EPSP_number.to_numpy() + (i - (len(conditions) - 1) / 2) * 0.12
        normalized_means = df_plot.mean_log.to_numpy()
        normalized_errors = np.array([
            normalized_means - df_plot.lower.to_numpy(),
            df_plot.upper.to_numpy() - normalized_means
        ])

        # Plot line with error bars
        color_key = f'{tr}{day}' if f'{tr}{day}' in colors else tr

        ax.errorbar(x, normalized_means, yerr=normalized_errors,
                    color=colors[color_key], marker='o', linewidth=1,
                    markersize=4, capsize=2.2, capthick=1, alpha = 0.6,
                    label=f'{tr} {day}' if day == 'D2' else 'RS0')

    ax.set_xticks(ticks=x_positions, labels=[f'EPSP {x}' for x in x_positions])
    ax.set_ylabel(f'Log({titles_dict[DV][0]})')

    # Reference line at log(1) = 0
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.7)

    # Save if directory provided
    if save_dir:
        plt.savefig(f'{save_dir}_normalized_{DV}_{data_type}_95CI.svg',
                    format='svg', bbox_inches='tight', dpi=300)

    plt.show()

def plot_con_facil(df, data_type, save_dir=False):
    '''
    df - facilitation counts dataframe exported from R
    data_type - 'slice' or 'connection_ID'
    '''
    w_cm, h_cm = 6, 4.5
    colors = figure_colors()[data_type]

    fig, ax = plt.subplots(1, 1, figsize=(w_cm/2.54, h_cm/2.54))

    conditions = {'RS0': ('Ctrl', 'D1'),
                  'Ctrl': ('Ctrl', 'D2'),
                  'high K': ('high K', 'D2')}

    for i, (condition, (tr, day)) in enumerate(conditions.items()):

        df_plot = df[df.condition == condition]

        if len(df_plot) == 0:
            continue

        mean = df_plot.percent.iloc[0]
        lower = df_plot.lower.iloc[0]
        upper = df_plot.upper.iloc[0]

        errors = [[mean - lower], [upper - mean]]
        color_key = f'{tr}{day}' if f'{tr}{day}' in colors else tr

        ax.errorbar(i, mean, yerr=errors,
                    color=colors[color_key], marker='o',
                    linewidth=1, markersize=4,
                    capsize=2.2, capthick=1)

        print(f'{condition}: {df_plot.n_facilitating.iloc[0]}/{df_plot.n.iloc[0]}')

    ax.set_xticks(ticks=np.arange(3), labels=list(conditions))
    ax.set_ylabel('Facilitating connections (%)')
    ax.set_ylim(0, 100)
    ax.set_yticks(np.arange(0, 101, 25))

    if save_dir:
        plt.savefig(f'{save_dir}_facilitation_{data_type}.svg',
                    format='svg', bbox_inches='tight', dpi=300)

    plt.show()

def plot_con_ppr_scatter_ci(df, df_ci, data_type, save_dir=False):
    '''
    df - original connection dataframe
    df_ci - summary dataframe exported from R (with mean_log, lower, upper)
    data_type - 'slice' or 'connection_ID'
    '''
    w_cm, h_cm = 7, 5
    colors = figure_colors()[data_type]
    stat_col = '#430704'

    fig, ax = plt.subplots(1, 1, figsize=(w_cm/2.54, h_cm/2.54))

    if 'TTX' in df_ci.condition.unique():
        conditions = {'RS0': ('Ctrl', 'D1'),
                        'Ctrl': ('Ctrl', 'D2'),
                        'high K': ('high K', 'D2'),
                        'TTX': ('TTX', 'D2')}
    else:
        conditions = {'RS0': ('Ctrl', 'D1'),
                    'Ctrl': ('Ctrl', 'D2'),
                    'high K': ('high K', 'D2')}

    for i, (condition, (tr, day)) in enumerate(conditions.items()):

        if day == 'D1':
            df_plot = df[df.day == day].copy()
        else:
            df_plot = df[(df.day == day) & (df.treatment == tr)].copy()

        df_plot = df_plot.dropna(subset=['Amp 1', 'Amp 2'])
        df_plot = df_plot[(df_plot[['Amp 1', 'Amp 2']] > 0).all(axis=1)]

        if len(df_plot) == 0:
            continue

        # raw individual PPR values
        values = (df_plot['Amp 2'] / df_plot['Amp 1']).to_numpy()

        # jitter the points slightly to the left of the group centre
        x_vals = np.linspace(0.7 + i, 1.3 + i, len(df_plot))

        color_key = f'{tr}{day}' if f'{tr}{day}' in colors else tr

        ax.scatter(x_vals, values,
                   color=colors[color_key],
                   alpha=0.7, s=12, edgecolors='none', zorder=2)

        # summary stats from R, back-transformed from log scale
        ci = df_ci[(df_ci.condition == condition) &
                   (df_ci.EPSP_number == 2)].iloc[0]

        mean = np.exp(float(ci['mean_log']))
        lower = np.exp(float(ci['lower']))
        upper = np.exp(float(ci['upper']))

        # manually draw the CI interval
        x_stat = i + 1
        cap_w = 0.1
        mean_w = 0.1

        # vertical CI line
        ax.vlines(x_stat, lower, upper, color=stat_col, linewidth=1.2, zorder=4)

        # CI caps
        ax.hlines([lower, upper], x_stat - cap_w, x_stat + cap_w,
                  color=stat_col, linewidth=1.2, zorder=4)

        # mean line
        ax.hlines(mean, x_stat - mean_w, x_stat + mean_w,
                  color=stat_col, linewidth=1.2, zorder=5)

        print(f'{condition}: n={len(values)}')

    ax.set_xticks(ticks=np.arange(len(conditions)+1), labels=[''] + list(conditions))
    ax.set_ylabel('PPR \n (EPSP2 / EPSP1)')

    # keep PPR axis with facilitation/depression reference at 1
    # ax.set_yscale('log')
    # ax.set_yticks([0.3, 0.5, 1, 2, 3])
    # ax.set_yticklabels(['0.3', '0.5', '1', '2', '3'])
    # ax.minorticks_off()

    ax.set_ylim(bottom=0)
    ax.set_yticks(ticks = np.linspace(0, 3, 7))
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.7)

    if save_dir:
        plt.savefig(f'{save_dir}_PPR_scatter_{data_type}_95CI.svg',
                    format='svg', bbox_inches='tight', dpi=300)

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

def plot_amp1_vs_failure_rate(df, df_fit, data_type, w_cm=18, h_cm=5.5,
                              model_col='failure_prop', save_dir=False):
    '''
    EPSP1 amplitude vs failure probability.
    One panel per treatment, comparing D1 and D2.
    Fits are taken from the beta-binomial GLMM in R.
    '''

    colors = figure_colors()[data_type]
    treatments = ['Ctrl', 'high K', 'TTX'] if 'TTX' in df.treatment.unique() else ['Ctrl', 'high K']

    fig, ax = plt.subplots(1, len(treatments), figsize=(w_cm/2.54, h_cm/2.54),
                           sharex=True, sharey=True)
    ax = np.atleast_1d(ax)

    for i, tr in enumerate(treatments):
        for day in ['D1', 'D2']:
            color = (colors.get('pre', 'darkgrey') if day == 'D1'
                     else colors.get(f'{tr}D2', colors.get(tr, 'grey')))

            df_plot = df[(df.treatment == tr) & (df.day == day)]
            # df_plot = df_plot[(df_plot['Amp 1'] > 0) & df_plot['failure_prop'].notna()]

            fitted = df_fit[(df_fit.treatment == tr) & (df_fit.day == day)]
            fitted = fitted.sort_values('Amp1_mV')

            # Individual connections
            ax[i].scatter(np.log(df_plot['Amp 1']), df_plot['failure_rate'] * 100,
                          color=color, alpha=0.45, s=16, zorder=2)

            # Beta-binomial GLMM predictions from R
            ax[i].plot(np.log(fitted['Amp1_mV']), fitted[model_col] * 100,
                       color=color, linewidth=1.7, label=day, zorder=3)

        ax[i].set_title(tr)
        ax[i].set_xlabel('Log(EPSP1 amplitude [mV])')
        ax[i].set_ylim(-5, 105)
        ax[i].legend(frameon=False, fontsize=8)

    ax[0].set_ylabel('Failure rate (%)')
    plt.tight_layout()

    if save_dir:
        fig.savefig(f'{save_dir}amp1_failure_{data_type}.svg',
                    format='svg', bbox_inches='tight')

    plt.show()
    return fig, ax

def plot_ppr_comparison(df_descriptive, df_model, data_type, DV='amp',
                        save_dir=False, w_cm=15, h_cm=5.7):
    """Side-by-side observed donor-bootstrap vs LMM means and 95% CIs (log scale)."""
    palette = figure_colors()['slice' if data_type in ('slice', 'slice_ID') else 'cell_ID_new']
    colors = {'RS0': palette['CtrlD1'], 'Ctrl': palette['CtrlD2'],
              'high K': palette['high KD2'], 'TTX': palette.get('TTXD2', '#0e62ff')}
    conditions = [c for c in ('RS0', 'Ctrl', 'high K', 'TTX')
                  if c in set(df_descriptive.condition) | set(df_model.condition)]
    offsets = dict(zip(conditions, np.linspace(-0.12, 0.12, len(conditions))))
    fig, axs = plt.subplots(1, 2, sharey=True, figsize=(w_cm/2.54, h_cm/2.54))

    for ax, (title, data) in zip(axs, [('Observed + donor bootstrap', df_descriptive),
                                       ('LMM-estimated', df_model)]):
        for condition in conditions:
            d = data.loc[data.condition.eq(condition)].sort_values('EPSP_number')
            if d.empty:
                continue
            x = d.EPSP_number.to_numpy(dtype=float) + offsets[condition]
            y = d.mean_log.to_numpy(dtype=float)
            lower, upper = d.lower.to_numpy(dtype=float), d.upper.to_numpy(dtype=float)
            if np.any((y < lower - 1e-8) | (y > upper + 1e-8)):
                raise ValueError(f'CI does not contain its mean: {title}, {condition}')
            ax.errorbar(x, y, yerr=np.vstack((y - lower, upper - y)),
                        color=colors[condition], marker='o', linewidth=1,
                        markersize=4, capsize=2.2, capthick=1, alpha=0.75,
                        label=condition)
        ax.axhline(0, color='grey', linestyle='--', alpha=0.7)
        ax.set_xticks(range(1, 5), [f'EPSP {i}' for i in range(1, 5)])
        ax.set_xlim(0.65, 4.35)
        ax.set_title(title)
        ax.set_xlabel('EPSP number')

    axs[0].set_ylabel('Log(normalized amplitude)' if DV == 'amp'
                      else 'Log(normalized latency)')
    axs[1].legend(frameon=False, loc='upper left', bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    if save_dir:
        fig.savefig(f'{save_dir}_PPR_comparison_{data_type}_{DV}.svg',
                    format='svg', bbox_inches='tight', dpi=300)
    return fig, axs

# #ANOHTE UGLY FUNTION

def plot_ppr_etypes(df, DV='amp', save_dir=False, w_cm=18, h_cm=5.7):
    """Mean log-normalized EPSPs grouped by postsynaptic e-type."""

    df = df.copy()
    cols = [f'Amp {i}' for i in range(1, 5)] if DV == 'amp' else [f'Lat{i}' for i in range(1, 5)]

    # Normalize to EPSP1 and log-transform
    baseline = df[cols[0]].where(df[cols[0]] > 0)
    for col in cols:
        df[col] = np.log(df[col].where(df[col] > 0) / baseline)

    df['condition'] = np.where(df.day == 'D1', 'RS0', df.treatment)

    df_long = df.melt(id_vars=['condition', 'e_type'], value_vars=cols,
                      var_name='EPSP', value_name='norm_response')
    df_long['EPSP_number'] = df_long.EPSP.str.extract(r'(\d+)$').astype(int)

    means = (df_long.groupby(['condition', 'e_type', 'EPSP_number'], as_index=False)
             ['norm_response'].mean())

    conditions = [c for c in ['RS0', 'Ctrl', 'high K', 'TTX'] if c in means.condition.unique()]
    etypes = sorted(means.e_type.dropna().unique())
    colors = dict(zip(etypes, plt.get_cmap('tab10').colors))

    fig, axs = plt.subplots(1, len(conditions), figsize=(w_cm/2.54, h_cm/2.54),
                            sharex=True, sharey=True)
    axs = np.atleast_1d(axs)

    for ax, condition in zip(axs, conditions):
        for etype in etypes:
            d = means[(means.condition == condition) & (means.e_type == etype)].sort_values('EPSP_number')
            if d.empty:
                continue

            ax.plot(d.EPSP_number, d.norm_response, '-o', color=colors[etype],
                    markersize=4, linewidth=1, alpha=0.8, label=str(etype))

        ax.axhline(0, color='grey', linestyle='--', alpha=0.6)
        ax.set_title(condition)
        ax.set_xticks(range(1, 5), [f'EPSP {i}' for i in range(1, 5)])
        ax.set_xlabel('EPSP number')

    axs[0].set_ylabel(f'Log(normalized {"amplitude" if DV == "amp" else "latency"})')
    axs[-1].legend(title='Postsynaptic e-type', frameon=False, bbox_to_anchor=(1.01, 1))
    fig.tight_layout()

    if save_dir:
        fig.savefig(f'{save_dir}PPR_etypes_{DV}.svg', format='svg', bbox_inches='tight')

    return fig, axs


def plot_con_ppr_scatter_etypes(df, save_dir=False, w_cm=11, h_cm=5.5):
    """Individual EPSP2/EPSP1 ratios colored by postsynaptic e-type."""

    df = df.copy()
    df['condition'] = np.where(df.day == 'D1', 'RS0', df.treatment)
    df = df.dropna(subset=['Amp 1', 'Amp 2', 'e_type_post'])
    df = df[(df[['Amp 1', 'Amp 2']] > 0).all(axis=1)]

    conditions = [c for c in ['RS0', 'Ctrl', 'high K', 'TTX'] if c in df.condition.unique()]
    etypes = sorted(df.e_type_post.unique())
    colors = dict(zip(etypes, plt.get_cmap('tab10').colors))

    fig, ax = plt.subplots(figsize=(w_cm/2.54, h_cm/2.54))

    for i, condition in enumerate(conditions):
        d = df[df.condition == condition]
        x = np.linspace(i + 0.7, i + 1.3, len(d))
        ppr = d['Amp 2'] / d['Amp 1']

        for etype in etypes:
            mask = d.e_type_post == etype
            ax.scatter(x[mask], ppr[mask], color=colors[etype],
                       alpha=0.7, s=15, label=f'e-type {etype}' if i == 0 else None)

    ax.axhline(1, color='grey', linestyle='--', alpha=0.7)
    ax.set_xticks(range(1, len(conditions)+1), conditions)
    ax.set_ylabel('PPR (EPSP2 / EPSP1)')
    ax.set_ylim(bottom=0)
    ax.legend(title='Postsynaptic e-type', frameon=False, bbox_to_anchor=(1.01, 1))

    fig.tight_layout()

    if save_dir:
        fig.savefig(f'{save_dir}_PPR_scatter_etypes.svg',
                    format='svg', bbox_inches='tight', dpi=300)

    return fig, ax


def plot_ppr_donor_heatmap(df_long, data_type, save_dir=False, EPSP_number=2,
                           w_cm=11, h_cm=None):
    """Donor × condition connection-count heatmap, calculated directly in pandas."""
    counts = get_ppr_donor_counts(df_long, EPSP_number=EPSP_number)
    conds = [c for c in ('RS0', 'TTX', 'high K', 'Ctrl') if c in counts.columns]
    matrix = counts.set_index('OP')[conds].sort_index(ascending=False)
    if matrix.empty:
        raise ValueError('No valid PPR observations found for the selected EPSP.')

    height = h_cm if h_cm is not None else max(6, 2 + 0.38 * len(matrix))
    fig, ax = plt.subplots(figsize=(w_cm/2.54, height/2.54))
    im = ax.imshow(matrix.to_numpy(), cmap='Blues', vmin=0,
                   aspect='auto', interpolation='nearest')
    ax.set_xticks(np.arange(len(conds)), conds)
    ax.set_yticks(np.arange(len(matrix)), matrix.index, fontsize=7)
    ax.set_ylabel('Donor')
    ax.set_title(f'Donor contributions to EPSP{EPSP_number} PPR')

    for row in range(len(matrix)):
        for col in range(len(conds)):
            n = int(matrix.iloc[row, col])
            if n:
                ax.text(col, row, str(n), ha='center', va='center', fontsize=7,
                        color='white' if n >= 0.65 * matrix.to_numpy().max() else 'black')

    fig.colorbar(im, ax=ax, label='Connections', shrink=0.65)
    ax.tick_params(length=0)

    if save_dir:
        fig.savefig(f'{save_dir}sup_PPR_donor_heatmap_{data_type}.svg',
                    format='svg', bbox_inches='tight', dpi=300)
    plt.show()

def get_ppr_donor_counts(df, EPSP_number=2):
    """Recompute each donor's connection counts from the individual long-format R export.

    Each valid row at this EPSP is counted once, exactly as R's `n()` does.
    """
    required = {'OP', 'condition', 'EPSP_number', 'norm_response'}
    if missing := required - set(df.columns):
        raise ValueError(f'Missing columns: {sorted(missing)}')

    d = df.loc[df.EPSP_number.eq(EPSP_number)].copy()
    d['norm_response'] = pd.to_numeric(d['norm_response'], errors='coerce')
    d = d.loc[np.isfinite(d.norm_response)].dropna(subset=['OP', 'condition'])
    d['OP'] = d.OP.astype(str)
    d['condition'] = d.condition.astype(str).str.strip()

    order = [c for c in ('RS0', 'TTX', 'high K', 'Ctrl') if c in d.condition.unique()]
    counts = (d.groupby(['OP', 'condition'], observed=True).size()
                .unstack(fill_value=0).reindex(columns=order, fill_value=0)
                .sort_index().astype(int))
    counts['total_connections'] = counts.sum(axis=1)
    return counts.reset_index()

# For traces?
def lowpass(trace, cutoff=1000, fs=20000):
    sos = butter(4, cutoff, btype='lowpass', fs=fs, output='sos')
    return sosfiltfilt(sos, trace)

def plot_connect_D1_D2(OP, patcher, file_index, active_channels, save_dir = False,
        human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):

    
    w_cm = len(active_channels[0]) * 3.2
    h_cm = len(active_channels[0]) * 1.5
    z1 = 1.5
    z2 = 40.5
    stim_window = stim_win.stim_window_con_screen

    linewidth_ = 0.8
    clrs = [['#332288',  # purple
        '#88CCEE',  # light blue
        '#44AA99',  # teal
        '#117733',  # dark green
        '#DDCC77',  # ochre
        '#CC6677',  # muted red
        '#882255',  # burgundy
        '#AA4499'],  # magenta
        ['#6699CC',  # medium blue
        '#999933',  # olive
        '#661100',  # dark brown
        '#E69F00',  # orange
        '#56B4E9',  # sky blue
        '#0072B2',  # deep blue
        '#D55E00']]   # vermillion


    #y_labels_D2 = [f're Ch{ch}' for ch in active_channels[0]]

    work_dir, filenames, _, _, _, _ = sort.get_OP_metadata(human_dir, OP, patcher)
    fn = work_dir + filenames[file_index[0]]
    fn2 = work_dir + filenames[file_index[1]]
    con_screen_data = [hcf.load_traces(fn), hcf.load_traces(fn2)]

    fig, ax = plt.subplots(len(active_channels[0]), 2 * len(active_channels[0]), sharex = True,
                           sharey = False, figsize = (w_cm/2.54, h_cm / 2.54),
                           gridspec_kw={'wspace': 0.02, 'hspace': 0.02})
    for f, data in enumerate(con_screen_data):
        for i, ch1 in enumerate(active_channels[f]):
            ch1_data = data[f'Ch{ch1}'][0]
            avg = np.mean(ch1_data, axis = 1)

            for j, ch2 in enumerate(active_channels[0]):
                if f==1:
                    j = j +len(active_channels[0])
                plotwin = avg[stim_window[f'Ch{ch2}'][0] -600 :stim_window[f'Ch{ch2}'][1]-1000] #from the average takes the signal for this stim_window
                plotwin = lowpass(plotwin, cutoff=1000)
                ax[i,j].plot(plotwin, c = clrs[f][i], lw = linewidth_)
                
                # if j == 0:
                #     ax[i, j].set_ylabel(f'Ch{ch1}')
                # if j == len(active_channels[f]):
                #     ax[i, j].set_ylabel(y_labels_D2[i])
                ax[i, j].yaxis.label.set_color(clrs[f][i])
                ax[i, j].set_xticks([])
                ax[i, j].set_yticks([])
                ax[i, j].spines['bottom'].set_visible(False)
                ax[i, j].spines['left'].set_visible(False)

                if plotwin.max()-plotwin.min() < 10:
                    ax[i,j].set_ylim([plotwin.min() - z1, plotwin.max() + z1])
                    # v1 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z1, color='k') # lw=0.7
                else:
                    ax[i,j].set_ylim([plotwin.min() - z1, plotwin.max() + z2])
                    # v2 = ax[i,j].vlines(0,plotwin.min() + z1, plotwin.max() + z2, color='k')

    if save_dir:
        plt.savefig(f'{save_dir}connect_trace_{OP}.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 500)
    plt.show()
from glob import glob
import os
from collections import Counter
# import importlib
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ephys_analysis.funcs_sorting as sort
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
from ephys_analysis.funcs_human_characterisation import cap_values_adjustment
from ephys_analysis.funcs_for_results_tables import collect_intrinsic_df, collect_IFF_dfs

#%% 
def main():
    '''
    plots intrinsic and time dependency data
    '''
    plot_time_dep()
    # plot_intrinsic_props()

def plot_time_dep():
    '''
    plot time dependency for 5 parameteres
    '''
    # plots over time
    param_list = ['resting_potential', 'Rin', 'cap_adj', 'TH', 'membra_time_constant_tau']

    adult_df_ctrl = load_intr_data('Ctrl')
    df_qc_hold = get_lab_book_info(adult_df_ctrl)
    df_qc_hold = get_match_cell_id(df_qc_hold)
    adult_df_ctrl = get_match_cell_id(adult_df_ctrl)

    # join QC and adult_df based on unique cell_ID
    adult_df_ctrl = pd.merge(adult_df_ctrl, df_qc_hold, on = 'cell_IDs_match', \
                        suffixes = (None, '_y')).reset_index(drop=True)

    adult_df_ctrl['hold_start'] = pd.to_numeric(adult_df_ctrl['hold_start'], errors='coerce')
    adult_df_ctrl = adult_df_ctrl[adult_df_ctrl.hold_start > -500].reset_index(drop = True)

    plot_param_over_time(adult_df_ctrl, param_list, 'Ctrl')

        # plots over time
    adult_df_hik = load_intr_data('high K')
    df_qc_hold_h = get_lab_book_info(adult_df_hik)
    df_qc_hold_h = get_match_cell_id(df_qc_hold_h)
    adult_df_hik = get_match_cell_id(adult_df_hik)

    # join QC and adult_df based on unique cell_ID
    adult_df_hik = pd.merge(adult_df_hik, df_qc_hold_h, on = 'cell_IDs_match', \
                        suffixes = (None, '_y')).reset_index(drop=True)

    adult_df_hik['hold_start'] = pd.to_numeric(adult_df_hik['hold_start'], errors='coerce')
    adult_df_hik = adult_df_hik[adult_df_hik.hold_start > -500].reset_index(drop = True)

    plot_param_over_time(adult_df_hik, param_list, 'high K')

    # mean plots
    mea_figure()
    mea_figure_normalized()

def plot_intrinsic_props():
    '''
    plot all the intrinsic properties for repatch data
    '''
    # df_repatch, df_slice = get_intrinsic_dfs()
    data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
        'paper_figs_collected_checked/data/'

    df_repatch = pd.read_excel(data_dir + 'repatch_data_complete.xlsx')
    df_slice = pd.read_excel(data_dir + 'slice_data_complete.xlsx')

    df_slice = df_slice[df_slice.treatment != 'fixed'].reset_index(drop = True)

    # plot distributions
    pl_intr.plot_hists(df_repatch, 'repatch')
    pl_intr.plot_hists(df_slice, 'slice')

    pl_intr.plot_param_for_days_no_15mM_slice(df_slice)
    pl_intr.plot_param_for_days_repatch_no_15mM(df_repatch)

    # get only firing cells
    df_repatch = pl_intr.remove_non_firing_cells_D1(df_repatch)
    df_slice = pl_intr.remove_non_firing_cells_D1(df_slice)

    # IFF and num_APs
    for dv in ['IFF', 'num_aps']:
        pl_intr.plot_iff_distribution(df_repatch, 'repatch', dv)
        pl_intr.plot_iff_avg_against_current_inj(df_repatch, 'firing cells repatch', dv, 17)
        # same for slice
        pl_intr.plot_iff_distribution(df_slice, 'slice', dv)
        pl_intr.plot_iff_avg_against_current_inj(df_slice, 'firing cells slice', dv, 17)

    # same but with 15 mM
    dest_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/' + \
    'results/human/paper_figs_collected_checked/plots/'

    df_repatch_15, df_slice_15 = get_intrinsic_dfs(True)
    df_slice_15 = df_slice_15[df_slice_15.treatment != 'fixed'].reset_index(drop = True)

    pl_intr.plot_param_for_days_slice(df_slice_15)
    pl_intr.plot_param_for_days_repatch(df_repatch_15)

    # get only firing cells
    df_repatch_15 = pl_intr.remove_non_firing_cells_D1(df_repatch_15)
    df_slice_15 = pl_intr.remove_non_firing_cells_D1(df_slice_15)

    # IFF and num_APs
    for dv in ['IFF', 'num_aps']:

        pl_intr.plot_iff_distribution(df_repatch_15, 'repatch_15', dv, dest_dir + 'iff_15/')
        pl_intr.plot_iff_avg_against_current_inj(df_repatch_15, 'firing cells repatch_15', dv, \
                                                 17, dest_dir + 'iff_15/')
        # same for slice
        pl_intr.plot_iff_distribution(df_slice_15, 'slice_15', dv, dest_dir + 'iff_15/')
        pl_intr.plot_iff_avg_against_current_inj(df_slice_15, 'firing cells slice_15', dv, \
                                                 17, dest_dir + 'iff_15/')

    # descriptive plots
    df_repatch.hold_start.hist()
    plt.xlabel('pA')
    plt.title('Holding start repatch')
    plt.savefig(dest_dir + 'qc_plots/hist/repatch_hold_start_repatch.png')
    plt.show()

    df_slice.hold_start.hist()
    plt.xlabel('pA')
    plt.title('Holding start slice')
    plt.savefig(dest_dir + '/qc_plots/hist/slice_hold_start_repatch.png')
    plt.show()

    # sets invalid (str) to NaN
    df_repatch.patient_age  = pd.to_numeric(df_repatch.patient_age, errors='coerce')
    df_repatch.patient_age.hist()
    plt.suptitle('Repatch')
    plt.xlabel('years')
    plt.savefig(dest_dir + 'age_distributions/repatch.png')
    plt.show()

    df_slice.patient_age  = pd.to_numeric(df_slice.patient_age, errors='coerce')
    df_slice.patient_age.hist()
    plt.suptitle('Slice')
    plt.xlabel('years')
    plt.savefig(dest_dir + 'age_distributions/slice.png')
    plt.show()

def load_intr_data(treatment):
    '''
    collects the intrinsic dfs 

    type (str): 'all', 'Ctrl', 'high K'
    '''
    df_intr_complete = collect_intrinsic_df() # for most updated version
    df_intr_complete = combine_columns(df_intr_complete)
    df_intr_complete = fix_cap_tau(df_intr_complete)
    # df_intr = pl_intr.get_column_RMPs_from_char(df_intr_complete)

    adult_df = pl_intr.filter_adult_hrs_incubation_data(df_intr_complete, min_age = 12, \
                                                        hrs_inc = 16, max_age = 151) # + QC criteria

    if treatment == 'all':
        return adult_df.reset_index(drop=True)

    return adult_df[adult_df.treatment == treatment].reset_index(drop=True)

def plot_param_over_time(adult_data, params, data_type):
    '''
    plots parameters over time
    args:
    params: list
    '''
    for param in params:
        adult_repatch = pl_intr.get_repatch_df(adult_data)
        adult_slice = adult_data.loc[adult_data['repatch'] == 'no'].dropna(subset=[param])
        dict_for_plotting = pl_intr.dict_for_plotting()

        fig, ax = plt.subplots(figsize = (8,6))
        ax.scatter(adult_repatch.hrs_after_OP, adult_repatch[param], \
                    c = 'orange', label = 'repatch', alpha = 0.45)
        ax.scatter(adult_slice.hrs_after_OP, adult_slice[param],\
                    c = 'green', label = 'slice', alpha = 0.3)

        # regression lines
        slope_r, intercept_r = np.polyfit(adult_repatch.hrs_after_OP, adult_repatch[param], 1)
        ax.plot(adult_repatch.hrs_after_OP, slope_r * adult_repatch.hrs_after_OP + intercept_r, \
                color = 'orange', linestyle = '--', linewidth = 2)

        slope, intercept = np.polyfit(adult_slice.hrs_after_OP, adult_slice[param], 1)
        ax.plot(adult_slice.hrs_after_OP, slope * adult_slice.hrs_after_OP + intercept, \
                color='green', linestyle = '--', linewidth = 2)

        # summary bins
        max_hours = int(adult_data['hrs_after_OP'].max())
        for start in range(10, max_hours + 1, 10):
            median_repatch = adult_repatch[(adult_repatch['hrs_after_OP'] >= start) & \
                                    (adult_repatch['hrs_after_OP'] < start + 10)][param].median()
            median_slice = adult_slice[(adult_slice['hrs_after_OP'] >= start) & \
                                        (adult_slice['hrs_after_OP'] < start + 10)][param].median()

            ax.plot([start], [median_repatch], 'orange', marker = 'D' , \
                markersize = 7, markeredgewidth = 4, alpha = 0.8)
            ax.plot([start], [median_slice], 'green', marker = 'D' , \
                markersize = 7, markeredgewidth = 4, alpha = 0.8)

        ax.set_xlabel('Hours after OP')
        ax.set_ylabel(dict_for_plotting[param][1])
        # handles, labels = ax.get_legend_handles_labels()
        fig.legend(ax.get_legend_handles_labels()[0], ax.get_legend_handles_labels()[1], \
                loc='upper right', fontsize='small')
        fig.suptitle(dict_for_plotting[param][0])

        # try:
        #     plt.show()
        # except (BrokenPipeError, IOError):
        #     pass
        dest_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
                    'paper_figs_collected_checked/plots/time_dependencies/'
        os.makedirs(dest_dir, exist_ok = True)
        fig.savefig(dest_dir + data_type + param + '.png')
        plt.close()

# MEA figure
def mea_figure():
    '''
    plot mea data
    '''
    mea_summary = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/' + \
                                'results/human/data/MEA/summary_pMEA.xlsx')

    k_incubation = mea_summary[mea_summary.Group == 'K']
    ctrl_incubation = mea_summary[mea_summary.Group == 'CTRL']

    fig, ax = plt.subplots(2, 1, figsize = (7,8))

    x_k = np.linspace(0.4, 0.8, len(k_incubation))
    ax[0].scatter(x_k, k_incubation.Baseline, c = 'magenta', label = 'aCSF + \n 8mM KCL', s = 45)
    ax[0].scatter(x_k + 1, k_incubation.Washout, marker = 'o', color ='darkblue', s = 45)
    ax[0].boxplot([k_incubation.Baseline, k_incubation.Washout], positions=[0.25, 2], widths=0.1)

    for i, val in enumerate(k_incubation.Baseline):
        ax[0].plot([ x_k[i], x_k[i]+1], [val, k_incubation.Washout.values[i]], \
                zorder = -1, c = 'darkgrey', alpha = 0.7)

    x_ctrl = np.linspace(0.4, 0.8, len(ctrl_incubation))
    ax[1].scatter(x_ctrl, ctrl_incubation.Baseline, c = 'cyan', label = 'aCSF', s = 45)
    ax[1].scatter(x_k + 1, ctrl_incubation.Washout, marker = 'o', color ='darkblue', s = 45)
    ax[1].boxplot([ctrl_incubation.Baseline, ctrl_incubation.Washout], positions=[0.25, 2], \
                  widths=0.1)

    for i, val in enumerate(ctrl_incubation.Baseline):
        ax[1].plot([x_ctrl[i], x_ctrl[i]+1], [val, ctrl_incubation.Washout.values[i]], \
                zorder = -1, c = 'darkgrey', alpha = 0.7)

    # ax[0].set_title('K Incubation')
    ax[0].set_xticks(ticks = [0.5, 1.5], labels = ['aCSF + \n 8mM KCL', 'aCSF'])
    ax[0].set_ylabel('# spikes / electrode')
    # Add secondary x-axis label
    secax = ax[0].secondary_xaxis('top') # try with 0 and top
    secax.set_xticks([0.5, 1.5])
    secax.set_xticklabels(['Incubation \n solution', 'Recording \n solution'])
    secax.spines['top'].set_visible(False)
    ax[0].tick_params(axis='x', which = 'minor', labelsize=45)

    ax[1].set_xticks(ticks = [0.5, 1.5], labels = ['aCSF', 'aCSF'])
    ax[1].set_yticks(ticks = [0, 200, 400, 600])
    # ax[1].set_title('Ctrl Incubation')
    ax[1].set_ylabel('# spikes / electrode')

    #handles, labels = ax[0].get_legend_handles_labels()
    legend = fig.legend(bbox_to_anchor = (0.97, 0.9), fontsize = 12, \
                        title = 'Incubation \n solution')
    plt.setp(legend.get_title(),fontsize = 14)
    plt.subplots_adjust(left = 0.2, hspace = 0.4)
    plt.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
                'paper_figs_collected_checked/plots/mea/mEA.png')
    # plt.show()

def mea_figure_normalized():
    '''
    plot mea data normalized to the recording solution values
    '''
    mea_summary = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/' + \
                                'results/human/data/MEA/summary_pMEA.xlsx')

    k_incubation = mea_summary[mea_summary.Group == 'K']
    ctrl_incubation = mea_summary[mea_summary.Group == 'CTRL']

    fig, ax = plt.subplots(2, 1, figsize = (7,8))

    x_k = np.linspace(0.4, 0.8, len(k_incubation))
    ax[0].scatter(x_k, k_incubation.Baseline / k_incubation.Washout , c = 'magenta', \
                  label = 'aCSF + \n 8mM KCL', s = 45)
    ax[0].scatter(x_k + 1, k_incubation.Washout /k_incubation.Washout , marker = 'o', \
                  color ='darkblue', s = 45)
    ax[0].boxplot([k_incubation.Baseline / k_incubation.Washout, \
                   k_incubation.Washout /k_incubation.Washout ], positions=[0.25, 2], widths=0.1)

    for i, val in enumerate(k_incubation.Baseline / k_incubation.Washout):
        ax[0].plot([ x_k[i], x_k[i]+1], [val, k_incubation.Washout.values[i] /  \
                    k_incubation.Washout.values[i]],zorder = -1, c = 'darkgrey', alpha = 0.7)

    x_ctrl = np.linspace(0.4, 0.8, len(ctrl_incubation))
    ax[1].scatter(x_ctrl, ctrl_incubation.Baseline / ctrl_incubation.Washout, \
                  c = 'cyan', label = 'aCSF', s = 45)
    ax[1].scatter(x_k + 1, ctrl_incubation.Washout / ctrl_incubation.Washout, \
                  marker = 'o', color ='darkblue', s = 45)
    ax[1].boxplot([ctrl_incubation.Baseline / ctrl_incubation.Washout, \
            ctrl_incubation.Washout / ctrl_incubation.Washout], positions=[0.25, 2], widths=0.1)

    for i, val in enumerate(ctrl_incubation.Baseline / ctrl_incubation.Washout):
        ax[1].plot([x_ctrl[i], x_ctrl[i]+1], [val, ctrl_incubation.Washout.values[i] / \
                    ctrl_incubation.Washout.values[i]], zorder = -1, c = 'darkgrey', alpha = 0.7)

    # ax[0].set_title('K Incubation')
    ax[0].set_xticks(ticks = [0.5, 1.5], labels = ['aCSF + \n 8mM KCL', 'aCSF'])
    ax[0].set_yticks(ticks = [5, 10, 15, 20])
    ax[0].set_ylabel('# spikes / electrode')
    # Add secondary x-axis label
    secax = ax[0].secondary_xaxis('top') # try with 0 and top
    secax.set_xticks([0.5, 1.5])
    secax.set_xticklabels(['Incubation \n solution', 'Recording \n solution'])
    secax.spines['top'].set_visible(False)
    ax[0].tick_params(axis='x', which = 'minor', labelsize=45)

    ax[1].set_xticks(ticks = [0.5, 1.5], labels = ['aCSF', 'aCSF'])
    ax[1].set_yticks(ticks = [5, 10, 15, 20])
    # ax[1].set_title('Ctrl Incubation')
    ax[1].set_ylabel('# spikes / electrode')

    #handles, labels = ax[0].get_legend_handles_labels()
    legend = fig.legend(bbox_to_anchor = (.93, 0.9), fontsize = 12, \
                        title = 'Incubation \n solution' )
    plt.subplots_adjust(left=0.15, hspace=0.4)
    plt.setp(legend.get_title(),fontsize = 14)
    # plt.show()
    plt.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
        'paper_figs_collected_checked/plots/mea/mea_normalized.png')

def get_lab_book_info(df, human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    '''
    collects info from the lab books
    '''
    patcher_dict = {
        'Verji': 'data_verji/',
        'Rosie': 'data_rosie/'
        }
    for patcher, dir_p in patcher_dict.items():

        for i, op in enumerate(df[df.patcher == patcher].OP.unique()):

            lab_book = sort.get_lab_book(human_dir + dir_p + op + '/')[0]

            if 'sl' not in lab_book.columns:
                print('fix lab book for', op)
            mask_depth = lab_book.columns.str.contains('depth', case=False, na = False)
            mask_seal = lab_book.columns.str.contains('seal', case=False, na = False)

            if 'Unnamed: 7' in lab_book.columns:
                holdings_end = lab_book['Unnamed: 7'][1:].tolist()
            else:
                holdings_end = []

            sl_indx = lab_book.index[lab_book['sl'].notnull()]
            all_slice_names = sort.fix_slice_names(lab_book['sl'][sl_indx].tolist(), sl_indx)

            slice_names = []
            for nam in all_slice_names:
                if 'd2' in nam:
                    slice_names.append(nam[:3] + 'D2')
                else:
                    slice_names.append(nam)
            slice_names = slice_names[:len(lab_book.channels[1:].tolist())]

            if len(slice_names) < len(lab_book.channels[1:].tolist()):
                slice_names.extend([float('nan')] * \
                                   (len(lab_book.channels[1:].tolist()) - len(slice_names)))

            hold_df = pd.DataFrame({'patcher': patcher, 'OP': op, 'slice':slice_names,
                'cell_ch':lab_book.channels[1:].tolist(), 
                'hold_start': lab_book['Holding (pA)'][1:].tolist(), 
                'depth': lab_book[lab_book.columns[mask_depth][0]][1:].tolist(), 
                'seal': lab_book[lab_book.columns[mask_seal][0]][1:].tolist(), 
                'bridge_balance': lab_book['Bridge balance (MΩ)'][1:].tolist(), 
                'capacitance_neutr': lab_book['Capacitance neutralization (pF) '][1:].tolist()})

            if len(holdings_end) > 0:
                hold_df.insert(5, 'hold_end', holdings_end)
            hold_df = hold_df.dropna(subset = ['cell_ch'])

            if i == 0 and patcher == 'Verji':
                df_qc_hold = hold_df
                continue

            df_qc_hold = pd.concat([df_qc_hold, hold_df], ignore_index = True)

    df_qc_hold = df_qc_hold.reset_index(drop = True)
    # df_qc_hold.to_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+ \
    #                     'data/summary_data_tables/metadata/holdings.xlsx')
    return df_qc_hold

def get_sample_size(adult_df):
    '''
    returns the sample size for slice and repatch in a dictionary
    '''
    adult_repatch = pl_intr.get_repatch_df(adult_df)
    adult_slice = adult_df.loc[adult_df['repatch'] == 'no']
    return {
            'repatch measurements (all)': len(adult_repatch),
            'repatched cells' : len(adult_repatch.cell_ID_new.unique()),
            'repatch Ctrl cells' : len(adult_repatch[adult_repatch['treatment'] == 'Ctrl'])/ 2,
            'repatch high K cells' : len(adult_repatch[adult_repatch['treatment'] == 'high K'])/2,
            'repatch patients' : len(adult_repatch.OP.unique()),
            'slice all cells' : len(adult_slice),
            'slice Ctrl cells': len(adult_slice[adult_slice['treatment'] == 'Ctrl']),
            'slice high K cells' : len(adult_slice[adult_slice['treatment'] == 'high K']),
            'slice Ctrl D1' : len(adult_slice[(adult_slice['treatment'] == 'Ctrl') & \
                                  (adult_slice['day'] == 'D1')]),
            'slice Ctrl D2' : len(adult_slice[(adult_slice['treatment'] == 'Ctrl') & \
                                  (adult_slice['day'] == 'D2')]),
            'slice high K D1' : len(adult_slice[(adult_slice['treatment'] == 'high K') & \
                                  (adult_slice['day'] == 'D1')]),
            'slice high K D2' : len(adult_slice[(adult_slice['treatment'] == 'high K') & \
                                  (adult_slice['day'] == 'D2')]),
            'slice patients' : len(adult_slice.OP.unique())
        }

def combine_columns(df):
    '''
    Combines columns in df that were measuring ramp_rheo into one column called rheo_ramp_c
    deletes other rheo ramp columns
    Parameters:
    df (pandas.DataFrame): The input DataFrame containing the columns to be combined.
    Returns:
    pandas.DataFrame: The DataFrame with the new 'rheo_ramp' column added.
    '''

    # combine ramp columns
    rheo_ramp_c = []
    for a in zip(df['Rheobse_ramp'], df['rheo_ramp'], df['rheos_ramp']):
        rheo_ramp_c.append(np.nanmean(a))

    col_delete = df.columns[df.columns.str.contains('unnamed', case=False)].tolist() + \
        df.columns[df.columns.str.contains('ramp', case=False)].tolist() + ['comment']
    df = df.drop(columns=col_delete)
    df.insert(len(df.columns), 'rheo_ramp_c', rheo_ramp_c)
    df = df.reset_index(drop=True)

    return df

def get_match_cell_id(df):
    '''
    cell_ID = patcher_initial + OP_date + slice + channel
    '''
    if 'cell_IDs_match' in df.columns:
        df = df.drop(columns = ['cell_IDs_match'])

    if 'cell_IDs_match_re' in df.columns:
        df = df.drop(columns = ['cell_IDs_match_re'])

    patcher_dict = {'Verji':'vm', 'Rosie': 'rs'}

    cell_IDs_match, cell_IDs_match_re = [], []
    for i in range(len(df)):
        chan = df.cell_ch[i]
        slic = df.slice[i]
        if isinstance(chan, (float, int, np.int64)):
            chan = str(int(chan))
            chan_re = chan
            slic_re = slic
        else:
            chan_re = str(chan)[-1] # ch repatched form D1
            slic_re = slic[:-3]
            chan = str(chan)[0] # active ch D2

        cell_IDs_match.append(patcher_dict[df.patcher[i]] +
                              df.OP[i][2:] +
                              slic +
                              'c' + chan)
        cell_IDs_match_re.append(patcher_dict[df.patcher[i]] +
                              df.OP[i][2:] +
                              slic_re +
                              'c' + chan_re)

    df.insert(len(df.columns), 'cell_IDs_match', cell_IDs_match)
    df.insert(len(df.columns), 'cell_IDs_match_re', cell_IDs_match_re)

    return df

def find_repeats(df):
    '''
    finds repeats
    '''
    if 'cell_IDs_match' not in df.columns:
        df = get_match_cell_id(df)
    counts = Counter(df.cell_IDs_match.tolist())
    repeats = [item for item, count in counts.items() if count > 1]

    if len(repeats) > 0:
        print('repeats adults_df', repeats)

def data_recordings_on_both_days(df):
    '''
    1. finds nidices for slices != 2 and != 5
    2. finds indices where there is no slice data for D2
    3. filters data based on the above
    returns filtered df
    '''
    # take care of special slices
    # where it's not the usual length of 2 in D1 and len of 5 in D2

    # remove OPs where no recording on second day
    op_d2 = df.OP[df.day.str.contains('D2', case=False).tolist()].unique()
    op_remove = list(set(df.OP.unique()) - set(op_d2))
    # addding OP210615 from Rosie, no repatches and difficult to analyze
    op_remove.append('OP210615')

    df = df[~df['OP'].isin(op_remove)] # some data for D2 available

    special_indx = []
    for sl in df.slice:
        if len(sl) != 2 and len(sl) != 5:
            if sl[-1] == '':
                continue
            # special slices to keep
            special_indx.append(df[df['slice'] == sl].index[0])

    # take the indices where for any slice on D1 there is a recording on D2
    true_indices = []
    for patcher in ['Rosie', 'Verji']:
        ops_patcher = df.OP[df['patcher'] == patcher].unique()
        for op in ops_patcher:
            df_p_op = df[(df.patcher == patcher) & (df.OP == op)]
            sl_d1 = df_p_op.slice[df.day == 'D1'].unique()
            for sl in sl_d1:
                if sl + '_D2' in df_p_op['slice'].values:
                    mask_keep = (df_p_op['slice'] == sl) | (df_p_op['slice'] == sl + '_D2')
                    true_indices = true_indices + df_p_op.index[mask_keep].tolist()

    true_indices = list(set(true_indices + special_indx))
    # false_indc = list(set(list(range(len(df)))) - set(true_indices))

    df = df.loc[true_indices].reset_index(drop = True)
    return df

def check_join_worked(df1, df2):
    '''
    checks the cell_IDs_match
    '''
    if 'cell_IDs_match' not in df1.columns or 'cell_IDs_match' not in df2.columns:
        print('fix your dfs, add cell_IDs_match columns')
        return
    # check if join operation worked
    # if the lists are not identics, aka if False in lists
    # print that something is wrong
    if False in list((df1.cell_IDs_match) == (df1.cell_IDs_match)):
        print('error')

# IFF functions
def check_no_iff_ops(human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    '''
    for all OP folders in human_dir
    checkl if IFF exists
    returns no_iff_dict for which no IFF found by patcher
    '''

    patcher_dict = {
        'Verji': 'data_verji/',
        'Rosie': 'data_rosie/'
        }

    no_iff_dict = {}
    for patcher, p_dir  in patcher_dict.items():
        patcher_dirs = glob(human_dir + p_dir + 'OP*')
        ops_no_iff = []
        for op_dir in patcher_dirs:
            if os.path.exists(os.path.join(op_dir, 'data_tables')):
                IFF = glob(os.path.join(op_dir, 'data_tables/') + '*IFF*')
                if len(IFF) == 0:
                    # print('no IFF df for', op_dir[op_dir.rfind('/')+1:], patcher)
                    ops_no_iff.append(op_dir[op_dir.rfind('/')+1:])
        no_iff_dict[patcher] = sorted(ops_no_iff)
    return no_iff_dict

def check_where_iff_needs_creating(df):
    '''
    based on the surgeries that are going into the slice 
    and repatch analysis
    '''

    df_r = df.OP[df['patcher'] == 'Rosie'].unique()
    df_v = df.OP[df['patcher'] == 'Verji'].unique()
    no_iff_dict = check_no_iff_ops()
    op_add_iff_r = sorted(list(set(df_r) & set(no_iff_dict['Rosie'])))
    op_add_iff_v = sorted(list(set(df_v) & set(no_iff_dict['Verji'])))

    if len(op_add_iff_r) > 0 or len(op_add_iff_v) > 0:
        print('Analyse IFF data for ', op_add_iff_r, 'or', op_add_iff_v)
    else:
        print('nothing to add, continue with peace of mind')

def qc_check_wrong_day(df):
    '''
    looks through the slices for each entry
    if the name of the slice is longer than 2, the day should be D2
    the entries for which this is not true are printed
    '''

    for i, sl in enumerate(df.slice):
        if len(sl) > 2 and df.day[i] != 'D2':
            print('Fix data for', df.OP[i], sl, df.day[i])

def copy_file_if_exists(src_path, dest_dir):
    '''
    Copies a file from src_path to dest_dir if it exists.

    Args:
        src_path (str): The source file path.
        dest_dir (str): The destination directory.
    '''
    # Check if the file exists
    if os.path.exists(src_path):
        # Copy the file to the destination directory
        shutil.copy(src_path, dest_dir)
    else:
        print(f"File '{src_path[src_path.rfind('/')+1:]}' does not exist")

def copy_plots_for_figures(df, data_type):
    '''
    Copies max spikes, onset and AP plots for each of the cells in the df to the qc_plots directory.
    data_type (str): slice, repatch
    '''
    human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
    patcher_dict = {'Verji': 'data_verji/',
    'Rosie': 'data_rosie/'}

    qc_plots = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human'+\
        '/paper_figs_collected_checked/plots/qc_plots/'

    for i, op in enumerate(df.OP):
        patcher = df.patcher[i]
        fn = df.filename[i]
        chan = str(df.cell_ch[i])

        op_plots = human_dir + patcher_dict[patcher] + op + '/plots/'

        src_path_spikes = op_plots + 'Max_Spikes/char_spikes_plot_' + fn[:-4] + '_Ch'+ chan + '.png'
        copy_file_if_exists(src_path_spikes, qc_plots + data_type + '/max_spikes/')

        # src_path_onset = op_plots + 'Onset/Char_onset_plot_' + fn[:-4] + '_' +chan + '.png'
        # copy_file_if_exists(src_path_onset, qc_plots + '/slice/onset/')

        # for one surgery those are the correct names
        src_path_onset = op_plots + 'Onset/Char_onset_plot_' + fn[:-4] + '_Ch' +chan + '.png'
        copy_file_if_exists(src_path_onset, qc_plots + data_type + '/onset/')

        src_path_ap = op_plots + 'AP_props/' + fn[:-4] + '_AP#' + str(2) + '_Ch'+ chan + '.png'
        copy_file_if_exists(src_path_ap, qc_plots + data_type +'/ap_props/')

def fix_cap_tau(df, human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    '''
    calcualtes the with the same way the capacitance and tau
    takes the 2nd value after the fit, not the first
    because the first sometimes falls on the artifact
    '''
    op_skip = ['OP210319']
    patcher_dict = {
        'Verji': 'data_verji/',
        'Rosie': 'data_rosie/'
        }
    taus, caps = [], []
    for i, fn in enumerate(df.filename):
        file_path = human_dir + patcher_dict[df.patcher[i]] + df.OP[i] + '/' + fn
        chan = [df.cell_ch[i]]
        op = df.OP[i]
        if op in op_skip:
            taus.append(np.nan)
            caps.append(np.nan)
            continue
        # currently function with no plotting
        tau, cap = cap_values_adjustment(file_path, chan) # when injection == - 200 pA
        if op == 'OP201029' and fn == '20o30008.abf':
            if chan ==[1]:
                tau = 7.95
                cap = 147.095200
            if chan ==[2]:
                tau = 8.10
                cap = 159.220600
        taus.append(tau)
        caps.append(cap)

    df.insert(len(df.columns), 'tau_adj', taus)
    df.insert(len(df.columns), 'cap_adj', caps)

    return df

def plots_comapre_old_new_tau_cap(df, data_type):
    '''
    plots the old and new tau and cap values
    '''
    for i, tau in enumerate(df.tau_adj):
        tau2 = df.membra_time_constant_tau[i]
        if abs(tau - tau2) > 0.15 * tau:
            print(df.OP[i], df.filename[i], df.cell_ch[i], tau, tau2)
        plt.scatter(0, tau, c = 'blue', alpha = 0.7)
        plt.scatter(1, tau2, c = 'red', alpha = 0.7)
        plt.plot([0, 1], [tau, tau2], c = 'blue', alpha = 0.6)
        plt.title(data_type + ' tau values')
    plt.show()

    for i, cap in enumerate(df.cap_adj):
        cap2 = df.capacitance[i]
        if abs(cap - cap2) > 0.15 * cap:
            print(df.OP[i], df.filename[i], df.cell_ch[i], cap, cap2)
        plt.scatter(0, cap, c = 'blue', alpha = 0.7)
        plt.scatter(1, cap2, c = 'red', alpha = 0.7)
        plt.plot([0, 1], [cap, cap2], c = 'blue', alpha = 0.6)
        plt.title(data_type + ' capacitance values')
    plt.show()

def get_intrinsic_dfs(_15mm = False):
    '''
    collect the QC_checked dataframes
    QC
    treatments: Ctrl, high K (only 8 mM)
    Age: 12 to max (75?)
    hours incubation > 16h
    columns from lab book during recording

    returns: adult_df_repatch, adult_df_slice with all info for plotting intrinsic properties
    '''
    adult_df = load_intr_data('all')
    # treatments to include/ exclude
    mask_exclude = (adult_df['treatment'] == 'TTX') | \
        ((adult_df['treatment'] == 'high K') & (adult_df['high K concentration'] == '15 mM'))
    if _15mm:
        mask_exclude = adult_df['treatment'] == 'TTX'
    adult_df = adult_df[~mask_exclude]
    adult_df = data_recordings_on_both_days(adult_df)

    check_where_iff_needs_creating(adult_df)
    df_qc_hold = get_lab_book_info(adult_df)
    df_qc_hold = get_match_cell_id(df_qc_hold)
    adult_df = get_match_cell_id(adult_df)

    iffs = collect_IFF_dfs()
    # remove ops not present in adult_df_qc
    ops_remove = set(iffs.OP.unique()) - set(adult_df.OP.unique())
    iffs = iffs[~iffs['OP'].isin(ops_remove)].reset_index(drop = True)
    iffs = get_match_cell_id(iffs)

    # find fuplicates
    # Filter the DataFrame to get the repeating values
    fix_op = sorted(iffs.OP[iffs['cell_IDs_match'].duplicated(keep=False)].unique())
    for op in fix_op:
        print(op)

    if len(adult_df.cell_IDs_match.unique()) < len(adult_df):
        print('smth wrong with unique cell_IDs adult')

    if len(df_qc_hold.cell_IDs_match.unique()) < len(df_qc_hold):
        print('smth wrong with unique cell_IDs QC')

    if len(iffs.cell_IDs_match.unique()) < len(iffs):
        print('smth wrong with unique cell_IDs iff')

    cells_not_in_qc=set(adult_df.cell_IDs_match.unique())-set(df_qc_hold.cell_IDs_match.unique())
    cells_not_in_iffs =set(adult_df.cell_IDs_match.unique()) - set(iffs.cell_IDs_match.unique())
    missing_cells = list(cells_not_in_qc) + list(cells_not_in_iffs)
    print('missign cells [QC], [IFF]', missing_cells)

    # remove the adult_df
    adult_df = adult_df[~adult_df['cell_IDs_match'].isin(missing_cells)]
    adult_df = adult_df.reset_index(drop=True)

    # join QC and adult_df based on unique cell_ID
    adult_df_qc = pd.merge(adult_df, df_qc_hold, on = 'cell_IDs_match', \
                        suffixes = (None, '_y')).reset_index(drop=True)

    adult_df_complete = pd.merge(adult_df_qc, iffs, on = 'cell_IDs_match', \
                        suffixes = (None, '_y_iff')).reset_index(drop=True)

    # if slices in new and old df identical --> remove the columns that contain '_y'
    # remove repeating columns, without the cell_IDs_match_re_y
    if False not in list((adult_df_qc.slice) == (adult_df_complete.slice_y)):
        col_remove = list(adult_df_complete.columns[adult_df_complete.columns.str.contains('_y')])
        col_remove.remove('cell_IDs_match_re_y')

        adult_df_complete = adult_df_complete.drop(columns = col_remove)
    adult_df_complete=adult_df_complete[adult_df_complete.hold_start>-500].reset_index(drop = True)
    qc_check_wrong_day(adult_df_complete)
    adult_df_repatch = pl_intr.get_repatch_df(adult_df_complete)
    adult_df_slice = adult_df_complete[adult_df_complete.repatch == 'no']
    adult_df_slice = data_recordings_on_both_days(adult_df_slice)

    # adult_df_qc.to_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/
    # paper_figs_collected_checked/adult_df_complete.xlsx')
    # save for stats
    # adult_repatch = pl_intr.get_repatch_df(adult_df)
    # adult_repatch.to_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+ \
    #                      'data/summary_data_tables/intrinsic_properties/2025-02-05_repatch.csv')
    return adult_df_repatch, adult_df_slice

# maybe can add a plot on the gradual effects of high K
# also some MEA data on that

# keep cells until they've reached maximum firing


# createa a specific dataframe for the rheos_ramp. for repatch adn slice
# summary stats patients


# plot the ones that don't fire on D2
# remove them from the IFF plots

if __name__ == "__main__":
    main()

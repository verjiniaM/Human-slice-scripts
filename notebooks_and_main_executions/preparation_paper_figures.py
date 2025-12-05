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
import ephys_analysis.funcs_for_results_tables as funcs_results
import ast

#%%
def main():
    '''
    plots intrinsic and time dependency data
    '''
    # plot_time_dep()
    # plot_intrinsic_props()

    data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
        'paper_figs_collected_checked/data/'

    df_repatch = pd.read_excel(data_dir + 'repatch_data_temporal.xlsx')
    df_slice = pd.read_excel(data_dir + 'slice_data_temporal.xlsx')

    # # get only cells where holding < -200 pA
    # df_repatch = df_repatch[df_repatch.hold_start > min_holding].reset_index(drop = True)
    # df_slice = df_slice[df_slice.hold_start > min_holding].reset_index(drop = True)

    df_slice = df_slice[df_slice.treatment != 'fixed'].reset_index(drop = True)

    cor1 = Correlations(df_slice, df_repatch)
    # cor1.sl_repatch_separate_means_confidence_ints()
    cor1.sl_repatch_separate_means_D1_all_ctrl()

def get_passing_meta_df(analysis_df, meta_df):
    # meta_complete = pd.read_excel("/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/meta_data.xlsx")

    ops, tissue_source = [], []
    for OP in analysis_df.OP.unique():
        ops.append(OP)
        tissue_source.append(analysis_df.tissue_source[analysis_df['OP'] == OP].to_list()[0])

    OPs_interest = pd.DataFrame({'OP': ops, 'tissue_source': tissue_source})
    relevant_meta = pd.merge(OPs_interest, meta_df)

    return relevant_meta
 
def plot_intrinsic_props(min_holding = -200 ):
    '''
    plot all the intrinsic properties for repatch data
    '''
    # df_repatch, df_slice = get_intrinsic_dfs()
    data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
        'paper_figs_collected_checked/data/'

    df_repatch = pd.read_excel(data_dir + 'repatch_data_temporal.xlsx')
    df_slice = pd.read_excel(data_dir + 'slice_data_temporal.xlsx')

    # # get only cells where holding < -200 pA
    # df_repatch = df_repatch[df_repatch.hold_start > min_holding].reset_index(drop = True)
    # df_slice = df_slice[df_slice.hold_start > min_holding].reset_index(drop = True)

    df_slice = df_slice[df_slice.treatment != 'fixed'].reset_index(drop = True)

    # plot distributions
    pl_intr.plot_hists(df_repatch, 'repatch')
    pl_intr.plot_hists(df_slice, 'slice')


    # pl_intr.plot_param_for_days_no_15mM_slice(df_slice)
    pl_intr.plot_param_for_days_repatch_no_15mM(df_repatch)
    pl_intr.plot_param_for_days_no_15mM_slice_v2(df_slice)

    # get only firing cells
    df_repatch = pl_intr.remove_non_firing_cells_D1(df_repatch)
    df_slice = pl_intr.remove_non_firing_cells_D1(df_slice)

    # IFF and num_APs
    for dv in ['IFF', 'num_aps']:
        # pl_intr.plot_iff_distribution(df_repatch, 'repatch', dv)
        pl_intr.plot_iff_avg_against_current_inj_plus_hist(df_repatch, 'firing cells repatch', dv, 13, True)

        # same for slice
        # pl_intr.plot_iff_distribution(df_slice, 'slice', dv)
        pl_intr.plot_iff_avg_against_current_inj_plus_hist(df_slice, 'firing cells slice', dv, 13, True)
        
    # plot non-firing cells
    pl_intr.plot_non_firing_cells_d2(df_repatch, 'repatch')
    pl_intr.plot_non_firing_cells_d2(df_slice, 'slice')

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
                                                 12, dest_dir + 'iff_15/')
        # same for slice
        pl_intr.plot_iff_distribution(df_slice_15, 'slice_15', dv, dest_dir + 'iff_15/')
        pl_intr.plot_iff_avg_against_current_inj(df_slice_15, 'firing cells slice_15', dv, \
                                                 12, dest_dir + 'iff_15/')

    pl_intr.plot_non_firing_cells_d2(df_repatch_15, 'repatch_15')
    pl_intr.plot_non_firing_cells_d2(df_slice_15, 'slice_15')

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
    df_intr_complete = funcs_results.collect_intrinsic_df() # for most updated version
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

        # regression liness
        slope_r, intercept_r = np.polyfit(adult_repatch.hrs_after_OP, adult_repatch[param], 1)
        ax.plot(adult_repatch.hrs_after_OP, slope_r * adult_repatch.hrs_after_OP + intercept_r, \
                color = 'orange', linestyle = '--', linewidth = 2)

        slope, intercept = np.polyfit(adult_slice.hrs_after_OP, adult_slice[param], 1)
        ax.plot(adult_slice.hrs_after_OP, slope * adult_slice.hrs_after_OP + intercept, \
                color='green', linestyle = '--', linewidth = 2)

        # summary bins
        max_hours = int(adult_data['hrs_after_OP'].max())

        medians_repatch, medians_slice = [], []
        for start in range(0,  50, 10):
            condition_repatch = (adult_repatch['hrs_after_OP'] >= start) \
                & (adult_repatch['hrs_after_OP'] < start + 10)
            condition_slice = (adult_slice['hrs_after_OP'] >= start) & \
                (adult_slice['hrs_after_OP'] < start + 10) # gives a bool

            medians_repatch.append(np.median(adult_repatch[param][condition_repatch]))
            medians_slice.append(np.median(adult_slice[param][condition_slice]))

        x_medians_r = np.arange(5, 55, 10)
        x_medians_s = np.arange(5, 55, 10)
        ax.plot(x_medians_r, medians_repatch, 'orange', marker = 'D' , \
                markersize = 7, markeredgewidth = 4, alpha = 0.8)
        ax.plot(x_medians_s, medians_slice, 'green', marker = 'D' , \
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
        # fig.savefig(dest_dir + data_type + param + '.png')
        fig.savefig(dest_dir + data_type + param + '.svg', \
            format = 'svg', bbox_inches = 'tight', dpi = 1000)
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
        if 'cell_ch' in df.columns:
            chan = df.cell_ch[i]
        elif 'Channels to use' in df.columns:
            chan = df['Channels to use'][i]
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
    
    return repeats

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
    adult_df = adult_df.reset_index(drop=True)
    # adult_df = data_recordings_on_both_days(adult_df)

    check_where_iff_needs_creating(adult_df)
    df_qc_hold = get_lab_book_info(adult_df)
    df_qc_hold = get_match_cell_id(df_qc_hold)
    adult_df = get_match_cell_id(adult_df)

    iffs = funcs_results.collect_IFF_dfs()
    # remove ops not present in adult_df_qc
    ops_remove = set(iffs.OP.unique()) - set(adult_df.OP.unique())
    iffs = iffs[~iffs['OP'].isin(ops_remove)].reset_index(drop = True)
    iffs = get_match_cell_id(iffs)

    # find duplicates
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
    adult_df_slice_all = adult_df_complete[adult_df_complete.repatch == 'no']
    adult_df_slice_no_rep = data_recordings_on_both_days(adult_df_slice_all)
    unused_df = adult_df_slice_all[~adult_df_slice_all['cell_IDs_match'].isin(adult_df_slice_no_rep['cell_IDs_match'])]

    # adult_df_qc.to_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/
    # paper_figs_collected_checked/adult_df_complete.xlsx')
    # save for stats
    # adult_repatch = pl_intr.get_repatch_df(adult_df)
    # adult_repatch.to_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+ \
    #                      'data/summary_data_tables/intrinsic_properties/2025-02-05_repatch.csv')
   #  return adult_df_repatch, adult_df_slice_no_rep

def get_spontan_df():
    '''
    collects all spontan dfs
    cross compares with alreadhy QC checked intrinsic data
    saves the results for statistical analysis
    '''
    # collecting and cross checking spontan / mini dataframes
    meta_df_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/manual_events_analysis/meta_dfs/'

    # collect event meta_dfs
    meta_spontan = funcs_results.collect_events_dfs('spontan')
    # remove cells where nothing to analyse and OP210615
    meta_spontan = meta_spontan[meta_spontan['OP'] != 'OP210615']
    meta_spontan = meta_spontan[meta_spontan['OP'] != 'OP210319']
    meta_spontan = meta_spontan[meta_spontan['OP'] != 'OP210323']
    meta_spontan = meta_spontan[meta_spontan['swps_to_analyse'] != '[]'].reset_index(drop=True)
    meta_spontan = get_match_cell_id(meta_spontan)
    repeats = find_repeats(meta_spontan)

    # load QC checked intrinsic data
    data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
    df_repatch = pd.read_excel(data_dir + 'repatch_data_temporal.xlsx')
    df_repatch.insert(len(df_repatch.columns), 'source', ['repatch'] * len(df_repatch))
    df_slice = pd.read_excel(data_dir + 'slice_data_temporal.xlsx')
    df_slice.insert(len(df_slice.columns), 'source', ['slice'] * len(df_slice))
    slice_inc = pd.read_excel(data_dir + '10.07.25_slice_incubation_only_no_move.xlsx')
    slice_inc.insert(len(slice_inc.columns), 'source', ['incubation'] * len(slice_inc))
    # combine in one df
    df_intr_complete = pd.concat([df_repatch, df_slice, slice_inc], ignore_index = True)
    df_intr_complete.reset_index(drop = True, inplace = True)

    df_intr_complete = get_match_cell_id(df_intr_complete)
    find_repeats(df_intr_complete) # Check if worked

    # with cell_ID from slice
    cells_not_in_spontan = set(df_intr_complete.cell_IDs_match.unique()) - \
        set(meta_spontan.cell_IDs_match.unique())
    print('removing ', str(len(cells_not_in_spontan)), 'cells,missing from spontan_meta')
    # remove from intr_df
    df_intr_complete = df_intr_complete[~df_intr_complete['cell_IDs_match'].isin(cells_not_in_spontan)]
    _remove = list(df_intr_complete.columns[df_intr_complete.columns.str.contains('_y')]) # all doubles
    df_intr_complete = df_intr_complete.drop(columns = _remove)
    df_intr_complete = df_intr_complete.reset_index(drop=True)

    # join QC and adult_df based on unique cell_ID
    spontan_qc = pd.merge(meta_spontan, df_intr_complete, on = 'cell_IDs_match', \
                        suffixes = (None, '_y')).reset_index(drop=True)

    # remove unnecessary columns
    col_remove = list(spontan_qc.columns[spontan_qc.columns.str.contains('_y')]) # all doubles
    col_remove = list(set(col_remove) - set(['cell_ID_y', 'treatment_y', 'hrs_incubation_y'])) # keep these
    col_remove += ['cell_ID', 'treatment', 'hrs_incubation']
    spontan_qc = spontan_qc.drop(columns = col_remove)
    spontan_qc = spontan_qc.loc[:, ~ spontan_qc.columns.str.contains('Unnamed')]
    # Rename the columns. to keep the ones from intr
    spontan_qc = spontan_qc.rename(columns={
        'cell_ID_y': 'cell_ID',
        'treatment_y': 'treatment', 
        'hrs_incubation_y': 'hrs_incubation'
    })
    qc_check_wrong_day(spontan_qc)

    spontan_repatch = pl_intr.get_repatch_df(spontan_qc)
    spontan_slice_all = spontan_qc[spontan_qc.repatch == 'no']
    spontan_slice_no_rep = data_recordings_on_both_days(spontan_slice_all)
    spontan_inc_only = spontan_qc[spontan_qc.source == 'incubation']

    save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/events'
    spontan_repatch.to_csv(save_dir + '/spontan_repatch.csv')
    spontan_slice_all.to_csv(save_dir + '/spontan_slice_all.csv')
    spontan_inc_only.to_csv(save_dir + '/spontan_inc_only.csv')
    
    spontan_repatch.to_excel(save_dir + '/spontan_repatch.xlsx')
    spontan_slice_all.to_excel(save_dir + '/spontan_slice_all.xlsx')
    spontan_inc_only.to_excel(save_dir + '/spontan_inc_only.xlsx')

def find_mini_xlsx_files(root_directory):
    """
    Search for .xlsx files containing 'connected' in filename within all subfolders
    """
    connected_files = []
    
    # Method 1: Using glob with recursive search
    pattern = os.path.join(root_directory, '**', '*connections*.xlsx')
    connected_files = glob.glob(pattern, recursive=True)
    
    return connected_files

def QC_mini_min_holding(mini_df, min_baseline = -300):
    '''
    goes through the mini_df
    if max val for all swps is < -300, excludes the fn
    '''
    mini_df = mini_df.reset_index(drop = True)

    indx_exclude = []
    for i, _ in enumerate(mini_df.min_vals_each_swp):
        list_max = ast.literal_eval(mini_df.max_vals_each_swp[i])

        if max(list_max) < min_baseline:
            indx_exclude.append(i)

    return mini_df.drop(index=indx_exclude).reset_index(drop=True)

def get_mini_df():
    '''
    collects all mini dfs
    cross compares and keeps only OPs with intrinsic data
    saves the results for statistical analysis
    '''
    # collecting and cross checking spontan / mini dataframes

    # collect event meta_dfs
    meta_mini, QC_df_all = funcs_results.collect_events_dfs('minis', QC = True)
    # remove cells where nothing to analyse and OP210615
    # meta_mini = meta_mini[meta_mini['OP'] != 'OP210615']
    meta_mini = meta_mini[meta_mini['swps_to_analyse'] != '[]'].reset_index(drop=True)
    QC_df_all = QC_df_all[QC_df_all['swps_to_analyse'] != '[]'].reset_index(drop=True)
    meta_mini = get_match_cell_id(meta_mini)
    QC_df_all = get_match_cell_id(QC_df_all)

    print(f'len mini_meta {len(meta_mini)}')
    print(f'len QC_df_all {len(QC_df_all)}')
    print(f'len QC_df_all cells unique {len(meta_mini.cell_IDs_match.unique())}')
    print(f'len QC_df_all cells unique {len(QC_df_all.cell_IDs_match.unique())}')

    set1 = set(QC_df_all.cell_IDs_match.unique()) - set(meta_mini.cell_IDs_match.unique())
    set2 = set(meta_mini.cell_IDs_match.unique()) - set(QC_df_all.cell_IDs_match.unique()) 

    if len(set1) > 0 or len(set2) >2:
        print('something wrong different cells in meta and QC')

    mini_all = pd.merge(meta_mini, QC_df_all, on = 'cell_IDs_match', \
                        suffixes = (None, '_y')).reset_index(drop=True)

    # load QC checked intrinsic data
    data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
 
    df_intr_complete = pd.read_excel(data_dir + 'adult_16min_hrs_inc_temporal.xlsx')
    df_intr_complete = get_match_cell_id(df_intr_complete)

    # with cell_ID from slice
    ops_exclude = set(meta_mini.OP.unique()) - set(df_intr_complete.OP.unique())
    print('removing ', str(len(ops_exclude)), 'OPs, missing from intr_df')
    
    mini_all = mini_all[~mini_all['OP'].isin(ops_exclude)].reset_index(drop = True)

    # remove unnecessary columns
    col_remove = list(mini_all.columns[mini_all.columns.str.contains('_y')]) # all doubles
    mini_all = mini_all.drop(columns = col_remove)
    mini_all = mini_all.loc[:, ~ mini_all.columns.str.contains('Unnamed')]

    mini_all = mini_all[mini_all.treatment != 'TTX']

    mini_all = QC_mini_min_holding(mini_all)
    save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
    mini_all.to_csv(save_dir + '/mini_QC_temporal_adult.csv')
    mini_all.to_excel(save_dir + '/mini_QC_temporal_adult.xlsx')


def get_counts_non_firing_cells():
    data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'

    df_repatch = pd.read_excel(data_dir+ 'repatch_data_temporal.xlsx')
    df_slice = pd.read_excel(data_dir+ 'slice_data_temporal.xlsx')
    df_slice_short_inc = pd.read_excel(data_dir + '10.07.25_slice_incubation_only_no_move.xlsx')
    df_slice_short_inc = df_slice_short_inc[df_slice_short_inc['hrs_after_OP'] < 35]
    df_slice_short_inc = df_slice_short_inc[df_slice_short_inc['hrs_incubation'] > 15.9]
    df_slice_short_inc.reset_index(drop = True, inplace = True)
    df_inc = df_slice_short_inc.loc[[i for i, sl in enumerate(df_slice_short_inc.slice) if len(sl) <= 2], :]

    df_dict = {'repatch': df_repatch,
            'slice': df_slice,
            'inc_only': df_inc}
    for data_type, df in df_dict.items():
        pl_intr.count_non_firing_cells_d2(df, data_dir, data_type)

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

        print(df.OP[i])
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

class Correlations:
    '''
    different versions of the correlation plots
    '''
    def __init__(self, df_slice, df_repatch):
        # remove rows where patient_age is unknown
        self.df_slice = df_slice[df_slice['patient_age'] != 'A']
        self.df_repatch = df_repatch[df_repatch['patient_age'] != 'A']
        self.scaling = 1.7
        self.ivs = {'patient_age': [12,
                        'Age (years)', 
                        [16.5, 26.5, 36.5, 46.5, 56.5, 66.5, 75],
                        ['12-21', '22-31', '32-42', '42-52', '53-62', '63-72','75']],
                    'hrs_after_OP': [7, 
                            'Time after tissue extraction (hrs)', 
                            [11.5, 21.5, 31.5, 41.5, 51.5],
                            ['7-16', '17-26', '27-36', '37-46', '47-57']],
                    'hrs_incubation': [0,
                                       'Time in incubstion solution (hrs)',
                                       [0, 5, 10, 15, 20, 25],
                                       ['0', '5', '10', '15', '20', '25']]}
        self.ivs_simple = {'patient_age': [[12, 32, 52, 76],
                'Age (years)', 
                [22, 42, 62],
                ['12-31', '32-51', '51-75']],
            'hrs_after_OP': [[7, 21, 41, 63], 
                    'Time after tissue extraction (hrs)', 
                    [15, 30, 50],
                    ['7-20', '21-40', '40-60']]}
        self.params = ['resting_potential', 'Rin', 'cap_adj', 'TH',
                'membra_time_constant_tau', 'rheo_ramp_c', 'sag',
                'max_depol', 'max_repol', 'Rs','AP_heigth', 'AP_halfwidth']
        self.dfs = {'slice': df_slice, 'repatch': df_repatch}
        self.colors = {'sliceCtrl': 'blue',
            'slicehigh K': 'green', 
            'repatchCtrl': 'orange',
            'repatchhigh K': 'red'}
        self.colors_2 = {'D1': 'darkgrey',
            'Ctrl': '#c875c4',
            'high K': '#4daf4a'}
        self.dict_for_plotting = pl_intr.dict_for_plotting()
        self.dest_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+ \
                    'paper_figs_collected_checked/plots/correlations/'
        os.makedirs(self.dest_dir, exist_ok = True)

    def boxplots_only(self):
        '''
        boxplots for each 10 units bin
        '''
        for iv, iv_items in self.ivs.items():
            for dv in self.params:
                fig, ax = plt.subplots(figsize = (9,6))
                for i, che in enumerate(self.dfs.items()):
                    name = che[0]
                    df  = che[1]
                    for t, treatment in enumerate(['Ctrl', 'high K']):
                        k = t + 2*i
                        df_tr = df[df['treatment'] == treatment]

                        df_tr = df_tr[pd.to_numeric(df_tr[iv], errors='coerce').notnull()]
                        df_tr = df_tr[pd.to_numeric(df_tr[dv], errors='coerce').notnull()]

                        max_iv = int(df_tr[iv].max())

                        data_plot, x_medians = [], []
                        for start in range(iv_items[0],  max_iv + 1, 10):
                            condition = (df_tr[iv] >= start) & (df_tr[iv] < start + 10)
                            x_medians.append(start + 3 + k)
                            data_plot.append(df_tr[dv][condition].tolist())
    
                        meds = [x/self.scaling for x in x_medians] # puts the boxplots closer to each other
                        ax.boxplot(data_plot, positions = meds, notch = True, patch_artist = True, \
                                    boxprops = dict(facecolor = self.colors[name+treatment], alpha = 0.35, color = self.colors[name+treatment],
                                                    linewidth = 2.3), \
                                    medianprops = dict(linewidth = 2.3, color = 'k'), \
                                    whiskerprops=dict(color = self.colors[name+treatment], linewidth=1.5), \
                                    capprops=dict(color = self.colors[name+treatment], linewidth=1.5), \
                                    label = name + ' ' + treatment)
                ax.set_xlabel(iv_items[1])
                ax.set_ylabel(self.dict_for_plotting[dv][1])
                ticks1 = [x / self.scaling for x in iv_items[2]]
                ax.set_xticks(ticks = ticks1, labels = iv_items[3], size = 20, rotation = 30)
                fig.legend(loc = 'upper right', fontsize = 'small', bbox_to_anchor=(1, 0.95))
                fig.suptitle(self.dict_for_plotting[dv][0])
                fig.savefig(self.dest_dir + 'boxplots/' + iv + dv + '.svg', \
                    format = 'svg', bbox_inches = 'tight', dpi = 500)
                plt.close()
    
    def boxplots_short(self):
        '''
        boxplots for rought overview
        '''
        for iv, iv_items in self.ivs_simple.items():
    
            for dv in self.params:
                fig, ax = plt.subplots(1,2, figsize = (9,6))
                for i, che in enumerate(self.dfs.items()):
                    name = che[0]
                    df  = che[1]

                    for t, treatment in enumerate(['D1', 'Ctrl', 'high K']):
                        if treatment == 'D1':
                            df_tr = df[df['day'] == 'D1']
                        else:
                            df_tr = df[(df['day'] == 'D2') & (df['treatment'] == treatment)]
                        df_tr = df_tr[pd.to_numeric(df_tr[iv], errors='coerce').notnull()]
                        df_tr = df_tr[pd.to_numeric(df_tr[dv], errors='coerce').notnull()]
                        
                        data_plot, x_medians = [], []
                        for start in iv_items[0]:
                            condition = (df_tr[iv] >= start) & (df_tr[iv] < start + 20)
                            if len(df_tr[condition]) < 0:
                                continue
                            x_medians.append(start + 10 + 1.5*t)
                            data_plot.append(df_tr[dv][condition].tolist())
    
                        ax[i].boxplot(data_plot, positions = x_medians, notch = True, patch_artist = True, \
                                    boxprops = dict(facecolor = self.colors_2[treatment], alpha = 0.35, color = self.colors_2[treatment],
                                                    linewidth = 2.3), \
                                    medianprops = dict(linewidth = 2.3, color = 'k'), \
                                    whiskerprops=dict(color = self.colors_2[treatment], linewidth=1.5), \
                                    capprops=dict(color = self.colors_2[treatment], linewidth=1.5), \
                                    label = treatment)
                        ax[i].set_xlabel(iv_items[1])
                        ax[i].set_ylabel(self.dict_for_plotting[dv][1])
                        ax[i].set_xticks(ticks = iv_items[2], labels = iv_items[3], size = 20, rotation = 30)
                fig.legend(loc = 'upper right', fontsize = 'small', bbox_to_anchor=(1, 0.95))
                fig.suptitle(self.dict_for_plotting[dv][0])
                dest_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/plots/correlations/'
                fig.savefig(dest_dir + 'boxplots/' + iv + dv + '.svg', \
                    format = 'svg', bbox_inches = 'tight', dpi = 500)
                plt.close()

    def boxplots_dots(self):
        '''
        boxplots for each 10 units // when <3 dots --> no boxplot
        '''

        for iv, iv_items in self.ivs.items():

            for dv in self.params:
                fig, ax = plt.subplots(figsize = (9,6))

                for i, che in enumerate(self.dfs.items()):
                    name = che[0]
                    df  = che[1]
                    for t, treatment in enumerate(['Ctrl', 'high K']):
                        k = t + 2*i
                        df_tr = df[df['treatment'] == treatment]

                        df_tr = df_tr[pd.to_numeric(df_tr[iv], errors='coerce').notnull()]
                        df_tr = df_tr[pd.to_numeric(df_tr[dv], errors='coerce').notnull()]
                        max_iv = int(df_tr[iv].max())

                        data_plot, x_medians = [], []
                        for start in range(iv_items[0],  max_iv + 1, 10):
                            condition = (df_tr[iv] >= start) & (df_tr[iv] < start + 10)
                            if len(df_tr[dv][condition]) < 3:
                                x_scatter = np.linspace(start + 3, start + 6, len(df_tr[dv][condition])) / self.scaling
                                ax.scatter(x_scatter, df_tr[dv][condition].tolist(),
                                           c = self.colors[name+treatment], alpha = 0.75)
                                continue
                            x_medians.append(start + 3 + k)
                            data_plot.append(df_tr[dv][condition].tolist())
   
                        meds = [x/self.scaling for x in x_medians]
                        ax.boxplot(data_plot, positions = meds, notch = True, patch_artist = True, \
                                # boxprops = dict(facecolor = colors[name+treatment], alpha = 0.75, color = 'black'), \
                                    boxprops = dict(facecolor = self.colors[name+treatment], alpha = 0.35, color = self.colors[name+treatment],
                                                    linewidth = 2.3), \
                                    medianprops = dict(linewidth = 2.3, color = 'k'), \
                                    whiskerprops=dict(color = self.colors[name+treatment], linewidth=1.5), \
                                    capprops=dict(color = self.colors[name+treatment], linewidth=1.5), \
                                    label = name + ' ' + treatment)
                ax.set_xlabel(iv_items[1])
                ax.set_ylabel(self.dict_for_plotting[dv][1])
                ticks1 = [x / self.scaling for x in iv_items[2]]
                ax.set_xticks(ticks = ticks1, labels = iv_items[3], size = 20, rotation = 30)
                fig.legend(loc = 'upper right', fontsize = 'small', bbox_to_anchor=(1, 0.93))
                fig.suptitle(self.dict_for_plotting[dv][0])

                fig.savefig(self.dest_dir + 'boxplots_dots/' + iv + dv + '.svg', \
                    format = 'svg', bbox_inches = 'tight', dpi = 500)
                plt.close()

    def sl_repatch_separate_scatter_median(self):
        '''
        separate scatterplots for slice and repatch
        medians for each 10 units bin
        '''
        for iv, iv_items in self.ivs.items():

            for dv in self.params:
                fig, ax = plt.subplots(1, 2,figsize = (16,6), sharex = True, sharey = True)

                for i, che in enumerate(self.dfs.items()):
                    name = che[0]
                    df  = che[1]
                    for t, treatment in enumerate(['Ctrl', 'high K']):
                        k = t + 2*i
                        df_tr = df[df['treatment'] == treatment]

                        df_tr = df_tr[pd.to_numeric(df_tr[iv], errors='coerce').notnull()]
                        df_tr = df_tr[pd.to_numeric(df_tr[dv], errors='coerce').notnull()]

                        ax[i].scatter(df_tr[iv], df_tr[dv], c = self.colors[name + treatment], \
                                      alpha = 0.45)

                        max_iv = int(df_tr[iv].max())

                        x_medians, medians = [], []
                        for start in range(iv_items[0],  max_iv + 1, 10):
                            condition = (df_tr[iv] >= start) & (df_tr[iv] < start + 10)
                            data_condition = df_tr[dv][condition].tolist()
                            if len(data_condition) < 1:
                                # not adding nans
                                continue
                            x_medians.append(start + 3 + k)
                            medians.append(np.nanmedian(data_condition))

                        ax[i].plot(x_medians, medians, self.colors[name + treatment], marker = 'D' , \
                         markersize = 7, markeredgewidth = 4, alpha = 0.8, label = name + ' ' + treatment)

                        ax[i].set_xlabel(iv_items[1])
                        ax[i].set_ylabel(self.dict_for_plotting[dv][1])
                        ax[i].set_xticks(ticks = iv_items[2], labels = iv_items[3], size = 20, rotation = 30)
                fig.legend(loc = 'upper right', fontsize = 'small', bbox_to_anchor=(1, 1))
                fig.suptitle(self.dict_for_plotting[dv][0])
                fig.savefig(self.dest_dir + 'slice_repatch_medians/' + iv + dv + '.svg', \
                    format = 'svg', bbox_inches = 'tight', dpi = 500)
                plt.close()

    def sl_repatch_separate_scatter_reg_lines(self):
        '''
        scatter with regression lines
        '''
        for iv, iv_items in self.ivs.items():

            for dv in self.params:
                fig, ax = plt.subplots(1, 2,figsize = (16,6), sharex = True, sharey = True)

                for i, (name, df) in enumerate(self.dfs.items()):
                    for t, treatment in enumerate(['Ctrl', 'high K']):
                        k = t + 2*i
                        df_tr = df[df['treatment'] == treatment]
                        
                        # removing non-numeric values
                        mask = (pd.to_numeric(df_tr[iv], errors='coerce').notnull()) & \
                            (pd.to_numeric(df_tr[dv], errors='coerce').notnull())
                        df_tr = df_tr.loc[mask]
                        x = df_tr[iv].to_numpy().astype(int)
                        y = df_tr[dv].to_numpy()

                        ax[i].scatter(df_tr[iv], df_tr[dv], c = self.colors[name + treatment], \
                                      alpha = 0.45)
                        
                        # regression lines
                        slope, intercept = np.polyfit(x, y, 1)
                        y_est = x * slope + intercept
                        y_err = x.std() * np.sqrt(1/len(x) + (x - x.mean())**2 / \
                                np.sum((x - x.mean())**2))
                    
                        ax[i].plot(x, slope * x + intercept, \
                               color = self.colors[name+treatment], linestyle = '--', linewidth = 2,
                               label = name + ' ' + treatment)
                        ax[i].fill_between(x, y_est - y_err, y_est + y_err, alpha = 0.2, \
                                           color = self.colors[name+treatment])

                        ax[i].set_xlabel(iv_items[1])
                        ax[i].set_ylabel(self.dict_for_plotting[dv][1])

                fig.legend(loc = 'upper right', fontsize = 'small', bbox_to_anchor=(1, 1))
                fig.suptitle(self.dict_for_plotting[dv][0])
                fig.savefig(self.dest_dir + 'slice_repatch_regression_line/' + iv + dv + '.svg', \
                    format = 'svg', bbox_inches = 'tight', dpi = 500)
                plt.close()
        
    def sl_repatch_separate_means_regression_lines(self):
        '''
        means for 10 unit bins and regression lines
        '''

        # for iv, iv_items in self.ivs.items():

        #     for dv in self.params:
        #         fig, ax = plt.subplots(1, 2,figsize = (16,6), sharex = True, sharey = True)

        #         for i, (name, df) in enumerate(self.dfs.items()):
       
        #             for t, treatment in enumerate(['Ctrl', 'high K']):
        #                 k = t + 2*i
        #                 df_tr = df[df['treatment'] == treatment]

        #                 df_plot = df_tr.groupby([iv])[dv].describe(exclude = 'nan')[["count", "mean", "std"]].reset_index()
        #                 mask2 = pd.to_numeric(df_plot['std'], errors='coerce').notnull()
        #                 df_plot = df_plot.loc[mask2].reset_index(drop = True)

        #                 df_plot["lower_ci"] = abs(df_plot["mean"] - 1.96*(df_plot["std"]/np.sqrt(df_plot["count"])))
        #                 df_plot["upper_ci"] = abs(df_plot["mean"] + 1.96*(df_plot["std"]/np.sqrt(df_plot["count"])))

        #                 df_plot = df_plot.sort_values(by=iv).reset_index(drop=True)
        #                 plt.plot(df_plot[iv], df_plot["mean"], label='Mean')
        #                 plt.plot(df_plot[iv], df_plot["lower_ci"],'-b', label='Upper CI')
        #                 plt.plot(df_plot[iv], df_plot["upper_ci"],'-r', label='Lower CI')
        #                 plt.fill_between(df_plot[iv], df_plot["lower_ci"], df_plot["upper_ci"],\
        #                                  alpha = 0.3, color = colors[name + treatment],
        #                                  interpolate = True)
        #                 plt.legend()
        #                 plt.show()



                        

        #                 # removing non-numeric values
        #                 mask = (pd.to_numeric(df_tr[iv], errors='coerce').notnull()) & \
        #                     (pd.to_numeric(df_tr[dv], errors='coerce').notnull())
        #                 df_tr = df_tr.loc[mask].reset_index(drop = True)


                    
                       
        #                 x = df_tr[iv].to_numpy().astype(int)
        #                 y = df_tr[dv].to_numpy()
                        
        #                 # regression lines
        #                 slope, intercept = np.polyfit(x, y, 1)
        #                 y_est = x * slope + intercept
        #                 y_err = x.std() * np.sqrt(1/len(x) + (x - x.mean())**2 / \
        #                         np.sum((x - x.mean())**2))
                    
        #                 ax[i].plot(x, slope * x + intercept, \
        #                        color = self.colors[name+treatment], linestyle = '--', linewidth = 2)
        #                 ax[i].fill_between(x, y_est - y_err, y_est + y_err, alpha = 0.2, \
        #                                    color = self.colors[name+treatment])

        #                 max_iv = int(x.max())
        #                 x_means, means= [], []
        #                 for start in range(iv_items[0],  max_iv + 1, 10):
        #                     condition = (df_tr[iv] >= start) & (df_tr[iv] < start + 10)
        #                     data_condition = df_tr[dv][condition].tolist()
        #                     if len(data_condition) < 1:
        #                         # not adding nans
        #                         continue
        #                     x_means.append(start + 4.3 + k)
        #                     means.append(np.nanmean(data_condition))

        #                 ax[i].scatter(x_means, means, s = 80, c = self.colors[name + treatment], \
        #                     marker = 'D', alpha = 0.8, label = name + ' ' + treatment)

        #                 ax[i].set_xlabel(iv_items[1])
        #                 ax[i].set_ylabel(self.dict_for_plotting[dv][1])
        #                 if iv == 'patient_age':
        #                     ax[i].set_xticks(ticks = iv_items[2][:-1], \
        #                     labels = iv_items[3][:-1], size = 20, rotation = 30)
        #                 else:
        #                     ax[i].set_xticks(ticks = iv_items[2], \
        #                     labels = iv_items[3], size = 20, rotation = 30)
                
        #         fig.legend(loc = 'upper right', fontsize = 'small', bbox_to_anchor=(1, 1))
        #         fig.suptitle(self.dict_for_plotting[dv][0])
        #         fig.savefig(self.dest_dir + 'slice_repatch_means_regression_line/' + iv + dv + '.svg', \
        #             format = 'svg', bbox_inches = 'tight', dpi = 500)
        #         plt.close()
                
    def sl_repatch_separate_means_confidence_ints(self):
        '''
        confidence intervals only
        not complete?
        '''

        for iv, iv_items in self.ivs.items():

            for dv in self.params:
                fig, ax = plt.subplots(3, 2,figsize = (16,6), sharex = True, sharey = True)

                for i, (name, df) in enumerate(self.dfs.items()):
       
                    for day in df.day:
                        for tr, treatment in enumerate(['Ctrl', 'high K']):
                            k = 0
                            # removing non-numeric values
                            df_day = df[df['day'] == day]
                            mask = (pd.to_numeric(df_day[iv], errors='coerce').notnull()) & \
                                (pd.to_numeric(df_day[dv], errors='coerce').notnull())
                            df_day = df_day.loc[mask]

                            if day == 'D2':
                                df_day = df_day[df_day['treatment'] == treatment]
                                k = 1 + tr

                            df_day = df_day.groupby([iv])[dv].describe()[["count", "mean", "std"]].reset_index()
                            df_day["lower_ci"] = df_day["mean"] - 1.96*(df_day["std"]/np.sqrt(df_day["count"]))
                            df_day["upper_ci"] = df_day["mean"] + 1.96*(df_day["std"]/np.sqrt(df_day["count"]))

                            df_day = df_day.sort_values(by = iv)
                            ax[k,i].plot(df_day[iv], df_day["mean"])
                            # ax[k,i].plot(df_day[iv], df_day["lower_ci"],'-b')
                            # ax[k,i].plot(df_day[iv], df_day["upper_ci"],'-r')
                            ax[k,i].fill_between(df_day[iv], df_day["mean"] - abs(df_day["lower_ci"]), df_day["mean"] + abs(df_day["upper_ci"]), alpha = 0.1)
                                

                ax[k,i].set_xlabel(iv_items[1])
                ax[2,i].set_ylabel(self.dict_for_plotting[dv][1])

                if iv == 'patient_age':
                    ax[2,i].set_xticks(ticks = iv_items[2][:-1], \
                    labels = iv_items[3][:-1], size = 20, rotation = 30)
                else:
                    ax[2,i].set_xticks(ticks = iv_items[2], \
                    labels = iv_items[3], size = 20, rotation = 30)
                    
                fig.legend(loc = 'upper right', fontsize = 'small', bbox_to_anchor=(1, 1))
                fig.suptitle(self.dict_for_plotting[dv][0])
            
                fig.savefig(self.dest_dir + 'slice_repatch_conf_int/' + iv + dv + '.svg', \
                    format = 'svg', bbox_inches = 'tight', dpi = 500)
                plt.close()
        
    def sl_repatch_separate_means_D1_all_ctrl(self):
        
        colors_D1_ctrl_highK = {'sliceD1': ['darkgrey', 'pre'],
            'sliceCtrl': ['#c875c4', 'slice post CTR'],
            'slicehigh K': ['#4daf4a', 'slice post HiK'], 
            'repatchD1': ['darkgrey','pre'],
            'repatchCtrl': ['darkblue', 'repatch post CTR'],
            'repatchhigh K': ['#ff7f00', 'repatch post HiK']}
        
        w_cm = 20
        h_cm = 7
        
        for iv, iv_items in self.ivs.items():

            for dv in self.params:
                fig, ax = plt.subplots(1, 2,figsize = (w_cm/2.54, h_cm / 2.54), sharex = True, sharey = True)

                for i, (name, df) in enumerate(self.dfs.items()):
                    df = df[pd.to_numeric(df[iv], errors='coerce').notnull()]
                    df = df[pd.to_numeric(df[dv], errors='coerce').notnull()]
       
                    for t, treatment in enumerate(['Ctrl', 'high K']):
                        k = t + 2*i

                        if t == 0:
                            df_day = df[df['day'] == 'D1']
                            ax[i].scatter(df_day[iv], df_day[dv], c = colors_D1_ctrl_highK[name + 'D1'][0], \
                                      alpha = 0.45)
                            
                            max_iv = int(df_day[iv].max())
                            day_x, day_medians =  [], []
                            for start in range(iv_items[0],  max_iv + 1, 10):
                                day_condition = (df_day[iv] >= start) & (df_day[iv] < start + 10)
                                day_data_condition = df_day[dv][day_condition].tolist()
                            
                                if len(day_data_condition) < 1:
                                    # not adding nans
                                    continue
                                day_x.append(start + 3 + k)
                                day_medians.append(np.nanmean(day_data_condition))

                            ax[i].plot(day_x, day_medians, colors_D1_ctrl_highK[name + 'D1'][0], marker = 'D' , \
                                    markersize = 7, markeredgewidth = 4, alpha = 0.8, \
                                        label = colors_D1_ctrl_highK[name + 'D1'][1])

                        
                        df_tr = df[(df['treatment'] == treatment) & (df['day'] == 'D2')]
                        ax[i].scatter(df_tr[iv], df_tr[dv], c = colors_D1_ctrl_highK[name + treatment][0], \
                                      alpha = 0.4)
            
                        max_iv = int(df_tr[iv].max())

                        x_medians, medians = [], []
                        for start in range(iv_items[0],  max_iv + 1, 10):
                            condition = (df_tr[iv] >= start) & (df_tr[iv] < start + 10)
                            data_condition = df_tr[dv][condition].tolist()
                            if len(data_condition) < 1:
                                # not adding nans
                                continue
                            x_medians.append(start + 3 + k)
                            medians.append(np.nanmean(data_condition))

                        ax[i].plot(x_medians, medians, colors_D1_ctrl_highK[name + treatment][0], marker = 'D' , \
                         markersize = 7, markeredgewidth = 4, alpha = 0.9, label = name + ' ' + treatment)

                        # ax[i].set_xlabel(iv_items[1])
                        ax[i].set_ylabel(self.dict_for_plotting[dv][0] + ' (' + self.dict_for_plotting[dv][1] + ')')
                        ax[i].set_xticks(ticks = iv_items[2], labels = iv_items[3], size = 20, rotation = 30)
                
                ax[i].set_xlabel(iv_items[1])       
                # fig.legend(loc = 'upper right', fontsize = 'small', bbox_to_anchor=(1, 1))
                # fig.suptitle(self.dict_for_plotting[dv][0])
                os.makedirs(self.dest_dir + 'slice_repatch_means_D1_ctrl_highK/', exist_ok = True)
                fig.savefig(self.dest_dir + 'slice_repatch_means_D1_ctrl_highK/' + iv + dv + '.svg', \
                    format = 'svg', bbox_inches = 'tight', dpi = 500)
                plt.close()
     
def _plot_time_dep(min_holding = -400):

    '''
    plot time dependency for 5 parameteres
    '''
    # plots over time
    param_list = ['resting_potential', 'Rin', 'cap_adj', 'TH',
                  'membra_time_constant_tau', 'Rheobase', 'rheo_ramp_c']

    adult_df_ctrl = load_intr_data('Ctrl')
    df_qc_hold = get_lab_book_info(adult_df_ctrl)
    df_qc_hold = get_match_cell_id(df_qc_hold)
    adult_df_ctrl = get_match_cell_id(adult_df_ctrl)

    # join QC and adult_df based on unique cell_ID
    adult_df_ctrl = pd.merge(adult_df_ctrl, df_qc_hold, on = 'cell_IDs_match', \
                        suffixes = (None, '_y')).reset_index(drop=True)

    adult_df_ctrl['hold_start'] = pd.to_numeric(adult_df_ctrl['hold_start'], errors='coerce')
    adult_df_ctrl = adult_df_ctrl[adult_df_ctrl.hold_start > min_holding].reset_index(drop = True)

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
    adult_df_hik = adult_df_hik[adult_df_hik.hold_start > min_holding].reset_index(drop = True)

    plot_param_over_time(adult_df_hik, param_list, 'high K')

    # mean plots
    mea_figure()
    mea_figure_normalized()

def add_treatment_r():
    '''
    adds treatment to the dataframe
    '''
    data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
        'paper_figs_collected_checked/data/'

    df_repatch = pd.read_excel(data_dir + 'repatch_data_temporal.xlsx')
    df_slice = pd.read_excel(data_dir + 'slice_data_temporal.xlsx')
    df_slice = df_slice[df_slice.treatment != 'fixed'].reset_index(drop = True)

    if 'treatment_r' in df_slice.columns:
        df_slice.drop(columns = ['treatment_r'], inplace = True)

    tr_r = []
    for i, hrs in enumerate(df_slice.hrs_incubation):
        tr = df_slice.treatment[i]
        if hrs == 0:
            tr_r.append('BL')
        elif hrs > 0:
            tr_r.append(tr)

    if len(tr_r) == len(df_slice):
        df_slice.insert(len(df_slice.columns), 'treatment_r', tr_r)
        df_slice.to_excel(data_dir + 'slice_data_temporal.xlsx')
        df_slice.to_csv(data_dir + 'slice_data_temporal.csv')

    if 'treatment_r' in df_repatch.columns:
        df_repatch.drop(columns = ['treatment_r'], inplace = True)

    tr_r = []
    for i, hrs in enumerate(df_repatch.hrs_incubation):
        tr = df_repatch.treatment[i]
        if hrs == 0:
            tr_r.append('BL')
        elif hrs > 0:
            tr_r.append(tr)

    if len(tr_r) == len(df_repatch):
        df_repatch.insert(len(df_repatch.columns), 'treatment_r', tr_r)
        df_repatch.to_excel(data_dir + 'repatch_data_temporal.xlsx')
        df_repatch.to_csv(data_dir + 'repatch_data_temporal.csv')

def get_meta_data_summary(meta_df):
    # meta_df = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/meta_data.xlsx')
    # filtered_meta_df = meta_df[meta_df['OP'].isin(used_ops)]
    counts = meta_df['sex'].value_counts()
    print(counts)
    num_male = counts.get('male', 0)
    num_female = counts.get('female', 0)
    print(f"male: {num_male}, female: {num_female}")
    print(f'median age {meta_df.age.median()}')
    print(f'min age {meta_df.age.min()}')
    print(f'max age {meta_df.age.max()}')

# if __name__ == "__main__":
#     main()


# %%
## FOR JUV DATA

# df_intr_complete = funcs_results.collect_intrinsic_df() # for most updated version
# df_intr_complete = combine_columns(df_intr_complete)

# # df_intr = pl_intr.get_column_RMPs_from_char(df_intr_complete)

def get_QC_data(df):
    #df['resting_potential'] = pd.to_numeric(df['resting_potential'],errors='coerce')
    mask = (df['Rs'] < 50) & \
        (df['resting_potential'] < -45 ) & (df['resting_potential'] > -90) & \
        (df['capacitance'] < 1000) & (df['membra_time_constant_tau'] < 1000)
                # (df['membra_time_constant_tau'] > -25 ) &\
                #     (df['capacitance'] < 900) &  (df['capacitance'] > 10)
    #(df['max_repol'] > -100)
    #(df['AP_heigth'] > 30) & \
    # (df['Rheobase'] > 0)

    #include filter for TH ?
    df = df.loc[mask, :]
    return df

def filter_adult_hrs_incubation_data(df_intr_props, min_age, hrs_inc, max_age = 151, max_hrs_incubation = 50):
    adult_df = pl_intr.patient_age_to_float(df_intr_props, min_age, max_age)
    adult_df = pl_intr.change_to_numeric(adult_df)
    adult_df = get_QC_data(adult_df)
    adult_df = pl_intr.filter_on_hrs_incubation(adult_df, hrs_inc, max_hrs_incubation)
    return adult_df

def QC_mini_min_holding(mini_df, min_baseline = -300):
    '''
    goes through the mini_df
    if max val for all swps is < -300, excludes the fn
    '''
    mini_df = mini_df.reset_index(drop = True)

    indx_exclude = []
    for i, _ in enumerate(mini_df.min_vals_each_swp):
        list_max = ast.literal_eval(mini_df.max_vals_each_swp[i])

        if max(list_max) < min_baseline:
            indx_exclude.append(i)

    return mini_df.drop(index=indx_exclude).reset_index(drop=True)




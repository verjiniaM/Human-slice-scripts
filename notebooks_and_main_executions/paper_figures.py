import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ephys_analysis.funcs_sorting as sort
from collections import Counter
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
from ephys_analysis.funcs_for_results_tables import collect_intrinsic_df

#%%

def main():
    '''
    collect the dault_df repatch for statistical analysis
    '''
    adult_df = load_intr_data('all')
    adult_df = adult_df[adult_df['treatment'] != 'TTX']
    adult_df = combine_columns(adult_df)
    df_qc_hold = get_lab_book_info(adult_df)
    df_qc_hold = get_match_cell_ID(df_qc_hold)
    adult_df = get_match_cell_ID(adult_df)

    if len(adult_df.cell_IDs_match) < len(adult_df):
        print('smth wrong with unique cell_IDs adult')

    if len(df_qc_hold.cell_IDs_match) < len(df_qc_hold):
        print('smth wrong with unique cell_IDs QC')

    # remove OPs where no recording on second day
    op_d2 = adult_df.OP[adult_df.day.str.contains('D2', case=False).tolist()].unique()
    op_remove = list(set(adult_df.OP.unique()) - set(op_d2))

    adult_df = adult_df[~adult_df['OP'].isin(op_remove)] # some data for D2 available

    missing_cell_ids = []
    _ = [missing_cell_ids.append(cell_id) for cell_id in adult_df['cell_IDs_match'].tolist() \
    if cell_id not in df_qc_hold['cell_IDs_match'].tolist()]

    # remove the adult_df
    adult_df = adult_df[~adult_df['cell_IDs_match'].isin(missing_cell_ids)].reset_index(drop=True)

    # join QC and adult_df based on unique cell_ID
    adult_df_qc = pd.merge(adult_df, df_qc_hold, on = 'cell_IDs_match', \
                        suffixes = (None, '_y')).reset_index(drop=True)

    # check if join operation worked
    if False in list((adult_df_qc.cell_IDs_match) == (adult_df.cell_IDs_match)):
        print('error')

    if False not in list((adult_df_qc.slice) == (adult_df_qc.slice_y)):
        col_remove = list(adult_df_qc.columns[adult_df_qc.columns.str.contains('_y')])
        col_remove.remove('cell_IDs_match_re_y')

        adult_df_qc = adult_df_qc.drop(columns = col_remove)


        # adult_df_qc.to_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/paper_figs_collected_checked/adult_df_complete.xlsx')

    adult_df_repatch = pl_intr.get_repatch_df(adult_df_qc)
    adult_df_slice = adult_df_qc[adult_df_qc.repatch == 'no']

    # save for stats
    # adult_repatch = pl_intr.get_repatch_df(adult_df)
    # adult_repatch.to_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+ \
    #                      'data/summary_data_tables/intrinsic_properties/2025-02-05_repatch.csv')
    return adult_df_repatch, adult_df_slice

def main_old():
    '''
    plot time dependency for 5 parameteres
    '''
    adult_df_ctrl = load_intr_data('Ctrl')
    plot_param_over_time(adult_df_ctrl, 'resting_potential')
    plot_param_over_time(adult_df_ctrl, 'Rin')
    plot_param_over_time(adult_df_ctrl, 'capacitance')
    plot_param_over_time(adult_df_ctrl, 'TH')
    plot_param_over_time(adult_df_ctrl, 'membra_time_constant_tau')
    mea_figure()
    mea_figure_normalized()

def load_intr_data(treatment):
    '''
    collects the intrinsic dfs 

    type (str): 'all', 'Ctrl', 'high K'
    '''
    # intrinsic df to OP241120
    df_intr_complete = collect_intrinsic_df() # for most updated version

    # df_intr = pl_intr.get_column_RMPs_from_char(df_intr_complete)

    adult_df = pl_intr.filter_adult_hrs_incubation_data(df_intr_complete, min_age = 12, \
                                                        hrs_inc = 16, max_age = 151) # + QC criteria

    if treatment == 'all':
        return adult_df.reset_index(drop=True)

    return adult_df[adult_df.treatment == treatment].reset_index(drop=True)

def plot_param_over_time(adult_data, param):
    '''
    plots parameter over time
    '''
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

    try:
        plt.show()
    except (BrokenPipeError, IOError):
        pass
    fig.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/time_dependencies/example.png')

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
    plt.show()

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
    plt.show()

def get_lab_book_info(df, human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    '''
    collects info from the lab books
    '''
    patcher_dict = {
        'Verji': 'data_verji/',
        'Rosie': 'data_rosie/'
        }
    for patcher, dir_p in patcher_dict.items():
        patcher_df = df[df.patcher == patcher]

        for i, op in enumerate(patcher_df.OP.unique()):

            lab_book = sort.get_lab_book(human_dir + dir_p + op + '/')[0]

            if 'sl' not in lab_book.columns:
                print('fix lab book for', op)
            chans = lab_book.channels[1:].tolist()
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
            slice_names = slice_names[:len(chans)]

            if len(slice_names) < len(chans):
                slice_names.extend([float('nan')] * (len(chans) - len(slice_names)))

            hold_df = pd.DataFrame({'patcher': patcher, 'OP': op, 'slice':slice_names,
                'cell_ch':chans, 
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

def get_match_cell_ID(df):
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
        if isinstance(chan, float) or isinstance(chan, int) \
        or isinstance(chan, np.int64):
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
        return 'no column for unique cell IDs'
    counts = Counter(df.cell_IDs_match.tolist())
    repeats = [item for item, count in counts.items() if count > 1]

    if len(repeats) > 0:
        print('repeats adults_df', repeats)

if __name__ == "__main__":
    main()


# see the ones that are not in the adult
# match on them
# check the repatch ones



# remove OPs where no recording on second day
op_d2 = adult_df.OP[adult_df.day.str.contains('D2', case=False).tolist()].unique()
op_remove = list(set(adult_df.OP.unique()) - set(op_d2))

adult_df = adult_df[~adult_df['OP'].isin(op_remove)] # some data for D2 available

missing_cell_ids = []
_ = [missing_cell_ids.append(cell_id) for cell_id in adult_df['cell_IDs_match'].tolist() \
 if cell_id not in df_qc_hold['cell_IDs_match'].tolist()]

# remove the adult_df
adult_df = adult_df[~adult_df['cell_IDs_match'].isin(missing_cell_ids)].reset_index(drop=True)

# join QC and adult_df based on unique cell_ID
adult_df_qc = pd.merge(adult_df, df_qc_hold, on = 'cell_IDs_match', \
                       suffixes = (None, '_y')).reset_index(drop=True)

# check if join operation worked
if False in list((adult_df_qc.cell_IDs_match) == (adult_df.cell_IDs_match)):
    print('error')

if False not in list((adult_df_qc.slice) == (adult_df_qc.slice_y)):
    col_remove = list(adult_df_qc.columns[adult_df_qc.columns.str.contains('_y')])
    col_remove.remove('cell_IDs_match_re_y')

    adult_df_qc = adult_df_qc.drop(columns = col_remove)


    # adult_df_qc.to_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/paper_figs_collected_checked/adult_df_complete.xlsx')

adult_df_repatch = pl_intr.get_repatch_df(adult_df_qc)
adult_df_slice = adult_df_qc[adult_df_qc.repatch == 'no']
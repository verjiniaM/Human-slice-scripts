import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
from ephys_analysis.funcs_for_results_tables import collect_intrinsic_df

def main():
    '''
    plot time dependency for 5 parameteres
    '''
    adult_df_ctrl = load_intr_data()
    plot_param_over_time(adult_df_ctrl, 'resting_potential')
    plot_param_over_time(adult_df_ctrl, 'Rin')
    plot_param_over_time(adult_df_ctrl, 'capacitance')
    plot_param_over_time(adult_df_ctrl, 'TH')
    plot_param_over_time(adult_df_ctrl, 'membra_time_constant_tau')
    mea_figure()
    mea_figure_normalized()

def load_intr_data():
    '''
    collects the intrinsic dfs 
    '''
    # intrinsic df to OP241120
    df_intr_complete = collect_intrinsic_df() # for most unpdated version

    # df_intr = pl_intr.get_column_RMPs_from_char(df_intr_complete)

    adult_df = pl_intr.filter_adult_hrs_incubation_data(df_intr_complete, min_age = 12, \
                                                        hrs_inc = 16, max_age = 151) # + QC criteria
    adult_ctrl = adult_df[adult_df.treatment == 'Ctrl']

    return adult_ctrl

def plot_param_over_time(adult_data, param):
    '''
    plots parameter over time
    '''
    adult_repatch = pl_intr.get_repatch_df(adult_data)
    adult_slice =  adult_data.loc[adult_data['repatch'] == 'no']
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
        end = start + 10
        median_repatch = adult_repatch[(adult_repatch['hrs_after_OP'] >= start) & \
                                      (adult_repatch['hrs_after_OP'] < end)][param].median()
        median_slice = adult_slice[(adult_slice['hrs_after_OP'] >= start) & \
                                    (adult_slice['hrs_after_OP'] < end)][param].median()

        ax.plot([start], [median_repatch], 'orange', marker = 'D' , \
            markersize = 7, markeredgewidth = 4, alpha = 0.8)
        ax.plot([start], [median_slice], 'green', marker = 'D' , \
            markersize = 7, markeredgewidth = 4, alpha = 0.8)

    ax.set_xlabel('Hours after OP')
    ax.set_ylabel(dict_for_plotting[param][1])
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', fontsize='small')
    fig.suptitle(dict_for_plotting[param][0])

    try:
        plt.show()
    except (BrokenPipeError, IOError):
        pass

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

if __name__ == "__main__":
    main()

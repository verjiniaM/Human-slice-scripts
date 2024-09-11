import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def QC_filter_for_plotting(results_df, holding, temp, K_concentr = 8, min_age = 0):
    '''
    holding [str] - 'yes' or 'no'
    temp [str] - 'yes' or 'no'
    '''
    results_df_plot = results_df[(results_df['holding_minus_70_y_o_n'] == holding) &\
        (results_df['K_concentration'] == K_concentr) &\
            (results_df['recording_in'] != 'puff high K') &\
                #(results_df['Average amplitude (pA)'] > min_avg_amp) &\ 
                    (results_df['temperature'] == temp) &
                    (results_df['patient_age'] > min_age)]

    amp_non_negative = [-i for i in results_df_plot['Average amplitude (pA)']]
    results_df_plot.insert(20, 'Average amp (positive)', amp_non_negative)
    return results_df_plot

def QC_RMP_Ctrl(df, max_allowed_RMP_Ctrl):
    '''
    checks that the RMP in Ctrl condition is not above max_allowed_RMP_Ctrl
    exludes cells from both conditions, where it is
    '''

    #
    not_repeated_cells = []
    for cell in df['cell_ID'].unique():
        #print(len(df[df['cell_ID'] == cell]))
        if len(df[df['cell_ID'] == cell]) == 1: # or len(df[df['cell_ID'] == cell]) == 3:
            not_repeated_cells.append(cell)
    for cell in not_repeated_cells:
        df = df.drop(df.index[df['cell_ID'] == cell])     
    df.reset_index(inplace = True, drop = True)

    cell_ID_keep = df['cell_ID'][(df['recording_in'] == 'Ctrl') & (df['resting_potential'] < max_allowed_RMP_Ctrl)].to_list()
    cells_to_delete = list(set(df['cell_ID'].unique().tolist()) - set(cell_ID_keep))
    for cell in cells_to_delete:
        df = df.drop(df.index[df['cell_ID'] == cell])     
    df.reset_index(inplace = True, drop = True)

    return df

#plotting funcs
#mpl.rcParams - for all parameter settings
# @mpl.rc_context({'axes.labelsize': 17, \
#     'axes.spines.right': False, \
#         'axes.spines.top': False, \
#             'axes.titlesize': 15, \
#                 'xtick.labelsize': 15, \
#                      'ytick.labelsize': 15, \
#                         'figure.titlesize': 20})
def plot_ (df, title, params = ['Average amp (positive)', 'Average interevent interval (ms)']):
    OP_colors = ['#dede00', '#ff7f00', '#4daf4a', '#984ea3', 'violet']
    op_color_dict = {'OP230914':'#dede00', 'OP231005':'#ff7f00', 'OP231109':'#4daf4a', 'OP231123':'#984ea3', 'OP231130':'violet'}
    df = df.sort_values(by = ['recording_in', 'cell_ID']).reset_index(drop= True)
    fig, ax = plt.subplots(len(params),1, sharex = False, figsize=(6,10))
    for p, param in enumerate(params):
        x_vals = []
        for i, rec_solution in enumerate(sorted(df.recording_in.unique())):
            df_plot = df[df['recording_in'] == rec_solution].reset_index(drop= True) 
            x = np.linspace(2*i, 1 + 2*i, len(df_plot))
            ax[p].scatter(x, df_plot[param], alpha = 0.6, label = 'Ctrl')
            ax[p].plot([0.4 + 2*i, 0.6 + 2*i], [np.nanmean(df_plot[param]), np.nanmean(df_plot[param])], c = 'k')
            ax[p].text(0.4 + 2*i, np.nanmean(df_plot[param]) + p +0.03,  str(round(np.nanmean(df_plot[param]),2)), size = 15, c = 'k', zorder = 10)
            #ax[p].set_title(param)

            x_vals.append(x)

            for j, OP in enumerate(sorted(df_plot.OP.unique())):
                indx = df_plot[df_plot['OP'] == OP].index
                x_op = x[indx]
                y_op = df_plot[param][indx]
                ax[p].scatter(x_op, y_op, c = op_color_dict[OP], s = 60, zorder = 5, label = OP)
                
        cell_IDs = df['cell_ID'][df['recording_in'] == 'Ctrl'].values
        for c, cell in enumerate(cell_IDs):
            #indx = index_s[c]
            #x_K = results_df_plot['x'].loc[results_df_plot['cell_ID'] == cell].tolist()[0]
            # if len(x_vals[1]) <= c :
            #     continue
            x = [x_vals[0][c], x_vals[1][c]]
            y = [df[param][df['recording_in'] == 'Ctrl'].tolist()[c], df[param][df['recording_in'] == 'high K'].tolist()[c]]
            #op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
            ax[p].plot(x, y, '-', color = 'darkgrey', alpha = 0.5, linewidth = 2, zorder = 1)
    
    ax[0].set_ylabel('Amplitude (pA)')
    ax[1].set_ylabel('IEI (ms)')
    ax[2].set_ylabel('RMP (mV)')
    # ax[0].set_ylabel('Input resistance (MΩ)')
    # ax[1].set_ylabel('Input resistance (MΩ)')

    #ax[2].set_title('Resting membrane potential')
    ax[0].set_xticks(ticks = [0.5, 2.5], labels = ['Ctrl', 'high K'],)
    ax[1].set_xticks(ticks = [0.5, 2.5], labels = ['Ctrl', 'high K'],)
    ax[2].set_xticks(ticks = [0.5, 2.5], labels = ['Ctrl', 'high K'],)

    fig.suptitle(title)
    fig.tight_layout()
    fig.patch.set_facecolor('white')
    #fig.legend()

    plt.show()

# @mpl.rc_context({'axes.labelsize': 17, \
#     'axes.spines.right': False, \
#         'axes.spines.top': False, \
#             'axes.titlesize': 15, \
#                 'xtick.labelsize': 15, \
#                      'ytick.labelsize': 15, \
#                         'figure.titlesize': 20})
def sns_plot_MEA_data(df, title):
    colors = ['darkblue','#4daf4a', 'violet']
    customPalette = sns.set_palette(sns.color_palette(colors))

    fig, ax1 = plt.subplots(1,1, sharex = False, figsize=(8,5))
    sns.lineplot(
        data = df, x = 'Condition', y = "value",
        hue="OP", palette = customPalette , ax = ax1)
    sns.scatterplot(
        data = df, x = 'Condition', y = "value", 
        hue="OP", palette = customPalette , ax = ax1)

    ax1.set_xticks(ax1.get_xticks(), df['Condition'].unique(), rotation=30)
    
    ax1.set_xlabel('')
    ax1.set_ylabel('Network Activity \n(spikes\electrode\second')

    fig.suptitle(title)
    fig.tight_layout()
    fig.patch.set_facecolor('white')
    #fig.legend()
    fig.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.show()


def plot_from_full_results_table(results_df_plot, title):
    median_amp_no_hold_no_temp_Ctrl = results_df_plot['Average amp (positive)'][results_df_plot['recording_in'] == 'Ctrl'].values
    median_amp_no_hold_no_temp_highK = results_df_plot['Average amp (positive)'][results_df_plot['recording_in'] == 'high K'].values
    freq_no_hold_no_temp_Ctrl = results_df_plot['Average interevent interval (ms)'][results_df_plot['recording_in'] == 'Ctrl'].values
    freq_no_hold_no_temp_highK = results_df_plot['Average interevent interval (ms)'][results_df_plot['recording_in'] == 'high K'].values

    fig, ax = plt.subplots(2,1, sharex = True, figsize=(12,10))
    x1 = np.linspace(0, 1,len(median_amp_no_hold_no_temp_Ctrl))
    x2 = np.linspace(2, 3,len(median_amp_no_hold_no_temp_highK))
    ax[0].scatter(x1, median_amp_no_hold_no_temp_Ctrl, alpha = 0.7, label = 'Ctrl')
    ax[0].scatter(x2, median_amp_no_hold_no_temp_highK, alpha = 0.7, label = 'high K')
    ax[0].plot([0.4, 0.6], [np.nanmean(median_amp_no_hold_no_temp_Ctrl), np.nanmean(median_amp_no_hold_no_temp_Ctrl)], c = 'k')
    ax[0].plot([2.4, 2.6], [np.nanmean(median_amp_no_hold_no_temp_highK), np.nanmean(median_amp_no_hold_no_temp_highK)], c = 'k')

    ax[0].text(0.5, np.nanmean(median_amp_no_hold_no_temp_Ctrl) + 0.07,  str(round(np.nanmean(median_amp_no_hold_no_temp_Ctrl),2)), size = 15, c = 'k')
    ax[0].text(2.5, np.nanmean(median_amp_no_hold_no_temp_highK) + 0.07,  str(round(np.nanmean(median_amp_no_hold_no_temp_highK),2)), size = 15, c = 'k')
    ax[0].set_title('Average EPSP amplitude for cell per condition, no APs')

    cell_IDs = results_df_plot['cell_ID'][results_df_plot['recording_in'] == 'Ctrl'].values
    for c, cell in enumerate(cell_IDs):
        #indx = index_s[c]
        #x_K = results_df_plot['x'].loc[results_df_plot['cell_ID'] == cell].tolist()[0]
        x = [x1[c],x2[c]]
        y = [median_amp_no_hold_no_temp_Ctrl[c], median_amp_no_hold_no_temp_highK[c]]
        #op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
        ax[0].plot(x, y, '-', color = 'darkgrey', alpha = 0.5, linewidth = 2, zorder = 1)

    x3 = np.linspace(0, 1,len(freq_no_hold_no_temp_Ctrl))
    x4 = np.linspace(2, 3,len(freq_no_hold_no_temp_highK))
    ax[1].scatter(x3, freq_no_hold_no_temp_Ctrl, alpha = 0.7, label = 'Ctrl')
    ax[1].scatter(x4, freq_no_hold_no_temp_highK, alpha = 0.7, label = 'high K')
    ax[1].plot([0.4, 0.6], [np.nanmean(freq_no_hold_no_temp_Ctrl), np.nanmean(freq_no_hold_no_temp_Ctrl)], c = 'k')
    ax[1].plot([2.4, 2.6], [np.nanmean(freq_no_hold_no_temp_highK), np.nanmean(freq_no_hold_no_temp_highK)], c = 'k')

    ax[1].text(0.5, np.nanmean(freq_no_hold_no_temp_Ctrl) + 0.1,  str(round(np.nanmean(freq_no_hold_no_temp_Ctrl),2)), size = 15, c = 'k')
    ax[1].text(2.5, np.nanmean(freq_no_hold_no_temp_highK) + 0.1,  str(round(np.nanmean(freq_no_hold_no_temp_highK),2)), size = 15, c = 'k')
    ax[1].set_title('Average interevent interval (ms), no APs')

    cell_IDs = results_df_plot['cell_ID'][results_df_plot['recording_in'] == 'Ctrl'].values
    for c, cell in enumerate(cell_IDs):
        #indx = index_s[c]
        #x_K = results_df_plot['x'].loc[results_df_plot['cell_ID'] == cell].tolist()[0]
        x = [x3[c],x4[c]]
        y = [freq_no_hold_no_temp_Ctrl[c], freq_no_hold_no_temp_highK[c]]
        #op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
        ax[1].plot(x, y, '-', color = 'darkgrey', alpha = 0.5, linewidth = 2, zorder = 1)

    ax[0].set_ylabel('Amplitude (pA)')
    ax[1].set_ylabel('IEI (ms)')

    ax[1].set_xticks(ticks = [0.5, 2.5], labels = ['Ctrl', 'high K'])

    fig.suptitle(title)
    fig.tight_layout()
    fig.patch.set_facecolor('white')
    fig.legend()

    plt.show()

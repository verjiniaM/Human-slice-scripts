# TEMPORARY SCRIPT TO CHECK THE s
# and the failure rate in firing of the repatch dataset
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ephys_analysis.funcs_plot_paper_figures as paper_plot

data_dir_single = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
     'paper_figs_collected_checked/data/'
df_repatch = pd.read_excel(data_dir_single + 'repatch_data_temporal.xlsx')
df_slice = pd.read_excel(data_dir_single + 'slice_data_temporal.xlsx')

data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/connectivity/'
df_repatch_cons_all = pd.read_excel(f'{data_dir}/repatched_connections_looser_QC.xlsx')
df_repatch_cons = df_repatch_cons_all[~df_repatch_cons_all['comment'].str.contains('appearing|disappearing', case=False, na=False)].reset_index(drop = True)

df_connections_QC = pd.read_excel(data_dir + 'all_QCed_connections.xlsx')
df_connections_QC = df_connections_QC[df_connections_QC.repatch == 'no'].reset_index(drop = True)


# ### Check for which df we are doing it
# df_cons = df_connections_QC
# df_intr = df_slice
# # check if all connections have the same length
# for i, con_name in enumerate(df_cons.connection_ID):
#     # if len(con_name) != len(df_cons.connection_ID[0]):
#     if len(con_name) != 11 and len(con_name) != 14:
#         print(i, con_name)

# not_found_cells, manual_check = [], []
# cell_IDs, days, treatments, Rs = [], [], [], []
# for i, con in enumerate(df_cons.connection_ID):
#     pre = int(con[-1])
#     post = int(con[-3])

#     slic = df_cons.slice[i]
#     day = df_cons.day[i]
#     tr = df_cons.treatment[i]

#     for cell in [pre, post]:
#         df = df_intr[(df_intr.cell_ch == cell) &\
#             (df_intr.day == day) &\
#             (df_intr.treatment == tr) &\
#                 (df_intr.slice == slic) &\
#                     (df_intr.cell_ID == f'{con[:-4]}c{cell}')]
#         #            (df_repatch.repatch == 'yes')]

#         if len(df) < 1:
#             not_found_cells.append(f'{con[:-4]}c{cell}')
#         elif len(df) > 1:
#             manual_check.append(f'{con[:-4]}c{cell}')
#         else:
#             cell_IDs.append(df.cell_ID.values[0])
#             days.append(df.day.values[0])
#             treatments.append(df.treatment.values[0])
#             Rs.append(df.Rs.values[0])

# print(len(not_found_cells))
# print(len(manual_check))

# df_Rs_check = pd.DataFrame({'cell_IDs': cell_IDs,
#                             'day': days,
#                             'treatment': treatments,
#                             'Rs': Rs})


# # PLOT
# order = ["Ctrl D1", "Ctrl D2", "high K D1", "high K D2"]
# plot_df = df_Rs_check
# plot_df["group"] = plot_df.apply(
#     lambda r: f"{r['treatment']} {r['day']}",
#     axis=1
# )
# plot_df["group"] = pd.Categorical(plot_df["group"], categories=order, ordered=True)

# fig, ax = plt.subplots(figsize=(7, 4))
# sns.stripplot(
#     data=plot_df, x="group", y="Rs", order=order,
#     jitter=0.18, size=6, alpha=0.75, ax=ax
# )

# means = plot_df.groupby("group", observed=True)["Rs"].mean()
# for i, g in enumerate(order):
#     if g in means.index and pd.notna(means[g]):
#         ax.scatter(i, means[g], marker="_", s=500, color="black", zorder=5)

# ax.set_xlabel("")
# ax.set_ylabel("Rs")
# ax.set_title("Rs by day and treatment")
# plt.tight_layout()
# plt.show()


# #### EXISTING SCRIPTS FOR PLOTTIG Rs
# paper_plot.intr_params_fig3(df_repatch, save_dir)

# importlib.reload(paper_plot)
# df_slice = pd.read_excel(data_dir_single + 'slice_data_temporal.xlsx')
# paper_plot.intr_params_fig2(df_slice, 'other', save_dir)

#%%
## FIRING FAILURE RATE
import PyQt5
import ephys_analysis.stimulation_windows_ms as stim_win
import ephys_analysis.funcs_con_screen as con_param
import ephys_analysis.funcs_sorting as sort

def plot_post_cell(con_screen_file, pre_cell_chan, post_cell_chan):
    
    end_fn = con_screen_file.rfind('/') + 1
    dir_connect = sort.make_dir_if_not_existing(con_screen_file[:end_fn] + 'plots/',  'connectivity_plots')

    stim_window = stim_win.stim_window_con_screen

    pre_sig  = con_param.presynaptic_screen_ALLOW_ALL(con_screen_file, pre_cell_chan)
    post_sig, vm0_post = con_param.postsynaptic_screen_ALLOW_ALL(con_screen_file, post_cell_chan, [])
    mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = con_param.get_analysis_window(pre_sig, post_sig)
    sampl_rate = 200_000

    fig, ax = plt.subplots(1, 2, figsize=(20, 10))

    xvals = preAPs[0][0]-100
    for i in range(len(post_sig[0])):
        ax[0].plot(post_sig[xvals:xvals+1000, i] - i * 4, color ='grey', lw = 0.5)
        ax[0].hlines((post_sig[xvals:xvals+1000, i] - i * 4)[0], 0, 1000, lw=0.2, color = 'k')
        ax[0].text(preAPs[0][0]-xvals - 100, post_sig[xvals, i] - i*4, i+1, size = 8)

    ax[0].vlines(preAPs[0][0]-xvals,-50, -250, lw=0.25, color = 'r')
    ax[0].set_title(f'Num swps {len(post_sig[0])}')

    y = mean_post[xvals:xvals+1000]+10
    times = np.linspace(0,len(y), len(y))/(sampl_rate*1000) #foe ms
    ax[1].plot(times, y, lw=1.5, color='k')
    ax[1].set_xlabel('ms')
    ax[1].set_ylabel('mV')

    fig.suptitle(f'{con_screen_file[end_fn:-4]} pre:{pre_cell_chan} post:{post_cell_chan}')
    fig.patch.set_facecolor('white')
    plt.show()

def get_swp_num_con_screen(con_screen_file, pre_cell_chan, post_cell_chan):
    '''
    returns the number of sweeps in the con scree that were included in the analysis
    '''
    post_sig, vm0_post = con_param.postsynaptic_screen_ALLOW_ALL(con_screen_file, post_cell_chan, [])
    
    return len(post_sig[0])


plt.switch_backend("QtAgg")
swps_nums = []
h_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
# df_repatch_cons = df_repatch_cons_all
for i, fn in enumerate(df_repatch_cons.fn):
    OP = df_repatch_cons.OP[i]
    patcher = df_repatch_cons.patcher[i]
    pre = df_repatch_cons.chan_pre[i]
    post = df_repatch_cons.chan_post[i]

    fn_path = f'{sort.get_work_dir(h_dir, OP, patcher)}{fn}'
    swps_nums.append(get_swp_num_con_screen(fn_path, pre, post))
    plot_post_cell(fn_path, pre, post)


df_repatch_cons.insert(len(df_repatch_cons.columns), 'num_swps', swps_nums)
# df_repatch_cons.to_excel(f'{data_dir}/repatched_connections_looser_QC_MARCH_2026_.xlsx')


def plot_failure_rates_repatch(df, dest_dir, w_cm = 18, h_cm = 8):
    '''
    plots the failure rates accross the repatch dataset
    '''
    if 'num_swps' not in df.columns:
        return print('reduqred columns not in dataframe')
    

    failure_rates = 100*(df.num_swps - df.num_swps_events) / df.num_swps
    df.insert(0, 'failure_rates', failure_rates)

    params = ['failure_rates', 'num_swps', 'num_swps_events']
    ax_elements = {'failure_rates' : (0,0),
                   'num_swps': (0,1),
                   'num_swps_events': (0,2)}
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
    if dest_dir:
        plt.savefig(dest_dir + date + '_plot_all_repatch_MEAN_sup.svg', \
        format = 'svg', bbox_inches = 'tight', dpi = 1000)
    plt.show(fig2)

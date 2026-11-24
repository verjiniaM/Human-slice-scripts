import importlib
import pandas as pd
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
import ephys_analysis.funcs_plot_paper_figures as paper_plot


data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/connectivity_paper/data/'
save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/connectivity_paper/plots/'


df_connections_QC = pd.read_excel(data_dir + 'all_slice_cons_with_failure_rates.xlsx')
df_connections_QC = df_connections_QC[df_connections_QC.repatch == 'no'].reset_index(drop = True)
df_con_percentage_all = pd.read_excel(data_dir + 'connectivity_percent_estimation.xlsx')


df_CTR_con = df_con_percentage_all[(df_con_percentage_all.treatment == 'Ctrl') | \
                               ((df_con_percentage_all.treatment == 'high K') & \
                                (df_con_percentage_all.day == 'D1'))]
df_CTR_con = df_CTR_con.dropna(subset=['hrs_after_OP', 'con_percentage_fixed'])
df_CTR_con = df_CTR_con[df_CTR_con.hrs_after_OP < 41.05] # THE OTHER CONNECTIONS ARE REPATCHED

paper_plot.plot_corr_single('hrs_after_OP', 'con_percentage_fixed', df_CTR_con, 'all', 'Ctrl', save_dir)



# Amp 1 over time in CTR
df_CTR = df_connections_QC[(df_connections_QC.treatment == 'Ctrl') \
                           | ((df_connections_QC.treatment == 'high K') & (df_connections_QC.day == 'D1'))]
df_CTR = df_CTR.dropna(subset=['hrs_after_OP', 'Lat1'])

paper_plot.plot_corr_single('hrs_after_OP', 'Amp 1', df_CTR, 'all', 'Ctrl', save_dir)
paper_plot.plot_corr_single('hrs_after_OP', 'Lat1', df_CTR, 'all', 'Ctrl', save_dir + 'lat')

## FOR PAPER
paper_plot.create_hist_sup_figures(df_CTR, 'Amp 1', 'Amplitude 1st EPSP (mV)', save_dir)
paper_plot.create_hist_sup_figures(df_CTR, 'Lat1', 'Latency 1st EPSP (ms)', save_dir)


# spontaneous
df_spontan = pd.read_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+\
                               'paper_figs_collected_checked/data/events/spontan/'+\
                                '_spontan_single_events_QCed.csv')
df_spontan_slice_CTR = df_spontan[df_spontan.source == 'slice']
df_spontan_slice_CTR = df_spontan_slice_CTR[(df_spontan_slice_CTR.treatment == 'Ctrl') | \
                               ((df_spontan_slice_CTR.treatment == 'high K') & \
                                (df_spontan_slice_CTR.day == 'D1'))]


paper_plot.plot_corr_single('hrs_after_OP', 'abs_amp_barb_pA', df_spontan_slice_CTR, 'all',\
                             'Ctrl', save_dir)


summary_df = pd.read_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/events/spontan/_slice_events_summary_QCed.csv')
summary_df = summary_df[summary_df.source == 'slice']
summary_df = summary_df[(summary_df.treatment == 'Ctrl') | \
                               ((summary_df.treatment == 'high K') & \
                                (summary_df.day == 'D1'))]

paper_plot.plot_corr_single('hrs_after_OP', 'avg_amp_adj', summary_df, 'all',\
                             'Ctrl', save_dir)
paper_plot.plot_corr_single('hrs_after_OP', 'avg_IEI_adj', summary_df, 'all',\
                             'Ctrl', save_dir)


# plotting correlations for fig 1
df_list = {'CTR_correlations_intr': [df_corrs, 'black']}
paper_plot.plot_corrs('hrs_after_OP', df_list, dir_read_corrs, destination_dir + 'mem', params_mem)

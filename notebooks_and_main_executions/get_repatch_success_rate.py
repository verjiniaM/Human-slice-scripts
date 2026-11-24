import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import connectivity_estimation as con_estm


#%%
work_dir = ''
save_plots_dir = False
#%%
#  ACQUIIRE THE NEEDED DATA


# VM_data_dir = f'{work_dir}results/human/connectivity_paper/data/'

# exp_overview = pd.read_excel(f'{work_dir}data/human/2025-03-08_experiments_overview.xlsx')

# # intr_all_VM = pd.read_csv(f'{VM_data_dir}all_intr_as_FM.csv')
# all_intr = pd.read_csv(f'{VM_data_dir}/archive/all_intr_props_all_treatments.csv')
# all_intr  = con_estm.create_unique_cell_IDs(all_intr, 'all')
# # remove the _D2 suffix
# all_intr['slice_name'] = all_intr['slice'].str[:2]

# df_repatch_non_QC = all_intr[all_intr.repatch == 'yes']
# df_repatch = df_repatch_non_QC[df_repatch_non_QC.groupby('cell_ID_new')['cell_ID_new'].transform('size') == 2].copy()

# # list of OPs to exclude
# exclude_OPs = exp_overview[exp_overview.repatch == 'no'].OP.tolist()

# # ops with technical problem
# al_ops = sorted(all_intr.OP.unique())
# indx_start = al_ops.index('OP230810')
# indx_end = al_ops.index('OP240201')
# ops_technical_problem = sorted(all_intr.OP.unique())[indx_start : indx_end]

# ops_keep = set(exp_overview[exp_overview.repatch == 'no'].OP) & set(df_repatch.OP.unique())

# list_OPs_exclude = list(set(exclude_OPs + ops_technical_problem) - ops_keep)
# df_intr_filtered = all_intr[~all_intr['OP'].isin(list_OPs_exclude)].copy()


# # Count cells recorded on D1
# n_D1 = df_intr_filtered[df_intr_filtered.day == 'D1'].groupby(['patcher', 'OP', 'slice_name', 'treatment'])['cell_ID_match'].nunique()

# # Count successfully repatched cells on D2
# n_D2 = df_repatch[df_repatch.day == 'D2'].groupby(['patcher', 'OP', 'slice_name', 'treatment'])['cell_ID_new'].nunique()

# # Combine and calculate success rate
# df_success = pd.concat([n_D1.rename('n_D1'), n_D2.rename('n_D2')], axis=1).fillna(0)
# df_success['success_rate'] = df_success['n_D2'] / df_success['n_D1']
# df_success = df_success.reset_index()
# df_success = df_success[df_success.treatment != 'fixed']




#%%
# LOAD DATA DIRECTLY

df_success = pd.read_csv(f'{work_dir}repatch_sucess_rate.csv')

#%%
# plot the figure

plot_colors = {'Ctrl': "#d383cf",
               'high K': "#3a9837",
               'TTX': "#4A74DD"}

fig, ax = plt.subplots(1,2, figsize = (4,2.5), sharey = True, width_ratios=(2,1))
for k, tr in enumerate(df_success.treatment.unique()):
    df_plot = df_success[df_success.treatment == tr]
    x = np.linspace(k +0.25, k+0.75, len(df_plot))
    ax[0].scatter(x, df_plot.success_rate, c = plot_colors[tr])
    ax[0].plot([k +0.41, k+0.59], [df_plot.success_rate.mean(), df_plot.success_rate.mean()],
               linewidth = 2, c = 'k', alpha = 0.8)

ax[0].set_xticks(ticks = [0.5, 1.5, 2.5], labels = df_success.treatment.unique())
ax[0].set_xlabel('treatment')
ax[0].set_ylabel('Success rate\n (num. repatch cells D2 / num.cells D1)')
ax[0].set_title('Colored by treatment')

# other subplot colored by the op
n = len(df_success['OP'].unique())  # Number of unique observations/groups
colors = plt.get_cmap('tab20', n)(np.arange(n))
color_map = dict(zip(df_success['OP'].unique(), colors))

x = np.linspace(k +0.25, k+0.75, len(df_success))
for g, OP in enumerate(df_success.OP.unique()):
    df_op = df_success[df_success.OP == OP]
    x_op = np.linspace(0.25+g/40, 0.75+g/40, len(df_op))
    ax[1].scatter(x_op, df_op.success_rate, c = color_map[OP], alpha = 0.75)

ax[1].plot([0.84, 1.16], [df_success.success_rate.mean(), df_success.success_rate.mean()],
           linewidth = 2, c = 'k', alpha = 0.8)
ax[1].text(0.82, 0.27, f'{round(df_success.success_rate.mean(), 2)}', fontsize = 9)

ax[1].set_xticks(ticks = [0, 1, 2], labels = ['', '', ''])
ax[1].set_title('Colored by patient')

if save_plots_dir:
    plt.savefig(f'{save_plots_dir}sup_repatch_success_rate.svg',
                format='svg')

plt.show()

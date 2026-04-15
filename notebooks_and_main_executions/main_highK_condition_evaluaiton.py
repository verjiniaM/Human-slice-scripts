# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast
import ephys_analysis.funcs_sorting as sort
import ephys_analysis.funcs_human_characterisation as hcf
from ephys_analysis.funcs_plotting_raw_traces import plot_trace

plt.style.use('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/Human-slice-scripts/style_plot_paper.mplstyle')

# FUNCTIONS

# GETTING DATA HiK wash in cellular
human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
# OPs where high K was washed in or out
OPs = ['OP221116', 'OP230126', 'OP231207', 'OP220413', 'OP230914', 'OP220615']
# cell attached: 'OP230330', 'OP230417'
OPs = sorted(OPs)

for i, OP in enumerate(OPs):
    work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, 'Verji', 'old')
    summary_df = pd.read_excel(f'{human_dir}/data_verji/{OP}/data_tables/RMP_high_K_{OP}.xlsx')
    if len(summary_df) == 0:
        continue
    
    # take one measurement per row, averages the RMP
    rmps, comment, condition = [], [], []
    for j, x in enumerate(summary_df.resting_potential):
        inc = summary_df.incubation_solution[j]
        rec_in = summary_df.recording_in[j]
        if inc == rec_in:
            condition.append('baseline')
        elif inc == 'Ctrl' and rec_in == 'high K':
            condition.append('wash in')
        elif inc == 'high K' and rec_in == 'Ctrl':
            condition.append('wash out')
        else:
            condition.append(rec_in)
            print(f'for {summary_df.filename[j]} chan {summary_df.cell_ch[j]} condition {rec_in}')

        rmp = ast.literal_eval(x)
        sd = np.std(rmp)
        if sd > 2:
            print(f'in depth analysis needed for {summary_df.filename[j]} chan {summary_df.cell_ch[j]}')
            comment.append('add more dps')
        else:
            comment.append('')
        rmps.append(np.mean(rmp))

    summary_df['RMP'] = rmps
    summary_df['comment'] = comment
    summary_df['condition'] = condition
    if i == 0:
        df = summary_df
    else:
        df = pd.concat([df, summary_df]).reset_index(drop = True)


df = df[(df.RMP > - 70) & (df.RMP < - 35)].reset_index(drop = True)
# add cell id

# 



df.to_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/high_K_wash_experiemnts.xlsx')
# %%
# plotting

cell_id_counts = df['cell_ID'].value_counts()
repeated_cell_ids = cell_id_counts[cell_id_counts > 1].index

# Filter dataframe to keep only rows with repeated cell_IDs
df_repeated = df[df['cell_ID'].isin(repeated_cell_ids)].reset_index(drop=True)

cells_skip = ['vm22413S3c3']
cells_long = []
# Create lists to store cell_IDs for wash_in and wash_out experiments
y_out_1, y_out_2, y_in_1, y_in_2, y_puff_1, y_puff_2 = [], [], [], [], [], []
for i, cell_ID in enumerate(df_repeated.cell_ID.unique()):
    if cell_ID in cells_skip:
        continue
    df_cell = df[df.cell_ID == cell_ID]
    if len(df_cell) == 1:
        # print(f'I am here {i}')
        if df_cell['comment'].iloc[0] != '':
            cells_long.append(i)
            rmps = ast.literal_eval(df_cell.resting_potential)
            if df_cell.condition == 'wash out':
                y_out_1.append(rmps[0])
                y_out_2.append(rmps[-1])
            elif df_cell.condition == 'wash in':
                y_in_1.append(rmps[0])
                y_in_2.append(rmps[-1])
            else:
                y_puff_1.append(rmps[0])
                y_puff_2.append(rmps[-1])
    else:
        # print(f'I am here {i}{i}')
        if len(df_cell.condition.unique()) > 1:
            if 'baseline' in df_cell.condition.unique():
                y = df_cell[df_cell.condition == 'baseline'].RMP
                # if len(y) > 1:
                # y_out_1.append(np.mean(y))
                if 'wash out' in df_cell.condition.unique():
                    y_out_1.append(np.mean(y))
                    y_2 = df_cell[df_cell.condition == 'wash out'].RMP
                    y_out_2.append(np.mean(y_2))
                if 'wash in' in df_cell.condition.unique():
                    y_in_1.append(np.mean(y))
                    y_2 = df_cell[df_cell.condition == 'wash in'].RMP
                    y_in_2.append(np.mean(y_2))

# Create the wash_in and wash_out dataframes
save_data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
df_wash_in = pd.DataFrame({'pre': y_in_1, 'post':y_in_2})
df_wash_out = pd.DataFrame({'pre': y_out_1, 'post':y_out_2})
df_wash_in.to_csv(save_data_dir + 'highKwashIN.csv')
df_wash_out.to_csv(save_data_dir + 'highKwashOUT.csv')


# plot
fig, ax = plt.subplots(1, 3, figsize=(18/2.54, 5/2.54),
                       gridspec_kw={'width_ratios': [2, 1, 1]})
plot_col = "#5274B6EA"
C = 0
df_long = df[df['comment'] != '']
rmps = []
for j, r in enumerate(df_long.resting_potential):
    rmp = ast.literal_eval(r)
    if len(rmp) < 10:
        continue
    rmps.append(rmp)
    x = np.linspace(0, 150, len(rmp))
    ax[C].plot(x, rmp, color = plot_col, alpha = 0.85)
ax[C].plot(x, np.mean(np.array(rmps), axis = 0), linewidth = 2.5, color = 'k', alpha = 0.8,
           linestyle = ':')

A = 1
excl = []
for i, a in enumerate(y_in_1):
    b = y_in_2[i]
    if abs(a - b) < 6 or b < a or a > -55:
        excl.append(i)
        continue
    # ax[A].scatter(0, a, color = plot_col)
    # ax[A].scatter(1, b, color = plot_col)
    ax[A].scatter(0, a, alpha = 0.7, facecolor = plot_col,
                            edgecolor = 'darkgrey')
    ax[A].scatter(1, b, alpha = 0.7, facecolor = plot_col,
                            edgecolor = 'darkgrey')
    ax[A].plot([0, 1], [a, b], color = plot_col, zorder = 0)

Y_1 = [it for i, it in enumerate(y_in_1) if i not in excl]
Y_2 = [it for i, it in enumerate(y_in_2) if i not in excl]
ax[A].scatter(0, np.mean(Y_1), marker = '_', s = 150, c = 'k')
ax[A].scatter(1, np.mean(Y_2), marker = '_', s = 150, c = 'k')
ax[A].set_xticks([-0.3, 0, 1, 1.3], ['', 'CTR', 'HIK', ''])
print(f'In {np.mean(Y_1)}')
print(f'In {np.mean(Y_2)}')

# plot the out
B = 2
excl = []
for i, a in enumerate(y_out_1):
    b = y_out_2[i]
    if abs(a - b) < 6:
        excl.append(i)
        continue
    # ax[B].scatter(0, a, color = plot_col)
    # ax[B].scatter(1, b, color = plot_col)
    ax[B].scatter(0, a, alpha = 0.7, facecolor = plot_col,
                            edgecolor = 'darkgrey')
    ax[B].scatter(1, b, alpha = 0.7, facecolor = plot_col,
                            edgecolor = 'darkgrey')
    ax[B].plot([0, 1], [a, b], color = plot_col, zorder = 0)
Y_1 = [it for i, it in enumerate(y_out_1) if i not in excl]
Y_2 = [it for i, it in enumerate(y_out_2) if i not in excl]
ax[B].scatter(0, np.mean(Y_1), marker = '_', s = 150, c = 'k')
ax[B].scatter(1, np.mean(Y_2), marker = '_', s = 150, c = 'k')
ax[B].set_xticks([-0.3, 0, 1, 1.3], ['', 'HiK', 'CTR', ''])
print(np.mean(Y_1))
print(np.mean(Y_2))

x_ticks = np.linspace(-80, -30, 5)
for g in [A, B, C]:
    ax[g].set_yticks(x_ticks)
ax[0].set_ylabel('Membrane potential (mV)')

# Save if directory provided
if save_dir:
    plt.savefig(f'{save_dir}highK_sum_figure.svg',
                format='svg', bbox_inches='tight', dpi=300)

# # check if there are APs
# aps = {} 
# # 18 / 84 cells fire APs
# for i, fn in enumerate(df.filename):
#     chan = df.cell_ch[i]
#     OP = df.OP[i]
#     work_dir, _, _, _ = sort.get_OP_metadata(human_dir, OP, 'Verji', 'old')
#     data = hcf.load_traces(work_dir + fn)[f'Ch{chan}'][0]
#     # data_f = data.flatten()
#     for swp in range(np.shape(data)[1]):
#         data_f = data[:, swp]
#         if len(data_f[data_f > 0]) > 0:
#             aps[f'{fn[:-4]}_ch_{chan}'] = [swp]
#             print(f'rec in {df.recording_in[i]}\n\
#                   inc {df.incubation_solution[i]}\n\
#                     condition {df.condition[i]}')
#             # plot_trace(work_dir + fn, [swp], chan)

# %%

import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
from ephys_analysis.funcs_for_results_tables import collect_intrinsic_df
import matplotlib.pyplot as plt
import numpy as np
import ephys_analysis.funcs_for_results_tables as get_results


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

DEST_DIR = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/Rosie_presentation_feb25/'

adult_df = load_intr_data('all')
adult_df = adult_df[adult_df['high K concentration'] == '8 mM']
adult_df = adult_df[adult_df.treatment != 'TTX']
adult_df, age_color_dict = pl_intr.get_age_color_dict(adult_df)

adult_df_slice_comparison = adult_df.loc[adult_df['repatch'] == 'no']
adult_df_repatch = pl_intr.get_repatch_df(adult_df)
adult_df_repatch = adult_df_repatch.sort_values(['cell_ID_new', 'treatment'])


x_labels =  ['D1','D2 \n after aCSF', 'D1', 'D2 \n after high K']
scatter_dot_size = 190
scatter_alpha = 0.8
ticks_and_text_size = 22
fig_x = 11
fig_y = 8

df = adult_df_repatch

titles_dict = pl_intr.dict_for_plotting()

params = ['TH']

if isinstance(params, list):
    keys_to_remove = [s for s in  list(titles_dict.keys()) if s not in params]
    titles_dict = {k: v for k, v in titles_dict.items() if k not in keys_to_remove}

data_type = 'repatch'
treatments = ['Ctrl', 'high K']
colors = ['gold', 'orange', 'gold', 'mediumseagreen']


for u, param in enumerate(titles_dict.keys()): 
    fig2 = plt.figure(figsize=(fig_x,fig_y))
    ax = plt.subplot(1,1,1)
    for i, tr in enumerate(treatments):
        x_plot =[]

        for j, day in enumerate(['D1', 'D2']):
            print(day, tr)
            k = j + 2*i 
            df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
            df_plot.reset_index(drop=True, inplace=True)
            median = df_plot[param].median()
            x = np.linspace(0.7+k, 1.3+k, len(df_plot))
            x_plot.append(x)       
            ax.scatter(x, df_plot[param], c = colors[int(k)], s = scatter_dot_size, zorder = 2, alpha = scatter_alpha)
            ax.scatter(1+k, median, color = 'k', marker = '_', s = 2000, zorder = 2)
            ax.text(0.85+k, (median + abs(.05*median)), str(round(median,2)), size = ticks_and_text_size)
            x_plot.append(x)
            if k in [1,3]:
                for c, cell in enumerate(df_plot['cell_ID_new']):
                    x1 = [x_plot[0][c], x[c]] 
                    y = df[param][df['cell_ID_new'] == cell]
                    plt.plot(x1, y, '-', color = colors[int(k)], alpha = 0.6, linewidth = 2, zorder = 1)

#             ax.text(1 + 2*i, int(np.max(df[param])), tr, size = 17, c = colors[2*i+1])
#             ax.text(1.7 + 2*i, int(np.max(df[param])), 'n = ' + str(len(df_plot)), size = 12, c = colors[2*i+1])

    ax.set_xticks(ticks = [1,2,3,4], labels = x_labels , size = ticks_and_text_size)
    ax.set_ylabel(titles_dict[param][1], size = ticks_and_text_size +2)
    ax.set_yticks(ticks = titles_dict[param][2])
    ax.tick_params(axis='y', labelsize = ticks_and_text_size)

    #plt.title(titles_dict[param][0], fontsize = 19, x = 0.5, y = 1)

    fig2.patch.set_facecolor('white')
    fig2.tight_layout()
    plt.savefig(DEST_DIR  + '_plot_' + param + '.png')

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


# def plot_param_for_days_slice(df, color_dict, params = None):
    
#     '''
#     params (list of params), for example ['TH']
#      '''
#     treatments = ['Ctrl', 'high K']
#     titles_dict = pl_intr.dict_for_plotting()

#     colors = ['#dede00', '#ff7f00', '#dede00', '#4daf4a','#dede00','#984ea3'] #Ctrl, 8mM Kcl, 15 mM KCl

#     if isinstance(params, list):
#         keys_to_remove = [s for s in  list(titles_dict.keys()) if s not in params]
#         titles_dict = {k: v for k, v in titles_dict.items() if k not in keys_to_remove}

#     for param in titles_dict.keys(): 
#         fig2 = plt.figure(figsize=(fig_x,fig_y))
#         ax = plt.subplot(1,1,1)
#         for i, tr in enumerate(treatments):
#             for j, day in enumerate(['D1', 'D2']):
#                 print(day, tr)
#                 k = j + 2*i 
#                 df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
#                 median = df_plot[param].median()
#                 x = np.linspace(0.65+k, 1.35+k, len(df_plot))

#                 _ = ax.boxplot(df_plot[param], positions = [k + 0.5])
#                 #, notch = True, patch_artist=True, boxprops=dict(facecolor=colors[k], alpha = 0.75),
#                 #medianprops = dict(linewidth=2.3, color = 'k'))    
#                 ax.scatter(x, df_plot[param], c = colors[int(k)], s = scatter_dot_size, alpha = scatter_alpha)
                
#                 ax.scatter(1+k, median, color = 'k', marker = '_', s = 2000, zorder = 2)
#                 ax.text(0.85+k, (median + abs(.05*median)), str(round(median,2)), size = ticks_and_text_size)

#     #         ax.text(1 + 2*i, int(np.max(df[param])), tr, size = 17, c = colors[2*i+1])
#     #         ax.text(1.7 + 2*i, int(np.max(df[param])), 'n = ' + str(len(df_plot)), size = ticks_and_text_size, c = colors[2*i+1])

#         ax.tick_params(axis='y', labelsize=22)
#         ax.set_xticks(ticks = [0.7,1.7,2.7,3.7], labels = x_labels, size = ticks_and_text_size)

#         #plt.title(titles_dict[param][0], fontsize = 19, x = 0.5, y = 1)
#         ax.set_ylabel(titles_dict[param][1], size = ticks_and_text_size +2)
#         ax.set_yticks(ticks = titles_dict[param][2])
#         ax.tick_params(axis='y', labelsize = ticks_and_text_size)

#         plt.subplots_adjust(hspace=0.35)
#         fig2.patch.set_facecolor('white')
#         fig2.tight_layout()
#         plt.show()
    
#         plt.savefig(DEST_DIR  + '_plot_' + param + '.pdf')
# plt.close(fig2)


# repatch


all_cons_IC, all_cons_VC =  get_results.collect_connections_df()
op_color_dict = pl_intr.get_op_color_dict(all_cons_IC)

all_cons_IC['Amp 1'] = all_cons_IC['Amp 1'].fillna(0)
all_cons_IC['Amp 2'] = all_cons_IC['Amp 2'].fillna(0)
all_cons_IC['Amp 3'] = all_cons_IC['Amp 3'].fillna(0)
all_cons_IC['Amp 4'] = all_cons_IC['Amp 4'].fillna(0)
all_cons_IC.columns[all_cons_IC.isna().any()].tolist()

connect_QC_passed = pl_intr.get_QC_connectivity_df(all_cons_IC)


df = df_connect_repatch
df = df.sort_values(by = ['day','connection_ID'], ascending = False)

data_type = 'repatch'

for u, param in enumerate(titles_dict_con.keys()): 
    fig2 = plt.figure(figsize=(8,fig_y))
    ax = plt.subplot(1,1,1)
    for i, tr in enumerate(treatments):
        x_plot =[]

        for j, day in enumerate(['D1', 'D2']):
            print(day, tr)
            k = j + 2*i 
            df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
            df_plot.reset_index(drop=True, inplace=True)
            median = df_plot[param].median()
            x = np.linspace(0.7+k, 1.3+k, len(df_plot))
            x_plot.append(x)       
            ax.scatter(x, df_plot[param], c = colors[int(k)], s = scatter_dot_size, zorder = 2, alpha = scatter_alpha)
            ax.scatter(1+k, median, color = 'k', marker = '_', s = 2000, zorder = 2)
            ax.text(0.9+k, (median + abs(.05*median)), str(round(median,1)), size = ticks_and_text_size)
            x_plot.append(x) 
            if k in [1,3]:
                for c, cell in enumerate(df_plot['connection_ID']):
                    x1 = [x[c], x_plot[0][c]] 
                    y = df[param][df['connection_ID'] == cell]
                    plt.plot(x1, y, '-', color = colors[int(k)], alpha = 0.6, linewidth = 2, zorder = 1)

#             ax.text(1 + 2*i, int(np.max(df[param])), tr, size = 17, c = colors[2*i+1])
#             ax.text(1.7 + 2*i, int(np.max(df[param])), 'n = ' + str(len(df_plot)), size = 12, c = colors[2*i+1])

    ax.set_xticks(ticks = [1,2,3,4], labels = x_labels , size = ticks_and_text_size)
    ax.set_ylabel(titles_dict_con[param][1], size = ticks_and_text_size +2)
    ax.set_yticks(ticks = titles_dict_con[param][2])
    ax.tick_params(axis='y', labelsize = ticks_and_text_size)

    #plt.title(titles_dict[param][0], fontsize = 19, x = 0.5, y = 1)

    fig2.patch.set_facecolor('white')
    fig2.tight_layout()

    date = str(datetime.date.today())
    plt.savefig(save_dir + 'con_screen/' + data_type + '_' + date + '_' + param + '.pdf')
    plt.close(fig2)



IFF_collected = get_results.collect_IFF_dfs(HUMAN_DIR = \
    '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/')

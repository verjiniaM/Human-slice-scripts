import glob
import funcs_sorting as sort
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
import math
import datetime
import funcs_for_results_tables as get_results
import funcs_plot_intrinsic_props as pl_intr

plt.style.use('./style_plot_intrinsic.mplstyle')

human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'

OP_dirs = glob.glob(human_dir + 'data_*/' + 'OP*')
OP_skip = ['OP210304', 'OP210319']

OPs, slices, days, treatments, num_chans, possible_cons, num_connections  = [], [], [], [], [], [], []


for i in range(15):
    path_json = glob.glob(OP_dirs[i] + '/*_con_screen_only.json')
    OP =  OP_dirs[i][OP_dirs[i].find('OP'):]

    if OP in OP_skip:
        continue

    if len(path_json) == 0:
        continue

    con_dict = sort.open_json(path_json[0])[0]
    for j in range(len(con_dict['active_chans'])):
        OPs.append(OP)
        slices.append(con_dict['slices'][j])
        if len(con_dict['slices'][j]) < 3:
            day = 'D1'
        else:
            day = 'D2'
        days.append(day)
        treatments.append(con_dict['treatment'][j])
        num_chans.append(len(con_dict['active_chans'][j]))
        possible_cons.append(len(con_dict['active_chans'][j])*(len(con_dict['active_chans'][j])-1))
        if len(con_dict['pre_chans']) == 0:
            num_connections.append(0)
            # print('For ',OP , 'slice ' , con_dict['slices'][j] , 'on ' + day + 'treatment ' + con_dict['treatment'][j], 'num chans' + str(len(con_dict['active_chans'][j])) + 'num possible connections' , str(len(con_dict['active_chans'][j])*(len(con_dict['active_chans'][j])-1)), con_dict['treatment'][j] ,  'no connectios')
        else:
            num_connections.append(len(con_dict['pre_chans'][j]))
            # print('For' , OP, 'slice ' , con_dict['slices'][j] , 'on ' , day , 'treatment ' ,con_dict['treatment'][j], 'num chans' , str(len(con_dict['active_chans'][j])) , 'num possible connections', str(len(con_dict['active_chans'][j])*(len(con_dict['active_chans'][j])-1)) ,con_dict['treatment'][j] , str(len(con_dict['pre_chans'][j])) ,'connectios found')

df_connections = pd.DataFrame({ 'OP':OPs, 'slices':slices, 'day':days, 'treatment':treatments,
    'num_channels':num_chans, 'num_possible_connections':possible_cons, 'found_connections':num_connections})
df_connections.loc[df_connections['day'] == 'D1']

# df_d1 = df_connections.loc[df_connections['day'] == 'D1']
# df_d2 = df_connections[df_connections.day == 'D2']

df_connections.insert(len(df_connections.columns), 'con_percentage', df_connections.num_possible_connections / df_connections.found_connections)

#calcualte connectivity percentage
con_percent_fixed = []
for con_percent in df_connections['con_percentage']:
    if math.isinf(con_percent):
        con_percent_fixed.append(0)
    else:  
        con_percent_fixed.append(con_percent)

df_connections.insert(len(df_connections.columns), 'con_percentage_fixed', con_percent_fixed)


#plot
fig, ax = plt.subplots(1, 1,figsize=(7,7))
colors = ['pink', 'green', 'orange']

mean_con = []

for j, day in enumerate(['D1', 'D2']):
    for i, tr in enumerate(['TTX', 'high K', 'Ctrl']):
        if day == 'D1':
            k_ = 1
        else:
            k_ = 2
        # k_ = j + 2*i 
        df_plot = df_connections.loc[(df_connections.treatment == tr) & (df_connections.day == day)]
        avg = np.mean(df_plot.con_percentage_fixed)
        #print(median)
        x = np.linspace(k_ - 0.3, k_ + 0.3, len(df_plot))
        ax.scatter(x, df_plot.con_percentage_fixed, c = colors[i])
        if k_ == 1:
            ax.scatter(k_, avg, marker = '_', s = 2000, c = colors[i])
        else:
            ax.scatter(k_, avg, marker = '_', s = 2000, c = colors[i], label = tr)
        ax.text(k_, 1.2*avg, str(round(avg,0)), size = 10)
        ax.set_xticks(ticks = [1,2], labels = ['D1', 'D2'])

    mean_all = np.mean(df_connections[df_connections.day == day].con_percentage_fixed)
    mean_con.append(mean_all)
    ax.scatter(k_, mean_all, marker = '_', s = 2000, c = 'k')

fig.suptitle('percentage connectivity on \n D1: ' + str(round(mean_con[0],2)) + '; on D2: ' + str(round(mean_con[1],2)))
plt.figlegend(loc = 'upper right',  bbox_to_anchor=(0.95, 0.85))
date = str(datetime.date.today())
plt.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/connectivity/connectivity_day1_day2/' + date + '_con_D1_D2.png')


#%% 
# estimating connectivity based on cells passed QC 

# loading data
# QC passed, age > 10, hrs incubation > 16
df_intr = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/summary_data_tables/intrinsic_properties/QC_passed_2024-07-04_collected.xlsx')
df_con, all_cons_VC = get_results.collect_connections_df()
df_con = pl_intr.filter_on_hrs_incubation(df_con, 16, 30)
df_con = pl_intr.change_to_numeric(df_con)
df_con = df_con.reset_index(drop = True)

#deleting cells that haven't passed the QC of intr_props
del_indx = []
for i, con in enumerate(df_con.connection_ID):
    cell_ID_pre = con[:-3] + str(df_con.chan_pre[i])
    cell_ID_post = con[:-3] + str(df_con.chan_post[i])

    if cell_ID_pre not in df_intr.cell_ID.tolist() or cell_ID_post not in df_intr.cell_ID.tolist():
        del_indx.append(i)
for cell in del_indx:
    df_con = df_con.drop(index = cell)     
df_con.reset_index(inplace = True, drop = True)


# estiamting connectivity
df_intr = df_intr.sort_values(['day'])
OPs, slices, day, tr, num_cells, num_possible_connections, num_connections, con_percentage = [], [], [], [], [], [], [], []

for OP in df_intr.OP.unique():
    df_OP_intr = df_intr[df_intr.OP == OP]

    if OP not in df_con.OP.unique():
        slice_num = len(df_OP_intr.slice.unique())
        OPs.extend([OP] * slice_num)
        for slic in df_OP_intr.slice.unique():
            tr.append(df_OP_intr['treatment'][df_OP_intr['slice'] == slic].values[0])
            slices.append(slic)
            if len(slic) > 2:
                day.append('D2')
            else:
                day.append('D1')
            cell_num = len(df_OP_intr[df_OP_intr['slice'] == slic])
            num_cells.append(cell_num)
            num_possible_connections.append(cell_num * (cell_num -1))
            num_connections.append(0)
            con_percentage.append(0)
        continue

    df_OP_con = df_con[df_con.OP == OP]
    for slic in df_OP_intr.slice.unique():
        if slic not in df_OP_con.slice.unique():
            OPs.append(OP)
            slices.append(slic)
            tr.append(df_OP_intr['treatment'][df_OP_intr['slice'] == slic].values[0])
            if len(slic) > 2:
                day.append('D2')
            else:
                day.append('D1')
            cell_num = len(df_OP_intr[df_OP_intr['slice'] == slic])
            num_cells.append(cell_num)
            num_possible_connections.append(cell_num * (cell_num -1))
            num_connections.append(0)
            con_percentage.append(0)
            continue
        
        OPs.append(OP)
        slices.append(slic)
        tr.append(df_OP_intr['treatment'][df_OP_intr['slice'] == slic].values[0])
        if len(slic) > 2:
                day.append('D2')
        else:
            day.append('D1')
        
        cell_num = len(df_OP_intr[df_OP_intr['slice'] == slic])
        num_cells.append(cell_num)
        num_possible_connections.append(cell_num * (cell_num -1))
        num_connections.append(len(df_OP_con[df_OP_con['slice'] == slic]))
        con_percentage.append((cell_num * (cell_num -1)) / (len(df_OP_con[df_OP_con['slice'] == slic])))
        

df_connec_intr_QC = pd.DataFrame({ 'OP':OPs, 'slices':slices, 'day':day, 'treatment':tr,
    'num_cells':num_cells, 'num_possible_connections':num_possible_connections, 
    'found_connections':num_connections, 'con_percentage':con_percentage})
df_connec_intr_QC = df_connec_intr_QC[df_connec_intr_QC['treatment'] != 'wash in high K']


# plot
fig, ax = plt.subplots(1, 1,figsize=(7,7))
colors = ['pink', 'green', 'orange']

df_connections = df_connec_intr_QC
mean_con = []
for j, day in enumerate(['D1', 'D2']):
    for i, tr in enumerate(['TTX', 'high K', 'Ctrl']):
        if day == 'D1':
            k_ = 1
        else:
            k_ = 2
        # k_ = j + 2*i 
        df_plot = df_connections.loc[(df_connections.treatment == tr) & (df_connections.day == day)]
        avg = np.mean(df_plot.con_percentage)
        #print(median)
        x = np.linspace(k_ - 0.3, k_ + 0.3, len(df_plot))
        ax.scatter(x, df_plot.con_percentage, c = colors[i])
        if k_ == 1:
            ax.scatter(k_, avg, marker = '_', s = 2000, c = colors[i])
        else:
            ax.scatter(k_, avg, marker = '_', s = 2000, c = colors[i], label = tr)
        ax.text(k_, 1.2*avg, str(round(avg,0)), size = 10)
        ax.set_xticks(ticks = [1,2], labels = ['D1', 'D2'])

    mean_all = np.mean(df_connections[df_connections.day == day].con_percentage)
    mean_con.append(mean_all)
    ax.scatter(k_, mean_all, marker = '_', s = 2000, c = 'k')

fig.suptitle('percentage connectivity on \n D1: ' + str(round(mean_con[0],2)) + '; on D2: ' + str(round(mean_con[1],2)))
plt.figlegend(loc = 'upper right',  bbox_to_anchor=(0.95, 0.85))
date = str(datetime.date.today())
plt.savefig('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/connectivity/connectivity_day1_day2/' + date + '_con_D1_D2_QC_passed_intr.png')



# %%

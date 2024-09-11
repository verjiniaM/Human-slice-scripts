import glob
import intrinsic_props_and_connectivity.funcs_sorting as sort
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
import math

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

for j, day in enumerate(['D1', 'D2']):
    for i, tr in enumerate(['TTX', 'high K', 'Ctrl']):
        if day == 'D1':
            k_ = 1
        else:
            k_ = 2
        # k_ = j + 2*i 
        df_plot = df_connections.loc[(df_connections.treatment == tr) & (df_connections.day == day)]
        median = np.mean(df_plot.con_percentage_fixed)
        #print(median)
        x = np.linspace(k_ - 0.3, k_ + 0.3, len(df_plot))
        ax.scatter(x, df_plot.con_percentage_fixed, c = colors[i])
        ax.scatter(k_, median, marker = '_', s = 2000, c = colors[i], label = tr)
        ax.text(k_, 1.2*median, str(round(median,0)), size = 10)
        ax.set_xticks(ticks = [1,2], labels = ['D1', 'D2'])


    median_all = np.mean(df_connections[df_connections.day == day].con_percentage_fixed)
    ax.scatter(k_, median_all, marker = '_', s = 2000, c = 'white')
    print('Connectivity on day ' + day + ' is ' + str(median_all))
plt.figlegend(loc = 'upper right',  bbox_to_anchor=(0.9, 0.9))
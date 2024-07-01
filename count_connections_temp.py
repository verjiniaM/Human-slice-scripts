import glob
import funcs_sorting as sort
import pandas as pd

human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'

OP_dirs = glob.glob(human_dir + 'data_*/' + 'OP*')
OP_skip = ['OP210304', 'OP210319']

OPs, slices, day, treatments, num_chans, possible_cons, num_connections  = [], [], [], [], [], [], []

for i in range(len(OP_dirs)):
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
            day == 'D1'
        else:
            day == 'D2'
        day.append(day)
        treatments.append(con_dict['treatment'][j])
        num_chans.append(len(con_dict['active_chans'][j]))
        possible_cons.append(len(con_dict['active_chans'][j])*(len(con_dict['active_chans'][j])-1))
        if len(con_dict['pre_chans']) == 0:
            num_connections.append(0)
        else:
            num_connections.append(len(con_dict['pre_chans'][j]))

    
df_connections = pd.DataFrame({ 'OP':OPs, 'slices':slices, 'day':day, 'treatment':treatments,
    'num_channels':num_chans, 'num_possible_connections':possible_cons, 'found_connections':num_connections})

df_d1 = df_connections.loc[df_connections['day'] == 'D1']
df_d2 = df_connections[df_connections.day == 'D2']

print(len(OPs))
print(len(df_connections))

df_connections

# df_d1

# con_d1 = df_d1['num_possible_connections'] / df_d1['found_connections']
# con_d2 = df_d2.num_possible_connections / df_d2.found_connections

# print('connectivity D1: ' + str(con_d1))
# print('connectivity D2: ' + str(con_d2))
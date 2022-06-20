import pandas as pd
import os
import numpy as np
import datetime


#%%

#this function opens the alerady existing experiments_overview file and adds the latest OP info
#patcher and op_folder have to be strings; so with " "
def update_op_list(op_folder, patcher):
    human_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
    exp_view = pd.read_excel(human_folder + 'experiemnts_overview.xlsx') 
    folder_list = np.array(sorted(os.listdir(human_folder + op_folder))) #lists all items in the op_folder
    OP_folders_indx = [i for i in range(len(folder_list)) if folder_list[i][0:2] == "OP"] #index of all items that start with "OP"
    OP_folders = folder_list[OP_folders_indx].tolist()

    existing_OPs = exp_view.loc[exp_view['patcher'] == patcher, 'OP'].tolist()

    #find non-added OPs
    all_OPs = OP_folders + existing_OPs
    not_added = [i for i in all_OPs if all_OPs.count(i)==1]

    not_added_df = pd.DataFrame({"OP": not_added, "patcher":[patcher]*len(not_added)})

    exp_view = pd.concat([exp_view, not_added_df],ignore_index=True, index = False)

    #for i in range(len(not_added)):
     #   exp_view = exp_view.append({'OP': not_added[i], 'patcher': patcher}, ignore_index=True)
    exp_view.to_excel(human_folder + 'experiemnts_overview.xlsx') 

    return print("Adding OPs:" + str(not_added) + " for " + patcher)

#%% 

#Adding the cortex_out info
def add_cortex_out_time():
    human_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
    exp_view = pd.read_excel(human_folder + 'experiemnts_overview.xlsx')

    #check where the OP time is not filled in
    for i in range(len(exp_view)):
        if pd.isnull(exp_view['cortex_out'][i]):
            #enter the data and time in %y%m%d %H:%M (e.g. 220525 = 25th May 2022);
            # if time is unknown add "  " (double space)
            cortex_out =  input('Cortex out for ' + exp_view['OP'][i] + '(yymmdd hh:mm)') 
            if cortex_out == "  ":
                exp_view.at[i,'cortex_out']= np.NaN
            else:
                date_time_cortex_out = datetime.datetime.strptime(cortex_out, '%y%m%d %H:%M')
                exp_view.at[i,'cortex_out'] = date_time_cortex_out
    
    exp_view.to_excel(human_folder + 'experiemnts_overview.xlsx', index = False) 

    
# %%

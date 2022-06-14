#%%
import datetime
import numpy as np
import pandas as pd
import os

#%%
def select_folders (dir):
    #selects the folders withing a given folder
    dirlist = []
    for filename in os.listdir(dir):
        if os.path.isdir(os.path.join(dir,filename)):
            dirlist.append(filename)
    dirlist = sorted(dirlist)
    return dirlist

def collect_tables (dir, folder):
    #selects the data frame with '.xlsx' ending in the /data_tables/ folder
    loc = dir + folder + "/"
    files = os.listdir(loc)
    for j in range(len(files)):
        file = files[j]
        if file == "data_tables":
            table_dir_content = os.listdir(loc + '/data_tables/')
            for i in range(len(table_dir_content)):
                if table_dir_content[i][-5:] == '.xlsx':
                    df = pd.read_excel(loc + '/data_tables/' + table_dir_content[i])
                    return df
            #if table_name != []:
    print('For ' + folder + " no data table found. Data not analyzed yet.")       


#%%%
human_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'

# verji = human_folder + 'data_verji/'
# rosie = human_folder + 'data_rosie/'

# folders_v = select_folders(verji)
# folders_r = select_folders(rosie)

# df_all = pd.DataFrame()
# for OPs in range(len(folders_v)):
#     folder = folders_v[OPs]
#     df1 = collect_tables(verji, folder)
#     if df1 is None:
#         continue
#     df_all = df_all.append(df1, ignore_index = True)

# for OPs in range(len(folders_r)):
#     folder = folders_r[OPs]
#     df2 = collect_tables(rosie, folder)
#     if df2 is None:
#         continue
#     df_all = df_all.append(df2, ignore_index = True)


#%%
#Append new data
folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/summary_data_tables/'
file_list = sorted(os.listdir(folder))

for i in range(len(file_list)):
    if file_list[i][-11:] == '+times.xlsx': 
        df_all_old_name = folder + file_list[i]

df_new_name  = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP220322/data_tables/OP220322_Intrinsic_and_synaptic_properties.xlsx'
df_new = pd.read_excel(df_new_name)

df_all_old = pd.read_excel(df_all_old_name)
df_all = df_all_old.append(df_new, ignore_index = True)
#%%
date  = str(datetime.date.today())
df_all.to_excel(folder + date + '_complete+times.xlsx')
df_all.to_csv(folder + date + '_complete+times.csv')

# %%

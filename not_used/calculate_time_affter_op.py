#%%
import neo
import datetime
import numpy as np
import pandas as pd

#%%%
human_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
df2 = pd.read_excel(human_folder + 'summary_data_tables/' + '2022-03-12_complete.xlsx')

OP_list = sorted(df2['OP'].unique().tolist())

#%%
#OP210929_t = datetime.datetime(2021, 9, 29) #minies only 
OP201027_t = datetime.datetime(2020, 10, 27, 10, 5)
OP201029_t = datetime.datetime(2020, 10, 29, 11)
OP210319_t = datetime.datetime(2021, 3, 19, 15, 5)
OP210323_t = datetime.datetime(2021, 3, 23, 17, 10)
OP210615_t = datetime.datetime(2021, 6, 15, 10, 15)
#OP211116_t = datetime.datetime(2021, 11, 16, 14, 0) #minies only
OP211123_t = datetime.datetime(2021, 11, 23, 10,30)
OP211209_t = datetime.datetime(2021, 12, 9, 12, 35)
OP211220_t = datetime.datetime(2021, 12, 20, 13, 40)
OP220111_t = datetime.datetime(2022, 1, 11, 10, 0)
OP220120_t = datetime.datetime(2022, 1, 20, 12, 45)
OP220127_t = datetime.datetime(2022, 1, 27, 12, 50)
OP220217_t = datetime.datetime(2022, 2, 17, 9, 25)
OP220228_t = datetime.datetime(2022, 2, 28, 10, 15)
#OP220308_t = datetime.datetime(2022, 3, 8, 10, 40)
OP220322_t = datetime.datetime(2022, 3, 22, 10,21 )

times = sorted([OP201027_t, OP201029_t, OP210319_t, OP210323_t, OP210615_t, OP211123_t, OP211209_t, OP211220_t, 
OP220111_t, OP220120_t, OP220127_t, OP220217_t, OP220228_t, OP220322_t])

if len(OP_list) != len(times): 
    print('Fix OP times')

#%%

crtx_out = {}
for i in range(len(OP_list)):
    crtx_out[OP_list[i]] = times[i]

#%% 

# for all vc files, start of recording
# files_list = df2['filename'].unique().tolist()
# for file in range(len(files_list)):
#     OP = df2['OP'][df2.index[df2['filename'] == files_list[file]][0]] 
# df_rec_time = pd.DataFrame(columns=['OP', 'filename','cortex_out', 'hrs_after_OP', '(hrs, min)'])
#     #df2 = pd.read_excel(human_folder + 'VM_04.02.2022_all_data.xlsx')
#     # df2.to_csv(human_folder + 'VM_04.02.2022_all_data.csv')

#     filename = human_folder + OP + '/' + files_list[file]

#     r = neo.io.AxonIO(filename=filename)
#     block=r.read(signal_group_mode = 'split-all')[0]
#     rec_time = block.rec_datetime
    
#     out_time = crtx_out[OP]

#     dt = rec_time - out_time #result in seconds

#     h_after_op = dt.seconds/3600
#     #(hrs, mins)
#     time_after_op = datetime.time(hour = int(dt.seconds /3600), minute = int((dt.seconds / 3600 - int(dt.seconds / 3600))*60))

#     df_rec_time = df_rec_time.append({'OP':OP, 'filename':files_list[file], 'cortex_out': crtx_out[OP],
#     'hrs_after_OP':h_after_op, '(hrs, min)': time_after_op}, ignore_index=True)
#     #df_rec_time.to_excel(human_folder + 'times_after_surgery.xlsx')


# %%
# for all files
df_rec_time = pd.DataFrame(columns=['OP', 'filename','cortex_out', 'rec_time', 'dt', 'hrs_after_OP', '(hrs, min)'])

all_files = df2['filename'].tolist()
patchers = df2['patcher'].tolist()
for file in range(len(all_files)):
    OP = df2['OP'][file]

    if patchers[file] == 'Rosie':
        filename = human_folder + 'data_rosie/' + OP + '/' + all_files[file]
    else:
        filename = human_folder + 'data_verji/' + OP + '/' + all_files[file]

    r = neo.io.AxonIO(filename=filename)
    block=r.read(signal_group_mode = 'split-all')[0]
    rec_time = block.rec_datetime
    
    out_time = crtx_out[OP]

    dt = rec_time - out_time #result in seconds

    h_after_op = dt.days*24 + dt.seconds/3600
    #(hrs, mins)
    time_after_op = datetime.time(hour = int(dt.seconds /3600), minute = int((dt.seconds / 3600 - int(dt.seconds / 3600))*60))

    df_rec_time = df_rec_time.append({'OP':OP, 'filename':all_files[file], 'cortex_out': crtx_out[OP],
    'rec_time': rec_time, 'dt': dt,
    'hrs_after_OP':h_after_op, '(hrs, min)': time_after_op}, ignore_index=True)
    
#%%
date  = str(datetime.date.today())
df_rec_time.to_excel(human_folder + 'summary_data_tables/' + date + '_times_after_OP.xlsx')

df2['hrs_after_op'] = df_rec_time['hrs_after_OP']

df2.to_excel(human_folder + 'summary_data_tables/' + date + '_complete+times.xlsx')
df2.to_csv(human_folder + 'summary_data_tables/' + date + '_complete+times.csv')

#%%
# filename1 = human_folder + 'data_rosie/OP210615/2021_06_15_0036.abf'
# r1 = neo.io.AxonIO(filename=filename1)
# block1=r1.read(signal_group_mode = 'split-all')[0]
# rec_time1 = block1.rec_datetime

# filename2 = human_folder + 'data_rosie/OP210615/2021_06_16_0009.abf'
# r2 = neo.io.AxonIO(filename=filename2)
# block2=r2.read(signal_group_mode = 'split-all')[0]
# rec_time2 = block2.rec_datetime

# inc = rec_time2 - rec_time1
# inc.days*24 + inc.seconds/3600
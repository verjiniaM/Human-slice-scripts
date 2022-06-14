#%%
import neo
import datetime
import numpy as np
import pandas as pd

#%%%
human_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/'
OP = 'OP220322/'

df = pd.read_excel(human_folder + OP + '/data_tables/'+ OP[:-1] + '_Intrinsic_and_synaptic_properties.xlsx')

# %%
time_out = datetime.datetime(2022, 3, 22, 10,21 )

#%%
crtx_out = {}
crtx_out[OP[:-1]] = time_out

# %%
# for all files
df_rec_time = pd.DataFrame(columns=['OP', 'filename','cortex_out', 'rec_time', 'dt', 'hrs_after_OP', '(hrs, min)'])

all_files = df['filename'].tolist()
patchers = df['patcher'].tolist()
for file in range(len(all_files)):
    OP = df['OP'][file]

    if patchers[file] == 'Rosie':
        filename = human_folder[:-11] + 'data_rosie/' + OP + '/' + all_files[file]
    else:
        filename = human_folder[:-11] + 'data_verji/' + OP + '/' + all_files[file]

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
  
# %%
#date  = str(datetime.date.today())
df_rec_time.to_excel(human_folder + OP + '/data_tables/'+ OP[:-1] + '_times_after_OP.xlsx')

df['hrs_after_op'] = df_rec_time['hrs_after_OP']

df.to_excel(human_folder + OP + '/data_tables/'+ OP[:-1] + '_Intrinsic_and_synaptic_properties+times.xlsx')

# %%

import human_characterisation_functions as hchf
import numpy as np
import pandas as pd
import os
import events_analysis_funcs as event


#%%
#sweeps for further analysis, not old recordings
human_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
op_folder = 'data_verji/'
patcher = 'Verji'
folders = np.array(sorted(os.listdir(human_folder + op_folder)))
OPs = []
for folder in range(len(folders)):
    if folders[folder][:2] == 'OP':
        OPs.append(folders[folder])

for op in range(len(OPs)):

    OP_folder = OPs[op] + '/'

    if patcher == 'Verji': work_dir = human_folder + 'data_verji/'+ OP_folder 
    if patcher == 'Rosie': work_dir = human_folder + 'data_rosie/'+ OP_folder 
    file_list = sorted(os.listdir(work_dir))

    filenames = []
    for file in range(len(file_list)):
        if  (file_list[file][:18] == 'recording_lab_book'): 
            df_rec = pd.read_excel(work_dir + file_list[file], header = 1)
            index_spontan = df_rec.index[df_rec['protocol'] == 'vm'].tolist()
            index_mini = df_rec.index[df_rec['protocol'] == 'minis'].tolist()
        if file_list[file][-4:] == '.abf': 
            filenames.append(file_list[file])
    
    print(OPs[op])
    analyse_y_n = input("Continue with analysis of " + OPs[op] + '?')
    if analyse_y_n == 'n':
        continue
    print('spontan files: ' + str(index_spontan))
    print('mini files: ' + str(index_mini))

    all_good = input("is the indexing correct?")
    if all_good == 'n': 
        print('correct the lab book entries')
        filenames = []
        for file in range(len(file_list)):
            if  (file_list[file][:18] == 'recording_lab_book'): 
                df_rec = pd.read_excel(work_dir + file_list[file], header = 1)
                index_spontan = df_rec.index[df_rec['protocol'] == 'vm'].tolist()
                index_mini = df_rec.index[df_rec['protocol'] == 'minis'].tolist()
            if file_list[file][-4:] == '.abf': 
                filenames.append(file_list[file])
        print('spontan files: ' + str(index_spontan))
        print('mini files: ' + str(index_mini))
        input()
    del df_rec

    print(filenames[1])   
    record_sorting = pd.DataFrame(columns = ['OP','patcher', 'filename', 'mini/spontan?', 'cell_ch','swps_to_analyse', 'swps_to_discard', 
    'min_vals_each_swp', 'max_vals_each_swp', 'drift'])

    for u in range(len(index_spontan)):
        spontan = index_spontan[u]
        filename_spontan = work_dir+filenames[spontan]
        
        active_channels = [int(item) for item in input('Channels used in ' + filenames[spontan] + ' ' + OP_folder[:-1]).split()]
        for q in range(len(active_channels)):
            ch = active_channels[q]
            min_vals, max_vals, good_swps, bad_swps, drift = events.rec_stability (filename_spontan, ch , 50)

            data_to_add = pd. DataFrame({'OP':OP_folder[:-1], 'patcher':patcher, 'filename':filenames[spontan], 'mini/spontan?':'spontan',
            'slice':slice, 'cell_ch':ch, 'swps_to_analyse':[good_swps], 'swps_to_discard':[bad_swps], 'min_vals_each_swp':[min_vals], 
            'max_vals_each_swp':[max_vals], 'drift':[drift]})
            record_sorting = pd.concat([record_sorting.loc[:], data_to_add]).reset_index(drop=True)

    for h in range(len(index_mini)):
        mini = index_mini[h]
        mini_file = work_dir+filenames[mini]

        active_channels = [int(item) for item in input('Channels used in ' + filenames[mini] + ' ' + OP_folder[:-1]).split()]
        for b in range(len(active_channels)):
            ch = active_channels[b]
            min_vals, max_vals, good_swps, bad_swps, drift = events.rec_stability (mini_file, ch , 50)

            data_to_add = pd. DataFrame({'OP':OP_folder[:-1], 'patcher':patcher, 'filename':filenames[mini], 'mini/spontan?':'minis',
            'slice':slice, 'cell_ch':ch, 'swps_to_analyse':[good_swps], 'swps_to_discard':[bad_swps], 'min_vals_each_swp':[min_vals], 
            'max_vals_each_swp':[max_vals], 'drift':[drift]})
            record_sorting = pd.concat([record_sorting.loc[:], data_to_add]).reset_index(drop=True)
            
    del index_mini, index_spontan
    record_sorting.to_excel(work_dir + 'data_tables/' + OP_folder[:-1] + '_QC_measures_events.xlsx') 
    record_sorting.to_csv(work_dir + 'data_tables/' + OP_folder[:-1] + '_QC_measures_events.csv')

   #print('The data for ' + OP + 'is fully analyzed and saved in ' + work_dir + 'data_tables/')

# %%

#Experimenting with the correct drift value
#collect average drifts for each cell --> plot

#accesing the QC tables

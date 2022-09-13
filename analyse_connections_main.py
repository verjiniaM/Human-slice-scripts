#%%
import os 
import pandas as pd
import numpy as np
import human_characterisation_functions as hchf
import connection_parameters as con_param

work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/'
OP =  'OP220615' #input("OP to analyse") + '/'


#collect the con screen files 

file_list = sorted(os.listdir(work_dir+OP+'/'))

filenames = []
for file in range(len(file_list)):
    #pclamp files
    if file_list[file][-4:] == '.abf': 
        filenames.append(file_list[file])
    #sorting the labbook entries and number of sweeps
    elif (file_list[file][-5:] == '.xlsx' and file_list[file][:2] == '~$'):
        print('Close the excel file')
        input()
    elif (file_list[file][-5:] == '.xlsx' and file_list[file][:2] != '~$'): 
        df_rec = pd.read_excel(work_dir + OP+'/' + file_list[file], header = 1)
        slice_indx = df_rec.index[df_rec['slice'].notnull()]
        slice_names = df_rec['slice'][slice_indx].tolist()
        index_con_screen = df_rec.index[df_rec['protocol'] == 'con screen'].tolist()

new_slice_names = []
for i in range(len(slice_names)):
    if i < len(slice_names)-1:
        new_slice_names.append([slice_names[i]]*(slice_indx[i+1]-slice_indx[i]))
    else :
        new_slice_names.append([slice_names[i]]*(15))
slice_names  = [x for xs in new_slice_names for x in xs]

con_data = pd.DataFrame(columns = ['fn', 'slice', 'chan_pre', 'chan_post', 'Vm pre', 'Vm post', 
    'Amp1',	'Amp2',	'Amp3',	'Amp4',	'Lat1',	'Lat2',	'Lat3',	'Lat4', 'num excluded swps', 'comments'])

#%%

for i in range(len(index_con_screen)):
    con_screen_file = work_dir + OP + '/' + filenames[index_con_screen[i]] 

    active_channels =  [int(item) for item in input('Active channels ' + filenames[index_con_screen[i]]).split()]
    con_param.plot_connect(con_screen_file, active_channels)

    pre_cells = [int(item) for item in input(filenames[index_con_screen[i]] + ' pre cells: ').split()]
    post_cells = [int(item) for item in input(filenames[index_con_screen[i]] + ' corresponding post cells: ').split()]

    # pre_cells = [2,4]
    # post_cells = [4,7]

    for cell in range(len(pre_cells)):
        pre_cell = pre_cells[cell]
        post_cell = post_cells[cell]
        slic = slice_names[post_cells[cell]]

        vm_pre = hchf.restingmembrane(con_screen_file, pre_cell)
        vm_post = hchf.restingmembrane(con_screen_file, post_cell)
        
        con_param.plot_post_cell(con_screen_file, pre_cell, post_cell) 
        pre_signal, es, exclude = con_param.presynaptic_screen(con_screen_file, pre_cell)
        post_signal = con_param.postsynaptic_screen(con_screen_file, post_cell, es)
        mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = con_param.get_analysis_window(pre_signal, post_signal)
        pre_signal, post_signal = con_param.remove_sweeps_with_POSTsynAPs(pre_signal, post_signal, preAPs)
        post_peaks = con_param.find_postsynaptic_peaks(post_window, preAPs_shifted)
        post_local_baseline = con_param.get_psp_baselines(post_window,preAPs_shifted)
        onsets = con_param.get_onsets(preAPs_shifted, post_window, post_peaks, post_local_baseline)
        latency = con_param.latencies(onsets, preAPs_shifted)

        con_param.plot_connection_window(con_screen_file, pre_cell, post_cell, pre_window, \
            post_window, preAPs_shifted, post_signal,onsets, preAPs, post_peaks, post_local_baseline)

        amps = []
        for u in range(4):
            amps.append(post_peaks[u][0] - post_local_baseline[0][0])

        # print(amps)
        # print(latency)

        df_add = pd.DataFrame({'fn': con_screen_file, 'slice': slic, 'chan_pre': pre_cell,
        'chan_post': post_cell, 'Vm pre' :vm_pre, 'Vm post': vm_post,
        'Amp1': amps[0], 'Amp2': amps[1],	'Amp3': amps[2],	'Amp4': amps[3],	
        'Lat1': latency[0][0],	'Lat2' : latency[1][0],	'Lat3': latency[2][0],	'Lat4': latency[3][0], 
        'num excluded swps': len(es), 'comments':exclude}, index=[0])
        con_data = pd.concat([con_data.loc[:], df_add]).reset_index(drop=True)

con_data.to_excel(work_dir + OP  + '/data_tables/' + OP + '_connected_cell_properties.xlsx') 


#%%
#analysis of the postsynaptic cell stimulated at different frequencies
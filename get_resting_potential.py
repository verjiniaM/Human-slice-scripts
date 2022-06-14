#%%

# Have to change this file to get the resring potential from the characterization file
import os
import human_characterisation_functions as hcf
import human_synaptic_functions1 as hsf
import trace_names_check as tn
import pandas as pd

#%%
human_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
H_all_files = sorted( os.listdir(human_folder))
OPs = [string for string in H_all_files if 'OP' in string]

dir_data = "resting_tables"
path = os.path.join(human_folder , dir_data)
if os.path.isdir(path) == False: os.mkdir(path) 

#%% 
for num in range(len(OPs)):
    OP = OPs[num] + '/'
    work_dir = human_folder + OP
    file_list = sorted(os.listdir(work_dir))

    filenames = []
    for i in range(len(file_list)):
        if file_list[i][-4:] == '.abf': 
            filenames.append(file_list[i])
        elif file_list[i][-5:] == '.xlsx': 
            df_rec = pd.read_excel(work_dir + file_list[i], header = 1)
            slice_indx = df_rec.index[df_rec['slice'].notnull()]
            slice_names = df_rec['slice'][slice_indx].tolist()
            index_vc = df_rec.index[df_rec['protocol'] == 'vc'].tolist()
            index_char = df_rec.index[df_rec['protocol'] == 'freq analyse'].tolist()
            index_vm = df_rec.index[df_rec['protocol'] == 'vm_mouse'].tolist()

    print(OP)
    print('vc:' + str(index_vc))
    print('vm:' + str(index_vm))
    print('fa:' + str(index_char))
    if len(index_vc) != len(index_vm): 
        print('Fix protocol names. Uneuqal number of VC and freq analyse protocold')
        continue

    #sorting the files based on lab journal entry and number of sweeps
    #plots and variables saved in indicated folders (see output)
    df = pd.DataFrame(columns=['filename','slice', 'cell_ch', 'cell_ID', 'resting_potential'])

    for i in range(len(index_vc)):
        vm = index_vm[i]
        slice = slice_names[i]
        char = index_char[i]

        filename_vm = work_dir+filenames[vm]

        #give inputs with space between entries
        active_channels = [int(item) for item in input('Which channels were active in' + filenames[vm-1] +'(look at plots)').split()]

        for j in range(len(active_channels)):
            ch = active_channels[j]
            resting_mem = hcf.restingmembrane(filename_vm, ch)

            cellID = filenames[char][:-7]+slice+'c'+str(ch)
            df = df.append({'filename':filenames[char],'slice':slice, 'cell_ch':ch, 'cell_ID':cellID,
            'resting_potential': resting_mem}, ignore_index=True)

    df2 = pd.read_excel(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.xlsx')

    #df.to_excel(human_folder + 'resting_tables/' + OP[:-1] + '_resting.xlsx') 
    #df2.to_csv(human_folder + 'data_tables/' + OP[:-1] + '_resting.csv')
# %%

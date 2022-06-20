
#%%%
import os
from re import A
import human_characterisation_functions as hcf
import human_synaptic_functions1 as hsf
import trace_names_check as tn
import pandas as pd
import glob

#loading the updated experiments_overview and the old summary_data_table
human_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'

exp_view = pd.read_excel(human_folder + 'experiemnts_overview.xlsx') 
latest_sum_data_table = glob.glob(human_folder + 'summary_data_tables/' + '*.xlsx')
summary_data_table = pd.read_excel(latest_sum_data_table[0])

#%%
#taking the unique OP IDs and updating the full OP list

last_update_op_list = summary_data_table.OP.unique().tolist()
newest_op_list = exp_view.OP.unique().tolist()
all_OPs = last_update_op_list + newest_op_list
op_to_analyse = [i for i in all_OPs if all_OPs.count(i)==1]

#%%
#Decide which OPs need to be analyzed 

for i in range(len(op_to_analyse)):
    OP = op_to_analyse[i]
    y_or_n = input('Do you want to start analysis of ' + op + '(y/n)?') 
    if y_or_n == "y":
        #input some functions
        OP_folder = OP + '/'
        tissue_source = input('Tissue source (Bielefeld, Mitte, Virchow): ')
        patcher = input('Patcher (Verji or Rosie): ')
        age = input('Patient age (if not known: A for adult, J for juvenile): ')
        print('starting analysis of '+ OP)
        work_dir = human_folder + OP
        file_list = sorted(os.listdir(work_dir))

        filenames = []
        for i in range(len(file_list)):
            #pclamp files
            if file_list[i][-4:] == '.abf': 
                filenames.append(file_list[i])
            #lab book
            elif file_list[i][-5:] == '.xlsx': 
                df_rec = pd.read_excel(work_dir + file_list[i], header = 1)
                slice_indx = df_rec.index[df_rec['slice'].notnull()]
                slice_names = df_rec['slice'][slice_indx].tolist()
                index_vc = df_rec.index[df_rec['protocol'] == 'vc'].tolist()
                index_char = df_rec.index[df_rec['protocol'] == 'freq analyse'].tolist()
                index_vm = df_rec.index[df_rec['protocol'] == 'vm_mouse'].tolist()

        #creating a dir to save plots (if not existing)
        dir_plots = "plots"
        path = os.path.join(work_dir, dir_plots)
        if os.path.isdir(path) == False: os.mkdir(path)

        #plotting the middle sweep of each filename
        for i in range(len(filenames)):
            filename = work_dir + filenames[i]
            tn.plot_traces(filename)

        print(index_vm)
        print(index_vc)
        print(index_char)
        print(slice_indx)
        if len(index_char) != len(index_vc): print('Fix protocol names. Unequal number of VC and freq analyse protocold')

        proceed_y_n = input("do all traces correspond to specified filenames in lab book (y/n)?")
        if proceed_y_n == 'y':
            dir_data = "data_tables"
            path = os.path.join(work_dir , dir_data)
            if os.path.isdir(path) == False: os.mkdir(path)

#%%


#sorting the files based on lab journal entry and number of sweeps
#plots and variables saved in indicated folders (see output)
df = pd.DataFrame(columns=['tissue_source','OP', 'patcher', 'patient_age', 
'filename','slice', 'cell_ch', 'day', 'treatment', 'hrs_incubation','cell_ID', 
'repatch', 'hrs_after_op','AP_heigth','Rheobase', 'TH', 'Vm', 
'capacitance', 'max_depol', 'max_repol', 'max_spikes','membra_time_constant_tau', 
'resting_potential', 'Rs', 'Rin', 'spon_freq', 'spon_ampl', 'mini_freq', 'mini_amp'])																	

#%%
# RUN for NO Vm file

for i in range(len(index_vc)):
    vc = index_vc[i]
    char = index_char[i]
    slice = slice_names[i][:2]

    filename_vc = work_dir+filenames[vc]
    filename_char = work_dir+filenames[char]
    #tn.plot_traces(filename_vc)

    #give inputs with space between entries
    active_channels = [int(item) for item in input('Which channels were active in' + filenames[vc] +'(look at plots)').split()]

    for j in range(len(active_channels)):
        ch = active_channels[j]
        hsf.get_holding_measures(filename_vc, ch)
        Rs, Rin = hcf.access_resistance(filename_vc, ch) #output ResRa, ResRi
        # cellID = filenames[vc][:-7]+slice+'c'+str(ch)
        # df = df.append({'filename':filenames[vc],'slice':slice,'cell_ch' :ch, 
        # 'cell_ID':cellID, 'Rs':Rs, 'Rin':Rin}, ignore_index=True)
        
        Vm, max_spikes, Rheobase, APheight, max_depol, max_repol, TH, capacitance, tau  = hcf.APproperties(filename_char, ch)
        cellID = filenames[char][:-7]+slice+'c'+str(ch)
        df = df.append({'tissue_source': tissue_source, 'OP':OP[:-1], 'patcher':patcher, 'patient_age':age, 
        'filename':filenames[char],'slice':slice, 'cell_ch':ch, 'cell_ID':cellID,
        'Rs':Rs, 'Rin':Rin, 'AP_heigth':APheight,'Rheobase':Rheobase, 'TH':TH, 'Vm' :Vm, 
        'capacitance':capacitance, 'max_depol':max_depol, 'max_repol':max_repol, 
        'max_spikes':max_spikes, 'membra_time_constant_tau':tau}, ignore_index=True)

df.to_excel(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.xlsx') 
# df.to_csv(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.csv')
# # %%

#%%
# RUN WHEN Vm file present

for i in range(len(index_vc)):
    vm = index_vm[i]
    vc = index_vc[i]
    char = index_char[i]
    slice = slice_names[i][:2]

    filename_vm = work_dir+filenames[vm]
    filename_vc = work_dir+filenames[vc]
    filename_char = work_dir+filenames[char]
    #tn.plot_traces(filename_vc)

    #give inputs with space between entries
    active_channels = [int(item) for item in input('Which channels were active in' + filenames[vc] +'(look at plots)').split()]

    for j in range(len(active_channels)):
        ch = active_channels[j]
        hsf.get_holding_measures(filename_vc, ch)
        Rs, Rin = hcf.access_resistance(filename_vc, ch) #output ResRa, ResRi
        # cellID = filenames[vc][:-7]+slice+'c'+str(ch)
        # df = df.append({'filename':filenames[vc],'slice':slice,'cell_ch' :ch, 
        # 'cell_ID':cellID, 'Rs':Rs, 'Rin':Rin}, ignore_index=True)
        
        resting_mem = hcf.restingmembrane(filename_vm, ch)
        Vm, max_spikes, Rheobase, APheight, max_depol, max_repol, TH, capacitance, tau  = hcf.APproperties(filename_char, ch)
        cellID = filenames[char][:-7]+slice+'c'+str(ch)
        df = df.append({'tissue_source': tissue_source, 'OP':OP[:-1], 'patcher':patcher, 'patient_age':age, 
        'filename':filenames[char],'slice':slice, 'cell_ch':ch, 'cell_ID':cellID,
        'Rs':Rs, 'Rin':Rin,'resting_potential': resting_mem,'AP_heigth':APheight,'Rheobase':Rheobase, 'TH':TH, 'Vm' :Vm, 
        'capacitance':capacitance, 'max_depol':max_depol, 'max_repol':max_repol, 
        'max_spikes':max_spikes, 'membra_time_constant_tau':tau}, ignore_index=True)

df.to_excel(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.xlsx') 
# df.to_csv(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.csv')
# %%

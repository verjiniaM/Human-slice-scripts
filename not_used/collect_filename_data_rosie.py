# %%
import os
import human_characterisation_functions as hcf
import human_synaptic_functions1 as hsf
import trace_names_check as tn
import pandas as pd

human_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_rosie/'

OP = 'OP210323/'
tissue_source = 'Virchow' #Bielefeld, Mitte, Virchow
patcher = 'Rosie' #Verji, Rosie
age = 'J' #A for adult, 'J' for juvenile

work_dir = human_folder + OP
file_list = sorted(os.listdir(work_dir))

# filenames = []
# for i in range(len(file_list)):
#     if file_list[i][-4:] == '.abf': 
#         filenames.append(file_list[i])
    
# for i in range(len(filenames)):
#     filename = work_dir + filenames[i]
#     tn.plot_traces(filename)

#%%

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

print(index_vc)
print(index_char)
print(slice_names)

# %%
dir_data = "data_tables"
path = os.path.join(work_dir , dir_data)
if os.path.isdir(path) == False: os.mkdir(path)

# %%
df = pd.DataFrame(columns=['tissue_source','OP', 'patcher', 'patient_age', 
'filename','slice', 'cell_ch', 'day', 'treatment', 'hrs_incubation','cell_ID', 
'repatch', 'hrs_after_op', 'repatch','AP_heigth','Rheobase', 'TH', 'Vm', 
'capacitance', 'max_depol', 'max_repol', 'max_spikes','membra_time_constant_tau', 
'resting_potential', 'Rs', 'Rin', 'spon_freq', 'spon_ampl', 'mini_freq', 'mini_amp'])

# %%
for i in range(len(index_vc)):
    vc = index_vc[i]
    char = index_char[i]
    slice = slice_names[i]

    filename_vc = work_dir+filenames[vc]
    filename_char = work_dir+filenames[char]

    active_channels = [int(item) for item in input('Which channels were active in' + filenames[vc] +'(look at plots)').split()]

    for j in range(len(active_channels)):
        ch = active_channels[j]
        hsf.get_holding_measures(filename_vc, ch)
        Rs, Rin = hcf.access_resistance(filename_vc, ch) #output ResRa, ResRi
        
        Vm, max_spikes, Rheobase, APheight, max_depol, max_repol, TH, capacitance, tau  = hcf.APproperties(filename_char, ch)
        cellID = filenames[char][:-7]+slice+'c'+str(ch)

        df = df.append({'tissue_source': tissue_source, 'OP':OP[:-1], 'patcher':patcher, 'patient_age':age, 
        'filename':filenames[char],'slice':slice, 'cell_ch':ch, 'cell_ID':cellID,
        'Rs':Rs, 'Rin':Rin,'AP_heigth':APheight,'Rheobase':Rheobase, 'TH':TH, 'Vm' :Vm, 
        'capacitance':capacitance, 'max_depol':max_depol, 'max_repol':max_repol, 
        'max_spikes':max_spikes, 'membra_time_constant_tau':tau}, ignore_index=True)


# %%
df

# %%
df.to_excel(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.xlsx') 


# %%

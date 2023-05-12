import os
import pandas as pd
import glob
import numpy as np
import json
import datetime
import human_characterisation_functions as hcf

#this function opens the alerady existing experiments_overview file and adds the latest OP info
#patcher and op_folder have to be strings; so with " "
# the op folder is the folder where all op folders are stored - *data_rosie or data_verji*
def update_op_list(human_dir, exp_view):
    '''
    opens the existing experiments_overview file; adds the latest OP and patcher
    '''
    folder_list_verji = np.array(sorted(os.listdir(human_dir + 'data_verji/'))) #lists all items in the op_folder
    folder_list_rosie = np.array(sorted(os.listdir(human_dir + 'data_rosie/')))
    folder_list = np.concatenate((folder_list_verji, folder_list_rosie))

    OP_folders_indx = [i for i in range(len(folder_list)) if folder_list[i][0:2] == "OP"] #index of all items that start with "OP"
    OP_folders = folder_list[OP_folders_indx].tolist()

    existing_OPs = exp_view['OP'].tolist()

    not_added = sorted(list(set(OP_folders) - set(existing_OPs)))

    patchers = []
    for x in not_added:
        patcher = input('Patcher for ' + x + '(Verji or Rosie): ')
        patchers.append(patcher)

    not_added_df = pd.DataFrame({"OP": not_added, "patcher": patchers})
    exp_view_new = pd.concat([exp_view, not_added_df],ignore_index=True).reset_index(drop=True) #maybe need to remove ignore_index = T

    date = str(datetime.date.today())
    exp_view_new.to_excel(human_dir + date + '_experiments_overview.xlsx', index = False) 
    exp_view.to_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/summary_data_tables/old_data/experiments_overview_before_'+ date + '.xlsx', index = False)

    return exp_view_new

def add_cortex_out_time(human_dir, exp_view):
    '''
    adds cortex_out time in experiments_overview table
    '''
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

    date = str(datetime.date.today())
    exp_view.to_excel(human_dir + date + '_experiments_overview.xlsx', index = False)

def get_sorted_file_list (dir):
    '''
    returns sorted file list
    '''
    file_list = sorted(os.listdir(dir))
    return file_list

def get_lab_book (OP_dir):
    '''
    returns a pandas data frame of the (.xlsx) file from the OP dir
    '''
    lab_book_path = glob.glob(OP_dir + '*.xlsx' )[0]
    if lab_book_path != '~$':
        df_rec = pd.read_excel(lab_book_path, header = 1)
    return df_rec

def get_abf_files (file_list):
    filenames = []
    for file in range(len(file_list)):
        #pclamp files
        if file_list[file][-4:] == '.abf': 
            filenames.append(file_list[file])
    return filenames

def sort_protocol_names (file_list, df_rec):
    '''
    takes a file list (contents of a folder) and pandas data frame (lab notebook)
    returns indices for all recording types
    '''

    slice_indx = df_rec.index[df_rec['slice'].notnull()]
    def_slice_names = df_rec['slice'][slice_indx].tolist()
    
    index_dict = {}
    dict_keys = ['vc', 'resting', 'freq analyse', 'ramp','con_screen', 'con_screen_VC',
    'spontan','vc_end', 'vc_mini', 'minis', 'vc_mini_end', 'resting long']
    for key in dict_keys:
        index = df_rec.index[df_rec['protocol'] == key].tolist()
        #df_rec['protocol'][index] = np.nan
        index_dict[key] = index
        

    ic_indices = []
    for i in range(len(df_rec['protocol'])):
        name = df_rec['protocol'][i]
        if type(name) == float:
            continue
        if name.find('IC') != -1 or name.find('ic') != -1:
            ic_indices.append(i)
    index_dict['IC_files'] = ic_indices

    #other_indices = df_rec.index[df_rec['protocol'].isnull() == 0].tolist()
    return slice_indx, def_slice_names, index_dict
#df_rec['protocol'][other_indices]

def fix_slice_names (def_slice_names, slice_indx):
    new_slice_names = []
    for i in range(len(def_slice_names)):
        if i < len(def_slice_names)-1:
            new_slice_names.append([def_slice_names[i]]*(slice_indx[i+1]-slice_indx[i]))
        else :
            new_slice_names.append([def_slice_names[i]]*(40))
    slice_names  = [x for xs in new_slice_names for x in xs]
    return slice_names

def get_work_dir(human_dir, OP, patcher):
    OP_folder = OP + '/'
    if patcher == 'Verji': 
        work_dir = human_dir + 'data_verji/'+ OP_folder
    else:
        work_dir = human_dir + 'data_rosie/'+ OP_folder 
    return work_dir

def get_OP_metadata (human_dir, OP, patcher):
    work_dir = get_work_dir(human_dir, OP, patcher)
    file_list = get_sorted_file_list(work_dir)
    df_rec = get_lab_book(work_dir)
    filenames = get_abf_files(file_list)
    slice_indx, def_slice_names, indices_dict = sort_protocol_names(file_list, df_rec)
    slice_names = fix_slice_names(def_slice_names, slice_indx)
    return work_dir, filenames, indices_dict, slice_names

def make_dir_if_not_existing(working_dir, new_dir):
    '''
    creates a dir new_dir in working_dir, if not already existing
    '''
    path = os.path.join(working_dir, new_dir)
    if os.path.isdir(path) == False: os.mkdir(path)
    return path

def to_json (work_dir, OP, file_out, out_data):
    '''
    returns a list with 3 dictionaries; 
    [1] slice names and active channels (for characterization)
    [0] con_screen_file indices and active channels
    [2] mini_meta - mini_file_indices and active chanels
    [3] treatment dictionary
    [4] connected cells in IC 
    '''
    fname = work_dir + OP + file_out
    with open(fname, "w") as outfile:
        json.dump(out_data , outfile)

def from_json (work_dir, OP, fn):
    fname = work_dir + OP + fn
    f = open(fname)
    return json.load(f)

def get_json_files (file_list):
    jsons = []
    for file in range(len(file_list)):
        #pclamp files
        if file_list[file][-5:] == '.json': 
            jsons.append(file_list[file])
    return jsons

def get_json_meta (human_dir, OP, patcher, file_out): # file_out = '_meta_active_chans.json'
    work_dir, filenames, indices_dict, slice_names = get_OP_metadata(human_dir, OP, patcher)
    exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0]) 

    file_list = get_sorted_file_list(work_dir)
    jsons = get_json_files(file_list)

    if OP + file_out in jsons:
        json_meta = from_json(work_dir, OP, file_out)
        return json_meta

    exp_view['cortex_out'] = exp_view['cortex_out'].dt.strftime('%Y-%m-%d %H:%M:%S')
    op_time = exp_view['cortex_out'][exp_view.index[exp_view['OP'] == OP]].tolist()
    #op_time = input('Cortex out [yyyy mm dd hh mm] ' + OP).split()

    treatment_dict = {}
    slices = np.unique(np.array(slice_names))
    unique_slices = []
    [unique_slices.append(S[:2]) for S in slices]
    unique_slices = np.unique(np.array(unique_slices))

    for i in unique_slices:
        treatment = input('Treatment (Ctrl, high K, or TTX) for ' + OP + ' ' + i)
        treatment_dict[i] = treatment
    
    active_chans_all, slice_names_dict, treatments = [], [], []
    for i in indices_dict['vc']:
        active_channels = [int(item) for item in input('Used channels in ' + OP + ' ' + filenames[i]).split()]
        active_chans_all.append(active_channels)
        slice_names_dict.append(slice_names[i])
        treatments.append(treatment_dict[slice_names[i][:2]])
    charact_meta = {
        'OP_time' : op_time,
        'slices' : slice_names_dict,
        'treatment' : treatments,
        'vc_files' : indices_dict['vc'],
        'active_chans': active_chans_all
    }

    pre_chans_all, post_chans_all, treatments = [], [], [] 
    for indx in indices_dict['con_screen']:
        #treatment = input('Treatment (Ctrl, high K, or TTX) for ' + OP + ' ' + slice_names[indx])
        pre_chans = [int(item) for item in input('Pre channels in ' + filenames[indx]).split()]
        post_chans = [int(item) for item in input('Post channels in ' + filenames[indx]).split()]
        pre_chans_all.append(pre_chans)
        post_chans_all.append(post_chans)
        treatments.append(treatment_dict[slice_names[indx][:2]])
    con_screen_meta = {
        'OP_time' : op_time,
        'con_screen_file_indices' : indices_dict['con_screen'],
        'treatment' : treatments,
        'pre_chans' : pre_chans_all,
        'post_chans' : post_chans_all
    }

    pre_chans_all_IC, post_chans_all_VC, treatments = [], [], []
    for i in indices_dict['IC_files']:
        pre_chans_IC = [int(item) for item in input('Pre channels in IC ' + filenames[i]).split()]
        post_chans_VC = [int(item) for item in input('Post channels in VC ' + filenames[i]).split()]
        pre_chans_all_IC.append(pre_chans_IC)
        post_chans_all_VC.append(post_chans_VC)
        treatments.append(treatment_dict[slice_names[i][:2]])
    con_screen_meta_IC = {
        'OP_time' : op_time,
        'con_screen_IC_file_indices' : indices_dict['IC_files'],
        'treatment' : treatments,
        'pre_chans_IC' : pre_chans_all_IC,
        'post_chans_VC' : post_chans_all_VC
    }

    chans_all, mini_slices, treatments = [],[], []
    for i in indices_dict['minis']:
        #treatment = input('Treatment (Ctrl, high K, or TTX) for ' + OP + ' ' + slice_names[i])
        active_channels = [int(item) for item in input('Used channels in ' + OP + ' ' + filenames[i]).split()]
        chans_all.append(active_channels)
        mini_slices.append(slice_names[i])
        treatments.append(treatment_dict[slice_names[i][:2]])
    mini_meta = {
        'OP_time' : op_time,
        'mini_slices' : mini_slices,
        'treatment' : treatments,
        'mini_files' : indices_dict['minis'],
        'mini_chans' : chans_all
    }
    
    to_json (work_dir, OP, file_out, [charact_meta, con_screen_meta, mini_meta, treatment_dict, con_screen_meta_IC])
    json_meta = from_json(work_dir, OP, file_out)
    return json_meta

def get_datetime_from_input (op_time):
    # year, month, day, hour, minute = map(int, op_time)
    # cortex_out_time = datetime.datetime(year, month, day, hour, minute)
    cortex_out_time = datetime.datetime.strptime(op_time, "%Y-%m-%d %H:%M:%S")
    return cortex_out_time

def get_time_after_OP (filename, cortex_out_time):

    rec_time = hcf.get_recording_time(filename)

    dt = rec_time - cortex_out_time
    h_after_op = dt.days*24 + dt.seconds/3600

    time_after_op = datetime.time(hour = int(dt.seconds /3600), minute = int((dt.seconds / 3600 - int(dt.seconds / 3600))*60))

    return h_after_op

def get_json_meta_high_K (human_dir, OP, patcher, file_out): # file_out = '_meta_active_chans.json'
    work_dir, filenames, indices_dict, slice_names = get_OP_metadata(human_dir, OP, patcher)
    exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0]) 

    file_list = get_sorted_file_list(work_dir)
    jsons = get_json_files(file_list)

    if OP + file_out in jsons:
        json_meta = from_json(work_dir, OP, file_out)
        return json_meta

    exp_view['cortex_out'] = exp_view['cortex_out'].dt.strftime('%Y-%m-%d %H:%M:%S')
    op_time = exp_view['cortex_out'][exp_view.index[exp_view['OP'] == OP]].tolist()

    active_chans_all, slice_names_dict = [], []
    for i in indices_dict['vc']:
        active_channels = [int(item) for item in input('Used channels in ' + OP + ' ' + filenames[i]).split()]
        active_chans_all.append(active_channels)
        slice_names_dict.append(slice_names[i])
    vc_meta = {
        'OP_time' : op_time,
        'slices' : slice_names_dict,
        'vc_files' : indices_dict['vc'],
        'active_chans': active_chans_all
    }

    active_chans_all, slice_names_dict = [], []
    for i in indices_dict['freq analyse']:
        active_channels = [int(item) for item in input('Used channels in ' + OP + ' ' + filenames[i]).split()]
        active_chans_all.append(active_channels)
        slice_names_dict.append(slice_names[i])
    charact_meta = {
        'OP_time' : op_time,
        'slices' : slice_names_dict,
        'char_files' : indices_dict['freq analyse'],
        'active_chans': active_chans_all
    }

    active_chans_all, slice_names_dict = [], []
    for i in indices_dict['resting']:
        active_channels = [int(item) for item in input('Used channels in ' + OP + ' ' + filenames[i]).split()]
        active_chans_all.append(active_channels)
        slice_names_dict.append(slice_names[i])
    resting_meta = {
        'slices' : slice_names_dict,
        'resting' : indices_dict['resting'],
        'active_chans': active_chans_all
    }

    active_chans_all, slice_names_dict = [], []
    for i in indices_dict['resting long']:
        active_channels = [int(item) for item in input('Used channels in ' + OP + ' ' + filenames[i]).split()]
        active_chans_all.append(active_channels)
        slice_names_dict.append(slice_names[i])
    resting_long_meta = {
        'slices' : slice_names_dict,
        'resting long' : indices_dict['resting long'],
        'active_chans': active_chans_all
    }

    to_json(work_dir, OP, file_out, [vc_meta, charact_meta, resting_meta, resting_long_meta])
    json_meta = from_json(work_dir, OP, file_out)
    return json_meta
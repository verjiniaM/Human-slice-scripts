import os
import pandas as pd
import glob
import numpy as np
import json
import datetime
import funcs_human_characterisation as hcf
import funcs_plotting as plotting_funcs
import pyabf
import re

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

def get_sorted_file_list(dir):
    '''
    returns sorted file list
    '''
    file_list = sorted(os.listdir(dir))
    return file_list

def get_lab_book(OP_dir):
    '''
    returns a pandas data frame of the (.xlsx) file from the OP dir
    '''
    lab_book_path = glob.glob(OP_dir + '*.xlsx' )[0]
    if lab_book_path != '~$':
        df_rec = pd.read_excel(lab_book_path, header = 1)
    return df_rec, lab_book_path

def get_abf_files(file_list):
    filenames = []
    for file in range(len(file_list)):
        #pclamp files
        if file_list[file][-4:] == '.abf': 
            filenames.append(file_list[file])
    return filenames

def sort_protocol_names (file_list, df_rec, lab_book_path, V = 'new'):
    '''
    takes a file list (contents of a folder) and pandas data frame (lab notebook)
    returns indices for all recording types
    '''

    slice_indx = df_rec.index[df_rec['slice'].notnull()]
    def_slice_names = df_rec['slice'][slice_indx].tolist()
    
    index_dict = {}
    dict_keys = ['vc', 'resting', 'freq analyse', 'characterization' ,'ramp',
    'con_screen', 'con_screen_VC',
    'spontan','vc_end', 'vc_mini', 'minis', 'vc_mini_end', 'resting_long']
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

    rec_file = pd.ExcelFile(lab_book_path)
    num_slices = rec_file.sheet_names[1:]
    
    if V == 'old':
        return slice_indx, def_slice_names, index_dict

    pre_chans, post_chans = [], []
    for slice_ in num_slices:
        con_screen_matrix = pd.read_excel(lab_book_path, slice_, index_col=0)

        pre_cells, post_cells = [], []
        for col in con_screen_matrix.columns:
            post_ = con_screen_matrix[col][con_screen_matrix[col] == 1].index.tolist()

            if len(post_) > 0:
                for i in post_:
                    pre_cells.append(i[-1])
                    post_cells.append(col[-1])
        pre_chans.append(pre_cells)
        post_chans.append(post_cells)
    #other_indices = df_rec.index[df_rec['protocol'].isnull() == 0].tolist()
    return slice_indx, def_slice_names, index_dict, pre_chans, post_chans
#df_rec['protocol'][other_indices]

def fix_slice_names (def_slice_names, slice_indx):
    '''
    makes a continuous list of the slice names
    indexing it with a file index shows the slices this files belongs to
    '''
    new_slice_names = []
    for i in range(len(def_slice_names)):
        if i < len(def_slice_names)-1:
            new_slice_names.append([def_slice_names[i]]*(slice_indx[i+1]-slice_indx[i]))
        else :
            new_slice_names.append([def_slice_names[i]]*(40))
    slice_names  = [x for xs in new_slice_names for x in xs]
    return slice_names

def get_work_dir(human_dir, OP, patcher: str):
    OP_folder = OP + '/'
    if patcher == 'Verji': 
        work_dir = human_dir + 'data_verji/'+ OP_folder
    elif patcher == '':
        work_dir = human_dir + OP_folder
    else:
        work_dir = human_dir + 'data_rosie/'+ OP_folder 
    return work_dir

def get_OP_metadata (human_dir, OP, patcher, V = 'new'):
    work_dir = get_work_dir(human_dir, OP, patcher)
    file_list = get_sorted_file_list(work_dir)
    jsons = get_json_files(file_list)
    df_rec, lab_book_path = get_lab_book(work_dir)
    filenames = get_abf_files(file_list)

    if V == 'old':
        slice_indx, def_slice_names, indices_dict = sort_protocol_names(file_list, df_rec, lab_book_path, V = 'old')
    else:
        slice_indx, def_slice_names, indices_dict, pre_chans, post_chans = sort_protocol_names(file_list, df_rec, lab_book_path)

    if OP + '_indices_dict.json' in jsons:
        indices_dict = from_json(work_dir, OP, '_indices_dict.json')
    else: 
        to_json(work_dir, OP, '_indices_dict.json', indices_dict)

    slice_names = fix_slice_names(def_slice_names, slice_indx)
    
    if V == 'old':
        return work_dir, filenames, indices_dict, slice_names

    return work_dir, filenames, indices_dict, slice_names, pre_chans, post_chans


def make_dir_if_not_existing(working_dir, new_dir):
    '''
    creates a dir new_dir in working_dir, if not already existing
    '''
    path = os.path.join(working_dir, new_dir)
    if os.path.isdir(path) == False: os.mkdir(path)
    return path

def plot_trace_if_not_done(work_dir, dir_plots, filenames):
    '''
    check if the traces dir is empty and
    only then plot the middle sweep for each filename
    '''
    traces_folder =  os.path.join(dir_plots, "traces/")
    if os.path.isdir(traces_folder) == 0 :
        for rec in range(len(filenames)):
            filename = work_dir + filenames[rec]
            plotting_funcs.plot_middle_sweep(filename)
    else:
            print("skipping plotting")

def to_json(work_dir, OP, file_out, out_data):
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

def open_json(dir):
    f = open(dir)
    return json.load(f)

def get_json_files (file_list):
    jsons = []
    for file in range(len(file_list)):
        #pclamp files
        if file_list[file][-5:] == '.json': 
            jsons.append(file_list[file])
    return jsons

def get_json_meta (human_dir, OP, patcher, file_out): # file_out = '_meta_active_chans.json'
    work_dir, filenames, indices_dict, slice_names, pre_chans, post_chans = get_OP_metadata(human_dir, OP, patcher)
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

    treatments = []
    for indx in indices_dict['con_screen']:
        #treatment = input('Treatment (Ctrl, high K, or TTX) for ' + OP + ' ' + slice_names[indx])
        # pre_chans = [int(item) for item in input('Pre channels in ' + filenames[indx]).split()]
        # post_chans = [int(item) for item in input('Post channels in ' + filenames[indx]).split()]
        # pre_chans_all.append(pre_chans)
        # post_chans_all.append(post_chans)
        treatments.append(treatment_dict[slice_names[indx][:2]])
    con_screen_meta = {
        'OP_time' : op_time,
        'con_screen_file_indices' : indices_dict['con_screen'],
        'treatment' : treatments,
        'pre_chans' : pre_chans,
        'post_chans' : post_chans
    }

    pre_chans_all_IC, post_chans_all_VC, treatments = [], [], []
    for i in indices_dict['IC_files']:
        # pre_chans_IC = [int(item) for item in input('Pre channels in IC ' + filenames[i]).split()]
        # post_chans_VC = [int(item) for item in input('Post channels in VC ' + filenames[i]).split()]
        pre_chans_IC = []
        post_chans_VC = []
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


def get_json_meta_connect_all (human_dir, OP, patcher, file_out = '_con_screen_only.json'): # file_out = '_con_screen_only.json'
    work_dir, filenames, indices_dict, slice_names = get_OP_metadata(human_dir, OP, patcher, 'old')
    
    if isinstance(indices_dict, list):
        indices_dict = indices_dict[0]

    exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0]) 
    op_time = exp_view['cortex_out'][exp_view.index[exp_view['OP'] == OP]].tolist()
    file_list = get_sorted_file_list(work_dir)
    
    jsons = get_json_files(file_list)   
    if OP + file_out in jsons:
        json_con_screen = from_json(work_dir, OP, file_out)
        if len(json_con_screen) == 1:

            work_dir, filenames, indices_dict, slice_names, pre_chans, post_chans = get_OP_metadata(human_dir, OP, patcher)

            treatment_dict = {}
            for a, s in enumerate(json_con_screen[0]['slices']):
                treatment_dict[s] = json_con_screen[0]['treatment'][a]

            if isinstance(indices_dict, list):
                indices_dict = indices_dict[0]
            json_con_screen[0]['pre_chans'] = pre_chans
            json_con_screen[0]['post_chans'] = post_chans

            if 'IC_files' in indices_dict.keys():
                pre_chans_all_IC, post_chans_all_VC, treatments = [], [], []
                df_rec, lab_book_path = get_lab_book(work_dir)
                for a, j in enumerate(indices_dict['IC_files']):
                    prot_name = df_rec['protocol'][j]
                    pre_chan_IC = int(re.findall(r'\d', prot_name)[0])
                    if slice_names[j] not in json_con_screen[0]['slices']:
                        post_chans_VC = [int(item) for item in input('Post channels in ' + filenames[j]).split()]
                    else:
                        indx_ = json_con_screen[0]['slices'].index(slice_names[j])

                        pre_chans_all = [eval(i) for i in pre_chans[indx_]]
                        post_chans_all = [eval(i) for i in post_chans[indx_]]

                        pre_chan_indx = [i for i, val in enumerate(pre_chans_all) if val==pre_chan_IC]
                        post_chans_VC = []
                        for i in pre_chan_indx:
                            post_chans_VC.append(post_chans_all[i])
        
                    pre_chans_IC = [pre_chan_IC] * len(post_chans_VC)

                    pre_chans_all_IC.append(pre_chans_IC)
                    post_chans_all_VC.append(post_chans_VC)
                    treatments.append(treatment_dict[slice_names[j][:2]])

                con_screen_meta_IC = {
                    'con_screen_IC_file_indices' : indices_dict['IC_files'],
                    'treatment' : treatments,
                    'pre_chans_IC' : pre_chans_all_IC,
                    'post_chans_VC' : post_chans_all_VC
                }
                to_json(work_dir, OP, file_out, [json_con_screen[0], con_screen_meta_IC])
                json_con_screen = from_json(work_dir, OP, file_out)
            else :
                to_json(work_dir, OP, file_out, [json_con_screen[0]])
                json_con_screen = from_json(work_dir, OP, file_out)
        return json_con_screen

    treatment_dict = {}
    slices = np.unique(np.array(slice_names))
    unique_slices = []
    [unique_slices.append(S[:2]) for S in slices]
    unique_slices = np.unique(np.array(unique_slices))

    for i in unique_slices:
        treatment = input('Treatment (Ctrl, high K, or TTX) for ' + OP + ' ' + i)
        treatment_dict[i] = treatment

    active_chans_con_screen, slice_names_dict, treatments = [], [], []
    for i in indices_dict['con_screen']:
        slice_names_dict.append(slice_names[i])
        treatments.append(treatment_dict[slice_names[i][:2]])
        active_channels = [int(item) for item in input('Used channels in ' + OP + ' ' + filenames[i]).split()]
        active_chans_con_screen.append(active_channels)
        con_screen_meta = {
            'OP_time' : op_time[0].strftime('%y%m%d %H:%M'),
            'slices' : slice_names_dict,
            'treatment' : treatments,
            'con_screen_file' : indices_dict['con_screen'],
            'active_chans': active_chans_con_screen
        }
    to_json(work_dir, OP, file_out, [con_screen_meta])
    json_con_screen = from_json(work_dir, OP, file_out)
    return json_con_screen


def get_datetime_from_input (op_time):
    ''' 
    op_time - 2020-10-20 10:00:00
    '''
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


#########
# used for EPSP manual analysis

def concat_dfs_in_folder(folder_path):
    '''
    reads all the excels files from the folder 
    example input(str): verji_dir + 'OP*' + '/data_tables/'  + '*RMP_high_K' + '*.xlsx'
    output: pandas datafframe with all excels concatinated

    '''
    if folder_path[-5:] == '.xlsx':
        df_paths = sorted(glob.glob(folder_path))
    else:
        df_paths = sorted(glob.glob(folder_path + '*.xlsx'))
    df_combined = pd.DataFrame()
    for path in df_paths:
        df = pd.read_excel(path)
        df_combined = pd.concat([df_combined[:], df]).reset_index(drop=True)
    return df_combined

def create_full_results_df(meta_keep_path, results_QC_path, verji_dir):
    '''
    takes one mate df where all files need to be kept
    takes a results df with summary sheet and a results sheet for each recording
    collects all RMP_high_K datatables from verji_dir
    adds metadata from them to the results_keep dataframe
    
    '''

    meta_keep = pd.read_excel(meta_keep_path)

    meta_keep_IDs = []
    for i in range(len(meta_keep)):
        meta_keep_IDs.append(meta_keep['Name of recording'][i][:-4] + '_' + str(meta_keep['Channels to use'][i]))

    #results
    results_df = pd.read_excel(results_QC_path + 'complete_QC_results.xlsx')
    results_file = pd.ExcelFile(results_QC_path + 'complete_QC_results.xlsx')

    results_IDs = []
    for j in range(len(results_df)):
        results_IDs.append(results_df['Recording filename'][j][:-4] + '_' + str(results_df['Channel'][j]))

    missing_file_lines = list(set(results_file.sheet_names) - set(results_IDs))
    in_meta_not_in_results = list(set(meta_keep_IDs) - set(results_IDs))
    results_exclude = list(set(results_IDs) - set(meta_keep_IDs))

    #if all of the above 3 lists are empy continue with a piece of heart

    #take more metadata from the tables in the data folders of the recordings
    RMP_high_K_all = concat_dfs_in_folder(verji_dir + 'OP*' + '/data_tables/'  + '*RMP_high_K' + '*.xlsx')

    sweeps, OPs, hrs_incub, cell_ID, recording_time, resting_potential, holding_minus_70_y_o_n, incubation_solution, \
        recording_in, tissue_source, patient_age, K_concentration,  temp = [],[],[], [], [], [], [], [], [], [], [], [], []
    #adding metadata columns 
    for i, fn in enumerate(results_df['Recording filename']):
        chan = results_df['Channel'][i]
        dat_to_add = RMP_high_K_all[(RMP_high_K_all['filename'] == fn) & 
                    (RMP_high_K_all['cell_ch'] == chan)]
        if len(dat_to_add) == 0:
            results_df =  results_df.drop(i)
            continue
        meta_keep_data =  meta_keep['swps_to_analyse'][(meta_keep['Name of recording'] == fn) & \
             (meta_keep['Channels to use'] == chan)].item()
        
        sweeps.append(meta_keep_data)
        OPs.append(dat_to_add['OP'].item())
        hrs_incub.append(dat_to_add['hrs_incubation'].item())
        cell_ID.append(dat_to_add['cell_ID'].item())
        recording_time.append(dat_to_add['recording_time'].item())
        list_restings = dat_to_add['resting_potential'].item()[1:-1].split(',') 
        RMP = np.mean([float(i) for i in list_restings ])
        resting_potential.append(RMP)
        holding_minus_70_y_o_n.append(dat_to_add['holding_minus_70_y_o_n'].item())
        incubation_solution.append(dat_to_add['incubation_solution'].item())
        recording_in.append(dat_to_add['recording_in'].item())
        tissue_source.append(dat_to_add['tissue_source'].item())
        patient_age.append(dat_to_add['patient_age'].item())
        K_concentration.append(dat_to_add['K concentration'].item())
        temp.append(dat_to_add['temperature'].item())

    list_of_lists = [sweeps, OPs, hrs_incub, cell_ID, recording_time, resting_potential, holding_minus_70_y_o_n, incubation_solution, \
        recording_in, tissue_source, patient_age, K_concentration, temp]
    list_col_names = ['analyzed_swps' ,'OP', 'hrs_incub', 'cell_ID', 'recording_time', 'resting_potential', 'holding_minus_70_y_o_n', 'incubation_solution', \
        'recording_in', 'tissue_source', 'patient_age', 'K_concentration', 'temperature']

    #columns_to_include = list(set(list(RMP_high_K_all.columns)) - set(list(results_df.columns)))
    for i, col in enumerate(list_of_lists):
        results_df.insert(i, list_col_names[i], col)
    
    return results_df


def get_h_m_s_from_timedelta(time_delta_obj):
    h = int(time_delta_obj.seconds/3600)
    m = int((time_delta_obj.seconds/3600 - int(time_delta_obj.seconds/3600))*60)
    s = time_delta_obj.seconds - h*3600 - m*60
    if s == 60:
        s = 0
        m = m+1
    return h, m, s


def get_time_of_recording_sweeps(results_df):
    '''
    takes the results df from the manually analyzed results
    adds columne with infor for the time  of the recordings
    relative to the starting time of recording for each slice
    and for the time withing the recording
    '''

    verji_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/'
    vc_times = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/EPSPs/vc_files_times.xlsx')

    slice_start_times_list, swp1_starts, swp2_starts, times_after_rec_slice_start_swp1, times_after_rec_slice_start_swp2 = [], [], [], [], []
    for i in range(len(results_df)):
        fn = verji_dir + results_df['OP'][i] + '/' + results_df['Recording filename'][i]
        rec_time = hcf.get_recording_time(fn)

        swps = [int(x) for x in results_df['analyzed_swps'][i][1:-1].split(',')]

        abf_fn = pyabf.ABF(fn)
        rec_time_time = datetime.timedelta(days = 0, seconds = (rec_time.hour*3600 + rec_time.minute*60 + rec_time.second))
        swp_time1 =  datetime.timedelta(days = 0, seconds = abf_fn.sweepLengthSec * (swps[0] - 1))
        swp_time2 =  datetime.timedelta(days = 0, seconds = abf_fn.sweepLengthSec * (swps[1] - 1))

        # h1, m1, s1 = get_h_m_s_from_timedelta(rec_time_time + swp_time1)
        # h2, m2, s2 = get_h_m_s_from_timedelta(rec_time_time + swp_time2)

        swp1_starts.append(rec_time_time.seconds + swp_time1.seconds)
        swp2_starts.append(rec_time_time.seconds + swp_time2.seconds)

        fn_results =  results_df['Recording filename'][i]
        for j, fn in enumerate(vc_times['fn'][vc_times['OP'] == results_df['OP'][i]]):
            if int(fn[5:8]) < int(fn_results[5:8])  and j + 1 == len(vc_times['fn'][vc_times['OP'] == results_df['OP'][i]]):
                slice_start_times_list.append(vc_times['rec_time'][vc_times['fn'] == fn].item())
            elif int(fn[5:8]) < int(fn_results[5:8]) and \
                int(vc_times['fn'][vc_times['OP'] == results_df['OP'][i]].tolist()[j+1][5:8]) > int(fn_results[5:8]): 
                slice_start_times_list.append(vc_times['rec_time'][vc_times['fn'] == fn].item())

        rec_timedelta = datetime.timedelta(days = 0, seconds = slice_start_times_list[i].time().hour*3600 + \
            slice_start_times_list[i].time().minute*60 + slice_start_times_list[i].time().second) 
        
        times_after_rec_slice_start_swp1.append(swp1_starts[i] - rec_timedelta.seconds)
        times_after_rec_slice_start_swp2.append(swp2_starts[i] - rec_timedelta.seconds)
        
    results_df.insert(len(results_df.columns), 'sweep1_time_within_rec (sec)', swp1_starts)
    results_df.insert(len(results_df.columns), 'sweep2_time_within_rec (sec)', swp2_starts)
    results_df.insert(len(results_df.columns), 'slice_rec_start', slice_start_times_list)
    results_df.insert(len(results_df.columns), 'time_after_rec_slice_start_swp1 (sec)', times_after_rec_slice_start_swp1)
    results_df.insert(len(results_df.columns), 'time_after_rec_slice_start_swp2 (sec)', times_after_rec_slice_start_swp2)
    #results_df.to_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/meta_events/EPSPs/DEL_240305-results_df_complete.xlsx')
    return results_df 

#%%

# to get the vc starting time - starting time of a recording 

# human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
# vc_times_all = pd.DataFrame(columns=['fn', 'rec_time'])
# for OP in results_df.OP.unique():

#     work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, 'Verji')

#     vc_times = {}
#     for j in indices_dict['vc']:
#         fn = human_dir + 'data_verji/' + OP + '/' + filenames[j]
#         vc_times[filenames[j]] = [hcf.get_recording_time(fn)]
#         vc_times_df = pd.DataFrame.from_dict(vc_times).T
#     vc_times_all = pd.concat([vc_times_all, vc_times_df], axis = 1)
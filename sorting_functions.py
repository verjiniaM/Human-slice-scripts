import os
import pandas as pd
import glob
import numpy as np
import json

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

def get_json_files (file_list):
    jsons = []
    for file in range(len(file_list)):
        #pclamp files
        if file_list[file][-5:] == '.json': 
            jsons.append(file_list[file])
    return jsons

def sort_protocol_names (file_list, df_rec):
    '''
    takes a file list (contents of a folder) and pandas data frame (lab notebook)
    returns indices for all recording types
    '''

    slice_indx = df_rec.index[df_rec['slice'].notnull()]
    def_slice_names = df_rec['slice'][slice_indx].tolist()
    
    index_dict = {}
    dict_keys = ['vc', 'vm_mouse', 'freq analyse', 'vc_end', 'vm', 'minis', 'con_screen', 'resting',]
    for key in dict_keys:
        index = df_rec.index[df_rec['protocol'] == key].tolist()
        #df_rec['protocol'][index] = np.nan
        index_dict[key] = index

    #other_indices = df_rec.index[df_rec['protocol'].isnull() == 0].tolist()
    return slice_indx, def_slice_names, index_dict
#df_rec['protocol'][other_indices]
def fix_slice_names (def_slice_names, slice_indx):
    new_slice_names = []
    for i in range(len(def_slice_names)):
        if i < len(def_slice_names)-1:
            new_slice_names.append([def_slice_names[i]]*(slice_indx[i+1]-slice_indx[i]))
        else :
            new_slice_names.append([def_slice_names[i]]*(15))
    slice_names  = [x for xs in new_slice_names for x in xs]
    return slice_names

def make_dir_if_not_existing(working_dir, new_dir):
    '''
    creates a dir new_dir in working_dir, if not already existing
    '''
    path = os.path.join(working_dir, new_dir)
    if os.path.isdir(path) == False: os.mkdir(path)
    return path


def to_json (work_dir, OP, fn, file_indices, pre_chans, post_chans, slices, active_chans):
    fname = work_dir + OP + fn
    dict1 = {
        'file_indices' : file_indices,
        'pre_chans' : pre_chans,
        'post_chans' : post_chans
    }
    dict2 = {
        'slices' : slices,
        'active_chans': active_chans
    }
    with open(fname, "w") as outfile:
        json.dump([dict1, dict2] , outfile)

def from_json (work_dir, OP, fn):
    fname = work_dir + OP + fn
    f = open(fname)
    return json.load(f)


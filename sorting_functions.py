import os
import pandas as pd
import glob

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
    slice_names = df_rec['slice'][slice_indx].tolist()
    
    index_dict = {}
    dict_keys = ['vc', 'vm_mouse', 'freq analyse', 'vc_end', 'vm', 'minis']
    for key in dict_keys:
        index = df_rec.index[df_rec['protocol'] == key].tolist()
        index_dict[key] = index

    # index_vc = df_rec.index[df_rec['protocol'] == 'vc'].tolist()
    # index_vm = df_rec.index[df_rec['protocol'] == 'vm_mouse'].tolist()
    # index_char = df_rec.index[df_rec['protocol'] == 'freq analyse'].tolist()
    # index_vc_end = df_rec.index[df_rec['protocol'] == 'vc_end'].tolist()
    # index_spontan = df_rec.index[df_rec['protocol'] == 'vm'].tolist()
    # index_mini = df_rec.index[df_rec['protocol'] == 'minis'].tolist()

    return slice_indx, def_slice_names, index_dict

    def fix_slice_names (def_slice_names, slice_indx):
        new_slice_names = []
        for i in def_slice_names:
            if i < len(def_slice_names)-1:
                new_slice_names.append([def_slice_names[i]]*(slice_indx[i+1]-slice_indx[i]))
            else :
                new_slice_names.append([def_slice_names[i]]*(15))
        slice_names  = [x for xs in new_slice_names for x in xs]
        return slice_names

    def make_dir_if_not_existing (working_dir, new_dir):
        '''
        creates a dir new_dir in working_dir, if not already existing
        '''
        path = os.path.join(working_dir, new_dir)
        if os.path.isdir(path) == False: os.mkdir(path)
        return path
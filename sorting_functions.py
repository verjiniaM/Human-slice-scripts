import os
import pandas as pd

def get_lab_book (OP_dir):
'''
returns a pandas data frame of the (.xlsx) file from the OP folder
and the sorted file list
'''
file_list = sorted(os.listdir(OP_dir))
if (file_list[file][-5:] == '.xlsx' and file_list[file][:2] != '~$'):
    df_rec = pd.read_excel(OP_dir + file_list[file], header = 1)
return df_rec, file_list

def get_abf_files (file_list):
    filenames = []
    for file in file_list:
        #pclamp files
        if file_list[file][-4:] == '.abf': 
            filenames.append(file_list[file])
    return filenames

def sort_protocol_names (file_list, df_rec):
    '''
    takes a file list (contents of a folder) and pandas data frame (lab notebook)
    returns indices for all recording types
    '''
    if (file_list[file][-5:] == '.xlsx' and file_list[file][:2] != '~$'): 
        slice_indx = df_rec.index[df_rec['slice'].notnull()]
        slice_names = df_rec['slice'][slice_indx].tolist()
        index_vc = df_rec.index[df_rec['protocol'] == 'vc'].tolist()
        index_vm = df_rec.index[df_rec['protocol'] == 'vm_mouse'].tolist()
        index_char = df_rec.index[df_rec['protocol'] == 'freq analyse'].tolist()
        index_vc_end = df_rec.index[df_rec['protocol'] == 'vc_end'].tolist()
        index_spontan = df_rec.index[df_rec['protocol'] == 'vm'].tolist()
        index_mini = df_rec.index[df_rec['protocol'] == 'minis'].tolist()

    return slice_indx, def_slice_names, index_vc, index_vm, index_char, index_vc_end, index_spontan, index_mini

    def fix_slice_names (def_slice_names, slice_indx):
        new_slice_names = []
        for i in def_slice_names:
            if i < len(def_slice_names)-1:
                new_slice_names.append([def_slice_names[i]]*(slice_indx[i+1]-slice_indx[i]))
            else :
                new_slice_names.append([def_slice_names[i]]*(15))
        slice_names  = [x for xs in new_slice_names for x in xs]
        return slice_names
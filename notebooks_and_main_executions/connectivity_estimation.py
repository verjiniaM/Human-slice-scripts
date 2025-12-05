import glob
import ephys_analysis.funcs_sorting as sort
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
import math
import os
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr
import ephys_analysis.funcs_con_screen as con_param
import ephys_analysis.funcs_human_characterisation as hcf
import ephys_analysis.funcs_plotting_raw_traces as funcs_plotting_raw_traces
import ephys_analysis.funcs_plotting_raw_traces as plotting_funcs
import ephys_analysis.funcs_for_results_tables as funcs_results

# %%
## FOR MONOGRAPH
def collect_con_screen_from_jsons():
    '''
    for every op in exp_dir, collects the data from the con_screen .json file
    '''
    human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
    exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0])

    # ops to skip because no data was collected
    ops_skip = ['OP230412', 'OP240124', 'OP240229', 'OP240318', 'OP240820', 'OP241203']

    df_con_screen_all = pd.DataFrame()
    ops_no_con_screen, op_no_pre_chans = [], []
    for i, op in enumerate(exp_view.OP):
        if op in ops_skip:
            continue
        patcher = exp_view.patcher[i]

        work_dir, _, _, _ = sort.get_OP_metadata(human_dir, op, patcher, 'old')
        try:
            con_screen_json = sort.from_json(work_dir, op, '_con_screen_only.json')
        except FileNotFoundError:
            print(f'no con_screen file for {op}, {patcher}')
            ops_no_con_screen.append(op)
            continue

        if 'pre_chans' in con_screen_json[0].keys():
            data_rows = []
            for i, pre_cells_list in enumerate(con_screen_json[0]['pre_chans']):
                post_cells_list = con_screen_json[0]['post_chans'][i]
                slice_name = con_screen_json[0]['slices'][i]
                treatment = con_screen_json[0]['treatment'][i]

                for j, pre_chan in enumerate(pre_cells_list):
                    post_chan = post_cells_list[j]

                    data_rows.append({
                        'OP': op,
                        'patcher': patcher, 
                        'slice': slice_name,
                        'treatment': treatment,
                        'pre_chan': pre_chan,
                        'post_chan': post_chan
                    })
            if len(df_con_screen_all) == 0:
                df_con_screen_all = pd.DataFrame(data_rows)
            else:
                df_con_screen_all = pd.concat([df_con_screen_all.loc[:], pd.DataFrame(data_rows)]).reset_index(drop=True)
        else:
            op_no_pre_chans.append(op)
    return df_con_screen_all,  ops_no_con_screen

def load_intr_adult_repatch_no_QC():
    '''
    collects the intrinsic dfs 

    type (str): 'all', 'Ctrl', 'high K'
    '''
    df_intr_complete = funcs_results.collect_intrinsic_df() # for most updated version
    # df_intr = pl_intr.get_column_RMPs_from_char(df_intr_complete)

    adult_df = pl_intr.filter_adult_hrs_incubation_data(df_intr_complete, min_age = 12, \
                                                        hrs_inc = 16, max_age = 151, QC = False) # + QC criteria
    adult_df_repatch = adult_df[adult_df.repatch == 'yes']

    return adult_df_repatch

def for_manual_check_and_edit(df_con_screen, df_repatch_no_QC):
    '''
    prints out if any slices have non-conventional names, ending in 'er' or 'd2'
    returns all op slice combinations that can't be found inn the collected repatched
    cells for intrinsic properties
    '''
    # manual check was completed here for all missing combinations
    for i, s in enumerate(df_con_screen.slice):
        if s[-2:] == 'd2' or s[-2:] == 'er':
            print(df_con_screen.iloc[i,:])

    # cehck for OP slice combinations, to exclude optential wrong namings of slices
    op_sl_con, op_sl_repatch = [], []
    for i, op in enumerate(df_con_screen.OP):
        op_sl_con.append(op + df_con_screen.slice[i])

    for i, op in enumerate(df_repatch_no_QC.OP):
        slic = df_repatch_no_QC.slice[i]
        op_sl_repatch.append(op + df_repatch_no_QC.slice[i])

    # manually check all missing combinations
    # some disappearing connections to be seen - TTX, OP220117, S4
    op_slice_not_in_repatch = sorted(set(op_sl_repatch) - set(op_sl_con))
    return op_slice_not_in_repatch

def create_unique_cell_IDs(df, data_type):
    '''
    if connectivity in data_type, creates also a connection ID
    '''
    # creating match cell_IDs
    pre_chan_ID, post_chan_ID, con_IDs = [], [], []
    if 'pre_chan' in df.columns:
        pre_chan = 'pre_chan'
        post_chan = 'post_chan'
    elif 'cell_ch' in df.columns:
        pre_chan = 'cell_ch'
    elif 'chan_pre' in df.columns:
        pre_chan = 'chan_pre'
        post_chan = 'chan_post'

    for i, pre_cell in enumerate(df[pre_chan]):
        op = df.OP[i][2:]
        slic = df.slice[i]
        patcher = df.patcher[i]
        pre_chan_ID.append(f'{patcher}_{op}_{slic}_{pre_cell}')
        if 'con' in data_type:
            post_cell = df[post_chan][i]
            post_chan_ID.append(f'{patcher}_{op}_{slic}_{post_cell}')
            con_IDs.append(f'{patcher}_{op}_{slic}_{pre_cell}_{post_cell}')
    if 'con' in data_type:
        df.insert(len(df.columns), 'pre_cell_ID', pre_chan_ID)
        df.insert(len(df.columns), 'post_cell_ID', post_chan_ID)
        df.insert(len(df.columns), 'con_ID', con_IDs)
        print(f'number of non-unique connections {data_type}: {len(df) - len(df['con_ID'].unique())}')
    else:
        df.insert(len(df.columns), 'cell_ID_match', pre_chan_ID)
        print(f'number of non-unique cells {data_type}: {len(df) - len(df['cell_ID_match'].unique())}')
    return df

def get_con_IDs_remove_manually_checked():
    return  ['Verji_221024_S3_4_7',
            'Verji_221024_S3_D2_4_7',
            'Verji_230208_S2_7_2',
            'Verji_230208_S2_D2_7_2',
            'Verji_230523_S2_1_2',
            'Verji_230523_S2_2_1_2',
            'Verji_230523_S2_2_3_1',
            'Verji_230523_S2_2_7_3',
            'Verji_230523_S2_2_7_8',
            'Verji_230523_S2_2_8_7',
            'Verji_230523_S2_6_1',
            'Verji_230523_S2_7_6',
            'Verji_230523_S2_7_8',
            'Verji_230523_S2_8_7',
            'Verji_230523_S3_2_8_3',
            'Verji_230523_S3_8_4',
            'Verji_231109_S2_4_8',
            'Verji_231109_S2_D2_4_8',
            'Verji_240411_S1_5_3',
            'Verji_240411_S1_6_5',
            'Verji_240411_S1_7_5',
            'Verji_240411_S1_D2_2_4',
            'Verji_240411_S1_D2_3_4',
            'Verji_240411_S1_D2_4_1',
            'Verji_240411_S2_3_4',
            'Verji_240411_S2_D2_3_4',
            'Verji_240417_S1_5_7',
            'Verji_240417_S1_D2_5_7',
            'Verji_240507_S3_7_5',
            'Verji_240507_S3_8_7',
            'Verji_240507_S3_D2_6_5',
            'Verji_240507_S3_D2_8_6']

def get_connection_IDs_remove_manually_checked_2():
    return ['23420S2c2#6',
            '23420S2c4#6',
            '23808S2c4#7',
            '24411S1c6#3',
            '24411S1c7#3',
            '24417S1c7#6',
            '24417S2c1#5']

def check_which_to_analyse(human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    # Lists to store folders without connection_analysis and with empty connection_analysis
    no_connection_analysis = []
    empty_connection_analysis = []

    # Iterate through the folders in human_dir
    for root, dirs, files in os.walk(human_dir):
        # Get the folder name
        folder_name = os.path.basename(root)
        
        # Check if the folder name starts with 'OP'
        if folder_name.startswith('OP'):
            # Check if 'connection_analysis' is not in the list of directories
            if 'connection_analysis' not in dirs:
                no_connection_analysis.append(folder_name)
            else:
                connection_analysis_path = os.path.join(root, 'connection_analysis')
                # Check if the connection_analysis folder is empty
                if not os.listdir(connection_analysis_path):
                    empty_connection_analysis.append(folder_name)

    # # Print the results
    # print("Folders without 'connection_analysis/' folder:")
    # for folder in sorted(no_connection_analysis):
    #     print(folder)

    # print("\nFolders with empty 'connection_analysis/' folder:")
    # for folder in sorted(empty_connection_analysis):
    #     print(folder)

    op_check = sorted(list(set(no_connection_analysis + empty_connection_analysis)))
    return op_check

def find_connected_xlsx_files(root_directory):
    """
    Search for .xlsx files containing 'connected' in filename within all subfolders
    """
    connected_files = []
    
    # Method 1: Using glob with recursive search
    pattern = os.path.join(root_directory, '**', '*connected*.xlsx')
    connected_files = glob.glob(pattern, recursive=True)
    
    return connected_files

def find_connections_xlsx_files(root_directory):
    """
    Search for .xlsx files containing 'connected' in filename within all subfolders
    """
    connected_files = []
    
    # Method 1: Using glob with recursive search
    pattern = os.path.join(root_directory, '**', '*connections*.xlsx')
    connected_files = glob.glob(pattern, recursive=True)
    
    return connected_files

def find_connected_xlsx_files_os_walk(root_directory):
    """
    Alternative method using os.walk for more control
    """
    connected_files = []
    
    # Method 2: Using os.walk
    for root, dirs, files in os.walk(root_directory):
        for file in files:
            if file.endswith('.xlsx') and 'connected' in file.lower():
                full_path = os.path.join(root, file)
                connected_files.append(full_path)
    
    return connected_files

def create_updated_con_screen_df(human_dir, patcher, op_add_data):
    exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0])

    for op in op_add_data:
        print(f'analyzing {op}')
        work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, op, patcher, 'old')
        dir_plots = sort.make_dir_if_not_existing(work_dir, 'connection_analysis')

        con_screen_json = sort.get_json_meta_connect_all(human_dir, op, 'Verji')

        con_data = pd.DataFrame(columns = ['OP', 'fn', 'slice', 'day', 'connection_ID', 'repatch', 'treatment',
        'hrs_after_OP', 'hrs_incubation','chan_pre', 'chan_post', 'Vm pre', 'Vm post', 
        'Amp 1','Amp 2','Amp 3', 'Amp 4',	'Lat1',	'Lat2',	'Lat3',	'Lat4', 'num excluded swps', 'comments'])

        for i, indx in enumerate(con_screen_json[0]['con_screen_file']):
            con_screen_file = work_dir + filenames[indx]

            if len(con_screen_json[0]['pre_chans']) == 0:
                con_data.to_excel(work_dir + '/connection_analysis/' + op + '_connections.xlsx', index=False)
                print('no connections in this dataset')
                continue

            pre_cells = con_screen_json[0]['pre_chans'][i]
            post_cells = con_screen_json[0]['post_chans'][i]
            
            time_after_op = sort.get_time_after_OP(con_screen_file, exp_view.cortex_out[exp_view['OP'] == op].tolist()[0])
            
            exclude = ''
            for j, pre_cell in enumerate(pre_cells):
                post_cell = post_cells[j]
                slic = slice_names[indx]
                if slic in con_screen_json[0]['slices']:
                    treatment = con_screen_json[0]['treatment'][con_screen_json[0]['slices'].index(slic)]
                else:
                    treatment  = ' '
                con_ID = hcf.get_connection_ID(con_screen_file, slic, pre_cell, post_cell)
                day = 'D1'
                if slic[-2:] == 'D2':
                    day = 'D2'

                pre_signal, es, vm_pre = con_param.presynaptic_screen(con_screen_file, pre_cell)
                post_signal, vm_post = con_param.postsynaptic_screen(con_screen_file, post_cell, es)
                if (np.array(vm_post)).size == 0:
                    print('QC not passed!!')
                    exclude = 'all'
                    es2 = []
                    post_signal, vm_post = con_param.postsynaptic_screen(con_screen_file, post_cell, es2)
                    continue
                mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = \
                    con_param.get_analysis_window(pre_signal, post_signal)
                pre_signal, post_signal = con_param.remove_sweeps_with_POSTsynAPs(pre_signal, post_signal, preAPs)
                post_peaks = con_param.find_postsynaptic_peaks(post_window, preAPs_shifted)
                post_local_baseline = con_param.get_psp_baselines(post_window,preAPs_shifted)
                onsets = con_param.get_onsets(preAPs_shifted, post_window, post_peaks, post_local_baseline)
                latency = con_param.latencies(onsets, preAPs_shifted)
                amps = con_param.get_amps(post_peaks, post_local_baseline)

                df_add = pd.DataFrame({'OP': op, 'fn': filenames[indx], 'slice': slic, 'day':day, 'treatment': treatment, 
                'connection_ID' : con_ID, 'hrs_after_OP' : time_after_op,
                'chan_pre': pre_cell, 'chan_post': post_cell, 'Vm pre' :vm_pre, 'Vm post': vm_post,
                'Amp 1': amps[0], 'Amp 2': amps[1],	'Amp 3': amps[2],	'Amp 4': amps[3],	
                    'Lat1': latency[0][0],	'Lat2' : latency[1][0],	'Lat3': latency[2][0],	'Lat4': latency[3][0], 
                'num excluded swps': len(es), 'comments': exclude}, index=[0])
                con_data = pd.concat([con_data.loc[:], df_add]).reset_index(drop=True)

                #plotting
                funcs_plotting_raw_traces.plot_connection_window(con_screen_file, pre_cell, post_cell, pre_window, \
                    post_window, preAPs_shifted, post_signal, onsets, preAPs, post_peaks, post_local_baseline)
                funcs_plotting_raw_traces.plot_post_cell(con_screen_file, pre_cell, post_cell)

        con_data.to_excel(work_dir + '/connection_analysis/' + op + '_connections.xlsx', index=False) 

def add_VC_con_df(human_dir, patcher, op_add_data):
    '''
    MOST PROBABLY NEEDS FIXING
    '''

    exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0])

    for op in op_add_data:
        print(f'analyzing {op}')
        work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, op, patcher, 'old')
        dir_plots = sort.make_dir_if_not_existing(work_dir, 'connection_analysis')

        
        con_data_VC = pd.DataFrame(columns = ['OP', 'fn', 'slice', 'day', 'connection_ID', 'repatch', 'treatment', 
        'hrs_after_OP', 'hrs_incubation','chan_pre', 'chan_post', 'Vm pre', 'Holding post', 
        'Amp 1',	'Amp 2',	'Amp 3',	'Amp 4',	'Lat1',	'Lat2',	'Lat3',	'Lat4', 'num excluded swps', 'comments'])

        for i, indx in enumerate(con_screen_json[1]['con_screen_IC_file_indices']):
            con_screen_file_IC = work_dir + filenames[indx]
            
            pre_cells = con_screen_json[1]['pre_chans_IC'][i]
            post_cells = con_screen_json[1]['post_chans_VC'][i]
            if pre_cells == [] or post_cells == []:
                continue
            time_after_op = sort.get_time_after_OP(con_screen_file_IC, exp_view.cortex_out[exp_view['OP'] == op].tolist()[0])
            
            exclude = ''
            for j, pre_cell in enumerate(pre_cells):
                post_cell = post_cells[j]
                slic = slice_names[indx]
                if slic in con_screen_json[0]['slices']:
                    treatment = con_screen_json[0]['treatment'][con_screen_json[0]['slices'].index(slic)]
                else:
                    treatment  = ' '

                con_ID = hcf.get_connection_ID (con_screen_file_IC, slic, pre_cell, post_cell)

                day = 'D1'
                if slic[-2:] == 'D2': 
                    day = 'D2'

                pre_sig, es, vm_pre = con_param.presynaptic_screen_IC(con_screen_file_IC, pre_cell)
                post_sig, holding_post = con_param.postsynaptic_screen_VC (con_screen_file_IC, post_cell, es)
                if (np.array(holding_post)).size == 0:
                    print('QC not passed!!')
                    exclude = 'all'
                    es2 = []
                    post_sig, holding_post = con_param.postsynaptic_screen_VC (con_screen_file_IC, post_cell, es)
                    continue
                mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = con_param.get_analysis_window_VC(pre_sig, post_sig)
                post_peaks = con_param.find_postsynaptic_peaks_VC(post_window, preAPs_shifted)
                bl = con_param.get_psp_baselines(post_window,preAPs_shifted)
                onsets = con_param.get_onsets_VC(preAPs_shifted, post_window, post_peaks, bl)
                latency = con_param.latencies(onsets, preAPs_shifted)
                amps = con_param.get_amps_VC(post_peaks, bl)

                df_add = pd.DataFrame({'OP': op, 'fn': filenames[indx], 'slice': slic, 'day':day, 'treatment': treatment, 
                'connection_ID' : con_ID, 'hrs_after_OP' : time_after_op,
                'chan_pre': pre_cell, 'chan_post': post_cell, 'Vm pre' :vm_pre, 'Holding post': holding_post,
                'Amp 1': amps[0], 'Amp 2': amps[1],	'Amp 3': amps[2],	'Amp 4': amps[3],	
                    'Lat1': latency[0][0],	'Lat2' : latency[1][0],	'Lat3': latency[2][0],	'Lat4': latency[3][0], 
                'num excluded swps': len(es), 'comments': exclude}, index=[0])
                con_data_VC = pd.concat([con_data_VC.loc[:], df_add]).reset_index(drop=True)

                #plotting
                funcs_plotting_raw_traces.plot_connection_window_VC(con_screen_file_IC, pre_cell, post_cell, pre_window, \
                        post_window, preAPs_shifted, post_sig, onsets, preAPs, post_peaks, bl)
                funcs_plotting_raw_traces.plot_post_cell_VC(con_screen_file_IC, pre_cell, post_cell)

        con_data_VC.to_excel(work_dir + '/connection_analysis/' + op + '_connected_cell_properties_post_in_VC.xlsx', index=False) 

def get_con_df_from_analyzed(human_dir, data_type):
    # Find all .xlsx files with 'connected' in the name
    if data_type == 'connected':
        all_connected = find_connected_xlsx_files(human_dir)
    elif data_type == 'connections':
        all_connected = find_connections_xlsx_files(human_dir)
    connected_IC = [f for f in all_connected if '_in_VC' not in f]

    con_df_complete = pd.DataFrame()
    for path in connected_IC:
        df_add = pd.read_excel(path)
        if df_add.empty or df_add.isna().all().all():
            print(f'no data skipping {path[path.rfind(('OP')): path.rfind(('OP'))+8]}')
            continue

        if 'rosie' in path:
            df_add.insert(0, 'patcher', 'Rosie')
        elif 'verji' in path:
            df_add.insert(0, 'patcher', 'Verji')
        else:
            print(f'no pather found for {path[path.rfind(('OP')): path.rfind(('OP'))+8]}')
        con_df_complete = pd.concat([con_df_complete.loc[:], df_add]).reset_index(drop=True)

    con_df_complete = con_df_complete[con_df_complete.treatment != 'TTX'].reset_index(drop=True)

    # remove OPs where no data
    con_df_complete = con_df_complete[~con_df_complete.OP.isna()].reset_index(drop=True)

    return con_df_complete

def get_match_cell_id_connectivity(df):
    '''
    cell_ID = patcher_initial + OP_date + slice + channel
    '''
    if 'cell_IDs_match_pre' in df.columns:
        df = df.drop(columns = ['cell_IDs_match_pre'])

    if 'cell_IDs_match_post' in df.columns:
        df = df.drop(columns = ['cell_IDs_match_post'])

    patcher_dict = {'Verji':'vm', 'Rosie': 'rs'}

    cell_IDs_match_pre, cell_IDs_match_post = [], []
    for i in range(len(df)):

        # print(df.OP[i])
        chan_pre = df.chan_pre[i]
        chan_post = df.chan_post[i]

        slic = df.slice[i]
        
        if isinstance(chan_pre, (float, int, np.int64)):
            chan_pre = str(int(chan_pre))
            chan_post = str(int(chan_post))
        else:
            chan_pre = str(chan_pre)[-1] # ch repatched form D1
            chan_post = str(chan_post)[-1]

        cell_IDs_match_pre.append(patcher_dict[df.patcher[i]] +
                              df.OP[i][2:] +
                              slic +
                              'c' + chan_pre)
        cell_IDs_match_post.append(patcher_dict[df.patcher[i]] +
                              df.OP[i][2:] +
                              slic +
                              'c' + chan_post)

    df.insert(len(df.columns), 'cell_IDs_match_pre', cell_IDs_match_pre)
    df.insert(len(df.columns), 'cell_IDs_match_post', cell_IDs_match_post)

    return df

def keep_cells_with_intr_data_only(con_df, intr_df):
    '''
    only connections with both pre and post cell in intr_df are kept
    '''
    # remove dead cells
    indx_keep = []
    for i, pre_cell in enumerate(con_df.cell_IDs_match_pre):
        post_cell = con_df.cell_IDs_match_post[i]

        if pre_cell in intr_df.cell_IDs_match.unique() and \
            post_cell in intr_df.cell_IDs_match.unique():
            indx_keep.append(i)

    return con_df.iloc[indx_keep].reset_index(drop=True)

def estimate_connectivity_percentage(con_df, intr_df):
    '''
    returns a df with connectivity estimation for each slice, only QC-ed cells included
    '''
    ops, slices, days, trs, hrs, num_cells, possible_cons, found_cons = [], [], [], [], [], [], [], []
    count_one_cell_per_slice = 0
    for _, op in enumerate(intr_df.OP.unique()):
        intr_op_df = intr_df[intr_df.OP == op]

        for slic in intr_op_df.slice.unique():
            df_slice = intr_op_df[intr_op_df.slice == slic]

            if len(df_slice) == 1:
                count_one_cell_per_slice += 1
                continue

            ops.append(op)
            slices.append(slic)
            days.append(df_slice[df_slice.slice == slic].day.unique()[0])
            trs.append(df_slice[df_slice.slice == slic].treatment.unique()[0])
            hrs.append(df_slice[df_slice.slice == slic].hrs_after_OP.unique()[0])
            num_cells.append(len(df_slice))
            possible_cons.append(len(df_slice) * (len(df_slice) - 1))
        
            if op in con_df.OP.unique():
                con_df_op = con_df[con_df.OP == op]
                
                if slic in con_df_op.slice.unique():
                    found_cons.append(len(con_df_op[con_df_op.slice == slic]))
                else:
                    found_cons.append(0)
            else:
                found_cons.append(0)

    con_percentage = np.array(possible_cons) / np.array(found_cons)

    #calcualte connectivity percentage, not dividing by 0
    con_percent_fixed = []
    for con_percent in con_percentage:
        if math.isinf(con_percent):
            con_percent_fixed.append(0)
        else:
            con_percent_fixed.append(con_percent)

    con_percentage_df = pd.DataFrame({'OP': ops, 'slices': slices, 'day': days,
                                    'treatment': trs, 'hrs_after_OP': hrs,
                                    'num_channels': num_cells,
                                    'num_possible_connections': possible_cons,
                                    'found_connections': found_cons,
                                    'con_percentage': con_percentage,
                                    'con_percentage_fixed': con_percent_fixed})
    return con_percentage_df, count_one_cell_per_slice

def get_connections_df_from_jsons(human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    OP_dirs = glob.glob(human_dir + 'data_*/' + 'OP*')
    #OP_skip = ['OP210304', 'OP210319']

    OPs, slices, days, treatments, num_chans, possible_cons, num_connections  = [], [], [], [], [], [], []

    # collect all connections enterted in the excel labbook during recording
    for i, OP_dir in enumerate(OP_dirs):
        path_json = glob.glob(OP_dir + '/*_con_screen_only.json')
        OP =  OP_dirs[i][OP_dirs[i].find('OP'):]

        if len(path_json) == 0:
            print(f'skipping {OP}')
            continue

        con_dict = sort.open_json(path_json[0])[0]
        for j in range(len(con_dict['active_chans'])):
            OPs.append(OP)
            slices.append(con_dict['slices'][j])
            if len(con_dict['slices'][j]) < 3:
                day = 'D1'
            else:
                day = 'D2'
            days.append(day)
            treatments.append(con_dict['treatment'][j])
            num_chans.append(len(con_dict['active_chans'][j]))
            possible_cons.append(len(con_dict['active_chans'][j])*(len(con_dict['active_chans'][j])-1))
            if 'pre_chans' not in con_dict.keys() or len(con_dict['pre_chans']) == 0:
                num_connections.append(0)
                print(f'skipping {OP}')
                # print('For ',OP , 'slice ' , con_dict['slices'][j] , 'on ' + day + 'treatment ' + con_dict['treatment'][j], 'num chans' + str(len(con_dict['active_chans'][j])) + 'num possible connections' , str(len(con_dict['active_chans'][j])*(len(con_dict['active_chans'][j])-1)), con_dict['treatment'][j] ,  'no connectios')
            else:
                num_connections.append(len(con_dict['pre_chans'][j]))
                # print('For' , OP, 'slice ' , con_dict['slices'][j] , 'on ' , day , 'treatment ' ,con_dict['treatment'][j], 'num chans' , str(len(con_dict['active_chans'][j])) , 'num possible connections', str(len(con_dict['active_chans'][j])*(len(con_dict['active_chans'][j])-1)) ,con_dict['treatment'][j] , str(len(con_dict['pre_chans'][j])) ,'connectios found')

    df_connections = pd.DataFrame({ 'OP':OPs, 'slices':slices, 'day':days, 'treatment':treatments,
        'num_channels':num_chans, 'num_possible_connections':possible_cons, 'found_connections':num_connections})
    return df_connections

def add_connectivity_percentage(df_connections):
    df_connections.insert(len(df_connections.columns), 'con_percentage',
                          df_connections.num_possible_connections / df_connections.found_connections)

    #calcualte connectivity percentage, not dividing by 0
    con_percent_fixed = []
    for con_percent in df_connections['con_percentage']:
        if math.isinf(con_percent):
            con_percent_fixed.append(0)
        else:
            con_percent_fixed.append(con_percent)

    df_connections.insert(len(df_connections.columns), 'con_percentage_fixed', con_percent_fixed)

    return df_connections

def cross_check_intr_con_df(intr_df, df_connections):
    # remove ops where no intrinsic
    ops_exclude = sorted(list(set(df_connections.OP) - set(intr_df.OP.unique())))
    df_connections = df_connections[~df_connections['OP'].isin(ops_exclude)].reset_index(drop = True)
    print(f'after removing ops not in intrinsic {len(df_connections)}')

    df_connections = get_match_cell_id_connectivity(df_connections)
    df_connections_QC = keep_cells_with_intr_data_only(df_connections, intr_df)
    print(f'after removing cells not in intrinsic {len(df_connections_QC)}')

    duplicates_connections = df_connections_QC[df_connections_QC.duplicated(subset=['OP', 'slice', 'day', 'chan_pre', 'chan_post'], keep='first')]
    if len(duplicates_connections) == 0:
        print('happy to continue, no duplicates in connections')

    ops_0_connectivity = sorted(set(intr_df.OP.unique()) - set(df_connections.OP.unique()))
    print(f'ops with no connectivity data {ops_0_connectivity}')

    combination_counts = intr_df.groupby(['OP', 'slice']).size().reset_index(name='count')
    print(f"Combinations OP x slice {len(combination_counts)}")

    df_con_percentage, num_alone_cell_in_slice = estimate_connectivity_percentage(df_connections_QC, intr_df)

    if num_alone_cell_in_slice + len(df_con_percentage) != len(combination_counts):
        print('something somewhere went wrong')
    else:
        print('happy baby')
        return df_con_percentage, num_alone_cell_in_slice

def quick_dirty_plot_connectivity(df_connections, title_, save_dir = False):

    if 'con_percentage_fixed' not in df_connections.columns:
        df_connections = add_connectivity_percentage(df_connections)
    fig, ax = plt.subplots(1, 1,figsize=(7,7))
    colors = ['pink', 'green', 'orange']

    mean_con = []

    for j, day in enumerate(['D1', 'D2']):
        for i, tr in enumerate(['Ctrl', 'high K']):
            # if day == 'D1':
            #     k_ = 1
            # else:
            #     k_ = 2
            k_ = j + 2*i
            df_plot = df_connections.loc[(df_connections.treatment == tr)\
                                         & (df_connections.day == day)]
            avg = np.mean(df_plot.con_percentage_fixed)
            #print(median)
            x = np.linspace(k_ - 0.3, k_ + 0.3, len(df_plot))
            ax.scatter(x, df_plot.con_percentage_fixed, c = colors[i])
            if k_ == 1:
                ax.scatter(k_, avg, marker = '_', s = 2000, c = colors[i])
            else:
                ax.scatter(k_, avg, marker = '_', s = 2000, c = colors[i], label = tr)
            ax.text(k_, 1.2*avg, str(round(avg,0)), size = 10)
            ax.set_xticks(ticks = [1,2], labels = ['D1', 'D2'])

        mean_all = np.mean(df_connections[df_connections.day == day].con_percentage_fixed)
        mean_con.append(mean_all)
        ax.scatter(k_, mean_all, marker = '_', s = 2000, c = 'k')
    fig.suptitle(title_)
    if save_dir:
        plt.savefig(f'{save_dir}{title_}connectivity_percentage.svg',
            format='svg', bbox_inches='tight', dpi=300)
    plt.show()

#%%

# what's happening
# I've analyzed the connections twice
# most up to date data is in /OP**/connection_analysis
# i've checked and all needed data is analyzed

human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'


# add data where needed, OPs found in both
# op_check = check_which_to_analyse()
# op_add_data = sorted(list(set(adult_df_temp_complete.OP.tolist()) & set(op_check)))
# patcher = 'Verji'
# create_updated_con_screen_df(human_dir, patcher, op_add_data)

# the connections ones, should be less and with a thorough quality check
# df_connected = get_con_df_from_analyzed(human_dir, 'connected')
df_connections = get_con_df_from_analyzed(human_dir, 'connections')
print(f'all collected connectivity dfs {len(df_connections)}')

# Remove all duplicate rows, keeping the first occurrence
df_connections = df_connections.drop_duplicates().reset_index(drop = True)
print(f'after removing duplicates {len(df_connections)}')

# all QCed cells
adult_df_temp_complete = pd.read_excel(data_dir + 'adult_16min_hrs_inc_temporal.xlsx')
df_con_percentage_all, num_alone_cell_all = cross_check_intr_con_df(adult_df_temp_complete, df_connections)
quick_dirty_plot_connectivity(df_con_percentage_all, 'all possible QC-ed')


# Very quck cross check, with all cells and connections. Almost no QC
# percentage connectivity estimation - no QC, just for comparison
# df_con_jsons = get_connections_df_from_jsons()
# df_con_jsons = add_connectivity_percentage(df_con_jsons)
# df_con_jsons = df_con_jsons.drop_duplicates(keep='first').reset_index(drop=True)
# ops_remove = sorted(set(df_con_jsons.OP.unique().tolist()) - set(adult_df_temp_complete.OP.unique()))
# df_con_jsons = df_con_jsons[~df_con_jsons['OP'].isin(ops_remove)].reset_index(drop=True)
# # plot
# quick_dirty_plot_connectivity(df_con_jsons, 'from json')


# single connections, amp and latency
df_connections = get_match_cell_id_connectivity(df_connections)
df_connections_QC = keep_cells_with_intr_data_only(df_connections, adult_df_temp_complete)
df_connections_QC = df_connections_QC[df_connections_QC.treatment != 'fixed']
df_connections_QC = df_connections_QC[(df_connections_QC['Vm pre'] < -45) & \
                                      (df_connections_QC['Vm post'] < -45)]

# save everything
# save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
# df_connections_QC.to_excel(save_dir + 'all_QCed_connections.xlsx')
# df_connections_QC.to_csv(save_dir + 'all_QCed_connections.csv')
# df_con_percentage_all.to_excel(save_dir + 'connectivity_percent_estimation.xlsx')
# df_con_percentage_all.to_csv(save_dir + 'connectivity_percent_estimation.csv')


# repatch
df_repatch = pd.read_excel(data_dir + 'repatch_data_temporal.xlsx')
df_con_percentage_repatch, num_alone_cell_repatch = cross_check_intr_con_df(df_repatch, df_connections)
quick_dirty_plot_connectivity(df_con_percentage_repatch, 'repatch')

df_connections_QC_repatch = keep_cells_with_intr_data_only(df_connections_QC, df_repatch)
df_connections_QC_repatch = df_connections_QC_repatch[df_connections_QC_repatch.treatment != 'fixed']

# op_color_dict = pl_intr.get_op_color_dict(df_connections_QC_repatch)
pl_intr.plot_connect_amplitude(df_connections_QC_repatch, 'repatch', op_color_dict, 'connection_ID', results_ = 'amp')
pl_intr.plot_connect_amplitude(df_connections_QC_repatch, 'repatch_lat', op_color_dict, 'connection_ID', results_ = 'lat')




#%%

# COLLECT DISAPPEARING CONNECTIONS

human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'

df_con_screen_all,  ops_no_con_screen = collect_con_screen_from_jsons()
adult_df_all = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human'+\
                             '/paper_figs_collected_checked/data/adult_complete_0_hrs_incubation.xlsx')
ops_worth_check = []
for op in ops_no_con_screen:
    if op in adult_df_all.OP.unique():
        ops_worth_check.append(op)

print(f'ops to check {ops_worth_check}')

# find all repatched cells, no QC
df_repatch_no_QC = load_intr_adult_repatch_no_QC()

# remove TTX
df_repatch_no_QC = df_repatch_no_QC[df_repatch_no_QC.treatment != 'TTX'].reset_index(drop = True)
df_con_screen = df_con_screen_all[df_con_screen_all.treatment != 'TTX'].reset_index(drop = True)

# remove OPs that are not in the collected repatched intrinsic data
ops_drop = []
for op in df_con_screen.OP.unique():
    if op not in df_repatch_no_QC.OP.unique():
        ops_drop.append(op)
# adding 'OP230810' because only 4 hrs incubation
ops_drop.append('OP230810')

df_repatch_no_QC = df_repatch_no_QC[~df_repatch_no_QC['OP'].isin(ops_drop)].reset_index(drop = True)
df_con_screen = df_con_screen[~df_con_screen['OP'].isin(ops_drop)].reset_index(drop = True)

op_slice_not_in_repatch = for_manual_check_and_edit(df_con_screen, df_repatch_no_QC)
print(f'missing op-slice combinations (in repatch intr; not in con_screen)\n{op_slice_not_in_repatch}')

# creating unique IDs and checking how unique they are
df_con_screen_json = create_unique_cell_IDs(df_con_screen, 'con_json')
df_repatch_no_QC_intr = create_unique_cell_IDs(df_repatch_no_QC, 'intr')

# collect repatched connections; maybe some excluded already
df_connections = get_con_df_from_analyzed(human_dir, 'connections')
df_con_repatch_collected_no_QC = df_connections[df_connections.repatch == 'yes'].reset_index(drop = True)
df_con_repatch_collected_no_QC = create_unique_cell_IDs(df_con_repatch_collected_no_QC, 'con')

# removing cells from json_collected that are not in the repatched intr
indx_keep = []
for i, pre_cell in enumerate(df_con_screen_json.pre_cell_ID):
    post_cell = df_con_screen_json.post_cell_ID[i]
    if pre_cell in df_repatch_no_QC_intr['cell_ID_match'].unique() and \
         post_cell in df_repatch_no_QC_intr['cell_ID_match'].unique():
        indx_keep.append(i)
df_con_screen_json_repatched_pre_post = df_con_screen_json[df_con_screen_json.index.isin(indx_keep)]

# looked at all 8 diff between repatched pre/post and collected cons repatch, no disappearing
sorted(set(df_con_screen_json_repatched_pre_post.con_ID) - set(df_con_repatch_collected_no_QC.con_ID))

# checking manually if all of those are indeed repatched connections
# creating a list of con_IDs to remove from the df_con_repatch_collected_no_QC
# after removing them --> the real df with repatched connections
sorted(set(df_con_repatch_collected_no_QC.con_ID) - set(df_con_screen_json_repatched_pre_post.con_ID))

# exclude connections that are in con_IDs_remove
con_IDs_remove = get_con_IDs_remove_manually_checked()
df_con_repatch_collected_QC = df_con_repatch_collected_no_QC[~df_con_repatch_collected_no_QC.con_ID.isin(con_IDs_remove)].reset_index(drop = True)

if len(df_con_repatch_collected_QC[df_con_repatch_collected_QC.duplicated(keep=False)]) == 0:
    print(f'all repatched cons are duplicated')

df_con_repatch_collected_QC[(df_con_repatch_collected_QC.hrs_incubation > 0) & (df_con_repatch_collected_QC.hrs_incubation < 16)]

# removing surgery with short incubation time
df_con_repatch_collected_QC = df_con_repatch_collected_QC[df_con_repatch_collected_QC.OP != 'OP230810']

# save and go over borderline cells again
df_con_repatch_collected_QC[(df_con_repatch_collected_QC['Vm pre'] > -45) |
                            (df_con_repatch_collected_QC['Vm post'] > -45) |
                            (df_con_repatch_collected_QC['Lat1'] < 0)].sort_values(['OP', 'slice'])# .to_excel('/Users/verjim/Desktop/temp_check_cells.xlsx')

connection_IDs_remove = get_connection_IDs_remove_manually_checked_2()
df_con_repatch_collected_QC = df_con_repatch_collected_QC[~df_con_repatch_collected_QC.connection_ID.isin(connection_IDs_remove)].reset_index(drop = True)

df_con_repatch_collected_QC[df_con_repatch_collected_QC['num excluded swps'] > 5]

# now the results can be plotted
plot_hist_params = ['hrs_incubation', 'Vm pre', 'Vm post', 'Amp 1', 'Lat1', 'num excluded swps']
for param in plot_hist_params:
    plt.hist(df_con_repatch_collected_QC[param].values)
    plt.title(param)
    plt.show()


# df_con_repatch_collected_QC.to_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/connectivity/repatched_connections_looser_QC.xlsx')

#%%

# adding connectivity plots if missing

op = 'OP240503'
con_screen_indx = 33
pre_cells = [8]
post_cells = [2]
patcher = 'Verji'
vc_file_indx = 0

for i, pre_cell in enumerate(pre_cells):
    post_cell = post_cells[i]
    work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, op, patcher, 'old')

    # for getting missing plot
    con_screen_file = work_dir + filenames[con_screen_indx]

    pre_signal = con_param.presynaptic_screen_ALLOW_ALL(con_screen_file, pre_cell)
    post_signal, vm_post = con_param.postsynaptic_screen_ALLOW_ALL(con_screen_file, post_cell, [])

    mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = \
        con_param.get_analysis_window(pre_signal, post_signal)
    pre_signal, post_signal = con_param.remove_sweeps_with_POSTsynAPs(pre_signal, post_signal, preAPs)
    post_peaks = con_param.find_postsynaptic_peaks(post_window, preAPs_shifted)
    post_local_baseline = con_param.get_psp_baselines(post_window,preAPs_shifted)
    onsets = con_param.get_onsets(preAPs_shifted, post_window, post_peaks, post_local_baseline)
    latency = con_param.latencies(onsets, preAPs_shifted)
    amps = con_param.get_amps(post_peaks, post_local_baseline)

    # plot from con_screen
    plotting_funcs.plot_connection_window(con_screen_file, pre_cell, post_cell, pre_window, \
        post_window, preAPs_shifted, post_signal, onsets, preAPs, post_peaks, post_local_baseline)




# FROM VC FILE
con_screen_file_IC = work_dir + filenames[vc_file_indx]

pre_sig, es, vm_pre = con_param.presynaptic_screen_IC(con_screen_file_IC, pre_cell)
post_sig, holding_post = con_param.postsynaptic_screen_VC (con_screen_file_IC, post_cell, es)
if (np.array(holding_post)).size == 0:
    print('QC not passed!!')
    exclude = 'all'
    es2 = []
    post_sig, holding_post = con_param.postsynaptic_screen_VC (con_screen_file_IC, post_cell, es)
    
mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs = con_param.get_analysis_window_VC(pre_sig, post_sig)
post_peaks = con_param.find_postsynaptic_peaks_VC(post_window, preAPs_shifted)
bl = con_param.get_psp_baselines(post_window,preAPs_shifted)
onsets = con_param.get_onsets_VC(preAPs_shifted, post_window, post_peaks, bl)
latency = con_param.latencies(onsets, preAPs_shifted)
amps = con_param.get_amps_VC(post_peaks, bl)

plotting_funcs.plot_connection_window_VC(con_screen_file_IC, pre_cell, post_cell, pre_window, \
        post_window, preAPs_shifted, post_sig, onsets, preAPs, post_peaks, bl)
plotting_funcs.plot_post_cell_VC(con_screen_file_IC, pre_cell, post_cell)

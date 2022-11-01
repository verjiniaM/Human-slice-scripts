


OP = 'OP220217'
patcher = 'Verji'

human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'

work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)

old_con_df = pd.read_excel(work_dir + '/data_tables/' + OP + '_connected_cell_properties.xlsx')

amps_1, amps_2, amps_3, amps_4, con_IDs = [], [], [], [], []
for j, fn in enumerate(old_con_df['fn']):
    con_screen_file = work_dir + fn

    pre_cell = old_con_df['chan_pre'][j]
    post_cell = old_con_df['chan_post'][j]
    slic = slice_names[filenames.index(fn)]

    con_ID = hcf.get_connection_ID (con_screen_file, slic, pre_cell, post_cell)

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
    amps = con_param.get_amps(post_peaks, post_local_baseline)

    #if old_con_df['Amp1'][j] == amps[0]:
    amps_1.append(amps[0])
    amps_2.append(amps[1])
    amps_3.append(amps[2])
    amps_4.append(amps[3])
    con_IDs.append(con_ID)

old_con_df.insert(14, 'Amp 1', amps_1)
old_con_df.insert(15, 'Amp 2', amps_2)
old_con_df.insert(16, 'Amp 3', amps_3)
old_con_df.insert(17, 'Amp 4', amps_4)
old_con_df.insert(5, 'connection_ID', con_IDs)

old_con_df = old_con_df.drop(['Unnamed: 0','Amp1', 'Amp2', 'Amp3', 'Amp4', 'cell_ID'], axis=1)

old_con_df.to_excel(work_dir + '/data_tables/' + OP + '_connected_cell_properties.xlsx') 

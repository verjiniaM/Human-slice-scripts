
import numpy as np
from scipy.signal import find_peaks
import stimulation_windows_ms as stim_win


def get_analysis_window_diff_freqs (pre_cell_chan, post_cell_chan, hz):
    con_screen_data = load_traces(con_screen_file)
    pre_ch, post_ch = 'Ch' + str(pre_cell_chan), 'Ch' + str(post_cell_chan)

    sweep_len = np.shape(con_screen_data[pre_ch][0])[0]

    stims = stim_win.stim_window_diff_freq
    win_start, win_end = stims[hz][0], stims[hz][1] 

    mean_pre = np.mean(con_screen_data[pre_ch][0], axis = 1)[win_start:win_end]
    vm0 = np.mean(con_screen_data[pre_ch][0][:,0][0:950]) #from column 0, 0:950
    preAPs = find_peaks(mean_pre, height=0)
    j = 0
    while len(preAPs[0]) == 0:
        j = j-5
        preAPs=find_peaks(mean_pre, height=j)
    print("Pre APs found with heigth " + str(j))
    num_aps = len(preAPs[0])

    pre_window = mean_pre[preAPs[0][0]-30:preAPs[0][num_aps-1]+300]

    mean_post = np.mean(con_screen_data[post_ch][0], axis = 1)
    post_window = mean_post[preAPs[0][0]-30:preAPs[0][num_aps-1]+3000] 
    preAPs_shifted = find_peaks(pre_window, height=j)
    pA0 = np.mean(con_screen_data[post_ch][0][:,0][0:950]) #from column 0, 0:950

    return vm0, pA0, mean_pre, mean_post, pre_window, post_window, preAPs_shifted, preAPs


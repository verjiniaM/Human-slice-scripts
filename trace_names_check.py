
import neo
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import human_characterisation_functions as hchf


#%%
# =============================================================================
# some fast plotting functions
# =============================================================================

# plots the middle sweep for each channel
# check the traces to see which channels were active or if the protocol names are enetered correctly
def plot_traces(filename): 
    r = neo.io.AxonIO(filename=filename)
    block=r.read(signal_group_mode = 'split-all')[0]
    sweeps = len(block.segments)
    sweep_len = len(block.segments[0].analogsignals[0])
    channels = len(block.segments[0].analogsignals) # no of channels in dep rec
    channel_dict = {}
    for ch in range(0, channels):
        name = block.segments[0].analogsignals[ch].name
        if name == '_Ipatch':
            name = 'Ch1'
            block.segments[0].analogsignals[ch].name = 'Ch1'
        channel_dict['AnalogSig%d' % (ch)] = name

    dir_plots = "plots"
    end_filename = filename.rfind('/')+1
    path1 = os.path.join(filename[:end_filename], dir_plots)
    if os.path.isdir(path1) == False: os.mkdir(path1)
    path = os.path.join(filename[:end_filename],dir_plots + '/traces')
    if os.path.isdir(path) == False: os.mkdir(path)

    plt.style.use(['fast'])
    x = math.ceil(np.sqrt(channels))
    fig, ax = plt.subplots(x,x,sharex=True, sharey=True,figsize=(6,6))

    for i in range(1,channels+1):
        ax = fig.add_subplot(x,x, i)
        cell_chan = i
        #ch1 = np.ndarray([1, sweeps])  #empty matrix
        # for key, val in channel_dict.items():
        #     if val == 'Ch'+ str(cell_chan): #-1
        #         cell_chan = int(key[-1])                  
        # for i in range(0,lens(block.segments)):
        middle_swp_num = int(sweeps/2)
        ch1 = block.segments[middle_swp_num].analogsignals[cell_chan-1].view(np.recarray).reshape(sweep_len).tolist()
        # plot_swp_num = int(ch1.shape[1]/2+1)
        ch_name = block.segments[0].analogsignals[cell_chan-1].name
        sampl_rate = block.segments[middle_swp_num].analogsignals[cell_chan-1].sampling_rate
        units = block.segments[middle_swp_num].analogsignals[cell_chan-1].units
        times = np.linspace(0,sweep_len,sweep_len)/sampl_rate
        ax.plot(times,ch1, lw=0.5)
        ax.set_xlabel('sec')
        ax.set_ylabel(str(units)[-2:])
        ax.title.set_text(ch_name)
        del ch1

    fig.tight_layout() 
    fig.patch.set_facecolor('white')
    plt.savefig(path + '/trace_plot_' + filename[end_filename:-4] + '.png')
    plt.close(fig)
    return 'Trace plots saved in'+path

#plots the sweeps with the min value and the 
def plot_mini_sweeps (filename, cell_chan, sweep):
    sweep = sweep - 1
    ch1, sweep_len, block = hchf.load_traces(filename, cell_chan) #ch1: each sweep is a column
    #no test pulse

    signal_no_test_pulse = ch1[4250:,sweep]
    
    min_val = np.amin(signal_no_test_pulse)
    max_val = np.amax(signal_no_test_pulse)
    #loc_max = np.where(signal_no_test_pulse == max_val)
    loc_min = np.where(signal_no_test_pulse == min_val)
    loc_max = np.where(signal_no_test_pulse == max_val)
    num_points = np.shape(signal_no_test_pulse)[0]
    times = np.linspace(0,num_points,num_points)/20000
    plt.plot(times,signal_no_test_pulse, lw=0.5)
    plt.plot(int(loc_min[0])/20000, np.amin(signal_no_test_pulse),'ro')
    plt.plot(int(loc_max[0])/20000, np.amax(signal_no_test_pulse),'k*')
    plt.xlabel('sec')
    plt.ylabel('pA')
    plt.show()

    #plot only min_interval
    min_interval = signal_no_test_pulse[int(loc_min[0][0])-2000:int(loc_min[0][0])+2000]
    min_val_int= np.amin(min_interval)
    loc_min_int = np.where(min_interval == min_val_int)

    num_points_int = np.shape(min_interval)[0]
    times_int = np.linspace(0,num_points_int,num_points_int)/20000
    plt.plot(times_int,min_interval)
    plt.plot(int(loc_min_int[0][0])/20000, min_val_int,'ro')
    plt.xlabel('sec')
    plt.ylabel('pA')
    plt.show()

    max_interval = signal_no_test_pulse[int(loc_max[0][0])-2000:int(loc_max[0][0])+2000]
    max_val_int= np.amin(max_interval)
    loc_max_int = np.where(max_interval == max_val_int)

    num_points_int = np.shape(max_interval)[0]
    plt.plot(times_int,max_interval)
    plt.plot(int(loc_max_int[0][0])/20000, max_val_int,'ro')
    plt.xlabel('sec')
    plt.ylabel('pA')
    plt.show()


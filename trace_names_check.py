
import neo
import numpy as np
import matplotlib.pyplot as plt
import math
import os


#%%
# =============================================================================
# Basic functions for characterisation of human neurons
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
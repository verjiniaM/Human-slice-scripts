#%%
import sorting_functions as sort
import plotting_funcs
import human_characterisation_functions as hcf
import pandas as pd

human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
OP = 'OP221116' # and also OP220615 and OP220413
patcher = 'Verji'

work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)
file_list = sort.get_sorted_file_list(work_dir)
jsons = sort.get_json_files(file_list)
if OP + '_indices_dict.json' in jsons:
    indices_dict = sort.from_json(work_dir, OP, '_indices_dict.json')
else: 
    sort.to_json(work_dir, OP, '_indices_dict.json', indices_dict)

#creating a dir to save plots and data_tables (if not existing)
dir_plots = sort.make_dir_if_not_existing (work_dir, 'plots')
sort.make_dir_if_not_existing (work_dir, 'data_tables')

#check if the traces dir is empty and only then plot the mimiddle sweep for each filename
traces_folder =  os.path.join(dir_plots, "traces/")
if os.path.isdir(traces_folder) == 0 :
    for rec in range(len(filenames)):
        filename = work_dir + filenames[rec]
        plotting_funcs.plot_middle_sweep(filename)
else:
        print("skipping plotting")

active_chans_meta = sort.get_json_meta_high_K(human_dir, OP, patcher, '_meta_active_chans_high_K.json')
cortex_out_time = sort.get_datetime_from_input(active_chans_meta[0]['OP_time'][0])

df_vc = pd.DataFrame(columns=['filename', 'slice', 'cell_ch', 'hrs_after_OP', 'recording_time',
'circumstance', 'Rs', 'Rin']) #circumstance: no high K, wash in, wash out

for i in range(len(indices_dict['vc'])):
    vc = indices_dict['vc'][i]
    slic = slice_names[vc]

    filename_vc = work_dir + filenames[vc]
    time_after_op = sort.get_time_after_OP(filename_vc, cortex_out_time)

    active_channels = active_chans_meta[0]['active_chans'][i]
    
    time_after_op = sort.get_time_after_OP(filename_vc, cortex_out_time)
    rec_time = str(hcf.get_recording_time(filename_vc))
    Rs, Rin = hcf.get_access_resistance(filename_vc, active_channels)

    add_df = pd.DataFrame({'filename': filenames[vc], 'slice' : slic, 'cell_ch': active_channels,
    'hrs_after_OP' : time_after_op, 'recording_time': rec_time, 'Rs' : Rs, 'Rin': Rin })

    plotting_funcs.plot_vc_holding (filename_vc, active_channels)

    df_vc = pd.concat([df_vc.loc[:], add_df]).reset_index(drop=True)
df_vc.to_excel(work_dir + 'data_tables/' + OP + '_changes_in_resistance.xlsx', index=False) 

df_AP_props = pd.DataFrame(columns=['filename', 'slice', 'cell_ch', 'hrs_after_OP', 'circumstance',
'recording_time', 'max_spikes', 'Rheobase', 'AP_heigth', 'TH', 
'max_depol', 'max_repol', 'membra_time_constant_tau', 'capacitance'])

for i in range(len(indices_dict['freq analyse'])):
    char = indices_dict['freq analyse'][i]
    slic = slice_names[char]

    filename_cahr = work_dir + filenames[char]
    time_after_op = sort.get_time_after_OP(filename_char, cortex_out_time)

    active_channels = active_chans_meta[1]['active_chans'][i]
    
    time_after_op = sort.get_time_after_OP(filename_char, cortex_out_time)
    rec_time = str(hcf.get_recording_time(filename_char))

    meta_df = pd.DataFrame({'filename': filenames[char], 'slice' : slic, 'cell_ch': active_channels,
    'hrs_after_OP' : time_after_op, 'recording_time': rec_time })

    charact_params  = hcf.all_chracterization_params(filename_char, active_channels, 'full')
    df_char = pd.DataFrame.from_dict(charact_params)

    df_to_add = pd.concat([meta_df, df_char], axis = 1)
    df_AP_props = pd.concat([df_AP_props.loc[:], df_to_add]).reset_index(drop=True)

    plotting_funcs.plots_for_charact_file(filename_char, active_channels, inj)

df_AP_props.to_excel(work_dir + 'data_tables/' + OP + '_changes_in_AP_props.xlsx', index=False) 



df_resting = pd.DataFrame(columns=['filename', 'slice', 'cell_ch', 'hrs_after_OP', 'circumstance',
'recording_time', 'resting_potential'])












def load_traces(filename,cell_chan): 
    r = neo.io.AxonIO(filename=filename)
    block=r.read(signal_group_mode = 'split-all')[0]
    rec_time = block.rec_datetime
    sweeps = len(block.segments)
    sweep_len = len(block.segments[0].analogsignals[0])
    channels = len(block.segments[0].analogsignals) # no of channels in dep rec
    channel_dict = {}
    for ch in range(0, channels):
        name = block.segments[0].analogsignals[ch].name
        if name == '_Ipatch' or name == 'IN0':
            name = 'Ch1'
            block.segments[0].analogsignals[ch].name = 'Ch1'
        channel_dict['AnalogSig%d' % (ch)] = name

    ch1 = np.ndarray([sweep_len, sweeps])  #empty matrix              
    for key, val in channel_dict.items():
        if val == 'Ch'+ str(cell_chan): #-1
            cell = int(key[-1])    
    for i in range(0,len(block.segments)):
        ch1[:,i] = block.segments[i].analogsignals[cell].view(np.recarray).reshape(sweep_len)
    # for i in range(0, ch1.shape[1]):
    #     plt.plot(ch1[:,i], lw=0.5)
    #     plt.title('Ch ' + str(cell_chan))
    #     plt.savefig(path + '\\' + filename[-12:-4]+'_'+str(cell_chan) + '_trace_plot.png')
    #     plt.close()

    return ch1, sweep_len, block

#from characterization file
def vm(ch1):
    Vm = np.median(ch1[:,5])
    return Vm

def APprops(filename, cell_chan, inj_len):
    ch1, sweep_len, block = load_traces(filename, cell_chan) #in ch1 each sweep is a column
    inj=[-300,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400,450,500,550,
         600,700,800,900,1000,1100,1200,1300,1400]
    inj = inj[0:inj_len]
    #start and end of the pulses
    onset=2624
    offset=22624
    mc = np.ndarray([5,3])

    # creating folder to save plots
    dir_plots = "Onset"
    end_filename = filename.rfind('/')+1
    path = os.path.join(filename[:end_filename]+'plots/', dir_plots)
    if os.path.isdir(path) == False: os.mkdir(path)

    # plotting
    clrs = ["b", "g", "r", "c", "m", "y", "#FF4500", "#800080"]
    fig = plt.figure()
    for i in range(0,5):
        I = inj[i]*1e-12 #check step size
        bl = np.median(ch1[0:onset-20,i])
        ss = np.median(ch1[offset-2000:offset-1,i]) #steady state, during the current step
        swp = list(ch1[:,i])
        Vdiff = bl-ss #voltage deflection size
        v65 = Vdiff*0.63
        V65 = bl-v65
        res = list(filter(lambda ii: ii < V65, swp))[0] #takes the first value in swp < V65
        tau65 = swp.index(res) #index of res
        R = (Vdiff/1000)/-I     
        tc = tau65 - onset
        mc[i,0] = tc*0.05 #membranec capacitance; tc - time constant
        mc[i,1] = R*1e-6 #resistance
        mc[i,2] = tc*5e-5/R #capacitance
        plt.plot(ch1[:,i], c=clrs[i])
        plt.scatter(onset+tc, V65, c=clrs[i])
        plt.annotate('V65  ', (onset+tc, V65), horizontalalignment='right')
        plt.scatter(onset,bl,c='r')

    fig.patch.set_facecolor('white')    
    plt.annotate('  Baseline', (onset,bl))
    plt.title('Ch ' + str(cell_chan))
    plt.savefig(path + '/Char_onset_plot_' + filename[end_filename:-4]+'_'+str(cell_chan) + '.png')
    plt.close()
    mc[:,2] = mc[:,2]/1e-12  
    capacitance=mc[1,2]
    tau=mc[1,0]    

    Vm = np.median(ch1[:,5]) #when the inj = 0mV
    max_spikes = 0
    for i, j in enumerate(range(0, len(ch1[0]))): #loop through all swps
        pks = detect_peaks(ch1[:,j], mph=0,mpd=50) # detects the peaks for each timepoint? 
        if len(pks)> max_spikes:
            max_spikes = len(pks) #find the max number of spikes (peaks)
    if max_spikes == 0: 
        print("No spikes found for" + filename[end_filename:-4] + ' Ch: ' + str(cell_chan))
        TH = math.nan
        max_depol = math.nan
        max_repol = math.nan
        APheight = math.nan
        Rheobase = math.nan
        return Vm, max_spikes, Rheobase,APheight, max_depol, max_repol, TH, capacitance, tau
    else:
        spike_counts = np.ndarray([len(inj),2])
        peaks = np.empty([len(inj), max_spikes, 3])
        peaks.fill(np.nan)
        for i, j in enumerate(range(0, len(ch1[0]))): #forr all swps
            pks = detect_peaks(ch1[:,j], mph=0,mpd=50) 
            peaks[i,0:len(pks), 0] = inj[i] #injected current
            peaks[i,0:len(pks), 1] = pks #
            peaks[i,0:len(pks), 2] = ch1[pks,j] #sweep number
        # number of spikes in each step
        for i, j in enumerate(inj):
            spike_counts[i,0] = j
            spike_counts[i,1] = (np.sum(np.isfinite(peaks[i,:,:])))/3
        
        spikes = np.where(np.isfinite(peaks)) 
        
        first_spike = spikes[0][0]
        if np.max(spike_counts[:,1]) == 1:
            print('MAX number of AP = 1 for ' + filename[end_filename:-4] + ' Ch ' + str(cell_chan))
            first_spiking_sweep = np.where(spike_counts[:,1]==1)[0][0]
        else:
            first_spiking_sweep=np.where(spike_counts[:,1]>1)[0][0] #where there is more than 1 AP
        
        #fig = plt.figure()
        #ax_spec = APpanelplot(first_spike,inj)
        # for i, spec in enumerate(ax_spec): #for this specific cell 
        #     ax = fig.add_subplot(plt.Subplot(fig, spec))
        #     ax.plot(ch1[:,first_spike+i], lw=0.5, c='grey')
        #     ax.scatter(peaks[first_spike+i, :, 1], 
        #                peaks[first_spike+i, :, 2], marker='+', c='r')
        #     ax.annotate('+'+str(inj[i+first_spike]), (25500,0), (25500,0), color='b', rotation=90)
        
        dir_plots = "Max_Spikes"
        path = os.path.join(filename[:end_filename]+'plots/', dir_plots)
        if os.path.isdir(path) == False: os.mkdir(path)

        #plots all of the sweeps with APs, each fo the detected spikes are marked with a cross
        win = len(inj) - first_spike
        x = math.ceil(np.sqrt(win))
        fig, ax = plt.subplots(x,x,sharex=True, sharey=False,figsize=(6,6))
        for i in range(0,win):
            ax = fig.add_subplot(x,x, i+1)
            if first_spike+i-1 < np.shape(ch1)[1]:
                ax.plot(ch1[:,first_spike+i-1], lw=0.5, c='grey')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.scatter(peaks[first_spike+i-1, :, 1], 
                        peaks[first_spike+i-1, :, 2], marker='+', c='r')
                ax.annotate('+'+str(inj[first_spike+i-1])+' pA', (25500,0), (25500,0), color='b', rotation=90)
        
        fig.patch.set_facecolor('white')
        fig.suptitle('Ch ' +str(cell_chan), fontsize=15)
        fig.tight_layout()
        plt.savefig(path + '/char_spikes_plot' + filename[end_filename:-4]+'_'+str(cell_chan) + '.png')
        plt.close(fig)

        dir_plots = "IV_curve"
        path = os.path.join(filename[:end_filename]+'plots/', dir_plots)
        if os.path.isdir(path) == False: os.mkdir(path)

        #IV cuve; input output function of the cell
        IOfit = np.polyfit(spike_counts[first_spike:,0], spike_counts[first_spike:,1],1)
        IO_slope = IOfit[0]
        plt.figure()
        plt.plot(spike_counts[:,0], spike_counts[:,1])
        plt.gca().set_xlabel('Injection current, pA')
        plt.gca().set_ylabel('Number of spikes')
        Rheobase = inj[first_spike]   
        
        fig.suptitle('Ch ' +str(cell_chan), fontsize=15)
        fig.patch.set_facecolor('white')
        plt.savefig(path + '/char_IV_curve_' + filename[end_filename:-4]+'_'+str(cell_chan) + '.png')
        plt.close()
        
        dir_plots = "AP_props"
        path = os.path.join(filename[:end_filename]+'plots/', dir_plots)
        if os.path.isdir(path) == False: os.mkdir(path)
        #calc TH from second AP in first spiking sweep. TH at vm when dV/dt > 10 mV/ms   
        # TH = threshold;
        
        if np.max(spike_counts[:,1]) == 1:
            peak_loc = np.where(ch1[:,first_spiking_sweep] == peaks[first_spiking_sweep,0,2])[0][0]
            AP1 = ch1[:, first_spiking_sweep][peak_loc - 200:peak_loc+200]
            d1_AP1 = np.diff(AP1)*20
            THloc_all = np.where(d1_AP1[:195] > 10)
            if THloc_all[0].size == 0:
                print("Only 1 SLOW AP found for " + filename[end_filename:-4] + ' Ch : ' + str(cell_chan))
                TH = math.nan
                max_depol = math.nan
                max_repol = math.nan
                APheight = math.nan
                Rheobase = math.nan
                return Vm, max_spikes, Rheobase,APheight, max_depol, max_repol, TH, capacitance, tau
            
            THloc = THloc_all[0][0]
            TH=AP1[THloc-1]
            APheight=peaks[first_spiking_sweep,0,2]-TH #amplitude of AP, the 2nd AP
            max_depol=np.max(d1_AP1[:200]) #how quickly is the voltage change occuring
            max_repol=np.min(d1_AP1[-200:])

            fig = plt.figure()
            plt.plot(AP1)
            plt.scatter(THloc,TH, color = 'red')
            plt.scatter(200,peaks[first_spiking_sweep,0,2], color = 'green')
            fig.suptitle('Ch: ' + str(cell_chan) + ', AP#1' +', TH = ' + str(round(TH,2)) + ', amp = ' + str(round(APheight,2)))
            fig.patch.set_facecolor('white') 
            plt.savefig(path + '/' + filename[end_filename:-4]+'_AP1_'+str(cell_chan) + '.png')
            plt.close()
        else:
            peak_loc = np.where(ch1[:,first_spiking_sweep] == peaks[first_spiking_sweep,1,2])[0][0]
            AP2 = ch1[:, first_spiking_sweep][peak_loc - 200:peak_loc+200]
            d1_AP2 = np.diff(AP2)*20
            THloc_all = np.where(d1_AP2[:195] > 10)
            if THloc_all[0].size == 0:
                print("No TH for SLOW AP in " + filename[end_filename:-4] + ' Ch : ' + str(cell_chan))
                TH = math.nan
                max_depol = math.nan
                max_repol = math.nan
                APheight = math.nan
                Rheobase = math.nan
                return Vm, max_spikes, Rheobase,APheight, max_depol, max_repol, TH, capacitance, tau
            THloc = THloc_all[0][0]
            TH=AP2[THloc-1]
            APheight=peaks[first_spiking_sweep,1,2]-TH #amplitude of AP, the 2nd AP
            max_depol=np.max(d1_AP2[:200]) #how quickly is the voltage change occuring
            max_repol=np.min(d1_AP2[-200:])

            fig = plt.figure()
            plt.plot(AP2)
            plt.scatter(THloc,TH, color = 'red')
            plt.scatter(200,peaks[first_spiking_sweep,1,2], color = 'green')
            fig.suptitle('Ch: ' + str(cell_chan) + ', AP#2' + ', TH = ' + str(round(TH,2)) + ', amp = ' + str(round(APheight,2)))
            fig.patch.set_facecolor('white') 
            plt.savefig(path + '/AP2_' + filename[end_filename:-4]+'_'+str(cell_chan) + '.png')
            plt.close()


        return Vm, max_spikes, Rheobase,APheight, max_depol, max_repol, TH, capacitance, tau



#%%

work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP220413/'

file_list = sorted(os.listdir(work_dir))

filenames = []
for i in range(len(file_list)):
    if file_list[i][-4:] == '.abf': 
        filenames.append(file_list[i])
    
# for i in range(len(filenames)):
#     filename = work_dir + filenames[i]
#     tn.plot_traces(filename)

#%%
#For constructing the filenames 

filenames = []
for i in range(len(file_list)):
    if file_list[i][-4:] == '.abf': 
        filenames.append(file_list[i])
    elif file_list[i][-5:] == '.xlsx': 
        df_rec = pd.read_excel(work_dir + file_list[i], header = 1)
        slice_indx = df_rec.index[df_rec['slice'].notnull()]
        slice_names = df_rec['slice'][slice_indx].tolist()
        index_vc = df_rec.index[df_rec['protocol'] == 'vc'].tolist()
        index_char = df_rec.index[df_rec['protocol'] == 'freq_analyse'].tolist()

index_vm_all=[]
for i in range(len(df_rec['protocol'])):  
    if df_rec['protocol'][i][0:3] == 'vm_':
        index_vm = i
        index_vm_all.append(index_vm)

index_vm_all = np.array(index_vm_all)

# %%
# run for all vm files

for slice in range(len(slice_indx)): #over each slice for a separate file
    slice_name = slice_names[slice]
    df = pd.DataFrame(columns=['file_name', 'channel', 'rec_start_time', 'sweep_num',
            'timepoint','RMP'])
    
    big = index_vm_all[index_vm_all >  slice_indx[slice]]
    if slice == len(slice_indx) -1: dat = big
    else: dat = big[big < slice_indx[slice+1]]

    for files in range(len(dat)):
        filename = filenames[dat[files]]
        

        active_channels = [int(item) for item in input('Which channels were active in' + filename +'(look at plots)').split()]
        for i in range(len(active_channels)): #over channels
            chan = active_channels[i]
            ch, rec_time = load_traces(work_dir+filename, chan)
            for j in range(np.shape(ch)[1]):  #over the sweeps
                sweep_data = ch[:,j]
                # swp_len = np.shape(ch)[0]
                # points = swp_len/1000
                # for k in range(int(points)):
                RMP = np.median(sweep_data) 
                df = df.append({'file_name': filename, 'channel': chan, 'rec_start_time':rec_time,
                    'sweep_num':j+1,'RMP':RMP}, ignore_index=True) 
    df.to_excel(work_dir  + 'data_tables/' + slice_name + '.xlsx')   

#%%

# FOR THE FREQ_ANALYSE FILES 

dir_data = "data_tables"
path = os.path.join(work_dir , dir_data)
if os.path.isdir(path) == False: os.mkdir(path)

#sorting the files based on lab journal entry and number of sweeps
#plots and variables saved in indicated folders (see output)
df = pd.DataFrame(columns=['tissue_source','OP', 'patcher', 'patient_age', 
'filename','slice', 'cell_ch', 'day', 'treatment', 'hrs_incubation','cell_ID', 
'repatch', 'hrs_after_op','AP_heigth','Rheobase', 'TH', 'Vm', 
'capacitance', 'max_depol', 'max_repol', 'max_spikes','membra_time_constant_tau', 
'resting_potential', 'Rs', 'Rin', 'spon_freq', 'spon_ampl', 'mini_freq', 'mini_amp'])																	

# %%

slice_names = ['S1', 'S1', 'S2', 'S2', 'S3', 'S3', 'S4', 'S1_d2']
inj_len = [26,26,21,26,25,26,26,26]


for i in range(len(index_char)):
    char = index_char[i]
    slice = slice_names[i]

    filename_char = work_dir+filenames[char]

    #give inputs with space between entries
    active_channels = [int(item) for item in input('Which channels were active in' + filename_char +'(look at plots)').split()]

    for j in range(len(active_channels)):
        ch = active_channels[j]
        inj = inj_len[i]

        Vm, max_spikes, Rheobase, APheight, max_depol, max_repol, TH, capacitance, tau  = APprops(filename = filename_char, cell_chan = ch, inj_len = inj)
        cellID = filenames[char][:-7]+slice+'c'+str(ch)
        df = df.append({'tissue_source': 'Bielefeld', 'OP':'OP220413', 'patcher':'Verji', 'patient_age':'A', 
        'filename':filenames[char],'slice':slice, 'cell_ch':ch, 'cell_ID':cellID,
        'AP_heigth':APheight,'Rheobase':Rheobase, 'TH':TH, 'Vm' :Vm, 
        'capacitance':capacitance, 'max_depol':max_depol, 'max_repol':max_repol, 
        'max_spikes':max_spikes, 'membra_time_constant_tau':tau}, ignore_index=True)

df.to_excel(work_dir + 'data_tables/'  + 'OP220413_Intrinsic_and_synaptic_properties.xlsx') 
df.to_csv(work_dir + 'data_tables/' + 'OP220413' + '_Intrinsic_and_synaptic_properties.csv')


#%%
#plotting g=funcs


vm_files = [25, 27, 28, 29, 30]
channels = [1,3,4,5,6,7,8]

start_washout = datetime.datetime(2022, 4, 14, 3, 40)

fig1 = plt.figure(figsize=(8,16))
plt.subplots_adjust(hspace=0.5)

for j, ch in enumerate(channels):
    chan = 'Ch' + str(ch)
    RMPs, rec_times = [], []
    for i in vm_files:
        vmfile =  work_dir + filenames[i]

        block = hcf.read_abf(vmfile)
        rec_time = block.rec_datetime

        dt = rec_time - start_washout
        time_after_washout =   datetime.time(hour = int(dt.seconds /3600), minute = int((dt.seconds / 3600 - int(dt.seconds / 3600))*60))

        rec_times.append(time_after_washout.minute)

        vm_dict = load_traces(vmfile)
        ch1 = vm_dict[chan][0]
        RMPs.append(np.median(ch1[:,0]))

    ax = plt.subplot(6, 1, j+1)
    ax.scatter(rec_times, RMPs)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    #ax.spines['left'].set_visible(False)

    ax.set_xticks(ticks = rec_times)
    ax.set_title(chan)
    ax.set_xlabel('Time after high K wash in (min)')
    ax.set_ylabel('RMP (mV)')

plt.show()


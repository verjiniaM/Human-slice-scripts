#%%
from detect_peaks import detect_peaks
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import pandas as pd
import intrinsic_props_and_connectivity.funcs_sorting as sort
import src.intrinsic_props_and_connectivity.funcs_human_characterisation as hcf
import glob


#%%
def rise_fall_speed_ap (filename, channels, inj, treatment, day, swp_step = 3, ap_step = 10):
# this function takes in the filename (.abf file), the cell channel, the injected current steps,
# and the ap_step. If ap_step = 3, then each 3rd action potential is analyzed
# the function outputs the AP peak value, threshold, the rise and the fall speeds (30 to 70 %) from threshold to peak
# additionally plots of each analyzed AP and detected peaks are saved in a folder [file_path**]/plots/APkinetics/
# the datatable with the results is saved in [file_path**]/data_tables/APkinetics/

    data_dict = hcf.load_traces (filename)
    inj = hcf.read_inj(inj)
    end_fn = filename.rfind('/')+1
    file_short = filename[end_fn:]
    OP = filename[filename.find('OP'):filename.find('OP')+8]
    dir_kinetics_plots = sort.make_dir_if_not_existing(filename[:end_fn] + 'plots/',  'AP_kinetics')
    dir_kinetics_data = sort.make_dir_if_not_existing(filename[:end_fn] + 'data_tables/',  'AP_kinetics')

    df = pd.DataFrame(columns=['OP', ' treatment', 'filename', 'day', 'ch', 'swp_num', 'injection', 'AP_num'
    'rise speed', 'fall speed'])
    for ch in channels: 
        key = 'Ch' + str(ch)
        ch1 = data_dict[key][0]

        if np.shape(ch1)[0] < 10_000:
            print('no characterization file found for ' + file_short )
            continue

        for i, j in enumerate(range(0, len(ch1[0]), swp_step)): #loop through all swps
            pks = detect_peaks(ch1[:,j], mph = 0,mpd = 50) # detects the peaks for each sweep
            if len(pks) <= 0:
                print('no peaks found for ' + file_short + ' Ch: ' + str(ch))
                continue

            num_plots =  int(math.ceil(np.sqrt(int(len(pks)/ap_step)+1)))

            fig1 = plt.figure(figsize=(12,12))
            plt.subplots_adjust(hspace=0.5)

            my_labels = {'l1' : 'TH', 'l2' : 'upstroke start', 'l3': 'upstroke end', 'l4': 'downstroke start', 'l5': 'downstroke end'}
            #ap_step = 3 #every 3rd ap is analyzed
            for h in range(0,len(pks),ap_step):

                win = 70
                ap_start = pks[h] - win
                ap_end = pks[h] + win
                AP = ch1[ap_start:ap_end,j]
                peak = AP[win]

                if AP[-1] > peak:
                    h =+ 1
                    ap_start = pks[h] - win
                    ap_end = pks[h] + win
                    AP = ch1[ap_start:ap_end,j]
                    peak = AP[win]

                d1_AP = np.diff(AP)*20 #multiply by 20 to get dv/dt in mV/ms, not for sample
                max_dv_dt = np.where(d1_AP == np.max(d1_AP))[0][0]
                THloc_all = np.where(d1_AP[:max_dv_dt] > 10) #find all instances of fast change before the max dv/dt
                if THloc_all[0].any() == False:
                    print('slow AP for ' + file_short + ' Ch: ' + str(ch))
                    continue
                TH_loc = THloc_all[len(THloc_all) - 1][0] #take the last time the dv/dt crossed 10 mV/ms before the maximum
                TH = AP[TH_loc-1] #-1 because of the length of the derivative
                amp = abs(TH) + peak

                y_30 = amp*0.3
                y_70 = amp*0.7

                AP_positive = AP + abs(TH)

                up_30 = np.where((AP_positive) > y_30)[0][0] #adding the TH so that we have only positive values; no change in the slope parameter
                up_70 = np.where((AP_positive) > y_70)[0][0] 
                x = range(155)
                x1 = x[up_30:up_70]
                fit1 = np.polyfit(x1, AP_positive[x1]*20, 1) # multiply by 20 to have mV/ms
                rise_speed = fit1[0]

                down_30 = np.where((AP_positive[win:] > y_30))[0][-1] + win
                down_70 = np.where((AP_positive[win:] > y_70))[0][-1] + win
                if down_30 == down_70:
                    continue
                x2 = x[down_70:down_30]
                fit2 = np.polyfit(x2, AP_positive[x2]*20, 1) # multiply by 20 to have mV/ms
                fall_speed = fit2[0]

                #add the data to the datafram
                data_to_add = pd.DataFrame({'OP' : OP, 'treatment' : treatment, 'filename': file_short, 'day' : day, "ch" : ch, "swp_num": i, 
                "injection" : inj[i], "AP_num" :  h , "rise speed" : rise_speed, "fall speed" : fall_speed}, index=[0])
                df = pd.concat([df.loc[:], data_to_add]).reset_index(drop=True)

                #plot the result to check
                ax = plt.subplot(num_plots, num_plots, int(h/ap_step+1))
                ax.plot(AP, c = 'k')
                ax.scatter(TH_loc, TH, marker = 'o', c = 'm', label = my_labels['l1'])

                ax.scatter(up_30, y_30-abs(TH), c = 'c', marker = 'p', label = my_labels['l2'])
                ax.scatter(up_70, y_70-abs(TH), c = 'c', marker = 'p', label = my_labels['l3'])
                
                ax.scatter(down_30, y_30-abs(TH), c = 'y', marker = 'p', label = my_labels['l4'])
                ax.scatter(down_70, y_70-abs(TH), c = 'y', marker = 'p', label = my_labels['l5'])
                
                my_labels['l1'], my_labels['l2'], my_labels['l3'], my_labels['l4'], my_labels['l5'] \
                = "_nolegend_", "_nolegend_", "_nolegend_", "_nolegend_", "_nolegend_"

                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                
                ax.set_ylabel("mV")
                ax.title.set_text("AP #" + str(h))
            
            fig1.patch.set_facecolor('white')
            plt.figlegend()
            fig1.suptitle('file ' + file_short + "; ch: " + str(ch)+  "; swp #" + str(i), fontsize=12)

            plt.savefig(dir_kinetics_plots + '/single_AP' + filename[end_fn:-4]+'_'+str(ch) + '.png')
            plt.close()

            fig = plt.figure(figsize=(6,6))
            plt.plot(ch1[:,j])
            for h in range(0,len(pks),ap_step):
                plt.scatter(pks[h], ch1[pks[h],j], marker='+', c='r')
            
            fig.patch.set_facecolor('white') 
            fig.suptitle('file ' + file_short + "; ch: " + str(ch)+  " swp #" + str(i), fontsize=12)

            plt.savefig(dir_kinetics_plots + '/Peaks_per_sweep_' + filename[end_fn:-4]+'_'+str(ch) + '.png')
            plt.close()

    df.to_excel(dir_kinetics_data + '/' + file_short[:-4] + '_AP_kinetics.xlsx', index = False) 


def ap_phase (filename, channels, inj, treatment, swp_step = 5, ap_step = 10):
# this function takes in the filename (.abf file), the cell channel, the injected current steps,
# and the ap_step. If ap_step = 3, then each 3rd action potential is analyzed
# the function outputs the AP peak value, threshold, the rise and the fall speeds (30 to 70 %) from threshold to peak
# additionally plots of each analyzed AP and detected peaks are saved in a folder [file_path**]/plots/APkinetics/
# the datatable with the results is saved in [file_path**]/data_tables/APkinetics/

    data_dict = hcf.load_traces (filename)
    inj = hcf.read_inj(inj)
    end_fn = filename.rfind('/')+1
    OP = filename[filename.find('OP'):filename.find('OP') + 8]
    dir_phase_plots = sort.make_dir_if_not_existing(filename[:end_fn] + 'plots/',  'AP_phase')

    for ch in channels: 
        key = 'Ch' + str(ch)
        ch1 = data_dict[key][0]

        if np.shape(ch1)[0] < 10_000:
            print('no characterization file found for ' + file_short )
            continue
        
        fig1 = plt.figure(figsize=(8,8))
        plt.subplots_adjust(hspace=0.5)

        for i, j in enumerate(range(0, len(ch1[0]), swp_step)): #loop through all swps
            pks = detect_peaks(ch1[:,j], mph = 0,mpd = 50) # detects the peaks for each sweep
            if len(pks) <= 0:
                print('no peaks found for ' + file_short + ' Ch: ' + str(ch))
                continue

            for h in range(0,len(pks),ap_step):
                win = 40
                ap_start = pks[h] - win
                ap_end = pks[h] + win
                AP = ch1[ap_start:ap_end,j]

                ap_change = np.diff(AP)
                ax = plt.subplot(1,1,1)
                ax.plot(AP[:-1], ap_change, label = "AP #" + str(h) + ' at ' + str(inj[i]) + 'pA')

                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                
                ax.set_xlabel("Membrane potential (mV)",fontsize = 16)
                ax.set_ylabel("mV/ ms", fontsize = 16)
                ax.title.set_text('AP phase plot for ' + filename[end_fn:-4] + 'Ch ' + str(ch))

            fig1.patch.set_facecolor('white')
            plt.figlegend()

            plt.savefig(dir_phase_plots + '/AP_phase_' + filename[end_fn:-4]  + '.png')
            plt.close()



# %%
human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0]) 
inj = 'full'

#start analysis from OP220322 
patchers = ['Verji', 'Rosie']
for patcher in patchers:
    data_op = exp_view['OP'][exp_view['patcher'] == patcher]
    for OP in data_op:
        work_dir = sort.get_work_dir(human_dir, OP, patcher)
        intrinsic_file = glob.glob(work_dir + 'data_tables/' + 'QC_passed_' + '*.xlsx')
        if intrinsic_file == []:
            continue
        intrinsic_df = pd.read_excel(intrinsic_file[0])
        for i, fn in enumerate(intrinsic_df['filename'].unique()):
            filename = work_dir + fn
            channels = intrinsic_df['cell_ch'][intrinsic_df['filename'] == intrinsic_df['filename'].unique()[i]]
            treatment = intrinsic_df['treatment'][intrinsic_df['filename'] == fn].tolist()[0]
            day = intrinsic_df['day'][intrinsic_df['filename'] == fn].tolist()[0]
            rise_fall_speed_ap(filename, channels, inj, treatment, day)
            #ap_phase (filename, channels, inj, treatment, swp_step = 5, ap_step = 10)



# %%
# collect AP kinetics folders and calcualte averages 
# correlate with wash out time for TTX
human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0]) 

AP_kinetics_all = pd.DataFrame()
patchers = ['Verji', 'Rosie']
for patcher in patchers:
    data_op = exp_view['OP'][exp_view['patcher'] == patcher]
    for OP in data_op:
        work_dir = sort.get_work_dir(human_dir, OP, patcher)
        AP_kin_files = glob.glob(work_dir + '/data_tables/AP_kinetics/' + '*AP_kinetics.xlsx')
        for fn in AP_kin_files:
            AP_kin = pd.read_excel(fn)
            AP_kin.insert(1, 'patcher', patcher)
            AP_kinetics_all = pd.concat([AP_kinetics_all.iloc[:], AP_kin]).reset_index(drop=True)

#%%

path_file = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/AP_kinetics/2022-11-07AP_kinetics_summary.xlsx'
AP_kinetics_all = pd.read_excel(path_file)

for i in reversed(range(len(AP_kinetics_all))):
    if AP_kinetics_all['rise speed'][i]  <= 0:
        AP_kinetics_all = AP_kinetics_all.drop([AP_kinetics_all.index[i]])
AP_kinetics_all.reset_index(inplace = True, drop = True)

for i in reversed(range(len(AP_kinetics_all))):
    if AP_kinetics_all['fall speed'][i]  < -250:
        AP_kinetics_all = AP_kinetics_all.drop([AP_kinetics_all.index[i]])
AP_kinetics_all.reset_index(inplace = True, drop = True)



#%%



param = 'rise speed'
# summay stats
fig1 = plt.figure(figsize=(13,5))
plt.subplots_adjust(hspace=0.5)
APs = [0, 10]
days = ['D1', 'D2']

colors = ['moccasin', 'moccasin', 'moccasin', 'red','cadetblue', 'mediumpurple']
for j, AP in enumerate(APs):
    ax = plt.subplot(1, 2, j+1)
    g = 0
    treatment = []
    for day in days:
        for i, tr in enumerate(AP_kinetics_all['treatment'].unique()[0:3]):
            df_plot = AP_kinetics_all[(AP_kinetics_all['treatment'] == tr) & (AP_kinetics_all['day'] == day) & (AP_kinetics_all['AP_num'] == AP)]
    
            y1 = df_plot[param].mean()
            ax.scatter(np.linspace(0.5+g, 1.4+g, len(df_plot)),  df_plot[param], c = colors[g], alpha = 0.5)
            ax.scatter(1+g, y1, color = 'k', marker = '_', s = 2000)
            ax.text(0.8+g, y1 + 10, str(round(y1,2)))
            g = g +1
            treatment.append(tr + ' ' +day)
    ax.set_title('AP #' + str(AP), size = 20)
    ax.set_xticks(ticks = [1, 2, 3, 4,5,6], labels = treatment) 
    ax.set_yticks(ticks = [0,100,200,300,400,500,600, 700])
    #ax.set_yticks(ticks = [0,-50, - 100,  -150, -200,-250])

fig1.suptitle('AP_' + param , fontsize = 24)
fig1.patch.set_facecolor('white')
plt.show()

# %%


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




# %%

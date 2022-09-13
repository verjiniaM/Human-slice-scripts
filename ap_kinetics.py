#%%
import human_characterisation_functions as hchf
from detect_peaks import detect_peaks
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import pandas as pd


#%%
def rise_fall_speed_ap (filename, cell_chan, inj, ap_step):
# this function takes in the filename (.abf file), the cell channel, the injected current steps,
# and the ap_step. If ap_step = 3, then each 3rd action potential is analyzed
# the function outputs the AP peak value, threshold, the rise and the fall speeds (30 to 70 %) from threshold to peak
# additionally plots of each analyzed AP and detected peaks are saved in a folder [file_path**]/plots/APkinetics/
# the datatable with the results is saved in [file_path**]/data_tables/APkinetics/

    ch1, sweep_len, block = hchf.load_traces(filename, cell_chan) #in ch1 each sweep is a column
    if inj == "full":
        inj=[-300,-200,-150,-100,-50,0,50,100,150,200,250,300,350,400,450,500,550,
         600,700,800,900,1000,1100,1200,1300,1400]
    else:
        inj = inj

    end_filename = filename.rfind('/')+1
    file_short = filename[end_filename:]

    df = pd.DataFrame(columns=["filename", "cell_chan", "swp_num", "injection", "rise speed", "fall speed"])
    for i, j in enumerate(range(0, len(ch1[0]))): #loop through all swps
        pks = detect_peaks(ch1[:,j], mph=0,mpd=50) # detects the peaks for each sweep
        if len(pks) > 0:
            num_plots =  int(math.ceil(np.sqrt(int(len(pks)/3)+1)))

            fig1 = plt.figure(figsize=(6,6))
            plt.subplots_adjust(hspace=0.5)

            #ap_step = 3 #every 3rd ap is analyzed; can be changes when executign the function
            for h in range(0,len(pks),ap_step):

                ap_start = pks[h] - 100
                ap_end = pks[h] + 100
                AP = ch1[ap_start:ap_end,j]
                d1_AP = np.diff(AP)*20 #multiply by 20 to get dv/dt in mV/ms, not for sample
                max_dv_dt = np.where(d1_AP == np.max(d1_AP))[0][0]
                THloc_all = np.where(d1_AP[:max_dv_dt] > 10) #find all instances of fast change before the max dv/dt
                TH_loc = THloc_all[len(THloc_all) - 1][0] #take the last time the dv/dt crossed 10 mV/ms before the maximum
                TH=AP[TH_loc-1] #-1 because of the length of the derivative

                peak = AP[100]
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

                down_30 = np.where((AP_positive[100:] > y_30))[0][-1]+100
                down_70 = np.where((AP_positive[100:] > y_70))[0][-1]+100
                x2 = x[down_70:down_30]
                fit2 = np.polyfit(x2, AP_positive[x2]*20, 1) # multiply by 20 to have mV/ms
                fall_speed = fit2[0]

                #add the data to the datafram
                data_to_add = pd.DataFrame({'filename': file_short, "cell_chan" : cell_chan, "swp_num": i, 
                "injection" : inj[i], "rise speed" : rise_speed, "fall speed" : fall_speed}, index=[0])
                df = pd.concat([df.loc[:], data_to_add]).reset_index(drop=True)

                #plot the result to check
                ax = plt.subplot(num_plots,num_plots, int(h/3+1))
                ax.plot(AP, c = 'k')
                ax.scatter(TH_loc, TH, marker = 'o', c = 'm')
                #ax.text(TH_loc-5, TH,'Threshold',horizontalalignment='right')

                ax.scatter(up_30, y_30-abs(TH), c = 'c', marker = 'p')
                #ax.text(up_30-5, y_30-abs(TH), 'rise start',horizontalalignment='right')
                ax.scatter(up_70, y_70-abs(TH), c = 'c', marker = 'p')
                #ax.text(up_70-5, y_70-abs(TH), 'rise end',horizontalalignment='right')
                
                ax.scatter(down_30, y_30-abs(TH), c = 'y', marker = 'p')
                #ax.text(down_30+5, y_30-abs(TH), 'down start')
                ax.scatter(down_70, y_70-abs(TH), c = 'y', marker = 'p')
                #ax.text(down_70+5, y_70-abs(TH), 'down end')
                
                ax.set_ylabel("mV")
                ax.title.set_text("AP #" + str(h))
        
            fig1.patch.set_facecolor('white')
            fig1.suptitle('file ' + file_short + "; ch: " + str(cell_chan)+  "; swp #" + str(i), fontsize=12)

            dir_plots = "APkinetics"
            path = os.path.join(filename[:end_filename]+'plots/', dir_plots)
            if os.path.isdir(path) == False: os.mkdir(path) #if the path doesn't exist, it creates it 
            plt.savefig(path + '/single_AP' + filename[end_filename:-4]+'_'+str(cell_chan) + '.png')
            plt.close()

            fig = plt.figure(figsize=(6,6))
            plt.plot(ch1[:,j])
            for h in range(0,len(pks),ap_step):
                plt.scatter(pks[h], ch1[pks[h],j], marker='+', c='r')
            
            fig.patch.set_facecolor('white') 
            fig.suptitle('file ' + file_short + "; ch: " + str(cell_chan)+  " swp #" + str(i), fontsize=12)

            plt.savefig(path + '/Peaks_per_sweep_' + filename[end_filename:-4]+'_'+str(cell_chan) + '.png')
            plt.close()

    path_data = os.path.join(filename[:end_filename]+'data_tables/', dir_plots)
    if os.path.isdir(path_data) == False: os.mkdir(path_data) #if the path doesn't exist, it creates it 
    df.to_excel(path_data + '/' + file_short[:-4] + '_AP_kinetics.xlsx') 
    df.to_csv(path_data +  '/' + file_short[:-4] + '_AP_kinetics.csv')



# %%

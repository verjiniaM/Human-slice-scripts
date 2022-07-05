
#%%
import os
import human_characterisation_functions as hcf
import human_synaptic_functions1 as hsf
import trace_names_check as tn
import pandas as pd
import glob

#loading the updated experiments_overview and the old summary_data_table
human_folder = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'

exp_view = pd.read_excel(human_folder + 'experiemnts_overview.xlsx') 
latest_sum_data_table = glob.glob(human_folder + 'summary_data_tables/' + '*.xlsx')
summary_data_table = pd.read_excel(latest_sum_data_table[0])

#%%
#taking the unique OP IDs and updating the full OP list

last_update_op_list = summary_data_table.OP.unique().tolist()
newest_op_list = exp_view.OP.unique().tolist()
all_OPs = last_update_op_list + newest_op_list
op_to_analyse = [i for i in all_OPs if all_OPs.count(i)==1]

#%%
#Decide which OPs need to be analyzed 

for i in range(len(op_to_analyse)):
    OP = op_to_analyse[i]
    y_or_n = input('Do you want to start analysis of ' + OP + '(y/n)?') 
    if y_or_n == "y":
        #input some functions
        OP_folder = OP + '/'
        tissue_source = input('Tissue source for ' + OP + '(Bielefeld, Mitte, Virchow): ')
        patcher = input('Patcher ' + OP + '(Verji or Rosie): ')
        age = input('Patient age ' + OP + '(if not known: A for adult, J for juvenile): ')
        
        print('starting analysis for '+ OP)
        
        if patcher == 'Verji': work_dir = human_folder + 'data_verji/'+ OP_folder 
        if patcher == 'Rosie': work_dir = human_folder + 'data_rosie/'+ OP_folder 
        file_list = sorted(os.listdir(work_dir))

        filenames = []
        for file in range(len(file_list)):
            #pclamp files
            if file_list[file][-4:] == '.abf': 
                filenames.append(file_list[file])
            #sorting the labbook entries and number of sweeps
            elif (file_list[file][-5:] == '.xlsx' and file_list[file][:2] == '~$'):
                print('Close the excel file')
                input()
            elif (file_list[file][-5:] == '.xlsx' and file_list[file][:2] != '~$'): 
                df_rec = pd.read_excel(work_dir + file_list[file], header = 1)
                slice_indx = df_rec.index[df_rec['slice'].notnull()]
                slice_names = df_rec['slice'][slice_indx].tolist()
                index_vc = df_rec.index[df_rec['protocol'] == 'vc'].tolist()
                index_vm = df_rec.index[df_rec['protocol'] == 'vm_mouse'].tolist()
                index_char = df_rec.index[df_rec['protocol'] == 'freq analyse'].tolist()
                index_vc_end = df_rec.index[df_rec['protocol'] == 'vc_end'].tolist()
                index_spontan = df_rec.index[df_rec['protocol'] == 'vm'].tolist()
                index_mini = df_rec.index[df_rec['protocol'] == 'minis'].tolist()
        
        #adjusting the slice names 
        new_slice_names = []
        for i in range(len(slice_names)):
            if i < len(slice_names)-1:
                new_slice_names.append([slice_names[i]]*(slice_indx[i+1]-slice_indx[i]))
            else :
                new_slice_names.append([slice_names[i]]*(15))
        slice_names = flat_list = [x for xs in new_slice_names for x in xs]

        #creating a dir to save plots (if not existing)
        dir_plots = "plots"
        path = os.path.join(work_dir, dir_plots)
        if os.path.isdir(path) == False: os.mkdir(path)

        #check if the traces dir is empty and only then plot the mimiddle sweep for each filename
        #if os.path.isdir(path) == False: os.mkdir(path)
        traces_folder =  os.path.join(path, "traces/")
        if os.path.isdir(traces_folder) == 1:
            print("skipping plotting")
        else: 
            for rec in range(len(filenames)):
                filename = work_dir + filenames[rec]
                tn.plot_traces(filename)

        print('number of vc files: ' + ' '.join(str(h) for h in index_vc))
        print('number of vm files: ' + ' '.join(str(h) for h in index_vm))
        print('number of characterization files: ' + ' '.join(str(g) for g in index_char))
        print('number of vc_end files: ' + str(index_vc_end))
       
        if len(index_char) != len(index_vc): 
            print('Fix protocol names. Unequal number of VC and freq analyse protocols')
            index_vc = [int(item) for item in input('Vc files corresponding to characterization files for ' + OP +' (input with spaces in between)').split()]
            index_char = [int(item) for item in input('corresponding characterization files for ' + OP +' (input with spaces in between)').split()]

        proceed_y_n = input("do all traces correspond to specified filenames in lab book for " + OP +  "(y/n)?")
        if proceed_y_n == 'n': 
            print('correct the lab book entries. Continuing to next OP')
            continue

        dir_data = "data_tables"
        path = os.path.join(work_dir , dir_data)
        if os.path.isdir(path) == False: os.mkdir(path)

        #creating the dataframe
        df = pd.DataFrame(columns=['tissue_source','OP', 'patcher', 'patient_age', 
        'filename','slice', 'cell_ch', 'day', 'treatment', 'hrs_incubation','cell_ID', 
        'repatch', 'hrs_after_op','AP_heigth','Rheobase', 'TH', 'Vm', 
        'capacitance', 'max_depol', 'max_repol', 'max_spikes','membra_time_constant_tau', 
        'resting_potential', 'Rs', 'Rin', 'spon_freq', 'spon_ampl', 'mini_freq', 'mini_amp'])	
 
        for i in range(len(index_vc)):
            vc = index_vc[i]
            vm = index_vm[i]
            char = index_char[i]
            slice = slice_names[vc]

            filename_vc = work_dir+filenames[vc]
            filename_vm = work_dir+filenames[vm]
            filename_char = work_dir+filenames[char]

            #give inputs with space between entries
            active_channels = [int(item) for item in input('Channels used in ' + filenames[vc] +'(input with spaces in between)').split()]

            for j in range(len(active_channels)):
                ch = active_channels[j]
                hsf.get_holding_measures(filename_vc, ch)
                Rs, Rin = hcf.access_resistance(filename_vc, ch) 
                cellID = filenames[vc][:-7]+slice+'c'+str(ch)
                
                resting_mem = hcf.restingmembrane(filename_vm, ch)
                Vm, max_spikes, Rheobase, APheight, max_depol, max_repol, TH, capacitance, tau  = hcf.APproperties(filename_char, ch, "full")
                cellID = filenames[char][:-7]+slice+'c'+str(ch)
                data_to_add = pd.DataFrame({'tissue_source': tissue_source, 'OP':OP[:-1], 'patcher':patcher, 'patient_age':age, 
                'filename':filenames[char],'slice':slice, 'cell_ch':ch, 'cell_ID':cellID,
                'Rs':Rs, 'Rin':Rin,'resting_potential': resting_mem,'AP_heigth':APheight,'Rheobase':Rheobase, 'TH':TH, 'Vm' :Vm, 
                'capacitance':capacitance, 'max_depol':max_depol, 'max_repol':max_repol, 
                'max_spikes':max_spikes, 'membra_time_constant_tau':tau}, index=[0])
                df = pd.concat([df.loc[:], data_to_add]).reset_index(drop=True)


        df.to_excel(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.xlsx') 
        df.to_csv(work_dir + 'data_tables/' + OP[:-1] + '_Intrinsic_and_synaptic_properties.csv')

        print("Analysis for " + OP + 'is complete. Double check the excel table and exclude recordings if necessary.')

        for file in range(len(file_list)):
            #pclamp files
            if  (file_list[file][-5:] == '.xlsx' and file_list[file][:2] != '~$'): 
                df_rec = pd.read_excel(work_dir + file_list[file], header = 1)
                index_vc = df_rec.index[df_rec['protocol'] == 'vc'].tolist()
                index_vc_end = df_rec.index[df_rec['protocol'] == 'vc_end'].tolist()
                index_spontan = df_rec.index[df_rec['protocol'] == 'vm'].tolist()
                index_min = df_rec.index[df_rec['protocol'] == 'minis'].tolist()

        print('VC files:     ' + ' '.join(str(h) for h in index_vc))
        print('VC_end files: ' + ' '.join(str(h) for h in index_vc_end))
        print('Spontaneous files: ' + str(index_spontan))
        print('Mini files: ' + str(index_mini))

        proceed = input('Proceeding with QC steps? All mvc, vc_end, mini, and spontaneous files detected? (y/n)')

        if proceed == 'n': 
            
            input('waiting to continue. Do not forget to close the excel file, if you change something')
            print('updating the excel file')

            print('Fix protocol names. Unequal number of VC and freq analyse protocols')
            index_vc = [int(item) for item in input('Vc files with appropriate vc_end files for ' + OP +' (input with spaces in between)').split()]
            index_vc_end = [int(item) for item in input('Vc_end files corresponding to characterization files for ' + OP +' (input with spaces in between)').split()]

            print('VC files:     ' + ' '.join(str(h) for h in index_vc))
            print('VC_end files: ' + ' '.join(str(h) for h in index_vc_end))
            print('Spontaneous files: ' + str(index_spontan))
            print('Mini files: ' + str(index_mini))

            if len(index_vc) != len(index_vc_end): 
                print('Fix protocol names. Unequal number of VC and vc_end files')
                index_vc = [int(item) for item in input('Vc files corresponding to vc_end files for ' + OP +' (input with spaces in between)').split()]
                index_vc_end = [int(item) for item in input('Vc_end files ' + OP +' (input with spaces in between)').split()]

        #date_frame for quality control - change in Rin, Rs 
        df_qc = pd.DataFrame(columns=['OP','patcher', 'filename', 'slice', 'cell_ch', 'Rs_start', 'Rin_start', 'Rs_end', 'Rin_end', 
        'chage_rs', 'change_rin'])

        for k in range(len(index_vc_end)):
            vc = index_vc[k]
            vc_end = index_vc_end[k]
            slice = slice_names[vc]

            filename_vc = work_dir+filenames[vc]
            filename_vc_end = work_dir+filenames[vc_end]

            active_channels = [int(item) for item in input('Channels used in ' + filenames[vc] +'(input with spaces in between)').split()]

            for l in range(len(active_channels)):
                ch = active_channels[l]
                Rs, Rin = hcf.access_resistance(filename_vc, ch) 
                Rs_end, Rin_end = hcf.access_resistance(filename_vc_end, ch) 
                cellID = filenames[vc][:-7]+slice+'c'+str(ch)

                data_to_add = pd.DataFrame({'OP':OP[:-1], 'patcher':patcher, 'filename':filenames[vc],'slice':slice, 
                'cell_ch':ch, 'Rs_start': Rs, 'Rin_start': Rin, 'Rs_end': Rs_end, 'Rin_end': Rin_end, 
                'chage_rs': Rs-Rs_end, 'change_rin':Rin-Rin_end}, index=[0])
                df_qc = pd.concat([df_qc.loc[:], data_to_add]).reset_index(drop=True)

        df_qc.to_excel(work_dir + 'data_tables/' + OP[:-1] + '_QC_measures_rs.xlsx') 
        df_qc.to_csv(work_dir + 'data_tables/' + OP[:-1] + '_QC_measures_rs.csv')

        #quality check for minis
        #calcualted the signal drift and sorts suitable/unsuitable recordings

        record_sorting = pd.DataFrame(columns = ['OP','patcher', 'filename', 'mini/spontan?', 'cell_ch','swps_to_analyse', 'swps_to_discard', 
        'min_vals_each_swp', 'max_vals_each_swp', 'drift'])

        for u in range(len(index_spontan)):
            spontan = index_spontan[u]
            slice = slice_names[spontan]

            filename_spontan = work_dir+filenames[spontan]
            
            active_channels = [int(item) for item in input('Channels used in ' + filenames[spontan] +'(input with spaces in between)').split()]

            for q in range(len(active_channels)):
                ch = active_channels[q]
                min_vals, max_vals, good_swps, bad_swps, drift = hcf.rec_stability (filename_spontan, ch , 50)

                data_to_add = pd. DataFrame({'OP':OP[:-1], 'patcher':patcher, 'filename':filenames[spontan], 'mini/spontan?':'spontan',
                'slice':slice, 'cell_ch':ch, 'swps_to_analyse':[good_swps], 'swps_to_discard':[bad_swps], 'min_vals_each_swp':[min_vals], 
                'max_vals_each_swp':[max_vals], 'drift':[drift]})
                record_sorting = pd.concat([record_sorting.loc[:], data_to_add]).reset_index(drop=True)



        for h in range(len(index_mini)):
            mini = index_mini[h]
            slice = slice_names[mini]

            mini_file = work_dir+filenames[mini]

            active_channels = [int(item) for item in input('Channels used in ' + filenames[mini] +'(input with spaces in between)').split()]

            for b in range(len(active_channels)):
                ch = active_channels[b]
                min_vals, max_vals, good_swps, bad_swps, drift = hcf.rec_stability (mini_file, ch , 50)

                data_to_add = pd. DataFrame({'OP':OP[:-1], 'patcher':patcher, 'filename':filenames[mini], 'mini/spontan?':'minis',
                'slice':slice, 'cell_ch':ch, 'swps_to_analyse':[good_swps], 'swps_to_discard':[bad_swps], 'min_vals_each_swp':[min_vals], 
                'max_vals_each_swp':[max_vals], 'drift':[drift]})
                record_sorting = pd.concat([record_sorting.loc[:], data_to_add]).reset_index(drop=True)

        record_sorting.to_excel(work_dir + 'data_tables/' + OP + '_QC_measures_events.xlsx') 
        record_sorting.to_csv(work_dir + 'data_tables/' + OP + '_QC_measures_events.csv')

        print('The data for ' + OP + 'is fully analyzed and saved in ' + work_dir + 'data_tables/')

    # %%

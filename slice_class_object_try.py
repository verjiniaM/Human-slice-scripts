import pandas as pd
import sorting_functions as sort
import human_characterisation_functions as hcf
import glob

class slice_meta :
    def __init__(self, exp_type, name, treatment, patcher, day, active_chans):

        self.exp_type = exp_type
        self.name = name
        self.treatment = treatment
        self.patcher = patcher #Verji, Rosie
        self.day = day
        self.channels = active_chans

        #incides of the filenames from the work_dir
        self.vc_files = []
        self.RMP = []
        self.char = []
        self.ramp = []
        self.con_screen = []
        self.con_screen_VC = []
        self.spontan = []
        self.minis = []
    
    def filenames_rec_type (self, df_slic_fn):
        self.vc_files.append(df_slic_fn['files'].iloc[indices_dict['vc']].loc[df_slic_fn['slice'] == self.name].tolist())
        self.RMP.append(df_slic_fn['files'].iloc[indices_dict['resting']].loc[df_slic_fn['slice'] == self.name].tolist())
        self.char.append(df_slic_fn['files'].iloc[indices_dict['freq analyse']].loc[df_slic_fn['slice'] == self.name].tolist())
        self.ramp.append(df_slic_fn['files'].iloc[indices_dict['ramp']].loc[df_slic_fn['slice'] == self.name].tolist())
        self.con_screen.append(df_slic_fn['files'].iloc[indices_dict['con_screen']].loc[df_slic_fn['slice'] == self.name].tolist())
        self.con_screen_VC.append(df_slic_fn['files'].iloc[indices_dict['con_screen_VC']].loc[df_slic_fn['slice'] == self.name].tolist())
        self.spontan.append(df_slic_fn['files'].iloc[indices_dict['spontan']].loc[df_slic_fn['slice'] == self.name].tolist())
        self.minis.append(df_slic_fn['files'].iloc[indices_dict['minis']].loc[df_slic_fn['slice'] == slice_name].tolist())
        self.vc_minis.append(df_slic_fn['files'].iloc[indices_dict['vc_mini']].loc[df_slic_fn['slice'] == slice_name].tolist())
        if self.vc_files == []:
            self.vc_files == df_slic_fn['files'].iloc[indices_dict['vc_mini']].loc[df_slic_fn['slice'] == slice_name].tolist()

    def new_cell_IDs(self, filenames):
        patcher_dict = {'Verji':'vm', 'Rosie': 'rs'}
        date_part = filenames[0][:-7]

        cell_IDs = []
        for ch in self.channels:
            cellID = patcher_dict[self.patcher] + date_part + \
                self.name + 'c' + str(ch)
            cell_IDs.append(cellID)
        return cell_IDs

#%%

#%%
OP = 'OP230417'
patcher = 'Verji'
tissue_source = 'Virchow'
inj = 'full'
age = '4.5'

#loading the updated experiments_overview and the old summary_data_table
human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
results_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/'

exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0]) 

#%%
work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)

df_slic_fn = pd.DataFrame({'slice': slice_names[:len(filenames)], 'files':filenames})

''' 
Each slice class is defined by one VC fiels and whatever else
'''

all_slices = sorted(list(set(slice_names)))

cortex_out_time = exp_view['cortex_out'].loc[exp_view['OP'] == OP].tolist()[0]

df_OP = pd.DataFrame(columns=['filename', 'slice', 'cell_ch', 'cell_ID', 'day', 'treatment',
    'hrs_incubation', 'repatch', 'hrs_after_OP', 'Rs', 'Rin', 'resting_potential', 'max_spikes',
    'Rheobase', 'AP_heigth', 'TH', 'max_depol', 'max_repol', 'membra_time_constant_tau', 'capacitance'])

for slice_name in all_slices:
    treatment = input('Treatment (Ctrl, high K, or TTX) for ' + slice_name + ' ' + OP)
    day = 'D1'
    if slice_name[-2:] == 'D2': 
        day = 'D2'

    vc_per_slice = df_slic_fn.iloc[indices_dict['vc']].loc[df_slic_fn['slice'] == slice_name]

    for vc_file in vc_per_slice['files']:
        active_chans = [int(item) for item in input('Used channels in ' + vc_file + ' ' + OP).split()]

        slice1 = slice_meta(slice_name, treatment, patcher, day, active_chans)
        slice1.filenames_rec_type (df_slic_fn)

        cell_IDs = slice1.new_cell_IDs(filenames)

        if len(slice1.vc_files[0]) != len(slice1.RMP[0]):

            time_after_op = sort.get_time_after_OP(work_dir + fn_vc, cortex_out_time)

            for i, fn_vc in enumerate(slice1.vc_files[0]):
                time_after_op = sort.get_time_after_OP(work_dir + fn_vc, cortex_out_time)
                Rs, Rin = hcf.get_access_resistance(work_dir + fn_vc, active_chans) 

                params1_df = pd.DataFrame({'filename': slice1.vc_files[0][i], 'slice' : slice1.name, 'cell_ch': slice1.channels,
                'hrs_after_OP' : time_after_op, 'cell_ID': cell_IDs, 'day' : slice1.day , 
                'treatment': slice1.treatment, 'Rs' : Rs, 'Rin': Rin})

                df_OP = pd.concat([df_OP.loc[:], params1_df]).reset_index(drop=True)
                df_name = '_VC_props.xlsx'
            print('Not usual repatch experiments. Saving only a vc dataframe at the end')

        else:
            #calculate parameters for repatch recordings
            for i, fn_vc in enumerate(slice1.vc_files[0]):
                time_after_op = sort.get_time_after_OP(work_dir + fn_vc, cortex_out_time)
                Rs, Rin = hcf.get_access_resistance(work_dir + fn_vc, active_chans) 

                RMPs = hcf.get_RMP(work_dir + slice1.RMP[0][i], active_channels)
                params1_df = pd.DataFrame({'filename': slice1.char[0], 'slice' : slice1.name, 'cell_ch': slice1.channels,
                    'hrs_after_OP' : time_after_op, 'cell_ID': cell_IDs, 'day' : slice1.day , 
                    'treatment': slice1.treatment, 'Rs' : Rs, 'Rin': Rin, 'resting_potential': RMPs })

                charact_params  = hcf.all_chracterization_params(work_dir + slice1.char[i], active_channels, inj)
                df_char = pd.DataFrame.from_dict(charact_params)

                df_to_add = pd.concat([params1_df, df_char], axis = 1)
                df_OP = pd.concat([df_OP.loc[:], df_to_add]).reset_index(drop=True)

                #plotting function
                plotting_funcs.plot_vc_holding (filename_vc, active_channels)
                plotting_funcs.plots_for_charact_file(filename_char, active_channels, inj)
                df_name = '_Intrinsic_and_synaptic_properties.xlsx'
                print('Intrinsic properties DataFrame for  ' + OP + ' saved successfully. ' + '\n' + 'Exclude recordings if necessary.')

tissue = pd.Series(tissue_source).repeat(len(df_OP))
OPs = pd.Series(OP).repeat(len(df_OP))
researcher = pd.Series(patcher).repeat(len(df_OP))
patient_age = pd.Series(age).repeat(len(df_OP))
series_df = pd.DataFrame({'tissue_source': tissue, 'OP': OPs, 'patcher': researcher, 'patient_age': patient_age}).reset_index(drop=True)

df_intrinsic = pd.concat([series_df, df_OP], axis = 1)

df_intrinsic.to_excel(work_dir + 'data_tables/TRY_' + OP + df_name, index=False) 



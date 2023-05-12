import pandas as pd
import sorting_functions as sort
import human_characterisation_functions as hcf
import glob

class slice_meta :
    def __init__(self, name, treatment, patcher, day, active_chans):

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

2
3
5
7
8
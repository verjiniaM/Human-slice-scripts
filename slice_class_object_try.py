import pandas as pd
import sorting_functions as sort

class slice_meta :
  def __init__(self, name, treatment, patcher, day, active_chans, 
  vc_files, RMP, char, ramp, con_screen, con_screen_VC,
  spontan, minis):
    self.name = name
    self.treatment = treatment
    self.patcher = patcher #Verji, Rosie
    self.day = day
    self.channels = active_chans
    #incides of the filenames from the work_dir
    self.vc_files = vc_files
    self.RMP = RMP
    self.char = char
    self.ramp = ramp
    self.con_screen = con_screen
    self.con_screen_VC = con_screen_VC
    self.spontan = spontan
    self.minis = minis

  def new_cell_IDs(self, filenames):
    patcher_dict = {'Verji':'vm', 'Rosie': 'rs'}
    date_part = filenames[0][:-7]

    cell_IDs = []
    for ch in self.channels[-1]:
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

#%%
work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)

df_slic_vc = pd.DataFrame({'slice': slice_names[:len(filenames)], 'files':filenames})


all_slices = sorted(list(set(slice_names)))
#for slice_name in all_slices:
slice_name = all_slices[0]
treatment = input('Treatment (Ctrl, high K, or TTX) for ' + slice_name + ' ' + OP)
day = 'D1'
if slice_name[-2:] == 'D2': 
    day = 'D2'

vc_per_slice = len(df_slic_vc.iloc[indices_dict['vc']].loc[df_slic_vc['slice'] == slice_name])

if vc_per_slice == 1:
    active_chans = [[int(item) for item in input('Used channels in ' + fn + ' ' + OP).split()]]
elif vc_per_slice > 1:
    active_chans = []
    for i, fn in enumerate(df_slic_vc['files'].iloc[indices_dict['vc']].loc[df_slic_vc['slice'] == slice_name]):
        active_chans.append([int(item) for item in input(slice_name + 'Used channels in ' + fn + ' ' + OP).split()])

vc_files = df_slic_vc['files'].iloc[indices_dict['vc']].loc[df_slic_vc['slice'] == slice_name].tolist()
RMP = df_slic_vc['files'].iloc[indices_dict['resting']].loc[df_slic_vc['slice'] == slice_name].tolist()
char = df_slic_vc['files'].iloc[indices_dict['freq analyse']].loc[df_slic_vc['slice'] == slice_name].tolist()
ramp = df_slic_vc['files'].iloc[indices_dict['ramp']].loc[df_slic_vc['slice'] == slice_name].tolist()
con_screen = df_slic_vc['files'].iloc[indices_dict['con_screen']].loc[df_slic_vc['slice'] == slice_name].tolist()
con_screen_VC = df_slic_vc['files'].iloc[indices_dict['con_screen_VC']].loc[df_slic_vc['slice'] == slice_name].tolist()
spontan = df_slic_vc['files'].iloc[indices_dict['spontan']].loc[df_slic_vc['slice'] == slice_name].tolist()
minis = df_slic_vc['files'].iloc[indices_dict['minis']].loc[df_slic_vc['slice'] == slice_name].tolist()

slice1 = slice_meta(slice_name, treatment, patcher, day, active_chans, 
vc_files, RMP, char, ramp, con_screen, con_screen_VC,
spontan, minis)










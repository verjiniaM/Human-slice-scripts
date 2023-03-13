import sorting_functions as sort
import human_characterisation_functions as hcf
import pandas as pd 
import numpy as np
import plot_intrinsic_props as plot_intr


human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
OP = 'OP221027'
patcher = 'Verji'
file_out = '_meta_active_chans.json'

work_dir, filenames, indices_dict, slice_names = sort.get_OP_metadata(human_dir, OP, patcher)
OP_meta = sort.from_json(work_dir, OP, file_out)

treatments = OP_meta[0]['treatment']
active_chans = OP_meta[0]['active_chans']

IFF_all = np.zeros((1,26))
for i, char in enumerate(indices_dict['freq analyse']):
    fn = work_dir + filenames[char]
    channels = active_chans[i]
    inj = 'full'
    IFF_dict = hcf.get_initial_firing_rate(fn, channels, inj = 'full')
    IFF_all = np.concatenate((IFF_all, IFF_dict['IFF']))

df = pd.DataFrame(IFF_all)
df.to_excel(work_dir + 'data_tables/IFF_all.xlsx',index=False)

plot_intr.plot_IFF_distribution(IFF_all)

plot_intr.plot_IFF_averages_for_I(IFF_all, IFF_dict['inj'])


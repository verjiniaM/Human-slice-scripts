import pandas as pd
import glob
import human_characterisation_functions as hcf
import funcs_for_results_tables as get_results
import plot_intrinsic_props as plot_intr


human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
exp_view = pd.read_excel(glob.glob(human_dir + '*experiments_overview.xlsx')[0]) 
exp_view_IFF = exp_view[exp_view['repatch'] == 'yes']

# OPs = exp_view_IFF['OP']
# patchers = exp_view_IFF['patcher'].tolist()
# for i, OP in enumerate(OPs):
#     patcher = patchers[i]
#     get_results.create_IFF_data_table(OP, patcher, file_out = '_meta_active_chans.json', inj = 'full', human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/')

IFF_collected = get_results.collect_IFF_dfs(human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/')
intrinsic_df = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/summary_data_tables/intrinsic_properties/2022-12-07_collected.xlsx')
IFF_all, IFF_repatch = get_results.get_QC_for_IFF_df(intrinsic_df, IFF_collected)

plot_intr.plot_IFF_distribution(IFF_repatch, 'repatch')
plot_intr.plot_IFF_distribution(IFF_all, 'all')


# plot_intr.plot_IFF_averages_for_I(IFF_all, IFF_dict['inj'])
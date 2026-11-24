import pandas as pd
import numpy as np



work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/'
save_plots_dir = f'{work_dir}results/human/connectivity_paper/plots/figs_draft2/'
VM_data_dir = f'{work_dir}results/human/connectivity_paper/data/'
stats_dir = f'{work_dir}results/human/connectivity_paper/stats/for_plotting/'


df_con_VM_etypes_OLD = pd.read_csv(f'{VM_data_dir}all_connect_cross_with_intr_with_etypes.csv')
all_event_props = pd.read_csv(f'{VM_data_dir}all_connections_QC_checked_complete.csv')
df_VM_etypes = pd.read_csv(f'{VM_data_dir}VM_intr_with_etype_QCed.csv')

#%%

# FIND CONNECTION PAIRS THAT WERE EXCLUDED IN THE LAST QC
keys = ['cell_IDs_match_pre', 'cell_IDs_match_post']

old_pairs = df_con_VM_etypes_OLD[keys].dropna().drop_duplicates()
all_pairs = all_event_props[keys].dropna().drop_duplicates()

missing_pairs = old_pairs.merge(all_pairs, on=keys, how='left', indicator=True)
missing_pairs = missing_pairs.loc[missing_pairs['_merge'] == 'left_only', keys]

print('Old pairs missing from all_event_props:', len(missing_pairs))
print(missing_pairs)


# USING MOST UPDATED DF, CREATE CONNECTIONS WITH ETYPES
df_con_VM_etypes = all_event_props.copy()

for side in ['pre', 'post']:
    df_con_VM_etypes = df_con_VM_etypes.merge(
        df_VM_etypes[['cell_IDs_match', 'e_type']].drop_duplicates('cell_IDs_match'),
        left_on=f'cell_IDs_match_{side}', right_on='cell_IDs_match',
        how='left', validate='m:1'
    ).rename(columns={'e_type': f'e_type_{side}'}).drop(columns='cell_IDs_match')

# keep only where post cell
df_con_post_etypes = df_con_VM_etypes.dropna(subset=['e_type_post']).copy()

# keep where both
df_con_pre_post_etypes = df_con_post_etypes.dropna(subset=['e_type_pre']).copy()





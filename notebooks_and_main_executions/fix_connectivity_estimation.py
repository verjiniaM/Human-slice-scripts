import pandas as pd
import preparation_paper_figures as paper_figs
import connectivity_estimation as con_estm

#%%
# intr_all = paper_figs.get_intrinsic_df_con_paper(min_holding = -250)
intr_all = pd.read_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+ \
                            '/connectivity_paper/data/con_analysis_clean_all/all_intr_props_all_treatments.csv')

df_all_connections = con_estm.get_all_connectivity_data(intr_all)
df_all_connections = con_estm.get_match_cell_id_connectivity(df_all_connections)

# analyzing duplicates
df_all_connections[df_all_connections.duplicated(keep = False, subset=['OP', 'slice', 'day', 'chan_pre', 'chan_post'])]


op_slice_day_combos = intr_all[['OP', 'patcher','slice', 'day']].drop_duplicates()
# in connect_dataframe keep only these
df_all_connections = df_all_connections.merge(
    op_slice_day_combos,
    on=['OP', 'slice', 'patcher', 'day'],
    how='inner'
)
df_all_connections = df_all_connections.dropna(subset=['Amp 1']).reset_index(drop=True)
df_all_connections.isna().sum()

# in the following, only cells that are in intr_all are included
df_con_percentage_only_intr, num_alone_cell_all = con_estm.cross_check_intr_con_df(intr_all, df_all_connections)
con_estm.quick_dirty_plot_connectivity(df_con_percentage_only_intr, 'none', save_dir = False)


# # save the results
    df_all_connections.to_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+ \
                                '/connectivity_paper/data/connection_properties_cross_checked_with_intr.csv')


    df_con_percentage_only_intr.to_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/'+ \
                                '/connectivity_paper/data/con_percentages_from_cross_checked_with_intr.csv')



#%%
# connectivitty estimation without intrinsic
# doesn't really work

df_con_jsons = con_estm.get_connections_df_from_jsons()
df_con_jsons = con_estm.add_connectivity_percentage(df_con_jsons)
df_con_jsons = df_con_jsons.drop_duplicates(keep='first').reset_index(drop=True)

df_con_jsons_found = df_con_jsons.merge(
    df_con_percentage_only_intr[['OP', 'patcher', 'slice', 'day', 'found_connections']],
    on=['OP', 'patcher', 'slice', 'day'],
    how='right'
)

df_con_jsons_found['found_connections'] = df_con_jsons_found['found_connections_y']

# removing columns and recalcualting
df_con_jsons_found_correct = df_con_jsons_found.copy()
df_con_jsons_found_correct = df_con_jsons_found_correct.drop(columns=['con_percentage', 'con_percentage_fixed'])

con_estm.quick_dirty_plot_connectivity(df_con_jsons_found_correct, 'none', save_dir = False)

import pandas as pd
import funcs_events_post_analysis as event_funcs
import funcs_plot_intrinsic_props as plot_intr
import funcs_for_results_tables as get_results
import datetime
import glob


# create metadata tables]
# !!!! This function needs rewriting but can be use interactivelly
# get_results.prapare_for_event_analysis(op_to_analyse = '')

# post analysis
spontan_new_df = event_funcs.post_events_analysis_add_metadata('spontan')
metadata_df = pd.read_excel('/Users/verjim/spontaneous-postsynaptic-currents-detection/metadata/metadata_S2_OP240215.xlsx')

spontan_QC_df = event_funcs.QC_post_event_analysis(spontan_new_df, metadata_df ,min_event_size = 3, min_hrs = 16)

n_nums_spontan = event_funcs.get_n_nums_per_day(spontan_QC_df, ['Ctrl', 'high K'])
event_funcs.plot_event_by_day_spontan(spontan_QC_df, n_nums_spontan, 
'repatch', save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/event_analysis/')

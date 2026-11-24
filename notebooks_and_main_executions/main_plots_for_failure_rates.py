import ephys_analysis.funcs_failure_rate as funcs_fail
import pandas as pd
work_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/'
VM_data_dir = f'{work_dir}results/human/connectivity_paper/data/'
df_repatch_newest = pd.read_csv(f'{VM_data_dir}repatched_connections_looser_QC_Nov_2026.csv')
df_analysis = df_repatch_newest[df_repatch_newest.fn == '21d10035.abf']
df_analysis.fn = '21d10036.abf'

funcs_fail.interactive_failure_rates(df_analysis.reset_index(drop=True), VM_data_dir)
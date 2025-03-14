import glob
import pandas as pd
import ephys_analysis.funcs_for_results_tables as get_results
import ephys_analysis.funcs_plot_intrinsic_props as pl_intr

DEST_DIR = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/Rosie_presentation_feb25/'

HUMAN_DIR = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'
exp_view = pd.read_excel(glob.glob(HUMAN_DIR + '*experiments_overview.xlsx')[0])

#collect all QC intrinsic datatables
df_intr_props = get_results.collect_intrinsic_df()
#df_intr_props = pl_intr.get_column_RMPs_from_char(df_intr_props)

# df_intr_props = df_intr_props[(df_intr_props.resting_potential > -100) & \
#                               (df_intr_props.resting_potential < 0)]

df_intr_props.columns[df_intr_props.isna().any()].tolist()

#check for unknown patient ages
df_intr_props.OP.iloc[[indx for indx in range(len(df_intr_props)) \
                       if isinstance(df_intr_props.patient_age.iloc[indx], str)]].unique()

#filter on age and hrs incubation
adult_df = pl_intr.filter_adult_hrs_incubation_data(df_intr_props, min_age = 12, \
                                                    hrs_inc = 16, max_age = 151) # QC criteria

#adult_df = adult_df[adult_df['area'] == 'temporal']
#adult_df = adult_df[adult_df['high K concentration'] == '8 mM']
#op_color_dict = pl_intr.get_op_color_dict(adult_df)
adult_df, age_color_dict = pl_intr.get_age_color_dict(adult_df)
adult_df = pl_intr.create_new_cell_IDs(adult_df)
#adult_df = pl_intr.get_precise_treatment(adult_df)

#create slice comparison and repatch dataframes
adult_df_slice_comparison = adult_df.loc[adult_df['repatch'] == 'no']
adult_df_repatch = pl_intr.get_repatch_df(adult_df)
adult_df_repatch = adult_df_repatch.sort_values(['cell_ID_new', 'treatment'])
#adult_repatch_norm = pl_intr.get_normalized_df(adult_df_repatch)

#plot save
pl_intr.plot_param_for_days_slice(adult_df_slice_comparison, \
                                  DEST_DIR + 'slice_age/')
pl_intr.plot_param_for_days_slice_TTX_incl(adult_df_slice_comparison, age_color_dict,\
                                            DEST_DIR + 'slice_age/')
#pl_intr.plot_param_for_days_repatch(adult_df_repatch, age_color_dict, DEST_DIR + 'repatch_age/')
pl_intr.plot_param_for_days_repatch_plus_TTX(adult_df_repatch, \
                                             age_color_dict, DEST_DIR + 'repatch_age/')
pl_intr.plot_param_for_days_repatch_all_params(adult_df_repatch, age_color_dict, DEST_DIR + 'age')
#pl_intr.plot_param_for_days_repatch_norm(adult_repatch_norm, age_color_dict)

#same but color on age not OP
#pl_intr.plot_param_for_days_repatch_no_15mM(adult_df_repatch, age_color_dict, \
# DEST_DIR + 'repatch_age/')

###########################################################################################
###### Plot 3 h incubation data

adult_df_short = pl_intr.filter_adult_hrs_incubation_data(df_intr_props, min_age = 0, \
                    hrs_inc = 0, max_age = 151, max_hrs_incubation = 5) # QC criteria
#adult_df_short = adult_df_short[adult_df_short['area'] == 'temporal']
adult_df_short, age_color_dict_short = pl_intr.get_age_color_dict(adult_df_short)
adult_df_short = pl_intr.create_new_cell_IDs(adult_df_short)
#adult_df = pl_intr.get_precise_treatment(adult_df)

#create slice comparison and repatch dataframes
adult_df_slice_comparison_short = adult_df_short.loc[adult_df_short['repatch'] == 'no']
adult_df_repatch_short = pl_intr.get_repatch_df(adult_df_short)
adult_df_repatch_short = adult_df_repatch_short.sort_values(['cell_ID_new', 'treatment'])
#adult_repatch_norm_short = pl_intr.get_normalized_df(adult_df_repatch_short)

#plot all parameters
pl_intr.plot_param_for_days_slice(adult_df_slice_comparison_short, age_color_dict_short, \
                                  DEST_DIR + 'short_incubation/slice_age/')
pl_intr.plot_param_for_days_repatch_plus_TTX(adult_df_repatch_short, age_color_dict_short, \
                                             DEST_DIR + 'short_incubation/repatch_age/')
pl_intr.plot_param_for_days_repatch_all_params(adult_df_repatch_short, age_color_dict_short, \
                                               DEST_DIR + 'short_incubation/')
#pl_intr.plot_param_for_days_repatch_norm(adult_repatch_norm_short, age_color_dict_short, \
#                                               DEST_DIR + 'short_incubation/3h_inc/')

#for connectivity summary plot

# all_cons_IC, all_cons_VC =  get_results.collect_connections_df()
all_cons_IC = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results' + \
                 '/human/data/summary_data_tables/connectivity/28_Nov_20224_all_cons_IC.xlsx')

op_color_dict = pl_intr.get_op_color_dict(all_cons_IC)
# exclude_not_adult = pl_intr.in_list1_not_in_list2(all_cons_IC['OP'].unique(), \
#    adult_df['OP'].unique())

all_cons_IC['Amp 1'] = all_cons_IC['Amp 1'].fillna(0)
all_cons_IC['Amp 2'] = all_cons_IC['Amp 2'].fillna(0)
all_cons_IC['Amp 3'] = all_cons_IC['Amp 3'].fillna(0)
all_cons_IC['Amp 4'] = all_cons_IC['Amp 4'].fillna(0)
all_cons_IC.columns[all_cons_IC.isna().any()].tolist()

# for OP in exclude_not_adult:
#     all_cons_IC = all_cons_IC.drop(all_cons_IC.index[all_cons_IC['OP'] == OP])
#     all_cons_VC = all_cons_IC.drop(all_cons_IC.index[all_cons_IC['OP'] == OP])
# all_cons_IC.reset_index(inplace = True, drop = True)
# all_cons_VC.reset_index(inplace = True, drop = True)

connect_QC_passed = pl_intr.get_QC_connectivity_df(all_cons_IC)
connect_slice = connect_QC_passed[connect_QC_passed['repatch'] == 'no']
pl_intr.plot_connect_amplitude(df = connect_slice, data_type = 'all_IC', con_ID_col='connection_ID', \
    op_color_dict = op_color_dict , results_ ='amp', destination_dir =DEST_DIR)

#repatch data
repatch_connect = pl_intr.get_repatch_connectivity_df(connect_QC_passed)
# how to deal with disappearing connections
repatch_connect.fillna(0, inplace = True) #fill missing values with 0
pl_intr.plot_connect_amplitude(repatch_connect, 'repatch_IC', \
    op_color_dict , 'connection_ID', 'amp', DEST_DIR)

#connect_QC_passed_VC = pl_intr.get_QC_connectivity_VC(all_cons_VC)
# pl_intr.plot_connect_amplitude(connect_QC_passed_VC, 'all_VC', op_color_dict,  \
#       '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results' + \
#       '/human/plots/for_lab_seminar_13.03.24/connectivity_plots/')

# repatch_connect_VC = pl_intr.get_repatch_connectivity_df(connect_QC_passed_VC)
# pl_intr.plot_connect_amplitude(repatch_connect_VC, 'repatch_VC')


# Plot Initial firing frequency and number of APs against injected current

IFF_collected = get_results.collect_IFF_dfs(HUMAN_DIR = \
    '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/')

IFF_QC, IFF_repatch = pl_intr.get_QC_for_IFF_df(adult_df, IFF_collected)

IFF_QC = pl_intr.remove_non_firing_cells_D1(IFF_QC)
IFF_QC = IFF_QC.loc[IFF_QC['treatment'] != 'TTX']
IFF_QC_slice = IFF_QC.loc[IFF_QC['repatch'] == 'no']
IFF_QC_slice_firing = pl_intr.remove_non_firing_cells_D1(IFF_QC_slice)
IFF_QC_slice_firing = pl_intr.get_precise_treatment(IFF_QC_slice_firing)


pl_intr.plot_IFF_avg_against_current_inj(IFF_QC_slice_firing, 'slice_firing', 'IFF',
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human' + \
    '/plots/for_lab_seminar_13.03.24/IFF/')


IFF_repatch_firing = pl_intr.remove_non_firing_cells_D1(IFF_repatch)
IFF_repatch_firing = IFF_repatch_firing.loc[IFF_repatch_firing['treatment'] != 'TTX']
IFF_repatch_firing = pl_intr.get_precise_treatment(IFF_repatch_firing)

pl_intr.plot_IFF_avg_against_current_inj(IFF_repatch_firing, 'repatch_firing', 'IFF',
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/' + \
    'results/human/plots/for_lab_seminar_13.03.24/IFF/')

#analysis high K

df_resting = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data' + \
            '/human/data_verji/OP221116/data_tables/OP221116_RMPs_over_time.xlsx')
SAVE_DIR_RMP = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/high K evaluation/'

df_resting = pl_intr.format_RMP_column(df_resting)
df_resting = pl_intr.get_hh_mm_from_rec_time(df_resting)
pl_intr.plot_RMP_time(df_resting, SAVE_DIR_RMP)

df_ap_props = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human' + \
                            '/data_verji/OP221116/data_tables/OP221116_changes_in_AP_props.xlsx')
df_ap_props = pl_intr.get_hh_mm_from_rec_time(df_ap_props)
ap_props_dict = pl_intr.dict_for_plotting()
pl_intr.plot_ap_props(df_ap_props, ap_props_dict, SAVE_DIR_RMP)


#plotting whole RMP changes of a Ch over recordings
FILE_FOLDER = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP221116/'
PLOTS_DEST = FILE_FOLDER + '/plots/'
files = ['22n16037.abf', '22n16039.abf', '22n16040.abf']
CHAN = 6

pl_intr.plot_full_RMP_trace(FILE_FOLDER, files, PLOTS_DEST, CHAN) #files = ['fn1', 'fn2']

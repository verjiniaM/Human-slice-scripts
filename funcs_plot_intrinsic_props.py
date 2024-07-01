
import matplotlib.pyplot as plt 
import pandas as pd
import datetime
import numpy as np
import math
import funcs_human_characterisation as hcf
import funcs_for_results_tables as get_results


#import seaborn as sns

plt.style.use('./style_plot_intrinsic.mplstyle')

#mpl.rcParams['axes.prop_cycle'] = cycler(color=['r', 'g', 'b', 'y']) 3 changes the default colors

def patient_age_to_float(df_intr_props, min_age, max_age=151):
    '''
    takes df intrinsic properties and min and max age
    for adult max_age has to be 151  
    filters df based on patient_age
    '''

    #change unknown adult age to 150 and unknown 'J' to -1
    mask_adult = df_intr_props['patient_age'] == 'A'
    df_intr_props.loc[mask_adult, 'patient_age'] = 150
    mask_juv = df_intr_props['patient_age'] == 'J'
    df_intr_props.loc[mask_juv, 'patient_age'] = 0

    #change the datatype to 
    df_intr_props['patient_age'].astype(float) 

    age_filtered_df = df_intr_props[(df_intr_props['patient_age'] >= min_age) &
     (df_intr_props['patient_age'] < max_age)]

    #return back 'A' or 'J'
    age_filtered_df = age_filtered_df.replace({'patient_age' : 0}, 'J')
    age_filtered_df = age_filtered_df.replace({'patient_age' : 150}, 'A')

    return age_filtered_df


def change_to_numeric(df):
    '''
    takes df intrinsic properties  
    changes datatype to infered ones (usully accurate)
    '''
    df = df.infer_objects()
    return df

def get_QC_data(df):
    mask = (df['Rs'] < 40) & \
        (df['resting_potential'] < -45 ) & (df['resting_potential'] > -90) & \
                (df['membra_time_constant_tau'] > -25 ) &\
                    (df['capacitance'] < 900) &  (df['capacitance'] > 10)
    #(df['max_repol'] > -100) 
    #(df['AP_heigth'] > 30) & \
    # (df['Rheobase'] > 0) & \                  

    #include filter for TH ?
    df = df.loc[mask, :]
    return df

def filter_on_hrs_incubation(df, min_inc_time, max_hrs_incubation):
    mask = (df['hrs_incubation'] == 0) | (df['hrs_incubation'] >= min_inc_time) 
    df = df.loc[mask, :]
    if max_hrs_incubation < 50:
        mask = (df['hrs_incubation'] <= max_hrs_incubation) 
        df = df.loc[mask, :]
    return df

def filter_adult_hrs_incubation_data(df_intr_props, min_age, hrs_inc, max_age = 151, max_hrs_incubation = 50):
    adult_df = patient_age_to_float(df_intr_props, min_age, max_age)
    adult_df = change_to_numeric(adult_df)
    adult_df = get_QC_data(adult_df)
    adult_df = filter_on_hrs_incubation(adult_df, hrs_inc, max_hrs_incubation)
    return adult_df

def create_new_cell_IDs(df):

    if 'cell_ID_new' in df.columns:
        return df

    patcher_dict = {'Verji':'vm', 'Rosie': 'rs'}
    cell_IDs_new = []
    for i in range(len(df)):
        cell_ID = patcher_dict[df['patcher'].tolist()[i]] + df['cell_ID'].tolist()[i]
        cell_IDs_new.append(cell_ID)
    df['cell_ID_new'] = cell_IDs_new
    return df

def get_repatch_df(df):
    '''
    creates new unique cell IDs for each patcher (with initials)
    removes cells which are reaptched but one day was removed (not up to quality standards)
    '''
    #take only repatched cells
    repatch_df = df[df['repatch'] == 'yes']
    repatch_df.reset_index(inplace = True, drop = True)

    #create unique cell ID
    repatch_df = create_new_cell_IDs(repatch_df)

    #check if cell IDs appear twice, if not remove
    not_repatched_cells = []
    for cell in repatch_df['cell_ID_new'].unique():
        if len(repatch_df[repatch_df['cell_ID_new'] == cell]) < 2:
            not_repatched_cells.append(cell)
    
    for cell in not_repatched_cells:
        repatch_df = repatch_df.drop(repatch_df.index[repatch_df['cell_ID_new'] == cell])     
    repatch_df.reset_index(inplace = True, drop = True)

    return repatch_df

def get_normalized_df(df_repatch):
    '''
    adds columns with normalized data from Rs:capacitance 
    '''
    df_repatch = df_repatch.infer_objects()

    sorted_D1 = df_repatch.loc[df_repatch['day'] == 'D1'].sort_values(by=['cell_ID_new'])
    sorted_D2 = df_repatch.loc[df_repatch['day'] == 'D2'].sort_values(by=['cell_ID_new'])
    df_sorted = pd.concat([sorted_D1.loc[:], sorted_D2]).reset_index(drop=True)

    cell_diffs = [print(cell) for i, cell in enumerate(sorted_D1['cell_ID_new'].tolist()) if cell != sorted_D2['cell_ID_new'].tolist()[i]] 
    if cell_diffs != []:
        print('not possible to normalize, check cell_IDs')
        return

    for col in (sorted_D1.dtypes == float)['Rs':'capacitance'].index:
        d1 = ((sorted_D1[col].values - sorted_D1[col].values)/sorted_D1[col].values)*100
        d2 = ((sorted_D2[col].values - sorted_D1[col].values)/sorted_D1[col].values)*100
        df_sorted[col + '_norm'] =  np.concatenate((d1, d2), axis=None)

    return df_sorted


def in_list1_not_in_list2(list1, list2):
    missing_in_list_2 = []
    [missing_in_list_2.append(item) for item in list1 if item not in list2]
    return missing_in_list_2

def fix_OP220322_cell_IDs(intrinsic_df):
    mask_OP = intrinsic_df['OP'] == 'OP220322'
    cell_IDs_fix = []
    for cell in intrinsic_df.loc[mask_OP, 'cell_ID']:
        cell_IDs_fix.append('22322' + cell[5:])

    intrinsic_df.loc[mask_OP, 'cell_ID'] = cell_IDs_fix
    return intrinsic_df

def find_repeating_cell_IDs(df):
    df = create_new_cell_IDs(df)
    repeating_cells, repeats = [], []
    for cell in df['cell_ID_new'].unique():
        cells = df['cell_ID_new'][df['cell_ID_new'] == cell]
        if len(cells) > 2:
            repeating_cells.append(cells)
            repeats.append(len(cells))
    return repeating_cells, repeats


def get_precise_treatment(df):
    high_K_concentr = df['high K concentration'].tolist()
    tr_precise = []
    for i, tr in enumerate(df['treatment']):
        if tr == 'Ctrl' or tr == 'TTX':
            tr_precise.append(tr)
        else:
            tr_precise.append(high_K_concentr[i])
    df['treatment'] = tr_precise
    return df


def dict_for_plotting():
    '''
    key is the param as in the data table column
    values[0] plot title
    values[1] y asix
    values[2] plot y ticks
    '''
    titles_dict = {'Rs': ['Series resistance', 'MΩ', [5,10,15,20,25,30]], 
    'AP_halfwidth': ['AP halfwidth', 'ms', [0.5, 1, 1.5, 2, 2.5]],
    'Rin': ['Input resistance', 'MΩ', [0,100,200,300,400]],
    'resting_potential': ['Resting membrane potential', 'mV',[-80, -70,-60,-50]],
    'max_spikes': ['Maximum number of spikes', 'count', [0,10,20,30,40,50]],
    'Rheobase': ['Rheobase', 'pA', [100,300,500,700,900]],
    'AP_heigth': ['AP amplitude', 'mV', [50,60,70,80,90,100]],    
    'TH': ['AP threshold', 'mV', [-50, -40, -30, -20, -10 ]],
    'max_depol': ['AP upstroke (30% to 70%)', 'mV/ms', [100,200,300,400,500]],
    'max_repol': ['AP downstroke (70% to 30%)', 'mV/ms', [-100, -80, -60, -40, -20]],
    'membra_time_constant_tau': ['Membrane time constant', 'ms', [10, 20, 30, 40, 50]],
    'capacitance': ['Capacitance', 'pF', [100, 300, 500, 700]]}
    return titles_dict

def dict_for_plotting_synaptic():
    '''
    key is the param as in the data table column
    values[0] plot title
    values[1] y asix
    values[2] plot y ticks
    '''
    titles_dict = {'amplitude mean': ['Mean amplitude', 'pA', [-50, -40, -30, -20, -10, 0]], 
    'frequency': ['Event frequency', 'Hz', [0,2,4,6,8,10]],
    'risetime mean': ['Mean risetime (10 to 90)', 'ms', [0.0002, 0.002, 0.01]],
    'decaytime mean': ['Mean half decay time', 'ms', [0.003, 0.004, 0.005, 0.006]]}
    return titles_dict

def dict_for_plotting_con_screen():
    titles_dict = {'Amp.1': ['EPSP amplitude', 'mV', [0,1,2,3,4,5]]}
    return titles_dict

def plot_time_after_OP_vs_param(df, param,
destintaion_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/intrinsic_properties/hrs_after_op_dependency/'): 
    titles_dict = dict_for_plotting()

    fig1 = plt.figure(figsize=(7,7))
    ax = plt.subplot(1, 1, 1)
    plt.subplots_adjust(hspace=0.5)
    plt.scatter(df['hrs_after_OP'], df[param])

    plt.title('AP ' + titles_dict[param][0] + ' time dependency', fontsize = 19)
    ax.set_xlabel('Time after cortex resection (hrs)', fontsize = 15)


    date = str(datetime.date.today())
    
    fig1.patch.set_facecolor('white')
    plt.savefig(destintaion_dir + date + param +'_vs_time_after_op.png')
    plt.close(fig1)

def get_op_color_dict(df):
    '''
    input should be adult df with all OPs
    assigns a color to each OP in the list
    '''
    cmap = plt.cm.get_cmap('tab20')
    op_color_dict = {}
    OP_list = sorted(list(df['OP'].unique()))
    for h, op in enumerate(OP_list):
        op_color_dict[op] = cmap((h+1)/len(OP_list))
    return op_color_dict


def get_age_color_dict(df):
    '''
    input should be adult df with all OPs
    assigns a color to each OP in the list
    '''
    cmap = plt.cm.get_cmap('copper')
    num_patient_age = []
    for i in range(len(df)):
        if df['patient_age'].tolist()[i] == 'A':
            num_patient_age.append(151)
        elif df['patient_age'].tolist()[i] == 'J':
            num_patient_age.append(0)
        else:
            num_patient_age.append(df['patient_age'].tolist()[i])
    df.insert(len(df.columns), 'patient_age_num', num_patient_age)

    age_color_dict = {}
    age_list = sorted(list(df.patient_age_num.unique()))
    for h, years in enumerate(age_list):
       age_color_dict[years] = cmap((h+1)/len(age_list))
    return df, age_color_dict

def plot_param_for_days_slice(df, op_color_dict,
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/intrinsic_properties/adult_slice/'):

    treatments = ['Ctrl', 'high K']
    #treatments = ['Ctrl', 'TTX', 'high K']
    titles_dict = dict_for_plotting()

    colors = ['#dede00', '#ff7f00', '#dede00', '#4daf4a','#dede00','#984ea3'] #Ctrl, 8mM Kcl, 15 mM KCl

    for u, param in enumerate(titles_dict.keys()): 
        fig2 = plt.figure(figsize=(7,7))
        ax = plt.subplot(1,1,1)
        day_label = []
        num_cels = {}
        data_boxplot = []
        for i, tr in enumerate(treatments):
            for j, day in enumerate(df['day'].unique()):
                k = j + 2*i 
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                median = df_plot[param].median()
                x = np.linspace(0.65+k, 1.35+k, len(df_plot))
                if tr == 'high K' and day == 'D2':
                    #ax.scatter(x, df_plot[param], alpha = 0.7, s = 40, c = df_plot['high K concentration'], cmap = colormap_K)
                    y_8mm = df_plot[param].loc[df_plot['high K concentration'] == '8 mM'] 
                    x_8 = np.linspace(0.8+k, 1.05+k, len(y_8mm))
                    median_8 = np.median(y_8mm)
                    y_15mm = df_plot[param].loc[df_plot['high K concentration'] == '15 mM'] 
                    x_15 = np.linspace(1.1+k, 1.35+k, len(y_15mm))
                    median_15 = np.median(y_15mm)
                    ax.scatter(x_8, y_8mm, c = colors[3], s = 80, alpha = 0.75)
                    ax.scatter(x_15, y_15mm, c = colors[5], s = 80, alpha = 0.75)
                    ax.boxplot(y_8mm, positions = [k+0.6], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[3], alpha = 0.75),
                    medianprops = dict(linewidth=2.3, color = 'k'))
                    ax.boxplot(y_15mm, positions = [k+1.55], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[5], alpha = 0.75),
                    medianprops = dict(linewidth=2.3, color = 'k'))
                else:
                    ax.boxplot(df_plot[param], positions = [k + 0.5], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[k], alpha = 0.75),
                    medianprops = dict(linewidth=2.3, color = 'k'))    
                    ax.scatter(x, df_plot[param], c = colors[int(k)], s = 80, alpha = 0.75)
                
                # yerr = 1.253*(df_plot[param].std()/(math.sqrt(len(df_plot))))
                # ax.scatter(1+k, median, color = 'k', marker = '_', s = 2000)
                # ax.text(0.9+k, (1.05*median), str(round(median,2)), size = 12)
                day_label.append(day)
                num_cels[tr + ' ' + day] = len(df_plot)
                data_boxplot.append(df_plot[param])
            ax.text(1 + 2*i, int(np.max(df[param])), tr, size = 17, c = colors[2*i+1])
            ax.text(1.7 + 2*i, int(np.max(df[param])), 'n = ' + str(len(df_plot)), size = 12, c = colors[2*i+1])

        #plt.boxplot(data_boxplot, showbox = False)
        ax.tick_params(axis='y', labelsize=22)
        ax.set_xticks(ticks = [0.6,1.6,2.6,3.55,4.6], labels = ['D1', 
        'D2 \n Ctrl', 'D1', 'D2 \n 8mM','D2 \n 15mM'], size = 20)
        #when plotting TTX as well
        # ax.set_xticks(ticks = [0.6,1.6,2.6,3.55,4.6, 5.6,6.6], labels = ['D1', 
        # 'D2 \n Ctrl','D1', 'D2 \n TTX' ,'D1', 'D2 \n 8mM','D2 \n 15mM'], size = 20)

        plt.title(titles_dict[param][0], fontsize = 19, x = 0.5, y = 1)
        ax.set_ylabel(titles_dict[param][1], fontsize = 24)
        ax.set_yticks(ticks = titles_dict[param][2])

        plt.subplots_adjust(hspace=0.35)
        fig2.patch.set_facecolor('white')
        fig2.tight_layout()

        date = str(datetime.date.today())
        plt.savefig(destination_dir  + date + '_plot_' + param + '.pdf')
        plt.close(fig2)

def plot_param_for_days_repatch(df, op_color_dict,
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/intrinsic_properties/adult_repatch/'):
    
    treatments = ['Ctrl', 'high K']
    #treatments = ['Ctrl', 'TTX', 'high K']
    titles_dict = dict_for_plotting()
        
    colors = ['#dede00', '#ff7f00', '#dede00', '#4daf4a','#dede00','#984ea3'] #Ctrl, 8mM Kcl, 15 mM Kcl
    for u, param in enumerate(titles_dict.keys()): 
        fig2 = plt.figure(figsize=(7,7))
        ax = plt.subplot(1,1,1)
        day_label = []
        num_cels = {}
        for i, tr in enumerate(treatments):
            x_plot =[]
            
            for j, day in enumerate(df['day'].unique()):
                k = j + 2*i 
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                df_plot.reset_index(drop=True, inplace=True)
                median = df_plot[param].median()
                x = np.linspace(0.7+k, 1.3+k, len(df_plot))
                x_plot.append(x)

                if tr == 'high K' and day == 'D2':
                    y_8mm = df_plot.loc[df_plot['high K concentration'] == '8 mM']
                    x_8 = np.linspace(0.8+k, 1.05+k, len(y_8mm))
                    #x_8 = x[:len(y_8mm)]
                    median_8 = np.median(y_8mm[param])
                    y_15mm = df_plot.loc[df_plot['high K concentration'] == '15 mM'] 
                    x_15 = np.linspace(1.4+k, 1.65+k, len(y_15mm))
                    #x_15 = x[-len(y_15mm):]
                    median_15 = np.median(y_15mm[param])
                    ax.scatter(x_8, y_8mm[param], c = colors[3], s = 80, zorder = 2)
                    ax.scatter(x_15, y_15mm[param], c = colors[5], s = 80, zorder = 2)
                    ax.plot([0.7+k, 1+k], [median_8, median_8], c = 'k', linestyle = 'solid')
                    ax.plot([1.3+k, 1.6+k], [median_15, median_15], c = 'k', linestyle = 'solid')
                    ax.text(k+0.65, median_8 + 0.5, str(round(median_8, 2)), size = 15)
                    ax.text(k+1.25, median_15 + 0.5, str(round(median_15, 2)), size = 15)
                    df_plot.insert(0, 'x', np.concatenate([x_15, x_8]))

                    for c, cell in enumerate(df_plot['cell_ID_new']):
                        #indx = index_s[c]
                        x_K = df_plot['x'].loc[df_plot['cell_ID_new'] == cell].tolist()[0]
                        x1 = [x_plot[0][c], x_K]
                        y = df[param][df['cell_ID_new'] == cell]
                        op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
                        age = df_plot['patient_age_num'][df_plot['cell_ID_new'] == cell].tolist()[0]
                        if type(list(op_color_dict.keys())[0]) is not str :
                            plt.plot(x1, y, '-', color = op_color_dict[age], alpha = 0.5, linewidth = 2, zorder = 1)
                        else:
                            plt.plot(x1, y, '-', color = op_color_dict[op], alpha = 0.5, linewidth = 2, zorder = 1)


                else:
                    ax.scatter(x, df_plot[param], c = colors[int(k)], s = 80, zorder = 2)
                    ax.scatter(1+k, median, color = 'k', marker = '_', s = 2000, zorder = 2)
                    ax.text(0.75+k, (1.05*median), str(round(median,2)), size = 15)
                    x_plot.append(x)
                    if k in [1,3,5]:
                        for c, cell in enumerate(df_plot['cell_ID_new']):
                            #indx = index_s[c]
                            x1 = [x_plot[0][c], x[c]] 
                            y = df[param][df['cell_ID_new'] == cell]
                            op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
                            age = df_plot['patient_age_num'][df_plot['cell_ID_new'] == cell].tolist()[0]
                            if type(list(op_color_dict.keys())[0]) is not str :
                                plt.plot(x1, y, '-', color = op_color_dict[age], alpha = 0.4, linewidth = 2, zorder = 1)
                            else:
                                plt.plot(x1, y, '-', color = op_color_dict[op], alpha = 0.4, linewidth = 2, zorder = 1)

                day_label.append(day)
                num_cels[tr + ' ' + day] = len(df_plot)

                # ax.text(1 + 2*i, int(np.max(df[param])), tr, size = 17, c = colors[2*i+1])
                # ax.text(1.7 + 2*i, int(np.max(df[param])), 'n = ' + str(len(y_8mm)), size = 12, c = colors[2*i+1])
    
        ax.tick_params(axis='y', labelsize=22) 
        ax.set_xticks(ticks = [1,2,3,3.8,4.55], labels = ['D1', 
        'D2 \n Ctrl', 'D1', 'D2 \n 8mM','D2 \n 15mM'], size = 20)
        #when plotting TTX as well
        # ax.set_xticks(ticks = [0.6,1.6,2.6,3.55,4.6, 5.6,6.6], labels = ['D1', 
        # 'D2 \n Ctrl','D1', 'D2 \n TTX' ,'D1', 'D2 \n 8mM','D2 \n 15mM'], size = 20)
        ax.set_ylabel(titles_dict[param][1], fontsize = 24)
        ax.set_yticks(ticks = titles_dict[param][2])
        #plt.title(titles_dict[param][0], fontsize = 19, x = 0.5, y = 1)

        fig2.patch.set_facecolor('white')
        fig2.tight_layout()

        date = str(datetime.date.today())
        plt.savefig(destination_dir + date + '_plot_' + param + '.pdf')
        plt.close(fig2)
    

def plot_param_for_days_repatch_no_15mM(df, op_color_dict,
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/intrinsic_properties/adult_repatch/'):
    
    treatments = ['Ctrl', 'high K']
    #treatments = ['Ctrl', 'TTX', 'high K']
    titles_dict = dict_for_plotting()
        
    colors = ['#dede00', '#ff7f00', '#dede00', '#4daf4a','#dede00','#984ea3'] #Ctrl, 8mM Kcl, 15 mM Kcl
    for u, param in enumerate(titles_dict.keys()): 
        fig2 = plt.figure(figsize=(7,7))
        ax = plt.subplot(1,1,1)
        day_label = []
        num_cels = {}
        for i, tr in enumerate(treatments):
            x_plot =[]
            
            for j, day in enumerate(df['day'].unique()):
                k = j + 2*i 
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                df_plot.reset_index(drop=True, inplace=True)
                median = df_plot[param].median()
                x = np.linspace(0.7+k, 1.3+k, len(df_plot))
                x_plot.append(x)

                if tr == 'high K' and day == 'D2':
                    y_8mm = df_plot.loc[df_plot['high K concentration'] == '8 mM']
                    x_8 = np.linspace(0.8+k, 1.3+k, len(y_8mm))
                    median_8 = np.median(y_8mm[param])
                    ax.scatter(x_8, y_8mm[param], c = colors[3], s = 80, zorder = 2)
                    ax.plot([0.7+k, 1.3+k], [median_8, median_8], c = 'k', linestyle = 'solid')
                    ax.text(k+0.75, median_8 + 0.5, str(round(median_8, 2)), size = 15)
                    df_plot.insert(0, 'x', x_8)

                    for c, cell in enumerate(df_plot['cell_ID_new']):
                        x_K = df_plot['x'].loc[df_plot['cell_ID_new'] == cell].tolist()[0]
                        x1 = [x_plot[0][c], x_K]
                        y = df[param][df['cell_ID_new'] == cell]
                        op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
                        age = df_plot['patient_age_num'][df_plot['cell_ID_new'] == cell].tolist()[0]
                        if type(list(op_color_dict.keys())[0]) is not str :
                            plt.plot(x1, y, '-', color = op_color_dict[age], alpha = 0.5, linewidth = 2, zorder = 1)
                        else:
                            plt.plot(x1, y, '-', color = op_color_dict[op], alpha = 0.5, linewidth = 2, zorder = 1)


                else:
                    ax.scatter(x, df_plot[param], c = colors[int(k)], s = 80, zorder = 2)
                    ax.scatter(1+k, median, color = 'k', marker = '_', s = 2000, zorder = 2)
                    ax.text(0.75+k, (1.05*median), str(round(median,2)), size = 15)
                    x_plot.append(x)
                    if k in [1,3,5]:
                        for c, cell in enumerate(df_plot['cell_ID_new']):
                            #indx = index_s[c]
                            x1 = [x_plot[0][c], x[c]] 
                            y = df[param][df['cell_ID_new'] == cell]
                            op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
                            age = df_plot['patient_age_num'][df_plot['cell_ID_new'] == cell].tolist()[0]
                            if type(list(op_color_dict.keys())[0]) is not str :
                                plt.plot(x1, y, '-', color = op_color_dict[age], alpha = 0.4, linewidth = 2, zorder = 1)
                            else:
                                plt.plot(x1, y, '-', color = op_color_dict[op], alpha = 0.4, linewidth = 2, zorder = 1)

                day_label.append(day)
                num_cels[tr + ' ' + day] = len(df_plot)

                # ax.text(1 + 2*i, int(np.max(df[param])), tr, size = 17, c = colors[2*i+1])
                # ax.text(1.7 + 2*i, int(np.max(df[param])), 'n = ' + str(len(y_8mm)), size = 12, c = colors[2*i+1])
    
        ax.tick_params(axis='y', labelsize=22) 
        ax.set_xticks(ticks = [1,2,3,4], labels = ['D1', 
        'D2 \n Ctrl', 'D1', 'D2 \n 8mM'], size = 20)
        ax.set_ylabel(titles_dict[param][1], fontsize = 24)
        ax.set_yticks(ticks = titles_dict[param][2])

        fig2.patch.set_facecolor('white')
        fig2.tight_layout()

        date = str(datetime.date.today())
        plt.savefig(destination_dir + date + '_plot_' + param + '.pdf')
        plt.close(fig2)


def plot_param_for_days_repatch_all_params(df, op_color_dict,
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/intrinsic_properties/adult_repatch/'):
    
    treatments = ['Ctrl', 'high K']
    #treatments = ['Ctrl', 'TTX', 'high K']
    titles_dict = dict_for_plotting()
        
    colors = ['#dede00', '#ff7f00', '#dede00', '#4daf4a','#dede00','#984ea3'] #Ctrl, 8mM Kcl, 15 mM Kcl
    nrows, ncols = 3, 4
    fig2, axs = plt.subplots(nrows, ncols,figsize=(22,30))

    for u, param in enumerate(titles_dict.keys()): 
        ax = axs[u // ncols, u % ncols]
        day_label = []
        num_cels = {}
        for i, tr in enumerate(treatments):
            x_plot =[]

            for j, day in enumerate(df['day'].unique()):
                k = j + 2*i 
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                df_plot.reset_index(drop=True, inplace=True)
                median = df_plot[param].median()
                x = np.linspace(0.7+k, 1.3+k, len(df_plot))
                x_plot.append(x)

                if tr == 'high K' and day == 'D2':
                    y_8mm = df_plot.loc[df_plot['high K concentration'] == '8 mM']
                    x_8 = np.linspace(0.8+k, 1.05+k, len(y_8mm))
                    #x_8 = x[:len(y_8mm)]
                    median_8 = np.median(y_8mm[param])
                    y_15mm = df_plot.loc[df_plot['high K concentration'] == '15 mM'] 
                    x_15 = np.linspace(1.4+k, 1.65+k, len(y_15mm))
                    #x_15 = x[-len(y_15mm):]
                    median_15 = np.median(y_15mm[param])
                    ax.scatter(x_8, y_8mm[param], c = colors[3], s = 80, zorder = 2)
                    ax.scatter(x_15, y_15mm[param], c = colors[5], s = 80, zorder = 2)
                    ax.plot([0.7+k, 1+k], [median_8, median_8], c = 'k', linestyle = 'solid')
                    ax.plot([1.3+k, 1.6+k], [median_15, median_15], c = 'k', linestyle = 'solid')
                    ax.text(k+0.65, median_8 + 0.5, str(round(median_8, 2)), size = 15)
                    ax.text(k+1.25, median_15 + 0.5, str(round(median_15, 2)), size = 15)
                    df_plot.insert(0, 'x', np.concatenate([x_15, x_8]))

                    for c, cell in enumerate(df_plot['cell_ID_new']):
                        #indx = index_s[c]
                        x_K = df_plot['x'].loc[df_plot['cell_ID_new'] == cell].tolist()[0]
                        x1 = [x_plot[0][c], x_K]
                        y = df[param][df['cell_ID_new'] == cell]
                        op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]                        
                        age = df_plot['patient_age_num'][df_plot['cell_ID_new'] == cell].tolist()[0]
                        if type(list(op_color_dict.keys())[0]) is not str:
                            ax.plot(x1, y, '-', color = op_color_dict[age], alpha = 0.5, linewidth = 2, zorder = 1)
                        else:
                            ax.plot(x1, y, '-', color = op_color_dict[op], alpha = 0.5, linewidth = 2, zorder = 1)


                else:
                    ax.scatter(x, df_plot[param], c = colors[int(k)], s = 80, zorder = 2)
                    ax.scatter(1+k, median, color = 'k', marker = '_', s = 2000, zorder = 2)
                    ax.text(0.75+k, (1.05*median), str(round(median,2)), size = 15)
                    x_plot.append(x)
                    if k in [1,3,5]:
                        for c, cell in enumerate(df_plot['cell_ID_new']):
                            #indx = index_s[c]
                            x1 = [x_plot[0][c], x[c]] 
                            y = df[param][df['cell_ID_new'] == cell]
                            op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
                            age = df_plot['patient_age_num'][df_plot['cell_ID_new'] == cell].tolist()[0]
                            if type(list(op_color_dict.keys())[0]) is not str :
                                ax.plot(x1, y, '-', color = op_color_dict[age], alpha = 0.4, linewidth = 2, zorder = 1)
                            else:
                                ax.plot(x1, y, '-', color = op_color_dict[op], alpha = 0.4, linewidth = 2, zorder = 1)


                day_label.append(day)
                num_cels[tr + ' ' + day] = len(df_plot)

            # ax.text(1 + 2*i, int(np.max(df[param])), tr, size = 17, c = colors[2*i+1])
            # ax.text(1.7 + 2*i, int(np.max(df[param])), 'n = ' + str(len(y_8mm)), size = 12, c = colors[2*i+1])

        ax.tick_params(axis='y', labelsize=22) 
        ax.set_xticks(ticks = [1,2,3,3.8,4.55], labels = ['D1', 
        'D2 \n Ctrl', 'D1', 'D2 \n 8mM','D2 \n 15mM'], size = 20)
        #when plotting TTX as well
        # ax.set_xticks(ticks = [0.6,1.6,2.6,3.55,4.6, 5.6,6.6], labels = ['D1', 
        # 'D2 \n Ctrl','D1', 'D2 \n TTX' ,'D1', 'D2 \n 8mM','D2 \n 15mM'], size = 20)
        ax.set_ylabel(titles_dict[param][1], fontsize = 24)
        ax.set_yticks(ticks = titles_dict[param][2])
        ax.set_title(titles_dict[param][0], fontsize = 14)

    fig2.patch.set_facecolor('white')
    fig2.tight_layout()

    date = str(datetime.date.today())
    plt.savefig(destination_dir + date + '_plot_all_params' + '.pdf')
    plt.close(fig2)


def plot_param_for_days_repatch_norm(df, op_color_dict,
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/intrinsic_properties/adult_repatch/'):
    
    treatments = ['Ctrl', 'high K']
    #treatments = ['Ctrl', 'TTX', 'high K']
    titles_dict = dict_for_plotting()   

    titles_dict_norm = {}
    for (key, value) in titles_dict.items():
        titles_dict_norm[key + '_norm'] = value
    titles_dict = titles_dict_norm

    colors = ['moccasin', '#ff7f00', 'moccasin', '#4daf4a','moccasin','#984ea3'] #Ctrl, 8mM Kcl, 15 mM Kcl

    for u, param in enumerate(titles_dict.keys()): 
        fig2 = plt.figure(figsize=(4,7))
        ax = plt.subplot(1,1,1)
        day_label = []
        num_cels = {}
        data_boxplot = []
        for i, tr in enumerate(treatments):
            x_plot =[]
            for j, day in enumerate(df['day'].unique()):
                k = j + 2*i 
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                x = np.repeat((k + 0.5), len(df_plot))
                x_plot.append(x)

            if k in [1,3,5]:
                for c, cell in enumerate(df_plot['cell_ID_new']):
                    x1 = [x_plot[0][c], x[c]] 
                    y = df[param][df['cell_ID_new'] == cell]
                    op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
                    age = df_plot['patient_age_num'][df_plot['cell_ID_new'] == cell].tolist()[0]
                    if type(list(op_color_dict.keys())[0]) is not str :
                        plt.plot(x1, y, '-', color = op_color_dict[age], alpha = 0.4, linewidth = 2, zorder = 1)
                    else:
                        plt.plot(x1, y, '-', color = op_color_dict[op], alpha = 0.4, linewidth = 2, zorder = 1)

        ax.tick_params(axis='y', labelsize=22)
        ax.set_xticks(ticks = [0.5,1.5,2.5,3.5], labels = ['D1', 
        'D2 \n Ctrl', 'D1', 'D2 \n high K'], size = 20)
        #ax.set_xticks(ticks = [0.5,1.5,2.5,3.5], labels = ['D1', 
        #'D2 \n Ctrl','D1', 'D2 \n TTX', 'D1', 'D2 \n high K'], size = 20)
        ax.set_ylabel('norm. % difference \n (D1-D2)/D1', fontsize = 24)

        plt.subplots_adjust(hspace=0.35)
        fig2.patch.set_facecolor('white')
        fig2.tight_layout()

        date = str(datetime.date.today())
        plt.savefig(destination_dir + date + param + '_plot.pdf')
        plt.close(fig2)


def plot_param_for_hrs_incubation(df, data_type, op_color_dict,
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/intrinsic_properties/time_dependencies/'):

    treatments = ['Ctrl', 'TTX', 'high K']
    titles_dict = dict_for_plotting()

    #destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/intrinsic_properties/'
    colors = [ 'red', 'cadetblue','mediumpurple']
    cmap = plt.cm.get_cmap('tab20')
    # op_color_dict = {}
    # for h, op in enumerate(df['OP'].unique()):
    #     op_color_dict[op] = cmap((h+1)/10)

    for param in titles_dict.keys(): 
        fig2 = plt.figure(figsize=(10,5))
        ax = plt.subplot(1,1,1)
        day_label = []
        num_cels = {}
        data_boxplot = []
        for i, tr in enumerate(treatments):
            x_plot =[]
            
            df_plot1 = df[(df['treatment'] == tr)]
            df_plot = df_plot1[df_plot1['hrs_after_OP'] > 0]
            x = df_plot['hrs_after_OP']
            ax.scatter(x, df_plot[param], c = colors[int(i)], s = 40, label = tr)

        #plt.boxplot(data_boxplot, showbox = False)
        plt.title(titles_dict[param][0] + ' over time', fontsize = 19, x = 0.5, y = 1)
        ax.set_xlabel('Hours after OP', fontsize = 15)
        ax.set_ylabel(titles_dict[param][1], fontsize = 15)

        plt.subplots_adjust(hspace=0.35)
        fig2.patch.set_facecolor('white')
        fig2.tight_layout()

        plt.figlegend(loc = 'upper right',  bbox_to_anchor=(1, 1))
        date = str(datetime.date.today())
        plt.savefig(destination_dir  + date + data_type + '_plot_' + param + '.png')
        plt.close(fig2)




#%%

#Funcs for analysis of intrinsic properties of testing high K condition

def format_RMP_column(df_resting):
    resting_format = []
    for j in range(len(df_resting)): 
        RMP_format = df_resting['resting_potential'].tolist()[j][1:-1]
        resting_format.append(RMP_format)
    df_resting['RMP_formatted'] = resting_format
    return df_resting

def get_hh_mm_from_rec_time(df):
    rec_time = []
    for j in range(len(df)): 
        rec_time.append(df['recording_time'].tolist()[j][11:16])
    df['rec_time_short'] = rec_time
    return df

def plot_RMP_time(df_resting, save_dir):
    for slic in df_resting['slice'].unique():
        df_slice = df_resting[df_resting['slice'] == slic]

        fig, ax = plt.subplots(len(df_slice['cell_ch'].unique()),1)
        fig.set_figheight(11)
        fig.set_figwidth(10)

        color_dict = {
                'puff high K': ['puff high K' , 'red'],
                'wash out' : ['wash out', 'blue'],
                'wash in high K': ['wash in high K', 'k'],
                'before' : ['before', 'white']
                }

        for i, chan in enumerate(df_slice['cell_ch'].unique()):
            k = -2
            df_chan = df_slice[df_slice['cell_ch'] == chan]

            time_label,x_all, y_all = [], [], []
            for fn in df_chan['filename'].unique():
                k = k + 1
                df_fn = df_chan[df_chan['filename'] == fn]
        
                RMPs = [ss.split(',') for ss in df_fn.RMP_formatted][0]
                num_vals = len(RMPs)
    
                x = np.linspace(1 + k, 1.9 + k, num_vals)
                y = [float(c) for c in RMPs]

                for val in range(num_vals):
                    x_all.append(x[val])
                    y_all.append(y[val])
                time_label.append(str(df_fn['recording_time'])[16:22])

            for condition in df_chan['circumstance'].unique():
                    x_1 = df_chan['rec_time_short'][df_chan['circumstance'] == condition]
                    
                    ax[i].vlines(x_1, ymin = max(y_all), ymax = min(y_all) - 5, lw = 0.45, linestyle = '-',
                    color = color_dict[condition][1], label = color_dict[condition][0])
                    color_dict[condition][0] = "_nolegend_"

            ax[i].scatter(x_all,y_all)
            ax[i].set_title('Ch' + str(chan))
            ax[i].set_xticks(ticks = list(range(2,len(time_label)+2)), labels = time_label,
            rotation = 45) 
            ax[i].set_ylabel('mV')
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)

        ax[i].set_xlabel('time')
        fig.patch.set_facecolor('white')
        plt.figlegend(loc = 'upper right',  bbox_to_anchor=(0.95, 0.95))
        fig.suptitle('Resting membrane potential', fontsize = 19)
        plt.subplots_adjust(hspace = 0.6)
        plt.savefig(save_dir + 'RMP_' + slic + '_'  + '.png')
        plt.close()

def plot_ap_props(df_ap_props, ap_props_dict, save_dir):
    for param in ap_props_dict.keys():
        if param not in df_ap_props.columns :
            continue
        
        for slic in df_ap_props['slice'].unique():
            df_slice = df_ap_props[df_ap_props['slice'] == slic]

            fig, ax = plt.subplots(len(df_slice['cell_ch'].unique()),1, sharex = True)
            fig.set_figheight(11)
            fig.set_figwidth(8)

            color_dict = {
                'puff high K': ['puff high K' , 'red'],
                'wash out' : ['wash out', 'blue'],
                'wash in high K': ['wash in high K', 'k'],
                'before' : ['before', 'white']
                }

            for i, chan in enumerate(df_slice['cell_ch'].unique()):
                k = 0
                df_chan = df_slice[df_slice['cell_ch'] == chan]

                rec_time = df_chan['rec_time_short']
                param_y = df_chan[param]

                ax[i].scatter(rec_time,param_y)
                ax[i].set_title('Ch' + str(chan))
                ax[i].set_ylabel(ap_props_dict[param][1])
                ax[i].spines['top'].set_visible(False)
                ax[i].spines['right'].set_visible(False)

                for condition in df_slice['circumstance'].unique():
                    x_1 = df_chan['rec_time_short'][df_chan['circumstance'] == condition]
                    
                    ax[i].vlines(x_1, ymin = min(param_y) -2 , ymax = max(param_y) + 2, lw = 0.3, linestyle = '-',
                    color = color_dict[condition][1], label = color_dict[condition][0])
                    color_dict[condition][0] = "_nolegend_"
            ax[i].set_xlabel('time')

            fig.patch.set_facecolor('white')
            plt.figlegend(loc = 'upper right',  bbox_to_anchor=(0.95, 0.95))
            fig.suptitle(ap_props_dict[param][0], fontsize = 19)
            plt.subplots_adjust(hspace=0.35)
            plt.savefig(save_dir + 'AP_' + slic + '_' + param + '.png')
            plt.close()


def plot_spiking_when_RMP_increases(filename, plots_destination):
    '''Finds where there is a jump in the RMP >30 mV 
    plots those sweeps '''
    data_dict = hcf.load_traces(filename)
    end_fn = filename.rfind('/') + 1

    for key in data_dict.keys():
        num_swps = np.shape(data_dict[key][0])[1]
        avg = np.mean(data_dict[key][0])

        for i in range(num_swps):
            if max(data_dict[key][0][:,i]) - min(data_dict[key][0][:,i]) > 40:
                #index = np.where(data_dict[key][0][:,i] > 40)[0][0] #finding the fifrst AP

                sampl_rate, units, times = hcf.get_abf_info(filename, int(key[-1]), num_swps, np.shape(data_dict[key][0])[0])

                fig = plt.figure(figsize=(15,7))
                ax = plt.subplot(1,1,1)
            
                ax.plot(times, data_dict[key][0][:,i])
                ax.set_xlabel('time [s]')
                ax.set_ylabel(str(units)[-2:])
                ax.set_title(filename[end_fn:-4] + key)
                fig.patch.set_facecolor('white')

                plt.savefig(plots_destination + 'RMP_big_change_' + filename[end_fn:-4] + '_sweep_' + str(i)  + '.png')
                plt.close()



def plot_full_RMP_trace(file_folder, files, plots_destination, channel, index_start_wash_in):
    '''
    files can be a list of RMP files to be plotted one after the other
    the first one should be control and then RMP recorded during wash in
    index_start_wash_in - index of the filename in files where the washin start; plotted red line
    '''
    key = 'Ch' + str(channel)
    data_plot, times_plot_all, len_times_plot = [], [], []
    for fn in files:
        filename = file_folder + fn
        data_dict = hcf.load_traces(filename)

        sweep_len = np.shape(data_dict[key][0])[0] 
        num_swps = np.shape(data_dict[key][0])[1]

        for i in range(num_swps):
            data_plot = np.append(data_plot, data_dict[key][0][:,i], 0)

        sampl_rate, units, times = hcf.get_abf_info(filename, int(key[-1]), num_swps, sweep_len)

        times_plot = np.linspace(len(times_plot_all), len(times_plot_all) + sweep_len*num_swps, sweep_len*num_swps)/sampl_rate
        len_times_plot.append(len(times_plot))
        times_plot_all = np.append(times_plot_all, times_plot)

    fig = plt.figure(figsize=(40,15))
    ax = plt.subplot(1,1,1)
    ax.plot(times_plot_all, data_plot)
    ax.vlines(times_plot_all[len_times_plot[index_start_wash_in]], ymin = min(data_plot)  - 5, 
    ymax = max(data_plot) + 5, lw = 1, linestyle = '-', color = 'red', label = 'start wash in')

    end_fn = filename.rfind('/') + 1

    ax.set_xlabel('time [s]', fontsize = 20)
    ax.set_ylabel(str(units)[-2:], fontsize = 20)
    ax.set_title(filename[end_fn:-4] + key)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.figlegend(loc = 'upper right',  bbox_to_anchor=(0.95, 0.95), fontsize = 20)
    fig.patch.set_facecolor('white')
    plt.xticks(fontsize = 20) 
    plt.yticks(fontsize = 20)

    plt.savefig(plots_destination + 'all_RMP_big_change_' + filename[end_fn:-4] + '_' + key + '.png')
    plt.close()



#%%
#Funcs for IFF and num aps agains current injection

def get_QC_for_IFF_df(intrinsic_df, IFF_df):
    ''' 
    takes cellIDs and repatch info from intrinsic df
    remove data for which no cell IDs are present (data didn't pass quality check)
    create a IFF excel sheet with data and sheet1 - all IFF data, sheet2 - repatch IFF
    '''
    
    cell_IDs, repatches, not_in_intrinsic = [], [], []
    for i in range(len(IFF_df)):
        OP = IFF_df['OP'][i]
        patcher = IFF_df['patcher'][i]
        fn = IFF_df['filename'][i]
        slic = IFF_df['slice'][i]
        day = IFF_df['day'][i]
        cell_ch = IFF_df['cell_ch'][i]

        cell_ID = intrinsic_df['cell_ID_new'][(intrinsic_df['OP'] == OP) & (intrinsic_df['patcher'] == patcher) & \
        (intrinsic_df['slice'] == slic) & (intrinsic_df['day'] == day) & \
        (intrinsic_df['cell_ch'] == cell_ch) & (intrinsic_df['filename'] == fn)]

        if len(cell_ID) == 0:
            cell_ID = hcf.get_new_cell_IDs_fn(fn, slic, [cell_ch], patcher)[0]
            
            not_in_intrinsic.append(cell_ID)
            cell_IDs.append(cell_ID)
            repatches.append('no')
            continue

        cell_IDs.append(cell_ID.tolist()[0])
        repatches.append(intrinsic_df['repatch'][intrinsic_df['cell_ID_new'] == cell_ID.tolist()[0]].tolist()[0])
    

    IFF_df.insert(7, 'cell_ID_new', cell_IDs)
    IFF_df.insert(8, 'repatch', repatches)

    #remove cells which are not in the intrinsic_df (not QC passed)
    not_QC_passed_cells = in_list1_not_in_list2(IFF_df['cell_ID_new'].tolist(), intrinsic_df['cell_ID_new'].tolist())
    for cell in not_QC_passed_cells :
        IFF_df = IFF_df.drop(IFF_df.index[IFF_df['cell_ID_new'] == cell])
    IFF_df.reset_index(inplace = True, drop = True)

    IFF_repatch = get_repatch_df(IFF_df)

    date = str(datetime.date.today())
    file_name = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/data/initial_firing_freqs/' + date + '_IFF_full_.xlsx'
    with pd.ExcelWriter(file_name) as writer:  
        IFF_df.to_excel(writer, sheet_name='Initial_firing_freq_all')
        IFF_repatch.to_excel(writer, sheet_name='Initial_firing_freq_repatch')
    
    return IFF_df, IFF_repatch

def remove_non_firing_cells_D1 (repatched_IFF_df):
    '''
    use this function on the repatch dataframe
    removes cells that are not firing any APs at any current injecion on day1
    '''
    num_aps_indx, IFF_indx = get_num_aps_and_IFF_data_culumns(repatched_IFF_df)
    
    not_firing_cells, weird_cells = [], []
    for i, cell in enumerate(repatched_IFF_df['cell_ID_new']):
        df = repatched_IFF_df[(repatched_IFF_df['cell_ID_new'] == cell) & (repatched_IFF_df['day'] == 'D1')]
        if len(df) == 0:
            weird_cells.append(cell)
            continue

        df2 = df.iloc[:, num_aps_indx].reset_index(drop=True)
        list_num_aps = df2.loc[0, :].values.flatten().tolist()
        if max(list_num_aps) <= 0:
            not_firing_cells.append(cell)

    for cell in not_firing_cells:
        repatched_IFF_df = repatched_IFF_df.drop(repatched_IFF_df.index[repatched_IFF_df['cell_ID_new'] == cell])
    repatched_IFF_df.reset_index(inplace = True, drop = True)

    return repatched_IFF_df

def get_num_aps_and_IFF_data_culumns(df):
    num_aps_indx, IFF_indx = [], []
    for i in range(len(df.columns)):
        if 'num_aps' in df.columns[i]:
            num_aps_indx.append(i)
        if 'IFF' in df.columns[i]:
            IFF_indx.append(i)
    return num_aps_indx, IFF_indx

def plot_IFF_distribution(IFF_df, data_type, DV):
    ''' 
    plots distribution of number of APs and initial firing frequencies at each current injection
    on D1 and D2
    data_type - all, repatch, repatch firing cells D1
    DV - (dependent variable): num_aps or IFF
    '''
    
    treatments = ['Ctrl', 'TTX', 'high K']
    num_aps_indx, IFF_indx = get_num_aps_and_IFF_data_culumns(IFF_df)
    date = str(datetime.date.today())
    destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/IFF/'
    kwargs = dict(alpha=0.6, bins=7)

    colors_dict = {'D1' : sns.color_palette("colorblind")[1], 'D2': sns.color_palette("colorblind")[2]}
    DV_dict = {'IFF' : [IFF_indx, 'Initial firing frequency (AP#1 to AP#2)',  'Firing frequency (Hz)'], 
    'num_aps' : [num_aps_indx, 'Number of fired action potentials', 'AP count']}

    for treatment in treatments:
        tr_df = IFF_df[IFF_df['treatment'] == treatment]

        fig, ax = plt.subplots(5,5,sharex = False, sharey = True ,figsize=(15,25))
        fig.subplots_adjust(hspace=0.3,  wspace = None)
        ax = ax.flatten()
        for day in tr_df['day'].unique():
            day_df = tr_df[tr_df['day'] == day]
            for i, col in enumerate(DV_dict[DV][0][1:]):
                plot_data = day_df.iloc[:,col]
                ax[i].hist(plot_data, **kwargs, label  = day)

                ax[i].set_ylabel('Frequency')
                ax[i].set_xlabel(DV_dict[DV][2])
                stim = day_df.columns[col][2:(day_df.columns[col]).rfind('pA')+2]
                ax[i].set_title(stim)
                ax[i].spines['top'].set_visible(False)
                ax[i].spines['right'].set_visible(False)

        plt.figlegend(loc = 'upper right', bbox_to_anchor=(0.95, 0.97), fontsize = 15)
        fig.patch.set_facecolor('white')

        fig.suptitle(DV_dict[DV][1] + ' ' +  data_type + ' ' + treatment, size = 20)

        plt.savefig(destination_dir + date + treatment + data_type + DV + '_distibution_plot.png')
        plt.close(fig)


def plot_IFF_avg_against_current_inj(IFF_df, data_type, DV, 
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/IFF/'):

    ''' 
    data_type - all, repatch, repatch firing cells D1
    DV - (dependent variable): num_aps or IFF
    for each current injection >0 plots the DV for D1 and D2 with
    mean with standard error around the mean
    '''
    
    num_aps_indx, IFF_indx = get_num_aps_and_IFF_data_culumns(IFF_df)
    treatments = ['Ctrl', '8 mM', '15 mM']
    #treatments = ['Ctrl', 'TTX', 'high K']
    date = str(datetime.date.today())

    colors_dict = {'D1' : ['#dede00','#dede00','#dede00'], 'D2': ['#ff7f00', '#4daf4a','#984ea3']}
    DV_dict = {'IFF' : [IFF_indx, 'Initial firing frequency (AP#1 to AP#2)',  'Instantaneous firing frequency (Hz)'], 
    'num_aps' : [num_aps_indx, 'Number of fired action potentials', 'AP count']}

    fig, ax = plt.subplots(1, len(treatments), figsize=(21,8), sharey = True)

    for k, treatment in enumerate(treatments):
        tr_df = IFF_df[IFF_df['treatment'] == treatment]

        for day in tr_df['day'].unique():
            day_df = tr_df[tr_df['day'] == day]
            avgs, sems, inj = [], [], []
            for i, col in enumerate(DV_dict[DV][0][5:]):
                data = day_df.iloc[:,col]
                x = np.linspace(0.75+i, 1.25+i, len(data))
                ax[k].scatter(x, data, alpha = 0.7, s = 7, c = colors_dict[day][k])

                avg = np.mean(data)
                sem = np.std(data, ddof=1) / np.sqrt(np.size(data))

                #yerr = np.linspace((avg - sem), (avg + sem), 5) 
                ax[k].errorbar(i + 1, avg, yerr = sem)

                avgs.append(avg)
                sems.append(sem) #standard error of the mean
                inj.append(day_df.columns[col][2:(day_df.columns[col]).rfind('pA')])
            
            ax[k].scatter(range(1, len(inj)+1), avgs, label = day, color = colors_dict[day][k], s = 30)
            ax[k].plot(range(1, len(inj)+1), avgs, label = day, color = colors_dict[day][k])
            ax[k].errorbar(range(1, len(inj)+1), avgs, yerr = sem, label = day, color = colors_dict[day][k])

            ax[k].set_title(treatment + ' ' + data_type, fontsize = 25)
            ax[k].set_xlabel('Current (pA)', fontsize = 24)
            ax[0].set_ylabel(DV_dict[DV][2], fontsize = 24      )
            ax[k].set_xticks(ticks = list(range(1,22,2)), labels = [0,100,200,300,400,500,600,800,1000,1200,1400],
            fontsize = 22, rotation = 45) 
            max_val = max(tr_df.iloc[:,DV_dict[DV][0]].max(numeric_only = True))
            ax[k].set_yticks(ticks = [0, 40, 80, 120, 160, 200], fontsize = 22) 
            if DV == 'num_aps':
                ax[k].set_yticks(ticks = [0, 10, 20, 30, 40, 50], fontsize = 22) 

            ax[k].tick_params(axis='y', labelsize = 22)
            ax[k].spines['top'].set_visible(False)
            ax[k].spines['right'].set_visible(False)

        #max_val = max(tr_df.max(numeric_only = True))
        if 'repatch' in data_type:
            ax[k].text(len(inj)-4, 200, 'n = ' + str(int(len(tr_df)/2)), size = 25, c = 'k')
        else: 
            ax[k].text(len(inj)-4, 200, 'n = ' + str(int(len(day_df))), size = 25, c = 'k')

            # add n numbers
            #ax.text(1.7 + 2*i, int(np.max(df[param])), 'n = ' + str(len(df_plot)), size = 10, c = colors[2*i+1])
    #plt.figlegend(loc = 'upper right', bbox_to_anchor=(1, 1), fontsize = 25)
    fig.patch.set_facecolor('white')
    fig.suptitle(DV_dict[DV][1], size = 30)

    fig.tight_layout()
    plt.savefig(destination_dir + date + data_type + DV + '_vs_curernt_inj_plot.png')
    plt.close(fig)



#%%
# connectivity ploting funcs (summary)

plt.style.use('./style_plot_intrinsic.mplstyle')

def get_QC_connectivity_df(df):
    mask = (df['Vm pre'] < -50 ) & (df['Vm post'] < -50 ) & \
        (df['Amp 1'] > 0) &  (df['Amp 2'] > 0 )   & \
            (df['Amp 3'] > 0 ) &  (df['Amp 4'] > 0)
    df = df.loc[mask, :]
    return df

def get_QC_connectivity_VC(df):
    mask = (df['Vm pre'] < -50 )
    df = df.loc[mask, :]
    return df

def get_repatch_connectivity_df(df):

    #take only repatched cells
    repatch_df = df[df['repatch'] == 'yes']
    repatch_df.reset_index(inplace = True, drop = True)

    #check if cell IDs appear twice, if not remove
    not_repatched_cells = []
    for cell in repatch_df['connection_ID'].unique():
        if len(repatch_df[repatch_df['connection_ID'] == cell]) != 2:
            not_repatched_cells.append(cell)

    for cell in not_repatched_cells:
        repatch_df = repatch_df.drop(repatch_df.index[repatch_df['connection_ID'] == cell])     
    repatch_df.reset_index(inplace = True, drop = True)

    return repatch_df

def get_amps_lats_columns(df):
    amps, lats = [], []
    for i in range(len(df.columns)):
        if 'Amp' in df.columns[i]:
            amps.append(i)
        if 'Lat' in df.columns[i]:
            lats.append(i)
    return amps, lats

def plot_connect_amplitude(df, data_type, op_color_dict, results_ = 'amp',
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/connectivity/'):

    treatments = ['Ctrl', 'TTX', 'high K']
    amps, lats = get_amps_lats_columns(df)
    plot_amps = df.columns[amps[0]:amps[3]+1]
    label_y = 'Amplitude (mV)'
    y_lim = 'Amp 1'
    if results_ == 'lat':
        plot_amps = df.columns[lats[0]:lats[3]+1]
        label_y = 'Delay (ms)'
        y_lim = 'Lat4'

    colors = ['#dede00', '#ff7f00', '#dede00', '#DA67BA', '#dede00', '#7FEE8F','#dede00']
   #colors = ['#dede00', '#ff7f00', '#dede00', '#4daf4a', '#dede00', 'mediumpurple']
    cmap = plt.cm.get_cmap('tab20')
    # op_color_dict = {}
    # for h, op in enumerate(df['OP'].unique()):
    #     op_color_dict[op] = cmap((h+1)/10)  

    fig, ax = plt.subplots(2,2, figsize = (16,16))
    ax = ax.flatten()
    for a, amp in enumerate(plot_amps): 
        
        day_label = []
        num_cels = {}
        
        for i, tr in enumerate(treatments):
            x_plot =[]
            for j, day in enumerate(sorted(df['day'].unique())):
                k = j + 2*i 
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                median = df_plot[amp].median()
                x = np.linspace(0.65+k, 1.35+k, len(df_plot))
                x_plot.append(x)
                ax[a].scatter(x, df_plot[amp], alpha = 0.9, c = colors[int(k)], s = 60)
                if len(df_plot) == 0:
                    continue
                yerr = 1.253*(df_plot[amp].std()/(math.sqrt(len(df_plot))))
                ax[a].scatter(1+k, median, color = 'k', marker = '_', s = 2000)
                ax[a].text(0.9+k, 1.2*median, str(round(median,2)), size = 10)
                day_label.append(day)
                num_cels[tr + ' ' + day] = len(df_plot)
                #data_boxplot.append(df_plot[amp)
            ax[a].text(1 + 2*i, int(np.max(df['Amp 1'])+0.6), tr, size = 17, c = colors[2*i+1])
            ax[a].text(1.7 + 2*i, int(np.max(df['Amp 1'])-0.3), 'n = ' + str(len(df_plot)), size = 10, c = colors[2*i+1])
            if k in [1,3,5] and 'repatch' in data_type:
                for c, cell in enumerate(df_plot['connection_ID']):
                    if k == 1 and len(df_plot['connection_ID'])-1 == c:
                        print(str(c+1) + 'Gencho')
                    x1 = [x_plot[0][c], x[c]] 
                    y = df[amp][df['connection_ID'] == cell]
                    op = df_plot['OP'][df_plot['connection_ID'] == cell].tolist()[0]
                    #age = df_plot['_new'][df_plot['connection_ID'] == cell].tolist()[0]
                    if type(list(op_color_dict.keys())[0]) is not str:
                        ax[a].plot(x1, y, '-', color = op_color_dict[age], alpha = 0.5, linewidth = 2, zorder = -1)
                    else:
                        ax[a].plot(x1, y, '-', color = op_color_dict[op], alpha = 0.5, linewidth = 2, zorder = -1)


        
        if 'VC' in data_type:
            ax[a].set_ylabel('Amplitude (pA)', fontsize = 24)
        else:
            ax[a].set_ylabel(label_y, fontsize = 20)
        if results_ == 'amp':
            ax[a].set_ylim([-0.5, int(np.max(df['Amp 1']))+1])
        else:
            max_ = np.max([np.max(df).Lat1, np.max(df).Lat2,np.max(df).Lat3,np.max(df).Lat4]) + 0.5
            ax[a].set_ylim([-0.5, int(max_)])
        #ax[a].set_xticks(ticks = list(range(1,len(day_label) +1)), labels = day_label) 
        ax[a].set_xticks(ticks = [1, 2, 3, 4, 5, 6], labels = ['D1', 
        'D2 \n Ctrl', 'D1', 'D2 \n TTX','D1', 'D2 \n high K'], size = 20)
        #ax[a].set_title('#' + str(a+1), fontsize = 15)
        ax[a].set_xlabel('Day', fontsize = 15)
        ax[a].tick_params(axis='y', labelsize=22) 
        
    fig.subplots_adjust(hspace = 0.4, wspace = 0.4)
    
    fig.suptitle('Postsynaptic responses', fontsize = 20)
    #fig.patch.set_facecolor('white')
    date = str(datetime.date.today())
    plt.savefig(destination_dir  + date + data_type +'_postsynaptic_responses.png')
    plt.close(fig)


###############
#%%
# Plots including TTX trreatment

def plot_param_for_days_slice_TTX_incl(df, op_color_dict,
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/intrinsic_properties/adult_slice/TTX_incl/'):

    treatments = ['Ctrl', 'TTX', 'high K']
    titles_dict = dict_for_plotting()

    colors = ['#dede00', '#ff7f00', '#dede00', '#DA67BA', '#dede00', '#7FEE8F','#dede00','#186F25'] #Ctrl, TTX, 8mM Kcl, 15 mM Kcl

    for u, param in enumerate(titles_dict.keys()): 
        fig2 = plt.figure(figsize=(11,7))
        ax = plt.subplot(1,1,1)
        day_label = []
        num_cels = {}
        data_boxplot = []
        for i, tr in enumerate(treatments):
            for j, day in enumerate(df['day'].unique()):
                k = j + 2*i 
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                median = df_plot[param].median()
                x = np.linspace(0.7+k, 1.2+k, len(df_plot))
                if tr == 'high K' and day == 'D2':
                    #ax.scatter(x, df_plot[param], alpha = 0.7, s = 40, c = df_plot['high K concentration'], cmap = colormap_K)
                    y_8mm = df_plot[param].loc[df_plot['high K concentration'] == '8 mM'] 
                    x_8 = np.linspace(0.7+k, 1.05+k, len(y_8mm))
                    median_8 = np.median(y_8mm)
                    y_15mm = df_plot[param].loc[df_plot['high K concentration'] == '15 mM'] 
                    x_15 = np.linspace(1.1+k, 1.25+k, len(y_15mm))
                    median_15 = np.median(y_15mm)
                    ax.scatter(x_8, y_8mm, c = colors[5], s = 60, alpha = 0.5)
                    ax.scatter(x_15, y_15mm, c = colors[7], s = 60, alpha = 0.5)
                    #ax.boxplot(y_8mm, positions = [k+0.6], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[5], alpha = 0.8),
                    #medianprops = dict(linewidth=2.3, color = 'k'))
                    #ax.boxplot(y_15mm, positions = [k+1.55], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[7], alpha = 0.8),
                    #medianprops = dict(linewidth=2.3, color = 'k'))
                else:
                    #ax.boxplot(df_plot[param], positions = [k + 0.5], notch = True, patch_artist=True, boxprops=dict(facecolor=colors[k], alpha = 0.8),
                    #medianprops = dict(linewidth=2.3, color = 'k'))    
                    ax.scatter(x, df_plot[param], c = colors[int(k)], s = 60, alpha = 0.5)
                
                # yerr = 1.253*(df_plot[param].std()/(math.sqrt(len(df_plot))))
                ax.scatter(1+k, median, color = 'k', marker = '_', s = 2000)
                ax.text(0.9+k, (1.05*median), str(round(median,2)), size = 20)
                day_label.append(day)
                num_cels[tr + ' ' + day] = len(df_plot)
                #data_boxplot.append(df_plot[param])
            #ax.text(1 + 2*i, int(np.max(df[param])), tr, size = 17, c = colors[2*i+1])
            #ax.text(1.7 + 2*i, int(np.max(df[param])), 'n = ' + str(len(df_plot)), size = 12, c = colors[2*i+1])

        #plt.boxplot(data_boxplot, showbox = False)
        ax.tick_params(axis='y', labelsize=22)
        ax.set_xticks(ticks = [0.7, 1.7, 2.7, 3.7, 4.7, 5.8, 6.55], labels = ['D1', 
        'D2 \n after Ctrl', 'D1', 'D2 \n after TTX','D1', 'D2 \n 8mM','D2 \n 15mM'], size = 20)

        #plt.title(titles_dict[param][0], fontsize = 19, x = 0.5, y = 1)
        ax.set_ylabel(titles_dict[param][1], fontsize = 24)
        ax.set_yticks(ticks = titles_dict[param][2])

        # plt.subplots_adjust(hspace=0.35)
        # fig2.patch.set_facecolor('white')
        # fig2.tight_layout()

        date = str(datetime.date.today())
        plt.savefig(destination_dir  + date + '_plot_' + param + '.pdf')
        plt.close(fig2)

def plot_param_for_days_repatch_plus_TTX(df, op_color_dict,
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/plots/intrinsic_properties/adult_repatch/TTX_included/'):
    
    treatments = ['Ctrl', 'TTX', 'high K']
    titles_dict = dict_for_plotting()
        
    colors = ['#dede00', '#ff7f00', '#dede00', '#DA67BA', '#dede00', '#7FEE8F','#dede00','#186F25'] #Ctrl, TTX, 8mM Kcl, 15 mM Kcl
    for u, param in enumerate(titles_dict.keys()): 
        fig2 = plt.figure(figsize=(11,7))
        ax = plt.subplot(1,1,1)
        day_label = []
        num_cels = {}
        for i, tr in enumerate(treatments):
            x_plot =[]
            
            for j, day in enumerate(df['day'].unique()):
                k = j + 2*i  
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                df_plot.reset_index(drop=True, inplace=True)
                median = df_plot[param].median()
                x = np.linspace(0.7+k, 1.3+k, len(df_plot))
                x_plot.append(x)

                if tr == 'high K' and day == 'D2':
                    y_8mm = df_plot.loc[df_plot['high K concentration'] == '8 mM']
                    x_8 = np.linspace(0.8+k, 1.05+k, len(y_8mm))
                    #x_8 = x[:len(y_8mm)]
                    median_8 = np.median(y_8mm[param])
                    y_15mm = df_plot.loc[df_plot['high K concentration'] == '15 mM'] 
                    x_15 = np.linspace(1.4+k, 1.65+k, len(y_15mm))
                    #x_15 = x[-len(y_15mm):]
                    median_15 = np.median(y_15mm[param])
                    ax.scatter(x_8, y_8mm[param], c = colors[5], s = 80, zorder = 2, alpha = 0.8)
                    ax.scatter(x_15, y_15mm[param], c = colors[7], s = 80, zorder = 2, alpha = 0.8)
                    ax.plot([0.7+k, 1+k], [median_8, median_8], c = 'k', linestyle = 'solid')
                    ax.plot([1.3+k, 1.6+k], [median_15, median_15], c = 'k', linestyle = 'solid')
                    ax.text(k+0.65, median_8 + 0.5, str(round(median_8, 2)), fontsize = 20)
                    ax.text(k+1.25, median_15 + 0.5, str(round(median_15, 2)), fontsize = 20)
                    df_plot.insert(0, 'x', np.concatenate([x_15, x_8]))

                    for c, cell in enumerate(df_plot['cell_ID_new']):
                        #indx = index_s[c]
                        x_K = df_plot['x'].loc[df_plot['cell_ID_new'] == cell].tolist()[0]
                        x1 = [x_plot[0][c], x_K]
                        y = df[param][df['cell_ID_new'] == cell]
                        op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
                        age = df_plot['patient_age_num'][df_plot['cell_ID_new'] == cell].tolist()[0]
                        if type(list(op_color_dict.keys())[0]) is not str :
                            plt.plot(x1, y, '-', colors[k], alpha = 0.3, linewidth = 2, zorder = 1)
                        else:
                            plt.plot(x1, y, '-', colors[k], alpha = 0.3, linewidth = 2, zorder = 1)


                else:
                    ax.scatter(x, df_plot[param], c = colors[int(k)], s = 80, zorder = 2, alpha = 0.8)
                    ax.scatter(1+k, median, color = 'k', marker = '_', s = 2000, zorder = 2, alpha = 0.8)
                    ax.text(0.75+k, (1.05*median), str(round(median,2)), size = 15)
                    x_plot.append(x)
                    if k in [1,3,5]:
                        for c, cell in enumerate(df_plot['cell_ID_new']):
                            #indx = index_s[c]
                            x1 = [x_plot[0][c], x[c]] 
                            y = df[param][df['cell_ID_new'] == cell]
                            op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
                            age = df_plot['patient_age_num'][df_plot['cell_ID_new'] == cell].tolist()[0]
                            if type(list(op_color_dict.keys())[0]) is not str :
                                plt.plot(x1, y, '-', color = colors[k], alpha = 0.3, linewidth = 2, zorder = -1)
                            else:
                                plt.plot(x1, y, '-', color = colors[k], alpha = 0.3, linewidth = 2, zorder = -1)


                day_label.append(day)
                num_cels[tr + ' ' + day] = len(df_plot)

                ax.text(1 + 2*i, int(np.max(df[param])), tr, size = 17, c = colors[2*i+1])
                ax.text(1.7 + 2*i, int(np.max(df[param])), 'n = ' + str(len(df_plot)), size = 12, c = colors[2*i+1])
    
        
        ax.set_xticks(ticks = [1, 2, 3, 4, 5, 5.8, 6.55], labels = ['D1', 
        'D2 \n after Ctrl', 'D1', 'D2 \n after TTX','D1', 'D2 \n 8mM','D2 \n 15mM'], size = 20)
        ax.set_ylabel(titles_dict[param][1], fontsize = 24)
        #ax.set_yticks(ticks = titles_dict[param][2])
        ax.tick_params(axis='y', labelsize=22) 
        #plt.title(titles_dict[param][0], fontsize = 19, x = 0.5, y = 1)
        plt.subplots_adjust(hspace=0.35)
        #fig2.tight_layout()

        date = str(datetime.date.today())
        plt.savefig(destination_dir + date + '_plot_' + param + '.pdf')
        plt.close(fig2)

 # %%
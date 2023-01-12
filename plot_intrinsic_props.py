#%%


import matplotlib.pyplot as plt 
import funcs_for_results_tables as get_results
import pandas as pd
import datetime
import numpy as np
import math

df_intr_props = get_results.collect_intrinsic_df()

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
    age_filtered_df.loc[mask_adult, 'patient_age'] = 'A'
    age_filtered_df.loc[mask_juv, 'patient_age'] = 'J'

    return age_filtered_df


def change_to_numeric(df):
    '''
    takes df intrinsic properties  
    changes datatype of columns 11 and 13:25 to numeric
    '''
    df.assign(hrs_incubation = pd.to_numeric(df['hrs_incubation'], errors = 'coerce'))
    for col in df.columns[13:25]:
        df.assign(col = pd.to_numeric(df[col], errors = 'coerce'))
    #check dtypes with display(df.dtypes)
    return df

def get_QC_data(df):
    mask = (df['Rs'] < 30) & \
        (df['resting_potential'] < -45 ) & (df['resting_potential'] > -90) & \
            (df['max_repol'] > -180) & \
                (df['membra_time_constant_tau'] > -25 ) &\
                    (df['capacitance'] < 900) &  (df['capacitance'] > 10) & \
                        (df['AP_heigth'] > 30) & \
                            (df['max_repol'] > -100) & \
                                (df['Rheobase'] > 0)

    #include filter for TH ?
    df = df.loc[mask, :]
    return df

def filter_on_hrs_incubation(df, min_inc_time):
    mask = (df['hrs_incubation'] == 0) | (df['hrs_incubation'] >= min_inc_time)
    df = df.loc[mask, :]
    return df

def get_adult_data(df_intr_props):
    adult_df = patient_age_to_float(df_intr_props, 18)
    adult_df = change_to_numeric(adult_df)
    adult_df = get_QC_data(adult_df)
    adult_df = filter_on_hrs_incubation(adult_df, 20)
    return adult_df

def get_repatch_df(df):
    repatch_df = df[df['repatch'] == 'yes']
    repatch_df.reset_index(inplace = True, drop = True)

    patcher_dict = {'Verji':'vm', 'Rosie': 'rs'}
    cell_IDs_new = []
    for i in range(len(repatch_df)):
        cell_ID = patcher_dict[repatch_df['patcher'][i]] + repatch_df['cell_ID'][i]
        cell_IDs_new.append(cell_ID)
    repatch_df['cell_ID_new'] = cell_IDs_new

    not_repatched_cells = []
    for cell in repatch_df['cell_ID_new'].unique():
        if len(repatch_df[repatch_df['cell_ID_new'] == cell]) < 2:
            not_repatched_cells.append(cell)

    for cell in not_repatched_cells:
        repatch_df = repatch_df.drop(repatch_df.index[repatch_df['cell_ID_new'] == cell])
        
    repatch_df.reset_index(inplace = True, drop = True)
    return repatch_df

def dict_for_plotting():
    titles_dict = {'Rs': ['Series resistance', 'Ω'], 
    'Rin': ['Input resistance', 'Ohm'],
    'resting_potential': ['Resting membrane potential', 'mV'],
    'max_spikes': ['Maximum number of spikes', 'count'],
    'Rheobase': ['Rheobase', 'pA'],
    'AP_heigth': ['AP amplitude', 'mV'],    
    'TH': ['AP threshold', 'mV'],
    'max_depol': ['AP upstroke (30% to 70%)', 'mV/ms'],
    'max_repol': ['AP downstroke (70% to 30%)', 'mV/ms'],
    'membra_time_constant_tau': ['Membrane time constant', 'tau'],
    'capacitance': ['Capacitance', 'F']}
    return titles_dict

#%%

def plot_time_after_OP_vs_param(df, param,
destintaion_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/intrinsic_properties/hrs_after_op_dependency/'): 
    titles_dict = dict_for_plotting()

    fig1 = plt.figure(figsize=(10,7))
    ax = plt.subplot(1, 1, 1)
    plt.subplots_adjust(hspace=0.5)
    plt.scatter(df['hrs_after_OP'], df[param])

    plt.title('AP ' + titles_dict[param][0] + ' time dependency', fontsize = 19)
    ax.set_xlabel('Time after cortex resection (hrs)', fontsize = 15)
    ax.set_ylabel(titles_dict[param][1], fontsize = 15)

    date = str(datetime.date.today())
    
    fig1.patch.set_facecolor('white')
    plt.savefig(destintaion_dir + date + param +'_vs_time_after_op.png')
    plt.close(fig1)


def plot_param_for_days(df,
destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/intrinsic_properties/repatch_all/'):

    treatments = ['Ctrl', 'TTX', 'high K']
    titles_dict = dict_for_plotting()

    #destination_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/intrinsic_properties/'
    colors = ['moccasin', 'red', 'moccasin', 'cadetblue', 'moccasin', 'mediumpurple']
    cmap = plt.cm.get_cmap('tab20')
    op_color_dict = {}
    for h, op in enumerate(df['OP'].unique()):
        op_color_dict[op] = cmap((h+1)/10)

    for param in titles_dict.keys(): 
        fig2 = plt.figure(figsize=(10,5))
        ax = plt.subplot(1,1,1)
        day_label = []
        num_cels = {}
        data_boxplot = []
        for i, tr in enumerate(treatments):
            x_plot =[]
            for j, day in enumerate(df['day'].unique()):
                k = j + 2*i 
                df_plot = df[(df['treatment'] == tr) & (df['day'] == day)]
                median = df_plot[param].median()
                x = np.linspace(0.65+k, 1.35+k, len(df_plot))
                x_plot.append(x)
                ax.scatter(x, df_plot[param], alpha = 0.8, c = colors[int(k)], s = 40)
                yerr = 1.253*(df_plot[param].std()/(math.sqrt(len(df_plot))))
                ax.scatter(1+k, median, color = 'k', marker = '_', s = 2000)
                #ax.text(0.9+k, median + 2, str(round(median,2)), size = 12)
                day_label.append(day)
                num_cels[tr + ' ' + day] = len(df_plot)
                data_boxplot.append(df_plot[param])
            ax.text(1 + 2*i, int(np.max(df[param])), tr, size = 17, c = colors[2*i+1])
            ax.text(1.7 + 2*i, int(np.max(df[param])), 'n = ' + str(len(df_plot)), size = 10, c = colors[2*i+1])
            if k in [1,3,5]:
                for c, cell in enumerate(df_plot['cell_ID_new']):
                    x1 = [x_plot[0][c], x[c]] 
                    y = df[param][df['cell_ID_new'] == cell]
                    op = df_plot['OP'][df_plot['cell_ID_new'] == cell].tolist()[0]
                    plt.plot(x1, y, '-', color = op_color_dict[op], alpha = 0.5)

        #plt.boxplot(data_boxplot, showbox = False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticks(ticks = list(range(1,7)), labels = day_label) 
        plt.title(titles_dict[param][0], fontsize = 19)
        ax.set_xlabel('Day', fontsize = 15)
        ax.set_ylabel(titles_dict[param][1], fontsize = 15)

        fig2.patch.set_facecolor('white')
        plt.savefig(destination_dir  + 'Repatch_plot_' + param + '.png')
        plt.close(fig2)

    #
# %%
#MAIN

adult_df = get_adult_data(df_intr_props)
adult_df = adult_df[adult_df['area'] == 'temporal']

repatch_df = get_repatch_df(adult_df)
repatch_df = repatch_df.sort_values(['cell_ID_new', 'treatment'])
plot_param_for_days(repatch_df) #for all params repatch!






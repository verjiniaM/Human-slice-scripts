import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#To run on the data table of all data
#cells that were repatched should be indicated
#%%
#%%

filename1 = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/OP211209/data_tables/Intrinsic_and_synaptic_properties copy.xlsx'
filename2 = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/OP211220/data_tables/Intrinsic_and_synaptic_properties copy.xlsx'
filename3 = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/OP211123/data_tables/Intrinsic_and_synaptic_properties copy.xlsx'
fn_RS = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/Intrinsic_and_synaptic_properties_rs_220113_Repatch.xlsx'

df1 = pd.read_excel(filename1, header = 0)
df2 = pd.read_excel(filename1, header = 0)
df3 = pd.read_excel(filename1, header = 0)
df_RS = pd.read_excel(fn_RS, header = 1)
df = pd.concat([df1, df2, df3], axis = 0, ignore_index=True)

# index_repatch = df.index[df['Repatch'] == 'y'].tolist()
# df_repatch = df.loc[index_repatch]

# df_D1 = df_repatch.iloc[:,4:17]
# df_D2 = df_repatch.iloc[:, 24:37]
# df_D2.columns = df_D1.columns
# df_combined = pd.concat([df_D1, df_D2], axis=0, ignore_index=True)

# day = []
# for i in range(0,int(df_combined.shape[0]/2)):
#     day.append('D1')
# for j in range(int(df_combined.shape[0]/2),int(df_combined.shape[0])):
#     day.append('D2')
# df_combined.insert(0, 'Day', day)
# treat = df_repatch['Treatment'].tolist() + df_repatch['Treatment'].tolist()
# df_combined['Treatment'] = treat
# df_RS = pd.concat([df_RS.iloc[:,4:18], df_RS['Treatment']], axis=1)

# df_all = pd.concat([df_combined, df_RS], axis=0, ignore_index=True)

df_all = df
treatments = 'TTX', 'high K', 'Ctrl'
for k in range(len(treatments)):
    tr = treatments[k]
    index_tr = df_all.index[df_all['treatment'] == tr].tolist()
    df_tr = df_all.loc[index_tr]
    #plt.figure()
    fig, ax = plt.subplots(3,4,sharex=True, sharey=True,figsize=(6,6))
    for i in range(0,12):
        ax = fig.add_subplot(3,4, i+1)
        ax.scatter(df_tr['day'], df_tr.iloc[:,i+8].tolist())
        avg_D1 = np.nanmean(df_tr.iloc[:,i+8].where(df_all['day'] == 'D1').tolist())
        avg_D2 = np.nanmean(df_tr.iloc[:,i+8].where(df_all['day'] == 'D2').tolist())
        ax.scatter(df_tr['day'][0], avg_D1, color='red', lw=6, ls='--', label="average plot")
        ax.scatter(df_tr['day'].iloc[-1], avg_D2, color='red', lw=6, ls='--', label="average plot")
        ax.title.set_text(df_tr.columns[i+8])
        # ax.set_xticks([])
        # ax.set_yticks([])
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        # # ax.spines['bottom'].set_visible(False)
        # # ax.spines['left'].set_visible(False)
    fig.patch.set_facecolor('white')
    fig.suptitle(tr, fontsize=15)
    fig.tight_layout()
    plt.show()
    
        


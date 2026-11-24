import pandas as pd

df_together = pd.read_csv('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/connectivity_paper/data/all_intrinsic.csv')
df_morpho = pd.read_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/connectivity_paper/data/morpho.xlsx')


# remove OPs that only in df_morpho but did not pass QC
OPs_remove = set(df_morpho.OP.unique()) - set(df_together.OP.unique())
df_morpho = df_morpho[~df_morpho['OP'].isin(OPs_remove)].reset_index(drop=True)

# see all unique combinations of cell_ID, opatient, day, slice 
# this should give 0 rows; if yes tagble ready for merge
df_morpho[df_morpho.duplicated(["cell_ID", "day", "slice", 'patcher'], keep=False)]


keys = ["cell_ID", "day", "slice", 'patcher']

df_join = df_together.merge(
    df_morpho,
    on=keys,
    how="left",
    suffixes=("", "_y")
)

print(len(df_join.columns))

shared = (set(df_together.columns) & set(df_morpho.columns)) - set(keys)

for col in shared:
    ycol = f"{col}_y"
    if ycol in df_join.columns and df_join[col].equals(df_join[ycol]):
        df_join = df_join.drop(columns=[ycol])
print(len(df_join.columns))
df_join = df_join.loc[:, ~ df_join.columns.str.contains('Unnamed')]
print(len(df_join.columns))

df_join.to_excel('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/connectivity_paper/data/complete_morpho_intr.xlsx')


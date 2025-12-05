import os
import glob
import pandas as pd

dir_to_search = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/imaging/H/'

def find_op_tiff_case_insensitive(directory):
    """Find .tiff files containing 'OP' (case-insensitive)"""
    op_tiff_files = []
    
    # Search for all .tiff and .TIFF files
    patterns = [
        os.path.join(directory, '**', '*.tif'),
        os.path.join(directory, '**', '*.lif'),
        os.path.join(directory, '*.tif'),
        os.path.join(directory, '*.lif'),
    ]
    
    all_files = []
    for pattern in patterns:
        all_files.extend(glob.glob(pattern, recursive=True))
    
    # Filter for files containing 'OP' (case-insensitive)
    for file_path in all_files:
        filename = os.path.basename(file_path)
        if 'OP' in filename.upper():  # Case-insensitive check
            op_tiff_files.append(file_path)
    
    return sorted(list(set(op_tiff_files)))

# Execute search
op_tiff_files = find_op_tiff_case_insensitive(dir_to_search)
print(f'found {len(op_tiff_files)} files')

imaged_op = []
for fn_path in op_tiff_files:
    imaged_op.append(fn_path[fn_path.rfind('OP'):(fn_path.rfind('OP')+8)])

imaged_op = sorted(list(set(imaged_op)))

# load ephys data
data_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/paper_figs_collected_checked/data/'
inc_df = pd.read_excel(data_dir + 'incubation_only_temporal.xlsx')
slice_df = pd.read_excel(data_dir + 'slice_data_temporal.xlsx')
repatch_df = pd.read_excel(data_dir + 'repatch_data_temporal.xlsx')
all_df = pd.concat([inc_df, slice_df, repatch_df])

op_list_ephys = []
for df in [inc_df, slice_df, repatch_df]:
    op_list_ephys.extend(df.OP.unique().tolist())

op_list_ephys = sorted(list(set(op_list_ephys)))

print(sorted(list(set(imaged_op) & set(op_list_ephys))))

ops_interest = ['OP211123', 'OP211209', 'OP220127', 'OP230808', 'OP240926']
for op in ops_interest:
    print(f'inc: {all_df.hrs_incubation[all_df.OP == op].tolist()[-1]}')
    print(f'after OP: {all_df.hrs_after_OP[all_df.OP == op].tolist()[-1]}')
    

df_qc = pd.DataFrame(columns=['OP','patcher', 'filename', 'slice', 'cell_ID', 'cell_ch', 'Rs_start', 'Rin_start', 'Rs_end', 'Rin_end', 
    'change_rs', 'change_rin'])

for i in all_vc_indices:
    slic = slic = slice_names[i]
    day = 'D2'
    active_chans = [int(item) for item in input('channels in ' + filenames[i]).split()]
    filename_vc = work_dir + filenames[i]
    cell_IDs = hcf.get_cell_IDs(filename_vc, slic, active_chans)
    Rs, Rin = hcf.get_access_resistance(filename_vc, active_chans) 
    params1_df = pd.DataFrame({'filename': filenames[i], 'slice' : slic, 'cell_ch': active_chans,
        cell_ID:'cell_IDs',  'Rs_start' : Rs, 'Rin_start': Rin})
    df_qc = pd.concat([df_qc.loc[:], data_to_add]).reset_index(drop=True)


df_qc.to_excel(work_dir + 'data_tables/' + OP + '_QC_measures_rs.xlsx') 



fname = work_dir + OP + '_indices_dict.json'
out_data = indices_dict

with open(fname, "w") as outfile:
    json.dump(out_data , outfile)

file_out = '_meta_active_chans.json'
sort.to_json (work_dir, OP, file_out, active_chans_meta[0], active_chans_meta[1], active_chans_meta[2])
active_chans_meta[2]
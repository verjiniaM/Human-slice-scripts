import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA
from itertools import permutations
from itertools import product
from scipy.optimize import curve_fit
from detecta import detect_peaks
import ephys_analysis.funcs_sorting as sort
import ephys_analysis.funcs_human_characterisation as hcf


def load_traces_FM(filename):
    '''
    returns a dictionary; key: channel name (e.g. 'Ch1')
    value: np.array(sweep_count, sweep_len), 
    each row is a sweep, len column - number of sweeps
    '''
    block = hcf.read_abf(filename)
    sweep_count = len(block.segments) #number of sweeps
    sweep_len = len(block.segments[0].analogsignals[0])
    channels = len(block.segments[0].analogsignals)

    data_dict = {}
    for ch in range(0, channels):
        name = block.segments[0].analogsignals[ch].name

        ch_data = np.ndarray([sweep_count, sweep_len])
        for i in range(0,sweep_count):
            ch_data[i] = block.segments[i].analogsignals[ch].view(np.recarray).reshape(sweep_len)

        data_dict[name]=[ch_data]
    return data_dict
    # len(data_dict) # number of channels
    # np.shape(data_dict1['Ch1'][0])[0] - sweep count
    # np.shape(data_dict1['Ch1'][0])[1] - sweep length

def get_RMP_FM(filename, active_chans):
    '''
    take average amplitudes of four 100-ms trace segments, 
    located before the –50, 50, 100 and 150-pA step current injections
    average to get RMP
    '''
    vals_for_RMP = [-50, 50, 100, 150]
    charact_data = load_traces_FM(filename)
    inj = hcf.get_inj_current_steps(filename)
    indexes = np.array([inj.index(value) for value in vals_for_RMP if value in inj])

    RMPs_all = []
    for ch in active_chans:
        key = 'Ch' + str(ch)
        chan_data = charact_data[key][0]
        rmp = np.mean(chan_data[indexes, 0:2000])

        RMPs_all.append(rmp)

    return RMPs_all

def get_Rin_FM(filename, active_chans):
    '''
    get steady state for -100, -50, 50, 100 injections
    '''
    vals_for_Rin = [-100, -50, 50, 100]
    charact_data = load_traces_FM(filename)
    inj = hcf.get_inj_current_steps(filename)
    indexes = np.array([inj.index(value) for value in vals_for_Rin if value in inj])

    rmp = get_RMP_FM(filename, active_chans)

    Rin_all = []
    for i, ch in enumerate(active_chans):
        key = 'Ch' + str(ch)
        chan_data = charact_data[key][0]

        SS = np.mean(chan_data[indexes, 15_000:20_000], axis=1)
        Rin = (SS - rmp[i]) / np.array(vals_for_Rin) * 1000

        Rin_all.append(np.mean(Rin))
    return Rin_all

def mono_exp(x, a, tau, b):
    return a * np.exp(-x / tau) + b

def get_tau_FM(filename, active_chans):
    '''
    fit mono-exponential function to initial voltage response
    for -100 and -50 pA injections
    average the two tau values
    '''
    vals_for_tau = [-100, -50]
    charact_data = load_traces_FM(filename)
    inj = hcf.get_inj_current_steps(filename)
    indexes = np.array([inj.index(value) for value in vals_for_tau if value in inj])
    onset, _ = hcf.find_charact_onset_offset(filename)

    tau_all = []
    for ch in active_chans:
        key = 'Ch' + str(ch)
        chan_data = charact_data[key][0]

        taus = []
        for index in indexes:
            trace = chan_data[index, onset:onset+2000]
            x = np.arange(len(trace)) / 20  # ms, for 20 kHz sampling

            p0 = [trace[0] - trace[-1], 20, trace[-1]]
            popt, _ = curve_fit(mono_exp, x, trace, p0=p0)

            taus.append(popt[1])

        tau_all.append(np.mean(taus))

    return tau_all

def get_tau_when_artifact_FM(filename, active_chans):
    '''
    fit mono-exponential function to initial voltage response
    for -100 and -50 pA injections
    average the two tau values
    '''
    vals_for_tau = [-100, -50]
    charact_data = load_traces_FM(filename)
    inj = hcf.get_inj_current_steps(filename)
    indexes = np.array([inj.index(value) for value in vals_for_tau if value in inj])
    onset, _ = hcf.find_charact_onset_offset(filename)

    tau_all = []
    for ch in active_chans:
        key = 'Ch' + str(ch)
        chan_data = charact_data[key][0]

        taus = []
        for index in indexes:
            trace = chan_data[index, onset+5:onset+2000]
            x = np.arange(len(trace)) / 20  # ms, for 20 kHz sampling

            p0 = [trace[0] - trace[-1], 20, trace[-1]]
            popt, _ = curve_fit(mono_exp, x, trace, p0=p0)

            taus.append(popt[1])

        tau_all.append(np.mean(taus))

    return tau_all

def sag_func(x, a, b, c, tau1, tau2, tau3, d):
    return (np.exp(-(x-c) / tau1)
            - (a * np.exp(-(x-c) / tau2)
            + b * np.exp(-(x-c) / tau3))
            - d)

def get_sag_FM(filename, active_chans):
    '''
    approximate sag from -100 pA injection
    sag ratio = peak voltage amplitude / steady-state voltage amplitude
    both measured from RMP
    '''
    val_for_sag = -100
    charact_data = load_traces_FM(filename)
    inj = hcf.get_inj_current_steps(filename)
    index = inj.index(val_for_sag)

    rmp = get_RMP_FM(filename, active_chans)

    sag_all = []
    for i, ch in enumerate(active_chans):
        key = 'Ch' + str(ch)
        chan_data = charact_data[key][0]

        trace = chan_data[index, 5000:20_000]
        x = np.arange(len(trace)) / 20  # ms

        p0 = [1, 1, 0, 20, 100, 500, -trace[-1]]
        
        try:
            popt, _ = curve_fit(sag_func, x, trace, p0=p0, maxfev=10000)
        except RuntimeError:
            sag_all.append(np.nan)
            print(f'Sag fit failed: {filename}, Ch{ch}')
            continue

        fitted_trace = sag_func(x, *popt)

        peak = np.min(fitted_trace)
        SS = np.mean(fitted_trace[-2000:])

        sag_ratio = (peak - rmp[i]) / (SS - rmp[i])

        sag_all.append(sag_ratio)

    return sag_all

def get_hyperpol_properties(filename, active_chans, dict_cell_ID):
    rmp = get_RMP_FM(filename, active_chans)
    Rin = get_Rin_FM(filename, active_chans)
    tau = get_tau_FM(filename, active_chans)
    sag = get_sag_FM(filename, active_chans)

    rows = []

    for i, chan in enumerate(active_chans):
        rows.append({
            'channel': f'Ch{chan}',
            'cell_ID': dict_cell_ID[chan],
            'RMP': rmp[i],
            'Rin': Rin[i],
            'tau': tau[i],
            'sag': sag[i]
        })

    return pd.DataFrame(rows)

def get_AP_HW(AP1, THloc, peak_loc, AP_TH, AP_amp, sampling_rate):
    '''
    estimated tehb AP Hw
    '''
    half_amp = AP_TH + AP_amp / 2
    t30 = AP_TH + AP_amp * 0.3
    t70 = AP_TH + AP_amp * 0.7

    rise_start = np.where(AP1[THloc:peak_loc] > t30)[0][0] + THloc
    rise_end = np.where(AP1[THloc:peak_loc] > t70)[0][0] + THloc

    x_rise = np.arange(rise_start, rise_end)
    y_rise = AP1[rise_start:rise_end]

    model_rise = np.polyfit(x_rise, y_rise, 1)
    hw1 = (half_amp - model_rise[1]) / model_rise[0]

    fall_start = np.where(AP1[peak_loc:] > t70)[0][-1] + peak_loc
    fall_end = np.where(AP1[peak_loc:] > t30)[0][-1] + peak_loc

    x_fall = np.arange(fall_start, fall_end)
    y_fall = AP1[fall_start:fall_end]

    model_fall = np.polyfit(x_fall, y_fall, 1)
    hw2 = (half_amp - model_fall[1]) / model_fall[0]

    return (hw2 - hw1) / sampling_rate

class AP_params_FM:

    def __init__(self, filename, active_chans, dict_cell_ID):
        self.sampling_rate = 20 # kHz
        self.filename = filename[-12:]
        self.active_chans = active_chans
        self.data = load_traces_FM(filename)
        self.inj = hcf.get_inj_current_steps(filename)
        self.rheo = {}
        self.peaks = self.find_AP_peak()
        self.cell_IDs = dict_cell_ID

    def find_AP_peak(self):
        '''
        finds the peak of the first AP
        returns the trace and the corresponding peak location
        '''
        peaks = {}
        first_ap_peak_loc = 150

        for chan in self.active_chans:
            first_ap_trace = None
            chan_data = self.data[f'Ch{chan}'][0]

            for swp in range(np.shape(chan_data)[0]):
                trace = chan_data[swp]
                pks = detect_peaks(trace, mph=2, mpd=1)

                if first_ap_trace is None and len(pks) > 0:
                    first_ap_trace = trace[pks[0]-first_ap_peak_loc:pks[0]+first_ap_peak_loc]
                    self.rheo[f'Ch{chan}'] = self.inj[swp]

                    peaks[f'Ch{chan}'] = first_ap_trace, first_ap_peak_loc
                    break

            if first_ap_trace is None:
                peaks[f'Ch{chan}'] = [np.nan] * 4
                self.rheo[f'Ch{chan}'] = np.nan
                print(f'No APs found for {self.filename} chan {chan}')

        return peaks

    def get_AP_params(self):
        '''
        defines the AP properties
        returns a dictionary with channels as keys
        and for eahc cha,t here is a dictionary with all AP params
        '''

        AP_keys = ['AP_TH', 'AP_up', 'AP_down', 'AP_amp',
           'AP_HW', 'AHP', 'AP_AHP_slope']
        AP_props = {}

        for chan in self.active_chans:
            AP1 = self.peaks[f'Ch{chan}'][0]
            peak_loc = self.peaks[f'Ch{chan}'][1]
            if len(AP1) > 0:
                d1_AP = np.diff(AP1)*20
                THlocs = np.where(d1_AP[:peak_loc-1] > 20)
                if THlocs[0].size != 0:
                    THloc = THlocs[0][0]
                    AP_TH = AP1[THloc-1]
                    AP_up = np.max(d1_AP[:peak_loc-1])
                    AP_down = np.min(d1_AP[peak_loc:])
                    AP_amp = AP1[peak_loc] - AP_TH

                    AP_HW = get_AP_HW(AP1, THloc, peak_loc, AP_TH, AP_amp, self.sampling_rate)

                    # if no min If no minimum was reached in a 10.0-ms time window,
                    #the amplitude at the point after the AP, where the slope exceeded −0.2 mV ms−1
                    min_ = np.min(AP1[peak_loc:])
                    AHP = AP_TH - min_

                    AHP_loc = peak_loc + np.argmin(AP1[peak_loc:]) + 20
                    AP_AHP_loc = AP1[AHP_loc:AHP_loc + 20]
                    AP_AHP_slope = np.mean(
                        np.diff(AP_AHP_loc) * self.sampling_rate)

                    AP_props[f'Ch{chan}'] = dict(zip(
                        AP_keys,
                        [AP_TH, AP_up, AP_down, AP_amp, AP_HW, AHP, AP_AHP_slope]))

                else:
                    AP_props[f'Ch{chan}'] = dict.fromkeys(AP_keys,np.nan)
            else:
                AP_props[f'Ch{chan}'] = dict.fromkeys(AP_keys,np.nan)

        return AP_props

    def get_firing_params(self):
        '''
        extracting rheobase, firing frequency at 700 pA, F/I slope
        and ISI adaptation
        '''

        firing_keys = ['rheobase', 'FF_700pA', 'F_I_slope', 'ISI_adaptation']
        firing_props = {}
        FF_700, ISI_adaptation = np.nan, np.nan

        for chan in self.active_chans:
            chan_data = self.data[f'Ch{chan}'][0]
            swp_num = np.shape(self.data[f'Ch{chan}'][0])[0]

            indx_0 = self.inj.index(0)
            indx_700 = self.inj.index(700)

            firing_freq = []
            count = 0
            ISI_calc = []
            for swp in range(indx_0, swp_num):

                trace = chan_data[swp]
                pks = detect_peaks(trace, mph=2, mpd=1)
                firing_freq.append(len(pks)) # becasuse 1000 ms stim window

                if swp == indx_700:
                    FF_700 = len(pks)

                if len(pks) >= 3 and count <= 5:
                    count += 1
                    ISI_1 = pks[1] - pks[0]
                    ISI_2 = pks[2] - pks[1]
                    ISI_calc.append(ISI_2 / ISI_1)

                if count == 5 or swp == len(self.inj) - 1:
                    ISI_adaptation = np.mean(ISI_calc) if len(ISI_calc) > 0 else np.nan

            if len(self.inj[indx_0:-1]) != len(firing_freq):
                print(f'Fix injection or fn for {self.filename}')
                FI_slope = np.nan
            else:
                FI_slope = np.polyfit(self.inj[indx_0:-1], firing_freq, 1)[0]

            firing_props[f'Ch{chan}'] = dict(zip(
                firing_keys,
                [self.rheo[f'Ch{chan}'], FF_700, FI_slope, ISI_adaptation]))

        return firing_props

    def get_dataframe(self):
        '''
        pack everything into a dataframe
        '''
        first_ap_params = self.get_AP_params()
        firing_params = self.get_firing_params()

        rows = []

        for chan in self.active_chans:
            rows.append({
                'filename': self.filename,
                'channel': f'Ch{chan}',
                'cell_ID': self.cell_IDs[chan],
                **first_ap_params[f'Ch{chan}'],
                **firing_params[f'Ch{chan}']
            })

        return pd.DataFrame(rows)


# if __name__ == '__main__':
#     filename = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/data_verji/OP230314/23314014.abf'
#     active_chans = [2]
#     AP_data = AP_params_FM(filename, active_chans)
#     df = AP_data.get_dataframe()

# TO GET analysis as FM of intr_poperties

def get_data_analyzed_as_FM(df_intr_VM,
                            human_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    '''
    for every filename of the VM ff, creates the equavalent data analyzed in FM way
    merges and keeps metadata columns
    '''
    dfs = []
    for fn in df_intr_VM.filename.unique():
        print(f'starting analysis of {fn}')
        fn_df = df_intr_VM[df_intr_VM.filename == fn]
        OP = fn_df.OP.values[0]
        patcher = fn_df.patcher.values[0]
        work_dir = sort.get_work_dir(human_dir, OP, patcher)
        active_chans = fn_df.cell_ch.to_list()
        dict_cell_ID =  dict(zip(fn_df['cell_ch'], fn_df['cell_IDs_match']))


        AP_data = AP_params_FM(filename=work_dir + fn, active_chans=active_chans, dict_cell_ID=dict_cell_ID)
        df_depol = AP_data.get_dataframe()

        df_hyperpol = get_hyperpol_properties(
            work_dir + fn,
            active_chans,
            dict_cell_ID)

        df = df_depol.merge(
            df_hyperpol,
            on=['channel', 'cell_ID'],
            how='left')

        dfs.append(df)

    df_AP = pd.concat(dfs, ignore_index=True)
    df_AP = df_AP.rename(columns={'cell_ID': 'cell_IDs_match'})

    if len(df_AP) == len(df_intr_VM):

        # merge dataframes
        cols_intr = [
        'cell_IDs_match',
        'resting_potential',
        'area',
        'tissue_source',
        'OP',
        'patcher',
        'patient_age',
        'filename',
        'slice',
        'cell_ch',
        'cell_ID',
        'day',
        'treatment',
        'hrs_incubation',
        'repatch',
        'hrs_after_OP',
        'dataset'
        ]

        df_intr_as_FM = df_AP.merge(
            df_intr_VM[cols_intr],
            on='cell_IDs_match',
            how='left'
        )
        
        # save result
        save_dir = '/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/results/human/connectivity_paper/data/'
        # df_intr_as_FM.to_csv(f'{save_dir}VM_data_analysed_as_FM.csv')
        return df_intr_as_FM


def get_fixed_tau(df, human_dir='/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/data/human/'):
    '''
    trying to compensate for the artifact
    taking the 5th point after onset
    '''
    df = df.copy()
    df['tau_fixed'] = df['tau']
    df_big_tau = df[df['tau'] > 40].reset_index()

    taus = []
    for i, fn in enumerate(df_big_tau.filename_x):
        chan = df_big_tau.cell_ch[i]
        OP = df_big_tau.OP[i]
        patcher = df_big_tau.patcher[i]

        work_dir = sort.get_work_dir(human_dir, OP, patcher)
        tau = get_tau_when_artifact_FM(work_dir + fn, [chan])

        taus.append(tau[0])

    df.loc[df_big_tau['index'], 'tau_fixed'] = taus
    return df

def apply_filter_properties_mask(df):
    '''
    mask to apply on VM cells to match as close as possible to FM
    converts AP down to positive
    combines the rheobase data from ramp and inj. renames to rheobase
    '''
    df = df.copy()
    if 'tau_fixed' not in df.columns:
        mask_exclude = (df['resting_potential'] > -60) | \
                        (df['tau'] > 40) | (df['tau'] < 0) |\
                        (df['AP_HW'] > 2) | (df['AP_HW'] < 0) | \
                        (df['sag'] < 1) | (df['sag'] > 1.5) | \
                        (df['rheobase'] < 100) | \
                        (df['F_I_slope'] > 0.3)
    else:
        mask_exclude = (df['resting_potential'] > -60) | \
                                (df['tau_fixed'] > 40) | (df['tau_fixed'] < 0) |\
                                (df['AP_HW'] > 2) | (df['AP_HW'] < 0) |\
                                (df['sag'] < 1) | (df['sag'] > 1.5) | \
                                (df['rheobase'] < 100) |\
                                (df['F_I_slope'] > 0.3)

    df = df.loc[~mask_exclude].copy()
    df.loc[:, 'AP_down'] = df['AP_down'].abs()
    # combine the rheobase columns
    df['rheobase_combined'] = (df['rheo_ramp_c'].fillna(df['rheobase']))
    print(f'remaining cell #: {len(df)}')
    return df


def get_dict_VM_FM():
    dict_intrinsic_VM_FM = {'resting_potential': 'V_m',
                        'Rin': 'R_in',
                        'tau_fixed': 'tau',
                        'AP_TH': 'AP_tamp',
                        'AP_amp': 'AP_amp',
                        'AP_HW': 'AP_hw',
                        'AP_up': 'AP_up',
                        'AP_down': 'AP_down',
                        'FF_700pA': 'freq_700',
                        'sag': 'sag',
                        'rheobase': 'rheo',
                        'rheobase_combined': 'rheo',
                        'AHP': 'AP_ahp',
                        'AP_AHP_slope': 'AP_adp',
                        'ISI_adaptation': 'isi_adapt',
                        'F_I_slope': 'f_i_slope'
                    }
    return dict_intrinsic_VM_FM

def remove_nans_from_dict(df):
    '''
    removes nans for KNN analysis
    '''
    df = df.copy()
    params_of_interest = [
    col for col in get_dict_VM_FM().keys()
    if col in df.columns]

    df = df.dropna(subset=params_of_interest).reset_index(drop=True)
    print(f'remaining cell #: {len(df)}')
    return df

def knn_analysis(dict_intrinsic_VM_FM, df_VM, df_paper):
    '''
    Standardize electrophysiological features, classify cells using 5-NN,
    and project paper and experimental data into the same 2D PCA space.
    '''
    # define the features
    if 'rheobase' in dict_intrinsic_VM_FM.keys():
        del dict_intrinsic_VM_FM['rheobase'] # looking at 'rheo_ramp_c' instead
    for key, val in dict_intrinsic_VM_FM.items():
        df_VM[val] =  df_VM[key]

    # creating an empty standardization tool
    scaler = StandardScaler()
    features = list(dict_intrinsic_VM_FM.values())

    # papaer mean and SD from papaer in in scaler
    X_paper = scaler.fit_transform(df_paper[features])
    # fitting the scaler to my data. so that they are in the same 'space'
    X_my = scaler.transform(df_VM[features])

    # train 5-nearest-neighbour classifier
    knn = KNeighborsClassifier(n_neighbors=5)
    knn.fit(X_paper, df_paper['etype_idx'])

    # classify your cells. predict their e-type
    # calcualtes distance --> finds 5 closes --> assigns majority e-type 
    df_VM['e_type'] = knn.predict(X_my)

    # plot the result

    # PCA: define space from the paper data
    pca = PCA(n_components=2)

    paper_pca = pca.fit_transform(X_paper)
    my_pca = pca.transform(X_my)

    return paper_pca, my_pca, df_VM.copy()

def create_df_plot(paper_pca, df_paper, my_pca, df_VM):
    '''
    create a combined dataframe for plotting
    '''
    df_paper_pca = pd.DataFrame(
        paper_pca,
        columns=['PC1', 'PC2']
    )

    df_paper_pca['e_type'] = df_paper['etype_idx'].values
    df_paper_pca['dataset'] = 'paper'

    df_my_pca = pd.DataFrame(
        my_pca,
        columns=['PC1', 'PC2']
    )

    df_my_pca['e_type'] = df_VM['e_type'].values
    df_my_pca['dataset'] = 'my data'

    df_plot = pd.concat(
        [df_paper_pca, df_my_pca],
        ignore_index=True
    )
    return df_plot


# FOR IMPUTATION
def prepare_imputation_data(df, df_paper, feature_dict):
    '''
    make AP down positive
    ocmbine rheobase
    rename columns
    include column of missing n_nums   
    '''

    df = df.copy()
    df['AP_down'] = df['AP_down'].abs()

    df['rheobase_combined'] = (
        df['rheo_ramp_c'].fillna(df['rheobase'])
    )

    features_my = list(feature_dict.keys())
    features_paper = list(feature_dict.values())

    X_my = df[features_my].rename(columns=feature_dict)
    X_paper = df_paper[features_paper]

    assert X_my.columns.equals(X_paper.columns)

    df['n_missing'] = X_my.isna().sum(axis=1)

    return df, X_my, X_paper

def fit_knn_model(X_paper, e_types, n_neighbors=5):
    """
    Standardize the reference ddataframe (df_FM) and 
    fit a k-nearest-neighbor classifier (n=5).

    Returns the fitted scaler and KNN model.
    """

    scaler = StandardScaler()
    X_paper_scaled = scaler.fit_transform(X_paper)

    knn = KNeighborsClassifier(n_neighbors=n_neighbors)
    knn.fit(X_paper_scaled, e_types)

    return scaler, knn

def fit_imputer(X_reference, scaler, n_neighbors=5):
    """
    Fit a KNN imputer on standardized reference-cell features.

    Returns the fitted imputer.
    """

    X_reference_scaled = scaler.transform(X_reference)

    imputer = KNNImputer(n_neighbors=n_neighbors)
    imputer.fit(X_reference_scaled)

    return imputer

def classify_cells(X_my, scaler, imputer, knn):
    """
    Impute missing features and classify cells using the fitted KNN model.

    Returns predicted e-types.
    """

    X_scaled = scaler.transform(X_my)
    X_imputed = imputer.transform(X_scaled)

    assert not np.isnan(X_imputed).any()

    return knn.predict(X_imputed)


def check_classification(df, df_e_types, pred_col):

    """
    BOT CREATED
    Compare predicted e-types with previously assigned e-types for complete cells.

    Returns matched cells, classification agreement, and a crosstab.
    """

    if pred_col not in df.columns:
        raise ValueError(f'{pred_col} not found in dataframe')

    df_complete = df[df['n_missing'] == 0].copy()

    df_complete = df_complete.merge(
        df_e_types[['cell_IDs_match', 'e_type']],
        on='cell_IDs_match',
        how='left'
    )

    df_compare = df_complete.dropna(subset=['e_type']).copy()

    agreement = (
        df_compare['e_type'].astype(int) ==
        df_compare[pred_col].astype(int)
    ).mean()

    table = pd.crosstab(
        df_compare['e_type'],
        df_compare[pred_col]
    )

    return df_compare, agreement, table

def compare_imputers(
    df_compare,
    features_my,
    feature_dict,
    scaler,
    imputer_paper,
    imputer_own,
    knn):
    """
    BOT CREATED
    Compare paper-based and own-data imputation by removing one feature at a time.

    Returns classification agreement for each removed feature.
    """
    results = []

    for feature in features_my:

        df_test = df_compare.copy()
        df_test[feature] = np.nan

        X_test = df_test[features_my].rename(columns=feature_dict)

        e_type_paper = classify_cells(
            X_test, scaler, imputer_paper, knn
        )

        e_type_own = classify_cells(
            X_test, scaler, imputer_own, knn
        )

        e_type_true = df_test['e_type'].astype(int).values

        results.append({
            'feature_removed': feature,
            'paper_imputer': (e_type_paper == e_type_true).mean(),
            'own_imputer': (e_type_own == e_type_true).mean()
        })

    return pd.DataFrame(results)

def impute_dataframe(df, X_my, feature_dict, scaler, imputer):
    """
    BOT CREATED
    Impute missing feature values and add them back to a copy of the dataframe.

    Returns the completed dataframe and standardized imputed feature matrix.
    """

    X_scaled = scaler.transform(X_my)
    X_imputed_scaled = imputer.transform(X_scaled)

    # back to original units
    X_imputed = scaler.inverse_transform(X_imputed_scaled)

    X_imputed = pd.DataFrame(
        X_imputed,
        index=X_my.index,
        columns=X_my.columns
    )

    # paper names -> your column names
    reverse_dict = {val: key for key, val in feature_dict.items()}
    X_imputed = X_imputed.rename(columns=reverse_dict)

    df_imputed = df.copy()
    features_my = list(feature_dict.keys())

    # only replace missing values; keep measured values untouched
    df_imputed[features_my] = (
        df_imputed[features_my]
        .fillna(X_imputed[features_my])
    )

    return df_imputed, X_imputed_scaled

def get_pca_projection(X_paper, X_my, scaler, imputer):
    """
    Fit PCA on standardized paper data and project imputed experimental cells into the same space.

    Returns paper PCA coordinates, experimental PCA coordinates, and the fitted PCA model.
    """

    X_paper_scaled = scaler.transform(X_paper)

    X_my_scaled = scaler.transform(X_my)
    X_my_imputed = imputer.transform(X_my_scaled)

    pca = PCA(n_components=2)

    paper_pca = pca.fit_transform(X_paper_scaled)
    my_pca = pca.transform(X_my_imputed)

    return paper_pca, my_pca

def compare_etype_distributions(
    df_original,
    df_complete,
    original_col='etype_idx',
    complete_col='e_type'):
    """
    Compare e-type counts and percentages across three datasets.
    """

    dfs = {
        'original': df_original[original_col],
        'own_complete': df_complete[complete_col]
        # 'imputed': df_imputed[imputed_col]
    }

    results = {}

    for name, data in dfs.items():
        results[f'{name}_n'] = data.value_counts()
        results[f'{name}_%'] = data.value_counts(normalize=True) * 100

    return pd.DataFrame(results).sort_index().round(2)

def impute_on_own_data_assign_e_types(df_VM, df_paper, dict_VM_FM):
    '''
    main impute function
    '''
    df_VM = df_VM.copy()
    df_paper = df_paper.copy()

    df_VM = df_VM.drop('tau', axis=1)
    dict_VM_FM.pop('rheobase')

    intr_all_impute, features_my, features_paper = prepare_imputation_data(
        df_VM, df_paper, dict_VM_FM)

    # create paper-based feature space
    scaler, knn = fit_knn_model(features_paper, df_paper['etype_idx'])

    # OWN-DATA IMPUTER
    X_my_complete = features_my.loc[intr_all_impute['n_missing'] == 0]

    # using imputer on my data because it showed better results
    imputer_own = fit_imputer(X_my_complete, scaler)

    intr_all_impute['e_type_own_imp'] = classify_cells(features_my, scaler, imputer_own,knn)

    # PCA space
    paper_pca, my_pca = get_pca_projection(features_paper, features_my, scaler, imputer_own)

    # prepare plot dataframe
    df_VM_plot = intr_all_impute.copy()
    df_VM_plot['e_type'] = df_VM_plot['e_type_own_imp'].astype(int)

    df_paper_VM = create_df_plot(paper_pca, df_paper, my_pca, df_VM_plot)

    return df_paper_VM, df_VM_plot


# FOR CONNECTION TYPES

def get_slice_connection_summary(df_cells, df_connections):
    rows = []

    for group, df_slice in df_cells.groupby(['OP', 'patcher', 'slice', 'day']):

        OP, patcher, slice_, day = group

        df_con_slice = df_connections[
            (df_connections['OP'] == OP) &
            (df_connections['patcher'] == patcher) &
            (df_connections['slice'] == slice_) &
            (df_connections['day'] == day)
        ]

        n_cells = len(df_slice)
        n_connections = len(df_con_slice)

        n_possible = n_cells * (n_cells - 1)

        connection_probability = (
            n_connections / n_possible * 100
            if n_possible > 0
            else 0
        )

        rows.append({
            'OP': OP,
            'patcher': patcher,
            'day' : day,
            'slice': slice_,
            'treatment' : df_slice['treatment'].values[0],
            'num_channels': n_cells,
            'found_connections': n_connections,
            'num_possible_connections': n_possible,
            'connection_probability': connection_probability
        })

    return pd.DataFrame(rows)


def get_connection_lookup(ETYPES = [2, 3, 4, 5]):

    connection_lookup = {
        (pre, post): i
        for i, (pre, post) in enumerate(product(ETYPES, repeat=2))
    }
    return connection_lookup

def get_all_possible_connection_from_intr(df_intr, df_cons):
    '''
    for each slices counts the number of possible connections
    and assigns a synapse types
    e_type column must be present in dataframe
    then assigns 1 if connection actually exists
    returns a dataframe with possible and found connections
    '''

    if 'e_type' not in df_intr.columns:
        print('No e_type column found in ointrinsic data')
        return

    con_pairs = set(zip(df_cons['cell_IDs_match_pre'], df_cons['cell_IDs_match_post']))
    etype_lookup = (df_intr.set_index('cell_IDs_match')['e_type'])
    connection_lookup = get_connection_lookup()

    all_pairs = []
    for group, df_slice in df_intr.groupby(['OP', 'patcher', 'day', 'slice', 'treatment']):
        OP, patcher, day, slic, tr = group

        cell_ids = df_slice['cell_IDs_match'].drop_duplicates()

        pairs = pd.DataFrame(
            permutations(cell_ids, 2),
            columns=['cell_IDs_match_pre', 'cell_IDs_match_post'])

        pairs['OP'] = OP
        pairs['patcher'] = patcher
        pairs['slice'] = slic
        pairs['day'] = day
        pairs['treatment'] = tr

        pairs['connected'] = [(pre, post) in con_pairs
                for pre, post in zip(pairs['cell_IDs_match_pre'],pairs['cell_IDs_match_post'])]

        pairs['pre_etype'] = pairs['cell_IDs_match_pre'].map(etype_lookup)
        pairs['post_etype'] = pairs['cell_IDs_match_post'].map(etype_lookup)

        pairs['connection_type'] = [connection_lookup[(pre, post)]
                    for pre, post in zip(pairs['pre_etype'], pairs['post_etype'])]

        all_pairs.append(pairs)

    df_pairs = pd.concat(all_pairs, ignore_index=True)

    return df_pairs

def get_first_param_matrix(df_event_props, df_FM, param, ETYPES):
    '''
    creates a param matrix based on the
    provided e_types

    shows hte average param per connection type
    '''
    stats = df_event_props.groupby(['e_type_pre', 'e_type_post'])[param].agg(['mean', 'count'])

    param_matrix = stats['mean'].unstack().reindex(index=ETYPES, columns=ETYPES)
    n_nums = stats['count'].unstack().reindex(index=ETYPES, columns=ETYPES).fillna(0).astype(int)
    
    etype_names = df_FM[['etype_idx', 'etype']].drop_duplicates().set_index('etype_idx')['etype'].to_dict()

    param_matrix = param_matrix.rename(index=etype_names, columns=etype_names)

    return param_matrix, n_nums


def get_connection_matrix(df_pairs, df_FM,  ETYPES = [2, 3, 4, 5]):
    '''
    creates a matrix based on the
    pre and post etypes

    shows connection probability and
    found/tested connections
    '''
    connection_lookup = get_connection_lookup(ETYPES)
    connection_matrix = pd.DataFrame(index=ETYPES, columns=ETYPES, dtype=object)

    n_found = dict.fromkeys(connection_lookup.values(), 0)
    n_tested = dict.fromkeys(connection_lookup.values(), 0)

    for _, row in df_pairs.iterrows():

        pre = int(row['pre_etype'])
        post = int(row['post_etype'])

        if (pre, post) in connection_lookup:
            con_type = connection_lookup[(pre, post)]
            n_tested[con_type] += 1

            if row['connected']:
                n_found[con_type] += 1

    connection_matrix = pd.DataFrame(index=ETYPES, columns=ETYPES, dtype=object)

    for pre in ETYPES:
        for post in ETYPES:

            con_type = connection_lookup[(pre, post)]

            found = n_found[con_type]
            tested = n_tested[con_type]

            if tested > 0:
                probability = found / tested * 100
                connection_matrix.loc[pre, post] = (f'{probability:.1f}% ({found}/{tested})')
            else:
                connection_matrix.loc[pre, post] = 'NA (0/0)'

    etype_names = (df_FM[['etype_idx', 'etype']].drop_duplicates().
                   set_index('etype_idx')['etype'].to_dict())

    connection_matrix = connection_matrix.rename(index=etype_names, columns=etype_names)

    connection_matrix.index.name = 'pre_etype'
    connection_matrix.columns.name = 'post_etype'

    return connection_matrix


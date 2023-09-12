import os
import numpy as np
import pandas as pd
import warnings
import pickle
import matplotlib.pyplot as plt
from pathlib import Path
from functools import reduce
from tqdm.auto import tqdm

from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

from util import MatfileIO, Bunch
from preprocessing import GaussianSmoothing


def preprocess_unitTimeBin(unitTimeBin, bin_id, g_sigma, g_length, g_decimal_out):
    """
    Preprocessing of the unit-by-timeBin neural data including smoothing.
    :param unitTimeBin:
    :param bin_id: timeBins (columns) to be included
    :param g_sigma: sigma of the gaussian kernel for smoothing
    :param g_length: length of the gaussian kernel for smoothing
    :param g_decimal_out: decimal points after convolution
    :return: ndarray with same dimensions as input after smoothing
    """
    utb = unitTimeBin[:, bin_id]
    sm_obj = GaussianSmoothing(utb, sigma=g_sigma, length=g_length, axis='row', decimal_out=g_decimal_out)
    conv = sm_obj.conv()
    return  conv  # conv.flatten()

def bin_step_size_win(step, width, window): 
    bin_start, bin_end = [], []
    start = 0
    while (window[0] + start + width) <= window[1]: 
        bin_start.append(window[0] + start)
        bin_end.append(window[0] + start + width)
        start += step    
    return np.array(bin_start)

def get_params():
    # Preprocessing - params
    p = Bunch()
    p.rerun = True
    p.bin_step = 0.25
    p.bin_width = 0.5
    p.time_win = [-1, 2]
    p.time_bins = bin_step_size_win(p.bin_step, p.bin_width, p.time_win)
    p.time_clip0 = p.time_win[0]  # no clip by default
    p.time_clip1 = p.time_win[1]  # no clip by default
    p.g_sigma = 3
    p.g_length = 15
    p.g_decimal_out = 3

    p.bl_I = Bunch()
    p.bl_I[1], p.bl_I[5] = 'lelo', 'lelo'  # left/low
    p.bl_I[2], p.bl_I[6] = 'lehi', 'lehi'  # left/high
    p.bl_I[3], p.bl_I[7] = 'rilo', 'rilo'  # right/low
    p.bl_I[4], p.bl_I[8] = 'rihi', 'rihi'  # right/high

    # params for classifier
    p.cf = Bunch()
    p.cf.resample_iter = 100
    p.cf.n = 1000  # n estimators
    p.cf.d = 100  # max depth
    p.cf.n_cell = 20  # number of cells to be resampled
    p.cf.to_scale = True
    return p

def train_classifier(X, y, classifier_type, resample_n_neuron=True,
                    resample_iter=100, n_tree=1000, d_tree=100, n_cell=20, cv_fold=5):
    # train_random_forest_classifier(X, y, resample_n_neuron=None, resample_iter=100):
    """
    :param X: 3-d array with neuron x timeBin x trial orientation
    :param y: numpy array with class labels
    :param resample_n_neuron: boolean for resampling
    :param resample_iter: the # of iterations to be run
    :param n_tree: n_estimators (the # of trees)
    :param d_tree: max depth
    :return: a dictionary containing classification report (mean, max, min)
    """
    X = np.array(X)
    y = np.array(y)
    assert X.shape[2] == y.shape[0]
   
    # Define the classifier
    if classifier_type.lower().find("r") != -1:
        clf = RandomForestClassifier(n_estimators=n_tree, max_depth=d_tree, random_state=42)
    elif classifier_type.lower().find("x") != -1:
        clf = xgb.XGBClassifier(n_estimators=n_tree, max_depth=d_tree, random_state=42)
    
    df_reports = []

    for i in tqdm(range(resample_iter)):
        if resample_n_neuron:
            X_idx = np.random.choice(X.shape[0], size=n_cell, replace=True)
            X_rs = X[X_idx, :, :]  # resample
        else:
            X_rs = X

        X_flat = np.reshape(X_rs, (1, -1, X_rs.shape[-1]))  # flatten
        X_in = np.transpose(np.squeeze(X_flat))  # transpose to trial by (neuron x time) 2d array
        assert X_in.shape[0] == y.shape[0]  # match the number of total trials

        skf = StratifiedKFold(n_splits=5)
        for train_idx, test_idx in skf.split(X_in, y):
            clf.fit(X_in[train_idx, :], y[train_idx])  # fit on test set
            y_pred = clf.predict(X_in[test_idx, :])
            print(np.around(clf.score(X_in[test_idx, :], y[test_idx]), 4))
            report = classification_report(y[test_idx], y_pred, output_dict=True, zero_division=0)
            df_reports.append(pd.DataFrame(report).reset_index())

    merged_df = reduce(lambda l, r: pd.concat([l, r]), df_reports)  # concatenate into a df  
    rez = Bunch()
    rez.mean = merged_df.groupby('index').agg('mean')  
    rez.med = merged_df.groupby('index').agg('median')
    rez.max = merged_df.groupby('index').agg('max')
    rez.min = merged_df.groupby('index').agg('min')
    return rez

def df_to_3d_arr(df, col):
    l = []
    for _, row in df.iterrows():
        l.append(row[col])
    return np.transpose(np.array(l), [1, 2, 0])  # transpose to neuron by time-bin by trial

def standard_scaler(neuron_time_trial):
    """
    :param neuron_time_trial: a 3d np array with neuron x time x trial orientation
    :return: a scaled array with same orientation/dimensions 
    """
    shp = neuron_time_trial.shape
    rs = neuron_time_trial.reshape(shp[0], shp[1] * shp[2])
    rs_t = rs.transpose([1, 0])  # reshape to sample by feature 
    
    scaler = StandardScaler()
    scaler.fit(rs_t)
    scaled = scaler.transform(rs_t)
    scaled_neuron_time_trial = scaled.transpose([1, 0]).reshape(shp[0], shp[1], shp[2])
    return scaled_neuron_time_trial

def to_run(dat_paths: list, param_setter={}):
    """
    :param dat_paths: list of filePaths (a string input will be converted to a list)  
    :param param_setter: a dictionary to overwrite specific parameters, e.g., param_setter={'resample_iter':10} to overwrite the 'resample_iter'. 
    :return: a dictionary containing classification report (mean, max, min)
    """
    _p = get_params()  # default params
    
    for key in param_setter.keys():  # update, if any, from param_setter
        if key in _p.keys():
            _p[key] = param_setter[key]
        if key in _p.cf.keys(): 
            _p.cf[key] = param_setter[key]
    
    if isinstance(dat_paths, str):
        dat_paths = list(dat_paths.split(" "))  # if a string input, convert to list
    
    for i in range(len(dat_paths)):  # iterate over paths
        if os.path.exists(dat_paths[i]): 
            out_dir = os.path.join(dat_paths[i], 'classifier_output')
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)     
                rez, _params = run(Path(dat_paths[i]), _p)
                with open(os.path.join(out_dir, 'trial_classifier_result_pre.pickle'), 'wb') as handle:
                    pickle.dump([rez, _params], handle)
                # with open(os.path.join('/Volumes/dudmanlab/junchol/js2p0/WR40_081919/Matfiles/classifier_output', 'trial_classifier_result.pickle'), 'rb') as handle:
                #    rez, p = pickle.load(handle)
            elif _p.rerun == True:  # run
                rez, _params = run(Path(dat_paths[i]), _p)     
                with open(os.path.join(out_dir, 'trial_classifier_result_pre.pickle'), 'wb') as handle:
                    pickle.dump([rez, _params], handle)
        else: 
            warnings.warn("The data directory does not exist!")
        
def run(dat_path, p):
    """
    :param dat_path: path to the data directory as a pathlib.PosixPath 
    :param p: relevant parameters and hyperparameters for preprocessing and training classifiers
    :return: a Bunch that contains all classifier performance metrics
    """
    # Load data with MatfileIO class
    io_cls = MatfileIO(Path(dat_path))

    df = io_cls.extract_dataframe()
    #df['blockNum'] = df.blockNum.apply(lambda x: str(x))

    # Use successful trials only
    success_I = df.trialType.str.contains("sp", case=False, na=False)
    df_s = df[success_I]
    """
    - df rows correspond to trials. 
    - df columns correspond to different variables. 
        - blockNum: block IDs. 
        - trialType: trial types.  
            - 'to': timeout. 
            - 'sp': successful pull. 
            - 'ps': push.  
            - 'pmpp': premature pull & push. 
    - df.ctx, df.str, df.cg contain binned spike counts per trial.
        - Each np.matrix is organized as # neurons by # time bins. 
        - By default, each np.matrix spans 3 s epoch (1 s pre and 2 s post event to which time bins are aligned) with the bin width of 50 ms (60 bins).
    """
    # define block type
    bl_type = list(df_s.blockNum.apply(lambda x: p.bl_I[x]))
    df_s.insert(0, 'bl_type', bl_type)
    y = df_s.bl_type.to_numpy()  # labels (block types)
    y_dir = np.frompyfunc(lambda a:a[:2],1,1)(y)
    y_trq = np.frompyfunc(lambda a:a[-2:],1,1)(y)    
    
    p.bin_I = (p.time_bins >= p.time_clip0) & (p.time_bins <= p.time_clip1) # to include pre-reach epoch only
    
    num_cells = []
    Xs = []
    if 'ctx' in df_s.keys():
        num_cells.append(df_s.ctx.iloc[0].shape[0])  # get the number of neurons
        ctx_proc = df_s.ctx.apply(lambda x:
                                preprocess_unitTimeBin(x, p.bin_I, p.g_sigma, p.g_length, p.g_decimal_out))
        df_s.insert(df_s.shape[1], 'ctx_proc', ctx_proc)
        X_ctx = df_to_3d_arr(df_s, 'ctx')  # ndarray neuron x time bin x trial
        if p.cf.to_scale:
            X_ctx = standard_scaler(X_ctx)
        Xs.append(X_ctx)
        
    if 'str' in df_s.keys():
        num_cells.append(df_s.str.iloc[0].shape[0]) 
        str_proc = df_s.str.apply(lambda x:
                                preprocess_unitTimeBin(x, p.bin_I, p.g_sigma, p.g_length, p.g_decimal_out))
        df_s.insert(df_s.shape[1], 'str_proc', str_proc)
        X_str = df_to_3d_arr(df_s, 'str')  # ndarray neuron x time bin x trial
        if p.cf.to_scale:
            X_str = standard_scaler(X_str)
        Xs.append(X_str)
        
    if 'cg' in df_s.keys():
        num_cells.append(df_s.cg.iloc[0].shape[0])
        cg_proc = df_s.cg.apply(lambda x:
                              preprocess_unitTimeBin(x, p.bin_I, p.g_sigma, p.g_length, p.g_decimal_out))
        df_s.insert(df_s.shape[1], 'cg_proc', cg_proc)
        X_cg = df_to_3d_arr(df_s, 'cg')  # ndarray neuron x time bin x trial
        if p.cf.to_scale:
            X_cg = standard_scaler(X_cg)
        Xs.append(X_cg)
        
    X_all = np.concatenate((Xs), axis=0)  # all units
    if p.cf.to_scale: 
        X_all = standard_scaler(X_all)
    resample_n = min(min(num_cells), p.cf.n_cell)
    rez = Bunch()  # container for the result    

    if 'ctx' in df_s.keys():
        # with resampling to match neuron numbers across regions
        rez.ctx_rs_block = train_classifier(X_ctx, y, 'random forest', resample_n_neuron=True, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d, n_cell=resample_n)
        rez.ctx_rs_dir = train_classifier(X_ctx, y_dir, 'random forest', resample_n_neuron=True, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d, n_cell=resample_n)
        rez.ctx_rs_trq = train_classifier(X_ctx, y_trq, 'random forest', resample_n_neuron=True, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d, n_cell=resample_n)
        # without resampling
        rez.ctx_block = train_classifier(X_ctx, y, 'random forest', resample_n_neuron=False, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d)
        rez.ctx_dir = train_classifier(X_ctx, y_dir, 'random forest', resample_n_neuron=False, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d)
        rez.ctx_trq = train_classifier(X_ctx, y_trq, 'random forest', resample_n_neuron=False, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d)
        
    if 'str' in df_s.keys():
        # with resampling to match neuron numbers across regions
        rez.str_rs_block = train_classifier(X_str, y, 'random forest', resample_n_neuron=True, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d, n_cell=resample_n)
        rez.str_rs_dir = train_classifier(X_str, y_dir, 'random forest', resample_n_neuron=True, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d, n_cell=resample_n)
        rez.str_rs_trq = train_classifier(X_str, y_trq, 'random forest', resample_n_neuron=True, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d, n_cell=resample_n)
        # without resampling
        rez.str_block = train_classifier(X_str, y, 'random forest', resample_n_neuron=False, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d)
        rez.str_dir = train_classifier(X_str, y_dir, 'random forest', resample_n_neuron=False, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d)
        rez.str_trq = train_classifier(X_str, y_trq, 'random forest', resample_n_neuron=False, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d)
        
    if 'cg' in df_s.keys():
        # with resampling to match neuron numbers across regions
        rez.cg_rs_block = train_classifier(X_cg, y, 'random forest', resample_n_neuron=True, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d, n_cell=resample_n)
        rez.cg_rs_dir = train_classifier(X_cg, y_dir, 'random forest', resample_n_neuron=True, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d, n_cell=resample_n)
        rez.cg_rs_trq = train_classifier(X_cg, y_trq, 'random forest', resample_n_neuron=True, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d, n_cell=resample_n)
        # without resampling
        rez.cg_block = train_classifier(X_cg, y, 'random forest', resample_n_neuron=False, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d)
        rez.cg_dir = train_classifier(X_cg, y_dir, 'random forest', resample_n_neuron=False, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d)
        rez.cg_trq = train_classifier(X_cg, y_trq, 'random forest', resample_n_neuron=False, resample_iter=p.cf.resample_iter, n_tree=p.cf.n, d_tree=p.cf.d)
    
    # without resampling concatenating all neurons
    rez.all_block = train_classifier(X_all, y, 'random forest', resample_n_neuron=False, resample_iter=1, n_tree=p.cf.n, d_tree=p.cf.d)
    rez.all_dir = train_classifier(X_all, y_dir, 'random forest', resample_n_neuron=False, resample_iter=1, n_tree=p.cf.n, d_tree=p.cf.d)
    rez.all_trq = train_classifier(X_all, y_trq, 'random forest', resample_n_neuron=False, resample_iter=1, n_tree=p.cf.n, d_tree=p.cf.d)
    return rez, p

    


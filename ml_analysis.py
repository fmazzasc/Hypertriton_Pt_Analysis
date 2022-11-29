#!/usr/bin/env python3

import sys
sys.path.append('utils')

import helpers as hp
from sklearn.model_selection import train_test_split
from hipe4ml.tree_handler import TreeHandler
from hipe4ml.model_handler import ModelHandler
from hipe4ml.analysis_utils import train_test_generator
from hipe4ml import analysis_utils, plot_utils
import yaml
import xgboost as xgb
import uproot
import ROOT
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import warnings
import pickle
import os


pd.options.mode.chained_assignment = None  # default='warn'


def rename_mc_df_columns(df):
    rename_dict = {}
    rename_dict['ptMC'] = 'gPt'
    rename_dict['etaMC'] = 'gEta'
    rename_dict['ctMC'] = 'gCt'
    rename_dict['yMC'] = 'gRapidity'
    df.rename(columns = rename_dict, inplace=True)
    df['gMatter'] = df['pdg'] > 0

parser = argparse.ArgumentParser(prog='ml_analysis', allow_abbrev=True)
parser.add_argument('-split', action='store_true')
parser.add_argument('-train', action='store_true')
parser.add_argument('-application', action='store_true')
parser.add_argument('config', help='Path to the YAML configuration file')
args = parser.parse_args()

with open(os.path.expandvars(args.config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
##################################

SPLIT = args.split

# training
PLOT_DIR = 'plots'
OPTIMIZE = False
OPTIMIZED = False
TRAIN = args.train
# application

APPLICATION = args.application


# avoid pandas warning
warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.gROOT.SetBatch()

##################################################################
# read configuration file
##################################################################

DATA_PATH = params['DATA_PATH']
AOD = params['AOD']

MC_PATH = params['MC_SIGNAL_PATH']
BKG_PATH = params['LS_BACKGROUND_PATH']
PT_BINS_CENT = params['PT_BINS_CENT']
CENTRALITY_LIST = params['CENTRALITY_LIST']
TRAINING_COLUMNS_LIST = params['TRAINING_COLUMNS']
RANDOM_STATE = params['RANDOM_STATE']
HYPERPARAMS = params['HYPERPARAMS']
HYPERPARAMS_RANGES = params['HYPERPARAMS_RANGES']
RESULTS_SUBDIR = params['RESULTS_SUBDIR']

USETOF = params['USETOF']
KINT7 = params['KINT7']
TRIGGERED = params['TRIGGERED']
MAX_EFF = params['MAX_EFF']




res_dir = 'results' + RESULTS_SUBDIR
if not os.path.isdir(res_dir):
    os.mkdir(res_dir)

model_dir = 'models' + RESULTS_SUBDIR
##################################################################
# split matter/antimatter
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']


if TRAIN:
    score_eff_arrays_dict = dict()
    if AOD:
        signal_tree_handler = TreeHandler(MC_PATH, "HyperTree")
        rename_mc_df_columns(signal_tree_handler._full_data_frame)
    else:
        signal_tree_handler = TreeHandler(MC_PATH, "SignalTable")

    if BKG_PATH=='':
        background_tree_handler = TreeHandler(DATA_PATH) if DATA_PATH[-1]!='t' else TreeHandler(DATA_PATH, 'DataTable' if not AOD else "HyperTree")
        background_tree_handler.apply_preselections("m>3.91 or m<2.975")
    else:
        background_tree_handler = TreeHandler(BKG_PATH, "DataTable")

    # make plot directory
    if not os.path.isdir(PLOT_DIR):
        os.mkdir(PLOT_DIR)

    # make dataframe directory
    if not os.path.isdir('df'):
        os.mkdir('df')

    root_file_presel_eff = ROOT.TFile(res_dir + "/PreselEff.root", "recreate")

    for i_cent_bins in range(len(CENTRALITY_LIST)):
        cent_bins = CENTRALITY_LIST[i_cent_bins]

        print('##########################################################################')
        print(f'CENTRALITY BIN: {cent_bins[0]}-{cent_bins[1]}\n')
        print('****************** Reweighting MC and computing Efficiency ***************')

        for split in SPLIT_LIST:
            split_ineq_sign = '> -0.1'
            if SPLIT:
                split_ineq_sign = '> 0.5'
            if split == 'antimatter':
                split_ineq_sign = '< 0.5'

            bw_file = ROOT.TFile.Open('utils/He3_fits.root')
            bw = bw_file.Get(f'BlastWave_{cent_bins[0]}_{cent_bins[1]}')
            print(f'Applying rejection with He3 spectra ...{cent_bins[0]}_{cent_bins[1]}')
            hp.apply_pt_rejection(signal_tree_handler, bw)
            bw_file.Close()              

            root_file_presel_eff.cd()
            if(cent_bins[0]==0 and USETOF):
                df_signal_cent = signal_tree_handler.apply_preselections(f'Matter {split_ineq_sign} and rej>0 and pt>0 and abs(Rapidity)<0.5 and ct < 35 and -4<TOFnSigmaHe3<4', inplace=False)
            else:
                df_signal_cent = signal_tree_handler.apply_preselections(f'Matter {split_ineq_sign} and rej>0 and pt>0 and abs(Rapidity)<0.5 and ct < 35', inplace=False)

            if AOD:
                df_signal_cent.apply_preselections('isReconstructed==True')

            df_generated_cent = signal_tree_handler.apply_preselections(f'gMatter {split_ineq_sign} and rej>0 and abs(gRapidity)<0.5', inplace=False)
            print(len(df_signal_cent), len(df_generated_cent))
            pt_bins_cent = np.sort(pd.unique([item for sublist in PT_BINS_CENT[i_cent_bins] for item in sublist]))  # flattening list
            hist_eff_pt = hp.presel_eff_hist([df_signal_cent, df_generated_cent], 'pt', split, cent_bins, pt_bins_cent)
            hist_eff_pt.Write()
            del df_generated_cent, df_signal_cent

            print('****************** Start BDT trainings ***************')
            pt_bins_cent = PT_BINS_CENT[i_cent_bins]
            for pt_bins in pt_bins_cent:
                print(f'Matter {split_ineq_sign} and ct < 35 and pt > {pt_bins[0]} and pt < {pt_bins[1]} and rej>0')

                if(cent_bins[0]==0 and USETOF):
                    signal_tree_handler_pt = signal_tree_handler.apply_preselections(f'Matter {split_ineq_sign} and ct < 35 and pt > {pt_bins[0]} and pt < {pt_bins[1]} and rej>0 and -4<TOFnSigmaHe3<4', inplace=False)
                else:
                    signal_tree_handler_pt = signal_tree_handler.apply_preselections(f'Matter {split_ineq_sign} and ct < 35 and pt > {pt_bins[0]} and pt < {pt_bins[1]} and rej>0', inplace=False)

                
                if AOD:
                    signal_tree_handler_pt.apply_preselections('isReconstructed==True')

                background_tree_handler_pt = background_tree_handler.apply_preselections(f'ct < 35 and pt > {pt_bins[0]} and pt < {pt_bins[1]} and {cent_bins[0]}<centrality<{cent_bins[1]}', inplace=False)

                print('Number of signal candidates: ',
                        len(signal_tree_handler_pt))
                print('Number of bkg candidates: ',
                        len(background_tree_handler_pt))

                if len(background_tree_handler_pt) > 5*len(signal_tree_handler_pt):
                    background_tree_handler_pt.shuffle_data_frame(
                        size=5*len(signal_tree_handler_pt))

                bin_results = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}'
                print("BIN RESULTS: ", bin_results)

                train_test_data = train_test_generator([signal_tree_handler_pt, background_tree_handler_pt], [1, 0], test_size=0.5, random_state=RANDOM_STATE)

                ##############################################################
                # TRAINING AND TEST SET PREPARATION
                ##############################################################
                # features plot
                leg_labels = ['background', 'signal']
                model_clf = xgb.XGBClassifier(use_label_encoder=False)
                model_hdl = ModelHandler(model_clf, TRAINING_COLUMNS_LIST)
                model_hdl.set_model_params(HYPERPARAMS)
                model_hdl.set_model_params({"n_jobs": 60})
                # hyperparameters optimization and model training

                if not os.path.isdir(model_dir):
                    os.mkdir(model_dir)

                if OPTIMIZE and TRAIN:
                    model_hdl.optimize_params_bayes(
                        train_test_data, HYPERPARAMS_RANGES, 'roc_auc', nfold=5, init_points=10, n_iter=10, njobs=40)

                model_hdl.train_test_model(train_test_data)
                model_file_name = str(f'{model_dir}/{bin_results}_trained')
                if OPTIMIZE:
                    model_file_name = str(
                        f'{model_dir}/{bin_results}_optimized_trained')
                model_hdl.dump_model_handler(model_file_name)

                # get predictions for training and test sets in a given centrality
                centrality_map = np.logical_and(train_test_data[2]['centrality']>cent_bins[0], train_test_data[2]['centrality']<cent_bins[1])
                test_y_score = model_hdl.predict(train_test_data[2])

                # get scores corresponding to BDT efficiencies using test set
                eff_array = np.arange(0.10, MAX_EFF, 0.01)
                score_array = analysis_utils.score_from_efficiency_array(
                    train_test_data[3][centrality_map], test_y_score[centrality_map], efficiency_selected=eff_array, keep_lower=False)
                score_eff_arrays_dict[bin_results] = [score_array, eff_array]

                train_test_data[2]['model_output'] = test_y_score
                mc_data_cent = train_test_data[2][train_test_data[3] > 0.5]
                mc_data_cent.to_parquet(f'df{RESULTS_SUBDIR}/mc_{bin_results}', compression='gzip')
                # save roc-auc
                del train_test_data
                    ##############################################################

    pickle.dump(score_eff_arrays_dict, open(
        res_dir + "/file_score_eff_dict", "wb"))

# apply model to data
if APPLICATION:    
    if not os.path.isdir('df'):
        os.mkdir('df')
    
    if not os.path.isdir(f'df{RESULTS_SUBDIR}'):
        os.mkdir(f'df{RESULTS_SUBDIR}')

    score_eff_arrays_dict = pickle.load(
        open(res_dir + "/file_score_eff_dict", "rb"))

    for split in SPLIT_LIST:
        df_data = TreeHandler(DATA_PATH) if DATA_PATH[-1]!='t' else TreeHandler(DATA_PATH, 'DataTable' if not AOD else "HyperTree")
        df_data = df_data.get_data_frame()

        split_ineq_sign = '> -0.1'
        if SPLIT:
            split_ineq_sign = '> 0.5'
            if split == 'antimatter':
                split_ineq_sign = '< 0.5'

        for i_cent_bins, pt_bins_cent in enumerate(PT_BINS_CENT):
            for pt_bins in pt_bins_cent:
                cent_bins = CENTRALITY_LIST[i_cent_bins]

                bin_results = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}'
                if cent_bins[0]==0 and USETOF:
                    df_data_cent = df_data.query(f'Matter {split_ineq_sign} and centrality >= {cent_bins[0]} and centrality < {cent_bins[1]} and pt > {pt_bins[0]} and pt < {pt_bins[1]} and ct < 35 and -4<TOFnSigmaHe3<4')
                else:
                    df_data_cent = df_data.query(f'Matter {split_ineq_sign} and centrality >= {cent_bins[0]} and centrality < {cent_bins[1]} and pt > {pt_bins[0]} and pt < {pt_bins[1]} and ct < 35')

                if KINT7:
                    print('selecting kINT7..')
                    df_data_cent = df_data_cent.query('trigger%2!=0')
                
                if TRIGGERED:
                    print('selecting TRIGGERED events..')
                    df_data_cent = df_data_cent.query('trigger%2==0')

                model_hdl = ModelHandler()

                if OPTIMIZED:
                    model_hdl.load_model_handler(
                        f'{model_dir}/{bin_results}_optimized_trained')
                else:
                    model_hdl.load_model_handler(
                        f'{model_dir}/{bin_results}_trained')

                eff_array = np.arange(0.10, MAX_EFF, 0.01)
                print(bin_results)
                data_y_score = model_hdl.predict(df_data_cent)
                df_data_cent['model_output'] = data_y_score

                df_data_cent = df_data_cent.query(
                    f'model_output > {score_eff_arrays_dict[bin_results][0][-1]}')
                df_data_cent.to_parquet(
                    f'df{RESULTS_SUBDIR}/{bin_results}', compression='gzip')

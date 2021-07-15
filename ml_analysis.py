#!/usr/bin/env python3
import os
import pickle
import warnings
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import uproot
import xgboost as xgb
import yaml
from hipe4ml import analysis_utils, plot_utils
from hipe4ml.analysis_utils import train_test_generator
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from sklearn.model_selection import train_test_split

pd.options.mode.chained_assignment = None  # default='warn'


def presel_eff_hist(df_list, col_name, split, cent_bins, bins):
    # fill histograms (vs. ct and vs. pt)
    bins_array = np.asarray(bins, dtype=float)
    hist_eff = ROOT.TH1F(
        f'fPreselEff_vs_{col_name}_{split}_{cent_bins[0]}_{cent_bins[1]}',
        f'Preselection Efficiency, {split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins) - 1, bins_array)
    hist_gen = ROOT.TH1F('fPreselGen_vs_{col_name}', 'Gen', len(bins)-1, bins_array)

    for val in df_list[0][col_name]:
        hist_eff.Fill(val)
    for val in df_list[1][col_name]:
        hist_gen.Fill(val)

    # compute efficiency and set properties
    hist_eff.Divide(hist_eff, hist_gen, 1, 1, "B")
    if col_name == 'ct':
        hist_eff.GetXaxis().SetTitle('#it{c}t (cm)')
    elif col_name == 'pt':
        hist_eff.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    hist_eff.GetYaxis().SetTitle('Efficiency')
    hist_eff.SetMinimum(0)
    hist_eff.SetDrawOption("histo")
    hist_eff.SetLineWidth(2)

    # return histogram
    return hist_eff


parser = argparse.ArgumentParser(prog='ml_analysis', allow_abbrev=True)
parser.add_argument('-split', action='store_true')
parser.add_argument('-eff', action='store_true')
parser.add_argument('-train', action='store_true')
parser.add_argument('-computescoreff', action='store_true')
parser.add_argument('-application', action='store_true')
args = parser.parse_args()

SPLIT = args.split
MAX_EFF = 0.9
DUMP_HYPERPARAMS = False

# training
TRAINING = not args.application
PLOT_DIR = 'plots'
MAKE_PRESELECTION_EFFICIENCY = args.eff
MAKE_FEATURES_PLOTS = False
MAKE_TRAIN_TEST_PLOT = False
OPTIMIZE = False
OPTIMIZED = False
TRAIN = args.train
COMPUTE_SCORES_FROM_EFF = args.computescoreff
# application
APPLICATION = args.application

# avoid pandas warning
warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.gROOT.SetBatch()

##################################################################
# read configuration file
##################################################################
config = 'config.yaml'
with open(os.path.expandvars(config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

DATA_PATH = params['DATA_PATH']
MC_PATH = params['MC_SIGNAL_PATH']
BKG_PATH = params['LS_BACKGROUND_PATH']
CT_BINS = params['CT_BINS']
PT_BINS_CENT = params['PT_BINS_CENT']
PT_BINS = params['PT_BINS']
CENTRALITY_LIST = params['CENTRALITY_LIST']
TRAINING_COLUMNS_LIST = params['TRAINING_COLUMNS']
RANDOM_STATE = params['RANDOM_STATE']
HYPERPARAMS = params['HYPERPARAMS']
HYPERPARAMS_RANGES = params['HYPERPARAMS_RANGES']
##################################################################
test_mode = False

# split matter/antimatter
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

if TRAIN:

    score_eff_arrays_dict = dict()

    signal_tree_handler = TreeHandler(MC_PATH, "SignalTable")
    background_tree_handler =  TreeHandler(BKG_PATH, "DataTable")

    signal_tree_handler.apply_preselections("ct<35 and 2<pt<9")
    background_tree_handler.apply_preselections("ct<35 and 2<pt<9")
    if test_mode:
        signal_tree_handler.shuffle_data_frame(size=int(1e4))
        background_tree_handler.shuffle_data_frame(size=int(1e4))

    # make plot directory
    if not os.path.isdir(PLOT_DIR):
        os.mkdir(PLOT_DIR)

    # make dataframe directory
    if not os.path.isdir('df'):
        os.mkdir('df')

    if not os.path.isdir(f'results'):
        os.mkdir(f'results')

    root_file_presel_eff = ROOT.TFile("results/PreselEff.root", "recreate")

    for split in SPLIT_LIST:

        split_ineq_sign = '> -0.1'
        if SPLIT:
            split_ineq_sign = '> 0.5'
            if split == 'antimatter':
                split_ineq_sign = '< 0.5'

        for i_cent_bins in range(len(CENTRALITY_LIST)):
            cent_bins = CENTRALITY_LIST[i_cent_bins]

            if MAKE_PRESELECTION_EFFICIENCY and not MAKE_FEATURES_PLOTS and not MAKE_TRAIN_TEST_PLOT:
                ##############################################################
                # PRESELECTION EFFICIENCY
                ##############################################################
                df_generated = uproot.open(os.path.expandvars(MC_PATH))['GenTable'].arrays(library="pd")
                df_signal_cent = signal_tree_handler.get_data_frame().query(
                    f'Matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]} and pt > 2 and pt < 9 and 0<ct<35 and -0.5<Rapidity<0.5')
                df_generated_cent = df_generated.query(
                    f'matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]} and -0.5<rapidity<0.5')
                del df_generated

                # print(PT_BINS_CENT[i_cent_bins])
                pt_bins_cent = np.sort(pd.unique([item for sublist in PT_BINS_CENT[i_cent_bins] for item in sublist])) #flattening list
                ct_bins_cent = np.sort(pd.unique([item for sublist in CT_BINS for item in sublist])) #flattening list
                print(pt_bins_cent)

                # fill histograms (vs. ct and vs. pt)
                hist_eff_ct = presel_eff_hist([df_signal_cent, df_generated_cent], 'ct', split, cent_bins, ct_bins_cent)
                hist_eff_pt = presel_eff_hist([df_signal_cent, df_generated_cent], 'pt', split, cent_bins, pt_bins_cent)

                # plot histograms
                if not os.path.isdir(f'{PLOT_DIR}/presel_eff'):
                    os.mkdir(f'{PLOT_DIR}/presel_eff')
                c1 = ROOT.TCanvas()
                ROOT.gStyle.SetOptStat(0)
                c1.cd()
                c1.SetTicks(1,1)
                hist_eff_ct.Draw("histo")
                c1.Print(f'{PLOT_DIR}/presel_eff/hPreselEffVsCt_{split}_{cent_bins[0]}_{cent_bins[1]}.png')
                hist_eff_pt.Draw("histo")
                c1.Print(f'{PLOT_DIR}/presel_eff/hPreselEffVsPt_{split}_{cent_bins[0]}_{cent_bins[1]}.png')

                hist_eff_ct.Write()
                hist_eff_pt.Write()

                del df_signal_cent
                del df_generated_cent
                ##############################################################


    root_file_presel_eff.Close()

    # second condition needed because of issue with Qt libraries
    if MAKE_FEATURES_PLOTS and not MAKE_PRESELECTION_EFFICIENCY and not TRAIN:
        ######################################################
        # PLOT FEATURES DISTRIBUTIONS AND CORRELATIONS
        ######################################################

        # define tree handlers

        if len(background_tree_handler)>5*len(signal_tree_handler):
            background_tree_handler.shuffle_data_frame(size=5*len(signal_tree_handler))


        if not os.path.isdir(f'{PLOT_DIR}/features'):
            os.mkdir(f'{PLOT_DIR}/features')

        leg_labels = ['background', 'signal']
        plot_utils.plot_distr(
            [background_tree_handler, signal_tree_handler],
            TRAINING_COLUMNS_LIST, bins=40, labels=leg_labels, log=True, density=True, figsize=(12, 7),
            alpha=0.3, grid=False)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.savefig(f'{PLOT_DIR}/features/FeaturePlots')
        plot_utils.plot_corr([background_tree_handler], TRAINING_COLUMNS_LIST, ['background'])
        plt.savefig(f'{PLOT_DIR}/features/BackgroundCorrelationMatrix')
        plot_utils.plot_corr([signal_tree_handler], TRAINING_COLUMNS_LIST, ['signal'])
        plt.savefig(f'{PLOT_DIR}/features/SignalCorrelationMatrix')
        plt.close('all')

        exit()
        ###########################################################


    for i_cent_bins, pt_bins_cent in enumerate(PT_BINS_CENT):
        cent_bins = CENTRALITY_LIST[i_cent_bins]
        for pt_bins in pt_bins_cent:
            

            signal_tree_handler_pt = signal_tree_handler.apply_preselections(f'ct < 35 and pt > {pt_bins[0]} and pt < {pt_bins[1]}', inplace=False)
            background_tree_handler_pt = background_tree_handler.apply_preselections(f'ct < 35 and pt > {pt_bins[0]} and pt < {pt_bins[1]}', inplace=False)

            if len(background_tree_handler_pt)>5*len(signal_tree_handler_pt):
                background_tree_handler_pt.shuffle_data_frame(size=5*len(signal_tree_handler_pt))

            train_test_data = train_test_generator([signal_tree_handler_pt, background_tree_handler_pt], [1, 0], test_size=0.5, random_state=RANDOM_STATE)

            for split in SPLIT_LIST:
                split_ineq_sign = '> -0.1'
                if SPLIT:
                    split_ineq_sign = '> 0.5'
                    if split == 'antimatter':
                        split_ineq_sign = '< 0.5'

                bin_results = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}'
                bin_model = f'{pt_bins[0]}_{pt_bins[1]}'

                ##############################################################
                # TRAINING AND TEST SET PREPARATION
                ##############################################################
                # features plot
                leg_labels = ['background', 'signal']
                model_clf = xgb.XGBClassifier(use_label_encoder=False)
                model_hdl = ModelHandler(model_clf, TRAINING_COLUMNS_LIST)
                model_hdl.set_model_params(HYPERPARAMS)
                model_hdl.set_model_params({"n_jobs":30})
                # hyperparameters optimization and model training
                if not os.path.isdir('models'):
                    os.mkdir('models')
                if OPTIMIZE and TRAIN:
                    model_hdl.optimize_params_bayes(train_test_data, HYPERPARAMS_RANGES,
                                                    'roc_auc', nfold=5, init_points=10, n_iter=10, njobs=40)
                isModelTrained = os.path.isfile(f'models/{bin_model}_trained')
                print(f'isModelTrained {bin_model}: {isModelTrained}')
                if TRAIN and not isModelTrained:
                    print(
                    f'Number of candidates ({split}) for training in {pt_bins[0]} <= pt < {pt_bins[1]} GeV/c: {len(train_test_data[0])}')
                    print(
                    f'signal candidates: {np.count_nonzero(train_test_data[1] == 1)}; background candidates: {np.count_nonzero(train_test_data[1] == 0)}; n_cand_bkg / n_cand_signal = {np.count_nonzero(train_test_data[1] == 0) / np.count_nonzero(train_test_data[1] == 1)}')

                    model_hdl.train_test_model(train_test_data)
                    model_file_name = str(f'models/{bin_model}_trained')
                    if OPTIMIZE:
                        model_file_name = str(f'models/{bin_model}_optimized_trained')
                    model_hdl.dump_model_handler(model_file_name)
                elif COMPUTE_SCORES_FROM_EFF:
                    if OPTIMIZED:
                        model_hdl.load_model_handler(f'models/{bin_model}_optimized_trained')
                    else:
                        model_hdl.load_model_handler(f'models/{bin_model}_trained')
                else:
                    continue


                    # get only centrality selected
                
                print("BIN RESULTS: ", bin_results)
                signal_tree_handler_cent = signal_tree_handler_pt.apply_preselections(f'Matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]}', inplace=False)
                background_tree_handler_cent = background_tree_handler_pt.apply_preselections(f'Matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]}', inplace=False)

                train_test_data_cent = train_test_generator([signal_tree_handler_cent, background_tree_handler_cent], [1, 0], test_size=0.5, random_state=RANDOM_STATE)


                # get predictions for training and test sets
                test_y_score = model_hdl.predict(train_test_data_cent[2])
                train_y_score = model_hdl.predict(train_test_data_cent[0])

                # second condition needed because of issue with Qt libraries
                if MAKE_TRAIN_TEST_PLOT:
                    if not os.path.isdir(f'{PLOT_DIR}/train_test_out'):
                        os.mkdir(f'{PLOT_DIR}/train_test_out')
                    plot_utils.plot_output_train_test(model_hdl, train_test_data_cent,
                                                        logscale=True, density=True, labels=leg_labels)
                    plt.savefig(f'{PLOT_DIR}/train_test_out/{bin_results}_out.png')

                    # plot_utils.plot_feature_imp(train_test_data_cent[0], train_test_data_cent[1], model_hdl)
                    # plt.savefig(f'{PLOT_DIR}/train_test_out/feature_imp_training_{bin_results}.png')
                    plot_utils.plot_roc_train_test(
                        train_test_data_cent[3],
                        test_y_score, train_test_data_cent[1],
                        train_y_score, labels=leg_labels)
                    plt.savefig(f'{PLOT_DIR}/train_test_out/roc_train_test_{bin_results}.png')
                    plt.close('all')

                if COMPUTE_SCORES_FROM_EFF:
                    # get scores corresponding to BDT efficiencies using test set
                    eff_array = np.arange(0.10, MAX_EFF, 0.01)
                    score_array = analysis_utils.score_from_efficiency_array(
                        train_test_data_cent[3], test_y_score, efficiency_selected=eff_array, keep_lower=False)
                    score_eff_arrays_dict[bin_results] = score_array

                    # write test set data frame
                    train_test_data_cent[2]['model_output'] = test_y_score
                    mc_data_cent = train_test_data_cent[2][train_test_data_cent[3]>0.5]
                    mc_data_cent.to_parquet(f'df/mc_{bin_results}', compression='gzip')

                    # get the model hyperparameters
                    if DUMP_HYPERPARAMS:
                        if not os.path.isdir('hyperparams'):
                            os.mkdir('hyperparams')
                        model_params_dict = model_hdl.get_model_params()
                        with open(f'hyperparams/model_params_{bin_model}.yml', 'w') as outfile:
                            yaml.dump(model_params_dict, outfile, default_flow_style=False)

                    # save roc-auc
                del train_test_data_cent
                ##############################################################

    pickle.dump(score_eff_arrays_dict, open("results/file_score_eff_dict", "wb"))

# apply model to data
if APPLICATION:
    if not os.path.isdir('df'):
        os.mkdir('df')
    score_eff_arrays_dict = pickle.load(open("results/file_score_eff_dict", "rb"))

    for split in SPLIT_LIST:
        df_data = uproot.open(os.path.expandvars(DATA_PATH))['DataTable'].arrays(library="pd")

        split_ineq_sign = '> -0.1'
        if SPLIT:
            split_ineq_sign = '> 0.5'
            if split == 'antimatter':
                split_ineq_sign = '< 0.5'

        for i_cent_bins, pt_bins_cent in enumerate(PT_BINS_CENT):
            for pt_bins in pt_bins_cent:
                cent_bins = CENTRALITY_LIST[i_cent_bins]

                bin_results = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}'
                bin_model = f'{pt_bins[0]}_{pt_bins[1]}'

                df_data_cent = df_data.query(
                    f'Matter {split_ineq_sign} and centrality > {cent_bins[0]} and centrality < {cent_bins[1]} and pt > {pt_bins[0]} and pt < {pt_bins[1]} and ct < 35')

                model_hdl = ModelHandler()

                if OPTIMIZED:
                    model_hdl.load_model_handler(f'models/{bin_model}_optimized_trained')
                else:
                    model_hdl.load_model_handler(f'models/{bin_model}_trained')

                eff_array = np.arange(0.10, MAX_EFF, 0.01)
                print(bin_results)
                data_y_score = model_hdl.predict(df_data_cent)
                df_data_cent['model_output'] = data_y_score

                df_data_cent = df_data_cent.query(f'model_output > {score_eff_arrays_dict[bin_results][len(eff_array)-1]}')
                df_data_cent.to_parquet(f'df/{bin_results}', compression='gzip')

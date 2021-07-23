#!/usr/bin/env python3
import sys
sys.path.append('utils')
import os
import pickle
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm
import ROOT
import uproot
import yaml


from helpers import significance_error, expected_signal_pt

SPLIT = True
MAX_EFF = 0.9

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

ANALYSIS_RESULTS_PATH = params['ANALYSIS_RESULTS_PATH']
PT_BINS_CENT = params['PT_BINS_CENT']
CENTRALITY_LIST = params['CENTRALITY_LIST']
RANDOM_STATE = params['RANDOM_STATE']
##################################################################

# split matter/antimatter
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']



RESULTS_SUBDIR = params['RESULTS_SUBDIR']
res_dir = 'results' + RESULTS_SUBDIR
if not os.path.isdir(res_dir):
    os.mkdir(res_dir)



score_eff_arrays_dict = pickle.load(open(res_dir + "/file_score_eff_dict", "rb"))
eff_array = np.arange(0.10, MAX_EFF, 0.01)
presel_eff_file = uproot.open(res_dir + '/PreselEff.root')
analysis_results_file = uproot.open(os.path.expandvars(ANALYSIS_RESULTS_PATH))

# get centrality selected histogram
cent_counts, cent_edges = analysis_results_file["AliAnalysisTaskHyperTriton2He3piML_custom_summary;1"][11].to_numpy()
cent_bin_centers = (cent_edges[:-1]+cent_edges[1:])/2

# cut dictionary
eff_cut_dict = dict()

for split in SPLIT_LIST:
    for i_cent_bins, pt_bins_cent in enumerate(PT_BINS_CENT):
        cent_bins = CENTRALITY_LIST[i_cent_bins]
        # get number of events
        cent_range_map = np.logical_and(cent_bin_centers > cent_bins[0], cent_bin_centers < cent_bins[1])
        counts_cent_range = cent_counts[cent_range_map]
        evts = np.sum(counts_cent_range)
        print(f'Number of events: {evts}')
        # get preselection efficiency histogram
        presel_eff_counts, presel_eff_edges = presel_eff_file[
            f'fPreselEff_vs_pt_{split}_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
        print(presel_eff_counts, presel_eff_edges)
        presel_eff_bin_centers = (presel_eff_edges[1:]+presel_eff_edges[:-1])/2

        for pt_bins in pt_bins_cent:

            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}'
            df_data = pd.read_parquet(f'df/{bin}')

            # plot directory
            if not os.path.isdir('plots/significance_scan'):
                os.mkdir('plots/significance_scan')

            # lists filled inside loop
            significance_list = []
            significance_err_list = []

            for eff_score in zip(eff_array, score_eff_arrays_dict[bin]):
                # if (ct_bins[0] > 0) and (eff_score[0] < 0.49):
                #     continue
                formatted_eff = "{:.2f}".format(eff_score[0])
                print(f'processing {bin}: eff = {eff_score[0]:.2f}, score = {eff_score[1]:.2f}...')

                # select data in the sidebands
                df_data_sel = df_data.query(f'model_output > {eff_score[1]}')

                # make histogram
                hyp_mass = 2.991
                sigma = 0.0015
                m_min = 2.98 #hyp_mass - 3*sigma
                m_max = 3.005 #hyp_mass + 3*sigma
                nBins = 26

                counts, bin_edges = np.histogram(np.array(df_data_sel['m']), nBins, range=(2.960, 3.025))
                bin_centers = (bin_edges[1:]+bin_edges[:-1])/2.
                side_map = np.logical_or(bin_centers < m_min, bin_centers > m_max)
                mass_map = np.logical_not(side_map)

                side_bins = bin_centers[side_map]
                side_counts = counts[side_map]
                side_errors = np.sqrt(side_counts)

                # polynomial fit to background
                pol = np.polynomial.Polynomial.fit(side_bins, side_counts, deg=1)
                xx_side, yy_side = pol.linspace()  # plot polynomial

                # compute background
                pol_integral = pol.integ()
                bin_size = (3.025-2.960)/nBins
                bkg = (pol_integral(m_max) - pol_integral(m_min))/bin_size

                # compute eff = presel_eff * BDT_eff
                presel_eff_map = np.logical_and(
                    presel_eff_bin_centers > pt_bins[0],
                    presel_eff_bin_centers < pt_bins[1])
                presel_eff = presel_eff_counts[presel_eff_map]
                eff = presel_eff * eff_score[0]
                print(eff_score)
                print(presel_eff)

                # compute expected signal
                sig = expected_signal_pt(cent_bins, pt_bins, eff, cent_counts)
                if not SPLIT:
                    sig *= 2
                mass_bins = bin_centers[mass_map]
                mass_counts = norm.pdf(mass_bins, hyp_mass, sigma)
                mass_counts = mass_counts*sig*bin_size
                pol_offset = pol(mass_bins)
                mass_counts = mass_counts+pol_offset

                xx_mass = np.linspace(m_min, m_max, 100)
                yy_mass = norm.pdf(xx_mass, hyp_mass, sigma)*sig*bin_size
                yy_offset = pol(xx_mass)  # plot polynomial
                yy_mass = yy_mass+yy_offset

                # compute significance
                significance_list.append(sig/np.sqrt(sig+bkg))
                significance_err_list.append(significance_error(sig, bkg))

                # plot histograms
                if not os.path.isdir(f'plots/significance_scan/{bin}'):
                    os.mkdir(f'plots/significance_scan/{bin}')
                fig = plt.figure()
                plt.errorbar(side_bins, side_counts, side_errors, fmt='o', label='Data', color='blue', ecolor='black')
                plt.errorbar(mass_bins, mass_counts, np.sqrt(mass_counts), fmt='o', label='Pseudodata', color='red', ecolor='black')
                plt.plot(xx_side, yy_side, label='Background fit', color='green')
                plt.plot(xx_mass, yy_mass, label='Gaussian model', color='orange')
                plt.xlabel(r'$M \ (^{3}He + \pi^{-}) \ (\mathrm{GeV}/\it{c}^{2})$')
                plt.ylabel(r'$\mathrm{Entries}/(2.5 \mathrm{MeV}/\it{c}^{2})$')
                handles, labels = fig.gca().get_legend_handles_labels()
                order = [2, 3, 0, 1]
                plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='upper right')
                plt.savefig(f'plots/significance_scan/{bin}/{formatted_eff}_{bin}.png')
                plt.close('all')

            eff_array_reduced = eff_array
            # if (ct_bins[0] == 0):
            #     eff_array_reduced = eff_array
            significance_array = np.asarray(significance_list)
            significance_err_array = np.asarray(significance_err_list)

            print(len(significance_array))
            significance_array = np.nan_to_num(significance_array*(eff_array_reduced))
            significance_err_array = np.nan_to_num(significance_err_array*(eff_array_reduced))



            low_limit = significance_array - significance_err_array
            up_limit = significance_array + significance_err_array

            # cut
            print(significance_array)
            cut_eff_index = significance_array == np.nanmax(significance_array)
            print(cut_eff_index)
            cut_eff = eff_array_reduced[cut_eff_index]

            fig = plt.figure()
            plt.plot(eff_array_reduced, significance_array, 'b', label='Expected Significance')
            plt.fill_between(eff_array_reduced, low_limit, up_limit,
                             facecolor='deepskyblue', label=r'$ \pm 1\sigma$', alpha=0.3)
            print("cut eff: ", cut_eff)
            plt.vlines(cut_eff[0], 0, 20, colors='red', linestyles='dashed',
                       label=f'cut = {"{:.2f}".format(cut_eff[0])}')

            handles, labels = fig.gca().get_legend_handles_labels()
            order = [0, 2, 1]
            plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower left')

            plt.xlabel("BDT Efficiency")
            plt.ylabel("Significance x BDT Efficiency")
            plt.xlim(0.5, MAX_EFF-0.01)
            # if ct_bins[0] == 0:
            #     plt.xlim(0.1, MAX_EFF-0.01)
            plt.ylim(0.3, up_limit.max()+0.3)

            plt.savefig(f'plots/significance_scan/{bin}.png')
            plt.close('all')

            eff_cut_dict[bin] = cut_eff[0]
            
            # get the systematic variation range
            i_eff = 1
            while i_eff < 11:
                tmp_eff_var_right = cut_eff[0] + i_eff/100
                tmp_eff_var_left = cut_eff[0] - i_eff/100
                
                # get indices in efficiency array
                i_eff_var_right = np.where(eff_array_reduced == tmp_eff_var_right)[0]
                i_eff_var_left = np.where(eff_array_reduced == tmp_eff_var_left)[0]
                
                # find corresponding significance
                tmp_signif_left = significance_array[i_eff_var_right] / tmp_eff_var_right
                tmp_signif_right = significance_array[i_eff_var_left] / tmp_eff_var_left
                
                # if one of the two is smaller than 3, exit
                if tmp_signif_left<3. or tmp_signif_right<3. :
                    break
                
                else:
                    i_eff += 1

            bin_range = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}_range'
            eff_cut_dict[bin_range] = i_eff
            print(f'BDT efficiency cut variation range: +/-{i_eff-1}%')

pickle.dump(eff_cut_dict, open(res_dir + "/file_eff_cut_dict", "wb"))

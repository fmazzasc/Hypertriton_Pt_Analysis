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
SPLIT=True
SPLIT_LIST = ['']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

if not os.path.isdir("plots/cpt_and_lifetime"):
    os.mkdir("plots/cpt_and_lifetime")

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

eff_cut_dict = pickle.load(open("results/file_eff_cut_dict", "rb"))
presel_eff_file = uproot.open('results/PreselEff.root')
signal_extraction_file = ROOT.TFile.Open('results/SignalExtraction.root')
signal_extraction_keys = uproot.open('results/SignalExtraction.root').keys()

abs_correction_file = ROOT.TFile.Open('results/He3_abs.root')
pt_spectra_file = ROOT.TFile.Open('results/pt_spectra.root', 'recreate')



for i_cent_bins, pt_bins_cent in enumerate(PT_BINS_CENT):
    cent_bins = CENTRALITY_LIST[i_cent_bins]
    for i_split,split in enumerate(SPLIT_LIST):
        print(f'{i_split} -> {split}')
        # get preselection efficiency and abs correction histograms
        presel_eff_counts, presel_eff_edges = presel_eff_file[f'fPreselEff_vs_pt_{split}_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
        presel_eff_bin_centers = (presel_eff_edges[1:]+presel_eff_edges[:-1])/2
        flat_pt_bins = [item for sublist in pt_bins_cent for item in sublist]
        bins = np.array(flat_pt_bins, dtype=float)
        # print(bins)
        h_corrected_yields[i_split] = ROOT.TH1D(f'fYields_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)
        for pt_bins in pt_bins_cent:

    
            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}'
            formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin])

            # look for plot with eff = eff_cut (or the nearest one)
            bkg_shape = 'pol1'
            eff_cut_increment = 0
            eff_cut_sign = -1
            while signal_extraction_keys.count(f'{bin}_{bkg_shape}/fInvMass_{formatted_eff_cut};1') == 0:
                if eff_cut_sign == -1:
                    eff_cut_increment += 0.01
                eff_cut_sign *= -1
                formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin]+eff_cut_increment*eff_cut_sign)

            # get signal
            h_raw_yield = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fRawYields;1')
            eff_index = h_raw_yield.FindBin(float(formatted_eff_cut))
            raw_yield = h_raw_yield.GetBinContent(eff_index)
            raw_yield_error = h_raw_yield.GetBinError(eff_index)
            print(f'eff_cut = {formatted_eff_cut}, raw_yield = {raw_yield}+{raw_yield_error}')


            presel_eff_map = np.logical_and(
                presel_eff_bin_centers > pt_bins[0],
                presel_eff_bin_centers < pt_bins[1])
            presel_eff = presel_eff_counts[presel_eff_map]
            bdt_eff = float(formatted_eff_cut)
            print(f'bin: {presel_eff_map}, presel_eff: {presel_eff}')
            eff = presel_eff * eff_cut_dict[bin]


            pt_bin_index = h_corrected_yields[i_split].FindBin(pt_bins[0]+0.5)
            h_corrected_yields[i_split].SetBinContent(pt_bin_index, raw_yield/eff[0])
            h_corrected_yields[i_split].SetBinError(pt_bin_index, raw_yield_error/eff[0])

        # set labels
        h_corrected_yields[i_split].GetXaxis().SetTitle("#it{p}_T (GeV/#it{c^2})")
        h_corrected_yields[i_split].GetYaxis().SetTitle("d#it{N}/d(#it{p}_T) (GeV/#it{c^2})")
        #h_corrected_yields[i_split].Scale(1, "width")
        for i_bin in range(len(bins))[2:]:
            bin_width = h_corrected_yields[i_split].GetBinWidth(i_bin)
            print(f"bin: {h_corrected_yields[i_split].GetBinLowEdge(i_bin)}; bin width: {bin_width}")
            bin_content = h_corrected_yields[i_split].GetBinContent(i_bin)
            bin_error = h_corrected_yields[i_split].GetBinError(i_bin)
            h_corrected_yields[i_split].SetBinContent(i_bin, bin_content/bin_width)
            h_corrected_yields[i_split].SetBinError(i_bin, bin_error/bin_width)
        h_corrected_yields[i_split].GetYaxis().SetRangeUser(1., 450.)
        h_corrected_yields[i_split].SetMarkerStyle(20)
        h_corrected_yields[i_split].SetMarkerSize(0.8)
        h_corrected_yields[i_split].Write()



pt_spectra_file.Close()



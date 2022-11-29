import sys
sys.path.append('utils')
import helpers as hp

import os
import pickle
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import ROOT
import uproot
import yaml


ROOT.gStyle.SetOptStat(0)
# ROOT.gStyle.SetOptFit(0)
ROOT.gROOT.SetBatch()

ROOT.gROOT.ProcessLine(".L utils/alipwgfunc/AliPWGFunc.cxx++")
from ROOT import AliPWGFunc

PT_BINS_CENT = [[[2, 4],[4, 9]], [[2, 3],[3,4],[4,5],[5,9]], [[2,4],[4,9]], [[2,9]]]
CENTRALITY_LIST = [[0,10],[10,30], [30, 50], [50,90]]
MAX_EFF = 0.95
##################################################################

res_dir_2015 = 'results/bins_offline_2015'
res_dir_2018 = 'results/bins_offline_2018_KINT7'
target_dir = 'results/bins_offline_merged'

if not os.path.isdir(target_dir):
    os.mkdir(target_dir)


# split matter/antimatter
SPLIT_LIST = ['antimatter', 'matter', 'all']

# bw file
#####################################################################
bw_file = ROOT.TFile('utils/BlastWaveFits.root')

## centrality
#####################################################################
analysis_results_2015 = uproot.open("/data/fmazzasc/PbPb_2body/2015/AnalysisResults_15o.root")
analysis_results_2018 = uproot.open("/data/fmazzasc/PbPb_2body/2018/AnalysisResults_18qr.root")
cent_counts, cent_edges = analysis_results_2015["AliAnalysisTaskHyperTriton2He3piML_custom_summary;1"][11].to_numpy()
cent_bin_centers = (cent_edges[:-1]+cent_edges[1:])/2


eff_cut_dict = pickle.load(open(res_dir_2015 + "/file_eff_cut_dict", "rb"))

## presel eff
#####################################################################
presel_eff_file_2015 = uproot.open(res_dir_2015 + '/PreselEff.root')
presel_eff_file_2018 = uproot.open(res_dir_2018 + '/PreselEff.root')


signal_extraction_file = ROOT.TFile.Open(target_dir + '/SignalExtraction.root')
signal_extraction_up = uproot.open(target_dir + '/SignalExtraction.root')
absorption_correction_file = uproot.open(f"utils/He3_abs_1.5.root")


pt_spectra_file = ROOT.TFile.Open(target_dir + '/systematics.root', 'recreate')
dire = pt_spectra_file.mkdir("trials")

N_TRIALS = 300

for i_cent_bins, pt_bins_cent in enumerate(PT_BINS_CENT):
    dire.cd()

    flat_pt_bins = [item for sublist in pt_bins_cent for item in sublist]
    bins = np.unique(np.array(flat_pt_bins, dtype=float))

    if cent_bins==[0,10]:
        bw = bw_file.Get('BlastWave_5_10')
    else:
        bw = bw_file.Get(f'BlastWave_{cent_bins[0]}_{cent_bins[1]}')

    if cent_bins!=[0,10]:
        continue

    cent_bins = CENTRALITY_LIST[i_cent_bins]
    cent_range_map = np.logical_and(cent_bin_centers > cent_bins[0], cent_bin_centers < cent_bins[1])
    counts_cent_range = cent_counts[cent_range_map]
    evts = np.sum(counts_cent_range)
    evts += hp.get_number_of_MB_ev(cent_bins, analysis_results_2018)  #2018 Events KINT7

    print("********************************************************")
    print(f'Number of events [{cent_bins[0]}, {cent_bins[1]}] : {evts}')


    h_corrected_yields = [ROOT.TH1D(), ROOT.TH1D(), ROOT.TH1D()]
    for i_split,split in enumerate(SPLIT_LIST):
        print("---------------------------")
        print(f'{i_split} -> {split}, cent: [{cent_bins[0]}, {cent_bins[1]}]')
        # get preselection efficiency and abs correction histograms
        if split=='all':
            presel_eff_counts_2015_matter, presel_eff_edges = presel_eff_file_2015[f'fPreselEff_vs_pt_matter_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
            presel_eff_counts_2018_matter, _ = presel_eff_file_2018[f'fPreselEff_vs_pt_matter_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
            presel_eff_counts_matter = 0.5*(presel_eff_counts_2015_matter + presel_eff_counts_2018_matter)
            presel_eff_counts_2015_antimatter, presel_eff_edges = presel_eff_file_2015[f'fPreselEff_vs_pt_antimatter_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
            presel_eff_counts_2018_antimatter, _ = presel_eff_file_2018[f'fPreselEff_vs_pt_antimatter_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
            presel_eff_counts_antimatter = 0.5*(presel_eff_counts_2015_antimatter + presel_eff_counts_2018_antimatter)
            presel_eff_counts = 0.5*(presel_eff_counts_matter + presel_eff_counts_antimatter)
        else:
            presel_eff_counts_2015, presel_eff_edges = presel_eff_file_2015[f'fPreselEff_vs_pt_{split}_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
            presel_eff_counts_2018, _ = presel_eff_file_2018[f'fPreselEff_vs_pt_{split}_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()
            presel_eff_counts = 0.5*(presel_eff_counts_2015 + presel_eff_counts_2018)    
        # get absorption correction
        func = "BlastWave" if cent_bins==[0,10] else "BGBW"
        absorption_counts, absorption_edges = absorption_correction_file[f'{cent_bins[0]}_{cent_bins[1]}'][f'fEffPt_antimatter_cent_{cent_bins[0]}_{cent_bins[1]}_func_{func};1'].to_numpy()
        presel_eff_bin_centers = (presel_eff_edges[1:]+presel_eff_edges[:-1])/2
        absorption_bin_centers = (absorption_edges[1:]+absorption_edges[:-1])/2


        if cent_bins[0]==0:
            hist_bins = [9e-6, 3.e-5]
        elif cent_bins[0]==10:
            hist_bins = [6e-6, 1.2e-5]
        else:
            hist_bins = [1e-6, 6e-6]
        trial = ROOT.TH1D(f'fParameterDistribution_{cent_bins[0]}_{cent_bins[1]}_{split}', f'{cent_bins[0]}-{cent_bins[1]}%_{split}', 60, hist_bins[0], hist_bins[1])
        trial.SetDirectory(0)
        i_trial=0
        h_corrected_yields[i_split] = ROOT.TH1D(f'fYields_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)

        while i_trial < N_TRIALS:
            print("trial n°:", i_trial)

            h_corrected_yields[i_split].Reset()
            for pt_bins in pt_bins_cent:
            
                bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}'
                formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin]
                )  
                # formatted_eff_cut = "0.80"

                # look for plot with eff = eff_cut (or the nearest one)
                bkg_shape = 'pol2'
                eff_cut_increment = 0
                eff_cut_sign = -1
                signal_extraction_keys = signal_extraction_up[f"{bin}_{bkg_shape}"].keys()
                # print('eff cut: ', formatted_eff_cut)
                if len(signal_extraction_keys)==0:
                    continue

                while signal_extraction_keys.count(f"fInvMass_{formatted_eff_cut};1")==0 and eff_cut_increment<0.20:
                    if eff_cut_sign == -1:
                        eff_cut_increment += 0.01
                    eff_cut_sign *= -1
                    formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin]+eff_cut_increment*eff_cut_sign)

                # print('eff cut: ', formatted_eff_cut)
                # random sample of cut (upper and lower limits from significance scan)
                bin_range = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}_range'
                eff_cut_range = eff_cut_dict[bin_range]/100 - 0.01
                #print(f"BDT efficiency cut variation range: {eff_cut_range}")
                lower_limit = float(formatted_eff_cut) - eff_cut_range
                upper_limit = float(formatted_eff_cut) + eff_cut_range
                if float(formatted_eff_cut) > MAX_EFF:
                    upper_limit = MAX_EFF
                random_cut = lower_limit + ROOT.gRandom.Rndm()*(upper_limit-lower_limit)
                formatted_eff_cut = "{:.2f}".format(random_cut)
                eff_cut_increment = 0
                eff_cut_sign = -1
                while signal_extraction_keys.count(f"fInvMass_{formatted_eff_cut};1")==0 and eff_cut_increment<0.20:
                    if eff_cut_sign == -1:
                        eff_cut_increment += 0.01
                    eff_cut_sign *= -1
                    formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin]+eff_cut_increment*eff_cut_sign)



                # get signal
                h_raw_yield = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fRawYields;1')
                eff_index = h_raw_yield.FindBin(float(formatted_eff_cut))
                raw_yield = h_raw_yield.GetBinContent(eff_index)
                raw_yield_error = h_raw_yield.GetBinError(eff_index)

                presel_eff_map = np.logical_and(presel_eff_bin_centers > pt_bins[0], presel_eff_bin_centers < pt_bins[1])
                absorption_map = np.logical_and(absorption_bin_centers > pt_bins[0], absorption_bin_centers < pt_bins[1])

                presel_eff = presel_eff_counts[presel_eff_map]
                absorption_corr = absorption_counts[absorption_map]
                absorption_corr = absorption_corr if split=='antimatter' else (0.98/0.95)*absorption_corr

                bdt_eff = float(formatted_eff_cut)
                eff = presel_eff * eff_cut_dict[bin]

                pt_bin_index = h_corrected_yields[i_split].FindBin(pt_bins[0]+0.2)
                h_corrected_yields[i_split].SetBinContent(pt_bin_index, raw_yield/eff[0]/absorption_corr)
                h_corrected_yields[i_split].SetBinError(pt_bin_index, raw_yield_error/eff[0]/absorption_corr)

        

                # print(f'bin: [{pt_bins[0]}, {pt_bins[1]}], BDT eff = {formatted_eff_cut}, presel_eff = {presel_eff}, raw_yield = {raw_yield}+{raw_yield_error}')


            # set labels
            h_corrected_yields[i_split].GetXaxis().SetTitle("#it{p}_T (GeV/#it{c^2})")
            h_corrected_yields[i_split].GetYaxis().SetTitle("d#it{N}/d(#it{p}_T) (GeV/#it{c^2})")

            for i_bin in range(len(bins))[1:]:
                bin_width = h_corrected_yields[i_split].GetBinWidth(i_bin)
                bin_content = h_corrected_yields[i_split].GetBinContent(i_bin)
                bin_error = h_corrected_yields[i_split].GetBinError(i_bin)
                h_corrected_yields[i_split].SetBinContent(i_bin, bin_content/bin_width/evts)
                h_corrected_yields[i_split].SetBinError(i_bin, bin_error/bin_width/evts)

                if split=='all':
                    h_corrected_yields[i_split].SetBinContent(i_bin, h_corrected_yields[i_split].GetBinContent(i_bin)/2)
                    h_corrected_yields[i_split].SetBinError(i_bin, h_corrected_yields[i_split].GetBinError(i_bin)/2)


            h_corrected_yields[i_split].SetMarkerStyle(20)
            h_corrected_yields[i_split].SetMarkerSize(0.8)

            pwg = AliPWGFunc()
            _, integral,integral_error,_ = hp.bw_fit(h_corrected_yields[i_split], bw, pwg, [2,9])

            trial.Fill(integral)

            if integral > 1.7e-5 or integral < 1.1e-5:
                h_corrected_yields[i_split].Write(f"trial_n_{i_trial}")


            i_trial += 1

        trial.GetXaxis().SetTitle("d#it{N}/dy")
        trial.GetYaxis().SetTitle("Counts")

        dire.Close()
        pt_spectra_file.cd()
        trial.Write()


pt_spectra_file.Close()



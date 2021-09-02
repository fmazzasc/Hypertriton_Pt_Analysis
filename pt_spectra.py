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

ROOT.gROOT.ProcessLine(".L utils/alipwgfunc/AliPWGFunc.cxx++")
from ROOT import AliPWGFunc


parser = argparse.ArgumentParser(prog='pt_spectra', allow_abbrev=True)
parser.add_argument('config', help='Path to the YAML configuration file')
args = parser.parse_args()

with open(os.path.expandvars(args.config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

HYP_HE_CROSS_SECT_SCALING = 96/98
ANALYSIS_RESULTS_PATH = params['ANALYSIS_RESULTS_PATH']
PT_BINS_CENT = params['PT_BINS_CENT']
CENTRALITY_LIST = params['CENTRALITY_LIST']
RESULTS_SUBDIR = params['RESULTS_SUBDIR']
##################################################################


res_dir = 'results' + RESULTS_SUBDIR
if not os.path.isdir(res_dir):
    os.mkdir(res_dir)

# split matter/antimatter
SPLIT_LIST = ['antimatter', 'matter']

ROOT.gStyle.SetOptStat(0)
# ROOT.gStyle.SetOptFit(0)
ROOT.gROOT.SetBatch()
#####################################################################
bw_file = ROOT.TFile('utils/BlastWaveFits.root')

#####################################################################
analysis_results_file = uproot.open(os.path.expandvars(ANALYSIS_RESULTS_PATH))
cent_counts, cent_edges = analysis_results_file["AliAnalysisTaskHyperTriton2He3piML_custom_summary;1"][11].to_numpy()  ##does not wprk with uproot3 




cent_bin_centers = (cent_edges[:-1]+cent_edges[1:])/2

print("Total number of events: ", np.sum(cent_counts))

eff_cut_dict = pickle.load(open(res_dir + "/file_eff_cut_dict", "rb"))
presel_eff_file = uproot.open(res_dir + '/PreselEff.root')
signal_extraction_file = ROOT.TFile.Open(res_dir + '/SignalExtraction.root')
signal_extraction_up = uproot.open(res_dir + '/SignalExtraction.root')
# absorption_correction_file = uproot.open(f"results{RESULTS_SUBDIR}/He3_abs.root")

pt_spectra_file = ROOT.TFile.Open(res_dir + '/pt_spectra.root', 'recreate')


for i_cent_bins, pt_bins_cent in enumerate(PT_BINS_CENT):


    flat_pt_bins = [item for sublist in pt_bins_cent for item in sublist]
    bins = np.unique(np.array(flat_pt_bins, dtype=float))
    cent_bins = CENTRALITY_LIST[i_cent_bins]

    cent_range_map = np.logical_and(cent_bin_centers > cent_bins[0], cent_bin_centers < cent_bins[1])
    counts_cent_range = cent_counts[cent_range_map]
    evts = np.sum(counts_cent_range)
    print("********************************************************")
    print(f'Number of events [{cent_bins[0]}, {cent_bins[1]}] : {evts}')

    h_corrected_yields = [ROOT.TH1D(), ROOT.TH1D()]
    h_corrected_ratio = ROOT.TH1D(f'fRatio_{cent_bins[0]}_{cent_bins[1]}', f'{cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)
    for i_split,split in enumerate(SPLIT_LIST):
        bw = bw_file.Get('BlastWave/BlastWave0')
        print("---------------------------")
        print(f'{i_split} -> {split}')
        # get preselection efficiency
        presel_eff_counts, presel_eff_edges = presel_eff_file[f'fPreselEff_vs_pt_{split}_{cent_bins[0]}_{cent_bins[1]};1'].to_numpy()        
        # get absorption correction
        func = "BlastWave" if cent_bins[0]<1 else "BGBW"
        # absorption_counts = absorption_correction_file[f'{cent_bins[0]}_{cent_bins[1]}'][f'fEffPt_{split}_cent_{cent_bins[0]}_{cent_bins[1]}_func_{func};1'].values
        # absorption_edges = absorption_correction_file[f'{cent_bins[0]}_{cent_bins[1]}'][f'fEffPt_{split}_cent_{cent_bins[0]}_{cent_bins[1]}_func_{func};1'].edges

        presel_eff_bin_centers = (presel_eff_edges[1:]+presel_eff_edges[:-1])/2
        # absorption_bin_centers = (absorption_edges[1:]+absorption_edges[:-1])/2


        h_corrected_yields[i_split] = ROOT.TH1D(f'fYields_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)
        for pt_bins in pt_bins_cent:
    
            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}'
            formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin])

            bkg_shape = 'pol1'
            eff_cut_increment = 0
            eff_cut_sign = -1
            signal_extraction_keys = signal_extraction_up[f"{bin}_{bkg_shape}"].keys()
            if len(signal_extraction_keys)==0:
                continue
            while signal_extraction_keys.count(f"fInvMass_{formatted_eff_cut};1")==0 and eff_cut_increment<0.20:
                if eff_cut_sign == -1:
                    eff_cut_increment += 0.01
                eff_cut_sign *= -1
                formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin]+eff_cut_increment*eff_cut_sign)

            


            # if signal_extraction_keys.count(f"fInvMass_{formatted_eff_cut};1".encode())==0:
            #     continue
            
            # print('ehi')



            # get signal
            h_raw_yield = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fRawYields;1')
            h_histo = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fInvMass_{formatted_eff_cut}')
            h_histo.SetName(f'fInvMass_{formatted_eff_cut}_pt_{pt_bins[0]}_{pt_bins[1]}')
            #########################
            eff_index = h_raw_yield.FindBin(float(formatted_eff_cut))
            raw_yield = h_raw_yield.GetBinContent(eff_index)
            raw_yield_error = h_raw_yield.GetBinError(eff_index)

            dire = f'{split}_{cent_bins[0]}_{cent_bins[1]}'
            pt_spectra_file.cd()
            if not pt_spectra_file.GetListOfKeys().Contains(dire):
                pt_spectra_file.mkdir(dire)
            pt_spectra_file.cd(dire)
            h_histo.Write()


            presel_eff_map = np.logical_and(presel_eff_bin_centers > pt_bins[0], presel_eff_bin_centers < pt_bins[1])
            # absorption_map = np.logical_and(absorption_bin_centers > pt_bins[0], absorption_bin_centers < pt_bins[1])

            presel_eff = presel_eff_counts[presel_eff_map]
            absorption_corr = 0.97*HYP_HE_CROSS_SECT_SCALING
            bdt_eff = float(formatted_eff_cut)
            eff = presel_eff * eff_cut_dict[bin]

            pt_bin_index = h_corrected_yields[i_split].FindBin(pt_bins[0]+0.2)
            h_corrected_yields[i_split].SetBinContent(pt_bin_index, raw_yield/eff[0]/absorption_corr)
            h_corrected_yields[i_split].SetBinError(pt_bin_index, raw_yield_error/eff[0]/absorption_corr)
            print(pt_bin_index, raw_yield, eff[0])

            if i_split==0:
                h_corrected_ratio.SetBinContent(pt_bin_index, raw_yield/eff[0]/absorption_corr)
                h_corrected_ratio.SetBinError(pt_bin_index, raw_yield_error/eff[0]/absorption_corr)
            else:
                bin_content = h_corrected_ratio.GetBinContent(pt_bin_index)
                bin_error = h_corrected_ratio.GetBinError(pt_bin_index)
                # ratio_error = np.sqrt((bin_error/bin_content)**2 + (raw_yield_error/raw_yield)**2)
                h_corrected_ratio.SetBinContent(pt_bin_index, bin_content/(raw_yield/eff[0]/absorption_corr))
                # h_corrected_ratio.SetBinError(pt_bin_index, ratio_error)           
    

            print(f'bin: [{pt_bins[0]}, {pt_bins[1]}], BDT eff = {formatted_eff_cut}, presel_eff = {presel_eff}, raw_yield = {raw_yield}+{raw_yield_error}')


        # set label
        h_corrected_yields[i_split].GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
        h_corrected_yields[i_split].GetYaxis().SetTitle("d#it{N}/d(y#it{p}_{T})")

        for i_bin in range(len(bins))[1:]:
            bin_width = h_corrected_yields[i_split].GetBinWidth(i_bin)
            bin_content = h_corrected_yields[i_split].GetBinContent(i_bin)
            bin_error = h_corrected_yields[i_split].GetBinError(i_bin)
            h_corrected_yields[i_split].SetBinContent(i_bin, bin_content/bin_width/evts)
            h_corrected_yields[i_split].SetBinError(i_bin, bin_error/bin_width/evts)

        h_corrected_yields[i_split].SetMarkerStyle(20)
        h_corrected_yields[i_split].SetMarkerSize(0.8)

        if not pt_spectra_file.GetListOfKeys().Contains('results'):
            pt_spectra_file.mkdir('results')
        pt_spectra_file.cd('results')
        print(h_corrected_yields[i_split].GetBinContent(2))
        pwg = AliPWGFunc()
        histo,Integral, integral_error, bw_fit = hp.bw_fit(h_corrected_yields[i_split], bw, pwg)
        bw = -1
        print("yield: ", Integral, ", error: ", integral_error)
        histo.Write()
        bw_fit.SetName(f'func_{cent_bins[0]}_{cent_bins[1]}')
        bw_fit.Write()
    
    h_corrected_ratio.Write()

pt_spectra_file.Close()
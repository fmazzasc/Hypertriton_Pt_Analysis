import sys
sys.path.append('utils')
import helpers as hp

import os
import pickle
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm
import ROOT
import uproot
import yaml


bw_file = ROOT.TFile('utils/BlastWaveFits.root')
bw_file2 = ROOT.TFile("utils/Anti_fits.root")

def create_tf1_dictionary():
    tf1_dic = {}
    tf1_dic['0_10'] = bw_file.Get('BlastWave/BlastWave0')
    tf1_dic['10_30'] = bw_file2.Get(f"BGBW/3/BGBW3")
    tf1_dic['30_50'] = bw_file2.Get(f"BGBW/4/BGBW4")
    tf1_dic['50_90'] = bw_file2.Get(f"BGBW/5/BGBW5")
    return tf1_dic


tf1_dic = create_tf1_dictionary()



HYP_HE_CROSS_SECT_SCALING = 96/98


parser = argparse.ArgumentParser(prog='ml_analysis', allow_abbrev=True)
parser.add_argument('config', help='Path to the YAML configuration file')
args = parser.parse_args()

with open(os.path.expandvars(args.config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
##################################

ANALYSIS_RESULTS_PATH = params['ANALYSIS_RESULTS_PATH']
PT_BINS_CENT = params['PT_BINS_CENT']
CENTRALITY_LIST = params['CENTRALITY_LIST']
RANDOM_STATE = params['RANDOM_STATE']
MC_PATH = params['MC_SIGNAL_PATH']
MERGE_SAMPLES = params['MERGE_SAMPLES']
##################################################################

RESULTS_SUBDIR = params['RESULTS_SUBDIR']
res_dir = 'results' + RESULTS_SUBDIR
if not os.path.isdir(res_dir):
    os.mkdir(res_dir)
# split matter/antimatter
FIT = True
SPLIT=True
SPLIT_LIST = ['']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

#####################################################################
#####################################################################
analysis_results_file = uproot.open(os.path.expandvars(ANALYSIS_RESULTS_PATH))
# cent_counts, cent_edges = analysis_results_file["AliAnalysisTaskHyperTriton2He3piML_custom_summary;1"][11].to_numpy()  ##does not wprk with uproot3 
cent_counts = analysis_results_file["AliAnalysisTaskHyperTriton2He3piML_custom_summary;1"][11].values ##does not wprk with uproot4
if MERGE_SAMPLES:
    cent_counts += uproot.open('/data/fmazzasc/PbPb_2body/AnalysisResults_2015.root')["AliAnalysisTaskHyperTriton2He3piML_custom_summary;1"][11].values
cent_edges = analysis_results_file["AliAnalysisTaskHyperTriton2He3piML_custom_summary;1"][11].edges ##does not wprk with uproot4


######################### compute the preselection efficiency #######################
df_rec = uproot.open("/data/fmazzasc/PbPb_2body/MC_tables/SignalTable_flat_pt.root")["SignalTable"].pandas.df().query('abs(Rapidity)<0.5')
df_gen = uproot.open("/data/fmazzasc/PbPb_2body/MC_tables/SignalTable_flat_pt.root")["GenTable"].pandas.df().query('abs(rapidity)<0.5')
pt_fine_bins = np.linspace(2,9,1000)
hist_eff_pt_fine = hp.presel_eff_hist([df_rec, df_gen], 'pt', "fine_matter", [0,90], pt_fine_bins)
####################################################################################################


cent_bin_centers = (cent_edges[:-1]+cent_edges[1:])/2

print("Total number of events: ", np.sum(cent_counts))

eff_cut_dict = pickle.load(open(res_dir + "/file_eff_cut_dict", "rb"))
signal_extraction_file = ROOT.TFile.Open(res_dir + '/SignalExtraction.root')
signal_extraction_up = uproot.open(res_dir + '/SignalExtraction.root')
absorption_correction_file = uproot.open(res_dir + "/He3_abs.root")

pt_spectra_file = ROOT.TFile.Open(res_dir + '/single_bin_yield.root', 'recreate')
hist_eff_pt_fine.Write()

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
        print("---------------------------")
        print(f'{i_split} -> {split}')

        # get preselection efficiency
        print('**********computing presel eff***********')
        tf1_cent = tf1_dic[f'{cent_bins[0]}_{cent_bins[1]}']
        tf1_cent.Write()
        hist_eff_pt = hp.presel_eff_hist([df_rec, df_gen], 'pt', split, cent_bins, bins)
        hp.reweight_efficiency(hist_eff_pt, hist_eff_pt_fine, tf1_cent)
        f_missing = 1 - tf1_cent.Integral(0,2)/tf1_cent.Integral(0,10)
        print('PRESEL EFF: ', hist_eff_pt.GetBinContent(1))
        print('f_missing: ', f_missing)
        print('******************************************')
        presel_eff = hist_eff_pt.GetBinContent(1)*f_missing

        # get absorption correction
        func = "BlastWave" if cent_bins[0]<1 else "BGBW"
        absorption_counts = absorption_correction_file[f'{cent_bins[0]}_{cent_bins[1]}'][f'fEffPt_{split}_cent_{cent_bins[0]}_{cent_bins[1]}_func_{func};1'].values
        absorption_edges = absorption_correction_file[f'{cent_bins[0]}_{cent_bins[1]}'][f'fEffPt_{split}_cent_{cent_bins[0]}_{cent_bins[1]}_func_{func};1'].edges
        absorption_bin_centers = (absorption_edges[1:]+absorption_edges[:-1])/2


        h_corrected_yields[i_split] = ROOT.TH1D(f'fYields_{split}_{cent_bins[0]}_{cent_bins[1]}', f'{split}, {cent_bins[0]}-{cent_bins[1]}%', len(bins)-1, bins)
        for pt_bins in pt_bins_cent:
            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}'
            formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin])
            bkg_shape = 'pol1'
            eff_cut_increment = 0
            eff_cut_sign = -1
            signal_extraction_keys = signal_extraction_up[f"{bin}_{bkg_shape}"].keys()
            while signal_extraction_keys.count(f"fInvMass_{formatted_eff_cut};1".encode())==0 and eff_cut_increment<20:
                if eff_cut_sign == -1:
                    eff_cut_increment += 0.01
                eff_cut_sign *= -1
                formatted_eff_cut = "{:.2f}".format(eff_cut_dict[bin]+eff_cut_increment*eff_cut_sign)

            # get signal
            h_raw_yield = signal_extraction_file.Get(f'{bin}_{bkg_shape}/fRawYields;1')
            eff_index = h_raw_yield.FindBin(float(formatted_eff_cut))
            raw_yield = h_raw_yield.GetBinContent(eff_index)
            raw_yield_error = h_raw_yield.GetBinError(eff_index)

            absorption_map = np.logical_and(absorption_bin_centers > pt_bins[0], absorption_bin_centers < pt_bins[1])
            absorption_corr = absorption_counts[absorption_map]

            bdt_eff = float(formatted_eff_cut)
            eff = float(presel_eff * eff_cut_dict[bin])
            absorption_corr = np.mean(absorption_corr)*HYP_HE_CROSS_SECT_SCALING


            pt_bin_index = h_corrected_yields[i_split].FindBin(pt_bins[0]+0.5)
            h_corrected_yields[i_split].SetBinContent(pt_bin_index, raw_yield/eff/absorption_corr)
            h_corrected_yields[i_split].SetBinError(pt_bin_index, raw_yield_error/eff/absorption_corr)
            if SPLIT:
                if i_split==0:
                    h_corrected_ratio.SetBinContent(pt_bin_index, raw_yield/eff/absorption_corr)
                    h_corrected_ratio.SetBinError(pt_bin_index, raw_yield_error/eff/absorption_corr)
                else:
                    bin_content = h_corrected_ratio.GetBinContent(pt_bin_index)
                    bin_error = h_corrected_ratio.GetBinError(pt_bin_index)
                    ratio_error = np.sqrt((bin_error/bin_content)**2 + (raw_yield_error/raw_yield)**2)
                    h_corrected_ratio.SetBinContent(pt_bin_index, bin_content/(raw_yield/eff/absorption_corr))
                    h_corrected_ratio.SetBinError(pt_bin_index, ratio_error)           
    

            print(f'bin: [{pt_bins[0]}, {pt_bins[1]}], BDT eff = {formatted_eff_cut}, presel_eff = {presel_eff}, raw_yield = {raw_yield}+{raw_yield_error}')


        # set labels
        h_corrected_yields[i_split].GetXaxis().SetTitle("#it{p}_T (GeV/#it{c^2})")
        h_corrected_yields[i_split].GetYaxis().SetTitle("d#it{N}/d(#it{p}_T) (GeV/#it{c^2})")

        for i_bin in range(len(bins))[1:]:
            bin_width = h_corrected_yields[i_split].GetBinWidth(i_bin)
            bin_content = h_corrected_yields[i_split].GetBinContent(i_bin)
            bin_error = h_corrected_yields[i_split].GetBinError(i_bin)
            h_corrected_yields[i_split].SetBinContent(i_bin, bin_content/evts)
            h_corrected_yields[i_split].SetBinError(i_bin, bin_error/evts)


        h_corrected_yields[i_split].SetMarkerStyle(20)
        h_corrected_yields[i_split].SetMarkerSize(0.8)
        print("Corrected yield: ", h_corrected_yields)
        h_corrected_yields[i_split].Write()
    
    h_corrected_ratio.Write()

pt_spectra_file.Close()



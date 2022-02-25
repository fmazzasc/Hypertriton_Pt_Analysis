#!/usr/bin/env python3
from numbers import Integral
import os
import pickle
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import uproot
import yaml


def generate_th1(counts, bins, name=''):
    th1 = ROOT.TH1D(f'{name}', f'{name}', len(bins) - 1, np.array(bins, 'double'))
    for index in range(0, len(counts)):
        th1.SetBinContent(index+1, counts[index])
    th1.SetDirectory(0)
    return th1

def presel_eff_hist(df_list, col_name, split, cent_bins, bins):

    counts_rec = np.histogram(df_list[0][col_name], bins=bins)
    counts_gen = np.histogram(df_list[1]['gPt'], bins=bins)
    print('------------------------------------')
    print('bins: ', bins)
    eff = counts_rec[0]/counts_gen[0]
    print('eff: ', eff)

    print('------------------------------------')

    eff = counts_rec[0]/counts_gen[0]
    hist_eff = generate_th1(eff, bins, f"fPreselEff_vs_{col_name}_{split}_{cent_bins[0]}_{cent_bins[1]}")
    hist_eff.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    hist_eff.GetYaxis().SetTitle('Efficiency')
    hist_eff.SetMinimum(0)
    hist_eff.SetDrawOption("histo")
    hist_eff.SetLineWidth(2)

    # return histogram
    return hist_eff



def ndarray2roo(ndarray, var):
    if isinstance(ndarray, ROOT.RooDataSet):
        print('Already a RooDataSet')
        return ndarray

    assert isinstance(ndarray, np.ndarray), 'Did not receive NumPy array'
    assert len(ndarray.shape) == 1, 'Can only handle 1d array'

    name = var.GetName()
    x = np.zeros(1, dtype=np.float64)

    tree = ROOT.TTree('tree', 'tree')
    tree.Branch(f'{name}', x, f'{name}/D')

    for i in ndarray:
        x[0] = i
        tree.Fill()

    array_roo = ROOT.RooDataSet(
        'data', 'dataset from tree', tree, ROOT.RooArgSet(var))
    return array_roo


def significance_error(signal, background):
    num_s_prop = 0.5*signal+background
    num_b_prop = 0.5*signal
    den_prop = np.sqrt((signal+background)*(signal+background)*(signal+background))
    return 1/den_prop*np.sqrt(signal*num_s_prop*num_s_prop+background*num_b_prop*num_b_prop)

def expo(x):
    return np.exp(-x / (252 * 0.029979245800)) #hyp tau taken from lifetime analysis

def expected_signal_ct(cent_class, ct_range, pt_range, eff, n_events):
    he3_yield_list = [2.35e-4, 2.03e-4, 6.58e-5]
    correction = 0.4  # he3/hyp ratio (Very optimistic, considering it constant with centrality)
    correction *= 0.25 # 2-body Branching ratio
    correction *= expo(ct_range[0])- expo(ct_range[1]) #selecting the correct ct bin
    correction *= eff
    cent_end_bin = [5., 10., 50.]
    for cent_bin, he3_yield in zip(cent_end_bin, he3_yield_list):
        if cent_bin==cent_class[1]:
            return he3_yield*correction*n_events

    # expected signal for 0-90% centrality
    he3_yield_0_90 = 2.74e4
    return correction*he3_yield_0_90


def expected_signal_pt(cent_range, pt_range, eff, nevents):
    # print(nevents)

    bw_file = ROOT.TFile('utils/BlastWaveFits.root', 'read')
    bw = [bw_file.Get('BlastWave/BlastWave{}'.format(i)) for i in [0, 1, 2]]
    bw_file.Close()

    correction = 0.4  # Very optimistic, considering it constant with centrality
    correction *= 0.25 #B.R. correction
    cent_bins = [10, 40, 90]
    signal = 0

    for cent in range(cent_range[0]+1, cent_range[1]):
        for index in range(0, 3):
            if cent < cent_bins[index]:
                signal = signal + \
                    nevents[cent] * \
                    bw[index].Integral(pt_range[0], pt_range[1], 1e-8)
                # print(signal)
                break
    print(eff)
    return int(signal * eff * correction)




def bw_fit(histo, bw, pwg, fit_range=[2,9]):
    params = bw.GetParameters()
    print("BW: ", bw)
    bw_fit = pwg.GetBGBW(params[0], params[1], params[2], params[3], params[4])   
    bw_fit.FixParameter(0,2.991)
    bw_fit.FixParameter(1, params[1])
    bw_fit.SetParLimits(2, 0, 5)
    bw_fit.SetParLimits(3, 0, 5)
    bw_fit.SetParLimits(4, params[4]/20, params[4])

    fit_result = histo.Fit(bw_fit, "SI", "", fit_range[0], fit_range[1])
    cov_matrix = fit_result.GetCovarianceMatrix()
    integral = bw_fit.Integral(0,10, 1e-4)
    integral_error = bw_fit.IntegralError(0,10, fit_result.GetParams(), cov_matrix.GetMatrixArray(), 1e-1)

    return histo, integral, integral_error, bw_fit



def he3_bw_fit(histo, bw, pwg, fit_range=[2,9]):
    params = bw.GetParameters()
    print("BW: ", bw)
    bw_fit = pwg.GetBGBW(params[0], params[1], params[2], params[3], params[4])   
    bw_fit.FixParameter(0, params[0])
    bw_fit.FixParameter(1, params[1])
    bw_fit.FixParameter(2, params[2])
    bw_fit.FixParameter(3, params[3])
    bw_fit.SetParLimits(4, params[4]/20, params[4])

    # bw_fit.SetParLimits(4, 0, 1e3)

    fit_result = histo.Fit(bw_fit, "SI", "", fit_range[0], fit_range[1])
    cov_matrix = fit_result.GetCovarianceMatrix()
    integral = bw_fit.Integral(0,10, 1e-4)
    integral_error = bw_fit.IntegralError(0,10, fit_result.GetParams(), cov_matrix.GetMatrixArray(), 1e-1)
    # integral = 1
    # integral_error = 1


    return histo, integral, integral_error, bw_fit


def mtexpo_fit(histo):

    fit_function = ROOT.TF1('mtexpo', '[0]*x*TMath::Exp(-TMath::Sqrt(x**2 + 2.990**2)/[1])', 0, 10)
    fit_function.SetParameter(0, 2e-5)
    fit_function.SetParameter(1, 3)

    fit_function.SetParLimits(0, 1e-6, 1e-4)
    fit_function.SetParLimits(1, 1, 5)
    # print(fit_function.GetParLimits(0))

    fit_result = histo.Fit(fit_function, "SMI+")
    cov_matrix = fit_result.GetCovarianceMatrix()
    par_after = fit_function.GetParameters()
    integral = fit_function.Integral(0,20)
    integral_error = fit_function.IntegralError(0,20, fit_result.GetParams(), cov_matrix.GetMatrixArray())

    print("yield: ", integral, ", error: ", integral_error)
    return histo, integral, integral_error



def reweight_efficiency(hEfficiencyAnalysisBins, hEfficiencySmallBins, fInputSpectrum):  
    epsilon = 1e-5
    for ibin in range(1,hEfficiencyAnalysisBins.GetNbinsX()+1):
        c = hEfficiencyAnalysisBins.GetBinCenter(ibin)
        w = hEfficiencyAnalysisBins.GetBinWidth(ibin)/2.0
        x1 = c-w
        x2 = c+w
        n1 = hEfficiencySmallBins.FindBin(x1+epsilon)
        n2 = hEfficiencySmallBins.FindBin(x2-epsilon)
        num=0
        den=0
        deltaP_squared = 0
        deltaG_squared = 0

        for i in range(n1,n2+1):
            x = hEfficiencySmallBins.GetBinCenter(i)
            eps_i = hEfficiencySmallBins.GetBinContent(i)
            gen_i = fInputSpectrum.Eval(x)
            error_eps_i = hEfficiencySmallBins.GetBinError(i)
            error_gen_i = 0  
            deltaP_squared = deltaP_squared + (eps_i*error_gen_i)**2 + (error_eps_i*gen_i)**2
            deltaG_squared = deltaG_squared + error_gen_i*error_gen_i
            num = num + eps_i*gen_i
            den = den + gen_i

        average_eff = num/den if den!=0 else 0   

        deltaP = np.sqrt(deltaP_squared)
        deltaG = np.sqrt(deltaG_squared)
        error_eff = np.sqrt((deltaP/den)**2 + (num*deltaG/(den*den))**2)

        ##Fill Average Efficiency
        hEfficiencyAnalysisBins.SetBinContent(ibin,average_eff)
        hEfficiencyAnalysisBins.SetBinError(ibin,error_eff)


def apply_pt_rejection(df, pt_shape):
    rej_flag = np.ones(len(df))
    random_arr = np.random.rand(len(df))
    max_bw = pt_shape.GetMaximum()

    for ind, (pt, rand) in enumerate(zip(df['gPt'],random_arr)):
        frac = pt_shape.Eval(pt)/max_bw
        if rand > frac:
            rej_flag[ind] = -1
    df._full_data_frame['rej'] = rej_flag



def get_number_of_MB_ev(cent_class, an_file):
    histo = an_file["AliAnalysisTaskHyperTriton2He3piML_custom_summary;1"][14].to_numpy()

    cent_bin_centers = (histo[1][:-1]+ histo[1][1:])/2
    cent_range_mask = np.logical_and(cent_bin_centers > cent_class[0], cent_bin_centers < cent_class[1])

    trigger_bin_centers = (histo[2][:-1]+ histo[2][1:])/2
    trigger_val = np.floor(trigger_bin_centers)
    trigger_mask = trigger_val%2 != 0

    filtered_histo = histo[0][cent_range_mask] #apply cent mask
    filtered_histo = filtered_histo[:,trigger_mask] #apply trigger mask

    return np.sum(filtered_histo)



def get_pt_shape_cent(cent_class):

    cent_bins_new =  [[0,5],[5,10],[10,30],[30,50],[50,90]]
    cent_bins_MB = [[0, 10], [10, 40], [40, 90]]
    # functions
    func = {}
    func_max = {}

    func_names = ["BGBW", "Boltzmann", "Mt-exp", "Pt-exp", "LevyTsallis"]
    func_names_MB = ["BlastWave", "Boltzmann", "LevyTsallis", "Mt-exp"]

    # functions input files
    input_func_file = ROOT.TFile("/home/fmazzasc/Hypertriton/Hypertriton_Pt_Analysis/utils/Anti_fits.root")
    input_func_file_MB = ROOT.TFile("/home/fmazzasc/Hypertriton/Hypertriton_Pt_Analysis/utils/BlastWaveFits.root")

    # get functions and maxima from file

    for i_cent in range(len(cent_bins_new)):
        for i_fun in range(len(func_names)):
            cent_bins = cent_bins_new[i_cent]
            key = f"cent_{cent_bins[0]}_{cent_bins[1]}_func_{func_names[i_fun]}"
            func[key] = input_func_file.Get(f"{func_names[i_fun]}/{i_cent + 1}/{func_names[i_fun]}{i_cent + 1}")
            func_max[key] = func[key].GetMaximum()

    for i_cent in range(len(cent_bins_MB)):
        for i_fun in range(len(func_names_MB)):
            cent_bins = cent_bins_MB[i_cent]
            key = f"cent_{cent_bins[0]}_{cent_bins[1]}_func_{func_names_MB[i_fun]}"
            func[key] = input_func_file_MB.Get(f"{func_names_MB[i_fun]}/{func_names_MB[i_fun]}{i_cent}")
            func_max[key] = func[key].GetMaximum()

    return func[f"cent_{cent_class[0]}_{cent_class[1]}_func_BGBW"]
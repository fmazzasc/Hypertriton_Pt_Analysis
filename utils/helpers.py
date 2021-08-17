#!/usr/bin/env python3
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
    counts_gen = np.histogram(df_list[1][col_name], bins=bins)
    print('------------------------------------')
    print(bins)
    print(counts_rec, counts_gen)
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




def bw_fit(histo, bw):
    params = bw.GetParameters()
    # params[0] = 2.991
    pwg = ROOT.AliPWGFunc()
    bw = pwg.GetBGBW(params[0], params[1], params[2], params[3], params[4])
    bw.FixParameter(0,2.991)
    bw.SetParLimits(1, 0, 2)
    bw.SetParLimits(2, 0, 2)
    bw.SetParLimits(3, 0, 2)
    bw.SetParLimits(4, 0, 1000)

    fit_result = histo.Fit(bw, "QSMI+", "", 0,100)
    cov_matrix = fit_result.GetCovarianceMatrix()
    integral = bw.Integral(0,20, 1e6)
    integral_error = bw.IntegralError(0,20, fit_result.GetParams(), cov_matrix.GetMatrixArray())

    return histo, integral, integral_error, bw
    


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


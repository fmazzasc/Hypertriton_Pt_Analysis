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
    bw.SetParLimits(0, 2.989, 2.992)
    bw.SetParLimits(1, 0, 2)
    bw.SetParLimits(2, 0, 2)
    bw.SetParLimits(3, 0, 2)
    bw.SetParLimits(4, 0, 1)

    fit_result = histo.Fit(bw, "SMI+")
    cov_matrix = fit_result.GetCovarianceMatrix()
    par_after = bw.GetParameters()
    integral = bw.Integral(0,20)
    print(cov_matrix[1][1])
    integral_error = bw.IntegralError(0,20, fit_result.GetParams(), cov_matrix.GetMatrixArray())

    print("yield: ", integral, ", error: ", integral_error)
    return histo, integral,integral_error
    

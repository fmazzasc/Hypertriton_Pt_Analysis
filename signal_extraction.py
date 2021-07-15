#!/usr/bin/env python3
import sys
sys.path.append('utils')

import os
import pickle
import warnings
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import uproot
import yaml
from helpers import significance_error, ndarray2roo

SPLIT = True
MAX_EFF = 1.00

# avoid pandas warning
warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.gROOT.SetBatch()

parser = argparse.ArgumentParser(prog='signal_extraction', allow_abbrev=True)
parser.add_argument('-bkgExpo', action='store_true')
parser.add_argument('-bkgPol2', action='store_true')

args = parser.parse_args()

BKG_EXPO = args.bkgExpo
BKG_POL2 =  args.bkgPol2

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
PT_BINS_CENT = params['PT_BINS_CENT']
CENTRALITY_LIST = params['CENTRALITY_LIST']
RANDOM_STATE = params['RANDOM_STATE']
##################################################################

# split matter/antimatter
SPLIT_LIST = ['all']
if SPLIT:
    SPLIT_LIST = ['antimatter', 'matter']

bkg_shape = 'pol1'
if BKG_EXPO:
    bkg_shape = 'expo'
if BKG_POL2:
    bkg_shape = 'pol2'


score_eff_arrays_dict = pickle.load(open("results/file_score_eff_dict", "rb"))
eff_array = np.arange(0.10, MAX_EFF, 0.01)

root_file_signal_extraction = ROOT.TFile("results/SignalExtraction.root", "recreate")
for split in SPLIT_LIST:
    for i_cent_bins, pt_bins_cent in enumerate(PT_BINS_CENT):
        cent_bins = CENTRALITY_LIST[i_cent_bins]
        for pt_bins in pt_bins_cent:
            bin = f'{split}_{cent_bins[0]}_{cent_bins[1]}_{pt_bins[0]}_{pt_bins[1]}'
            df_data = pd.read_parquet(f'df/{bin}')
            df_signal = pd.read_parquet(f'df/mc_{bin}')

            # ROOT.Math.MinimizerOptions.SetDefaultTolerance(1e-2)
            root_file_signal_extraction.mkdir(f'{bin}_{bkg_shape}')

            # raw yileds histogram
            h_raw_yields = ROOT.TH1D("fRawYields", "fRawYields", 101, -0.005, 1.005)

            # significance histogram
            h_significance = ROOT.TH1D("fSignificance", "fSignificance", 101, -0.005, 1.005)

            for eff_score in zip(eff_array, score_eff_arrays_dict[bin]):
                if eff_score[0] < 0.50:
                    continue
                formatted_eff = "{:.2f}".format(eff_score[0])
                print(f'processing {bin}: eff = {eff_score[0]:.2f}, score = {eff_score[1]:.2f}...')

                df_data_sel = df_data.query(f'model_output > {eff_score[1]}')
                df_signal_sel = df_signal.query(f'model_output > {eff_score[1]}')
                if len(df_signal_sel) > 10000:
                    print('Sampling 10000 events...')
                    df_signal_sel = df_signal_sel.sample(10000)

                # get invariant mass distribution (data and mc)
                roo_m = ROOT.RooRealVar("m", "#it{M} (^{3}He + #pi^{-})", 2.960, 3.025, "GeV/#it{c}^{2}")
                roo_data = ndarray2roo(np.array(df_data_sel['m']), roo_m)
                roo_mc_signal = ndarray2roo(np.array(df_signal_sel['m']), roo_m)

                # declare fit model
                # kde
                roo_n_signal = ROOT.RooRealVar('N_{signal}', 'Nsignal', 0., 1.e3)
                delta_mass = ROOT.RooRealVar("#Deltam", 'deltaM', -0.004, 0.004, 'GeV/c^{2}')
                shifted_mass = ROOT.RooAddition("mPrime", "m + #Deltam", ROOT.RooArgList(roo_m, delta_mass))
                roo_signal = ROOT.RooKeysPdf("signal", "signal", shifted_mass, roo_m,
                                             roo_mc_signal, ROOT.RooKeysPdf.NoMirror, 2)
                roo_signal_plot = ROOT.RooKeysPdf(roo_signal)

                # background
                roo_n_background = ROOT.RooRealVar('N_{bkg}', 'Nbackground', 0., 1.e4)
                roo_slope = ROOT.RooRealVar('slope', 'slope', -20., 20.)
                roo_slope2 = ROOT.RooRealVar('slope2', 'slope2', -20., 20.)


                roo_bkg = ROOT.RooRealVar()

                if bkg_shape=="pol1":
                    roo_bkg = ROOT.RooPolynomial('background', 'background', roo_m, ROOT.RooArgList(roo_slope))
                
                elif bkg_shape=="pol2":
                    roo_bkg = ROOT.RooPolynomial('background', 'background', roo_m, ROOT.RooArgList(roo_slope, roo_slope2))

                else:
                    roo_bkg = ROOT.RooExponential('background', 'background', roo_m, roo_slope)

                # model
                roo_model = ROOT.RooAddPdf(
                    'model', 'model', ROOT.RooArgList(roo_signal, roo_bkg),
                    ROOT.RooArgList(roo_n_signal, roo_n_background))

                # fit
                ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
                ROOT.RooMsgService.instance().setSilentMode(ROOT.kTRUE)
                ROOT.gErrorIgnoreLevel = ROOT.kError
                r = roo_model.fitTo(roo_data, ROOT.RooFit.Save(), ROOT.RooFit.Extended(ROOT.kTRUE))

                print(f'fit status: {r.status()}')
                if r.status() == 0 and delta_mass.getError() > 1.e-6:

                    # plot
                    nBins = 26
                    xframe = roo_m.frame(2.96, 3.025, nBins)
                    xframe.SetTitle(
                        str(pt_bins[0]) + '#leq #it{c}t<' + str(pt_bins[1]) + ' cm, ' + str(cent_bins[0]) + '-' +
                        str(cent_bins[1]) + '%, BDT efficiency = ' + str(formatted_eff))
                    xframe.SetName(f'fInvMass_{formatted_eff}')
                    roo_data.plotOn(xframe, ROOT.RooFit.Name('data'))
                    roo_model.plotOn(
                        xframe, ROOT.RooFit.Components('background'),
                        ROOT.RooFit.Name('background'),
                        ROOT.RooFit.LineStyle(ROOT.kDashed),
                        ROOT.RooFit.LineColor(ROOT.kGreen))
                    roo_model.plotOn(xframe, ROOT.RooFit.Components('signal'), ROOT.RooFit.Name('signal'),
                                     ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))
                    roo_model.plotOn(xframe, ROOT.RooFit.Name('model'), ROOT.RooFit.LineColor(ROOT.kBlue))

                    formatted_chi2 = "{:.2f}".format(xframe.chiSquare('model', 'data'))
                    roo_model.paramOn(xframe, ROOT.RooFit.Label(
                        '#chi^{2}/NDF = '+formatted_chi2),
                        ROOT.RooFit.Layout(0.55, 0.85, 0.88))
                    xframe.getAttText().SetTextSize(0.035)
                    xframe.getAttLine().SetLineWidth(0)

                    print(f'chi2/NDF: {formatted_chi2}, edm: {r.edm()}')
                    if float(formatted_chi2) < 2 and r.edm() < 1:

                        # fit mc distribution to get sigma and mass
                        roo_mean_mc = ROOT.RooRealVar("mean", "mean", 2.98, 3.0)
                        roo_sigma_mc = ROOT.RooRealVar("sigma", "sigma", 0.0005, 0.0040)
                        gaus = ROOT.RooGaussian('gaus', 'gaus', roo_m, roo_mean_mc, roo_sigma_mc)
                        gaus.fitTo(roo_mc_signal)

                        # mass
                        mass_val = roo_mean_mc.getVal()-delta_mass.getVal()

                        # significance
                        m_set = ROOT.RooArgSet(roo_m)
                        normSet = ROOT.RooFit.NormSet(m_set)
                        roo_m.setRange(
                            'signalRange', mass_val - 3 * roo_sigma_mc.getVal(),
                            mass_val + 3 * roo_sigma_mc.getVal())
                        signal_int = (roo_model.pdfList().at(0).createIntegral(
                            m_set, normSet, ROOT.RooFit.Range("signalRange"))).getVal()
                        print(f'signal integral = {signal_int}')
                        bkg_int = (roo_model.pdfList().at(1).createIntegral(
                            m_set, normSet, ROOT.RooFit.Range("signalRange"))).getVal()
                        print(f'background integral = {bkg_int}')
                        sig = signal_int*roo_n_signal.getVal()
                        bkg = bkg_int*roo_n_background.getVal()
                        significance_val = sig/np.sqrt(sig+bkg)
                        significance_err = significance_error(sig, bkg)

                        # fill significance histogram
                        eff_index = h_significance.FindBin(float(formatted_eff))
                        h_significance.SetBinContent(eff_index, significance_val)
                        h_significance.SetBinError(eff_index, significance_err)

                        if significance_val > 2:
                            # fill raw yields histogram
                            h_raw_yields.SetBinContent(eff_index, roo_n_signal.getVal())
                            h_raw_yields.SetBinError(eff_index, roo_n_signal.getError())

                            # write to file
                            root_file_signal_extraction.cd(f'{bin}_{bkg_shape}')
                            xframe.Write()

                            # draw on canvas and save plots
                            canv = ROOT.TCanvas()
                            canv.cd()
                            text_mass = ROOT.TLatex(
                                2.965, 0.74 * xframe.GetMaximum(),
                                "#it{m}_{^{3}_{#Lambda}H} = " + "{:.6f}".format(mass_val) + " GeV/#it{c^{2}}")
                            text_mass.SetTextSize(0.035)
                            text_signif = ROOT.TLatex(2.965, 0.91 * xframe.GetMaximum(),
                                                    "S/#sqrt{S+B} (3#sigma) = " + "{:.3f}".format(significance_val) + " #pm " +
                                                    "{:.3f}".format(significance_err))
                            text_signif.SetTextSize(0.035)
                            text_sig = ROOT.TLatex(2.965, 0.84 * xframe.GetMaximum(), "S (3#sigma) = " + "{:.1f}".format(sig) + " #pm " + "{:.1f}".format(signal_int*roo_n_signal.getError()))
                            text_sig.SetTextSize(0.035)
                            text_bkg = ROOT.TLatex(2.965, 0.77 * xframe.GetMaximum(), "B (3#sigma) = " + "{:.1f}".format(bkg) + " #pm" + "{:.1f}".format(bkg_int*roo_n_background.getError()))
                            text_bkg.SetTextSize(0.035)
                            xframe.Draw("")
                            # text_mass.Draw("same")
                            text_signif.Draw("same")
                            text_sig.Draw("same")
                            text_bkg.Draw("same")
                            print(
                                f'significance = {"{:.3f}".format(significance_val)} +/- {"{:.3f}".format(significance_err)}')
                            if not os.path.isdir('plots/signal_extraction'):
                                os.mkdir('plots/signal_extraction')
                            if not os.path.isdir(f'plots/signal_extraction/{bin}_{bkg_shape}'):
                                os.mkdir(f'plots/signal_extraction/{bin}_{bkg_shape}')
                            canv.Print(f'plots/signal_extraction/{bin}_{bkg_shape}/{eff_score[0]:.2f}_{bin}.png')

                            # plot kde and mc
                            frame = roo_m.frame(2.96, 3.025, 130)
                            frame.SetTitle(str(cent_bins[0])+"-"+str(cent_bins[1])+"%, "+str(pt_bins[0])+"#leq #it{c}t<"+str(pt_bins[1])+" cm, BDT efficiency = "+str(formatted_eff))
                            roo_mc_signal.plotOn(frame)
                            roo_signal_plot.plotOn(frame, ROOT.RooFit.Name("KDE"))
                            gaus.plotOn(frame, ROOT.RooFit.Name("gaussian"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed))
                            cc = ROOT.TCanvas("cc", "cc")
                            if not os.path.isdir('plots/kde_signal'):
                                os.mkdir('plots/kde_signal')
                            if not os.path.isdir(f'plots/kde_signal/{bin}'):
                                os.mkdir(f'plots/kde_signal/{bin}')
                            frame.Draw()
                            leg_mc = ROOT.TLegend(0.6, 0.8, 0.85, 0.7)
                            leg_mc.AddEntry(frame.findObject("KDE"), "KDE")
                            leg_mc.AddEntry(frame.findObject("gaussian"), "Gaussian")
                            leg_mc.SetBorderSize(0)
                            leg_mc.Draw("same")
                            cc.SetLogy(ROOT.kTRUE)
                            cc.Print(f'plots/kde_signal/{bin}/{formatted_eff}_{bin}.png')

            h_raw_yields.GetXaxis().SetTitle("BDT efficiency")
            h_raw_yields.GetYaxis().SetTitle("#it{N_{raw}}")
            h_raw_yields.Write()

            h_significance.GetXaxis().SetTitle("BDT efficiency")
            h_significance.GetYaxis().SetTitle("S / #sqrt{S + B}")
            h_significance.Write()

root_file_signal_extraction.Close()

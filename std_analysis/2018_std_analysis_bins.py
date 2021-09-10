import sys
sys.path.append('../utils')
import helpers as hp

import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
import helpers_std as hp_std
import ROOT
from hipe4ml.tree_handler import TreeHandler
ROOT.gROOT.SetBatch()

ROOT.gROOT.ProcessLine(".L ../utils/alipwgfunc/AliPWGFunc.cxx++")

from ROOT import AliPWGFunc


isAOD = False


def rename_mc_df_columns(df):
    rename_dict = {}
    rename_dict['ptMC'] = 'gPt'
    rename_dict['etaMC'] = 'gEta'
    rename_dict['ctMC'] = 'gCt'
    rename_dict['yMC'] = 'gRapidity'
    
    df.rename(columns = rename_dict, inplace=True)

def apply_pt_rejection(df, pt_shape):
    rej_flag = -1*np.ones(len(df))
    random_arr = np.random.rand(len(df))
    max_bw = pt_shape.GetMaximum()

    for ind, (pt, rand) in enumerate(zip(df['gPt'],random_arr)):
        frac = pt_shape.Eval(pt)/max_bw
        if rand < frac:
            rej_flag[ind] = True
    df['rej'] = rej_flag


cosPAcuts = np.linspace(0.9999, 0.9999, 1)
pidHe3cuts = np.linspace(100, 121, 1)
he3PtCuts = np.linspace(1.7, 1.9, 1)
piPtCuts = np.linspace(0.15, 0.18, 1)
prongsDCA = np.linspace(1,2, 1)

cent_class = [0,10]
n_ev = 98e6

mass_range = [2.96,3.04]
n_bins = 40

if isAOD:
    df_mc = uproot.open("/data/fmazzasc/PbPb_2body/AOD/HyperTritonTree_MC.root")["HyperTree"].arrays(library="pd")
    rename_mc_df_columns(df_mc)
    df = pd.read_parquet("/data/fmazzasc/PbPb_2body/AOD/HyperTritonTree_18_red.parquet.gzip")


else:
    hdl = TreeHandler('/data/fmazzasc/PbPb_2body/SignalTable_20g7_flat_pt.root', 'SignalTable')
    df_mc = hdl.get_data_frame()
    df = uproot.open("/data/fmazzasc/PbPb_2body/old_stuff/DataTable_18qr_pass3.root")["DataTable"].arrays(library="pd")


# print(np.min(df.query('pt>0')))

bw_file = ROOT.TFile('../utils/BlastWaveFits.root')
bw = bw_file.Get("BlastWave/BlastWave0")

apply_pt_rejection(df_mc, bw)




pt_bins = [[2,3],[3,3.5],[3.5,4],[4,4.5], [4.5,5], [5,6], [6,9]]
# pt_bins = [[0,9]]

flat_pt_bins = [item for sublist in pt_bins for item in sublist]
bins = np.unique(np.array(flat_pt_bins, dtype=float))



if isAOD:
    df_rec = df_mc.query('rej>0 and isReconstructed==True and abs(Rapidity)<0.5')
    df_gen = df_mc.query('rej>0 and abs(gRapidity)<0.5')
    ffile = ROOT.TFile(f"std_analysis_2018_{cent_class[0]}_{cent_class[1]}_AOD.root", "recreate")
else:
    ffile = ROOT.TFile(f"std_analysis_2018_{cent_class[0]}_{cent_class[1]}.root", "recreate")
    df_rec = df_mc.query('rej>0 and pt>0 and abs(Rapidity)<0.5')
    df_gen = df_mc.query('rej>0 and abs(gRapidity)<0.5')


ffile.mkdir(f'{cent_class[0]}_{cent_class[1]}')


pt_spectrum =  ROOT.TH1D(f'Yield', 'Yield', len(bins)-1, bins)
efficiency =  ROOT.TH1D(f'Efficiency', 'Efficiency', len(bins)-1, bins)
print(pt_spectrum)


for ind,pt_bin in enumerate(pt_bins):
    ffile.cd()
    for cosPA in cosPAcuts:
        for pidHe3 in pidHe3cuts:
            for he3Pt in he3PtCuts:
                for piPt in piPtCuts:
                    for pDCA in prongsDCA:
                        print(pt_bin)
                        print('********************************************************************')
                        cut = f'2<ct<35 and V0CosPA > {cosPA} and NpidClustersHe3 > {pidHe3} and pt > {pt_bin[0]} and pt < {pt_bin[1]} and He3ProngPvDCA > 0.05 and PiProngPvDCA > 0.2 and abs(TPCnSigmaHe3) < 3.5 and ProngsDCA < {pDCA}'
                        cut_cent = cut + f" and {cent_class[0]}<centrality<{cent_class[1]}"
                        print("#################################")
                        print("CUT: ", cut_cent)
                        print("#################################")
                        print("Efficiency before cuts")
                        print(len(df_rec.query(f'pt > {pt_bin[0]} and pt < {pt_bin[1]}'))/len(df_gen.query(f'gPt > {pt_bin[0]} and gPt < {pt_bin[1]}')))
                        print("#################################")
                        print("Efficiency after cuts")
                        eff = len(df_rec.query(cut))/len(df_gen.query(f'gPt > {pt_bin[0]} and gPt < {pt_bin[1]}'))
                        print(eff)

                        df_sel = df.query(cut_cent)
                        selected_data_hist = np.histogram(np.array(df_sel['m']), bins=n_bins, range=mass_range)                        
                        selected_data_hist = hp_std.h1_invmass(selected_data_hist[0], mass_range=mass_range, bins=n_bins, name=f'{cut}')
                        fit_result = hp_std.fit_hist(selected_data_hist, cent_class =cent_class, pt_range = pt_bin, ct_range = [2,35])
                        bin_width = pt_bin[1] - pt_bin[0]
                        pt_spectrum.SetBinContent(ind + 1, fit_result[0]/eff/n_ev/bin_width/2)
                        pt_spectrum.SetBinError(ind + 1, fit_result[1]/eff/n_ev/bin_width/2)
                        efficiency.SetBinContent(ind + 1, eff)
pwg = AliPWGFunc()
histo,Integral, integral_error, bw_fit = hp.bw_fit(pt_spectrum, bw, pwg)
print("yield: ", Integral, ", error: ", integral_error)
histo.Write()
bw_fit.SetName(f'func_{cent_class[0]}_{cent_class[1]}')
bw_fit.Write()
pt_spectrum.Write()
efficiency.Write()



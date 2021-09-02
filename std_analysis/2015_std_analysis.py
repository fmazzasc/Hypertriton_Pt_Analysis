import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
import helpers_std as hp
import ROOT
from hipe4ml.tree_handler import TreeHandler
ROOT.gROOT.SetBatch()



def apply_pt_rejection(df, pt_shape):
    rej_flag = -1*np.ones(len(df))
    random_arr = np.random.rand(len(df))
    max_bw = pt_shape.GetMaximum()

    for ind, (pt, rand) in enumerate(zip(df['gPt'],random_arr)):
        frac = pt_shape.Eval(pt)/max_bw
        if rand < frac:
            rej_flag[ind] = True
    df['rej'] = rej_flag


cosPAcuts = np.linspace(0.9995, 0.9998, 1)
pidHe3cuts = np.linspace(120, 121, 1)
he3PtCuts = np.linspace(1.7, 1.9, 1)
piPtCuts = np.linspace(0.15, 0.18, 1)
prongsDCA = np.linspace(1,2, 1)

cent_classes = [[0,10], [10,30], [30,50]]
n_ev = [99e6, 11e6*2, 11e6*2]

mass_range = [2.96,3.04]
n_bins = 40

hdl = TreeHandler('/data/fmazzasc/PbPb_2body/SignalTable_20g7_flat_pt.root', 'SignalTable')
df = hdl.get_data_frame()
print(np.min(df.query('pt>0')))

bw_file = ROOT.TFile('../utils/BlastWaveFits.root')
bw = bw_file.Get("BlastWave/BlastWave0")

apply_pt_rejection(df, bw)

presel_eff = len(df.query('rej>0 and pt>0 and abs(Rapidity)<0.5'))/len(df.query('rej>0 and abs(gRapidity)<0.5'))
df_rec = df.query('rej>0 and pt>0 and abs(Rapidity)<0.5')


df = uproot.open("/data/fmazzasc/PbPb_2body/old_stuff/DataTable_18qr_pass3.root")["DataTable"].arrays(library="pd")
ffile = ROOT.TFile("std_analysis.root", "recreate")

for ind, cent in enumerate(cent_classes):
    ffile.mkdir(f'{cent[0]}_{cent[1]}')
    ffile.cd()
    for cosPA in cosPAcuts:
        for pidHe3 in pidHe3cuts:
            for he3Pt in he3PtCuts:
                for piPt in piPtCuts:
                    for pDCA in prongsDCA:
                        print(cosPA)
                        print('********************************************************************')
                        cut = ('2<ct<35 and V0CosPA > {} and NpidClustersHe3 > {} and He3ProngPt > {} and pt > 2 and pt < 9 and PiProngPt > {} and He3ProngPvDCA > 0.05 and PiProngPvDCA > 0.2 and abs(TPCnSigmaHe3) < 3.5 and ProngsDCA < {}'.format(cosPA, pidHe3, he3Pt, piPt, pDCA))
                        cut_cent = cut + f" and {cent[0]}<centrality<{cent[1]}"
                        print("#################################")
                        print("CUT: ", cut_cent)
                        print("#################################")
                        print("Efficiency before cuts")
                        print(presel_eff)
                        print("#################################")
                        print("Efficiency after cuts")
                        eff = presel_eff*len(df_rec.query(cut))/len(df_rec)
                        print(eff)

                        df_sel = df.query(cut_cent)
                        selected_data_hist = np.histogram(np.array(df_sel['m']), bins=n_bins, range=mass_range)                        
                        selected_data_hist = hp.h1_invmass(selected_data_hist[0], mass_range=mass_range, bins=n_bins, name=f'{cut}')
                        fit_result = hp.fit_hist(selected_data_hist, cent_class =cent, pt_range = [2,9], ct_range = [2,35])
                        print(fit_result)
                        print('yield: ', fit_result[0]/eff/n_ev[ind]/2) #matter + antimatter /2


ffile.Close()
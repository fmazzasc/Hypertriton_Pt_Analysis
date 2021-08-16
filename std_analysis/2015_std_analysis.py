import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt

df = uproot.open("/data/fmazzasc/PbPb_2body/DataTable_15o.root")["DataTable"].arrays(library="pd")
# print(df)

df = df.query('Matter==1 and 0<centrality<10 and V0CosPA > 0.9995 and He3ProngPt > 1.8 and pt > 2 and pt < 10 and PiProngPt > 0.15 and He3ProngPvDCA > 0.05 and PiProngPvDCA > 0.2 and abs(TPCnSigmaHe3) < 3.5 and ProngsDCA < 1')

plt.hist(df["m"], bins=40, range=(2.96,3.04))
plt.savefig("inv.png")







import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
import helpers_std as hp
import ROOT
ROOT.gROOT.SetBatch()

cosPAcuts = np.linspace(0.9996, 0.99995, 1)
pidHe3cuts = np.linspace(80, 101, 1)
he3PtCuts = np.linspace(1.7, 1.9, 1)
piPtCuts = np.linspace(0.15, 0.18, 1)
prongsDCA = np.linspace(1,2, 1)

cent_classes = [[0,10], [10,30], [30,50]]
n_ev = [11e6, 11e6*2, 11e6*2]

mass_range = [2.96,3.04]
n_bins = 40


mc_file = uproot.open("/data/fmazzasc/PbPb_2body/SignalTable_16h7abc.root")
df = uproot.open("/data/fmazzasc/PbPb_2body/DataTable_15o.root")["DataTable"].arrays(library="pd")
rec_df = mc_file["SignalTable"].arrays(library="pd")
gen_df = mc_file["GenTable"].arrays(library="pd")

ffile = ROOT.TFile("std_analysis.root", "recreate")

for ind, cent in enumerate(cent_classes):
    ffile.mkdir(f'{cent[0]}_{cent[1]}')
    ffile.cd()
    for cosPA in cosPAcuts:
        for pidHe3 in pidHe3cuts:
            for he3Pt in he3PtCuts:
                for piPt in piPtCuts:
                    for pDCA in prongsDCA:
                        
                        cut = ('Matter==1 and 2<ct<35 and V0CosPA > {} and NpidClustersHe3 > {} and He3ProngPt > {} and pt > 2 and pt < 9 and PiProngPt > {} and He3ProngPvDCA > 0.05 and PiProngPvDCA > 0.2 and abs(TPCnSigmaHe3) < 3.5 and ProngsDCA < {}'.format(
                                cosPA, pidHe3, he3Pt, piPt, pDCA))
                        cut = cut + f" and {cent[0]}<centrality<{cent[1]}"
                        print("#################################")
                        print("CUT: ", cut)
                        eff = len(rec_df.query(cut + "and abs(Rapidity)<0.5"))/len(gen_df.query("abs(rapidity)<0.5" + f" and {cent[0]}<centrality<{cent[1]} and matter==1"))
                        print("EFFICIENCY: ", eff)
                        df_sel = df.query(cut)
                        selected_data_hist = np.histogram(np.array(df_sel['m']), bins=n_bins, range=mass_range)                        
                        selected_data_hist = hp.h1_invmass(selected_data_hist[0], mass_range=mass_range, bins=n_bins, name=f'{cut}')
                        fit_result = hp.fit_hist(selected_data_hist, cent_class =cent, pt_range = [2,9], ct_range = [2,35])
                        print(fit_result)
                        print('yield: ', fit_result[0]/eff/n_ev[ind])


ffile.Close()
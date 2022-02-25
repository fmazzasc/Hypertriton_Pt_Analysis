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

#######################################################################
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
#######################################################################
cent_class = [0, 10]

isAOD = True
is2015 = False
minimumBias = False
is_2018_and_2015 = False

cosPAcuts = np.linspace(0.9999, 0.99995, 1)
pidHe3cuts = np.linspace(110, 101, 1)
he3PtCuts = np.linspace(1.7, 1.9, 1)
piPtCuts = np.linspace(0.15, 0.18, 1)
prongsDCA = np.linspace(1,2, 1)

absorption_corr = 0.98
mass_range = [2.96,3.04]
n_bins = 40

if isAOD:
    if is2015:
        print('No AOD tree for 2015 sample')
        exit()

    df_mc = uproot.open("/data/fmazzasc/PbPb_2body/AOD/HyperTritonTree_MC.root")["HyperTree"].arrays(library="pd")
    rename_mc_df_columns(df_mc)
    df = pd.read_parquet("/data/fmazzasc/PbPb_2body/AOD/HyperTritonTree_15_red.parquet.gzip")
    an_file = uproot.open("/data/fmazzasc/PbPb_2body/AOD/HyperTritonTree_2015_2.root")

else:
    if is2015:
        hdl = TreeHandler('/data/fmazzasc/PbPb_2body/2015/SignalTable_16h7abc_flat_pt.root', 'SignalTable')
        df_mc = hdl.get_data_frame()
        df = uproot.open("/data/fmazzasc/PbPb_2body/2015/DataTable_15o.root")["DataTable"].arrays(library="pd")
        an_file = uproot.open('/data/fmazzasc/PbPb_2body/2015/AnalysisResults_15o.root')
    else:
        hdl = TreeHandler('/data/fmazzasc/PbPb_2body/2018/SignalTable_20g7_flat_pt.root', 'SignalTable')
        df_mc = hdl.get_data_frame()
        df = pd.read_parquet('/data/fmazzasc/PbPb_2body/2018/DataTable_18qr.parquet.gzip') # df = uproot.open("/data/fmazzasc/PbPb_2body/old_stuff/DataTable_18qr_pass3.root")["DataTable"].arrays(library="pd")
        an_file = uproot.open("/data/fmazzasc/PbPb_2body/2018/AnalysisResults_18qr.root")

        if is_2018_and_2015:
            df_2015 = uproot.open("/data/fmazzasc/PbPb_2body/DataTable_15o.root")["DataTable"].arrays(library="pd")
            an_2015 = uproot.open('/data/fmazzasc/PbPb_2body/AnalysisResults_15o.root')
        


cent_counts, cent_edges = an_file["AliAnalysisTaskHyperTriton2He3piML_custom_summary;1" if not isAOD else 'Hypertriton_summary;1'][11].to_numpy()
cent_bin_centers = (cent_edges[:-1]+cent_edges[1:])/2
cent_range_map = np.logical_and(cent_bin_centers > cent_class[0], cent_bin_centers < cent_class[1])
counts_cent_range = cent_counts[cent_range_map]
n_ev = np.sum(counts_cent_range)


if is2015==False and minimumBias:
    an_file = uproot.open("/data/fmazzasc/PbPb_2body/2018/AnalysisResults_18qr.root")
    n_ev = hp.get_number_of_MB_ev(cent_class, an_file)


if is_2018_and_2015:
    cent_counts, cent_edges = an_2015["AliAnalysisTaskHyperTriton2He3piML_custom_summary;1"][11].to_numpy()  ##does not wprk with uproot3 
    cent_bin_centers = (cent_edges[:-1]+cent_edges[1:])/2
    cent_range_map = np.logical_and(cent_bin_centers > cent_class[0], cent_bin_centers < cent_class[1])
    counts_cent_range = cent_counts[cent_range_map]
    n_ev += np.sum(counts_cent_range)
    df = pd.concat([df, df_2015])

print(f'Number of events: {n_ev}')

if minimumBias: 
    df = df.query('trigger%2!=0 and trigger!=0')



# pt_bins = [[3, 3.5],[3.5,4],[4,4.5],[4.5,5],[5,6],[6,9]]
pt_bins = [[2,9]]

bw_index = 0
if cent_class[0]==10 and len(pt_bins)==1:
    bw_index = 1

if cent_class[0]==30 and len(pt_bins)==1:
    bw_index = 2

if cent_class[0]==50 and len(pt_bins)==1:
    bw_index = 3


bw_file = ROOT.TFile('../utils/BlastWaveFits.root')
bw = bw_file.Get(f"BlastWave/BlastWave{bw_index}")

# bw = hp.get_pt_shape_cent(cent_class)

apply_pt_rejection(df_mc, bw)


# df = df.query('Matter==1')

flat_pt_bins = [item for sublist in pt_bins for item in sublist]
bins = np.unique(np.array(flat_pt_bins, dtype=float))


AOD_string = "AOD" if isAOD else "ESD"
MB_string = "_MB" if minimumBias else ""
year_string = "2015" if is2015 else "2018"
year_string = "2015_2018" if is_2018_and_2015==True else year_string
directory = "bins" if len(pt_bins)>1 else "single_bin"

ffile = ROOT.TFile(f"results/{directory}/{AOD_string}_{year_string}_{cent_class[0]}_{cent_class[1]}{MB_string}.root", "recreate")


if isAOD:
    df_rec = df_mc.query('rej>0 and isReconstructed==True and abs(Rapidity)<0.5')
    df_gen = df_mc.query('rej>0 and abs(gRapidity)<0.5')
  
else:
    # df_rec = df_mc.query('rej>0 and pt>0 and abs(Rapidity)<0.5 and Matter==0')
    # df_gen = df_mc.query('rej>0 and abs(gRapidity)<0.5 and gMatter==0')
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
                        # if len(pt_bins)>1 and pt_bin[0]<3:
                        #     cosPA = 0.9998
                        #     pidHe3 = 110
                        print('********************************************************************')
                        cut = f'2<ct<35 and V0CosPA > {cosPA} and NpidClustersHe3 > {pidHe3} and pt > {pt_bin[0]} and pt < {pt_bin[1]} and He3ProngPvDCA > 0.05 and PiProngPvDCA > 0.2 and abs(TPCnSigmaHe3) < 3.5 and ProngsDCA < {pDCA}'
                        cut_cent = cut + f" and {cent_class[0]}<=centrality<{cent_class[1]}"
                        print("#################################")
                        print("CUT: ", cut_cent)
                        print("#################################")
                        print("Efficiency before cuts")
                        print(len(df_rec.query(f'pt > {pt_bin[0]} and pt < {pt_bin[1]}'))/len(df_gen.query(f'gPt > {pt_bin[0]} and gPt < {pt_bin[1]}')))
                        print("#################################")
                        print("Efficiency after cuts")
                        if len(pt_bins)>1:
                            eff = len(df_rec.query(cut))/len(df_gen.query(f'gPt > {pt_bin[0]} and gPt < {pt_bin[1]}'))
                        else:
                            print("SINGLE BIN EFFICIENCY")
                            eff = len(df_rec.query(cut))/len(df_gen)
                        print(eff)

                        df_sel = df.query(cut_cent)
                        selected_data_hist = np.histogram(np.array(df_sel['m']), bins=n_bins, range=mass_range)                        
                        selected_data_hist = hp_std.h1_invmass(selected_data_hist[0], mass_range=mass_range, bins=n_bins, name=f'{cut}')
                        fit_result = hp_std.fit_hist(selected_data_hist, cent_class=cent_class, pt_range=pt_bin, ct_range=[1,35])
                        bin_width = pt_bin[1] - pt_bin[0]
                        if len(pt_bins)>1:
                            pt_spectrum.SetBinContent(ind + 1, fit_result[0]/eff/n_ev/2/bin_width/absorption_corr)
                            pt_spectrum.SetBinError(ind + 1, fit_result[1]/eff/n_ev/2/bin_width/absorption_corr)
                            efficiency.SetBinContent(ind + 1, eff)
                        else:
                            pt_spectrum.SetBinContent(ind + 1, fit_result[0]/eff/n_ev/2/absorption_corr)
                            pt_spectrum.SetBinError(ind + 1, fit_result[1]/eff/n_ev/2/absorption_corr)
                            efficiency.SetBinContent(ind + 1, eff)


pwg = AliPWGFunc()
histo,Integral, integral_error, bw_fit = hp.bw_fit(pt_spectrum, bw, pwg, fit_range=[2,9])
print("yield: ", Integral, ", error: ", integral_error)
histo.SetTitle(";#it{p}_{T} (GeV/#it{c});1/N_{ev} d#it{N}/(dyd#it{p}_{T})")

histo.Write()
bw_fit.SetName(f'func_{cent_class[0]}_{cent_class[1]}')
bw_fit.Write()
pt_spectrum.Write()
efficiency.Write()
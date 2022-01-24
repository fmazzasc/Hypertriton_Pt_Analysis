import sys
sys.path.append('../utils')
import helpers as hp
import uproot
import numpy as np
import ROOT
ROOT.gROOT.SetBatch()

def rename_mc_df_columns(df):
    rename_dict = {}
    rename_dict['ptMC'] = 'gPt'
    rename_dict['etaMC'] = 'gEta'
    rename_dict['ctMC'] = 'gCt'
    rename_dict['yMC'] = 'gRapidity'
    
    df.rename(columns = rename_dict, inplace=True)

def generate_th1(counts, bins, name=''):
    th1 = ROOT.TH1D(f'{name}', f'{name}', len(bins) - 1, np.array(bins, 'double'))
    for index in range(0, len(counts)):
        th1.SetBinContent(index+1, counts[index])
    th1.SetDirectory(0)
    return th1

def presel_eff_hist(df_list, bins):

    counts_rec = np.histogram(df_list[0]["pt"], bins=bins)
    counts_gen = np.histogram(df_list[1]["gPt"], bins=bins)
    eff = counts_rec[0]/counts_gen[0]
    hist_eff = generate_th1(eff, bins, f"fPreselEff_vs_Pt")
    hist_eff.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    hist_eff.GetYaxis().SetTitle('Efficiency')
    hist_eff.SetMinimum(0)
    hist_eff.SetDrawOption("histo")
    hist_eff.SetLineWidth(2)

    # return histogram
    return hist_eff

def apply_pt_rejection(df, pt_shape):
    rej_flag = -1*np.ones(len(df))
    random_arr = np.random.rand(len(df))
    max_bw = pt_shape.GetMaximum()

    for ind, (pt, rand) in enumerate(zip(df['gPt'],random_arr)):
        frac = pt_shape.Eval(pt)/max_bw
        if rand < frac:
            rej_flag[ind] = True
    df['rej'] = rej_flag




bins = np.array([1, 9])

df_mc_esd = uproot.open("/data/fmazzasc/PbPb_2body/SignalTable_20g7_flat_pt.root")["SignalTable"].arrays(library="pd")
df_mc_aod = uproot.open("/data/fmazzasc/PbPb_2body/AOD/HyperTritonTree_MC.root")["HyperTree"].arrays(library="pd")


bw_file = ROOT.TFile('../utils/BlastWaveFits.root')
bw = bw_file.Get(f"BlastWave/BlastWave2")

rename_mc_df_columns(df_mc_aod)
apply_pt_rejection(df_mc_esd, bw)
apply_pt_rejection(df_mc_aod, bw)





df_rec_aod = df_mc_aod.query('rej>0 and isReconstructed==True and abs(Rapidity)<0.5')
df_gen_aod = df_mc_aod.query('rej>0 and abs(gRapidity)<0.5')

df_rec_esd = df_mc_esd.query('rej>0 and pt>0 and abs(Rapidity)<0.5')
df_gen_esd = df_mc_esd.query('rej>0 and abs(gRapidity)<0.5')


presel_eff_esd = presel_eff_hist([df_rec_esd,df_gen_esd], bins=bins)
presel_eff_aod = presel_eff_hist([df_rec_aod,df_gen_aod], bins=bins)


cv = ROOT.TCanvas("")
presel_eff_esd.Draw('')
presel_eff_esd.SetMarkerColor(ROOT.kRed)
presel_eff_esd.SetLineColor(ROOT.kRed)
presel_eff_esd.SetStats(0)
presel_eff_esd.SetTitle("")


presel_eff_aod.Draw('same')

leg = ROOT.TLegend()
leg.AddEntry(presel_eff_esd, 'ESD')
leg.AddEntry(presel_eff_aod, 'AOD')
leg.Draw()

cv.SaveAs('eff.png')

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




bins = np.array([1,2,3,4,5,6,7,8,9])

df_mc_2018 = uproot.open("/data/fmazzasc/PbPb_2body/SignalTable_20g7_flat_pt.root")["SignalTable"].arrays(library="pd")
df_mc_2015 =  uproot.open("/data/fmazzasc/PbPb_2body/SignalTable_16h7abc_flat_pt.root")["SignalTable"].arrays(library="pd")


bw_file = ROOT.TFile('../utils/BlastWaveFits.root')
bw = bw_file.Get(f"BlastWave/BlastWave2")

# rename_mc_df_columns(df_mc_2015)
apply_pt_rejection(df_mc_2018, bw)
apply_pt_rejection(df_mc_2015, bw)





df_rec_2015 = df_mc_2015.query('rej>-2 and pt>0 and abs(Rapidity)<0.5')
df_gen_2015 = df_mc_2015.query('rej>-2 and abs(gRapidity)<0.5')

df_rec_2018 = df_mc_2018.query('rej>-2 and pt>0 and abs(Rapidity)<0.5')
df_gen_2018 = df_mc_2018.query('rej>-2 and abs(gRapidity)<0.5')


presel_eff_2018 = presel_eff_hist([df_rec_2018,df_gen_2018], bins=bins)
presel_eff_2015 = presel_eff_hist([df_rec_2015,df_gen_2015], bins=bins)


cv = ROOT.TCanvas("")
presel_eff_2018.Draw('')
presel_eff_2018.SetMarkerColor(ROOT.kRed)
presel_eff_2018.SetLineColor(ROOT.kRed)
presel_eff_2018.SetStats(0)
presel_eff_2018.SetTitle("")


presel_eff_2015.Draw('same')

leg = ROOT.TLegend()
leg.AddEntry(presel_eff_2018, '2018')
leg.AddEntry(presel_eff_2015, '2015')
leg.Draw()

cv.SaveAs('eff.png')

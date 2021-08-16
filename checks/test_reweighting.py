import helpers as hp
import uproot
import numpy as np
import ROOT
import aghast


def generate_th1(counts, bins, name=''):
    th1 = ROOT.TH1D(f'{name}', f'{name}', len(bins) - 1, np.array(bins, 'double'))
    for index in range(0, len(counts)):
        th1.SetBinContent(index+1, counts[index])
    th1.SetDirectory(0)
    return th1

def presel_eff_hist(df_list, col_name, split, cent_bins, bins):

    counts_rec = np.histogram(df_list[0]["pt"], bins=bins)
    counts_gen = np.histogram(df_list[1]["pt"], bins=bins)
    eff = counts_rec[0]/counts_gen[0]
    hist_eff = generate_th1(eff, bins, f"fPreselEff_vs_{col_name}_{split}_{cent_bins[0]}_{cent_bins[1]}")
    hist_eff.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    hist_eff.GetYaxis().SetTitle('Efficiency')
    hist_eff.SetMinimum(0)
    hist_eff.SetDrawOption("histo")
    hist_eff.SetLineWidth(2)

    # return histogram
    return hist_eff


df_rec = uproot.open("/data/fmazzasc/PbPb_2body/MC_tables/SignalTable_flat_pt.root")["SignalTable"].arrays(library="pd").query("abs(Rapidity)<0.5 and pt>2")
df_gen = uproot.open("/data/fmazzasc/PbPb_2body/MC_tables/SignalTable_flat_pt.root")["GenTable"].arrays(library="pd").query("abs(rapidity)<0.5")

cent_bins = [0,10]
pt_bins = [2,3,4,5,6,7,8,9]
pt_fine_bins = np.linspace(2,9,1000)

hist_eff_pt = presel_eff_hist([df_rec, df_gen], 'pt', "matter", cent_bins, pt_bins)
for ibin in range(1,hist_eff_pt.GetNbinsX() + 1):
    print(hist_eff_pt.GetBinContent(ibin))

fromroot = aghast.from_root(hist_eff_pt)
np_hist = aghast.to_numpy(fromroot)

hist_eff_pt_fine = presel_eff_hist([df_rec, df_gen], 'pt', "fine_eff_reweighting", cent_bins, pt_fine_bins)

bw_file = ROOT.TFile("BlastWaveFits.root")
hypPtShape = bw_file.Get("BlastWave/BlastWave0")

hp.reweight_efficiency(hist_eff_pt, hist_eff_pt_fine, hypPtShape)
print(hist_eff_pt.GetBinContent(ibin))


df_rec = uproot.open("/data/fmazzasc/PbPb_2body/MC_tables/SignalTable_bw_0_10.root")["SignalTable"].arrays(library="pd").query("abs(Rapidity)<0.5 and pt>2")
df_gen = uproot.open("/data/fmazzasc/PbPb_2body/MC_tables/SignalTable_bw_0_10.root")["GenTable"].arrays(library="pd").query("abs(rapidity)<0.5")
hist_eff_pt_distr = presel_eff_hist([df_rec, df_gen], 'pt', "distr_reweighting", cent_bins, pt_bins)

cv = ROOT.TCanvas('hist_comp')
hist_eff_pt.Draw()
hist_eff_pt_distr.Draw('same')



ffile = ROOT.TFile("prova.root", "recreate")
hist_eff_pt.Write()
hist_eff_pt_fine.Write()
hypPtShape.Write()
cv.Write()
ffile.Close()

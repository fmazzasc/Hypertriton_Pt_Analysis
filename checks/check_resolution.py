import uproot
import ROOT


def fill_2D_hist(hist, arr1, arr2):
    for x, y in zip(arr1, arr2):
        hist.Fill(x, y)


df = uproot.open('/data/fmazzasc/PbPb_2body/SignalTable_20g7_flat_pt.root')[
    'SignalTable'].arrays(library='pd').query('pt>0')

hist_pt = ROOT.TH2D("res pt", "; pT gen; (pT rec - pT gen)/pT_gen",
                    100, 1, 10, 100, -0.2, 0.2)

hist_ct = ROOT.TH2D("res ct", "; ct gen; (ct rec - ct gen)/ct_gen", 100, 1, 35, 100, -0.15, 0.15)

hist_rapidity = ROOT.TH2D(
    "res rap", "; rap gen; (rap rec - rap gen)", 100, -1, 1, 100, -0.01, 0.01)


res_pt = (df['gPt'] - df['pt'])/df['gPt']
res_rap = (df['gRapidity'] - df['Rapidity'])
res_ct = (df['gCt'] - df['ct'])/df['gCt']


fill_2D_hist(hist_pt, df['gPt'], res_pt)
fill_2D_hist(hist_rapidity, df['gRapidity'], res_rap)
fill_2D_hist(hist_ct, df['gCt'], res_ct)


outfile = ROOT.TFile('resolutions.root', "RECREATE")
hist_pt.Write()
hist_ct.Write()
hist_rapidity.Write()
outfile.Close()

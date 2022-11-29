import uproot
import ROOT


def fill_2D_hist(hist, arr1, arr2):
    for x, y in zip(arr1, arr2):
        hist.Fill(x, y)

def rename_mc_df_columns(df):
    rename_dict = {}
    rename_dict['ptMC'] = 'gPt'
    rename_dict['etaMC'] = 'gEta'
    rename_dict['ctMC'] = 'gCt'
    rename_dict['yMC'] = 'gRapidity'
    df.rename(columns = rename_dict, inplace=True)


df = uproot.open('/data/fmazzasc/PbPb_2body/AOD/HyperTritonTree_MC.root')[
    'HyperTree'].arrays(library='pd').query('isReconstructed==1')

rename_mc_df_columns(df)

hist_pt = ROOT.TH2D("res pt", "; pT gen; (pT rec - pT gen)/pT_gen",
                    100, 1, 10, 100, -0.2, 0.2)

hist_ct = ROOT.TH2D("res ct", "; ct gen; (ct gen - ct rec)/ct_gen", 100, 1, 35, 300, -0.6, 0.6)

hist_rapidity = ROOT.TH2D(
    "res rap", "; rap gen; (rap rec - rap gen)", 100, -1, 1, 100, -0.01, 0.01)

hist_cpa = ROOT.TH2D("res cpa", "; pT gen; CosPA rec", 100, 1, 10, 100, 0.998, 1)
hist_cpa_ct = ROOT.TH2D("res cpa ct", "; ct gen; CosPA rec", 100, 0, 20, 100, -1, 1)


print(df.columns)

res_pt = (df['gPt'] - df['pt'])/df['gPt']
res_rap = (df['gRapidity'] - df['Rapidity'])
res_ct = (df['gCt'] - df['ct'])/df['gCt']


res_cpa = (df['V0CosPA'])




fill_2D_hist(hist_pt, df['gPt'], res_pt)
fill_2D_hist(hist_rapidity, df['gRapidity'], res_rap)
fill_2D_hist(hist_ct, df['gCt'], res_ct)
fill_2D_hist(hist_cpa, df['gPt'], res_cpa)
fill_2D_hist(hist_cpa_ct, df['gCt'], res_cpa)




outfile = ROOT.TFile('resolutions.root', "RECREATE")
hist_pt.Write()
hist_ct.Write()
hist_rapidity.Write()
hist_cpa.Write()
hist_cpa_ct.Write()


outfile.Close()

import sys
import os
import yaml
import argparse
import numpy as np
sys.path.append('utils')
import helpers as hp

import ROOT
ROOT.gROOT.SetBatch()


ROOT.gROOT.ProcessLine(".L utils/alipwgfunc/AliPWGFunc.cxx++")
from ROOT import AliPWGFunc


yields = ROOT.TH1D('yields', 'yields', 3, 0, 2)
stats = ROOT.TH1D('stats', 'stats', 3, 0, 2 )
systs = ROOT.TH1D('systs', 'systs', 3, 0, 2 )

pwg = AliPWGFunc()
bw_file = ROOT.TFile('utils/He3_fits.root')
corrected_yields_file = ROOT.TFile.Open('results/bins_offline_2018_TOF' + '/corrected_yields.root')




hyp_lambda_ratio = {}
histos = []
integrals = []

#### 0-10
histo_matter = corrected_yields_file.Get('matter_0_10/fYields_matter_0_10')
histo_matter.SetDirectory(0)
histo_antimatter = corrected_yields_file.Get('antimatter_0_10/fYields_antimatter_0_10')
histo_antimatter.SetDirectory(0)
bw = bw_file.Get('BlastWave_5_10')
histo, integral, integral_error, bw_fit = hp.he3_bw_fit(histo_matter, bw, pwg, fit_range=[3,9])
histo_anti, integral_anti, integral_error_anti, bw_fit = hp.he3_bw_fit(histo_antimatter, bw, pwg, fit_range=[3,9])
print(integral, integral_error, integral_anti, integral_error_anti)

hyp_lambda_ratio['010'] = [integral/19.575, integral_error/19.575]
histos.append(histo)
yields.SetBinContent(1, integral)
stats.SetBinContent(1, integral_error)
systs.SetBinContent(1, integral_error)




#### 10-30
corrected_yields_file = ROOT.TFile.Open('results/bins_offline_merged' + '/corrected_yields.root')
histo_matter = corrected_yields_file.Get('matter_10_30/fYields_matter_10_30')
histo_matter.SetDirectory(0)
histo_antimatter = corrected_yields_file.Get('antimatter_10_30/fYields_antimatter_10_30')
histo_antimatter.SetDirectory(0)
bw = bw_file.Get('BlastWave_10_30')
histo, integral, integral_error, bw_fit = hp.he3_bw_fit(histo_matter, bw, pwg, fit_range=[2,9])
histo_anti, integral_anti, integral_error_anti, bw_fit = hp.he3_bw_fit(histo_antimatter, bw, pwg, fit_range=[2,9])
hyp_lambda_ratio['1030'] = [integral/12, integral_error/12]

histos.append(histo)
yields.SetBinContent(2, integral)
stats.SetBinContent(2, integral_error)
systs.SetBinContent(2, integral_error)



#### 30-50
corrected_yields_file = ROOT.TFile.Open('results/bins_offline_2018' + '/corrected_yields.root')
histo_matter = corrected_yields_file.Get('matter_30_50/fYields_matter_30_50')
histo_matter.SetDirectory(0)
histo_antimatter = corrected_yields_file.Get('antimatter_30_50/fYields_antimatter_30_50')
histo_antimatter.SetDirectory(0)
bw = bw_file.Get('BlastWave_30_50')
histo, integral, integral_error, bw_fit = hp.he3_bw_fit(histo_matter, bw, pwg, fit_range=[2,9])
histo_anti, integral_anti, integral_error_anti, bw_fit = hp.he3_bw_fit(histo_antimatter, bw, pwg, fit_range=[2,9])
hyp_lambda_ratio['3050'] = [integral/5.3, integral_error/5.3]

histos.append(histo)
yields.SetBinContent(3, integral)
stats.SetBinContent(3, integral_error)
systs.SetBinContent(3, integral_error)




outfile = ROOT.TFile('res.root', 'recreate')
for histo in histos:
    print(histo)
    histo.Write()

yields.Write()
stats.Write()
systs.Write()

outfile.Close()


# hyp/lambda curve



# d$N$/d$\eta$ obtained by simple weighted average of the values published in https://arxiv.org/pdf/1910.14401.pdf


kBlueC  = ROOT.TColor.GetColor("#2077b4");
kRedC  = ROOT.TColor.GetColor("#d62827");
kGreenC  = ROOT.TColor.GetColor("#2ba02b");
kOrangeC  = ROOT.TColor.GetColor("#ff7f0f");
kVioletC  = ROOT.TColor.GetColor("#9467bd");
kPinkC  = ROOT.TColor.GetColor("#e377c1");
kGreyC  = ROOT.TColor.GetColor("#7f7f7f");
kBrownC  = ROOT.TColor.GetColor("#8c564c");
kAzureC  = ROOT.TColor.GetColor("#18becf");
kGreenBC  = ROOT.TColor.GetColor("#bcbd21");

hp_ratio_csm_3 = ROOT.TGraphErrors("utils/ProdModels/csm_models/VanillaCSM.S3.Vc.eq.3dVdy.dat","%*s %*s %*s %lg %*s %*s %*s %*s %lg %*s")
hp_ratio_csm_1 = ROOT.TGraphErrors("utils/ProdModels/csm_models/VanillaCSM.S3.Vc.eq.dVdy.dat","%*s %*s %*s %lg %*s %*s %*s %*s %lg %*s")


hp_ratio_csm_1.SetLineColor(922)
hp_ratio_csm_1.SetLineWidth(2)
hp_ratio_csm_1.SetTitle("SHM, #it{Vc} = d#it{V}/d#it{y}")
hp_ratio_csm_1.SetMarkerSize(0)


hp_ratio_csm_3.SetLineColor(922)
hp_ratio_csm_3.SetLineWidth(2)
hp_ratio_csm_3.SetLineStyle(2)
hp_ratio_csm_3.SetMarkerSize(0)
hp_ratio_csm_3.SetTitle("SHM, #it{Vc} = 3d#it{V}/d#it{y}")

n = hp_ratio_csm_1.GetN()
grshade = ROOT.TGraph(2*n)
for i in range(n) : 
   grshade.SetPoint(i, hp_ratio_csm_3.GetPointX(i), hp_ratio_csm_3.GetPointY(i))
   grshade.SetPoint(n + i, hp_ratio_csm_1.GetPointX(n - i -1), hp_ratio_csm_1.GetPointY(n - i - 1))
   
grshade.SetFillColorAlpha(16, 0.571)
# grshade.SetFillStyle(3013)




hp_2body = ROOT.TGraphErrors("utils/ProdModels/coalescence/hp_ratio_2body_coal.csv","%lg %lg %lg")
hp_2body.SetLineColor(kBlueC)
hp_2body.SetMarkerColor(kBlueC)
hp_2body.SetTitle("2-body coalescence")
hp_2body.SetFillStyle(3145)
hp_2body.SetMarkerSize(0)
hp_2body.SetLineWidth(2)

hp_3body = ROOT.TGraphErrors("utils/ProdModels/coalescence/hp_ratio_3body_coal.csv","%lg %lg %lg")
hp_3body.SetLineColor(kAzureC)
hp_3body.SetMarkerColor(kAzureC)
hp_3body.SetTitle("3-body coalescence")
hp_3body.SetFillStyle(3014)
hp_3body.SetMarkerSize(0)
hp_3body.SetLineWidth(2)


hp_3body.SetFillColorAlpha(kAzureC, 0.571)
hp_2body.SetFillColorAlpha(kBlueC, 0.571)

mg = ROOT.TMultiGraph()
mg.Add(hp_2body)
mg.Add(hp_3body)



cv = ROOT.TCanvas("cv", "cv", 700,700)
cv.SetBottomMargin(0.145)
cv.SetLeftMargin(0.17)
cv.SetTopMargin(0.01)
cv.SetRightMargin(0.01)
frame=cv.DrawFrame(5., 1e-7, 200, 6e-6,";#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}; {}_{#Lambda}^{3}H/p")

cv.SetLogx()
cv.SetLogy()
mg.Draw("4al same")
# grshade.Draw("f same")
hp_ratio_csm_1.Draw("L same")
hp_ratio_csm_3.Draw("L same")
# hp_ratio_csm_van.Draw("L same")
mg.GetYaxis().SetRangeUser(2e-7, 2e-5)
mg.GetXaxis().SetRangeUser(5, 3e3)
mg.GetXaxis().SetTitle('#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}')
mg.GetYaxis().SetTitle('{}_{#Lambda}^{3}H/#Lambda')
mg.GetXaxis().SetTitleOffset(1.1)
mg.GetYaxis().SetTitleOffset(1.16)
mg.GetYaxis().SetTitleSize(0.06)
mg.GetXaxis().SetTitleSize(0.06)
mg.GetYaxis().SetLabelSize(0.04)
mg.GetXaxis().SetLabelSize(0.04)



zero = np.array([0], dtype=np.float64)
x=np.array([29.4], dtype=np.float64)
ex=np.array([0.6], dtype=np.float64)
y = np.array([1.19e-6], dtype=np.float64)
ey = np.array([3.49e-7], dtype=np.float64)
eys = np.array([2.2e-7], dtype=np.float64)
ppb_stat040 = ROOT.TGraphErrors(1, x, y,zero, ey)
ppb_stat040.SetLineColor(kRedC)
ppb_stat040.SetMarkerColor(kRedC)
ppb_stat040.SetMarkerStyle(20)
ppb_stat040.SetMarkerSize(1)
ppb_stat040.SetLineWidth(1)
ppb_syst040 = ROOT.TGraphErrors(1, x, y, ex, eys)
ppb_syst040.SetTitle("ALICE p#font[122]{-}Pb, 0#font[122]{-}40%, #sqrt{#it{s}_{NN}} = 5.02 TeV")
ppb_syst040.SetLineColor(kRedC)
ppb_syst040.SetMarkerColor(kRedC)
ppb_syst040.SetFillStyle(0)
ppb_syst040.SetMarkerStyle(20)
ppb_syst040.SetLineWidth(1)
ppb_syst040.SetMarkerSize(1)



x = np.array([30.81], dtype=np.float64)
ex = np.array([0.44], dtype=np.float64)
y = np.array([5.49*10**(-7)], dtype=np.float64)
ey = np.array([1.34*10**(-7)], dtype=np.float64)
eys = np.array([0.72e-7], dtype=np.float64)
zero = np.array([0], dtype=np.float64)
pp_stat = ROOT.TGraphErrors(1,x,y,zero,ey)
pp_stat.SetLineColor(kOrangeC)
pp_stat.SetMarkerColor(kOrangeC)
pp_stat.SetMarkerStyle(21)
pp_stat.SetMarkerSize(1)
pp_stat.SetLineWidth(1)
pp_syst = ROOT.TGraphErrors(1,x,y,ex,eys)
pp_syst.SetTitle("ALICE Preliminary pp, HM trigger, #sqrt{#it{s}} = 13 TeV")
pp_syst.SetLineColor(kOrangeC)
pp_syst.SetMarkerColor(kOrangeC)
pp_syst.SetFillStyle(0)
pp_syst.SetMarkerStyle(21)
pp_syst.SetMarkerSize(1)
pp_syst.SetLineWidth(1)



x = np.array([1447], dtype=np.float64)
ex = np.array([39], dtype=np.float64)
y = np.array([4*hyp_lambda_ratio['010'][0]], dtype=np.float64)
ey = np.array([4*hyp_lambda_ratio['010'][1]], dtype=np.float64)
eys = np.array([4*hyp_lambda_ratio['010'][1]], dtype=np.float64)
zero = np.array([0], dtype=np.float64)
pbpb_stat = ROOT.TGraphErrors(1,x,y,zero,ey)
pbpb_stat.SetLineColor(ROOT.kBlue + 3)
pbpb_stat.SetMarkerColor(ROOT.kBlue +3)
pbpb_stat.SetMarkerStyle(ROOT.kFullDiamond)
pbpb_stat.SetMarkerSize(1.5)
pbpb_stat.SetLineWidth(1)
pbpb_syst = ROOT.TGraphErrors(1,x,y,ex,eys)
pbpb_syst.SetTitle("ALICE Pb#font[122]{-}Pb Preliminary")
pbpb_syst.SetLineColor(ROOT.kBlue +3)
pbpb_syst.SetMarkerColor(ROOT.kBlue +3)
pbpb_syst.SetFillStyle(0)
pbpb_syst.SetMarkerStyle(ROOT.kFullDiamond)
pbpb_syst.SetMarkerSize(1.5)
pbpb_syst.SetLineWidth(1)


x = np.array([980], dtype=np.float64)
ex = np.array([30], dtype=np.float64)
y = np.array([4*hyp_lambda_ratio['1030'][0]], dtype=np.float64)
ey = np.array([4*hyp_lambda_ratio['1030'][1]], dtype=np.float64)
eys = np.array([4*hyp_lambda_ratio['1030'][1]], dtype=np.float64)
zero = np.array([0], dtype=np.float64)
pbpb_1030_stat = ROOT.TGraphErrors(1,x,y,zero,ey)
pbpb_1030_stat.SetLineColor(ROOT.kBlue + 3)
pbpb_1030_stat.SetMarkerColor(ROOT.kBlue +3)
pbpb_1030_stat.SetMarkerStyle(ROOT.kFullDiamond)
pbpb_1030_stat.SetMarkerSize(1.5)
pbpb_1030_stat.SetLineWidth(1)
pbpb_1030_syst = ROOT.TGraphErrors(1,x,y,ex,eys)
pbpb_1030_syst.SetLineColor(ROOT.kBlue +3)
pbpb_1030_syst.SetMarkerColor(ROOT.kBlue +3)
pbpb_1030_syst.SetFillStyle(0)
pbpb_1030_syst.SetMarkerStyle(ROOT.kFullDiamond)
pbpb_1030_syst.SetMarkerSize(1.5)
pbpb_1030_syst.SetLineWidth(1)


x = np.array([400], dtype=np.float64)
ex = np.array([15], dtype=np.float64)
y = np.array([4*hyp_lambda_ratio['3050'][0]], dtype=np.float64)
ey = np.array([4*hyp_lambda_ratio['3050'][1]], dtype=np.float64)
eys = np.array([4*hyp_lambda_ratio['3050'][1]], dtype=np.float64)
zero = np.array([0], dtype=np.float64)
pbpb_3050_stat = ROOT.TGraphErrors(1,x,y,zero,ey)
pbpb_3050_stat.SetLineColor(ROOT.kBlue + 3)
pbpb_3050_stat.SetMarkerColor(ROOT.kBlue +3)
pbpb_3050_stat.SetMarkerStyle(ROOT.kFullDiamond)
pbpb_3050_stat.SetMarkerSize(1.5)
pbpb_3050_stat.SetLineWidth(1)
pbpb_3050_syst = ROOT.TGraphErrors(1,x,y,ex,eys)
pbpb_3050_syst.SetLineColor(ROOT.kBlue +3)
pbpb_3050_syst.SetMarkerColor(ROOT.kBlue +3)
pbpb_3050_syst.SetFillStyle(0)
pbpb_3050_syst.SetMarkerStyle(ROOT.kFullDiamond)
pbpb_3050_syst.SetMarkerSize(1.5)
pbpb_3050_syst.SetLineWidth(1)



leg = ROOT.TLegend(0.221,0.81,0.87,0.96)
leg.SetMargin(0.12)
leg.SetNColumns(1)

pp_stat.Draw("Pz")
pp_syst.Draw("P2")


ppb_stat040.Draw("Pz")
ppb_syst040.Draw("P2")

pbpb_stat.Draw("Pz")
pbpb_syst.Draw("P2")

pbpb_1030_stat.Draw("Pz")
pbpb_1030_syst.Draw("P2")

pbpb_3050_stat.Draw("Pz")
pbpb_3050_syst.Draw("P2")

leg.AddEntry(ppb_syst040,"","pf")
leg.AddEntry(pp_syst,"","pf")
leg.AddEntry(pbpb_syst,"","pf")




leg.SetEntrySeparation(0.2)
legT = ROOT.TLegend(0.62,0.17,0.97,0.44)
legT.SetMargin(0.14)
legT.SetBorderSize(0)
legT.AddEntry(hp_3body, hp_3body.GetTitle(), "LF")
legT.AddEntry(hp_2body,  hp_2body.GetTitle(), "LF")
legT.AddEntry(hp_ratio_csm_1,  hp_ratio_csm_1.GetTitle(), "LF")
legT.AddEntry(hp_ratio_csm_3, hp_ratio_csm_3.GetTitle(), "LF")
# legT.AddEntry(hp_ratio_csm_van)
leg.SetFillStyle(0)
legT.SetFillStyle(0)
leg.Draw()
legT.Draw()

pinfo = ROOT.TPaveText(0.2208,0.75, 0.47, 0.8, 'NDC')
pinfo.SetBorderSize(0)
pinfo.SetFillStyle(0)
# pinfo.SetTextAlign(30+3)
pinfo.SetTextFont(42)
pinfo.AddText('B.R. = 0.25 #pm 0.02')
pinfo.Draw()


cv.Draw()

cv.SaveAs("hl_ratio.png")
import ROOT 
ROOT.gROOT.SetBatch()

aod_file = ROOT.TFile('std_analysis_2018_0_10_AOD.root')
aod_spectrum = aod_file.Get('Yield')

esd_file = ROOT.TFile('std_analysis_2018_0_10.root')
esd_spectrum = esd_file.Get('Yield')

esd_spectrum.GetFunction('fBGBW').SetBit(ROOT.TF1.kNotDraw)
aod_spectrum.GetFunction('fBGBW').SetBit(ROOT.TF1.kNotDraw)


cv = ROOT.TCanvas()
esd_spectrum.Draw('')
esd_spectrum.SetMarkerColor(ROOT.kRed)
esd_spectrum.SetLineColor(ROOT.kRed)
esd_spectrum.SetStats(0)

aod_spectrum.Draw('same')

leg = ROOT.TLegend()
leg.AddEntry(esd_spectrum, 'ESD')
leg.AddEntry(aod_spectrum, 'AOD')
leg.Draw()

cv.SaveAs('comp.png')

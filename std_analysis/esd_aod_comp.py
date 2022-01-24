import ROOT 
ROOT.gROOT.SetBatch()

aod_file = ROOT.TFile('results/bins/AOD_2018_0_10.root')
aod_spectrum = aod_file.Get('Yield')

esd_file = ROOT.TFile('results/bins/ESD_2018_0_10.root')
esd_spectrum = esd_file.Get('Yield')

esd_spectrum.GetFunction('fBGBW').SetBit(ROOT.TF1.kNotDraw)
aod_spectrum.GetFunction('fBGBW').SetBit(ROOT.TF1.kNotDraw)


cv = ROOT.TCanvas()
esd_spectrum.Draw('')
esd_spectrum.SetMarkerColor(ROOT.kRed)
esd_spectrum.SetLineColor(ROOT.kRed)
esd_spectrum.SetStats(0)
esd_spectrum.SetTitle(";#it{p}_{T} (GeV/#it{c});1/N_{ev} d#it{N}/(dyd#it{p}_{T})")
aod_spectrum.Draw('same')

leg = ROOT.TLegend()
leg.AddEntry(esd_spectrum, 'ESD')
leg.AddEntry(aod_spectrum, 'AOD')
leg.Draw()

cv.SaveAs('comp.png')

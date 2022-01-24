import ROOT 
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine(".L ../utils/alipwgfunc/AliPWGFunc.cxx++")
from ROOT import AliPWGFunc
import numpy as np


cent_dict_he = {'0_5':1, '5_10':2, '10_30':3, '30_50':4}
cent_class = [5,10]

he_bw = cent_dict_he[f'{cent_class[0]}_{cent_class[1]}'] if cent_class != [0,10] else 2

def get_bin_integral(histo, ibin):
    bin_content = histo.GetBinContent(ibin)
    bin_error =  histo.GetBinError(ibin)
    bin_width = histo.GetBinWidth(ibin)
    return bin_content*bin_width, bin_error*bin_width

def get_histo_average(hist1, hist2):
    histo = hist1.Clone()
    for ibin in range(1, histo.GetNbinsX() + 1):
        histo.SetBinContent(ibin, 0.5*(hist1.GetBinContent(ibin) + hist2.GetBinContent(ibin)))
        histo.SetBinError(ibin, 0.5*(hist1.GetBinError(ibin) + hist2.GetBinError(ibin)))
    return histo



he3_file = ROOT.TFile('../utils/Anti_fits.root')

if cent_class == [0,10]:
    he3_spectrum1 = he3_file.Get('data/1/stat1')
    he3_spectrum2 = he3_file.Get('data/2/stat2')
    he3_spectrum = get_histo_average(he3_spectrum1, he3_spectrum2)

else:
    he3_spectrum = he3_file.Get(f'data/{he_bw}/stat{he_bw}')

hyp_file = ROOT.TFile(f'results/bins/ESD_2018_{cent_class[0]}_{cent_class[1]}.root')
hyp_spectrum = hyp_file.Get('Yield')
hyp_spectrum.GetFunction('fBGBW').SetBit(ROOT.TF1.kNotDraw)
he3_func = he3_file.Get(f'BGBW/{he_bw}/BGBW{he_bw}')



he3_spectrum_norm = he3_spectrum.Clone('he3_spectrum_norm')
he3_spectrum_norm.Scale(hyp_spectrum.Integral(3,6)/he3_spectrum_norm.Integral(3,6))

out_file = ROOT.TFile(f'results/comp_hyp_he3_{cent_class[0]}_{cent_class[1]}.root', 'RECREATE')

cv = ROOT.TCanvas()
hyp_spectrum.Draw('')
hyp_spectrum.SetStats(0)

he3_spectrum_norm.SetMarkerColor(ROOT.kRed)
he3_spectrum_norm.SetLineColor(ROOT.kRed)
he3_spectrum_norm.Draw('same')
leg = ROOT.TLegend()
leg.AddEntry(hyp_spectrum, 'Hypertriton')
leg.AddEntry(he3_spectrum, 'He3')
leg.Draw()
cv.Write('he3_hyp')



ratio_spectra = hyp_spectrum.Clone('Ratio')
hyp_integral_arr = []
he3_integral_arr = []

for ibin_hyp in range(1, ratio_spectra.GetNbinsX() + 1):
    hyp_integral, hyp_error = get_bin_integral(hyp_spectrum, ibin_hyp) ##branching ratio
    hyp_integral *= 4
    hyp_error *= 4

    lower_edge = hyp_spectrum.GetBinLowEdge(ibin_hyp)
    upper_edge = hyp_spectrum.GetBinLowEdge(ibin_hyp) + hyp_spectrum.GetBinWidth(ibin_hyp)


    if ibin_hyp==ratio_spectra.GetNbinsX():
        ratio_spectra.SetBinContent(ibin_hyp, 0)
        continue

    he3_integral = 0
    he3_error = 0

    for ibin_he3 in range(1, he3_spectrum.GetNbinsX() + 1):
        lower_edge_he3 = he3_spectrum.GetBinLowEdge(ibin_he3)
        if lower_edge_he3==lower_edge:
            for ibin_up in range(ibin_he3, he3_spectrum.GetNbinsX() + 1):
                upper_edge_he3 = he3_spectrum.GetBinLowEdge(ibin_up) + he3_spectrum.GetBinWidth(ibin_up)
                if upper_edge_he3<=upper_edge:
                    he3_integral_incr, he3_error_incr = get_bin_integral(he3_spectrum, ibin_up)
                    he3_integral += he3_integral_incr
                    he3_error += he3_error_incr


                    print(lower_edge, upper_edge, lower_edge_he3, upper_edge_he3, he3_integral, he3_error)


    print(lower_edge, upper_edge, hyp_integral, hyp_error) 

    ratio = hyp_integral/he3_integral
    ratio_error = ratio*np.sqrt((hyp_error/hyp_integral)**2 + (he3_error/he3_integral)**2) 
    ratio_spectra.SetBinContent(ibin_hyp, ratio)
    ratio_spectra.SetBinError(ibin_hyp, ratio_error)
    ratio_spectra.SetTitle(";#it{p}_{T} (GeV/#it{c});{}_{#Lambda}^{3}H/^{3}He;")
    hyp_integral_arr.append(hyp_integral)
    he3_integral_arr.append(he3_integral)



thermal_pred = ROOT.TLine(2,0.335, 6,0.335)
thermal_pred.SetLineStyle(ROOT.kDashed)
thermal_pred.SetLineColor(ROOT.kRed)

# hyp_spectrum.Fit(he3_func)
hyp_spectrum.GetXaxis().SetRangeUser(2,6)
ratio_spectra.GetXaxis().SetRangeUser(2,6)


cv = ROOT.TCanvas('ratio')
ratio_spectra.Draw()
thermal_pred.Draw('same')
leg = ROOT.TLegend(0.6,0.3, 0.8,0.35)
leg.AddEntry(thermal_pred, 'SHM, V_{c} = dV/dy', "L")
leg.Draw()
cv.Write()

print(hyp_integral_arr, he3_integral_arr)


pwg = AliPWGFunc()
hyp_integral_arr = np.sum(hyp_integral_arr[1: -1])
params = he3_func.GetParameters()

# scaling = he3_func.Integral(2,6)/hyp_spectrum.Integral(2,6)
scaling = he3_func.Integral(3,6)/(hyp_integral_arr/4)
bw = pwg.GetBGBW(params[0], params[1], params[2], params[3], params[4]/16)   



cv = ROOT.TCanvas('he3_bw')
hyp_spectrum.Draw()
bw.Draw('same')
leg = ROOT.TLegend()
leg.AddEntry(bw, 'He3 Blast Wave')
leg.Draw()

cv.Write()



out_file.Close()

print('yield ratio: ', 4/(scaling))
import ROOT 
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine(".L ../utils/alipwgfunc/AliPWGFunc.cxx++")
from ROOT import AliPWGFunc
import numpy as np


cent_dict_he = {'0_5':1, '5_10':2, '10_30':3, '30_50':4}
cent_class = [0,10]

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

hyp_file = ROOT.TFile(f'../results/bins_offline_2018_AOD/corrected_yields.root')
hyp_spectrum = hyp_file.Get('all_0_10/fYields_all_0_10')
# hyp_spectrum.GetFunction('fBGBW').SetBit(ROOT.TF1.kNotDraw)
he3_func = he3_file.Get(f'BGBW/{he_bw}/BGBW{he_bw}')


out_file = ROOT.TFile(f'comp_hyp_he3_{cent_class[0]}_{cent_class[1]}.root', 'RECREATE')



ratio_spectra = hyp_spectrum.Clone('Ratio')
hyp_integral_arr = []
he3_integral_arr = []

for ibin_hyp in range(1, ratio_spectra.GetNbinsX() + 1):
    hyp_integral, hyp_error = get_bin_integral(hyp_spectrum, ibin_hyp) ##branching ratio
    hyp_integral *= 2
    hyp_error *= 2

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



thermal_pred = ROOT.TLine(3,0.335, 6,0.335)
thermal_pred.SetLineStyle(ROOT.kDashed)
thermal_pred.SetLineColor(ROOT.kRed)

# hyp_spectrum.Fit(he3_func)
hyp_spectrum.GetXaxis().SetRangeUser(3,6)


ratio_spectra.GetXaxis().SetRangeUser(3,6)

ratio_spectra.Fit('pol0')

cv = ROOT.TCanvas('ratio')
ratio_spectra.Draw()
thermal_pred.Draw('same')
leg = ROOT.TLegend(0.6,0.3, 0.8,0.35)
leg.AddEntry(thermal_pred, 'SHM, V_{c} = dV/dy', "L")
leg.Draw()

ratio_spectra.SetMinimum(0)
ratio_spectra.SetMaximum(0.4)

cv.Write()

out_file.Close()

import ROOT 
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine(".L ../utils/alipwgfunc/AliPWGFunc.cxx++")
from ROOT import AliPWGFunc
import numpy as np

cent_dict_he = {'0_5':1, '5_10':2, '10_30':3, '30_50':4}
cent_class = [5,10]
he_bw = cent_dict_he[f'{cent_class[0]}_{cent_class[1]}']

def get_bin_integral(histo, ibin):
    bin_content = histo.GetBinContent(ibin)
    bin_width = histo.GetBinWidth(ibin)
    return bin_content*bin_width



he3_file = ROOT.TFile('../utils/Anti_fits.root')
he3_spectrum = he3_file.Get(f'data/{he_bw}/stat{he_bw}')

hyp_file = ROOT.TFile(f'results/ESD_2018_{cent_class[0]}_{cent_class[1]}.root')
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

for ibin_hyp in range(1, ratio_spectra.GetNbinsX() + 1):
    hyp_integral = get_bin_integral(hyp_spectrum, ibin_hyp)*4 ##branching ratio
    lower_edge = hyp_spectrum.GetBinLowEdge(ibin_hyp)
    upper_edge = hyp_spectrum.GetBinLowEdge(ibin_hyp) + hyp_spectrum.GetBinWidth(ibin_hyp)


    if ibin_hyp==ratio_spectra.GetNbinsX():
        ratio_spectra.SetBinContent(ibin_hyp, 0)
        continue

    he3_integral = 0

    for ibin_he3 in range(1, he3_spectrum.GetNbinsX() + 1):
        lower_edge_he3 = he3_spectrum.GetBinLowEdge(ibin_he3)
        if lower_edge_he3==lower_edge:
            for ibin_up in range(ibin_he3, he3_spectrum.GetNbinsX() + 1):
                upper_edge_he3 = he3_spectrum.GetBinLowEdge(ibin_up) + he3_spectrum.GetBinWidth(ibin_up)
                if upper_edge_he3<=upper_edge:
                    he3_integral += get_bin_integral(he3_spectrum, ibin_up)
                    # print(lower_edge, upper_edge, lower_edge_he3, upper_edge_he3, he3_integral)


    print(lower_edge, upper_edge, hyp_integral)    

    ratio_spectra.SetBinContent(ibin_hyp, hyp_integral/he3_integral)
    hyp_integral_arr.append(hyp_integral)

ratio_spectra.Write()

pwg = AliPWGFunc()

hyp_integral_arr = np.sum(hyp_integral_arr[1: -1])

params = he3_func.GetParameters()

# scaling = he3_func.Integral(2,6)/hyp_spectrum.Integral(2,6)
scaling = he3_func.Integral(3,6)/(hyp_integral_arr/4)
bw = pwg.GetBGBW(params[0], params[1], params[2], params[3], params[4]/15)   

# hyp_spectrum.Fit(he3_func)

cv = ROOT.TCanvas('he3_bw')
hyp_spectrum.Draw()
bw.Draw('same')

leg = ROOT.TLegend()
leg.AddEntry(bw, 'He3 Blast Wave')
leg.Draw()

cv.Write()



out_file.Close()

print('yield ratio: ', 4/(scaling))







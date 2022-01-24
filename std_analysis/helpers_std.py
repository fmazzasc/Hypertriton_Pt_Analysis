
import ROOT
import numpy as np


def fit_hist(histo, cent_class, pt_range, ct_range, nsigma=3, model="pol2", fixsigma=-1, sigma_limits=None, mode=3, split =''):
    # canvas for plotting the invariant mass distribution
    cv = ROOT.TCanvas(f'cv_{histo.GetName()}')

    # define the number of parameters depending on the bkg model
    if 'pol' in str(model):
        n_bkgpars = int(model[3]) + 1
    elif 'expo' in str(model):
        n_bkgpars = 2
    else:
        print(f'Unsupported model {model}')

    # define the fit function bkg_model + gauss
    fit_tpl = ROOT.TF1('fitTpl', f'{model}(0)+gausn({n_bkgpars})', 0, 5)

    # redefine parameter names for the bkg_model
    for i in range(n_bkgpars):
        fit_tpl.SetParName(i, f'B_{i}')

    # define parameter names for the signal fit
    fit_tpl.SetParName(n_bkgpars, 'N_{sig}')
    fit_tpl.SetParName(n_bkgpars + 1, '#mu')
    fit_tpl.SetParName(n_bkgpars + 2, '#sigma')
    # define parameter values and limits
    fit_tpl.SetParameter(n_bkgpars, 40)
    fit_tpl.SetParLimits(n_bkgpars, 0.001, 10000)
    fit_tpl.SetParameter(n_bkgpars + 1, 2.991)
    fit_tpl.SetParLimits(n_bkgpars + 1, 2.986, 3)

    # define signal and bkg_model TF1 separately
    sigTpl = ROOT.TF1('fitTpl', 'gausn(0)', 0, 5)
    bkg_tpl = ROOT.TF1('fitTpl', f'{model}(0)', 0, 5)

    # plotting stuff for fit_tpl
    fit_tpl.SetNpx(300)
    fit_tpl.SetLineWidth(2)
    fit_tpl.SetLineColor(2)
    # plotting stuff for bkg model
    bkg_tpl.SetNpx(300)
    bkg_tpl.SetLineWidth(2)
    bkg_tpl.SetLineStyle(2)
    bkg_tpl.SetLineColor(2)

    # define limits for the sigma if provided
    if sigma_limits != None:
        fit_tpl.SetParameter(n_bkgpars + 2, 0.5 *
                             (sigma_limits[0] + sigma_limits[1]))
        fit_tpl.SetParLimits(n_bkgpars + 2, sigma_limits[0], sigma_limits[1])
    # if the mc sigma is provided set the sigma to that value
    elif fixsigma > 0:
        fit_tpl.FixParameter(n_bkgpars + 2, fixsigma)
    # otherwise set sigma limits reasonably
    else:
        fit_tpl.SetParameter(n_bkgpars + 2, 0.002)
        fit_tpl.SetParLimits(n_bkgpars + 2, 0.001, 0.003)

    ########################################
    # plotting the fits
    ax_titles = ';m (^{3}He + #pi) (GeV/#it{c})^{2};Counts' + f' / {round(1000 * histo.GetBinWidth(1), 2)} MeV'

    # invariant mass distribution histo and fit
    histo.UseCurrentStyle()
    histo.SetLineColor(1)
    histo.SetMarkerStyle(20)
    histo.SetMarkerColor(1)
    histo.SetTitle(ax_titles)
    histo.SetMaximum(1.5 * histo.GetMaximum())
    histo.Fit(fit_tpl, "QRL", "", 2.96, 3.04)
    histo.Fit(fit_tpl, "QRL", "", 2.96, 3.04)
    histo.SetDrawOption("e")
    histo.GetXaxis().SetRangeUser(2.96, 3.04)
    # represent the bkg_model separately
    bkg_tpl.SetParameters(fit_tpl.GetParameters())
    bkg_tpl.SetLineColor(600)
    bkg_tpl.SetLineStyle(2)
    bkg_tpl.Draw("same")
    # represent the signal model separately
    sigTpl.SetParameter(0, fit_tpl.GetParameter(n_bkgpars))
    sigTpl.SetParameter(1, fit_tpl.GetParameter(n_bkgpars+1))
    sigTpl.SetParameter(2, fit_tpl.GetParameter(n_bkgpars+2))
    sigTpl.SetLineColor(600)
    # sigTpl.Draw("same")

    # get the fit parameters
    mu = fit_tpl.GetParameter(n_bkgpars+1)
    muErr = fit_tpl.GetParError(n_bkgpars+1)
    sigma = fit_tpl.GetParameter(n_bkgpars+2)
    sigmaErr = fit_tpl.GetParError(n_bkgpars+2)
    signal = fit_tpl.GetParameter(n_bkgpars) / histo.GetBinWidth(1)
    errsignal = fit_tpl.GetParError(n_bkgpars) / histo.GetBinWidth(1)
    bkg = bkg_tpl.Integral(mu - nsigma * sigma, mu +
                           nsigma * sigma) / histo.GetBinWidth(1)

    if bkg > 0:
        errbkg = np.sqrt(bkg)
    else:
        errbkg = 0
    # compute the significance
    if signal+bkg > 0:
        signif = signal/np.sqrt(signal+bkg)
        deriv_sig = 1/np.sqrt(signal+bkg)-signif/(2*(signal+bkg))
        deriv_bkg = -signal/(2*((signal+bkg)**(1.5)))
        errsignif = np.sqrt((errsignal*deriv_sig)**2+(errbkg*deriv_bkg)**2)
    else:
        signif = 0
        errsignif = 0

    # print fit info on the canvas
    pinfo2 = ROOT.TPaveText(0.5, 0.5, 0.91, 0.9, "NDC")
    pinfo2.SetBorderSize(0)
    pinfo2.SetFillStyle(0)
    pinfo2.SetTextAlign(30+3)
    pinfo2.SetTextFont(42)

    string = f'ALICE Internal, Pb-Pbs {cent_class[0]}-{cent_class[1]}%'
    pinfo2.AddText(string)
    
    decay_label = {
        "": '{}^{3}_{#Lambda}H#rightarrow ^{3}He#pi^{-} + c.c.',
        "_matter": '{}^{3}_{#Lambda}H#rightarrow ^{3}He#pi^{-}',
        "_antimatter": '{}^{3}_{#bar{#Lambda}}#bar{H}#rightarrow ^{3}#bar{He}#pi^{+}',
    }

    string = decay_label[split] + ', %i #leq #it{ct} < %i cm %.1f #leq #it{p}_{T} < %.1f GeV/#it{c} ' % (ct_range[0], ct_range[1], pt_range[0], pt_range[1])
    pinfo2.AddText(string)

    string = f'Significance ({nsigma:.0f}#sigma) {signif:.1f} #pm {errsignif:.1f} '
    pinfo2.AddText(string)

    string = f'S ({nsigma:.0f}#sigma) {signal:.0f} #pm {errsignal:.0f}'
    pinfo2.AddText(string)

    string = f'B ({nsigma:.0f}#sigma) {bkg:.0f} #pm {errbkg:.0f}'
    pinfo2.AddText(string)

    if bkg > 0:
        ratio = signal/bkg
        string = f'S/B ({nsigma:.0f}#sigma) {ratio:.4f}'

    pinfo2.AddText(string)
    pinfo2.Draw()
    ROOT.gStyle.SetOptStat(0)

    st = histo.FindObject('stats')
    if isinstance(st, ROOT.TPaveStats):
        st.SetX1NDC(0.12)
        st.SetY1NDC(0.62)
        st.SetX2NDC(0.40)
        st.SetY2NDC(0.90)
        st.SetOptStat(0)

    histo.Write()
    cv.Write()

    return (signal, errsignal, signif, errsignif, mu, muErr, sigma, sigmaErr)



def h1_invmass(counts, mass_range=[2.96, 3.04] , bins=34, name=''):
    th1 = ROOT.TH1D(f'{name}', f'{name}_x', int(bins), mass_range[0], mass_range[1])
    for index in range(0, len(counts)):
        th1.SetBinContent(index+1, counts[index])
        # th1.SetBinError(index + 1, np.sqrt(counts[index]))
    th1.SetDirectory(0)
    return th1
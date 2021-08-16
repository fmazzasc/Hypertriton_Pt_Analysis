import numpy as np
import ROOT

fitted_file = ROOT.TFile("../results/bins_offline_2018/pt_spectra_iter.root")
fitted_func = fitted_file.Get('results/fYields_matter_0_10;1').GetFunction('fBGBW')

print(fitted_func.Integral(0,2)/fitted_func.Integral(0,10))


bw_file = ROOT.TFile("../utils/BlastWaveFits.root")
bw_shape = bw_file.Get('BlastWave/BlastWave0')

print(bw_shape.Integral(0,2)/bw_shape.Integral(0,10))

fact = 1+(bw_shape.Integral(0,2)/bw_shape.Integral(0,10))

print(fitted_func.Integral(2,10))
print(fitted_func.Integral(2,10)*fact)
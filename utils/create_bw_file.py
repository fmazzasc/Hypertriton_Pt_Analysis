import ROOT 

out_file = ROOT.TFile('He3_fits.root', "recreate")

bw_old = ROOT.TFile('BlastWaveFits.root')
bw_new = ROOT.TFile('Anti_fits.root')


cent_old = [[0,10]]
cent_new = [[0,5], [5,10], [10,30], [30,50], [50,90]]

bw_0 = bw_new.Get('BGBW/1/BGBW1')
bw_0.SetName('BlastWave_0_5')

bw_1 = bw_new.Get('BGBW/2/BGBW2')
bw_1.SetName('BlastWave_5_10')

bw_2 = bw_new.Get('BGBW/3/BGBW3')
bw_2.SetName('BlastWave_10_30')

bw_3 = bw_new.Get('BGBW/4/BGBW4')
bw_3.SetName('BlastWave_30_50')

bw_4 = bw_new.Get('BGBW/5/BGBW5')
bw_4.SetName('BlastWave_50_90')

bw_5 = bw_old.Get('BlastWave/BlastWave0')
bw_5.SetName('BlastWave_0_10')

out_file.cd()
bw_0.Write()
bw_1.Write()
bw_2.Write()
bw_3.Write()
bw_4.Write()
bw_5.Write()

out_file.Close()

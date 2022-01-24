import ROOT
input_file0 = ROOT.TFile('ESD_2018_0_10.root')
bw_0 = input_file0.Get('func_0_10')
# input_file1 = ROOT.TFile('ESD_2018_5_10.root')
# bw_1 = input_file1.Get('func_5_10')


input_file2 = ROOT.TFile('ESD_2015_2018_10_30_MB.root')
bw_2 = input_file2.Get('func_10_30')



output_file = ROOT.TFile('hyp_blast_waves.root', 'recreate')
bw_0.Write()
# bw_1.Write()
bw_2.Write()



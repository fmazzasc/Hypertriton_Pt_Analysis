import os
import sys
import yaml
import argparse
import ROOT



parser = argparse.ArgumentParser(prog='iter_analysis', allow_abbrev=True)
parser.add_argument('config', help='Path to the YAML configuration file')
args = parser.parse_args()
##################################
with open(os.path.expandvars(args.config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


iter = 0
while(iter<4):
    os.system(f"python3 ml_analysis.py {args.config} -split -train -n_iter {iter}")
    os.system(f"python3 ml_analysis.py {args.config} -split  -application -n_iter {iter}")
    os.system(f"python3 signal_extraction.py {args.config}")
    os.system(f"python3 significance_scan.py {args.config}")
    os.system(f"python3 pt_spectra.py {args.config}")
    iter += 1



# ##################################
# bw_file = ROOT.TFile('utils/BlastWaveFits.root')
# bw_file2 = ROOT.TFile("utils/Anti_fits.root")

# def create_tf1_dictionary():
#     tf1_dic = {}
#     tf1_dic['0_10'] = bw_file.Get('BlastWave/BlastWave0')
#     tf1_dic['10_30'] = bw_file2.Get(f"BGBW/3/BGBW3")
#     tf1_dic['30_50'] = bw_file2.Get(f"BGBW/4/BGBW4")
#     tf1_dic['50_90'] = bw_file2.Get(f"BGBW/5/BGBW5")
#     return tf1_dic
# ###################################
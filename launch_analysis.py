import os
import argparse

parser = argparse.ArgumentParser(prog='ml_analysis', allow_abbrev=True)
parser.add_argument('config', help='Path to the YAML configuration file')
args = parser.parse_args()
##################################



os.system(f"python3 ml_analysis.py {args.config} -split -train -n_iter 1")
os.system(f"python3 ml_analysis.py {args.config} -split  -application -n_iter 1")
os.system(f"python3 signal_extraction.py {args.config}")
os.system(f"python3 significance_scan.py {args.config}")
os.system(f"python3 pt_spectra.py {args.config}")






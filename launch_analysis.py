import os
import argparse

parser = argparse.ArgumentParser(prog='ml_analysis', allow_abbrev=True)
parser.add_argument('config', help='Path to the YAML configuration file')
args = parser.parse_args()
##################################



os.system(f"python3 ml_analysis.py {args.config} -split -eff -train -computescore")
os.system(f"python3 ml_analysis.py {args.config} -split  -application")
os.system(f"python3 signal_extraction.py {args.config}")
os.system(f"python3 significance_scan.py {args.config}")
os.system(f"python3 he3_absorption_analysis.py {args.config}")






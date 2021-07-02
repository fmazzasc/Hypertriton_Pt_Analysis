import os


os.system("python3 ml_analysis.py -split -eff -train -computescore")
os.system("python3 ml_analysis.py -split  -application")
os.system("python3 signal_extraction.py")




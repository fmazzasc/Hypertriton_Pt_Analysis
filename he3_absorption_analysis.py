#!/usr/bin/env python3
import os
import re
import pickle
import warnings
import numpy as np
import ROOT
import yaml

HE_3_MASS = 2.809230089

##################################################################
# read configuration file
##################################################################
config = 'config.yaml'
with open(os.path.expandvars(config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

PT_BINS_CENT = params['PT_BINS_CENT']
CENT_BINS_ANALYSIS = params['CENTRALITY_LIST']
RESULTS_SUBDIR = params['RESULTS_SUBDIR']


##################################################################
# define helper function

def pt_bins_generator(key):   #returns the pT bin corresponding to a given centrality. If there is no correspondance, it returns a None
    for cent_ind, cent_bin in enumerate(CENT_BINS_ANALYSIS):
        key_cent = re.findall(r"[-+]?\d*\.\d+|\d+", key) 
        if str(cent_bin[0]) == key_cent[0] and str(cent_bin[1]) == key_cent[1] :
            pt_bins = PT_BINS_CENT[cent_ind]
            return pt_bins
    pt_bins = None
    return pt_bins
################################################################s

cent_bins_new =  [[0,5],[5,10],[10,30],[30,50],[50,90]]
cent_bins_MB = [[0, 10], [10, 40], [40, 90]]
##################################################################

if not os.path.isdir(f'plots/absorption_correction'):
    os.mkdir(f'plots/absorption_correction')

split_list = ['antimatter', 'matter']

# mc input file
mc_file = '/data/fmazzasc/PbPb_2body/AnalysisResults_He3_abs.root'
outfile = ROOT.TFile(f"results{RESULTS_SUBDIR}/He3_abs.root", "recreate")
centfile = ROOT.TFile("/data/fmazzasc/PbPb_2body/AnalysisResults_18.root")

# get event centrality distribution
cent_dist = centfile.Get("AliAnalysisTaskHyperTriton2He3piML_custom_summary")[11]
cent_dist_max = cent_dist.GetMaximum()

##################################################################
# functions
func = {}
func_max = {}

func_names = ["BGBW", "Boltzmann", "Mt-exp", "Pt-exp", "LevyTsallis"]
func_names_MB = ["BlastWave", "Boltzmann", "LevyTsallis", "Mt-exp"]

# functions input files
input_func_file = ROOT.TFile("utils/Anti_fits.root")
input_func_file_MB = ROOT.TFile("utils/BlastWaveFits.root")

# get functions and maxima from file

for i_cent in range(len(cent_bins_new)):
    for i_fun in range(len(func_names)):
        cent_bins = cent_bins_new[i_cent]
        key = f"cent_{cent_bins[0]}_{cent_bins[1]}_func_{func_names[i_fun]}"
        func[key] = input_func_file.Get(f"{func_names[i_fun]}/{i_cent + 1}/{func_names[i_fun]}{i_cent + 1}")
        func_max[key] = func[key].GetMaximum()

for i_cent in range(len(cent_bins_MB)):
    for i_fun in range(len(func_names_MB)):
        cent_bins = cent_bins_MB[i_cent]
        key = f"cent_{cent_bins[0]}_{cent_bins[1]}_func_{func_names_MB[i_fun]}"
        func[key] = input_func_file_MB.Get(f"{func_names_MB[i_fun]}/{func_names_MB[i_fun]}{i_cent}")
        func_max[key] = func[key].GetMaximum()

# book histograms
h_abs_radius = {}
h_abs_ct = {}
h_gen_radius = {}
h_gen_ct = {}
h_gen_pt = {}
h_rec_radius = {}
h_rec_ct = {}
h_rec_pt = {}

for key in func.keys():
        for split in split_list:
            h_abs_radius[f"{split}_" + key] = ROOT.TH1D(f"fAbsRadius_{split}_{cent_bins[0]}_" + key, ";#it{R}_{#it{abs}} (cm);Entries", 1000, 0, 1000)
            h_abs_ct[f"{split}_" + key] = ROOT.TH1D(f"fAbsCt_{split}_" + key, ";#it{c}t (cm);Entries",50, 1, 35)
            h_gen_radius[f"{split}_" + key] = ROOT.TH1D(f"fGenRadius_{split}_" + key, ";#it{R}_{#it{abs}} (cm);Entries", 1000, 0, 1000)
            h_gen_ct[f"{split}_" + key] = ROOT.TH1D(f"fGenCt_{split}_" + key, ";#it{c}t (cm);Entries", 50, 1, 35)
            h_rec_radius[f"{split}_" + key] = ROOT.TH1D(f"fRecRadius_{split}_" + key, ";#it{R}_{#it{abs}} (cm);Entries", 1000, 0, 1000)
            h_rec_ct[f"{split}_" + key] = ROOT.TH1D(f"fRecCt_{split}_" + key, ";#it{c}t (cm);Entries", 50, 1, 35)

### allocate pt histograms
for key in h_rec_ct.keys():
    pt_bins = pt_bins_generator(key)
    if type(pt_bins) != list:
            pt_bins = PT_BINS_CENT[0]
            pt_bins = np.unique(np.array(pt_bins, dtype=float))
    pt_bins = np.unique(np.array(pt_bins, dtype=float))
    h_rec_pt[key] = ROOT.TH1D(f"fRecPt_" + key, ";#it{p}_{T} (GeV/#it{c});Entries", len(pt_bins) - 1, pt_bins)
    h_gen_pt[key] = ROOT.TH1D(f"fGenPt_" + key, ";#it{p}_{T} (GeV/#it{c});Entries", len(pt_bins) - 1, pt_bins)

# read tree
data_frame_he3 = ROOT.RDataFrame('STree', mc_file)
data_frame_he3 = data_frame_he3.Filter('pt > 2. and pt < 10. and (flag & 1)==1')
np_he3 = data_frame_he3.AsNumpy(["pt", "pdg", "absCt", "eta"])

# analysis in centrality classes
counter = 0
print_step_index = 0
num_entries = len(np_he3["pt"])
print_steps = num_entries*np.arange(0,1,0.01)


for he3 in zip(np_he3['pt'], np_he3['pdg'], np_he3['absCt'], np_he3['eta']):

    # if counter > 10000:
    #     break
    if np.floor(counter/num_entries*100) < 99:
        if counter > print_steps[print_step_index]:
            print("Loading.... : ", np.floor(counter/num_entries*100), " %")
            print_step_index +=1

    split = "antimatter"
    if he3[1] == 1000020030:
        split = "matter"
    absCt = he3[2]

    for key in func.keys():
        # rejection sampling to reweight pt
        if ROOT.gRandom.Rndm()*func_max[key] > func[key].Eval(he3[0]):
            continue
        # sample decay ct and ckeck for absorption
        decCt = ROOT.gRandom.Exp(7.6)
        # polar angle from eta
        tmp = abs(he3[3])
        tmp = ROOT.TMath.Exp(tmp)
        theta = 2*ROOT.TMath.ATan(tmp) # eta = -log[tan(theta/2)]
        # momentum from transverse momentum and angle
        mom = he3[0]/ROOT.TMath.Sin(theta)
        # absorption radius
        abs_radius = absCt*mom/HE_3_MASS
        # decay radius
        dec_radius = decCt*mom/HE_3_MASS

        h_abs_ct[f"{split}_" + key].Fill(absCt)
        h_abs_radius[f"{split}_" + key].Fill(abs_radius)
        h_gen_radius[f"{split}_" + key].Fill(dec_radius)
        h_gen_ct[f"{split}_" + key].Fill(decCt)
        h_gen_pt[f"{split}_" + key].Fill(he3[0])
        if not (decCt > absCt):  # decCt < absCt
            h_rec_radius[f"{split}_" + key].Fill(dec_radius)
            h_rec_ct[f"{split}_" + key].Fill(decCt)
            h_rec_pt[f"{split}_" + key].Fill(he3[0])
        if (absCt < -0.5):  # decCt < absCt
            h_rec_radius[f"{split}_" + key].Fill(dec_radius)
            h_rec_ct[f"{split}_" + key].Fill(decCt)
            h_rec_pt[f"{split}_" + key].Fill(he3[0])

        # sample decay ct and check for absorption
        decCt = ROOT.gRandom.Exp(7.6)
        h_gen_ct[f"{split}_" + key].Fill(decCt)
        h_gen_pt[f"{split}_" + key].Fill(he3[0])
        if (decCt < absCt) or (absCt < -0.5):  # decCt < absCt
            h_rec_ct[f"{split}_" + key].Fill(decCt)
            h_rec_pt[f"{split}_" + key].Fill(he3[0])
    counter +=1
for key in h_rec_ct.keys():
    key_cent = re.findall(r"[-+]?\d*\.\d+|\d+", key)
    if not outfile.GetDirectory(f"{key_cent[0]}_{key_cent[1]}"):
        outfile.mkdir(f"{key_cent[0]}_{key_cent[1]}")
    outfile.cd(f"{key_cent[0]}_{key_cent[1]}")

    ############### eff radius
    h_rec_radius[key].Divide(h_gen_radius[key])
    h_rec_radius[key].GetXaxis().SetTitle("#it{R}_{#it{abs}} (cm)")
    h_rec_radius[key].GetYaxis().SetTitle("1 - #it{f}_{abs}")
    h_rec_radius[key].Write(f"fEffRadius_" + key)  
    ############### eff ct
    h_rec_ct[key].Divide(h_gen_ct[key])
    h_rec_ct[key].GetXaxis().SetTitle("#it{c}t (cm)")
    h_rec_ct[key].GetYaxis().SetTitle("1 - #it{f}_{abs}")
    h_rec_ct[key].Write(f"fEffCt_" + key)
    ############### eff pT
    h_rec_pt[key].Divide(h_gen_pt[key])
    h_rec_pt[key].GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
    h_rec_pt[key].GetYaxis().SetTitle("1 - #it{f}_{abs}")
    h_rec_pt[key].Write(f"fEffPt_" + key)

outfile.Close()

from os import PRIO_USER
import uproot 
import numpy as np
import pandas as pd
from hipe4ml.tree_handler import TreeHandler
import matplotlib.pyplot as plt
import ROOT

import pickle


def apply_pt_rejection(df, pt_shape):
    rej_flag = -1*np.ones(len(df))
    random_arr = np.random.rand(len(df))
    max_bw = pt_shape.GetMaximum()

    for ind, (pt, rand) in enumerate(zip(df['gPt'],random_arr)):
        frac = pt_shape.Eval(pt)/max_bw
        if rand < frac:
            rej_flag[ind] = True
    df['rej'] = rej_flag


hdl = TreeHandler('/data/fmazzasc/PbPb_2body/SignalTable_20g7_flat_pt.root', 'SignalTable')
df = hdl.get_data_frame()
print(np.min(df.query('pt>0')))

bw_file = ROOT.TFile('../utils/BlastWaveFits.root')
bw = bw_file.Get("BlastWave/BlastWave0")

apply_pt_rejection(df, bw)

plt.hist(df.query('rej>0')['gPt'], bins = 200, alpha = 0.5)
plt.hist(df['gPt'], bins = 200, alpha = 0.5)

plt.savefig('rej.png')

print(df.query('rej>0'))

print(len(df.query('pt>0 and abs(Rapidity)<0.5'))/len(df.query('abs(gRapidity)<0.5')))
print(len(df.query('rej>0 and pt>0 and abs(Rapidity)<0.5'))/len(df.query('rej>0 and abs(gRapidity)<0.5')))


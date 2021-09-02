import uproot 
import numpy as np
import pandas as pd
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
import matplotlib.pyplot as plt

import pickle

df = TreeHandler('/data/fmazzasc/PbPb_2body/MC_tables/SignalTable_20g7_flat_pt.root', 'SignalTable')
model_hdl = ModelHandler()
model_hdl.load_model_handler(f'../models/2018_single_bin/2_9_trained')
dict = pickle.load(open('../results/2018_single_bin/file_score_eff_dict','rb'))

print(dict)
score_cut = 3.48671889

df.apply_preselections('centrality<10 and 2<pt<9 and Matter==1', inplace=True)
df.apply_model_handler(model_hdl)

print(len(df.apply_preselections('model_output>3.48671889', inplace=False))/len(df))

df_sel = df.apply_preselections('model_output>5', inplace=False)


plt.hist(df_sel['V0CosPA'], bins = 200, alpha=0.6, range = [0.999,1], density = True)
plt.hist(df['V0CosPA'], bins = 200, alpha=0.6, range = [0.999,1], density = True)
plt.yscale('log')
plt.savefig('cospa.png')
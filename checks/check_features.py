import sys
sys.path.append('../utils')
import helpers as hp
import uproot 
import numpy as np
import pandas as pd
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
import hipe4ml.plot_utils as pu

import matplotlib.pyplot as plt

training_columns = ['V0CosPA', 'ct', 'ProngsDCA','PiProngPvDCAXY', 'He3ProngPvDCAXY', 'He3ProngPvDCA', 'PiProngPvDCA', 'NpidClustersHe3', 'TPCnSigmaHe3','TPCnSigmaPi']
df = TreeHandler('/data/fmazzasc/PbPb_2body/SignalTable_20g7_flat_pt.root', 'SignalTable')
df_old = TreeHandler('/data/fmazzasc/PbPb_2body/SignalTable_16h7abc_flat_pt.root', 'SignalTable')


pu.plot_distr([df.apply_preselections('centrality < 5', inplace=False), df.apply_preselections('centrality > 50', inplace=False)], column=training_columns,labels=["cent < 5", "cent > 50"], log=True, density=True, figsize=(20, 20), alpha=0.3, grid=False)
plt.savefig('check_centrality.png')
pu.plot_distr([df, df_old], column=training_columns,labels=["20g7", "16h7abc"], log=True, density=True, figsize=(20, 20), alpha=0.3, grid=False)
plt.savefig('mc_comp.png')

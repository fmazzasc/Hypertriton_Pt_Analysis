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
df_sig = TreeHandler('/data/fmazzasc/PbPb_2body/KF/AnalysisResults_MC.root', 'HyperTree', columns_names=training_columns)
print('Signal loaded')
df_dat = TreeHandler('/data/fmazzasc/PbPb_2body/KF/AnalysisResults.root', 'HyperTree', columns_names=training_columns)
print('Data loaded')

df_sig.apply_preselections('V0CosPA>0.99 and ct<50 and NpidClustersHe3>100')
df_dat.apply_preselections('V0CosPA>0.99 and ct<50 and NpidClustersHe3>100')



pu.plot_distr([df_sig, df_dat], column=training_columns,labels=["Signal", "Data"], log=True, density=True, figsize=(20, 20), alpha=0.3, grid=False)
plt.savefig('mc_comp.png')

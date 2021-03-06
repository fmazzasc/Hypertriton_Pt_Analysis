
import sys
sys.path.append('../utils')
import helpers as hp
import uproot 
import numpy as np
import pandas as pd
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
import matplotlib.pyplot as plt
import pickle
import ROOT


dict = pickle.load(open('../results/bins_offline_2018/file_score_eff_dict','rb'))
score = dict['matter_0_5_3_4'][0]
eff_arr = dict['matter_0_5_3_4'][1]

df = pd.read_parquet('TreeDataFrame.parquet.gzip')

model_hdl = ModelHandler()
model_hdl.load_model_handler('../models/bins_offline_2018/matter_0_5_3_4_trained')

print(len(df))
df.query('centrality < 5 and 3<pt<4 and rej>0 and Matter==1 and ct<35', inplace=True)
print(len(df))
model_output = np.array(model_hdl.predict(df))

cut = score[60]
eff = eff_arr[60]

print(np.sum(model_output>cut)/len(df), " vs ", eff)

import pandas as pd
import matplotlib.pyplot as plt
import uproot

df = uproot.open('/data/fmazzasc/PbPb_2body/KF/AnalysisResults.root')['HyperTree'].arrays(library='pd')
df = df.query('2<pt<3 and NpidClustersHe3>130 and Matter>0')
plt.hist(df['TPCnSigmaHe3'], bins = 100)
plt.savefig('prova.png')
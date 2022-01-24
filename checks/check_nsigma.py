import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_parquet('/data/fmazzasc/PbPb_2body/no_pt_cut/DataTable_18qr.parquet.gzip')
df = df.query('2<pt<3 and NpidClustersHe3>130 and Matter>=0')
plt.hist(df['TPCnSigmaHe3'], bins = 100)
plt.savefig('prova.png')
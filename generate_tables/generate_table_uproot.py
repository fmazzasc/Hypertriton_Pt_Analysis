import uproot
import awkward as ak
import numpy as np
import pandas as pd

from concurrent.futures import ThreadPoolExecutor
executor = ThreadPoolExecutor(20)

kHyperMass = 2.99131
kHe3Mass = 2.80839160743
kPiMass = 0.13957


def remove_col_dot(df):
    new_col = []
    for col in df.columns:
        new_col.append(col[col.find('.') + 1:])
    return new_col


def Hypote(a, b, c=0):
    return np.sqrt(a**2 + b**2 + c**2)


cols = ['fCent', 'fTrigger', 'RHyperTriton.fPxHe3', 'RHyperTriton.fPyHe3', 'RHyperTriton.fPzHe3', 'RHyperTriton.fPxPi',
        'RHyperTriton.fPyPi', 'RHyperTriton.fPzPi', 'RHyperTriton.fDecayX', 'RHyperTriton.fDecayY', 'RHyperTriton.fDecayZ',
        'RHyperTriton.fDcaV0daughters', 'RHyperTriton.fDcaHe32PrimaryVertex', 'RHyperTriton.fDcaPi2PrimaryVertex',
        'RHyperTriton.fDcaHe32PrimaryVertexXY', 'RHyperTriton.fDcaPi2PrimaryVertexXY', 'RHyperTriton.fNpidClustersHe3',
        'RHyperTriton.fTPCnSigmaHe3', 'RHyperTriton.fTPCnSigmaPi']


df = uproot.open('HyperTritonTree_18q.root', decompression_executor=executor)[
    '_custom']['fTreeV0'].arrays(cols, library="pd", decompression_executor=executor)

df.columns = remove_col_dot(df)

he3_en = np.sqrt(df['fPxHe3']**2 + df['fPyHe3'] **
                 2 + df['fPzHe3']**2 + kHe3Mass**2)
pi_en = np.sqrt(df['fPxPi']**2 + df['fPyPi']**2 + df['fPzPi']**2 + kPiMass**2)
hyp_en = he3_en + pi_en
hyp_p = [df['fPxHe3'] + df['fPxPi'], df['fPyHe3'] +
         df['fPyPi'], df['fPzHe3'] + df['fPzPi']]
hyp_dec = [df["fDecayX"], df["fDecayY"], df["fDecayZ"]]


hyp_pt = np.sqrt(hyp_p[0]**2 + hyp_p[1]**2)
hyp_mass = np.sqrt(hyp_en**2 - hyp_p[0]**2 - hyp_p[1]**2 - hyp_p[2]**2)
hyp_ct = kHyperMass * \
    (Hypote(hyp_dec[0], hyp_dec[1], hyp_dec[2]) /
     Hypote(hyp_p[0], hyp_p[1], hyp_p[2]))
hyp_rap = 0.5*np.log((hyp_en+hyp_p[2])/(hyp_en - hyp_p[2]))

p_dec_prod = hyp_p[0]*df["fDecayX"] + hyp_p[1] * \
    df["fDecayY"] + hyp_p[2]*df["fDecayZ"]
cpa = p_dec_prod/(Hypote(hyp_p[0], hyp_p[1], hyp_p[2])
                  * Hypote(hyp_dec[0], hyp_dec[1], hyp_dec[2]))

new_df = pd.DataFrame()
new_df['centrality'] = df['fCent']
new_df['trigger'] = df['fTrigger']
new_df['PiProngPt'] = Hypote(df['fPxPi'], df['fPyPi'])
new_df['He3ProngPt'] = Hypote(df['fPxHe3'], df['fPyHe3'])
new_df['ProngsDCA'] = df['fDcaV0daughters']
new_df['PiProngPvDCA'] = df['fDcaPi2PrimaryVertex']
new_df['He3ProngPvDCA'] = df['fDcaHe32PrimaryVertex']
new_df['PiProngPvDCAXY'] = df['fDcaPi2PrimaryVertexXY']
new_df['He3ProngPvDCAXY'] = df['fDcaHe32PrimaryVertexXY']
new_df['NpidClustersHe3'] = df['fNpidClustersHe3']
new_df['TPCnSigmaPi'] = df['fTPCnSigmaPi']
new_df['TPCnSigmaHe3'] = df['fTPCnSigmaHe3']
new_df['pt'] = hyp_pt
new_df['ct'] = hyp_ct
new_df['m'] = hyp_mass
new_df['V0CosPA'] = cpa
new_df['Rapidity'] = hyp_rap


new_df.reset_index(drop=True, inplace=True)
new_df.to_parquet('DataTable_18q.parquet.gzip')
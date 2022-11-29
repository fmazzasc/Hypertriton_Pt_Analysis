import uproot 

df_mc = uproot.open('/data/fmazzasc/PbPb_2body/KF/AnalysisResults_MC.root')['HyperTree'].arrays(library='pd')
df_mc = df_mc.query('3<pt<4 and isReconstructed==True')

df_mc_2 = uproot.open('/data/fmazzasc/PbPb_2body/AOD/HyperTritonTree_MC.root')['HyperTree'].arrays(library='pd')
df_mc_2 = df_mc_2.query('3<pt<4 and isReconstructed==True')

#compare masses of MC and MC KF
import matplotlib.pyplot as plt
plt.hist(df_mc['m'], bins = 100, alpha=0.5, range=[2.96, 3.02], label='MC KF')
plt.hist(df_mc_2['m'], bins = 100, alpha=0.5, range=[2.96, 3.02], label='MC AOD')
plt.legend()
plt.savefig('mass.png')


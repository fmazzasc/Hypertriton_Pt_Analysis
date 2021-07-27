

import os
from ROOT import gROOT

gROOT.SetBatch(True)

gROOT.LoadMacro("GenerateTableFromMC.cc")
gROOT.LoadMacro("GenerateTableFromData.cc")
from ROOT import GenerateTableFromMC, GenerateTableFromData


print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate Data Table")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromData(False)
print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate Like-Sign Background Table")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromData(True)

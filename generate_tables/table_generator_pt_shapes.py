#!/usr/bin/env python3

import os
from ROOT import gROOT

gROOT.SetBatch(True)

gROOT.LoadMacro("GenerateTableFromMC.cc")
gROOT.LoadMacro("GenerateTableFromData.cc")
from ROOT import GenerateTableFromMC, GenerateTableFromData

GenerateTableFromMC(True, "bw", "0_10")
print("++++++++++++++++++++++++++++++++++++++++++")
# GenerateTableFromMC(True, "mtexp", "10_30")
# print("++++++++++++++++++++++++++++++++++++++++++")
# GenerateTableFromMC(True, "mtexp", "30_50")
# print("++++++++++++++++++++++++++++++++++++++++++")
# GenerateTableFromMC(True, "mtexp", "50_90")
# print("++++++++++++++++++++++++++++++++++++++++++")


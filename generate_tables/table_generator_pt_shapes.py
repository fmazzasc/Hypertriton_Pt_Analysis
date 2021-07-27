#!/usr/bin/env python3

import os
from ROOT import gROOT

gROOT.SetBatch(True)

gROOT.LoadMacro("GenerateTableFromMC.cc")
gROOT.LoadMacro("GenerateTableFromData.cc")
from ROOT import GenerateTableFromMC, GenerateTableFromData

GenerateTableFromMC(True, "levy", "0_10")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromMC(True, "levy", "10_30")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromMC(True, "levy", "30_50")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromMC(True, "levy", "50_90")
print("++++++++++++++++++++++++++++++++++++++++++")


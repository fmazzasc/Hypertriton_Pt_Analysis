import sys
sys.path.append('../utils')
import helpers as hp

import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
import helpers_std as hp_std
import ROOT
from hipe4ml.tree_handler import TreeHandler
ROOT.gROOT.SetBatch()

ROOT.gROOT.ProcessLine(".L ../utils/alipwgfunc/AliPWGFunc.cxx++")
from ROOT import AliPWGFunc


multiplicity = [1756, 983, 415]

our_yield = 

stefano_yield = 

he3_yield = 
#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".L tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

import os
objs = []
#mass = 651
mass = 1599
samples = [
    "ntuple/DYJetsToLL_M-50_noPU.root",
    "ntuple/HSCPppstau_m%d_LGW25.root" % mass,
    "ntuple/DYJetsToLL_M-50_PU200*.root",
    "ntuple/HSCPppstau_M_%d_PU200*.root" % mass,
]

gSystem.CompileMacro("TreeAnalyzer.C", "k");

for sample in samples:
    foutName = "hist_" + os.path.basename(sample.replace('*',''))

    tree = TChain("HSCPTree/tree")
    tree.Add(sample)
    ana = TreeAnalyzer(tree)
    fout = TFile(foutName, "RECREATE")
    ana.Loop(fout)
    ana = None


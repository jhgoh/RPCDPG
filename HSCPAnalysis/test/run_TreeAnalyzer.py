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
#fName1 = "ntuple/DYJetsToLL_M-50_noPU.root"
#fName2 = "ntuple/HSCPppstau_m%d_LGW25.root" % mass
fName1 = "ntuple/DYJetsToLL_M-50_PU200*.root"
fName2 = "ntuple/HSCPppstau_M_%d_PU200*.root" % mass
title1 = "Z#rightarrow#mu#mu"
title2 = "#tilde{#tau}^{-} (M=%d GeV)" % mass

gSystem.CompileMacro("TreeAnalyzer.C", "k");

for sample in [fName1, fName2]:
    foutName = "hist_" + os.path.basename(sample.replace('*',''))

    tree = TChain("HSCPTree/tree")
    tree.Add(sample)
    ana = TreeAnalyzer(tree)
    fout = TFile(foutName, "RECREATE")
    ana.Loop(fout)
    ana = None


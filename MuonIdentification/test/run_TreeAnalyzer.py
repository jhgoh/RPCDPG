#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".L tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

import os
objs = []
fName = "RelValZMM_14_NoPU.root"
#fName = "ntuple/ZMM_14_D17NoPU.root"
#fName = "ntuple/ZMM_14_D17PU140.root"
#fName = "DYJetsToLL_M-50_NoPU.root"

def save(c): c.Print("%s.png" % c.GetName())

gSystem.CompileMacro("TreeAnalyzer.C", "k");

foutName = "hist_" + os.path.basename(fName)

if not os.path.exists(foutName):
    fin = TFile(fName)
    tree = fin.Get("HSCPTree/tree")
    ana = TreeAnalyzer(tree)
    fout = TFile(foutName, "RECREATE")
    ana.Loop(fout)
    ana = None
fout = TFile(foutName)

objs = []
for var in ["all", "loose", "tight", "RPC"]:
    cEff_abseta = TCanvas("c_abseta_%s" % var, "Eff. #abseta "+var, 500, 500)
    eff_abseta = fout.Get("%s/eff_abseta" % var)
    eff_abseta.Draw()

    cEff_phi = TCanvas("c_phi_%s" % var, "Eff. #phi "+var, 500, 500)
    eff_phi = fout.Get("%s/eff_phi" % var)
    eff_phi.Draw()

    #cdR = TCanvas("cDR_%s" % var, "#Delta{}R "+var+" dR", 500, 500)
    #hDR = fout.Get("%s/h_minDR" % var)
    #hDR.Draw()

    #objs.extend([cEff_abseta, eff_abseta, cEff_phi, eff_phi, cdR, hDR])
    objs.extend([cEff_abseta, eff_abseta, cEff_phi, eff_phi])


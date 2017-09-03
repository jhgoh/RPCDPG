#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".L tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

import os
objs = []
mass = 651
#mass = 1599
#res = [0,1,2,3,4,5,10]
res = 0
fName1 = "ntuple/DYJetsToLL_M-50_PU200*.root"
fName2 = "ntuple/HSCPppstau_M_%d_PU200*.root" % mass
title1 = "Z#rightarrow#mu#mu"
title2 = "#tilde{#tau}^{-} (M=%d GeV)" % mass

def normalize(h):
    if h.Integral() == 0: return
    h.Scale(1./h.Integral())

def save(c): c.Print("%s.png" % c.GetName())

def drawH1(varName, dir1, dir2, options):
    options = options.split(',')

    h1 = dir1.Get("h_%s" % varName)
    h2 = dir2.Get("h_%s" % varName)

    if 'normalize' in options:
        yTitle = "Normalized"
        normalize(h1)
        normalize(h2)
    else:
        yTitle = h1.GetYaxis().GetTitle()

    if 'stack' in options:
        h1.SetFillColor(kRed)
        h2.SetFillColor(kBlue+1)
    else:
        h1.SetLineColor(kRed)
        h2.SetLineColor(kBlue+1)

    h1.SetOption("hist")
    h2.SetOption("hist")

    hs = THStack("hs_%s" % varName, "%s;%s;%s" % (varName, h1.GetXaxis().GetTitle(), yTitle))
    hs.Add(h2)
    hs.Add(h1)

    leg = TLegend(0.2, 0.75, 0.5, 0.85, "", "NDC")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    leg.AddEntry(h1, title1, "l")
    leg.AddEntry(h2, title2, "l")

    c = TCanvas("c_%s" % varName, varName, 500, 500)
    if 'stack' in options: hs.Draw("hist")
    else: hs.Draw("nostack,hist")

    leg.Draw()

    save(c)

    objs.append([c, leg, hs, h1, h2])


gSystem.CompileMacro("TreeAnalyzer.C", "k");

for sample in [fName1, fName2]:
    foutName = "hist_" + os.path.basename(sample.replace('*',''))
    if os.path.exists(foutName): continue

    tree = TChain("HSCPTree/tree")
    tree.Add(sample)
    ana = TreeAnalyzer(tree)
    fout = TFile(foutName, "RECREATE")
    ana.Loop(fout)
    ana = None

f1 = TFile("hist_"+os.path.basename(fName1.replace('*','')))
f2 = TFile("hist_"+os.path.basename(fName2.replace('*','')))

dirGen1 = f1.Get("gen")
dirGen2 = f2.Get("gen")

drawH1("muon1_beta_res%03d" % res, dirGen1, dirGen2, "normalize")
drawH1("muon1_iRPC_beta_res%03d" % res, dirGen1, dirGen2, "normalize")
drawH1("muon1_cRPC_beta_res%03d" % res, dirGen1, dirGen2, "normalize")

##########################


dirMuon1 = f1.Get("muon")
dirMuon2 = f2.Get("muon")

for leg in ("muon1", "muon2", "iRPC_muon1", "iRPC_muon2", "cRPC_muon1", "cRPC_muon2"):
    for var in ("beta", "m"):
        drawH1("%s_%s" % (leg, var), dirMuon1, dirMuon2, "normalize")

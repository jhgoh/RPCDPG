#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".L tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

import os
objs = []

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

    leg.AddEntry(h1, "Z#rightarrow#mu#mu", "l")
    leg.AddEntry(h2, "#tilde{#tau}^{-} (M=%d GeV)" % mass, "l")

    c = TCanvas("c_%s" % varName, varName, 500, 500)
    if 'stack' in options: hs.Draw("hist")
    else: hs.Draw("nostack,hist")

    leg.Draw()

    save(c)

    objs.append([c, leg, hs, h1, h2])


gSystem.CompileMacro("TreeAnalyzer.C", "k");
#mass = 651
mass = 1599

fNameDY = "ntuple/PhaseIIFall16DR82/DYJetsToLL_M-50_NoPU.root"
fNameHSCP = "ntuple/PhaseIIFall16DR82/HSCPppstau_M_%d_NoPU.root" % mass

for sample in [fNameDY, fNameHSCP]:
    foutName = "hist_" + os.path.basename(sample)
    if os.path.exists(foutName): continue

    fin = TFile(sample)
    tree = fin.Get("HSCPTree/tree")
    ana = TreeAnalyzer(tree)
    fout = TFile(foutName, "RECREATE")
    ana.Loop(fout)
    ana = None

#res = [0,1,2,3,4,5,10]
res = 0

fDY = TFile("hist_"+os.path.basename(fNameDY))
fHSCP = TFile("hist_"+os.path.basename(fNameHSCP))

dirGenDY = fDY.Get("gen")
dirGenHSCP = fHSCP.Get("gen")

drawH1("muon1_beta_res%03d" % res, dirGenDY, dirGenHSCP, "normalize")
drawH1("muon1_iRPC_beta_res%03d" % res, dirGenDY, dirGenHSCP, "normalize")
drawH1("muon1_cRPC_beta_res%03d" % res, dirGenDY, dirGenHSCP, "normalize")

##########################


dirMuonDY = fDY.Get("muon")
dirMuonHSCP = fHSCP.Get("muon")

drawH1("muon1_beta", dirMuonDY, dirMuonHSCP, "normalize")
drawH1("muon2_beta", dirMuonDY, dirMuonHSCP, "normalize")
drawH1("muon1_m", dirMuonDY, dirMuonHSCP, "normalize")
drawH1("muon2_m", dirMuonDY, dirMuonHSCP, "normalize")

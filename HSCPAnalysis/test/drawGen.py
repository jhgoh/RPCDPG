#!/usr/bin/env python

from ROOT import *

gROOT.ProcessLine("#include \"tdrstyle.C\"")
gROOT.ProcessLine("setTDRStyle();")
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

f = TFile("HSCPppstau_M_1599_NoPU.root")
tree = f.Get("HSCPTree/tree")

h_gen1_eta = TH1D("h_gen1_eta", "gen1_eta;#tilde{#tau}^{-} pseudorapidity #eta;Events / 0.1", 60, -3, 3)
h_gen2_eta = TH1D("h_gen2_eta", "gen2_eta;#tilde{#tau}^{+} pseudorapidity #eta;Events / 0.1", 60, -3, 3)
h_gen1_eta_gen2_eta = TH2D("h_gen1_eta_gen2_eta", "gen1_eta_gen2_eta;#tilde{#tau}^{-} pseudorapidity #eta;#tilde{#tau}^{+} pseudorapidity #eta", 60, -3, 3, 60, -3, 3)

tree.Draw("gen1_eta>>h_gen1_eta", "", "goff")
tree.Draw("gen2_eta>>h_gen2_eta", "", "goff")
tree.Draw("gen2_eta:gen1_eta>>h_gen1_eta_gen2_eta", "", "goff")

h_gen1_eta.SetMinimum(0)
h_gen2_eta.SetMinimum(0)
h_gen1_eta.SetMaximum(1.2*max(h_gen1_eta.GetMaximum(), h_gen2_eta.GetMaximum()))
h_gen2_eta.SetMaximum(h_gen1_eta.GetMaximum())

c_gen1_eta = TCanvas("c_gen1_eta", "gen1_eta", 500, 500)
lines_gen1_eta = [
    TLine(-2.4,0,-2.4,h_gen1_eta.GetMaximum()), TLine(2.4,0,2.4,h_gen1_eta.GetMaximum()),
    TLine(-2.1,0,-2.1,h_gen1_eta.GetMaximum()), TLine(2.1,0,2.1,h_gen1_eta.GetMaximum()),
]
h_gen1_eta.Draw()
for l in lines_gen1_eta:
    l.SetLineColor(kRed)
    l.SetLineStyle(kDashed)
    l.Draw()

c_gen2_eta = TCanvas("c_gen2_eta", "gen2_eta", 500, 500)
h_gen2_eta.Draw()
for l in lines_gen1_eta:
    l.SetLineColor(kRed)
    l.SetLineStyle(kDashed)
    l.Draw()

c_gen1_eta_gen2_eta = TCanvas("c_gen1_eta_gen2_eta", "gen1_eta_gen2_eta", 500, 500)
lines_gen1_eta_gen2_eta = [
    TLine(-2.4,-3,-2.4,3), TLine(2.4,-3,2.4,3), TLine(-3,-2.4,3,-2.4), TLine(-3,2.4,3,2.4),
    TLine(-2.1,-3,-2.1,3), TLine(2.1,-3,2.1,3), TLine(-3,-2.1,3,-2.1), TLine(-3,2.1,3,2.1),
]
h_gen1_eta_gen2_eta.Draw("COLZ")
for l in lines_gen1_eta_gen2_eta:
    l.SetLineColor(kRed)
    l.Draw()

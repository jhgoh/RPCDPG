#!/usr/bin/env python

from ROOT import *

def save(can):
    can.Print(can.GetName()+".png")

def normalize(h):
    if h.Integral() == 0: return
    h.Scale(1./h.Integral())

def buildCut(cut):
    if cut == "": return basecut
    return "(%s) && (%s)" % (basecut, cut)

gROOT.ProcessLine("#include \"tdrstyle.C\"")
gROOT.ProcessLine("setTDRStyle();")
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

#f = TFile("HSCPppstau_M_651_NoPU.root")
#f = TFile("HSCPppstau_M_1599_NoPU.root")
f = TFile("DYJetsToLL_M-50_NoPU.root")
tree = f.Get("HSCPTree/tree")

basecut = TCut("gen1_pdgId != 0 && gen2_pdgId != 0")

### fraction by eta regions
h_regions = TH1D("h_regions", "regions;;Events", 3, 0, 3)
h_regions.GetXaxis().SetBinLabel(1, "C+C")
h_regions.GetXaxis().SetBinLabel(2, "C+I")
h_regions.GetXaxis().SetBinLabel(3, "I+I")
tree.Draw("(abs(gen1_eta)>=1.8)+(abs(gen2_eta)>=1.8)>>h_regions", buildCut("abs(gen1_eta)<2.4 && abs(gen2_eta)<2.4"), "goff")
normalize(h_regions)

c_regions = TCanvas("c_regions", "c_regions", 500, 500)
c_regions.SetLogy()
h_regions.Draw("hist")
save(c_regions)

### gen1 and gen2 eta plots
h_gen1_eta = TH1D("h_gen1_eta", "gen1_eta;#tilde{#tau}^{-} pseudorapidity #eta;Events / 0.1", 60, -3, 3)
h_gen2_eta = TH1D("h_gen2_eta", "gen2_eta;#tilde{#tau}^{+} pseudorapidity #eta;Events / 0.1", 60, -3, 3)
h_gen1_eta_gen2_eta = TH2D("h_gen1_eta_gen2_eta", "gen1_eta_gen2_eta;#tilde{#tau}^{-} pseudorapidity #eta;#tilde{#tau}^{+} pseudorapidity #eta", 60, -3, 3, 60, -3, 3)

tree.Draw("gen1_eta>>h_gen1_eta", buildCut(""), "goff")
tree.Draw("gen2_eta>>h_gen2_eta", buildCut(""), "goff")
tree.Draw("gen2_eta:gen1_eta>>h_gen1_eta_gen2_eta", buildCut(""), "goff")

h_gen1_eta.SetMinimum(0)
h_gen2_eta.SetMinimum(0)
h_gen1_eta.SetMaximum(1.2*max(h_gen1_eta.GetMaximum(), h_gen2_eta.GetMaximum()))
h_gen2_eta.SetMaximum(h_gen1_eta.GetMaximum())

c_gen1_eta = TCanvas("c_gen1_eta", "gen1_eta", 500, 500)
lines_gen1_eta = [
    TLine(-2.4,0,-2.4,h_gen1_eta.GetMaximum()), TLine(2.4,0,2.4,h_gen1_eta.GetMaximum()),
    TLine(-1.8,0,-1.8,h_gen1_eta.GetMaximum()), TLine(1.8,0,1.8,h_gen1_eta.GetMaximum()),
]
h_gen1_eta.Draw()
for l in lines_gen1_eta:
    l.SetLineColor(kRed)
    l.SetLineStyle(kDashed)
    l.Draw()
save(c_gen1_eta)

c_gen2_eta = TCanvas("c_gen2_eta", "gen2_eta", 500, 500)
h_gen2_eta.Draw()
for l in lines_gen1_eta:
    l.SetLineColor(kRed)
    l.SetLineStyle(kDashed)
    l.Draw()
save(c_gen2_eta)

c_gen1_eta_gen2_eta = TCanvas("c_gen1_eta_gen2_eta", "gen1_eta_gen2_eta", 500, 500)
lines_gen1_eta_gen2_eta = [
    TLine(-2.4,-3,-2.4,3), TLine(2.4,-3,2.4,3), TLine(-3,-2.4,3,-2.4), TLine(-3,2.4,3,2.4),
    TLine(-1.8,-3,-1.8,3), TLine(1.8,-3,1.8,3), TLine(-3,-1.8,3,-1.8), TLine(-3,1.8,3,1.8),
]
h_gen1_eta_gen2_eta.Draw("COLZ")
for l in lines_gen1_eta_gen2_eta:
    l.SetLineColor(kRed)
    l.Draw()
save(c_gen1_eta_gen2_eta)

### gen1 and gen2 beta plots
h_gen1_beta = TH1D("h_gen1_beta", "gen1_beta;#tilde{#tau}^{-} #beta;Events / 0.02", 50, 0, 1)
h_gen2_beta = TH1D("h_gen2_beta", "gen2_beta;#tilde{#tau}^{+} #beta;Events / 0.02", 50, 0, 1)
h_gen1_beta_gen2_beta = TH2D("h_gen1_beta_gen2_beta", "gen1_beta_gen2_beta;#tilde{#tau}^{-} #beta;#tilde{#tau}^{+} #beta", 50, 0, 1, 50, 0, 1)

tree.Draw("gen1_beta>>h_gen1_beta", buildCut(""), "goff")
tree.Draw("gen2_beta>>h_gen2_beta", buildCut(""), "goff")
tree.Draw("gen2_beta:gen1_beta>>h_gen1_beta_gen2_beta", buildCut(""), "goff")

h_gen1_beta.SetMinimum(0)
h_gen2_beta.SetMinimum(0)
h_gen1_beta.SetMaximum(1.2*max(h_gen1_beta.GetMaximum(), h_gen2_beta.GetMaximum()))
h_gen2_beta.SetMaximum(h_gen1_beta.GetMaximum())

c_gen1_beta = TCanvas("c_gen1_beta", "gen1_beta", 500, 500)
h_gen1_beta.Draw()
save(c_gen1_beta)

c_gen2_beta = TCanvas("c_gen2_beta", "gen2_beta", 500, 500)
h_gen2_beta.Draw()
save(c_gen2_beta)

c_gen1_beta_gen2_beta = TCanvas("c_gen1_beta_gen2_beta", "gen1_beta_gen2_beta", 500, 500)
h_gen1_beta_gen2_beta.Draw("COLZ")
save(c_gen1_beta_gen2_beta)

### gen1 and gen2 beta by regions
h_gen1_beta_iRPC = TH1D("h_gen1_beta_iRPC", "gen1_beta_iRPC;#tilde{#tau}^{-} #beta_iRPC;Events / 0.02", 50, 0, 1)
h_gen2_beta_iRPC = TH1D("h_gen2_beta_iRPC", "gen2_beta_iRPC;#tilde{#tau}^{+} #beta_iRPC;Events / 0.02", 50, 0, 1)
h_gen1_beta_cRPC = TH1D("h_gen1_beta_cRPC", "gen1_beta_cRPC;#tilde{#tau}^{-} #beta_cRPC;Events / 0.02", 50, 0, 1)
h_gen2_beta_cRPC = TH1D("h_gen2_beta_cRPC", "gen2_beta_cRPC;#tilde{#tau}^{+} #beta_cRPC;Events / 0.02", 50, 0, 1)
h_gen1_beta_nRPC = TH1D("h_gen1_beta_nRPC", "gen1_beta_nRPC;#tilde{#tau}^{-} #beta_cRPC;Events / 0.02", 50, 0, 1)
h_gen2_beta_nRPC = TH1D("h_gen2_beta_nRPC", "gen2_beta_nRPC;#tilde{#tau}^{+} #beta_cRPC;Events / 0.02", 50, 0, 1)

tree.Draw("gen1_beta>>h_gen1_beta_iRPC", buildCut("abs(gen1_eta)>=1.8 && abs(gen1_eta)<2.4"), "goff")
tree.Draw("gen2_beta>>h_gen2_beta_iRPC", buildCut("abs(gen2_eta)>=1.8 && abs(gen2_eta)<2.4"), "goff")
tree.Draw("gen1_beta>>h_gen1_beta_cRPC", buildCut("abs(gen1_eta)<1.8"), "goff")
tree.Draw("gen2_beta>>h_gen2_beta_cRPC", buildCut("abs(gen2_eta)<1.8"), "goff")
tree.Draw("gen1_beta>>h_gen1_beta_nRPC", buildCut("abs(gen1_eta)>=2.4"), "goff")
tree.Draw("gen2_beta>>h_gen2_beta_nRPC", buildCut("abs(gen2_eta)>=2.4"), "goff")

h_gen1_beta_iRPC_norm = h_gen1_beta_iRPC.Clone()
h_gen2_beta_iRPC_norm = h_gen2_beta_iRPC.Clone()
h_gen1_beta_cRPC_norm = h_gen1_beta_cRPC.Clone()
h_gen2_beta_cRPC_norm = h_gen2_beta_cRPC.Clone()
h_gen1_beta_nRPC_norm = h_gen1_beta_nRPC.Clone()
h_gen2_beta_nRPC_norm = h_gen2_beta_nRPC.Clone()

normalize(h_gen1_beta_iRPC_norm)
normalize(h_gen2_beta_iRPC_norm)
normalize(h_gen1_beta_cRPC_norm)
normalize(h_gen2_beta_cRPC_norm)
normalize(h_gen1_beta_nRPC_norm)
normalize(h_gen2_beta_nRPC_norm)

h_gen1_beta_iRPC.SetFillColor(kRed)
h_gen2_beta_iRPC.SetFillColor(kRed)
h_gen1_beta_cRPC.SetFillColor(kAzure+2)
h_gen2_beta_cRPC.SetFillColor(kAzure+2)
h_gen1_beta_nRPC.SetFillColor(kBlack)
h_gen2_beta_nRPC.SetFillColor(kBlack)

h_gen1_beta_iRPC_norm.SetLineColor(kRed)
h_gen2_beta_iRPC_norm.SetLineColor(kRed)
h_gen1_beta_cRPC_norm.SetLineColor(kAzure+2)
h_gen2_beta_cRPC_norm.SetLineColor(kAzure+2)
h_gen1_beta_nRPC_norm.SetLineColor(kBlack)
h_gen2_beta_nRPC_norm.SetLineColor(kBlack)

leg_gen1_beta_ciRPC = TLegend(0.2, 0.7, 0.5, 0.85, "", "NDC")
leg_gen1_beta_ciRPC.SetBorderSize(0)
leg_gen1_beta_ciRPC.SetFillStyle(0)

leg_gen2_beta_ciRPC = leg_gen1_beta_ciRPC.Clone()
leg_gen1_beta_ciRPC_norm = leg_gen1_beta_ciRPC.Clone()
leg_gen2_beta_ciRPC_norm = leg_gen1_beta_ciRPC.Clone()

c_gen1_beta_ciRPC = TCanvas("c_gen1_beta_ciRPC", "gen1_beta_ciRPC", 500, 500)
hs_gen1_beta_ciRPC = THStack("hs_gen1_beta_ciRPC", "hs_gen1_beta_ciRPC;#tilde{#tau}^{-} #beta;Events / 0.02")
hs_gen1_beta_ciRPC.Add(h_gen1_beta_nRPC)
hs_gen1_beta_ciRPC.Add(h_gen1_beta_iRPC)
hs_gen1_beta_ciRPC.Add(h_gen1_beta_cRPC)
leg_gen1_beta_ciRPC.AddEntry(h_gen1_beta_cRPC, "|#eta|<1.8", "f")
leg_gen1_beta_ciRPC.AddEntry(h_gen1_beta_iRPC, "1.8#leq|#eta|<2.4", "f")
leg_gen1_beta_ciRPC.AddEntry(h_gen1_beta_nRPC, "|#eta|>2.4", "f")
hs_gen1_beta_ciRPC.Draw()
leg_gen1_beta_ciRPC.Draw()
save(c_gen1_beta_ciRPC)

c_gen2_beta_ciRPC = TCanvas("c_gen2_beta_ciRPC", "gen2_beta_ciRPC", 500, 500)
hs_gen2_beta_ciRPC = THStack("hs_gen2_beta_ciRPC", "hs_gen2_beta_ciRPC;#tilde{#tau}^{+} #beta;Events / 0.02")
hs_gen2_beta_ciRPC.Add(h_gen2_beta_nRPC)
hs_gen2_beta_ciRPC.Add(h_gen2_beta_iRPC)
hs_gen2_beta_ciRPC.Add(h_gen2_beta_cRPC)
leg_gen2_beta_ciRPC.AddEntry(h_gen2_beta_cRPC, "|#eta|<1.8", "f")
leg_gen2_beta_ciRPC.AddEntry(h_gen2_beta_iRPC, "1.8#leq|#eta|<2.4", "f")
leg_gen2_beta_ciRPC.AddEntry(h_gen2_beta_nRPC, "|#eta|>2.4", "f")
hs_gen2_beta_ciRPC.Draw()
leg_gen2_beta_ciRPC.Draw()
save(c_gen2_beta_ciRPC)

c_gen1_beta_ciRPC_norm = TCanvas("c_gen1_beta_ciRPC_norm", "gen1_beta_ciRPC_norm", 500, 500)
hs_gen1_beta_ciRPC_norm = THStack("hs_gen1_beta_ciRPC_norm", "hs_gen1_beta_ciRPC_norm;#tilde{#tau}^{-} #beta;Normalized")
hs_gen1_beta_ciRPC_norm.Add(h_gen1_beta_nRPC_norm)
hs_gen1_beta_ciRPC_norm.Add(h_gen1_beta_iRPC_norm)
hs_gen1_beta_ciRPC_norm.Add(h_gen1_beta_cRPC_norm)
leg_gen1_beta_ciRPC_norm.AddEntry(h_gen1_beta_cRPC_norm, "|#eta|<1.8", "l")
leg_gen1_beta_ciRPC_norm.AddEntry(h_gen1_beta_iRPC_norm, "1.8#leq|#eta|<2.4", "l")
leg_gen1_beta_ciRPC_norm.AddEntry(h_gen1_beta_nRPC_norm, "|#eta|>2.4", "l")
hs_gen1_beta_ciRPC_norm.Draw("nostackhist")
leg_gen1_beta_ciRPC_norm.Draw()
save(c_gen1_beta_ciRPC_norm)

c_gen2_beta_ciRPC_norm = TCanvas("c_gen2_beta_ciRPC_norm", "gen2_beta_ciRPC_norm", 500, 500)
hs_gen2_beta_ciRPC_norm = THStack("hs_gen2_beta_ciRPC_norm", "hs_gen2_beta_ciRPC_norm;#tilde{#tau}^{+} #beta;Normalized")
hs_gen2_beta_ciRPC_norm.Add(h_gen2_beta_nRPC_norm)
hs_gen2_beta_ciRPC_norm.Add(h_gen2_beta_iRPC_norm)
hs_gen2_beta_ciRPC_norm.Add(h_gen2_beta_cRPC_norm)
leg_gen2_beta_ciRPC_norm.AddEntry(h_gen2_beta_cRPC_norm, "|#eta|<1.8", "l")
leg_gen2_beta_ciRPC_norm.AddEntry(h_gen2_beta_iRPC_norm, "1.8#leq|#eta|<2.4", "l")
leg_gen2_beta_ciRPC_norm.AddEntry(h_gen2_beta_nRPC_norm, "|#eta|>2.4", "l")
hs_gen2_beta_ciRPC_norm.Draw("nostackhist")
leg_gen2_beta_ciRPC_norm.Draw()
save(c_gen2_beta_ciRPC_norm)


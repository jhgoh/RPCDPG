#!/usr/bin/env python
from ROOT import *
gROOT.ProcessLine(".L tdrstyle.C")
gROOT.ProcessLine("setTDRStyle()")
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

def buildLegend(x1, y1, x2, y2):
    leg = TLegend(x1, y1, x2, y2)
    leg.SetLineWidth(0)
    leg.SetFillStyle(0)
    return leg

#fB = TFile("hist_DYJetsToLL_M-50_noPU.root")
#fS = TFile("hist_HSCPppstau_m1599_LGW25.root")
fB = TFile("hist_DYJetsToLL_M-50_PU200.root")
fS = TFile("hist_HSCPppstau_M_1599_14TeV_PU200.root")

tB = fB.Get("tree")
tS = fS.Get("tree")

###
c_beta = TCanvas("c_beta", "c_beta", 500, 500)

hS_beta = TH1D("hS_beta", "hS_beta;Reconstructed #beta", 100, -1, 2)
tS.Draw("fit_beta1>>hS_beta", "gen1_pdgId!=0 && (fit_qual1 < 1e9 && fit_qual2 < 1e9)", "goff")
hB_beta = TH1D("hB_beta", "hB_beta;Reconstructed #beta", 100, -1, 2)
tB.Draw("fit_beta1>>hB_beta", "gen1_pdgId!=0 && (fit_qual1 < 1e9 && fit_qual2 < 1e9)", "goff")

hS_beta.Scale(1./hS_beta.GetEntries())
hB_beta.Scale(1./hB_beta.GetEntries())

hS_beta.SetLineColor(kBlue)
hB_beta.SetLineColor(kRed)

hstack_beta = THStack("hstack_beta", "hstack_beta;%s;Arbitrary unit" % hS_beta.GetXaxis().GetTitle())
hstack_beta.Add(hS_beta)
hstack_beta.Add(hB_beta)
hstack_beta.Draw("nostackhist")
leg = buildLegend(0.23, 0.6, 0.47, 0.8)
leg.AddEntry(hS_beta, "HSCP", "l")
leg.AddEntry(hB_beta, "Z+Jets", "l")
leg.Draw()

###
c_beta_eta = TCanvas("c_beta_eta", "c_beta_eta", 1000, 500)

hS_beta_eta = TH2D("hS_beta_eta", "hS_beta_eta;Generated #eta;Reconstructed #beta", 100, -2.5, 2.5, 100, -1, 2)
hB_beta_eta = TH2D("hB_beta_eta", "hB_beta_eta;Generated #eta;Reconstructed #beta", 100, -2.5, 2.5, 100, -1, 2)
tS.Draw("fit_beta1:gen1_p4.Eta()>>hS_beta_eta", "gen1_pdgId!=0 && (fit_qual1 < 1e9 && fit_qual2 < 1e9)", "goff")
tB.Draw("fit_beta1:gen1_p4.Eta()>>hB_beta_eta", "gen1_pdgId!=0 && (fit_qual1 < 1e9 && fit_qual2 < 1e9)", "goff")
c_beta_eta.Divide(2,1)
c_beta_eta.cd(1)
hS_beta_eta.Draw("COLZ")
c_beta_eta.cd(2)
hB_beta_eta.Draw("COLZ")

###
c_beta_Eff = TCanvas("c_beta_Eff", "c_beta_Eff", 1000, 500)
gStyle.SetOptTitle(1)

hS_beta_genTot = TH1D("hS_beta_genTot", "hS_beta_genTot;Total Generated #beta", 100, -1, 2) 
hB_beta_genTot = TH1D("hB_beta_genTot", "hB_beta_genTot;Total Generated #beta", 100, -1, 2)
hS_beta_gen = TH1D("hS_beta_gen", "hS_beta_gen;Generated #beta", 100, -1, 2)
hB_beta_gen = TH1D("hB_beta_gen", "hB_beta_gen;Generated #beta", 100, -1, 2)

#tS.Draw("gen1_p4.Beta()>>hS_beta_gen", "gen1_pdgId!=0 && (fit_qual1 < 1e9 && fit_qual2 < 1e9)", "goff")
#tB.Draw("gen1_p4.Beta()>>hB_beta_gen", "gen1_pdgId!=0 && (fit_qual1 < 1e9 && fit_qual2 < 1e9)", "goff")
tS.Draw("gen1_p4.Beta()>>hS_beta_gen", "gen1_pdgId!=0 && fit_qual1 < 1e9 && fit_nhit1 >= 3 && fit_qual1/fit_beta1 < 0.3 && fit_beta1 < 0.7", "goff")
tB.Draw("gen1_p4.Beta()>>hB_beta_gen", "gen1_pdgId!=0 && fit_qual1 < 1e9 && fit_nhit1 >= 3 && fit_qual1/fit_beta1 < 0.3 && fit_beta1 < 0.7", "goff")
#tS.Draw("gen1_p4.Beta()>>hS_beta_gen", "gen1_pdgId!=0 && fit_qual1 < 1e9 && fit_nhit1 >= 4", "goff")
#tB.Draw("gen1_p4.Beta()>>hB_beta_gen", "gen1_pdgId!=0 && fit_qual1 < 1e9 && fit_nhit1 >= 4", "goff")
tS.Draw("gen1_p4.Beta()>>hS_beta_genTot", "gen1_pdgId!=0", "goff")
tB.Draw("gen1_p4.Beta()>>hB_beta_genTot", "gen1_pdgId!=0", "goff")
c_beta_Eff.Divide(2,1)

Sig_betaEffCurve = TEfficiency(hS_beta_gen, hS_beta_genTot)
Sig_betaEffCurve.SetTitle("HSCP Generated #beta efficiency;#beta;#epsilon")
c_beta_Eff.cd(1)
Sig_betaEffCurve.Draw()

Bg_betaEffCurve = TEfficiency(hB_beta_gen, hB_beta_genTot)
Bg_betaEffCurve.SetTitle("DYJets Generated #beta efficiency;#beta;#epsilon")
c_beta_Eff.cd(2)
Bg_betaEffCurve.Draw()

###
c_2nd_beta_Eff = TCanvas("c_2nd_beta_Eff", "c_2nd_beta_Eff", 1000, 500)
gStyle.SetOptTitle(1)

hS2_beta_genTot = TH1D("hS2_beta_genTot", "hS2_beta_genTot;Total Generated #beta", 100, -1, 2)
hB2_beta_genTot = TH1D("hB2_beta_genTot", "hB2_beta_genTot;Total Generated #beta", 100, -1, 2)
hS2_beta_gen = TH1D("hS2_beta_gen", "hS2_beta_gen;Generated #beta", 100, -1, 2)
hB2_beta_gen = TH1D("hB2_beta_gen", "hB2_beta_gen;Generated #beta", 100, -1, 2)

#tS.Draw("gen1_p4.Beta()>>hS2_beta_gen", "gen1_pdgId!=0 && (fit_qual2 < 1e9 && fit_qual2 < 1e9)", "goff")
#tB.Draw("gen1_p4.Beta()>>hB2_beta_gen", "gen1_pdgId!=0 && (fit_qual2 < 1e9 && fit_qual2 < 1e9)", "goff")
tS.Draw("gen2_p4.Beta()>>hS2_beta_gen", "gen2_pdgId!=0 && fit_qual2 < 1e9 && fit_nhit2 >= 3 && fit_qual2/fit_beta2 < 0.3 && fit_beta2 < 0.7", "goff")
tB.Draw("gen2_p4.Beta()>>hB2_beta_gen", "gen2_pdgId!=0 && fit_qual2 < 1e9 && fit_nhit2 >= 3 && fit_qual2/fit_beta2 < 0.3 && fit_beta2 < 0.7", "goff")
#tS.Draw("gen2_p4.Beta()>>hS2_beta_gen", "gen2_pdgId!=0 && fit_qual2 < 1e9 && fit_nhit2 >= 4", "goff")
#tB.Draw("gen2_p4.Beta()>>hB2_beta_gen", "gen2_pdgId!=0 && fit_qual2 < 1e9 && fit_nhit2 >= 4", "goff")
tS.Draw("gen2_p4.Beta()>>hS2_beta_genTot", "gen2_pdgId!=0", "goff")
tB.Draw("gen2_p4.Beta()>>hB2_beta_genTot", "gen2_pdgId!=0", "goff")
c_2nd_beta_Eff.Divide(2,1)

Sig2_betaEffCurve = TEfficiency(hS2_beta_gen, hS2_beta_genTot)
Sig2_betaEffCurve.SetTitle("HSCP Generated #beta efficiency;#beta;#epsilon")
c_2nd_beta_Eff.cd(1)
Sig2_betaEffCurve.Draw()

Bg2_betaEffCurve = TEfficiency(hB2_beta_gen, hB2_beta_genTot)
Bg2_betaEffCurve.SetTitle("DYJets Generated #beta efficiency;#beta;#epsilon")
c_2nd_beta_Eff.cd(2)
Bg2_betaEffCurve.Draw()

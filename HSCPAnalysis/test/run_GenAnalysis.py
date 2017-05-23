#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".L tdrstyle.C")
setTDRStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

import os

def normalize(h):
    if h.Integral() == 0: return
    h.Scale(1./h.Integral())

def save(c): c.Print("%s.png" % c.GetName())

gSystem.CompileMacro("GenAnalysis.C", "k");
mass = 651

fNameDY = "DYJetsToLL_M-50_NoPU.root"
#fNameHSCP = "HSCPppstau_M_1599_NoPU.root"
fNameHSCP = "HSCPppstau_M_%d_NoPU.root" % mass

for sample in [fNameDY, fNameHSCP]:
    foutName = "hist_" + sample
    if os.path.exists(foutName): continue

    fin = TFile(sample)
    tree = fin.Get("HSCPTree/tree")
    ana = GenAnalysis(tree)
    fout = TFile(foutName, "RECREATE")
    ana.Loop(fout)
    ana = None

#res = [0,1,2,3,4,5,10]
res = 2

fDY = TFile("hist_"+fNameDY)
fHSCP = TFile("hist_"+fNameHSCP)

hDY_beta1_all = fDY.Get("h_beta1_res%03d" % res)
hHSCP_beta1_all = fHSCP.Get("h_beta1_res%03d" % res)
hDY_beta1_iRPC = fDY.Get("h_iRPC_beta1_res%03d" % res)
hHSCP_beta1_iRPC = fHSCP.Get("h_iRPC_beta1_res%03d" % res)
hDY_beta1_cRPC = fDY.Get("h_cRPC_beta1_res%03d" % res)
hHSCP_beta1_cRPC = fHSCP.Get("h_cRPC_beta1_res%03d" % res)

hDY_beta1_all.SetLineColor(kRed)
hHSCP_beta1_all.SetLineColor(kBlue+1)
hDY_beta1_iRPC.SetLineColor(kRed)
hHSCP_beta1_iRPC.SetLineColor(kBlue+1)
hDY_beta1_cRPC.SetLineColor(kRed)
hHSCP_beta1_cRPC.SetLineColor(kBlue+1)

leg_beta1_all = TLegend(0.2, 0.75, 0.5, 0.85, "", "NDC")
leg_beta1_all.SetBorderSize(0)
leg_beta1_all.SetFillStyle(0)

leg_beta1_iRPC = leg_beta1_all.Clone()
leg_beta1_cRPC = leg_beta1_all.Clone()

leg_beta1_all.AddEntry(hDY_beta1_all, "Z#rightarrow#mu#mu", "l")
leg_beta1_all.AddEntry(hHSCP_beta1_all, "#tilde{#tau}^{-} (M=%d GeV)" % mass, "l")
leg_beta1_iRPC.AddEntry(hDY_beta1_iRPC, "Z#rightarrow#mu#mu", "l")
leg_beta1_iRPC.AddEntry(hHSCP_beta1_iRPC, "#tilde{#tau}^{-} (M=%d GeV)" % mass, "l")
leg_beta1_cRPC.AddEntry(hDY_beta1_cRPC, "Z#rightarrow#mu#mu", "l")
leg_beta1_cRPC.AddEntry(hHSCP_beta1_cRPC, "#tilde{#tau}^{-} (M=%d GeV)" % mass, "l")

normalize(hDY_beta1_all)
normalize(hHSCP_beta1_all)
normalize(hDY_beta1_iRPC)
normalize(hHSCP_beta1_iRPC)
normalize(hDY_beta1_cRPC)
normalize(hHSCP_beta1_cRPC)

c_beta1_all = TCanvas("c_beta_all_res%03d" % res, "all beta", 500, 500)
hDY_beta1_all.Draw("hist")
hHSCP_beta1_all.Draw("samehist")
leg_beta1_all.Draw()
save(c_beta1_all)

c_beta1_iRPC = TCanvas("c_beta_iRPC_res%03d" % res, "iRPC beta", 500, 500)
hDY_beta1_iRPC.Draw("hist")
hHSCP_beta1_iRPC.Draw("samehist")
leg_beta1_iRPC.Draw()
save(c_beta1_iRPC)

c_beta1_cRPC = TCanvas("c_beta_cRPC_res%03d" % res, "cRPC beta", 500, 500)
hDY_beta1_cRPC.Draw("hist")
hHSCP_beta1_cRPC.Draw("samehist")
leg_beta1_cRPC.Draw()
save(c_beta1_cRPC)


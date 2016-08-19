#!/usr/bin/env python

from ROOT import *
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

modules = [
    "SingleMuPt100",
    "HSCPppstau_M_156", "HSCPppstau_M_200", "HSCPppstau_M_308", "HSCPppstau_M_494",
    "HSCPppstau_M_651", "HSCPppstau_M_1029", "HSCPppstau_M_1218", "HSCPppstau_M_1599",
    #"HSCPstop_M_100", "HSCPstop_M_200", "HSCPstop_M_400", "HSCPstop_M_600",
    #"HSCPstop_M_1000", "HSCPstop_M_1200", "HSCPstop_M_1600", "HSCPstop_M_2600",
]
colors = [
    kBlack,
    kRed+1, kOrange+1, kGreen+1, kAzure+1,
    kBlue+1, kMagenta-2, kMagenta+1, kRed+3
]

c = TCanvas("cM", "cM", 500, 500)
leg = TLegend(0.5, 0.5, 0.9, 0.9)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
files = []
hists = []
opt = ""
for i, module in enumerate(modules):
    f = TFile("hist_%s.root" % module)
    h = f.Get("tofAnalysis/hMuonM")
    h.SetLineColor(colors[i])
    h.Rebin(20)
    #h.Scale(1./h.GetEntries())
    #h.SetMinimum(1e-3)
    h.SetMaximum(200)

    h.Draw(opt)
    leg.AddEntry(h, module, "l")
                                       
    opt = "same"

    files.append(f)
    hists.append(h)

leg.Draw()
c.Print(c.GetName()+".png")

#!/usr/bin/env python

import sys, os
from ROOT import *

gROOT.ProcessLine(".L %s/src/SUSYBSMAnalysis/HSCP/test/ICHEP_Analysis/tdrstyle.C" % os.environ['CMSSW_RELEASE_BASE'])
gROOT.ProcessLine("setTDRStyle()")
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetLabelSize(.037, "XY")
gStyle.SetTitleSize(.05, "XY")
gStyle.SetTitleOffset(1.1, "X")
gStyle.SetPadTopMargin(1)

headerTitle = TLatex(gStyle.GetPadLeftMargin(), 0.92, "CMS Phase-II Simulation Work in Progress")
headerTitle.SetNDC()
headerTitle.SetTextSize(0.04)
headerTitle.SetTextFont(42)
histDir = "20161009/TOFAnalysis"
modules = ["HSCPppstau_M_156", "HSCPppstau_M_1029",]
timeReses = ["100ps", "1ns", "2ns", "4ns", "25ns"]
colors = [kRed+1, kBlack, kOrange+1, kGreen+1, kAzure+1, kBlue+1, kMagenta-2, kMagenta+1, kRed+3]

canvases = []
legs = []
files = []

for module in modules:
    c = TCanvas("c_%s" % module, module, 500, 500)
    leg = TLegend(0.7, 0.67, 0.95, 0.85)
    opt = ""
    for i, timeRes in enumerate(timeReses):
        f = TFile("%s/hist_%s__timeRes%s.root" % (histDir, module, timeRes))
        h = f.Get("tofAnalysis/hBetaRes")
        h.SetTitle(";#beta resolution (#beta-#beta_{SIM})/#beta_{SIM};Arbitrary Unit")

        h.Rebin(2)
        h.SetLineColor(colors[i])
        h.Scale(1./h.GetEntries())
        h.SetMaximum(0.6)
        h.SetMinimum(0)

        h.Draw(opt+"hist")
        leg.AddEntry(h, timeRes, "l")

        opt = "same"

        files.append(f)

    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.Draw()

    headerTitle.Draw()

    legs.append(leg)
    canvases.append(c)
    files.append(f)

    c.Print(c.GetName()+".png")
    c.Print(c.GetName()+".C")
    c.Print(c.GetName()+".pdf")

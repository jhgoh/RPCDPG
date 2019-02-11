#!/usr/bin/env python

from ROOT import *
from array import array
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetPalette(1)

f1 = TFile("HGCalStainless_default_me0_singleNu_RECO/ntuple.root")
f2 = TFile("HGCalStainless_modified_me0_singleNu_RECO/ntuple.root")

hArea = f1.Get("rpcHit/hArea")
h1 = f1.Get("rpcHit/hCounts")
h2 = f2.Get("rpcHit/hCounts")
nEvent1 = f1.Get("rpcHit/hEvents").GetBinContent(1)
nEvent2 = f1.Get("rpcHit/hEvents").GetBinContent(1)

bins = [1e-4]
while len(bins)<100:
    bins.append(bins[-1]*0.9)
bins = array('d', reversed(bins))
nbins = len(bins)-1
hRateRB1 = TH1F("hRateRB1", "hRateRB1;/cm^{2};Number of Rolls", nbins, bins)
hRateRB2 = TH1F("hRateRB2", "hRateRB2;/cm^{2};Number of Rolls", nbins, bins)
hRateRE1 = TH1F("hRateRE1", "hRateRE1;/cm^{2};Number of Rolls", nbins, bins)
hRateRE2 = TH1F("hRateRE2", "hRateRE2;/cm^{2};Number of Rolls", nbins, bins)
hRateRX1 = TH1F("hRateRX1", "hRateRX1;/cm^{2};Number of Rolls", nbins, bins)
hRateRX2 = TH1F("hRateRX2", "hRateRX2;/cm^{2};Number of Rolls", nbins, bins)

hRateVsRate = TH2D("hRateVsRate", "hRateVsRate;Rate in default geom. [/cm^{2}];Rate in modified geom. [/cm^{2}]", nbins, bins, nbins, bins)
hRateVsRateRing3 = TH2D("hRateVsRateRing3", "hRateVsRateRing3;Rate in default geom. [/cm^{2}];Rate in modified geom. [/cm^{2}]", nbins, bins, nbins, bins)

for b in range(1, hArea.GetNbinsX()+1):
    rollName = hArea.GetXaxis().GetBinLabel(b)
    area = hArea.GetBinContent(b)
    if area == 0: break
    rate1 = h1.GetBinContent(b)/area/nEvent1
    rate2 = h2.GetBinContent(b)/area/nEvent2
    if rollName.startswith("W"):
        hRateRB1.Fill(rate1)
        hRateRB2.Fill(rate2)
    elif "_R1_" in rollName:
        hRateRX1.Fill(rate1)
        hRateRX2.Fill(rate2)
    else:
        hRateRE1.Fill(rate1)
        hRateRE2.Fill(rate2)

    hRateVsRate.Fill(rate1, rate2)
    if "_R3_" in rollName:
        hRateVsRateRing3.Fill(rate1, rate2)

cRateRB = TCanvas("cRateRB", "cRateRB", 500, 500)
cRateRB.SetLogx()
hRateRB1.SetMinimum(0)
#hRateRB1.SetMaximum(100)
hRateRB1.SetLineColor(kBlue)
hRateRB2.SetLineColor(kRed)
hRateRB1.Draw()
hRateRB2.Draw("same")

cRateRE = TCanvas("cRateRE", "cRateRE", 500, 500)
cRateRE.SetLogx()
hRateRE1.SetMinimum(0)
#hRateRE1.SetMaximum(100)
hRateRE1.SetLineColor(kBlue)
hRateRE2.SetLineColor(kRed)
hRateRE1.Draw()
hRateRE2.Draw("same")

cRateRX = TCanvas("cRateRX", "cRateRX", 500, 500)
cRateRX.SetLogx()
hRateRX1.SetMinimum(0)
#hRateRX1.SetMaximum(100)
hRateRX1.SetLineColor(kBlue)
hRateRX2.SetLineColor(kRed)
hRateRX1.Draw()
hRateRX2.Draw("same")

cRateVsRate = TCanvas("cRateVsRate", "cRateVsRate", 500, 500)
cRateVsRate.SetLogx()
cRateVsRate.SetLogy()
hRateVsRate.Draw("COLZ")
fx = TF1("fx", "x", bins[0], bins[-1])
fx.SetLineWidth(1)
fx.SetLineColor(kBlack)
fx.Draw("same")

cRateVsRateRing3 = TCanvas("cRateVsRateRing3", "cRateVsRateRing3", 500, 500)
cRateVsRateRing3.SetLogx()
cRateVsRateRing3.SetLogy()
hRateVsRateRing3.Draw("COLZ")
#fx = TF1("fx", "x", bins[0], bins[-1])
#fx.SetLineWidth(1)
#fx.SetLineColor(kBlack)
fx.Draw("same")

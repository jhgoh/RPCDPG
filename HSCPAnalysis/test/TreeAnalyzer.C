#define TreeAnalyzer_cxx
#include "TreeAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include <iostream>

void TreeAnalyzer::Loop(TFile* fout)
{
  fout->cd();

  TDirectory* dirGenBeta = fout->mkdir("gen");
  dirGenBeta->cd();
  std::vector<TH1D*> h_all_beta1, h_iRPC_beta1, h_cRPC_beta1;
  std::vector<TH1D*> h_all_beta2, h_iRPC_beta2, h_cRPC_beta2;
  std::vector<int> betaRes = {0,1,2,3,4,5,10};
  for ( auto r : betaRes ) {
    TH1D* hh_beta1 = new TH1D(Form("h_muon1_beta_res%03d", r), Form("Beta distribution res=%d%%;#tilde{#tau}^{-} #beta;Events / 0.02", r), 50, 0, 1);
    TH1D* hh_iRPC_beta1 = new TH1D(Form("h_muon1_iRPC_beta_res%03d", r), Form("iRPC Beta distribution res=%d%%;#tilde{#tau}^{-} #beta;Events / 0.02", r), 50, 0, 1);
    TH1D* hh_cRPC_beta1 = new TH1D(Form("h_muon1_cRPC_beta_res%03d", r), Form("cRPC Beta distribution res=%d%%;#tilde{#tau}^{-} #beta;Events / 0.02", r), 50, 0, 1);
    TH1D* hh_beta2 = new TH1D(Form("h_muon2_beta_res%03d", r), Form("Beta distribution res=%d%%;#tilde{#tau}^{+} #beta;Events / 0.02", r), 50, 0, 1);
    TH1D* hh_iRPC_beta2 = new TH1D(Form("h_muon2_iRPC_beta_res%03d", r), Form("iRPC Beta distribution res=%d%%;#tilde{#tau}^{+} #beta;Events / 0.02", r), 50, 0, 1);
    TH1D* hh_cRPC_beta2 = new TH1D(Form("h_muon2_cRPC_beta_res%03d", r), Form("cRPC Beta distribution res=%d%%;#tilde{#tau}^{+} #beta;Events / 0.02", r), 50, 0, 1);
    h_all_beta1.push_back(hh_beta1);
    h_iRPC_beta1.push_back(hh_iRPC_beta1);
    h_cRPC_beta1.push_back(hh_cRPC_beta1);
    h_all_beta2.push_back(hh_beta2);
    h_iRPC_beta2.push_back(hh_iRPC_beta2);
    h_cRPC_beta2.push_back(hh_cRPC_beta2);
  }

  TDirectory* dirReco = fout->mkdir("muon");
  dirReco->cd();
  TH1D* h_muon1_m = new TH1D("h_muon1_m", "TOF mass;Leading #tilde{#tau} mass (GeV);Events / 10GeV", 200, 0, 2000);
  TH1D* h_muon2_m = new TH1D("h_muon2_m", "TOF mass;2nd leading #tilde{#tau} mass (GeV);Events / 10GeV", 200, 0, 2000);
  TH1D* h_muon1_beta = new TH1D("h_muon1_beta", "Beta distribution;Leading #tilde{#tau} #beta;Events / 0.02", 50, 0, 1);
  TH1D* h_muon2_beta = new TH1D("h_muon2_beta", "Beta distribution;2nd leading #tilde{#tau} #beta;Events / 0.02", 50, 0, 1);
  TH2D* h_muon1_beta__muon2_beta = new TH2D("h_muon1_beta__muon2_beta", "Beta distribution;Leading #tilde{#tau} #beta;2nd leading #tilde{#tau} #beta", 50, 0, 1, 50, 0, 1);
  TH2D* h_muon1_m__muon2_m = new TH2D("h_muon1_m__muon2_m", "TOF mass;Leading #tilde{#tau} mass (GeV);2nd leading #tilde{#tau} mass (GeV)", 50, 0, 1, 50, 0, 1);

  TH1D* h_iRPC_muon1_m = new TH1D("h_iRPC_muon1_m", "TOF mass;Leading #tilde{#tau} mass (GeV);Events / 10GeV", 200, 0, 2000);
  TH1D* h_iRPC_muon2_m = new TH1D("h_iRPC_muon2_m", "TOF mass;2nd leading #tilde{#tau} mass (GeV);Events / 10GeV", 200, 0, 2000);
  TH1D* h_iRPC_muon1_beta = new TH1D("h_iRPC_muon1_beta", "Beta distribution;Leading #tilde{#tau} #beta;Events / 0.02", 50, 0, 1);
  TH1D* h_iRPC_muon2_beta = new TH1D("h_iRPC_muon2_beta", "Beta distribution;2nd leading #tilde{#tau} #beta;Events / 0.02", 50, 0, 1);

  TH1D* h_cRPC_muon1_m = new TH1D("h_cRPC_muon1_m", "TOF mass;Leading #tilde{#tau} mass (GeV);Events / 10GeV", 200, 0, 2000);
  TH1D* h_cRPC_muon2_m = new TH1D("h_cRPC_muon2_m", "TOF mass;2nd leading #tilde{#tau} mass (GeV);Events / 10GeV", 200, 0, 2000);
  TH1D* h_cRPC_muon1_beta = new TH1D("h_cRPC_muon1_beta", "Beta distribution;Leading #tilde{#tau} #beta;Events / 0.02", 50, 0, 1);
  TH1D* h_cRPC_muon2_beta = new TH1D("h_cRPC_muon2_beta", "Beta distribution;2nd leading #tilde{#tau} #beta;Events / 0.02", 50, 0, 1);

  h_muon1_m->GetXaxis()->SetNdivisions(505);
  h_muon2_m->GetXaxis()->SetNdivisions(505);
  h_iRPC_muon1_m->GetXaxis()->SetNdivisions(505);
  h_iRPC_muon2_m->GetXaxis()->SetNdivisions(505);
  h_cRPC_muon1_m->GetXaxis()->SetNdivisions(505);
  h_cRPC_muon2_m->GetXaxis()->SetNdivisions(505);
  h_muon1_m__muon2_m->GetXaxis()->SetNdivisions(505);
  h_muon1_m__muon2_m->GetYaxis()->SetNdivisions(505);

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if ( gen1_pdgId != 0 and gen2_pdgId != 0 ) {

      const int region1 = std::abs(gen1_eta) < 1.8 ? 0 : std::abs(gen1_eta) < 2.4 ? 1 : 2;
      const int region2 = std::abs(gen2_eta) < 1.8 ? 0 : std::abs(gen2_eta) < 2.4 ? 1 : 2;

      h_all_beta1[0]->Fill(gen1_beta);
      if      ( region1 == 0 ) h_cRPC_beta1[0]->Fill(gen1_beta);
      else if ( region1 == 1 ) h_iRPC_beta1[0]->Fill(gen1_beta);
      h_all_beta2[0]->Fill(gen2_beta);
      if      ( region2 == 0 ) h_cRPC_beta2[0]->Fill(gen2_beta);
      else if ( region2 == 1 ) h_iRPC_beta2[0]->Fill(gen2_beta);
      for ( size_t i=1; i<betaRes.size(); ++i ) {
        const double newBeta1 = gRandom->Gaus(gen1_beta, gen1_beta*betaRes[i]/100.);
        h_all_beta1[i]->Fill(newBeta1);
        if      ( region1 == 0 ) h_cRPC_beta1[i]->Fill(newBeta1);
        else if ( region1 == 1 ) h_iRPC_beta1[i]->Fill(newBeta1);

        const double newBeta2 = gRandom->Gaus(gen2_beta, gen2_beta*betaRes[i]/100.);
        h_all_beta2[i]->Fill(newBeta2);
        if      ( region2 == 0 ) h_cRPC_beta2[i]->Fill(newBeta2);
        else if ( region2 == 1 ) h_iRPC_beta2[i]->Fill(newBeta2);
      }
    }

    if ( muon_n > 0 ) {
      const double pt1 = muon_pt[0], eta1 = muon_eta[0];
      const double beta1 = muon_RPCBeta[0];
      const double m1 = beta1 == 0 ? 0 : pt1*cosh(abs(eta1))*sqrt(1./beta1/beta1-1);

      h_muon1_beta->Fill(beta1);
      h_muon1_m->Fill(m1);

      if ( muon_nIRPC[0] > 0 ) {
        h_iRPC_muon1_beta->Fill(beta1);
        h_iRPC_muon1_m->Fill(m1);
      }
      else {
        h_cRPC_muon1_beta->Fill(beta1);
        h_cRPC_muon1_m->Fill(m1);
      }
    }
    if ( muon_n > 1 ) {
      const double pt1 = muon_pt[0], eta1 = muon_eta[0];
      const double beta1 = muon_RPCBeta[0];
      const double m1 = beta1 == 0 ? 0 : pt1*cosh(abs(eta1))*sqrt(1./beta1/beta1-1);
      const double pt2 = muon_pt[1], eta2 = muon_eta[1];
      const double beta2 = muon_RPCBeta[1];
      const double m2 = beta2 == 0 ? 0 : pt2*cosh(abs(eta2))*sqrt(1./beta2/beta2-1);

      h_muon2_beta->Fill(beta2);
      h_muon2_m->Fill(m2);

      h_muon1_beta__muon2_beta->Fill(beta1, beta2);
      h_muon1_m__muon2_m->Fill(m1, m2);

      if ( muon_nIRPC[1] > 0 ) {
        h_iRPC_muon2_beta->Fill(beta2);
        h_iRPC_muon2_m->Fill(m2);
      }
      else {
        h_cRPC_muon2_beta->Fill(beta2);
        h_cRPC_muon2_m->Fill(m2);
      }
    }
  }

  fout->Write();
}


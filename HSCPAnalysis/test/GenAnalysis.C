#define GenAnalysis_cxx
#include "GenAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include <iostream>

void GenAnalysis::Loop(TFile* fout)
{
  //   In a ROOT session, you can do:
  //      root> .L GenAnalysis.C
  //      root> GenAnalysis t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  fout->cd();

  std::vector<TH1D*> h_all_beta1, h_iRPC_beta1, h_cRPC_beta1;
  std::vector<TH1D*> h_all_beta2, h_iRPC_beta2, h_cRPC_beta2;
  std::vector<int> betaRes = {0,1,2,3,4,5,10};
  for ( auto r : betaRes ) {
    TH1D* hh_beta1 = new TH1D(Form("h_beta1_res%03d", r), Form("Beta distribution res=%d%%;#tilde{#tau}^{-} #beta;Events / 0.02", r), 50, 0, 1);
    TH1D* hh_iRPC_beta1 = new TH1D(Form("h_iRPC_beta1_res%03d", r), Form("iRPC Beta distribution res=%d%%;#tilde{#tau}^{-} #beta;Events / 0.02", r), 50, 0, 1);
    TH1D* hh_cRPC_beta1 = new TH1D(Form("h_cRPC_beta1_res%03d", r), Form("cRPC Beta distribution res=%d%%;#tilde{#tau}^{-} #beta;Events / 0.02", r), 50, 0, 1);
    TH1D* hh_beta2 = new TH1D(Form("h_beta2_res%03d", r), Form("Beta distribution res=%d%%;#tilde{#tau}^{+} #beta;Events / 0.02", r), 50, 0, 1);
    TH1D* hh_iRPC_beta2 = new TH1D(Form("h_iRPC_beta2_res%03d", r), Form("iRPC Beta distribution res=%d%%;#tilde{#tau}^{+} #beta;Events / 0.02", r), 50, 0, 1);
    TH1D* hh_cRPC_beta2 = new TH1D(Form("h_cRPC_beta2_res%03d", r), Form("cRPC Beta distribution res=%d%%;#tilde{#tau}^{+} #beta;Events / 0.02", r), 50, 0, 1);
    h_all_beta1.push_back(hh_beta1);
    h_iRPC_beta1.push_back(hh_iRPC_beta1);
    h_cRPC_beta1.push_back(hh_cRPC_beta1);
    h_all_beta2.push_back(hh_beta2);
    h_iRPC_beta2.push_back(hh_iRPC_beta2);
    h_cRPC_beta2.push_back(hh_cRPC_beta2);
  }

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if ( gen1_pdgId == 0 or gen2_pdgId == 0 ) continue;

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

  fout->Write();
}


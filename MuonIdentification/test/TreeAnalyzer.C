#define TreeAnalyzer_cxx
#include "../../HSCPAnalysis/test/TreeAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TEfficiency.h>

#include <iostream>

void TreeAnalyzer::Loop(TFile* fout)
{
  const double minPt = 20, maxAbsEta = 2.5;

  fout->cd();

  auto dir_gen = fout->mkdir("gen"); dir_gen->cd();
  TH1D* h_gen_abseta = new TH1D("h_abseta", "muon |#eta|;Muon |#eta|;Events / 0.05", 50, -maxAbsEta, maxAbsEta);
  TH1D* h_gen_phi = new TH1D("h_phi", "muon #phi;Muon #phi;Events / 0.05", 50, -3.15, 3.15);
  auto dir_all = fout->mkdir("all"); dir_all->cd();
  TH1D* h_all_abseta = new TH1D("h_abseta", "muon |#eta|;Muon |#eta|;Events / 0.05", 50, -maxAbsEta, maxAbsEta);
  TH1D* h_all_phi = new TH1D("h_phi", "muon #phi;Muon #phi;Events / 0.05", 50, -3.15, 3.15);
  TH1D* h_all_minDR = new TH1D("h_minDR", "#Delta{}R(muon, gen);Events / 0.005", 60, 0, 0.1);
  auto dir_loose = fout->mkdir("loose"); dir_loose->cd();
  TH1D* h_loose_abseta = new TH1D("h_abseta", "muon |#eta|;Muon |#eta|;Events / 0.05", 50, -maxAbsEta, maxAbsEta);
  TH1D* h_loose_phi = new TH1D("h_phi", "muon #phi;Muon #phi;Events / 0.05", 50, -3.15, 3.15);
  TH1D* h_loose_minDR = new TH1D("h_minDR", "#Delta{}R(muon, gen);Events / 0.005", 60, 0, 0.1);
  auto dir_tight = fout->mkdir("tight"); dir_tight->cd();
  TH1D* h_tight_abseta = new TH1D("h_abseta", "muon |#eta|;Muon |#eta|;Events / 0.05", 50, -maxAbsEta, maxAbsEta);
  TH1D* h_tight_phi = new TH1D("h_phi", "muon #phi;Muon #phi;Events / 0.05", 50, -3.15, 3.15);
  TH1D* h_tight_minDR = new TH1D("h_minDR", "#Delta{}R(muon, gen);Events / 0.005", 60, 0, 0.1);
  auto dir_RPC = fout->mkdir("RPC"); dir_RPC->cd();
  TH1D* h_RPC_abseta = new TH1D("h_abseta", "muon |#eta|;Muon |#eta|;Events / 0.05", 50, -maxAbsEta, maxAbsEta);
  TH1D* h_RPC_phi = new TH1D("h_phi", "muon #phi;Muon #phi;Events / 0.05", 50, -3.15, 3.15);
  TH1D* h_RPC_minDR = new TH1D("h_minDR", "#Delta{}R(muon, gen);Events / 0.005", 60, 0, 0.1);

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    printf("Processing %lldth event (%.0f%%)\r", jentry, 100.*jentry/nentries);

    std::vector<TLorentzVector> muons_p4;
    for ( unsigned int i=0; i<muon_n; ++i ) {
      if ( muon_pt[i] < minPt or std::abs(muon_eta[i]) > maxAbsEta ) continue;
      TLorentzVector muon_p4;
      muon_p4.SetPtEtaPhiM(muon_pt[i], muon_eta[i], muon_phi[i], 0);
  
      muons_p4.push_back(muon_p4);
    }

    if ( gen1_pdgId != 0 and gen1_pt > minPt and abs(gen1_eta) < maxAbsEta ) {
      const double gen1_abseta = gen1_eta;//abs(gen1_eta);
      h_gen_abseta->Fill(gen1_abseta);
      h_gen_phi->Fill(gen1_phi);

      TLorentzVector gen_p4;
      gen_p4.SetPtEtaPhiM(gen1_pt, gen1_eta, gen1_phi, 0);

      int matchedIdx = -1;
      double minDR = 0.1;
      for ( unsigned int i=0; i<muons_p4.size(); ++i ) {
        const double dR = gen_p4.DeltaR(muons_p4[i]);
        if ( dR < minDR ) {
          minDR = dR;
          matchedIdx = i;
        }
      }

      if ( matchedIdx >= 0 ) {
        const auto muon_p4 = muons_p4[matchedIdx];
        const double dR = gen_p4.DeltaR(muons_p4[matchedIdx]);
        h_all_abseta->Fill(gen1_abseta);
        h_all_phi->Fill(gen1_phi);
        h_all_minDR->Fill(dR);
        if ( muon_isLoose[matchedIdx] ) {
          h_loose_abseta->Fill(gen1_abseta);
          h_loose_phi->Fill(gen1_phi);
          h_loose_minDR->Fill(dR);
        }
        if ( muon_isTight[matchedIdx] ) {
          h_tight_abseta->Fill(gen1_abseta);
          h_tight_phi->Fill(gen1_phi);
          h_tight_minDR->Fill(dR);
        }
        if ( muon_isRPC[matchedIdx] ) {
          h_RPC_abseta->Fill(gen1_abseta);
          h_RPC_phi->Fill(gen1_phi);
          h_RPC_minDR->Fill(dR);
        }

      }
    }

    if ( gen2_pdgId != 0 and gen2_pt > minPt and abs(gen2_eta) < maxAbsEta ) {
      const double gen2_abseta = gen2_eta; // abs(gen2_eta);
      h_gen_abseta->Fill(gen2_abseta);
      h_gen_phi->Fill(gen2_phi);

      TLorentzVector gen_p4;
      gen_p4.SetPtEtaPhiM(gen2_pt, gen2_eta, gen2_phi, 0);

      int matchedIdx = -1;
      double minDR = 0.1;
      for ( unsigned int i=0; i<muons_p4.size(); ++i ) {
        const double dR = gen_p4.DeltaR(muons_p4[i]);
        if ( dR < minDR ) {
          minDR = dR;
          matchedIdx = i;
        }
      }

      if ( matchedIdx >= 0 ) {
        const auto muon_p4 = muons_p4[matchedIdx];
        const double dR = gen_p4.DeltaR(muons_p4[matchedIdx]);
        h_all_abseta->Fill(gen2_abseta);
        h_all_phi->Fill(gen2_phi);
        h_all_minDR->Fill(dR);
        if ( muon_isLoose[matchedIdx] ) {
          h_loose_abseta->Fill(gen2_abseta);
          h_loose_phi->Fill(gen2_phi);
          h_loose_minDR->Fill(dR);
        }
        if ( muon_isTight[matchedIdx] ) {
          h_tight_abseta->Fill(gen2_abseta);
          h_tight_phi->Fill(gen2_phi);
          h_tight_minDR->Fill(dR);
        }
        if ( muon_isTight[matchedIdx] ) {
          h_RPC_abseta->Fill(gen2_abseta);
          h_RPC_phi->Fill(gen2_phi);
          h_RPC_minDR->Fill(dR);
        }
      }
    }
  }

  dir_all->cd();
  TEfficiency eff_all_abseta(*h_all_abseta, *h_gen_abseta);
  eff_all_abseta.SetName("eff_abseta"); eff_all_abseta.Write();
  TEfficiency eff_all_phi(*h_all_phi, *h_gen_phi);
  eff_all_phi.SetName("eff_phi"); eff_all_phi.Write();
  dir_loose->cd();
  TEfficiency eff_loose_abseta(*h_loose_abseta, *h_gen_abseta);
  eff_loose_abseta.SetName("eff_abseta"); eff_loose_abseta.Write();
  TEfficiency eff_loose_phi(*h_loose_phi, *h_gen_phi);
  eff_loose_phi.SetName("eff_phi"); eff_loose_phi.Write();
  dir_tight->cd();
  TEfficiency eff_tight_abseta(*h_tight_abseta, *h_gen_abseta);
  eff_tight_abseta.SetName("eff_abseta"); eff_tight_abseta.Write();
  TEfficiency eff_tight_phi(*h_tight_phi, *h_gen_phi);
  eff_tight_phi.SetName("eff_phi"); eff_tight_phi.Write();
  dir_RPC->cd();
  TEfficiency eff_RPC_abseta(*h_RPC_abseta, *h_gen_abseta);
  eff_RPC_abseta.SetName("eff_abseta"); eff_RPC_abseta.Write();
  TEfficiency eff_RPC_phi(*h_RPC_phi, *h_gen_phi);
  eff_RPC_phi.SetName("eff_phi"); eff_RPC_phi.Write();

  fout->Write();
}


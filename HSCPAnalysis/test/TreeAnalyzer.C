#define TreeAnalyzer_cxx
#include "TreeAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>

#include <vector>
#include <iostream>
#include <algorithm>

const double muonMass = 0.1057;
const double speedOfLight = 29.979; // 30cm/ns
const int bxLo = -8, nBx = 17;
const unsigned minNhitCluster = 2;

double deltaPhi(const double phi1, const double phi2)
{
  double dphi = phi1 - phi2;
  if ( dphi < -TMath::TwoPi() ) dphi += TMath::TwoPi();
  else if ( dphi > TMath::TwoPi() ) dphi -= TMath::TwoPi();
  return dphi;
}

struct HEtaPhiBeta
{
  constexpr static int nbinsEta = 25, nbinsPhi = 50, nbinsBeta = 15;
  constexpr static double minEta = -2.5, maxEta = 2.5;
  constexpr static double minPhi = 0, maxPhi = 2*3.14159265358979323846L+1e-9;
  constexpr static double minBeta = -1, maxBeta = 2;
  constexpr static double binwEta  = (maxEta-minEta)/nbinsEta;
  constexpr static double binwPhi  = (maxPhi-minPhi)/nbinsPhi;
  constexpr static double binwBeta = (maxBeta-minBeta)/nbinsBeta;

  unsigned findBin(const double eta, const double phi, const double beta) const
  {
    const int ieta  = std::max(-1, std::min(nbinsEta, int(floor((eta-minEta)/binwEta))));
    const int iphi  = std::max(-1, std::min(nbinsPhi, int(floor((phi-minPhi)/binwPhi))));
    const int ibeta = std::max(-1, std::min(nbinsBeta, int(floor((beta-minBeta)/binwBeta))));

    const unsigned ibin = (ieta+1)
                        + (iphi+1)*(nbinsEta+2);
                        //+ (ibeta+1)*(nbinsEta+2)*(nbinsPhi+2);
    return ibin;
  }

  std::vector<double> findBinLowEdge(unsigned ibin) const
  {
    const int ieta = ibin % (nbinsEta+2) - 1;
    ibin /= nbinsEta+2;
    const int iphi = ibin % (nbinsPhi+2) - 1;
    //ibin /= nbinsPhi+2;
    //const int ibeta = ibin - 1;

    const double eta  = ieta*binwEta+minEta;
    const double phi  = iphi*binwPhi+minPhi;
    //const double beta = ibeta*binwBeta+minBeta;

    std::vector<double> res = {eta, phi};//, beta};
    return res;
  }

  void fill(const double eta, const double phi, const double beta, const unsigned idx)
  {
    const unsigned ibin = findBin(eta, phi, beta);
    auto itr = contents.find(ibin);
    if ( itr == contents.end() ) itr = contents.insert(std::make_pair(ibin, std::vector<unsigned>())).first;
    itr->second.push_back(idx);
  }

  std::map<unsigned, std::vector<unsigned>> contents;
};

std::vector<std::vector<unsigned>> TreeAnalyzer::clusterHitsByGenP4s(const TLorentzVector p4s[]) const
{
  // Cluster hits along the GenParticle's four momentum, just for initial testing.
  std::vector<std::vector<unsigned>> clusters;
  clusters.resize(2);

  for ( unsigned i=0; i<rpcHit_n; ++i ) {
    const TVector3 pos(rpcHit_x[i], rpcHit_y[i], rpcHit_z[i]);
    double minDR = 0.3;
    int match = -1;
    for ( unsigned j=0; j<2; ++j ) {
      const double dR = p4s[j].Vect().DeltaR(pos);
      if ( dR < minDR ) {
        minDR = dR;
        match = j;
      }
    }
    if ( match >= 0 ) clusters.at(match).push_back(i);
  }

  return clusters;
}

std::vector<std::vector<unsigned>> TreeAnalyzer::clusterHitsByEtaPhi() const
{
  // Hit clustering in 3D, eta-phi-beta space with vx=vy=vz=0 & bx=0 hypothesis
  std::vector<std::vector<unsigned>> clusters;
  // First step to fill "histograms"
  HEtaPhiBeta h;
  for ( unsigned i=0; i<rpcHit_n; ++i ) {
    const TVector3 pos(rpcHit_x[i], rpcHit_y[i], rpcHit_z[i]);
    const double eta = pos.Eta(), phi = pos.Phi();
    const double ct = speedOfLight*rpcHit_time[i];
    const double beta = 1./(1+ct/pos.Mag());

    h.fill(eta, phi, beta, i);
  }

  for ( auto item : h.contents ) {
    if ( item.second.size() < minNhitCluster ) continue;
    clusters.push_back(item.second);
  }

  return clusters;
}

std::vector<double> TreeAnalyzer::fitTrackBxConstrained(const std::vector<unsigned>& hits) const
{
  // result: qual, beta, t0
  std::vector<double> result = {1e9, 0, 0, 0, 0};
  const unsigned n = hits.size();
  if ( n == 0 ) return result;

  double sumBeta = 0;
  double sumTimeDiff2 = 0;
  double firstPhi = 0; // this is necessary to shift all phi's to the new origin. without shifting, +2pi and -2pi will result unphysical large variance
  double sumEta = 0, sumEta2 = 0, sumDphi = 0, sumDphi2 = 0;

  unsigned nValidBeta = 0;
  for ( auto i : hits ) {
    const TVector3 pos(rpcHit_x[i], rpcHit_y[i], rpcHit_z[i]);
    const double ct = speedOfLight*rpcHit_time[i];
    const double ibeta = 1./(1+ct/pos.Mag());
    if ( ibeta > 2 or ibeta < -1 ) continue;
    sumBeta += ibeta;
    sumTimeDiff2 += rpcHit_time[i]*rpcHit_time[i];

    sumEta += pos.Eta();
    sumEta2 += pos.Eta()*pos.Eta();
    if ( nValidBeta == 0 ) firstPhi = pos.Phi();
    else {
      const double dphi = pos.Phi()-firstPhi;
      sumDphi += dphi;
      sumDphi2 += dphi*dphi;
    }

    ++nValidBeta;
  }
  const double dRErr2 = ((sumDphi2-sumDphi*sumDphi/nValidBeta) +
                         (sumEta2-sumEta*sumEta/nValidBeta))/nValidBeta;

  result = {dRErr2, sumBeta/nValidBeta, sumTimeDiff2/nValidBeta};

  return result;
}

std::vector<double> TreeAnalyzer::fitTrackSlope(const std::vector<unsigned>& hits) const
{
  std::vector<double> result = {1e9, 0, 0, 0, 0};
  const unsigned n = hits.size();
  if ( n <= 2 ) return result;

  double sx = 0, sy = 0, sxy = 0, sxx = 0, syy = 0;
  for ( auto i : hits ) {
    const double r2 = rpcHit_x[i]*rpcHit_x[i] + rpcHit_y[i]*rpcHit_y[i] + rpcHit_z[i]*rpcHit_z[i];
    const double r = std::sqrt(r2);
    const double t = rpcHit_time[i];
    sx  += r;
    sy  += t;
    sxy += r*t;
    sxx += r2;
    syy += t*t;
  }
  const double ssxy = sxy-sx*sy/n;
  const double ssxx = sxx-sx*sx/n;
  const double ssyy = syy-sy*sy/n;
  const double b = ssxy/ssxx;
  const double s = std::sqrt((ssyy - b*ssxy)/(n-2));
  const double bStdErr = s/std::sqrt(ssxx);

  const double t0 = (sy-b*sx)/n; // is "a" in the original code
  const double t0Err = s*std::sqrt(1./n + sx*sx*ssxx/n/n); // is "aStdErr" in the original code
  const double beta = 1./(b*speedOfLight+1.);
  const double betaErr = speedOfLight*bStdErr/((b*speedOfLight+1)*(b*speedOfLight+1));

  const int nbx = std::round(t0/25.);
  const double bxPull = t0 - nbx*25;
  //const double bxPull = t0/25.;
  result = {bxPull, beta, betaErr, t0, t0Err};

  return result;
}

void TreeAnalyzer::Loop(TFile* fout)
{
  fout->cd();
  TTree* tree = new TTree("tree", "tree");

  TLorentzVector out_gens_p4[2];
  int out_gens_pdgId[2];
  tree->Branch("gen1_p4", "TLorentzVector", &out_gens_p4[0]);
  tree->Branch("gen2_p4", "TLorentzVector", &out_gens_p4[1]);
  tree->Branch("gen1_pdgId", &out_gens_pdgId[0], "gen1_pdgId/I");
  tree->Branch("gen2_pdgId", &out_gens_pdgId[1], "gen2_pdgId/I");

  unsigned out_muons_n;
  TLorentzVector out_muons_p4[3];
  int out_muons_q[3];
  tree->Branch("muon1_p4", "TLorentzVector", &out_muons_p4[0]);
  tree->Branch("muon2_p4", "TLorentzVector", &out_muons_p4[1]);
  tree->Branch("muon3_p4", "TLorentzVector", &out_muons_p4[2]);
  tree->Branch("muon1_q", &out_muons_q[0], "muon1_q/I");
  tree->Branch("muon2_q", &out_muons_q[1], "muon2_q/I");
  tree->Branch("muon3_q", &out_muons_q[2], "muon3_q/I");

  float out_fit_quals[2];
  float out_fit_betas[2];
  tree->Branch("fit_qual1", &out_fit_quals[0], "fit_qual1/F");
  tree->Branch("fit_qual2", &out_fit_quals[1], "fit_qual2/F");
  tree->Branch("fit_beta1", &out_fit_betas[0], "fit_beta1/F");
  tree->Branch("fit_beta2", &out_fit_betas[1], "fit_beta2/F");

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( (10000*(jentry+1)/nentries) % 100 == 0 ) {
      printf("Processing %lld/%lld, %.0f %% done...\r",  jentry+1, nentries, 100.*jentry/nentries);
    }

    // Initialize variables
    for ( unsigned i=0; i<2; ++i ) {
      out_gens_p4[i].SetXYZT(0,0,0,0);
      out_gens_pdgId[i] = 0;
    }
    out_muons_n = 0;
    for ( unsigned i=0; i<3; ++i ) {
      out_muons_p4[i].SetXYZT(0,0,0,0);
      out_muons_q[i] = 0;
    }
    for ( unsigned i=0; i<2; ++i ) {
      out_fit_quals[i] = 1e9;
      out_fit_betas[i] = 0;
    }

    // Fill particles
    if ( gen1_pdgId != 0 ) {
      out_gens_p4[0].SetPtEtaPhiM(gen1_pt, gen1_eta, gen1_phi, gen1_m);
      out_gens_pdgId[0] = gen1_pdgId;
    }
    if ( gen2_pdgId != 0 ) {
      out_gens_p4[1].SetPtEtaPhiM(gen2_pt, gen2_eta, gen2_phi, gen2_m);
      out_gens_pdgId[1] = gen2_pdgId;
    }

    std::vector<unsigned> muonIdxs;
    for ( unsigned i=0; i<muon_n; ++i ) {
      if ( muon_pt[i] < 30 or std::abs(muon_eta[i]) > 2.4 ) continue;
      muonIdxs.push_back(i);
    }
    std::sort(muonIdxs.begin(), muonIdxs.end(), [&](unsigned i, unsigned j){return muon_pt[i] > muon_pt[j];});
    for ( unsigned i=0, n=std::min(3ul, muonIdxs.size()); i<n; ++i ) {
      const auto ii = muonIdxs[i];
      out_muons_p4[i].SetPtEtaPhiM(muon_pt[ii], muon_eta[ii], muon_phi[ii], muonMass);
      out_muons_q[i] = muon_q[i];
    }

    // Cluster hits and do the fitting
    //const auto hitClusters = clusterHitsByGenP4s(out_gens_p4);
    auto hitClusters = clusterHitsByEtaPhi();
    std::sort(hitClusters.begin(), hitClusters.end(),
              [](const std::vector<unsigned>& a, const std::vector<unsigned>& b){return a.size() > b.size();});
    for ( unsigned i=0, n=std::min(2ul, hitClusters.size()); i<n; ++i ) {
      //const auto res = fitTrackBxConstrained(hitClusters[i]);
      const auto res = fitTrackSlope(hitClusters[i]);
      out_fit_quals[i] = res[0];
      out_fit_betas[i] = res[1];
    }

/*
      const double ct = speedOfLight*rpcHit_time[i];
      for ( int bx = bxLo; bx <= bxLo+nBx; ++bx ) {
        const double ibeta = 1./(1+speedOfLight*(rpcHit_time[i]-bx*25)/r);
        hitsByBeta[bx-bxLo][ibeta] = i;
      }
*/

    if ( out_fit_quals[0] >= 1e9 or out_fit_quals[1] >= 1e9 ) continue;

    tree->Fill();
  }
  cout << "Processing " << nentries << "/" << nentries << "\n"; // Just to print last event

  fout->Write();
}


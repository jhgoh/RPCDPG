//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct  9 02:04:nHit 2018 by ROOT version 6.10/05
// from TTree tree/tree
// found on file: ntuple/DYJetsToLL_M-50_noPU.root
//////////////////////////////////////////////////////////

#ifndef TreeAnalyzer_h
#define TreeAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>

// Header file for the classes stored in the TTree if any.

class TreeAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   const static unsigned int nHit = 1000;
   const static unsigned int nMuon = 100;

   // Declaration of leaf types
   Int_t           gen1_pdgId;
   Double_t        gen1_pt;
   Double_t        gen1_eta;
   Double_t        gen1_phi;
   Double_t        gen1_m;
   Double_t        gen1_beta;
   Int_t           gen2_pdgId;
   Double_t        gen2_pt;
   Double_t        gen2_eta;
   Double_t        gen2_phi;
   Double_t        gen2_m;
   Double_t        gen2_beta;
   UShort_t        rpcHit_n;
   Bool_t          rpcHit_isBarrel[nHit];   //[rpcHit_n]
   Bool_t          rpcHit_isIRPC[nHit];   //[rpcHit_n]
   Short_t         rpcHit_sector[nHit];   //[rpcHit_n]
   Short_t         rpcHit_station[nHit];   //[rpcHit_n]
   Short_t         rpcHit_wheel[nHit];   //[rpcHit_n]
   Short_t         rpcHit_layer[nHit];   //[rpcHit_n]
   Short_t         rpcHit_disk[nHit];   //[rpcHit_n]
   Short_t         rpcHit_ring[nHit];   //[rpcHit_n]
   Double_t        rpcHit_x[nHit];   //[rpcHit_n]
   Double_t        rpcHit_y[nHit];   //[rpcHit_n]
   Double_t        rpcHit_z[nHit];   //[rpcHit_n]
   Double_t        rpcHit_time[nHit];   //[rpcHit_n]
   Double_t        rpcHit_timeErr[nHit];   //[rpcHit_n]
   Double_t        rpcHit_lx[nHit];   //[rpcHit_n]
   Double_t        rpcHit_ly[nHit];   //[rpcHit_n]
   Short_t         rpcHit_bx[nHit];   //[rpcHit_n]
   UShort_t        dtSegment_n;
   Double_t        dtSegment_x[nHit];   //[dtSegment_n]
   Double_t        dtSegment_y[nHit];   //[dtSegment_n]
   Double_t        dtSegment_z[nHit];   //[dtSegment_n]
   Double_t        dtSegment_dx[nHit];   //[dtSegment_n]
   Double_t        dtSegment_dy[nHit];   //[dtSegment_n]
   Double_t        dtSegment_dz[nHit];   //[dtSegment_n]
   Double_t        dtSegment_time[nHit];   //[dtSegment_n]
   Double_t        dtSegment_lx[nHit];   //[dtSegment_n]
   Double_t        dtSegment_ly[nHit];   //[dtSegment_n]
   UShort_t        cscSegment_n;
   Double_t        cscSegment_x[nHit];   //[cscSegment_n]
   Double_t        cscSegment_y[nHit];   //[cscSegment_n]
   Double_t        cscSegment_z[nHit];   //[cscSegment_n]
   Double_t        cscSegment_dx[nHit];   //[cscSegment_n]
   Double_t        cscSegment_dy[nHit];   //[cscSegment_n]
   Double_t        cscSegment_dz[nHit];   //[cscSegment_n]
   Double_t        cscSegment_time[nHit];   //[cscSegment_n]
   Double_t        cscSegment_lx[nHit];   //[cscSegment_n]
   Double_t        cscSegment_ly[nHit];   //[cscSegment_n]
   UShort_t        gemSegment_n;
   Double_t        gemSegment_x[nHit];   //[gemSegment_n]
   Double_t        gemSegment_y[nHit];   //[gemSegment_n]
   Double_t        gemSegment_z[nHit];   //[gemSegment_n]
   Double_t        gemSegment_dx[nHit];   //[gemSegment_n]
   Double_t        gemSegment_dy[nHit];   //[gemSegment_n]
   Double_t        gemSegment_dz[nHit];   //[gemSegment_n]
   Double_t        gemSegment_time[nHit];   //[gemSegment_n]
   Double_t        gemSegment_lx[nHit];   //[gemSegment_n]
   Double_t        gemSegment_ly[nHit];   //[gemSegment_n]
   UShort_t        muon_n;
   Double_t        muon_pt[nMuon];   //[muon_n]
   Double_t        muon_eta[nMuon];   //[muon_n]
   Double_t        muon_phi[nMuon];   //[muon_n]
   Short_t         muon_q[nMuon];   //[muon_n]
   Bool_t          muon_isLoose[nMuon];   //[muon_n]
   Bool_t          muon_isTight[nMuon];   //[muon_n]
   Bool_t          muon_isRPC[nMuon];   //[muon_n]
   Double_t        muon_time[nMuon];   //[muon_n]
   Double_t        muon_RPCTime[nMuon];   //[muon_n]
   Double_t        muon_RPCTimeNew[nMuon];   //[muon_n]
   Double_t        muon_RPCBeta[nMuon];   //[muon_n]
   UShort_t        muon_nRPC[nMuon];   //[muon_n]
   UShort_t        muon_nIRPC[nMuon];   //[muon_n]
   Double_t        muon_genDR[nMuon];   //[muon_n]
   Short_t         muon_genPdgId[nMuon];   //[muon_n]

   // List of branches
   TBranch        *b_gen1_pdgId;   //!
   TBranch        *b_gen1_pt;   //!
   TBranch        *b_gen1_eta;   //!
   TBranch        *b_gen1_phi;   //!
   TBranch        *b_gen1_m;   //!
   TBranch        *b_gen1_beta;   //!
   TBranch        *b_gen2_pdgId;   //!
   TBranch        *b_gen2_pt;   //!
   TBranch        *b_gen2_eta;   //!
   TBranch        *b_gen2_phi;   //!
   TBranch        *b_gen2_m;   //!
   TBranch        *b_gen2_beta;   //!
   TBranch        *b_rpcHit_n;   //!
   TBranch        *b_rpcHit_isBarrel;   //!
   TBranch        *b_rpcHit_isIRPC;   //!
   TBranch        *b_rpcHit_sector;   //!
   TBranch        *b_rpcHit_station;   //!
   TBranch        *b_rpcHit_wheel;   //!
   TBranch        *b_rpcHit_layer;   //!
   TBranch        *b_rpcHit_disk;   //!
   TBranch        *b_rpcHit_ring;   //!
   TBranch        *b_rpcHit_x;   //!
   TBranch        *b_rpcHit_y;   //!
   TBranch        *b_rpcHit_z;   //!
   TBranch        *b_rpcHit_time;   //!
   TBranch        *b_rpcHit_timeErr;   //!
   TBranch        *b_rpcHit_lx;   //!
   TBranch        *b_rpcHit_ly;   //!
   TBranch        *b_rpcHit_bx;   //!
   TBranch        *b_dtSegment_n;   //!
   TBranch        *b_dtSegment_x;   //!
   TBranch        *b_dtSegment_y;   //!
   TBranch        *b_dtSegment_z;   //!
   TBranch        *b_dtSegment_dx;   //!
   TBranch        *b_dtSegment_dy;   //!
   TBranch        *b_dtSegment_dz;   //!
   TBranch        *b_dtSegment_time;   //!
   TBranch        *b_dtSegment_lx;   //!
   TBranch        *b_dtSegment_ly;   //!
   TBranch        *b_cscSegment_n;   //!
   TBranch        *b_cscSegment_x;   //!
   TBranch        *b_cscSegment_y;   //!
   TBranch        *b_cscSegment_z;   //!
   TBranch        *b_cscSegment_dx;   //!
   TBranch        *b_cscSegment_dy;   //!
   TBranch        *b_cscSegment_dz;   //!
   TBranch        *b_cscSegment_time;   //!
   TBranch        *b_cscSegment_lx;   //!
   TBranch        *b_cscSegment_ly;   //!
   TBranch        *b_gemSegment_n;   //!
   TBranch        *b_gemSegment_x;   //!
   TBranch        *b_gemSegment_y;   //!
   TBranch        *b_gemSegment_z;   //!
   TBranch        *b_gemSegment_dx;   //!
   TBranch        *b_gemSegment_dy;   //!
   TBranch        *b_gemSegment_dz;   //!
   TBranch        *b_gemSegment_time;   //!
   TBranch        *b_gemSegment_lx;   //!
   TBranch        *b_gemSegment_ly;   //!
   TBranch        *b_muon_n;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_q;   //!
   TBranch        *b_muon_isLoose;   //!
   TBranch        *b_muon_isTight;   //!
   TBranch        *b_muon_isRPC;   //!
   TBranch        *b_muon_time;   //!
   TBranch        *b_muon_RPCTime;   //!
   TBranch        *b_muon_RPCTimeNew;   //!
   TBranch        *b_muon_RPCBeta;   //!
   TBranch        *b_muon_nRPC;   //!
   TBranch        *b_muon_nIRPC;   //!
   TBranch        *b_muon_genDR;   //!
   TBranch        *b_muon_genPdgId;   //!

   TreeAnalyzer(TTree *tree=0);
   virtual ~TreeAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TFile* fout);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   std::vector<std::vector<unsigned>> clusterHitsByGenP4s(const TLorentzVector p4s[]) const;
   std::vector<double> fitTrackBxConstrained(const std::vector<unsigned>& hits) const;
   std::vector<double> fitTrackSlope(const std::vector<unsigned>& hits) const;


};

#endif

#ifdef TreeAnalyzer_cxx
TreeAnalyzer::TreeAnalyzer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ntuple/DYJetsToLL_M-50_noPU.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ntuple/DYJetsToLL_M-50_noPU.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ntuple/DYJetsToLL_M-50_noPU.root:/HSCPTree");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

TreeAnalyzer::~TreeAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("gen1_pdgId", &gen1_pdgId, &b_gen1_pdgId);
   fChain->SetBranchAddress("gen1_pt", &gen1_pt, &b_gen1_pt);
   fChain->SetBranchAddress("gen1_eta", &gen1_eta, &b_gen1_eta);
   fChain->SetBranchAddress("gen1_phi", &gen1_phi, &b_gen1_phi);
   fChain->SetBranchAddress("gen1_m", &gen1_m, &b_gen1_m);
   fChain->SetBranchAddress("gen1_beta", &gen1_beta, &b_gen1_beta);
   fChain->SetBranchAddress("gen2_pdgId", &gen2_pdgId, &b_gen2_pdgId);
   fChain->SetBranchAddress("gen2_pt", &gen2_pt, &b_gen2_pt);
   fChain->SetBranchAddress("gen2_eta", &gen2_eta, &b_gen2_eta);
   fChain->SetBranchAddress("gen2_phi", &gen2_phi, &b_gen2_phi);
   fChain->SetBranchAddress("gen2_m", &gen2_m, &b_gen2_m);
   fChain->SetBranchAddress("gen2_beta", &gen2_beta, &b_gen2_beta);
   fChain->SetBranchAddress("rpcHit_n", &rpcHit_n, &b_rpcHit_n);
   fChain->SetBranchAddress("rpcHit_isBarrel", rpcHit_isBarrel, &b_rpcHit_isBarrel);
   fChain->SetBranchAddress("rpcHit_isIRPC", rpcHit_isIRPC, &b_rpcHit_isIRPC);
   fChain->SetBranchAddress("rpcHit_sector", rpcHit_sector, &b_rpcHit_sector);
   fChain->SetBranchAddress("rpcHit_station", rpcHit_station, &b_rpcHit_station);
   fChain->SetBranchAddress("rpcHit_wheel", rpcHit_wheel, &b_rpcHit_wheel);
   fChain->SetBranchAddress("rpcHit_layer", rpcHit_layer, &b_rpcHit_layer);
   fChain->SetBranchAddress("rpcHit_disk", rpcHit_disk, &b_rpcHit_disk);
   fChain->SetBranchAddress("rpcHit_ring", rpcHit_ring, &b_rpcHit_ring);
   fChain->SetBranchAddress("rpcHit_x", rpcHit_x, &b_rpcHit_x);
   fChain->SetBranchAddress("rpcHit_y", rpcHit_y, &b_rpcHit_y);
   fChain->SetBranchAddress("rpcHit_z", rpcHit_z, &b_rpcHit_z);
   fChain->SetBranchAddress("rpcHit_time", rpcHit_time, &b_rpcHit_time);
   fChain->SetBranchAddress("rpcHit_timeErr", rpcHit_timeErr, &b_rpcHit_timeErr);
   fChain->SetBranchAddress("rpcHit_lx", rpcHit_lx, &b_rpcHit_lx);
   fChain->SetBranchAddress("rpcHit_ly", rpcHit_ly, &b_rpcHit_ly);
   fChain->SetBranchAddress("rpcHit_bx", rpcHit_bx, &b_rpcHit_bx);
   fChain->SetBranchAddress("dtSegment_n", &dtSegment_n, &b_dtSegment_n);
   fChain->SetBranchAddress("dtSegment_x", dtSegment_x, &b_dtSegment_x);
   fChain->SetBranchAddress("dtSegment_y", dtSegment_y, &b_dtSegment_y);
   fChain->SetBranchAddress("dtSegment_z", dtSegment_z, &b_dtSegment_z);
   fChain->SetBranchAddress("dtSegment_dx", dtSegment_dx, &b_dtSegment_dx);
   fChain->SetBranchAddress("dtSegment_dy", dtSegment_dy, &b_dtSegment_dy);
   fChain->SetBranchAddress("dtSegment_dz", dtSegment_dz, &b_dtSegment_dz);
   fChain->SetBranchAddress("dtSegment_time", dtSegment_time, &b_dtSegment_time);
   fChain->SetBranchAddress("dtSegment_lx", dtSegment_lx, &b_dtSegment_lx);
   fChain->SetBranchAddress("dtSegment_ly", dtSegment_ly, &b_dtSegment_ly);
   fChain->SetBranchAddress("cscSegment_n", &cscSegment_n, &b_cscSegment_n);
   fChain->SetBranchAddress("cscSegment_x", cscSegment_x, &b_cscSegment_x);
   fChain->SetBranchAddress("cscSegment_y", cscSegment_y, &b_cscSegment_y);
   fChain->SetBranchAddress("cscSegment_z", cscSegment_z, &b_cscSegment_z);
   fChain->SetBranchAddress("cscSegment_dx", cscSegment_dx, &b_cscSegment_dx);
   fChain->SetBranchAddress("cscSegment_dy", cscSegment_dy, &b_cscSegment_dy);
   fChain->SetBranchAddress("cscSegment_dz", cscSegment_dz, &b_cscSegment_dz);
   fChain->SetBranchAddress("cscSegment_time", cscSegment_time, &b_cscSegment_time);
   fChain->SetBranchAddress("cscSegment_lx", cscSegment_lx, &b_cscSegment_lx);
   fChain->SetBranchAddress("cscSegment_ly", cscSegment_ly, &b_cscSegment_ly);
   fChain->SetBranchAddress("gemSegment_n", &gemSegment_n, &b_gemSegment_n);
   fChain->SetBranchAddress("gemSegment_x", gemSegment_x, &b_gemSegment_x);
   fChain->SetBranchAddress("gemSegment_y", gemSegment_y, &b_gemSegment_y);
   fChain->SetBranchAddress("gemSegment_z", gemSegment_z, &b_gemSegment_z);
   fChain->SetBranchAddress("gemSegment_dx", gemSegment_dx, &b_gemSegment_dx);
   fChain->SetBranchAddress("gemSegment_dy", gemSegment_dy, &b_gemSegment_dy);
   fChain->SetBranchAddress("gemSegment_dz", gemSegment_dz, &b_gemSegment_dz);
   fChain->SetBranchAddress("gemSegment_time", gemSegment_time, &b_gemSegment_time);
   fChain->SetBranchAddress("gemSegment_lx", gemSegment_lx, &b_gemSegment_lx);
   fChain->SetBranchAddress("gemSegment_ly", gemSegment_ly, &b_gemSegment_ly);
   fChain->SetBranchAddress("muon_n", &muon_n, &b_muon_n);
   fChain->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_q", muon_q, &b_muon_q);
   fChain->SetBranchAddress("muon_isLoose", muon_isLoose, &b_muon_isLoose);
   fChain->SetBranchAddress("muon_isTight", muon_isTight, &b_muon_isTight);
   fChain->SetBranchAddress("muon_isRPC", muon_isRPC, &b_muon_isRPC);
   fChain->SetBranchAddress("muon_time", muon_time, &b_muon_time);
   fChain->SetBranchAddress("muon_RPCTime", muon_RPCTime, &b_muon_RPCTime);
   fChain->SetBranchAddress("muon_RPCTimeNew", muon_RPCTimeNew, &b_muon_RPCTimeNew);
   fChain->SetBranchAddress("muon_RPCBeta", muon_RPCBeta, &b_muon_RPCBeta);
   fChain->SetBranchAddress("muon_nRPC", muon_nRPC, &b_muon_nRPC);
   fChain->SetBranchAddress("muon_nIRPC", muon_nIRPC, &b_muon_nIRPC);
   fChain->SetBranchAddress("muon_genDR", muon_genDR, &b_muon_genDR);
   fChain->SetBranchAddress("muon_genPdgId", muon_genPdgId, &b_muon_genPdgId);
   Notify();
}

Bool_t TreeAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeAnalyzer_cxx

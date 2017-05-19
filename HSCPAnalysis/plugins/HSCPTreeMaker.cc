#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//#include "SimDataFormats/Vertex/interface/SimVertex.h"
//#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
//#include "DataFormats/RPCDigi/interface/RPCDigi.h"
//#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLinkfwd.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class HSCPTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  HSCPTreeMaker(const edm::ParameterSet& pset);
  ~HSCPTreeMaker() {};

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;
  std::vector<std::pair<double, int>> getRPCTimes(const reco::Muon& mu,
                                                  const RPCRecHitCollection& rpcHits, edm::ESHandle<RPCGeometry>& rpcGeom,
                                                  const bool doOnlyWithTime) const;

private:
  const double speedOfLight_;
  const double signalSpeed_;
  const int signalPdgId_;

  const bool doSimHit_;
  const bool doSimDigi_;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  //edm::EDGetTokenT<edm::SimVertexContainer> simVertexToken_;
  edm::EDGetTokenT<edm::PSimHitContainer> rpcSimHitsToken_;
  //edm::EDGetTokenT<RPCDigiCollection> rpcDigisToken_;
  edm::EDGetTokenT<edm::DetSetVector<RPCDigiSimLink>> rpcSimDigisToken_;

  edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

  TTree* tree_;

  int b_gen1_pdgId;
  double b_gen1_pt, b_gen1_eta, b_gen1_phi, b_gen1_m, b_gen1_beta;
  int b_gen2_pdgId;
  double b_gen2_pt, b_gen2_eta, b_gen2_phi, b_gen2_m, b_gen2_beta;

  const static unsigned short simHit1_N = 1000;
  unsigned short b_simHit1_n;
  bool b_simHit1_isBarrel[simHit1_N], b_simHit1_isIRPC[simHit1_N];
  short b_simHit1_sector[simHit1_N], b_simHit1_station[simHit1_N], b_simHit1_roll[simHit1_N];
  short b_simHit1_wheel[simHit1_N], b_simHit1_layer[simHit1_N];
  short b_simHit1_disk[simHit1_N], b_simHit1_ring[simHit1_N];
  double b_simHit1_x[simHit1_N], b_simHit1_y[simHit1_N], b_simHit1_z[simHit1_N];
  double b_simHit1_tof[simHit1_N];
  double b_simHit1_lx[simHit1_N], b_simHit1_ly[simHit1_N];

  const static unsigned short simHit2_N = 1000;
  unsigned short b_simHit2_n;
  bool b_simHit2_isBarrel[simHit2_N], b_simHit2_isIRPC[simHit2_N];
  short b_simHit2_sector[simHit2_N], b_simHit2_station[simHit2_N], b_simHit2_roll[simHit2_N];
  short b_simHit2_wheel[simHit2_N], b_simHit2_layer[simHit2_N];
  short b_simHit2_disk[simHit2_N], b_simHit2_ring[simHit2_N];
  double b_simHit2_x[simHit2_N], b_simHit2_y[simHit2_N], b_simHit2_z[simHit2_N];
  double b_simHit2_tof[simHit2_N];
  double b_simHit2_lx[simHit2_N], b_simHit2_ly[simHit2_N];

  const static unsigned short simDigi1_N = 1000;
  unsigned short b_simDigi1_n;
  bool b_simDigi1_isBarrel[simDigi1_N], b_simDigi1_isIRPC[simDigi1_N];
  short b_simDigi1_sector[simDigi1_N], b_simDigi1_station[simDigi1_N], b_simDigi1_roll[simDigi1_N];
  short b_simDigi1_wheel[simDigi1_N], b_simDigi1_layer[simDigi1_N];
  short b_simDigi1_disk[simDigi1_N], b_simDigi1_ring[simDigi1_N];
  double b_simDigi1_x[simDigi1_N], b_simDigi1_y[simDigi1_N], b_simDigi1_z[simDigi1_N];
  double b_simDigi1_tof[simDigi1_N], b_simDigi1_t0[simDigi1_N];
  double b_simDigi1_lx[simDigi1_N], b_simDigi1_ly[simDigi1_N];
  short b_simDigi1_bx[simDigi1_N];

  const static unsigned short simDigi2_N = 1000;
  unsigned short b_simDigi2_n;
  bool b_simDigi2_isBarrel[simDigi2_N], b_simDigi2_isIRPC[simDigi2_N];
  short b_simDigi2_sector[simDigi2_N], b_simDigi2_station[simDigi2_N], b_simDigi2_roll[simDigi2_N];
  short b_simDigi2_wheel[simDigi2_N], b_simDigi2_layer[simDigi2_N];
  short b_simDigi2_disk[simDigi2_N], b_simDigi2_ring[simDigi2_N];
  double b_simDigi2_x[simDigi2_N], b_simDigi2_y[simDigi2_N], b_simDigi2_z[simDigi2_N];
  double b_simDigi2_tof[simDigi2_N], b_simDigi2_t0[simDigi2_N];
  double b_simDigi2_lx[simDigi2_N], b_simDigi2_ly[simDigi2_N];
  short b_simDigi2_bx[simDigi2_N];

  const static unsigned short rpcHit_N = 1000;
  unsigned short b_rpcHit_n;
  bool b_rpcHit_isBarrel[rpcHit_N], b_rpcHit_isIRPC[rpcHit_N];
  short b_rpcHit_sector[rpcHit_N], b_rpcHit_station[rpcHit_N], b_rpcHit_roll[rpcHit_N];
  short b_rpcHit_wheel[rpcHit_N], b_rpcHit_layer[rpcHit_N];
  short b_rpcHit_disk[rpcHit_N], b_rpcHit_ring[rpcHit_N];
  double b_rpcHit_x[rpcHit_N], b_rpcHit_y[rpcHit_N], b_rpcHit_z[rpcHit_N];
  double b_rpcHit_time[rpcHit_N], b_rpcHit_timeErr[rpcHit_N], b_rpcHit_t0[rpcHit_N];
  double b_rpcHit_lx[rpcHit_N], b_rpcHit_ly[rpcHit_N];
  short b_rpcHit_bx[rpcHit_N];

  const static unsigned short muon_N = 10;
  unsigned short b_muon_n;
  double b_muon_pt[muon_N], b_muon_eta[muon_N], b_muon_phi[muon_N];
  bool b_muon_isLoose[muon_N], b_muon_isTight[muon_N];
  double b_muon_time[muon_N], b_muon_RPCTime[muon_N], b_muon_RPCTimeNew[muon_N];
  unsigned short b_muon_nRPC[muon_N], b_muon_nIRPC[muon_N];
  double b_muon_genDR[muon_N];
  short b_muon_genPdgId[muon_N];

};

HSCPTreeMaker::HSCPTreeMaker(const edm::ParameterSet& pset):
  speedOfLight_(29.9792458), // unit in (cm/ns)
  signalSpeed_(speedOfLight_*pset.getParameter<double>("signalPropagationSpeed")),
  signalPdgId_(int(pset.getParameter<unsigned int>("signalPdgId"))),
  doSimHit_(pset.getParameter<bool>("doSimHit")),
  doSimDigi_(pset.getParameter<bool>("doSimDigi")),
  genParticleToken_(consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genParticle"))),
  //simVertexToken_(consumes<edm::SimVertexContainer>(pset.getParameter<edm::InputTag>("simVertex"))),
  rpcSimHitsToken_(consumes<edm::PSimHitContainer>(pset.getParameter<edm::InputTag>("rpcSimHits"))),
  //rpcDigisToken_(consumes<RPCDigiCollection>(pset.getParameter<edm::InputTag>("rpcDigis"))),
  rpcSimDigisToken_(consumes<edm::DetSetVector<RPCDigiSimLink>>(pset.getParameter<edm::InputTag>("rpcSimDigis"))),
  rpcRecHitsToken_(consumes<RPCRecHitCollection>(pset.getParameter<edm::InputTag>("rpcRecHits"))),
  muonToken_(consumes<reco::MuonCollection>(pset.getParameter<edm::InputTag>("muons"))),
  vertexToken_(consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertex")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("tree", "tree");

  tree_->Branch("gen1_pdgId", &b_gen1_pdgId, "gen1_pdgId/I");
  tree_->Branch("gen1_pt", &b_gen1_pt, "gen1_pt/D");
  tree_->Branch("gen1_eta", &b_gen1_eta, "gen1_eta/D");
  tree_->Branch("gen1_phi", &b_gen1_phi, "gen1_phi/D");
  tree_->Branch("gen1_m", &b_gen1_m, "gen1_m/D");
  tree_->Branch("gen1_beta", &b_gen1_beta, "gen1_beta/D");

  tree_->Branch("gen2_pdgId", &b_gen2_pdgId, "gen2_pdgId/I");
  tree_->Branch("gen2_pt", &b_gen2_pt, "gen2_pt/D");
  tree_->Branch("gen2_eta", &b_gen2_eta, "gen2_eta/D");
  tree_->Branch("gen2_phi", &b_gen2_phi, "gen2_phi/D");
  tree_->Branch("gen2_m", &b_gen2_m, "gen2_m/D");
  tree_->Branch("gen2_beta", &b_gen2_beta, "gen2_beta/D");

  if ( doSimHit_ ) {
    tree_->Branch("simHit1_n", &b_simHit1_n, "simHit1_n/s");
    tree_->Branch("simHit1_isBarrel", b_simHit1_isBarrel, "simHit1_isBarrel[simHit1_n]/O");
    tree_->Branch("simHit1_isIRPC", b_simHit1_isIRPC, "simHit1_isIRPC[simHit1_n]/O");
    tree_->Branch("simHit1_sector", b_simHit1_sector, "simHit1_sector[simHit1_n]/S");
    tree_->Branch("simHit1_station", b_simHit1_station, "simHit1_station[simHit1_n]/S");
    tree_->Branch("simHit1_wheel", b_simHit1_wheel, "simHit1_wheel[simHit1_n]/S");
    tree_->Branch("simHit1_layer", b_simHit1_layer, "simHit1_layer[simHit1_n]/S");
    tree_->Branch("simHit1_disk", b_simHit1_disk, "simHit1_disk[simHit1_n]/S");
    tree_->Branch("simHit1_ring", b_simHit1_ring, "simHit1_ring[simHit1_n]/S");
    tree_->Branch("simHit1_x", b_simHit1_x, "simHit1_x[simHit1_n]/D");
    tree_->Branch("simHit1_y", b_simHit1_y, "simHit1_y[simHit1_n]/D");
    tree_->Branch("simHit1_z", b_simHit1_z, "simHit1_z[simHit1_n]/D");
    tree_->Branch("simHit1_tof", b_simHit1_tof, "simHit1_tof[simHit1_n]/D");
    tree_->Branch("simHit1_lx", b_simHit1_lx, "simHit1_lx[simHit1_n]/D");
    tree_->Branch("simHit1_ly", b_simHit1_ly, "simHit1_ly[simHit1_n]/D");

    tree_->Branch("simHit2_n", &b_simHit2_n, "simHit2_n/s");
    tree_->Branch("simHit2_isBarrel", b_simHit2_isBarrel, "simHit2_isBarrel[simHit2_n]/O");
    tree_->Branch("simHit2_isIRPC", b_simHit2_isIRPC, "simHit2_isIRPC[simHit2_n]/O");
    tree_->Branch("simHit2_sector", b_simHit2_sector, "simHit2_sector[simHit2_n]/S");
    tree_->Branch("simHit2_station", b_simHit2_station, "simHit2_station[simHit2_n]/S");
    tree_->Branch("simHit2_wheel", b_simHit2_wheel, "simHit2_wheel[simHit2_n]/S");
    tree_->Branch("simHit2_layer", b_simHit2_layer, "simHit2_layer[simHit2_n]/S");
    tree_->Branch("simHit2_disk", b_simHit2_disk, "simHit2_disk[simHit2_n]/S");
    tree_->Branch("simHit2_ring", b_simHit2_ring, "simHit2_ring[simHit2_n]/S");
    tree_->Branch("simHit2_x", b_simHit2_x, "simHit2_x[simHit2_n]/D");
    tree_->Branch("simHit2_y", b_simHit2_y, "simHit2_y[simHit2_n]/D");
    tree_->Branch("simHit2_z", b_simHit2_z, "simHit2_z[simHit2_n]/D");
    tree_->Branch("simHit2_tof", b_simHit2_tof, "simHit2_tof[simHit2_n]/D");
    tree_->Branch("simHit2_lx", b_simHit2_lx, "simHit2_lx[simHit2_n]/D");
    tree_->Branch("simHit2_ly", b_simHit2_ly, "simHit2_ly[simHit2_n]/D");
  }

  if ( doSimDigi_ ) {
    tree_->Branch("simDigi1_n", &b_simDigi1_n, "simDigi1_n/s");
    tree_->Branch("simDigi1_isBarrel", b_simDigi1_isBarrel, "simDigi1_isBarrel[simDigi1_n]/O");
    tree_->Branch("simDigi1_isIRPC", b_simDigi1_isIRPC, "simDigi1_isIRPC[simDigi1_n]/O");
    tree_->Branch("simDigi1_sector", b_simDigi1_sector, "simDigi1_sector[simDigi1_n]/S");
    tree_->Branch("simDigi1_station", b_simDigi1_station, "simDigi1_station[simDigi1_n]/S");
    tree_->Branch("simDigi1_wheel", b_simDigi1_wheel, "simDigi1_wheel[simDigi1_n]/S");
    tree_->Branch("simDigi1_layer", b_simDigi1_layer, "simDigi1_layer[simDigi1_n]/S");
    tree_->Branch("simDigi1_disk", b_simDigi1_disk, "simDigi1_disk[simDigi1_n]/S");
    tree_->Branch("simDigi1_ring", b_simDigi1_ring, "simDigi1_ring[simDigi1_n]/S");
    tree_->Branch("simDigi1_x", b_simDigi1_x, "simDigi1_x[simDigi1_n]/D");
    tree_->Branch("simDigi1_y", b_simDigi1_y, "simDigi1_y[simDigi1_n]/D");
    tree_->Branch("simDigi1_z", b_simDigi1_z, "simDigi1_z[simDigi1_n]/D");
    tree_->Branch("simDigi1_tof", b_simDigi1_tof, "simDigi1_tof[simDigi1_n]/D");
    tree_->Branch("simDigi1_t0", b_simDigi1_t0, "simDigi1_t0[simDigi1_n]/D");
    tree_->Branch("simDigi1_lx", b_simDigi1_lx, "simDigi1_lx[simDigi1_n]/D");
    tree_->Branch("simDigi1_ly", b_simDigi1_ly, "simDigi1_ly[simDigi1_n]/D");
    tree_->Branch("simDigi1_bx", b_simDigi1_bx, "simDigi1_bx[simDigi1_n]/S");

    tree_->Branch("simDigi2_n", &b_simDigi2_n, "simDigi2_n/s");
    tree_->Branch("simDigi2_isBarrel", b_simDigi2_isBarrel, "simDigi2_isBarrel[simDigi2_n]/O");
    tree_->Branch("simDigi2_isIRPC", b_simDigi2_isIRPC, "simDigi2_isIRPC[simDigi2_n]/O");
    tree_->Branch("simDigi2_sector", b_simDigi2_sector, "simDigi2_sector[simDigi2_n]/S");
    tree_->Branch("simDigi2_station", b_simDigi2_station, "simDigi2_station[simDigi2_n]/S");
    tree_->Branch("simDigi2_wheel", b_simDigi2_wheel, "simDigi2_wheel[simDigi2_n]/S");
    tree_->Branch("simDigi2_layer", b_simDigi2_layer, "simDigi2_layer[simDigi2_n]/S");
    tree_->Branch("simDigi2_disk", b_simDigi2_disk, "simDigi2_disk[simDigi2_n]/S");
    tree_->Branch("simDigi2_ring", b_simDigi2_ring, "simDigi2_ring[simDigi2_n]/S");
    tree_->Branch("simDigi2_x", b_simDigi2_x, "simDigi2_x[simDigi2_n]/D");
    tree_->Branch("simDigi2_y", b_simDigi2_y, "simDigi2_y[simDigi2_n]/D");
    tree_->Branch("simDigi2_z", b_simDigi2_z, "simDigi2_z[simDigi2_n]/D");
    tree_->Branch("simDigi2_tof", b_simDigi2_tof, "simDigi2_tof[simDigi2_n]/D");
    tree_->Branch("simDigi2_t0", b_simDigi2_t0, "simDigi2_t0[simDigi2_n]/D");
    tree_->Branch("simDigi2_lx", b_simDigi2_lx, "simDigi2_lx[simDigi2_n]/D");
    tree_->Branch("simDigi2_ly", b_simDigi2_ly, "simDigi2_ly[simDigi2_n]/D");
    tree_->Branch("simDigi2_bx", b_simDigi2_bx, "simDigi2_bx[simDigi2_n]/S");
  }

  tree_->Branch("rpcHit_n", &b_rpcHit_n, "rpcHit_n/s");
  tree_->Branch("rpcHit_isBarrel", b_rpcHit_isBarrel, "rpcHit_isBarrel[rpcHit_n]/O");
  tree_->Branch("rpcHit_isIRPC", b_rpcHit_isIRPC, "rpcHit_isIRPC[rpcHit_n]/O");
  tree_->Branch("rpcHit_sector", b_rpcHit_sector, "rpcHit_sector[rpcHit_n]/S");
  tree_->Branch("rpcHit_station", b_rpcHit_station, "rpcHit_station[rpcHit_n]/S");
  tree_->Branch("rpcHit_wheel", b_rpcHit_wheel, "rpcHit_wheel[rpcHit_n]/S");
  tree_->Branch("rpcHit_layer", b_rpcHit_layer, "rpcHit_layer[rpcHit_n]/S");
  tree_->Branch("rpcHit_disk", b_rpcHit_disk, "rpcHit_disk[rpcHit_n]/S");
  tree_->Branch("rpcHit_ring", b_rpcHit_ring, "rpcHit_ring[rpcHit_n]/S");
  tree_->Branch("rpcHit_x", b_rpcHit_x, "rpcHit_x[rpcHit_n]/D");
  tree_->Branch("rpcHit_y", b_rpcHit_y, "rpcHit_y[rpcHit_n]/D");
  tree_->Branch("rpcHit_z", b_rpcHit_z, "rpcHit_z[rpcHit_n]/D");
  tree_->Branch("rpcHit_time", b_rpcHit_time, "rpcHit_time[rpcHit_n]/D");
  tree_->Branch("rpcHit_timeErr", b_rpcHit_timeErr, "rpcHit_timeErr[rpcHit_n]/D");
  tree_->Branch("rpcHit_t0", b_rpcHit_t0, "rpcHit_t0[rpcHit_n]/D");
  tree_->Branch("rpcHit_lx", b_rpcHit_lx, "rpcHit_lx[rpcHit_n]/D");
  tree_->Branch("rpcHit_ly", b_rpcHit_ly, "rpcHit_ly[rpcHit_n]/D");
  tree_->Branch("rpcHit_bx", b_rpcHit_bx, "rpcHit_bx[rpcHit_n]/S");

  tree_->Branch("muon_n", &b_muon_n, "muon_n/s");
  tree_->Branch("muon_pt", b_muon_pt, "muon_pt[muon_n]/D");
  tree_->Branch("muon_eta", b_muon_eta, "muon_eta[muon_n]/D");
  tree_->Branch("muon_phi", b_muon_phi, "muon_phi[muon_n]/D");
  tree_->Branch("muon_isLoose", b_muon_isLoose, "muon_isLoose[muon_n]/O");
  tree_->Branch("muon_isTight", b_muon_isTight, "muon_isTight[muon_n]/O");
  tree_->Branch("muon_time", b_muon_time, "muon_time[muon_n]/D");
  tree_->Branch("muon_RPCTime", b_muon_RPCTime, "muon_RPCTime[muon_n]/D");
  tree_->Branch("muon_RPCTimeNew", b_muon_RPCTimeNew, "muon_RPCTimeNew[muon_n]/D");
  tree_->Branch("muon_nRPC", b_muon_nRPC, "muon_nRPC[muon_n]/s");
  tree_->Branch("muon_nIRPC", b_muon_nIRPC, "muon_nIRPC[muon_n]/s");
  tree_->Branch("muon_genDR", b_muon_genDR, "muon_genDR[muon_n]/D");
  tree_->Branch("muon_genPdgId", b_muon_genPdgId, "muon_genPdgId[muon_n]/S");

}

template<typename T1, typename T2>
double dist(const T1& a, const T2& b)
{
  const double dx = a.x()-b.x();
  const double dy = a.y()-b.y();
  const double dz = a.z()-b.z();
  return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void HSCPTreeMaker::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  b_simHit1_n = b_simHit2_n = b_simDigi1_n = b_simDigi2_n = 0;
  b_rpcHit_n = 0;
  b_muon_n = 0;

  b_gen1_pdgId = b_gen1_pt = b_gen1_eta = b_gen1_phi = b_gen1_m = b_gen1_beta = 0;
  b_gen2_pdgId = b_gen2_pt = b_gen2_eta = b_gen2_phi = b_gen2_m = b_gen2_beta = 0;

  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  event.getByToken(genParticleToken_, genParticleHandle);

  //edm::Handle<edm::SimVertexContainer> simVertexHandle;
  //event.getByToken(simVertexToken_, simVertexHandle);

  edm::Handle<edm::PSimHitContainer> rpcSimHitsHandle;
  if ( doSimHit_ ) event.getByToken(rpcSimHitsToken_, rpcSimHitsHandle);

  edm::Handle<edm::DetSetVector<RPCDigiSimLink>> rpcSimDigisHandle;
  if ( doSimDigi_ ) event.getByToken(rpcSimDigisToken_, rpcSimDigisHandle);

  edm::Handle<RPCRecHitCollection> rpcRecHitsHandle;
  event.getByToken(rpcRecHitsToken_, rpcRecHitsHandle);

  edm::Handle<reco::MuonCollection> muonHandle;
  event.getByToken(muonToken_, muonHandle);

  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByToken(vertexToken_, vertexHandle);

  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  /*const math::XYZTLorentzVector simPV = [&](){
    if ( !simVertexHandle.isValid() or
         simVertexHandle->empty() ) return math::XYZTLorentzVector();
    return simVertexHandle->at(0).position();
  }();*/

  const reco::GenParticle* genParticle1 = 0, * genParticle2 = 0;
  if ( genParticleHandle.isValid() ) {
    for ( auto& p : *genParticleHandle ) {
      if ( p.status() != 1 ) continue;

      if ( p.pdgId() == signalPdgId_ ) {
        if ( !genParticle1 or genParticle1->pt() < p.pt() ) genParticle1 = &p;
      }
      else if ( p.pdgId() == -signalPdgId_ ) {
        if ( !genParticle2 or genParticle2->pt() < p.pt() ) genParticle2 = &p;
      }
    }
  }
  if ( genParticle1 ) {
    const auto& p = *genParticle1;
    b_gen1_pdgId = p.pdgId();
    b_gen1_pt = p.pt();
    b_gen1_eta = p.eta();
    b_gen1_phi = p.phi();
    b_gen1_m = p.mass();
    b_gen1_beta = 1/hypot(1, p.mass()/p.p());
  }
  if ( genParticle2 ) {
    const auto& p = *genParticle2;
    b_gen2_pdgId = p.pdgId();
    b_gen2_pt = p.pt();
    b_gen2_eta = p.eta();
    b_gen2_phi = p.phi();
    b_gen2_m = p.mass();
    b_gen2_beta = 1/hypot(1, p.mass()/p.p());
  }

  if ( doSimHit_ and rpcSimHitsHandle.isValid() ) {
    for ( auto& simHit : *rpcSimHitsHandle ) {
      const int pid = simHit.particleType();
      if ( std::abs(pid) == signalPdgId_ ) continue;

      const RPCDetId rpcId(simHit.detUnitId());
      const RPCRoll* roll = rpcGeom->roll(rpcId);
      const auto simHitGPos = roll->toGlobal(simHit.localPosition());
      //const double r = dist(simPV, simHitGPos);
      //const double r0 = simHitGPos.mag();
      //const double beta = r/simHit.tof()/speedOfLight;

      if ( pid == signalPdgId_ and b_simHit1_n < simHit1_N ) {
        b_simHit1_isBarrel[b_simHit1_n] = (rpcId.region() == 0);
        b_simHit1_isIRPC[b_simHit1_n] = roll->isIRPC();
        b_simHit1_sector[b_simHit1_n] = rpcId.sector();
        b_simHit1_layer[b_simHit1_n] = rpcId.layer();
        b_simHit1_roll[b_simHit1_n] = rpcId.roll();
        b_simHit1_x[b_simHit1_n] = simHitGPos.x();
        b_simHit1_y[b_simHit1_n] = simHitGPos.y();
        b_simHit1_z[b_simHit1_n] = simHitGPos.z();
        b_simHit1_tof[b_simHit1_n] = simHit.tof();
        b_simHit1_lx[b_simHit1_n] = simHit.localPosition().x();
        b_simHit1_ly[b_simHit1_n] = simHit.localPosition().y();

        if ( rpcId.region() == 0 ) {
          b_simHit1_disk[b_simHit1_n] = b_simHit1_ring[b_simHit1_n] = 0;
          b_simHit1_wheel[b_simHit1_n] = rpcId.ring();
          b_simHit1_station[b_simHit1_n] = rpcId.station();
        }
        else {
          b_simHit1_wheel[b_simHit1_n] = b_simHit1_station[b_simHit1_n] = 0;
          b_simHit1_disk[b_simHit1_n] = rpcId.region()*rpcId.station();
          b_simHit1_ring[b_simHit1_n] = rpcId.ring();
        }

        ++b_simHit1_n;
      }
      else if ( pid == -signalPdgId_ and b_simHit2_n < simHit2_N ) {
        b_simHit2_isBarrel[b_simHit2_n] = (rpcId.region() == 0);
        b_simHit2_isIRPC[b_simHit2_n] = roll->isIRPC();
        b_simHit2_sector[b_simHit2_n] = rpcId.sector();
        b_simHit2_layer[b_simHit2_n] = rpcId.layer();
        b_simHit2_roll[b_simHit2_n] = rpcId.roll();
        b_simHit2_x[b_simHit2_n] = simHitGPos.x();
        b_simHit2_y[b_simHit2_n] = simHitGPos.y();
        b_simHit2_z[b_simHit2_n] = simHitGPos.z();
        b_simHit2_tof[b_simHit2_n] = simHit.tof();
        b_simHit2_lx[b_simHit2_n] = simHit.localPosition().x();
        b_simHit2_ly[b_simHit2_n] = simHit.localPosition().y();

        if ( rpcId.region() == 0 ) {
          b_simHit2_disk[b_simHit2_n] = b_simHit2_ring[b_simHit2_n] = 0;
          b_simHit2_wheel[b_simHit2_n] = rpcId.ring();
          b_simHit2_station[b_simHit2_n] = rpcId.station();
        }
        else {
          b_simHit2_wheel[b_simHit2_n] = b_simHit2_station[b_simHit2_n] = 0;
          b_simHit2_disk[b_simHit2_n] = rpcId.station();
          b_simHit2_ring[b_simHit2_n] = rpcId.ring();
        }

        ++b_simHit2_n;
      }
    }
  }

  if ( doSimDigi_ and rpcSimDigisHandle.isValid() ) {
    for ( const auto& detSet : *rpcSimDigisHandle ) {
      const unsigned int rawId = detSet.detId();
      if ( DetId(rawId).det() != DetId::Muon or
           DetId(rawId).subdetId() != 3 ) continue;

      const RPCDetId rpcId(rawId);
      const RPCRoll* roll = rpcGeom->roll(rpcId);
      if ( !roll ) continue;
      const double r0 = roll->toGlobal(LocalPoint(0,0,0)).mag();

      for ( const auto& simDigi : detSet ) {
        const int pid = simDigi.getParticleType();
        if ( std::abs(pid) != signalPdgId_ ) continue;

        const auto lp = simDigi.getEntryPoint();
        const auto gp = roll->toGlobal(lp);
        const double tof = simDigi.getTimeOfFlight();

        if ( pid == signalPdgId_ and b_simDigi1_n < simDigi1_N ) {
          b_simDigi1_isBarrel[b_simDigi1_n] = (rpcId.region() == 0);
          b_simDigi1_isIRPC[b_simDigi1_n] = roll->isIRPC();
          b_simDigi1_sector[b_simDigi1_n] = rpcId.sector();
          b_simDigi1_layer[b_simDigi1_n] = rpcId.layer();
          b_simDigi1_roll[b_simDigi1_n] = rpcId.roll();
          b_simDigi1_x[b_simDigi1_n] = gp.x();
          b_simDigi1_y[b_simDigi1_n] = gp.y();
          b_simDigi1_z[b_simDigi1_n] = gp.z();
          b_simDigi1_tof[b_simDigi1_n] = tof;
          b_simDigi1_t0[b_simDigi1_n] = tof - r0/speedOfLight_;
          b_simDigi1_lx[b_simDigi1_n] = lp.x();
          b_simDigi1_ly[b_simDigi1_n] = lp.y();
          b_simDigi1_bx[b_simDigi1_n] = simDigi.getBx();

          if ( rpcId.region() == 0 ) {
            b_simDigi1_disk[b_simDigi1_n] = b_simDigi1_ring[b_simDigi1_n] = 0;
            b_simDigi1_wheel[b_simDigi1_n] = rpcId.ring();
            b_simDigi1_station[b_simDigi1_n] = rpcId.station();
          }
          else {
            b_simDigi1_wheel[b_simDigi1_n] = b_simDigi1_station[b_simDigi1_n] = 0;
            b_simDigi1_disk[b_simDigi1_n] = rpcId.region()*rpcId.station();
            b_simDigi1_ring[b_simDigi1_n] = rpcId.ring();
          }

          ++b_simDigi1_n;
        }
        else if ( pid == -signalPdgId_ and b_simDigi2_n < simDigi2_N ) {
          b_simDigi2_isBarrel[b_simDigi2_n] = (rpcId.region() == 0);
          b_simDigi2_isIRPC[b_simDigi2_n] = roll->isIRPC();
          b_simDigi2_sector[b_simDigi2_n] = rpcId.sector();
          b_simDigi2_layer[b_simDigi2_n] = rpcId.layer();
          b_simDigi2_roll[b_simDigi2_n] = rpcId.roll();
          b_simDigi2_x[b_simDigi2_n] = gp.x();
          b_simDigi2_y[b_simDigi2_n] = gp.y();
          b_simDigi2_z[b_simDigi2_n] = gp.z();
          b_simDigi2_tof[b_simDigi2_n] = tof;
          b_simDigi2_t0[b_simDigi2_n] = tof - r0/speedOfLight_;
          b_simDigi2_lx[b_simDigi2_n] = lp.x();
          b_simDigi2_ly[b_simDigi2_n] = lp.y();
          b_simDigi2_bx[b_simDigi2_n] = simDigi.getBx();

          if ( rpcId.region() == 0 ) {
            b_simDigi2_disk[b_simDigi2_n] = b_simDigi2_ring[b_simDigi2_n] = 0;
            b_simDigi2_wheel[b_simDigi2_n] = rpcId.ring();
            b_simDigi2_station[b_simDigi2_n] = rpcId.station();
          }
          else {
            b_simDigi2_wheel[b_simDigi2_n] = b_simDigi2_station[b_simDigi2_n] = 0;
            b_simDigi2_disk[b_simDigi2_n] = rpcId.region()*rpcId.station();
            b_simDigi2_ring[b_simDigi2_n] = rpcId.ring();
          }

          ++b_simDigi2_n;
        }
      }
    }
  }

  if ( muonHandle.isValid() and vertexHandle.isValid() and rpcRecHitsHandle.isValid() ) {
    const reco::Vertex* pv = 0;
    for ( auto& vtx : *vertexHandle ) {
      if ( vtx.isFake() ) continue;
      if ( vtx.ndof() <= 4 ) continue;
      if ( std::abs(vtx.z()) > 24 ) continue;
      if ( std::abs(vtx.position().rho()) > 2 ) continue;

      pv = &vtx;
      break;
    }

    for ( auto& mu : *muonHandle ) {
      if ( std::abs(mu.eta()) >= 2.5 or mu.pt() < 20 ) continue;

      b_muon_pt[b_muon_n] = mu.pt();
      b_muon_eta[b_muon_n] = mu.eta();
      b_muon_phi[b_muon_n] = mu.phi();
      b_muon_isLoose[b_muon_n] = muon::isLooseMuon(mu);
      b_muon_isTight[b_muon_n] = !pv ? false : muon::isTightMuon(mu, *pv);

      b_muon_time[b_muon_n] = mu.time().timeAtIpInOut;
      b_muon_RPCTime[b_muon_n] = mu.rpcTime().timeAtIpInOut;

      const std::vector<std::pair<double, int>> hitTimes = getRPCTimes(mu, *rpcRecHitsHandle, rpcGeom, true);
      const double sumTime = std::accumulate(hitTimes.begin(), hitTimes.end(), 0.0, [](const double& a, const std::pair<double, int>& b){return a+b.first;});
      b_muon_RPCTimeNew[b_muon_n] = hitTimes.empty() ? 0 : sumTime/hitTimes.size();
      b_muon_nRPC[b_muon_n] = hitTimes.size();
      b_muon_nIRPC[b_muon_n] = std::count_if(hitTimes.begin(), hitTimes.end(), [](const std::pair<double, int>& p){return p.second == 2;});

      const double dR1 = !genParticle1 ? 999 : deltaR(mu, *genParticle1);
      const double dR2 = !genParticle2 ? 999 : deltaR(mu, *genParticle2);
      if ( dR1 < 999 and dR1 < dR2 ) {
        b_muon_genDR[b_muon_n] = dR1;
        b_muon_genPdgId[b_muon_n] = genParticle1->pdgId();
      }
      else if ( dR2 < 999 and dR2 < dR1 ) {
        b_muon_genDR[b_muon_n] = dR2;
        b_muon_genPdgId[b_muon_n] = genParticle2->pdgId();
      }

      if ( ++b_muon_n >= muon_N ) break;
    }
  }

  if ( rpcRecHitsHandle.isValid() ) {
    for ( auto& rpcHit : *rpcRecHitsHandle ) {
      const RPCDetId rpcId = rpcHit.rpcId();
      const RPCRoll* roll = rpcGeom->roll(rpcId);
      const auto gp = roll->toGlobal(rpcHit.localPosition());

      b_rpcHit_isBarrel[b_rpcHit_n] = (rpcId.region() == 0);
      b_rpcHit_isIRPC[b_rpcHit_n] = roll->isIRPC();
      b_rpcHit_sector[b_rpcHit_n] = rpcId.sector();
      b_rpcHit_layer[b_rpcHit_n] = rpcId.layer();
      b_rpcHit_roll[b_rpcHit_n] = rpcId.roll();
      b_rpcHit_x[b_rpcHit_n] = gp.x();
      b_rpcHit_y[b_rpcHit_n] = gp.y();
      b_rpcHit_z[b_rpcHit_n] = gp.z();
      b_rpcHit_time[b_rpcHit_n] = rpcHit.time();
      b_rpcHit_timeErr[b_rpcHit_n] = rpcHit.timeError();
      b_rpcHit_lx[b_rpcHit_n] = rpcHit.localPosition().x();
      b_rpcHit_ly[b_rpcHit_n] = rpcHit.localPosition().y();
      b_rpcHit_bx[b_rpcHit_n] = rpcHit.BunchX();

      if ( rpcId.region() == 0 ) {
        b_rpcHit_disk[b_rpcHit_n] = b_rpcHit_ring[b_rpcHit_n] = 0;
        b_rpcHit_wheel[b_rpcHit_n] = rpcId.ring();
        b_rpcHit_station[b_rpcHit_n] = rpcId.station();
      }
      else {
        b_rpcHit_wheel[b_rpcHit_n] = b_rpcHit_station[b_rpcHit_n] = 0;
        b_rpcHit_disk[b_rpcHit_n] = rpcId.region()*rpcId.station();
        b_rpcHit_ring[b_rpcHit_n] = rpcId.ring();
      }

      if ( ++b_rpcHit_n >= rpcHit_N ) break;
    }
  }

  tree_->Fill();
}

std::vector<std::pair<double, int>> HSCPTreeMaker::getRPCTimes(const reco::Muon& mu, 
                                                               const RPCRecHitCollection& rpcHits, edm::ESHandle<RPCGeometry>& rpcGeom,
                                                               const bool doOnlyWithTime) const
{
  std::vector<std::pair<double, int>> hitTimes;

  reco::TrackRef trackRef;
  if ( mu.isGlobalMuon() ) trackRef = mu.globalTrack();
  else if ( mu.isStandAloneMuon() ) trackRef = mu.standAloneMuon();

  if ( trackRef.isNull() ) return hitTimes;

  for ( auto ihit = trackRef->recHitsBegin(); ihit != trackRef->recHitsEnd(); ++ihit ) {
    if ( !(*ihit)->isValid() ) continue;
    if ( (*ihit)->geographicalId().det() != DetId::Muon or (*ihit)->geographicalId().subdetId() != 3 ) continue;

    const RPCDetId detId((*ihit)->rawId());
    const RPCRoll* roll = rpcGeom->roll(detId);

    const auto hitRange = rpcHits.get(detId);
    for ( auto rpcHit = hitRange.first; rpcHit != hitRange.second; ++rpcHit ) {
      // "detType" is 0 if no time info, 1 with time info after LB upgrade, 2 for iRPC
      const int detType = rpcHit->timeError() >= 0 ? roll->isIRPC() ? 2 : 1 : 0;
      if ( doOnlyWithTime and detType == 0 ) continue; 

      const double time = detType != 0 ? rpcHit->time() : 25.*rpcHit->BunchX();
      hitTimes.emplace_back(time, detType);
    }
  }

  return hitTimes;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HSCPTreeMaker);


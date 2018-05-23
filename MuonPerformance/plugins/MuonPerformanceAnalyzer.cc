#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "SimMuon/MCTruth/interface/MuonToSimAssociatorByHits.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"

#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0Chamber.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMSuperChamber.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCChamber.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TTree.h>

struct RPCHitInfo
{
  RPCHitInfo(float _x, float _y, float _t,
             float _gx, float _gy, float _gz, float _x0, bool _isIRPC):
    x(_x), y(_y), t(_t), gx(_gx), gy(_gy), gz(_gz), x0(_x0), isIRPC(_isIRPC) {};
  float x, y, t;
  float gx, gy, gz, x0;
  bool isIRPC;
};

class MuonPerformanceAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  MuonPerformanceAnalyzer(const edm::ParameterSet& pset);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  std::vector<RPCHitInfo> collectRPCHits(const RPCGeometry* rpcGeom, const reco::TrackRef track, const RPCRecHitCollection* rpcHits) const;

private:
  edm::EDGetTokenT<edm::View<reco::Muon>> muonsToken_;
  edm::EDGetTokenT<RPCRecHitCollection> rpcHitsToken_;
  edm::EDGetTokenT<TrackingParticleCollection> simTPsToken_;
  edm::EDGetTokenT<reco::MuonToTrackingParticleAssociator> muAssocToken_;
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;

  TrackingParticleSelector tpSelector_;

private:
  TTree* tree_;

  const static unsigned short muons_N = 100;
  unsigned short b_muons_n;
  float b_muons_pt[muons_N], b_muons_eta[muons_N], b_muons_phi[muons_N];
  float b_muons_time[muons_N], b_muons_rpcTime[muons_N], b_muons_rpcTime1[muons_N];
  float b_muons_sta_pt[muons_N], b_muons_sta_eta[muons_N], b_muons_sta_phi[muons_N];
  float b_muons_glb_pt[muons_N], b_muons_glb_eta[muons_N], b_muons_glb_phi[muons_N];
  float b_muons_gen_pt[muons_N], b_muons_gen_eta[muons_N], b_muons_gen_phi[muons_N];
  float b_muons_glb_ptErr[muons_N], b_muons_sta_ptErr[muons_N];

  unsigned int b_muons_hits_n[muons_N], b_muons_muonHits_n[muons_N];
  unsigned int b_muons_DTHits_n[muons_N], b_muons_CSCHits_n[muons_N];
  unsigned int b_muons_RPCHits_n[muons_N];
  unsigned int b_muons_GEMHits_n[muons_N], b_muons_ME0Hits_n[muons_N];

  unsigned int b_muons_glb_hits_n[muons_N], b_muons_glb_muonHits_n[muons_N];
  unsigned int b_muons_glb_DTHits_n[muons_N], b_muons_glb_CSCHits_n[muons_N];
  unsigned int b_muons_glb_RPCHits_n[muons_N], b_muons_glb_iRPCHits_n[muons_N];
  unsigned int b_muons_glb_GEMHits_n[muons_N], b_muons_glb_ME0Hits_n[muons_N];

  unsigned int b_muons_sta_hits_n[muons_N], b_muons_sta_muonHits_n[muons_N];
  unsigned int b_muons_sta_DTHits_n[muons_N], b_muons_sta_CSCHits_n[muons_N];
  unsigned int b_muons_sta_RPCHits_n[muons_N], b_muons_sta_iRPCHits_n[muons_N];
  unsigned int b_muons_sta_GEMHits_n[muons_N], b_muons_sta_ME0Hits_n[muons_N];

  bool b_muons_isGen[muons_N];
  bool b_muons_isGlb[muons_N], b_muons_isSta[muons_N], b_muons_isTrk[muons_N], b_muons_isRPC[muons_N];
  bool b_muons_isLoose[muons_N], b_muons_isTight[muons_N];

  unsigned int b_run;
  unsigned long b_event;
  unsigned short b_vertices_n;

}; 

MuonPerformanceAnalyzer::MuonPerformanceAnalyzer(const edm::ParameterSet& pset)
{
  muonsToken_ = consumes<edm::View<reco::Muon>>(pset.getParameter<edm::InputTag>("muons"));
  rpcHitsToken_ = consumes<RPCRecHitCollection>(pset.getParameter<edm::InputTag>("rpcHits"));
  simTPsToken_ = consumes<TrackingParticleCollection>(pset.getParameter<edm::InputTag>("simTPs"));
  muAssocToken_ = consumes<reco::MuonToTrackingParticleAssociator>(pset.getParameter<edm::InputTag>("muAssoc"));
  verticesToken_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertices"));

  auto tpSet = pset.getParameter<edm::ParameterSet>("tpSelector");
  tpSelector_ = TrackingParticleSelector(tpSet.getParameter<double>("ptMin"), 1e9,
                                         tpSet.getParameter<double>("minRapidity"),
                                         tpSet.getParameter<double>("maxRapidity"),
                                         tpSet.getParameter<double>("tip"),
                                         tpSet.getParameter<double>("lip"),
                                         tpSet.getParameter<int>("minHit"),
                                         tpSet.getParameter<bool>("signalOnly"),
                                         tpSet.getParameter<bool>("intimeOnly"),
                                         tpSet.getParameter<bool>("chargedOnly"),
                                         tpSet.getParameter<bool>("stableOnly"),
                                         tpSet.getParameter<std::vector<int> >("pdgId"));

  // Book TTree
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "tree");

  tree_->Branch("run", &b_run, "run/i");
  tree_->Branch("event", &b_event, "event/l");
  tree_->Branch("muons_n", &b_muons_n, "muons_n/s");

  tree_->Branch("muons_isGlb", b_muons_isGlb, "muons_isGlb[muons_n]/O");
  tree_->Branch("muons_isGen", b_muons_isGen, "muons_isGen[muons_n]/O");
  tree_->Branch("muons_isSta", b_muons_isSta, "muons_isSta[muons_n]/O");
  tree_->Branch("muons_isTrk", b_muons_isTrk, "muons_isTrk[muons_n]/O");
  tree_->Branch("muons_isRPC", b_muons_isRPC, "muons_isRPC[muons_n]/O");
  tree_->Branch("muons_isLoose", b_muons_isLoose, "muons_isLoose[muons_n]/O");
  tree_->Branch("muons_isTight", b_muons_isTight, "muons_isTight[muons_n]/O");

  tree_->Branch("muons_pt", b_muons_pt, "muons_pt[muons_n]/F");
  tree_->Branch("muons_eta", b_muons_eta, "muons_eta[muons_n]/F");
  tree_->Branch("muons_phi", b_muons_phi, "muons_phi[muons_n]/F");

  tree_->Branch("muons_time", b_muons_time, "muons_time[muons_n]/F");
  tree_->Branch("muons_rpcTime", b_muons_rpcTime, "muons_rpcTime[muons_n]/F");
  tree_->Branch("muons_rpcTime1", b_muons_rpcTime1, "muons_rpcTime1[muons_n]/F");

  tree_->Branch("muons_sta_pt" , b_muons_sta_pt, "muons_sta_pt[muons_n]/F");
  tree_->Branch("muons_sta_eta", b_muons_sta_eta, "muons_sta_eta[muons_n]/F");
  tree_->Branch("muons_sta_phi", b_muons_sta_phi, "muons_sta_phi[muons_n]/F");
  tree_->Branch("muons_sta_ptErr", b_muons_sta_ptErr, "muons_sta_ptErr[muons_n]/F");

  tree_->Branch("muons_glb_pt" , b_muons_glb_pt, "muons_glb_pt[muons_n]/F");
  tree_->Branch("muons_glb_eta", b_muons_glb_eta, "muons_glb_eta[muons_n]/F");
  tree_->Branch("muons_glb_phi", b_muons_glb_phi, "muons_glb_phi[muons_n]/F");
  tree_->Branch("muons_glb_ptErr", b_muons_glb_ptErr, "muons_glb_ptErr[muons_n]/F");

  tree_->Branch("muons_gen_pt" , b_muons_gen_pt, "muons_gen_pt[muons_n]/F");
  tree_->Branch("muons_gen_eta", b_muons_gen_eta, "muons_gen_eta[muons_n]/F");
  tree_->Branch("muons_gen_phi", b_muons_gen_phi, "muons_gen_phi[muons_n]/F");

  tree_->Branch("muons_hits_n"    , b_muons_hits_n    , "muons_hits_n[muons_n]/s"    );
  tree_->Branch("muons_muonHits_n", b_muons_muonHits_n, "muons_muonHits_n[muons_n]/s");
  tree_->Branch("muons_DTHits_n"  , b_muons_DTHits_n  , "muons_DTHits_n[muons_n]/s"  );
  tree_->Branch("muons_CSCHits_n" , b_muons_CSCHits_n , "muons_CSCHits_n[muons_n]/s" );
  tree_->Branch("muons_RPCHits_n" , b_muons_RPCHits_n , "muons_RPCHits_n[muons_n]/s" );
  tree_->Branch("muons_GEMHits_n" , b_muons_GEMHits_n , "muons_GEMHits_n[muons_n]/s" );
  tree_->Branch("muons_ME0Hits_n" , b_muons_ME0Hits_n , "muons_ME0Hits_n[muons_n]/s" );

  tree_->Branch("muons_glb_hits_n"    , b_muons_glb_hits_n    , "muons_glb_hits_n[muons_n]/s"    );
  tree_->Branch("muons_glb_muonHits_n", b_muons_glb_muonHits_n, "muons_glb_muonHits_n[muons_n]/s");
  tree_->Branch("muons_glb_DTHits_n"  , b_muons_glb_DTHits_n  , "muons_glb_DTHits_n[muons_n]/s"  );
  tree_->Branch("muons_glb_CSCHits_n" , b_muons_glb_CSCHits_n , "muons_glb_CSCHits_n[muons_n]/s" );
  tree_->Branch("muons_glb_RPCHits_n" , b_muons_glb_RPCHits_n , "muons_glb_RPCHits_n[muons_n]/s" );
  tree_->Branch("muons_glb_iRPCHits_n", b_muons_glb_iRPCHits_n, "muons_glb_iRPCHits_n[muons_n]/s");
  tree_->Branch("muons_glb_GEMHits_n" , b_muons_glb_GEMHits_n , "muons_glb_GEMHits_n[muons_n]/s" );
  tree_->Branch("muons_glb_ME0Hits_n" , b_muons_glb_ME0Hits_n , "muons_glb_ME0Hits_n[muons_n]/s" );

  tree_->Branch("muons_sta_hits_n"    , b_muons_sta_hits_n    , "muons_sta_hits_n[muons_n]/s"    );
  tree_->Branch("muons_sta_muonHits_n", b_muons_sta_muonHits_n, "muons_sta_muonHits_n[muons_n]/s");
  tree_->Branch("muons_sta_DTHits_n"  , b_muons_sta_DTHits_n  , "muons_sta_DTHits_n[muons_n]/s"  );
  tree_->Branch("muons_sta_CSCHits_n" , b_muons_sta_CSCHits_n , "muons_sta_CSCHits_n[muons_n]/s" );
  tree_->Branch("muons_sta_RPCHits_n" , b_muons_sta_RPCHits_n , "muons_sta_RPCHits_n[muons_n]/s" );
  tree_->Branch("muons_sta_iRPCHits_n", b_muons_sta_iRPCHits_n, "muons_sta_iRPCHits_n[muons_n]/s");
  tree_->Branch("muons_sta_GEMHits_n" , b_muons_sta_GEMHits_n , "muons_sta_GEMHits_n[muons_n]/s" );
  tree_->Branch("muons_sta_ME0Hits_n" , b_muons_sta_ME0Hits_n , "muons_sta_ME0Hits_n[muons_n]/s" );

  tree_->Branch("vertices_n", &b_vertices_n, "vertices_n/s");
}

void MuonPerformanceAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  b_muons_n = 0;
  b_run = event.id().run();
  b_event = event.id().event();

  edm::Handle<reco::VertexCollection> verticesHandle;
  event.getByToken(verticesToken_, verticesHandle);
  b_vertices_n = verticesHandle->size();
  const reco::Vertex* pv = verticesHandle->empty() ? 0 : &verticesHandle->at(0);

  edm::Handle<TrackingParticleCollection> simTPsHandle;
  event.getByToken(simTPsToken_, simTPsHandle);

  edm::Handle<edm::View<reco::Muon>> muonsHandle;
  event.getByToken(muonsToken_, muonsHandle);

  edm::Handle<RPCRecHitCollection> rpcHitsHandle;
  event.getByToken(rpcHitsToken_, rpcHitsHandle);

  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  //edm::ESHandle<RPCGeometry> gemGeom;
  //eventSetup.get<MuonGeometryRecord>().get(gemGeom);

  edm::Handle<reco::MuonToTrackingParticleAssociator> muAssocHandle;
  event.getByToken(muAssocToken_, muAssocHandle);
  reco::MuonToSimCollection glbToSimColl, staToSimColl, trkToSimColl;
  reco::SimToMuonCollection simToGlbColl, simToStaColl, simToTrkColl;
  muAssocHandle->associateMuons(glbToSimColl, simToGlbColl, muonsHandle, reco::GlobalTk, simTPsHandle);
  muAssocHandle->associateMuons(staToSimColl, simToStaColl, muonsHandle, reco::OuterTk, simTPsHandle);
  muAssocHandle->associateMuons(trkToSimColl, simToTrkColl, muonsHandle, reco::InnerTk, simTPsHandle);

  for ( size_t i=0, n=simTPsHandle->size(); i<n; ++i ) {
    TrackingParticleRef simTP(simTPsHandle, i);
    if ( std::abs(simTP->pdgId()) != 13 or simTP->genParticles().empty() ) continue;
    if ( simTP->eventId().event() != 0 or simTP->eventId().bunchCrossing() != 0 ) continue;

    edm::RefToBase<reco::Muon> muonRef;
    if ( simToGlbColl.find(simTP) != simToGlbColl.end() ) {
      auto matches = simToGlbColl[simTP];
      if ( !matches.empty() ) muonRef = matches.begin()->first;
    }
    else if ( simToStaColl.find(simTP) != simToStaColl.end() ) {
      auto matches = simToStaColl[simTP];
      if ( !matches.empty() ) muonRef = matches.begin()->first;
    }
    else if ( simToTrkColl.find(simTP) != simToTrkColl.end() ) {
      auto matches = simToTrkColl[simTP];
      if ( !matches.empty() ) muonRef = matches.begin()->first;
    }
    
    b_muons_isGen[b_muons_n] = true;
    b_muons_isSta[b_muons_n] = b_muons_isGlb[b_muons_n] = b_muons_isTrk[b_muons_n] = b_muons_isRPC[b_muons_n] = false;
    b_muons_isLoose[b_muons_n] = b_muons_isTight[b_muons_n] = false;

    b_muons_pt[b_muons_n] = b_muons_eta[b_muons_n] = b_muons_phi[b_muons_n] = -1e5;
    b_muons_time[b_muons_n] = b_muons_rpcTime[b_muons_n] = b_muons_rpcTime1[b_muons_n] = -1e5;
    b_muons_glb_pt[b_muons_n] = b_muons_glb_eta[b_muons_n] = b_muons_glb_phi[b_muons_n] = b_muons_glb_ptErr[b_muons_n] = -1e5;
    b_muons_sta_pt[b_muons_n] = b_muons_sta_eta[b_muons_n] = b_muons_sta_phi[b_muons_n] = b_muons_sta_ptErr[b_muons_n] = -1e5;

    b_muons_hits_n[b_muons_n] = b_muons_muonHits_n[b_muons_n] = 0;
    b_muons_DTHits_n[b_muons_n] = b_muons_CSCHits_n[b_muons_n] = 0;
    b_muons_RPCHits_n[b_muons_n] = 0;
    b_muons_GEMHits_n[b_muons_n] = b_muons_ME0Hits_n[b_muons_n] = 0;

    b_muons_glb_hits_n[b_muons_n] = b_muons_glb_muonHits_n[b_muons_n] = 0;
    b_muons_glb_DTHits_n[b_muons_n] = b_muons_glb_CSCHits_n[b_muons_n] = 0;
    b_muons_glb_RPCHits_n[b_muons_n] = b_muons_glb_iRPCHits_n[b_muons_n] = 0;
    b_muons_glb_GEMHits_n[b_muons_n] = b_muons_glb_ME0Hits_n[b_muons_n] = 0;

    b_muons_sta_hits_n[b_muons_n] = b_muons_sta_muonHits_n[b_muons_n] = 0;
    b_muons_sta_DTHits_n[b_muons_n] = b_muons_sta_CSCHits_n[b_muons_n] = 0;
    b_muons_sta_RPCHits_n[b_muons_n] = b_muons_sta_iRPCHits_n[b_muons_n] = 0;
    b_muons_sta_GEMHits_n[b_muons_n] = b_muons_sta_ME0Hits_n[b_muons_n] = 0;

    b_muons_gen_pt[b_muons_n] = simTP->pt();
    b_muons_gen_eta[b_muons_n] = simTP->eta();
    b_muons_gen_phi[b_muons_n] = simTP->phi();

    if ( muonRef.isNonnull() ) {
      if ( true ) {
        b_muons_isSta[b_muons_n] = muonRef->isStandAloneMuon();
        b_muons_isGlb[b_muons_n] = muonRef->isGlobalMuon();
        b_muons_isTrk[b_muons_n] = muonRef->isTrackerMuon();
        b_muons_isRPC[b_muons_n] = muonRef->isRPCMuon();
        b_muons_isLoose[b_muons_n] = muon::isLooseMuon(*muonRef);
        b_muons_isTight[b_muons_n] = pv ? muon::isTightMuon(*muonRef, *pv) : false;

        auto track = muonRef->muonBestTrack();
        b_muons_pt[b_muons_n] = track->pt();
        b_muons_eta[b_muons_n] = track->eta();
        b_muons_phi[b_muons_n] = track->phi();

        b_muons_time[b_muons_n] = muonRef->time().timeAtIpInOut;
        b_muons_rpcTime[b_muons_n] = muonRef->rpcTime().timeAtIpInOut;

        auto hitPattern = track->hitPattern();
        b_muons_hits_n[b_muons_n] = hitPattern.numberOfValidHits();
        b_muons_muonHits_n[b_muons_n] = hitPattern.numberOfValidMuonHits();
        b_muons_DTHits_n[b_muons_n] = hitPattern.numberOfValidMuonDTHits();
        b_muons_CSCHits_n[b_muons_n] = hitPattern.numberOfValidMuonCSCHits();
        b_muons_RPCHits_n[b_muons_n] = hitPattern.numberOfValidMuonRPCHits();
        b_muons_GEMHits_n[b_muons_n] = hitPattern.numberOfValidMuonGEMHits();
        b_muons_ME0Hits_n[b_muons_n] = hitPattern.numberOfValidMuonME0Hits();
      }

      if ( muonRef->isGlobalMuon() ) {
        auto track = muonRef->combinedMuon();
        b_muons_glb_pt[b_muons_n] = track->pt();
        b_muons_glb_eta[b_muons_n] = track->eta();
        b_muons_glb_phi[b_muons_n] = track->phi();
        b_muons_glb_ptErr[b_muons_n] = track->ptError();

        std::vector<RPCHitInfo> rpcHits = collectRPCHits(rpcGeom.product(), track, rpcHitsHandle.product());
        const int nIRPCHits = std::accumulate(rpcHits.begin(), rpcHits.end(), 0, [](int n, const RPCHitInfo h){return h.isIRPC ? n+1 : n;});
        const double sumTime = std::accumulate(rpcHits.begin(), rpcHits.end(), 0.0,
                                               [](double t, const RPCHitInfo h){return h.t == 0.0 ? t : t+h.t;}); // to be removed for >= 911p1
        b_muons_rpcTime1[b_muons_n] = rpcHits.empty() ? 0 : sumTime/rpcHits.size();

        auto hitPattern = track->hitPattern();
        b_muons_glb_hits_n[b_muons_n] = hitPattern.numberOfValidHits();
        b_muons_glb_muonHits_n[b_muons_n] = hitPattern.numberOfValidMuonHits();
        b_muons_glb_DTHits_n[b_muons_n] = hitPattern.numberOfValidMuonDTHits();
        b_muons_glb_CSCHits_n[b_muons_n] = hitPattern.numberOfValidMuonCSCHits();
        b_muons_glb_RPCHits_n[b_muons_n] = hitPattern.numberOfValidMuonRPCHits();
        b_muons_glb_iRPCHits_n[b_muons_n] = nIRPCHits;
        b_muons_glb_GEMHits_n[b_muons_n] = hitPattern.numberOfValidMuonGEMHits();
        b_muons_glb_ME0Hits_n[b_muons_n] = hitPattern.numberOfValidMuonME0Hits();
      }

      if ( muonRef->isStandAloneMuon() ) {
        auto track = muonRef->standAloneMuon();
        b_muons_sta_pt[b_muons_n] = track->pt();
        b_muons_sta_eta[b_muons_n] = track->eta();
        b_muons_sta_phi[b_muons_n] = track->phi();
        b_muons_sta_ptErr[b_muons_n] = track->ptError();

        std::vector<RPCHitInfo> rpcHits = collectRPCHits(rpcGeom.product(), track, rpcHitsHandle.product());
        const int nIRPCHits = std::accumulate(rpcHits.begin(), rpcHits.end(), 0, [](int n, const RPCHitInfo h){return h.isIRPC ? n+1 : n;});
        if ( b_muons_rpcTime1[b_muons_n] == -1e5 ) {
          const double sumTime = std::accumulate(rpcHits.begin(), rpcHits.end(), 0.0,
              [](double t, const RPCHitInfo h){return h.t == 0.0 ? t : t+h.t;}); // to be removed for >= 911p1
          b_muons_rpcTime1[b_muons_n] = rpcHits.empty() ? 0 : sumTime/rpcHits.size();
        }

        auto hitPattern = track->hitPattern();
        b_muons_sta_hits_n[b_muons_n] = hitPattern.numberOfValidHits();
        b_muons_sta_muonHits_n[b_muons_n] = hitPattern.numberOfValidMuonHits();
        b_muons_sta_DTHits_n[b_muons_n] = hitPattern.numberOfValidMuonDTHits();
        b_muons_sta_CSCHits_n[b_muons_n] = hitPattern.numberOfValidMuonCSCHits();
        b_muons_sta_RPCHits_n[b_muons_n] = hitPattern.numberOfValidMuonRPCHits();
        b_muons_sta_iRPCHits_n[b_muons_n] = nIRPCHits;
        b_muons_sta_GEMHits_n[b_muons_n] = hitPattern.numberOfValidMuonGEMHits();
        b_muons_sta_ME0Hits_n[b_muons_n] = hitPattern.numberOfValidMuonME0Hits();
      }
    }

    ++b_muons_n;
  }

  // For the fake muons
  for ( size_t i=0, n=muonsHandle->size(); i<n; ++i ) {
    edm::RefToBase<reco::Muon> muonRef(muonsHandle, i);
    edm::Ref<TrackingParticleCollection> simTP;
    if ( glbToSimColl.find(muonRef) != glbToSimColl.end() ) {
      auto matches = glbToSimColl[muonRef];
      if ( !matches.empty() ) simTP = matches.begin()->first;
    }
    else if ( staToSimColl.find(muonRef) != staToSimColl.end() ) {
      auto matches = staToSimColl[muonRef];
      if ( !matches.empty() ) simTP = matches.begin()->first;
    }
    else if ( trkToSimColl.find(muonRef) != trkToSimColl.end() ) {
      auto matches = trkToSimColl[muonRef];
      if ( !matches.empty() ) simTP = matches.begin()->first;
    }
    if ( simTP.isNonnull() ) continue;

    b_muons_isGen[b_muons_n] = false;
    b_muons_isSta[b_muons_n] = muonRef->isStandAloneMuon();
    b_muons_isGlb[b_muons_n] = muonRef->isGlobalMuon();
    b_muons_isTrk[b_muons_n] = muonRef->isTrackerMuon();
    b_muons_isRPC[b_muons_n] = muonRef->isRPCMuon();
    b_muons_isLoose[b_muons_n] = muon::isLooseMuon(*muonRef);
    b_muons_isTight[b_muons_n] = pv ? muon::isTightMuon(*muonRef, *pv) : false;
    b_muons_gen_pt[b_muons_n] = b_muons_gen_eta[b_muons_n] = b_muons_gen_phi[b_muons_n] = -1e5;
    b_muons_time[b_muons_n] = b_muons_rpcTime[b_muons_n] = b_muons_rpcTime1[b_muons_n] = -1e5;

    b_muons_glb_pt[b_muons_n] = b_muons_glb_eta[b_muons_n] = b_muons_glb_phi[b_muons_n] = b_muons_glb_ptErr[b_muons_n] = 0;
    b_muons_sta_pt[b_muons_n] = b_muons_sta_eta[b_muons_n] = b_muons_sta_phi[b_muons_n] = b_muons_sta_ptErr[b_muons_n] = 0;

    b_muons_glb_hits_n[b_muons_n] = b_muons_glb_muonHits_n[b_muons_n] = 0;
    b_muons_glb_DTHits_n[b_muons_n] = b_muons_glb_CSCHits_n[b_muons_n] = 0;
    b_muons_glb_RPCHits_n[b_muons_n] = b_muons_glb_muonHits_n[b_muons_n] = 0;
    b_muons_glb_GEMHits_n[b_muons_n] = b_muons_glb_muonHits_n[b_muons_n] = 0;

    b_muons_sta_hits_n[b_muons_n] = b_muons_sta_muonHits_n[b_muons_n] = 0;
    b_muons_sta_DTHits_n[b_muons_n] = b_muons_sta_CSCHits_n[b_muons_n] = 0;
    b_muons_sta_RPCHits_n[b_muons_n] = b_muons_sta_muonHits_n[b_muons_n] = 0;
    b_muons_sta_GEMHits_n[b_muons_n] = b_muons_sta_muonHits_n[b_muons_n] = 0;

    if ( true ) {
      auto track = muonRef->muonBestTrack();
      b_muons_pt[b_muons_n] = track->pt();
      b_muons_eta[b_muons_n] = track->eta();
      b_muons_phi[b_muons_n] = track->phi();

      b_muons_time[b_muons_n] = muonRef->time().timeAtIpInOut;
      b_muons_rpcTime[b_muons_n] = muonRef->rpcTime().timeAtIpInOut;

      auto hitPattern = track->hitPattern();
      b_muons_hits_n[b_muons_n] = hitPattern.numberOfValidHits();
      b_muons_muonHits_n[b_muons_n] = hitPattern.numberOfValidMuonHits();
      b_muons_DTHits_n[b_muons_n] = hitPattern.numberOfValidMuonDTHits();
      b_muons_CSCHits_n[b_muons_n] = hitPattern.numberOfValidMuonCSCHits();
      b_muons_RPCHits_n[b_muons_n] = hitPattern.numberOfValidMuonRPCHits();
      b_muons_GEMHits_n[b_muons_n] = hitPattern.numberOfValidMuonGEMHits();
      b_muons_ME0Hits_n[b_muons_n] = hitPattern.numberOfValidMuonME0Hits();
    }
 
    if ( muonRef->isGlobalMuon() ) {
      auto track = muonRef->combinedMuon();
      b_muons_glb_pt[b_muons_n] = track->pt();
      b_muons_glb_eta[b_muons_n] = track->eta();
      b_muons_glb_phi[b_muons_n] = track->phi();
      b_muons_glb_ptErr[b_muons_n] = track->ptError();

      std::vector<RPCHitInfo> rpcHits = collectRPCHits(rpcGeom.product(), track, rpcHitsHandle.product());
      const int nIRPCHits = std::accumulate(rpcHits.begin(), rpcHits.end(), 0, [](int n, const RPCHitInfo h){return h.isIRPC ? n+1 : n;});
      const double sumTime = std::accumulate(rpcHits.begin(), rpcHits.end(), 0.0,
                                             [](double t, const RPCHitInfo h){return h.t == 0.0 ? t : t+h.t;}); // to be removed for >= 911p1
      b_muons_rpcTime1[b_muons_n] = rpcHits.empty() ? 0 : sumTime/rpcHits.size();

      auto hitPattern = track->hitPattern();
      b_muons_glb_hits_n[b_muons_n] = hitPattern.numberOfValidHits();
      b_muons_glb_muonHits_n[b_muons_n] = hitPattern.numberOfValidMuonHits();
      b_muons_glb_DTHits_n[b_muons_n] = hitPattern.numberOfValidMuonDTHits();
      b_muons_glb_CSCHits_n[b_muons_n] = hitPattern.numberOfValidMuonCSCHits();
      b_muons_glb_RPCHits_n[b_muons_n] = hitPattern.numberOfValidMuonRPCHits();
      b_muons_glb_iRPCHits_n[b_muons_n] = nIRPCHits;
      b_muons_glb_GEMHits_n[b_muons_n] = hitPattern.numberOfValidMuonGEMHits();
      b_muons_glb_ME0Hits_n[b_muons_n] = hitPattern.numberOfValidMuonME0Hits();
    }

    if ( muonRef->isStandAloneMuon() ) {
      auto track = muonRef->standAloneMuon();
      b_muons_sta_pt[b_muons_n] = track->pt();
      b_muons_sta_eta[b_muons_n] = track->eta();
      b_muons_sta_phi[b_muons_n] = track->phi();
      b_muons_sta_ptErr[b_muons_n] = track->ptError();

      std::vector<RPCHitInfo> rpcHits = collectRPCHits(rpcGeom.product(), track, rpcHitsHandle.product());
      const int nIRPCHits = std::accumulate(rpcHits.begin(), rpcHits.end(), 0, [](int n, const RPCHitInfo h){return h.isIRPC ? n+1 : n;});
      if ( b_muons_rpcTime1[b_muons_n] == -1e5 ) {
        const double sumTime = std::accumulate(rpcHits.begin(), rpcHits.end(), 0.0,
            [](double t, const RPCHitInfo h){return h.t == 0.0 ? t : t+h.t;}); // to be removed for >= 911p1
        b_muons_rpcTime1[b_muons_n] = rpcHits.empty() ? 0 : sumTime/rpcHits.size();
      }

      auto hitPattern = track->hitPattern();
      b_muons_sta_hits_n[b_muons_n] = hitPattern.numberOfValidHits();
      b_muons_sta_muonHits_n[b_muons_n] = hitPattern.numberOfValidMuonHits();
      b_muons_sta_DTHits_n[b_muons_n] = hitPattern.numberOfValidMuonDTHits();
      b_muons_sta_CSCHits_n[b_muons_n] = hitPattern.numberOfValidMuonCSCHits();
      b_muons_sta_RPCHits_n[b_muons_n] = hitPattern.numberOfValidMuonRPCHits();
      b_muons_sta_iRPCHits_n[b_muons_n] = nIRPCHits;
      b_muons_sta_GEMHits_n[b_muons_n] = hitPattern.numberOfValidMuonGEMHits();
      b_muons_sta_ME0Hits_n[b_muons_n] = hitPattern.numberOfValidMuonME0Hits();
    }

    ++b_muons_n;
  }

  tree_->Fill();
}

std::vector<RPCHitInfo> MuonPerformanceAnalyzer::collectRPCHits(const RPCGeometry* rpcGeom, const reco::TrackRef track, 
                                                                const RPCRecHitCollection* rpcHits) const
{
  std::vector<RPCHitInfo> out;
  if ( track.isNull() ) return out;  

  for ( auto ihit = track->recHitsBegin(); ihit != track->recHitsEnd(); ++ihit ) {
    if ( !(*ihit)->isValid() ) continue;
    if ( (*ihit)->geographicalId().det() != DetId::Muon or (*ihit)->geographicalId().subdetId() != 3 ) continue;

    const RPCDetId detId((*ihit)->rawId());
    const auto roll = rpcGeom->roll(detId);
    //const double r0 = roll->toGlobal(LocalPoint(0,0)).mag();

    const auto hitRange = rpcHits->get(detId);
    for ( auto rpcHit = hitRange.first; rpcHit != hitRange.second; ++rpcHit ) {
      const auto lp = rpcHit->localPosition();
      const auto gp = roll->toGlobal(lp);

      // redo the position without y-position (effective for iRPC)
      const float fstrip = (roll->centreOfStrip(rpcHit->firstClusterStrip())).x();
      const float lstrip = (roll->centreOfStrip(rpcHit->firstClusterStrip()+rpcHit->clusterSize())).x();
      const float x0 = (fstrip + lstrip)/2;

      //out.push_back(RPCHitInfo(lp.x(), lp.y(), rpcHit->time(), gp.x(), gp.y(), gp.z(), x0, roll->isIRPC()));
      out.emplace_back(lp.x(), lp.y(), rpcHit->time(), gp.x(), gp.y(), gp.z(), x0, roll->isIRPC());
    }
  }

  return out;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonPerformanceAnalyzer);


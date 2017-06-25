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
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLinkfwd.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4D.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegment.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class MuonRPCAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  MuonRPCAnalyzer(const edm::ParameterSet& pset);
  ~MuonRPCAnalyzer() {};

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  //edm::EDGetTokenT<edm::SimVertexContainer> simVertexToken_;
  //edm::EDGetTokenT<RPCDigiCollection> rpcDigisToken_;
  edm::EDGetTokenT<edm::DetSetVector<RPCDigiSimLink>> rpcSimDigisToken_;

  edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;
  //edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentsToken_;
  //edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsToken_;
  //edm::EDGetTokenT<GEMSegmentCollection> gemSegmentsToken_;

  edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;

public:
  TTree* tree_;

  const static unsigned int rpcDet_N = 5000;
  unsigned int b_rpcDet_n;
  bool b_rpcDet_isBarrel[rpcDet_N], b_rpcDet_isIRPC[rpcDet_N];
  short b_rpcDet_wheel[rpcDet_N], b_rpcDet_disk[rpcDet_N];
  unsigned short b_rpcDet_station[rpcDet_N];
  short b_rpcDet_ring[rpcDet_N];
  unsigned short b_rpcDet_sector[rpcDet_N], b_rpcDet_layer[rpcDet_N], b_rpcDet_roll[rpcDet_N];
  unsigned short b_rpcDet_nHit[rpcDet_N], b_rpcDet_nSimHit[rpcDet_N];
  unsigned short b_rpcDet_nOverlap[rpcDet_N];

};

MuonRPCAnalyzer::MuonRPCAnalyzer(const edm::ParameterSet& pset):
  genParticlesToken_(consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genParticles"))),
  //simVertexToken_(consumes<edm::SimVertexContainer>(pset.getParameter<edm::InputTag>("simVertex"))),
  //rpcSimHitsToken_(consumes<edm::PSimHitContainer>(pset.getParameter<edm::InputTag>("rpcSimHits"))),
  //rpcDigisToken_(consumes<RPCDigiCollection>(pset.getParameter<edm::InputTag>("rpcDigis"))),
  rpcSimDigisToken_(consumes<edm::DetSetVector<RPCDigiSimLink>>(pset.getParameter<edm::InputTag>("rpcSimDigis"))),
  rpcRecHitsToken_(consumes<RPCRecHitCollection>(pset.getParameter<edm::InputTag>("rpcRecHits"))),
  //dtSegmentsToken_(consumes<DTRecSegment4DCollection>(pset.getParameter<edm::InputTag>("dtSegments"))),
  //cscSegmentsToken_(consumes<CSCSegmentCollection>(pset.getParameter<edm::InputTag>("cscSegments"))),
  //gemSegmentsToken_(consumes<GEMSegmentCollection>(pset.getParameter<edm::InputTag>("gemSegments"))),
  muonsToken_(consumes<reco::MuonCollection>(pset.getParameter<edm::InputTag>("muons"))),
  verticesToken_(consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertices")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("tree", "tree");
  tree_->Branch("rpcDet_n", &b_rpcDet_n, "rpcDet_n/s");
  tree_->Branch("rpcDet_isBarrel", b_rpcDet_isBarrel, "rpcDet_isBarrel[rpcDet_n]/O");
  tree_->Branch("rpcDet_isIRPC"  , b_rpcDet_isIRPC  , "rpcDet_isIRPC[rpcDet_n]/O"  );
  tree_->Branch("rpcDet_wheel"   , b_rpcDet_wheel   , "rpcDet_wheel[rpcDet_n]/S"   );
  tree_->Branch("rpcDet_disk"    , b_rpcDet_disk    , "rpcDet_disk[rpcDet_n]/S"    );
  tree_->Branch("rpcDet_station" , b_rpcDet_station , "rpcDet_station[rpcDet_n]/s" );
  tree_->Branch("rpcDet_ring"    , b_rpcDet_ring    , "rpcDet_ring[rpcDet_n]/s"    );
  tree_->Branch("rpcDet_sector"  , b_rpcDet_sector  , "rpcDet_sector[rpcDet_n]/s"  );
  tree_->Branch("rpcDet_layer"   , b_rpcDet_layer   , "rpcDet_layer[rpcDet_n]/s"   );
  tree_->Branch("rpcDet_roll"    , b_rpcDet_roll    , "rpcDet_roll[rpcDet_n]/s"    );
  tree_->Branch("rpcDet_nHit"    , b_rpcDet_nHit    , "rpcDet_nHit[rpcDet_n]/s"    );
  tree_->Branch("rpcDet_nSimHit" , b_rpcDet_nSimHit , "rpcDet_nSimHit[rpcDet_n]/s" );
  tree_->Branch("rpcDet_nOverlap", b_rpcDet_nOverlap, "rpcDet_nOverlap[rpcDet_n]/s");

}

void MuonRPCAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  b_rpcDet_n = 0;

  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  event.getByToken(genParticlesToken_, genParticlesHandle);

  edm::Handle<RPCRecHitCollection> rpcRecHitsHandle;
  event.getByToken(rpcRecHitsToken_, rpcRecHitsHandle);

  edm::Handle<reco::VertexCollection> verticesHandle;
  event.getByToken(verticesToken_, verticesHandle);

  edm::Handle<edm::DetSetVector<RPCDigiSimLink>> rpcSimDigisHandle;
  event.getByToken(rpcSimDigisToken_, rpcSimDigisHandle);

  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  std::map<unsigned int, std::vector<RPCDigiSimLink>> detTrackToSimDigiMap;
  for ( auto detSet : *rpcSimDigisHandle ) {
    const auto detId = detSet.detId();
    if ( detTrackToSimDigiMap.find(detId) == detTrackToSimDigiMap.end() ) {
      detTrackToSimDigiMap[detId] = std::vector<RPCDigiSimLink>();
    }
    for ( auto idigi = detSet.begin(); idigi != detSet.end(); ++idigi ) {
      detTrackToSimDigiMap[detId].emplace_back(*idigi);
    }
  }

  for ( const auto roll : rpcGeom->rolls() ) {
    const RPCDetId rpcId = roll->id();

    b_rpcDet_isBarrel[b_rpcDet_n] = (rpcId.region() == 0);
    b_rpcDet_isIRPC[b_rpcDet_n] = roll->isIRPC();
    if ( b_rpcDet_isBarrel[b_rpcDet_n] ) {
      b_rpcDet_wheel[b_rpcDet_n] = rpcId.ring();
      b_rpcDet_station[b_rpcDet_n] = rpcId.station();
      b_rpcDet_disk[b_rpcDet_n] = b_rpcDet_ring[b_rpcDet_n] = 0;
    }
    else {
      b_rpcDet_wheel[b_rpcDet_n] = b_rpcDet_station[b_rpcDet_n] = 0;
      b_rpcDet_disk[b_rpcDet_n] = rpcId.region()*rpcId.station();
      b_rpcDet_ring[b_rpcDet_n] = rpcId.ring();
    }
    b_rpcDet_sector[b_rpcDet_n] = rpcId.sector();
    b_rpcDet_layer[b_rpcDet_n] = rpcId.layer();
    b_rpcDet_roll[b_rpcDet_n] = rpcId.roll();
    
    const auto hitRange = rpcRecHitsHandle->get(roll->id());
    const unsigned int nRecHit = hitRange.second-hitRange.first;
    b_rpcDet_nHit[b_rpcDet_n] = nRecHit;
    //for ( auto ihit = hitRange.first; ihit != hitRange.second; ++ihit ) {
    //  hCLs_->Fill(ihit->clusterSize());
    //  hXerr_->Fill(sqrt(max(0.F, ihit->localPositionError().xx())));
    //  hYerr_->Fill(sqrt(max(0.F, ihit->localPositionError().yy())));
    //}

    b_rpcDet_nSimHit[b_rpcDet_n] = b_rpcDet_nOverlap[b_rpcDet_n] = 0;
    if ( detTrackToSimDigiMap.find(roll->id()) != detTrackToSimDigiMap.end() ) {
      std::set<unsigned int> simTrackIds;
      const auto& rpcSimDigis = detTrackToSimDigiMap[rpcId.rawId()];
      std::vector<unsigned int> stripProfile(roll->nstrips());
      for ( auto& simDigi : rpcSimDigis ) {
        const unsigned int trackId = simDigi.getTrackId();
        simTrackIds.insert(trackId);
        ++stripProfile[simDigi.getStrip()-1];
      }
      for ( auto n : stripProfile ) {
        if ( n > 1 ) ++b_rpcDet_nOverlap[b_rpcDet_n];
      }
      b_rpcDet_nSimHit[b_rpcDet_n] = simTrackIds.size();
    }

    ++b_rpcDet_n;
  }

/*
  std::vector<const reco::GenParticle*> genMuons;
  for ( const auto& p : *genParticleHandle ) {
    if ( std::abs(p.pdgId()) != 13 ) continue;
    if ( p.status() != 1 ) continue;

    const double pt = p.pt();
    const double abseta = std::abs(p.eta());
    if ( abseta < 1.8 or abseta > 2.4 ) continue;
    if ( pt < 20 ) continue;

    genMuons.push_back(&p);
  }
  std::sort(genMuons.begin(), genMuons.end(), 
            [](const reco::GenParticle* a, const reco::GenParticle* b){return a->pt()>b->pt();});
  if ( !genMuons.empty() ) {
    const reco::GenParticle* genMu = genMuons[0];
    const reco::Muon* recMu = 0;
    double minDR = 0.1;
    for ( auto& mu : *muonHandle ) {
      const double dR = deltaR(mu, *genMu);
      if ( dR > minDR ) continue;

      minDR = dR;
      recMu = &mu;
    }

    if ( recMu ) {
      const double dPtPt = (recMu->pt()-genMu->pt())/genMu->pt();
      hDR_->Fill(deltaR(*recMu, *genMu));
      hDPt_->Fill(dPtPt);
      hDPt_abseta_->Fill(std::abs(genMu->eta()), dPtPt);

      if ( recMu->isStandAloneMuon() ) {
        hDPtSta_->Fill( (recMu->standAloneMuon()->pt()-genMu->pt())/genMu->pt() );
      }
      if ( recMu->isGlobalMuon() ) {
        hDPtGlb_->Fill( (recMu->globalTrack()->pt()-genMu->pt())/genMu->pt() );
      }
    }
  }
*/

/*
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

    reco::MuonCollection muons;
    for ( auto& mu : *muonHandle ) {
      if ( std::abs(mu.eta()) >= 2.5 or mu.pt() < 20 ) continue;
      muons.push_back(mu);
    }
    std::sort(muons.begin(), muons.end(), [](const reco::Muon& a, const reco::Muon& b){return a.pt() > b.pt();});

    for ( auto& mu : muons ) {
      b_muon_pt[b_muon_n] = mu.pt();
      b_muon_eta[b_muon_n] = mu.eta();
      b_muon_phi[b_muon_n] = mu.phi();
      b_muon_q[b_muon_n] = mu.charge();
      b_muon_isLoose[b_muon_n] = muon::isLooseMuon(mu);
      b_muon_isTight[b_muon_n] = !pv ? false : muon::isTightMuon(mu, *pv);
      b_muon_isRPC[b_muon_n] = muon::isGoodMuon(mu, muon::RPCMuLoose, reco::Muon::RPCHitAndTrackArbitration);

      b_muon_time[b_muon_n] = mu.time().timeAtIpInOut;
      b_muon_RPCTime[b_muon_n] = mu.rpcTime().timeAtIpInOut;

      const std::vector<HitInfo> hitTimes = getRPCTimes(mu, *rpcRecHitsHandle, rpcGeom, true);
      const double sumTime = std::accumulate(hitTimes.begin(), hitTimes.end(), 0.0, [](const double& a, const HitInfo& b){return a+b.first[0];});
      b_muon_RPCTimeNew[b_muon_n] = hitTimes.empty() ? 0 : sumTime/hitTimes.size();
      b_muon_nRPC[b_muon_n] = hitTimes.size();
      b_muon_nIRPC[b_muon_n] = std::count_if(hitTimes.begin(), hitTimes.end(), [](const HitInfo& p){return p.second == 2;});
      const double sumBeta = std::accumulate(hitTimes.begin(), hitTimes.end(), 0.0, [&](const double& sum, const HitInfo& p){
        if ( p.second == 0 ) return sum;
        const double t = p.first[0];
        const double r = p.first[1];
        const double r0 = p.first[2];
        const double beta = r/(r0+speedOfLight_*t);
        return sum+beta;
      });
      b_muon_RPCBeta[b_muon_n] = hitTimes.empty() ? 0 : sumBeta/hitTimes.size();

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
*/

  tree_->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonRPCAnalyzer);


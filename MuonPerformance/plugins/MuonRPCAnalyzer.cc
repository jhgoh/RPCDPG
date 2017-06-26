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
  TTree* detTree_, * hitTree_;

  const static unsigned int rpcDet_N = 5000;
  unsigned int b_rpcDet_n;
  bool b_rpcDet_isBarrel[rpcDet_N], b_rpcDet_isIRPC[rpcDet_N];
  short b_rpcDet_wheel[rpcDet_N], b_rpcDet_disk[rpcDet_N];
  unsigned short b_rpcDet_station[rpcDet_N];
  short b_rpcDet_ring[rpcDet_N];
  unsigned short b_rpcDet_sector[rpcDet_N], b_rpcDet_layer[rpcDet_N], b_rpcDet_roll[rpcDet_N];
  unsigned short b_rpcDet_nHit[rpcDet_N], b_rpcDet_nSimHit[rpcDet_N];
  unsigned short b_rpcDet_nOverlap[rpcDet_N];

  const static unsigned int rpcHit_N = 1000;
  unsigned int b_rpcHit_n;
  bool b_rpcHit_isBarrel[rpcHit_N], b_rpcHit_isIRPC[rpcHit_N];
  short b_rpcHit_wheel[rpcHit_N], b_rpcHit_disk[rpcHit_N];
  unsigned short b_rpcHit_station[rpcHit_N];
  short b_rpcHit_ring[rpcHit_N];
  unsigned short b_rpcHit_sector[rpcHit_N], b_rpcHit_layer[rpcHit_N], b_rpcHit_roll[rpcHit_N];
  unsigned short b_rpcHit_cls[rpcHit_N], b_rpcHit_bx[rpcHit_N];
  float b_rpcHit_time[rpcHit_N];
  float b_rpcHit_lx[rpcHit_N], b_rpcHit_ly[rpcHit_N];
  float b_rpcHit_gx[rpcHit_N], b_rpcHit_gy[rpcHit_N], b_rpcHit_gz[rpcHit_N];
  float b_simHit_dx[rpcHit_N], b_simHit_dy[rpcHit_N], b_simHit_dt[rpcHit_N];

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

  detTree_ = fs->make<TTree>("rpcDet", "rpcDet");
  detTree_->Branch("n", &b_rpcDet_n, "n/s");
  detTree_->Branch("isBarrel", b_rpcDet_isBarrel, "isBarrel[n]/O");
  detTree_->Branch("isIRPC"  , b_rpcDet_isIRPC  , "isIRPC[n]/O"  );
  detTree_->Branch("wheel"   , b_rpcDet_wheel   , "wheel[n]/S"   );
  detTree_->Branch("disk"    , b_rpcDet_disk    , "disk[n]/S"    );
  detTree_->Branch("station" , b_rpcDet_station , "station[n]/s" );
  detTree_->Branch("ring"    , b_rpcDet_ring    , "ring[n]/s"    );
  detTree_->Branch("sector"  , b_rpcDet_sector  , "sector[n]/s"  );
  detTree_->Branch("layer"   , b_rpcDet_layer   , "layer[n]/s"   );
  detTree_->Branch("roll"    , b_rpcDet_roll    , "roll[n]/s"    );
  detTree_->Branch("nHit"    , b_rpcDet_nHit    , "nHit[n]/s"    );
  detTree_->Branch("nSimHit" , b_rpcDet_nSimHit , "nSimHit[n]/s" );
  detTree_->Branch("nOverlap", b_rpcDet_nOverlap, "nOverlap[n]/s");

  hitTree_ = fs->make<TTree>("rpcHit", "rpcHit");
  hitTree_->Branch("n", &b_rpcHit_n, "n/s");
  hitTree_->Branch("isBarrel", b_rpcHit_isBarrel, "isBarrel[n]/O");
  hitTree_->Branch("isIRPC"  , b_rpcHit_isIRPC  , "isIRPC[n]/O"  );
  hitTree_->Branch("wheel"   , b_rpcHit_wheel   , "wheel[n]/S"   );
  hitTree_->Branch("disk"    , b_rpcHit_disk    , "disk[n]/S"    );
  hitTree_->Branch("station" , b_rpcHit_station , "station[n]/s" );
  hitTree_->Branch("ring"    , b_rpcHit_ring    , "ring[n]/s"    );
  hitTree_->Branch("sector"  , b_rpcHit_sector  , "sector[n]/s"  );
  hitTree_->Branch("layer"   , b_rpcHit_layer   , "layer[n]/s"   );
  hitTree_->Branch("roll"    , b_rpcHit_roll    , "roll[n]/s"    );
  hitTree_->Branch("cls", b_rpcHit_cls, "cls[n]/s");
  hitTree_->Branch("bx", b_rpcHit_bx, "bx[n]/s");
  hitTree_->Branch("time", b_rpcHit_time, "time[n]/F");
  hitTree_->Branch("lx" , b_rpcHit_lx, "lx[n]/F" );
  hitTree_->Branch("ly" , b_rpcHit_ly, "ly[n]/F" );
  hitTree_->Branch("gx" , b_rpcHit_gx, "gx[n]/F" );
  hitTree_->Branch("gy" , b_rpcHit_gy, "gy[n]/F" );
  hitTree_->Branch("gz" , b_rpcHit_gz, "gz[n]/F" );
  hitTree_->Branch("dx" , b_simHit_dx, "dx[n]/F" );
  hitTree_->Branch("dy" , b_simHit_dy, "dy[n]/F" );
  hitTree_->Branch("dt" , b_simHit_dt, "dt[n]/F" );

}

void MuonRPCAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  b_rpcDet_n = b_rpcHit_n = 0;

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
    const double tof0 = roll->toGlobal(LocalPoint(0,0,0)).mag()/30;

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
    
    // Fill simDigis and find cluster first
    // Find overlaping simHit and fill them
    b_rpcDet_nSimHit[b_rpcDet_n] = b_rpcDet_nOverlap[b_rpcDet_n] = 0;
    struct HitInfo { unsigned int n; unsigned int strip1, strip2; float sumX, sumY, sumToF; };
    std::map<unsigned int, HitInfo> simTrackIdToStripRange;
    if ( detTrackToSimDigiMap.find(roll->id()) != detTrackToSimDigiMap.end() ) {
      const auto& rpcSimDigis = detTrackToSimDigiMap[rpcId.rawId()];
      std::vector<unsigned int> stripProfile(roll->nstrips());
      for ( auto& simDigi : rpcSimDigis ) {
        const unsigned int trackId = simDigi.getTrackId();
        const unsigned int strip = simDigi.getStrip();
        auto key = simTrackIdToStripRange.find(trackId);
        if ( key == simTrackIdToStripRange.end() ) {
          simTrackIdToStripRange[trackId] = {
            1, strip, strip,
            simDigi.getEntryPoint().x(), simDigi.getEntryPoint().y(), simDigi.getTimeOfFlight()
          };
        }
        else {
          ++key->second.n;
          key->second.strip1 = std::min(key->second.strip1, strip);
          key->second.strip2 = std::max(key->second.strip2, strip);
          key->second.sumX += simDigi.getEntryPoint().x();
          key->second.sumY += simDigi.getEntryPoint().y();
          key->second.sumToF += simDigi.getTimeOfFlight();
        }
        ++stripProfile[simDigi.getStrip()-1];
      }
      for ( auto n : stripProfile ) {
        if ( n > 1 ) ++b_rpcDet_nOverlap[b_rpcDet_n];
      }
      b_rpcDet_nSimHit[b_rpcDet_n] = simTrackIdToStripRange.size();
    }

    const auto hitRange = rpcRecHitsHandle->get(roll->id());
    const unsigned int nRecHit = hitRange.second-hitRange.first;
    b_rpcDet_nHit[b_rpcDet_n] = nRecHit;
    for ( auto ihit = hitRange.first; ihit != hitRange.second; ++ihit ) {
      b_rpcHit_isBarrel[b_rpcHit_n] = (rpcId.region() == 0);
      b_rpcHit_isIRPC[b_rpcHit_n] = roll->isIRPC();
      if ( b_rpcHit_isBarrel[b_rpcHit_n] ) {
        b_rpcHit_wheel[b_rpcHit_n] = rpcId.ring();
        b_rpcHit_station[b_rpcHit_n] = rpcId.station();
        b_rpcHit_disk[b_rpcHit_n] = b_rpcHit_ring[b_rpcHit_n] = 0;
      }
      else {
        b_rpcHit_wheel[b_rpcHit_n] = b_rpcHit_station[b_rpcHit_n] = 0;
        b_rpcHit_disk[b_rpcHit_n] = rpcId.region()*rpcId.station();
        b_rpcHit_ring[b_rpcHit_n] = rpcId.ring();
      }
      b_rpcHit_sector[b_rpcHit_n] = rpcId.sector();
      b_rpcHit_layer[b_rpcHit_n] = rpcId.layer();
      b_rpcHit_roll[b_rpcHit_n] = rpcId.roll();

      b_rpcHit_cls[b_rpcHit_n] = ihit->clusterSize();
      b_rpcHit_bx[b_rpcHit_n] = ihit->BunchX();
      b_rpcHit_time[b_rpcHit_n] = ihit->time();
      const auto lp = ihit->localPosition();
      b_rpcHit_lx[b_rpcHit_n] = lp.x();
      b_rpcHit_ly[b_rpcHit_n] = lp.y();
      const auto gp = roll->toGlobal(lp);
      b_rpcHit_gx[b_rpcHit_n] = gp.x();
      b_rpcHit_gy[b_rpcHit_n] = gp.y();
      b_rpcHit_gz[b_rpcHit_n] = gp.z();

      b_simHit_dx[b_rpcHit_n] = b_simHit_dy[b_rpcHit_n] = b_simHit_dt[b_rpcHit_n] = -1e5;
      for ( auto itr = simTrackIdToStripRange.begin(); itr != simTrackIdToStripRange.end(); ++itr ) {
        //const unsigned int trkId = itr->first;
        const int strip1 = itr->second.strip1;
        const int strip2 = itr->second.strip2;
        if ( strip1 > ihit->firstClusterStrip()+ihit->clusterSize() ) continue;
        if ( strip2 < ihit->firstClusterStrip() ) continue;

        b_simHit_dx[b_rpcHit_n] = lp.x() - itr->second.sumX/itr->second.n;
        b_simHit_dy[b_rpcHit_n] = lp.y() - itr->second.sumY/itr->second.n;
        b_simHit_dt[b_rpcHit_n] = ihit->time() - (itr->second.sumToF/itr->second.n - tof0);

        break;
      }
      if ( detTrackToSimDigiMap.find(roll->id()) != detTrackToSimDigiMap.end() ) {
        std::set<unsigned int> simTrackIds;
        const auto& rpcSimDigis = detTrackToSimDigiMap[rpcId.rawId()];
        for ( auto& simDigi : rpcSimDigis ) {
          const unsigned int trackId = simDigi.getTrackId();
          simTrackIds.insert(trackId);
        }
      }

      ++b_rpcHit_n;
    }

    ++b_rpcDet_n;
  }

  hitTree_->Fill();
  detTree_->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonRPCAnalyzer);


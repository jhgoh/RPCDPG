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

typedef std::pair<std::array<double, 3>, int> HitInfo;

class MuonRPCAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  MuonRPCAnalyzer(const edm::ParameterSet& pset);
  ~MuonRPCAnalyzer() {};

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  //edm::EDGetTokenT<edm::SimVertexContainer> simVertexToken_;
  edm::EDGetTokenT<RPCDigiCollection> rpcDigisToken_;
  edm::EDGetTokenT<edm::DetSetVector<RPCDigiSimLink>> rpcSimDigisToken_;

  edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;
  //edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentsToken_;
  //edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsToken_;
  //edm::EDGetTokenT<GEMSegmentCollection> gemSegmentsToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

public:
  TH1D* hStripProfile_;
  TH1D* hCLs_;
  TH1D* hXerr_, * hYerr_;
  TH1D* hDR_, * hDPt_;
  TH2D* hDPt_abseta_;

  TH1D* hDPtSta_, * hDPtGlb_;

  TH1D* hNSimHitInRoll_, * hNDigiInRoll_, * hNRecHitInRoll_;
  TH2D* hNSimHitVsNRecHitInRoll_;
};

MuonRPCAnalyzer::MuonRPCAnalyzer(const edm::ParameterSet& pset):
  genParticleToken_(consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genParticle"))),
  //simVertexToken_(consumes<edm::SimVertexContainer>(pset.getParameter<edm::InputTag>("simVertex"))),
  //rpcSimHitsToken_(consumes<edm::PSimHitContainer>(pset.getParameter<edm::InputTag>("rpcSimHits"))),
  rpcDigisToken_(consumes<RPCDigiCollection>(pset.getParameter<edm::InputTag>("rpcDigis"))),
  rpcSimDigisToken_(consumes<edm::DetSetVector<RPCDigiSimLink>>(pset.getParameter<edm::InputTag>("rpcSimDigis"))),
  rpcRecHitsToken_(consumes<RPCRecHitCollection>(pset.getParameter<edm::InputTag>("rpcRecHits"))),
  //dtSegmentsToken_(consumes<DTRecSegment4DCollection>(pset.getParameter<edm::InputTag>("dtSegments"))),
  //cscSegmentsToken_(consumes<CSCSegmentCollection>(pset.getParameter<edm::InputTag>("cscSegments"))),
  //gemSegmentsToken_(consumes<GEMSegmentCollection>(pset.getParameter<edm::InputTag>("gemSegments"))),
  muonToken_(consumes<reco::MuonCollection>(pset.getParameter<edm::InputTag>("muons"))),
  vertexToken_(consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertex")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  hStripProfile_ = fs->make<TH1D>("hStripProfile", "Strip profile;Strip;Entries", 200, 0, 200);
  hCLs_ = fs->make<TH1D>("hCLs", "Cluster size;Cluster size;Entries", 10, 0, 10);
  hXerr_ = fs->make<TH1D>("hXerr", "x error;x error (cm);Entries", 100, -10, 10);
  hYerr_ = fs->make<TH1D>("hYerr", "y error;y error (cm);Entries", 100, -10, 10);

  hDR_ = fs->make<TH1D>("hDR", "#DeltaR;#DeltaR;Entries", 100, 0, 0.5);

  hDPt_ = fs->make<TH1D>("hDPt", "#Deltap_{T};#Deltap_{T}/p_{T};Entries", 100, -0.2, 0.2);
  hDPtSta_ = fs->make<TH1D>("hDPtSta", "Standalone muon #Deltap_{T};#Deltap_{T}/p_{T};Entries", 100, -1, 1);
  hDPtGlb_ = fs->make<TH1D>("hDPtGlb", "Global muon #Deltap_{T};#Deltap_{T}/p_{T};Entries", 100, -0.2, 0.2);
  hDPt_abseta_ = fs->make<TH2D>("hDPt_abseta", "#Deltap_{T}/p_{T};|#eta|;Entries", 100, 0, 2.5, 100, -0.2, 0.2);

  hNSimHitInRoll_ = fs->make<TH1D>("hNSimHitInRoll", "SimHit multiplicity in a roll (unique trackId)", 20, 0, 20);
  hNDigiInRoll_ = fs->make<TH1D>("hNDigiInRoll", "Digi multiplicity in a roll", 20, 0, 20);
  hNRecHitInRoll_ = fs->make<TH1D>("hNRecHitInRoll", "RecHit multiplicity in a roll", 20, 0, 20);

  hNSimHitVsNRecHitInRoll_ = fs->make<TH2D>("hNSimHitVsNRecHitInRoll", "SimHit multiplicity vs RecHit multiplicity in a roll;SimHit multiplicity;RecHit multiplicity", 10, 0, 10, 10, 0, 10);
}

void MuonRPCAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  event.getByToken(genParticleToken_, genParticleHandle);

  edm::Handle<RPCRecHitCollection> rpcRecHitsHandle;
  event.getByToken(rpcRecHitsToken_, rpcRecHitsHandle);

  edm::Handle<RPCDigiCollection> rpcDigisHandle;
  event.getByToken(rpcDigisToken_, rpcDigisHandle);

  edm::Handle<reco::MuonCollection> muonHandle;
  event.getByToken(muonToken_, muonHandle);

  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByToken(vertexToken_, vertexHandle);

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
    if ( !roll->isIRPC() ) continue;

    const auto digiRange = rpcDigisHandle->get(roll->id());
    const int nDigi = digiRange.second-digiRange.first;
    for ( auto idigi = digiRange.first; idigi != digiRange.second; ++idigi ) {
      const int strip = idigi->strip();
      hStripProfile_->Fill(strip);
    }
    if ( nDigi > 0 ) hNDigiInRoll_->Fill(nDigi);

    const auto hitRange = rpcRecHitsHandle->get(roll->id());
    const int nRecHit = hitRange.second-hitRange.first;
    for ( auto ihit = hitRange.first; ihit != hitRange.second; ++ihit ) {
      hCLs_->Fill(ihit->clusterSize());
      hXerr_->Fill(sqrt(max(0.F, ihit->localPositionError().xx())));
      hYerr_->Fill(sqrt(max(0.F, ihit->localPositionError().yy())));
    }
    if ( nRecHit > 0 ) hNRecHitInRoll_->Fill(nRecHit);

    std::map<unsigned int, unsigned int> nSimDigiByTrackId;
    if ( detTrackToSimDigiMap.find(roll->id()) != detTrackToSimDigiMap.end() ) {
      const auto& rpcSimDigis = detTrackToSimDigiMap[roll->id()];
      for ( auto& simDigi : rpcSimDigis ) {
        const unsigned int trackId = simDigi.getTrackId();
        if ( nSimDigiByTrackId.find(trackId) == nSimDigiByTrackId.end() ) {
          nSimDigiByTrackId[trackId] = 0;
        }
        else {
          ++nSimDigiByTrackId[trackId];
        }
      }
    }
    const int nSimDigi = nSimDigiByTrackId.size();
    if ( nSimDigi > 0 ) hNSimHitInRoll_->Fill(nSimDigi);

    hNSimHitVsNRecHitInRoll_->Fill(nSimDigi, nRecHit);
  }

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

//  tree_->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonRPCAnalyzer);


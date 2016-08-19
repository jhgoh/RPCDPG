#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class HSCPTimeOfFlightAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  HSCPTimeOfFlightAnalyzer(const edm::ParameterSet& pset);
  ~HSCPTimeOfFlightAnalyzer() {};

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  const double speedOfLight_;
  const double signalSpeed_;

  edm::EDGetTokenT<edm::SimVertexContainer> simVertexToken_;
  edm::EDGetTokenT<edm::PSimHitContainer> rpcSimHitsToken_;
  edm::EDGetTokenT<RPCDigiCollection> rpcDigisToken_;
  edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonsToken_;

  TH1F* hSimTime_;
  TH1F* hDigiTime_;
  TH1F* hRecHitTime_;
  TH1F* hMuonHitTime_;
  TH1F* hMuonHitCorrTime_;
};

HSCPTimeOfFlightAnalyzer::HSCPTimeOfFlightAnalyzer(const edm::ParameterSet& pset):
  speedOfLight_(29.9792458), // unit in (cm/ns)
  signalSpeed_(speedOfLight_*pset.getParameter<double>("signalPropagationSpeed")),
  simVertexToken_(consumes<edm::SimVertexContainer>(pset.getParameter<edm::InputTag>("simVertex"))),
  rpcSimHitsToken_(consumes<edm::PSimHitContainer>(pset.getParameter<edm::InputTag>("rpcSimHits"))),
  rpcDigisToken_(consumes<RPCDigiCollection>(pset.getParameter<edm::InputTag>("rpcDigis"))),
  rpcRecHitsToken_(consumes<RPCRecHitCollection>(pset.getParameter<edm::InputTag>("rpcRecHits"))),
  muonsToken_(consumes<reco::MuonCollection>(pset.getParameter<edm::InputTag>("muons")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  const double nbins = 100;
  const double xmax = 2;
  hSimTime_ = fs->make<TH1F>("hSimTime", "SimTime;Time (ns)", nbins, -xmax, xmax);
  hDigiTime_ = fs->make<TH1F>("hDigiTime", "DigiTime;Time (ns)", nbins, -xmax, xmax);
  hRecHitTime_ = fs->make<TH1F>("hRecHitTime", "RecHitTime;Time (ns)", nbins, -xmax, xmax);
  hMuonHitTime_ = fs->make<TH1F>("hMuonHitTime", "MuonHitTime;Time (ns)", nbins, -xmax, xmax);
  hMuonHitCorrTime_ = fs->make<TH1F>("hMuonHitCorrTime", "MuonHitCorrTime;Time (ns)", nbins, -xmax, xmax);
}

template<typename T1, typename T2>
double dist(const T1& a, const T2& b)
{
  const double dx = a.x()-b.x();
  const double dy = a.y()-b.y();
  const double dz = a.z()-b.z();
  return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void HSCPTimeOfFlightAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::SimVertexContainer> simVertexHandle;
  event.getByToken(simVertexToken_, simVertexHandle);

  edm::Handle<edm::PSimHitContainer> rpcSimHitsHandle;
  event.getByToken(rpcSimHitsToken_, rpcSimHitsHandle);

  edm::Handle<RPCDigiCollection> rpcDigisHandle;
  event.getByToken(rpcDigisToken_, rpcDigisHandle);

  edm::Handle<RPCRecHitCollection> rpcRecHitsHandle;
  event.getByToken(rpcRecHitsToken_, rpcRecHitsHandle);

  edm::Handle<reco::MuonCollection> muonsHandle;
  event.getByToken(muonsToken_, muonsHandle);

  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  /*const math::XYZTLorentzVector simPV = [&](){
    if ( !simVertexHandle.isValid() or
         simVertexHandle->empty() ) return math::XYZTLorentzVector();
    return simVertexHandle->at(0).position();
  }();*/

  if ( rpcSimHitsHandle.isValid() ) {
    for ( auto& simHit : *rpcSimHitsHandle ) {
      const RPCDetId detId(simHit.detUnitId());
      const RPCRoll* roll = rpcGeom->roll(detId);
      if ( !roll->isIRPC() ) continue;

      const auto simHitGPos = roll->toGlobal(simHit.localPosition());
      //const double r = dist(simPV, simHitGPos);
      const double r0 = simHitGPos.mag();
      //const double beta = r/simHit.tof()/speedOfLight;

      hSimTime_->Fill(simHit.tof()-r0/speedOfLight_);
    }
  }

  if ( rpcDigisHandle.isValid() ) {
    for ( auto rpcDetDigi : *rpcDigisHandle ) {
      const RPCDetId detId = rpcDetDigi.first;
      const RPCRoll* roll = rpcGeom->roll(detId);
      if ( !roll->isIRPC() ) continue;

      for ( auto digi = rpcDetDigi.second.first;
            digi != rpcDetDigi.second.second; ++digi ) {
        hDigiTime_->Fill(digi->time());
      }
    }
  }

  if ( rpcRecHitsHandle.isValid() ) {
    for ( auto& recHit : *rpcRecHitsHandle ) {
      const RPCDetId detId = recHit.rpcId();
      const RPCRoll* roll = rpcGeom->roll(detId);
      if ( !roll->isIRPC() ) continue;
      //if ( recHit.BunchX() != 0 ) continue;
      if ( recHit.time() == 0 ) continue;

      //const auto recHitGPos = roll->toGlobal(LocalPoint(0,0,0));
      //const double r0 = recHitGPos.mag();
      hRecHitTime_->Fill(recHit.time());
    }
  }

  if ( muonsHandle.isValid() and rpcRecHitsHandle.isValid() ) {
    for ( auto& mu : *muonsHandle ) {
      if ( !mu.isRPCMuon() ) continue;

      for ( auto& match : mu.matches() ) {
        const DetId& detId = match.id;
        if ( detId.det() != DetId::Muon or detId.subdetId() != MuonSubdetId::RPC ) continue;
        const RPCRoll* roll = rpcGeom->roll(detId);
        if ( !roll->isIRPC() ) continue;

        const auto rpcMatches = match.rpcMatches;
        if ( rpcMatches.empty() ) continue;
        const auto& rpcHitsRange = rpcRecHitsHandle->get(detId); 
        if ( rpcHitsRange.second-rpcHitsRange.first == 0 ) continue;

        auto rpcMatch = rpcMatches.begin();
        if ( true ) {
          double bestDx = 1e9;
          for ( auto m = rpcMatch; m != rpcMatches.end(); ++m ) {
            const double dx = std::abs(m->x-match.x);
            if ( dx < bestDx ) {
              rpcMatch = m;
              bestDx = dx;
            }
          }
        }

        auto rpcHit = rpcHitsRange.first;
        if ( true ) {
          double bestDx = 1e9;
          for ( auto hit = rpcHit; hit != rpcHitsRange.second; ++hit ) {
            const double dx = std::abs(hit->localPosition().x()-rpcMatch->x);
            if ( dx < bestDx ) {
              rpcHit = hit;
              bestDx = dx;
            }
          }
        }
        //if ( rpcHit->BunchX() != 0 ) continue;
        if ( rpcHit->time() == 0 ) continue;
        //if ( rpcHit->timeError() < 0 ) continue;

        const auto& bounds = roll->surface().bounds();
        const double timeCorr = [&](){
          const double x1 = bounds.width()/2;
          const double x0 = bounds.widthAtHalfLength()/2;
          const double y1 = bounds.length()/2;
          const double smax = (x1-x0)/y1;
          const double dx = std::abs(match.y)*std::abs(match.x)/(x0+smax*match.y)*smax;

          return -(match.y > 0 ? 1 : -1)*std::hypot(dx, match.y)/signalSpeed_;
        }();

        hMuonHitTime_->Fill(rpcHit->time());
        hMuonHitCorrTime_->Fill(rpcHit->time()-timeCorr);
      }
    }
  }

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HSCPTimeOfFlightAnalyzer);


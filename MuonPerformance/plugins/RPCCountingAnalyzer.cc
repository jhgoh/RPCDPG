#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"

#include "TH1D.h"
#include "TProfile.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class RPCHitCounterAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources>
{
public:
  RPCHitCounterAnalyzer(const edm::ParameterSet& pset);
  virtual ~RPCHitCounterAnalyzer() = default;

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;
  void beginRun(const edm::Run& run, const edm::EventSetup& eventSetup) override;
  void endRun(const edm::Run& run, const edm::EventSetup&) override {};

private:
  edm::EDGetTokenT<RPCRecHitCollection> rpcHitToken_;

  TProfile* hArea_;
  TH1D* hCounts_;
  TH1D* hEvents_;
  
};

RPCHitCounterAnalyzer::RPCHitCounterAnalyzer(const edm::ParameterSet& pset)
{
  rpcHitToken_ = consumes<RPCRecHitCollection>(pset.getParameter<edm::InputTag>("rpcRecHits"));
  hArea_ = nullptr;
}

void RPCHitCounterAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  if ( hArea_ != nullptr ) return;

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  // Set the roll names
  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  hArea_ = fs->make<TProfile>("hArea", "Roll area;Roll Name;Area [cm^{2}]", 5000, 1, 5001, 0, 100000);
  hCounts_ = fs->make<TH1D>("hCounts", "Counts;Roll Index;Number of RecHits", 5000, 1, 5001);
  hEvents_ = fs->make<TH1D>("hEvents", "hEvents;Types;Number of Events", 1, 1, 2);

  int i = 0;
  for ( const RPCRoll* roll : rpcGeom->rolls() ) {
    ++i;
    const auto detId = roll->id();
    const string rollName = RPCGeomServ(detId).name();
    const double width = roll->surface().bounds().width();
    const double height = roll->surface().bounds().length();

    hArea_->GetXaxis()->SetBinLabel(i, rollName.c_str());
    hArea_->Fill(i, width*height);
  }
}

void RPCHitCounterAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace std;

  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  edm::Handle<RPCRecHitCollection> rpcHitHandle;
  event.getByToken(rpcHitToken_, rpcHitHandle);

  hEvents_->Fill(1);
  for ( auto rpcHit : *rpcHitHandle ) {
    const auto detId = rpcHit.rawId();
    const string rollName = RPCGeomServ(detId).name();

    const int idx = hArea_->GetXaxis()->FindBin(rollName.c_str());
    hCounts_->Fill(idx);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RPCHitCounterAnalyzer);

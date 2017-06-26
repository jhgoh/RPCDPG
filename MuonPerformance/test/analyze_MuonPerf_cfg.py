import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process("MuonAnalyser",eras.Phase2C2)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(),)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('SimMuon.MCTruth.muonAssociatorByHitsHelper_cfi')
process.muonAssociatorByHitsHelper.usePhase2Tracker = cms.bool(True)
process.muonAssociatorByHitsHelper.useGEMs = cms.bool(True)
process.muonAssociatorByHitsHelper.pixelSimLinkSrc = cms.InputTag("simSiPixelDigis:Pixel")
process.muonAssociatorByHitsHelper.stripSimLinkSrc = cms.InputTag("simSiPixelDigis:Tracker")

from Validation.RecoMuon.selectors_cff import muonTPSet
process.muonPerf = cms.EDAnalyzer("MuonPerformanceAnalyzer",
    rpcHits = cms.InputTag("rpcRecHits"),
    muons = cms.InputTag("muons"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    muAssoc =  cms.InputTag("muonAssociatorByHitsHelper"),
    tpSelector = muonTPSet,
    simTPs = cms.InputTag("mix","MergedTrackTruth"),
)
process.muonPerf.tpSelector.maxRapidity = 3.0
process.muonPerf.tpSelector.minRapidity = -3.0
process.muonPerf.tpSelector.pdgId = cms.vint32(13,-13)

process.rpcDet = cms.EDAnalyzer("MuonRPCAnalyzer",
    rpcRecHits = cms.InputTag("rpcRecHits"),
    rpcSimDigis = cms.InputTag("simMuonRPCReDigis", "RPCDigiSimLink"),
    muons = cms.InputTag("muons"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("ntuple.root"))

process.p = cms.Path(
    process.muonAssociatorByHitsHelper + process.muonPerf
  + process.rpcDet
)

process.source.fileNames = [
    '/store/user/jhgoh/RPCUpgrade/20170623_1/step3_redigi_SingleMuPt100_n96_no2D/step3_000.root',
]

process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')

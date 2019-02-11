import FWCore.ParameterSet.Config as cms

import sys
process = cms.Process("ANA")

process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
)

process.source.fileNames.append('file:/xrootd/store/user/yekang/CRAB_PrivateMC/HGCalStainless_default_me0_singleNu_RECO/190129_065738/0000/singleNu_GEN-SIM-DIGI_98.root')

#process.rpcInfo = cms.EDAnalyzer("MuonRPCAnalyzer",
#    genParticles = cms.InputTag("genParticles"),
#    rpcDigis = cms.InputTag("simMuonRPCDigis"),
#    rpcRecHits = cms.InputTag("rpcRecHits"),
#    rpcSimDigis = cms.InputTag("simMuonRPCDigis:RPCDigiSimLink"),
#    muons = cms.InputTag("muons"),
#    vertices = cms.InputTag("offlinePrimaryVerticies"),
#)

process.rpcHit = cms.EDAnalyzer("RPCHitCounterAnalyzer",
    rpcRecHits = cms.InputTag("rpcRecHits"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)
process.p = cms.Path(
    process.rpcHit
)


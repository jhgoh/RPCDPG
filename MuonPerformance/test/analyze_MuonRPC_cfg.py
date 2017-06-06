import FWCore.ParameterSet.Config as cms

import sys
try:
    n = int(sys.argv[-1])
except:
    n = 192
#n = 192
#n = 96
#n = 32

process = cms.Process("ANA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
)

for i in range(99):
    process.source.fileNames.append('/store/user/jhgoh/RPCUpgrade/20170605_1/SingleMuPt100HighEta/step3_n%d/step3_%03d.root' % (n, i))

process.rpcInfo = cms.EDAnalyzer("MuonRPCAnalyzer",
    genParticle = cms.InputTag("genParticles"),
    rpcDigis = cms.InputTag("simMuonRPCDigis"),
    rpcRecHits = cms.InputTag("rpcRecHits"),
    rpcSimDigis = cms.InputTag("simMuonRPCDigis:RPCDigiSimLink"),
    muons = cms.InputTag("muons"),
    vertex = cms.InputTag("offlinePrimaryVerticies"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist_%d.root" % n),
)
process.p = cms.Path(
    process.rpcInfo
)


import FWCore.ParameterSet.Config as cms

HSCPTree = cms.EDAnalyzer("HSCPTreeMaker",
    genParticle = cms.InputTag("genParticles"),
    doSimHit = cms.bool(False),
    simVertex = cms.InputTag("g4SimHits"),
    rpcSimHits = cms.InputTag("g4SimHits:MuonRPCHits"),
    rpcDigis = cms.InputTag("rpcDigis"),
    rpcRecHits = cms.InputTag("rpcRecHits"),
    dtSegments = cms.InputTag("dt4DSegments"),
    cscSegments = cms.InputTag("cscSegments"),
    gemSegments = cms.InputTag("gemSegments"),
    muons = cms.InputTag("muons"),
    vertex = cms.InputTag("offlinePrimaryVertices"),
    signalPropagationSpeed = cms.double(0.66), ## To be overridden by SimMuon.RPCDigitizer.muonRPCDigis_cfi.simMuonRPCDigis, this may depend on eras
    signalPdgId = cms.uint32(13),
    #signalPdgId = cms.uint32(1000015),
    doSimDigi = cms.bool(True),
    rpcSimDigis = cms.InputTag("simMuonRPCDigis:RPCDigiSimLink"),
)

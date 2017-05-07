import FWCore.ParameterSet.Config as cms

HSCPL1TAnalyzer = cms.EDAnalyzer("HSCPL1TAnalyzer",
    genParticle = cms.InputTag("genParticles"),
    simVertex = cms.InputTag("g4SimHits"),
    rpcSimHits = cms.InputTag("g4SimHits:MuonRPCHits"),
    rpcDigis = cms.InputTag("rpcDigis"),
    rpcRecHits = cms.InputTag("rpcRecHits"),
    muons = cms.InputTag("muons"),
    signalPropagationSpeed = cms.double(0.66), ## To be overridden by SimMuon.RPCDigitizer.muonRPCDigis_cfi.simMuonRPCDigis, this may depend on eras
    #signalPdgId = cms.uint32(13),
    signalPdgId = cms.uint32(1000015),
    rpcSimDigis = cms.InputTag("simMuonRPCDigis:RPCDigiSimLink"),
)

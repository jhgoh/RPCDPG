import FWCore.ParameterSet.Config as cms

process = cms.Process("COPY")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("PhysicsTools.PatAlgos.slimming.prunedGenParticles_cfi")
process.load("PhysicsTools.PatAlgos.slimming.primaryVertexAssociation_cfi")
process.load("PhysicsTools.PatAlgos.slimming.offlineSlimmedPrimaryVertices_cfi")
process.load("PhysicsTools.PatAlgos.slimming.slimmedAddPileupInfo_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/xrootd/store/user/jhgoh/RPCUpgrade/20180208_1/step3_RAW2DIGI_L1Reco_RECO_PU_1.root',
    ),
    secondaryFileNames = cms.untracked.vstring(
    )
)

process.copyAll = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('copy.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_TriggerResults_*_*',
        'keep *_simMuonRPCDigis_*_*',
        'keep *_gemRecHits_*_*',
        'keep *_me0RecHits_*_*',
        'keep *_rpcRecHits_*_*',
        'keep *_dtSegments_*_*',
        'keep *_cscSegments_*_*',
        'keep *_gemSegments_*_*',
        'keep recoTrack*_globalMuons_*_*',
        'keep recoTrack*_standAloneMuons_*_*',
        'keep recoMuons_muons_*_*',
        'drop *_*_*_COPY',
        'keep *_prunedGenParticles_*_*',
        'keep *_slimmedAddPileupInfo_*_*',
        'keep *_offlineSlimmedPrimaryVertices_*_*',
    ),
)

process.offlineSlimmedPrimaryVertices.score = "offlinePrimaryVertices"

process.p = cms.Path(
    process.prunedGenParticles
  #+ process.primaryVertexAssociation * process.offlineSlimmedPrimaryVertices
  + process.offlineSlimmedPrimaryVertices
  + process.slimmedAddPileupInfo
)
process.out = cms.EndPath(process.copyAll)


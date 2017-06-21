import FWCore.ParameterSet.Config as cms

#from SimMuon.RPCDigitizer.customizeRPCDigi import customise_rpcRedigi
#from SimMuon.Configuration.customizeRPCDigi import customise_rpcRedigi # not in CMSSW yet (910)
def customise_rpcRedigi(process):
    process.load('Configuration.StandardSequences.Digi_cff')
    process.simMuonRPCReDigis = process.simMuonRPCDigis.clone()
    process.simMuonRPCReDigis.digiModelConfig = process.simMuonRPCDigis.digiModelConfig.clone(
        IRPC_time_resolution = cms.double(1.5),
        IRPC_electronics_jitter = cms.double(0.1),
        linkGateWidth = cms.double(25.0),
    )
    process.simMuonRPCReDigis.digiIRPCModelConfig = process.simMuonRPCReDigis.digiModelConfig.clone(
        IRPC_time_resolution = cms.double(1.0),
        linkGateWidth = cms.double(25.0),
        do_Y_coordinate = cms.bool(True),
    )
    process.simMuonRPCReDigis.digiModel = cms.string('RPCSimModelTiming')
    process.RandomNumberGeneratorService.simMuonRPCReDigis = cms.PSet(
        initialSeed = cms.untracked.uint32(13579),
        engineName = cms.untracked.string('TRandom3')
    )
    process.rpcRecHits.rpcDigiLabel = cms.InputTag("simMuonRPCReDigis")
    process.validationMuonRPCDigis.rpcDigiTag = cms.untracked.InputTag("simMuonRPCReDigis")
    process.reconstruction_step.replace(
        process.rpcRecHits,
        cms.Sequence(process.simMuonRPCReDigis+process.rpcRecHits)
    )
    return process

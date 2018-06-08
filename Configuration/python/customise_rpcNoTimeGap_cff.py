import FWCore.ParameterSet.Config as cms

def customise_rpcNoTimeGap(process):
    if not hasattr(process, 'simMuonRPCDigis'):
        process.load('Configuration.StandardSequences.Digi_cff')
    process.simMuonRPCDigis.digiModelConfig.linkGateWidth = cms.double(25.0)
    process.simMuonRPCDigis.digiIRPCModelConfig.linkGateWidth = cms.double(25.0)

    return process

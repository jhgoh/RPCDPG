import FWCore.ParameterSet.Config as cms

def customise_eventContent(process):
    process.setName_('RERECO')

    process.AODSIMoutput.outputCommands = [
        "drop *",
        "keep *_genParticles_*_*",
        "keep *_g4SimHits_MuonRPCHits_*",
        "keep *_simMuonRPCDigis_*_*",
        "keep *_simMuonRPCReDigis_*_*",
        "keep *_muons__*",
        "keep *_rpcRecHits_*_*",
    ]

    return process

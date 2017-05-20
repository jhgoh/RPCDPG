import FWCore.ParameterSet.Config as cms

def customise_eventContent(process):
    process.setName_('RERECO')

    process.AODSIMoutput.outputCommands += [
        "keep *_simMuonRPCDigis_*_*",
        "keep *_g4SimHits_MuonRPCHits_*",
    ]

    return process

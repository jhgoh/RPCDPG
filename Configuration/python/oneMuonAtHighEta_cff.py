import FWCore.ParameterSet.Config as cms

def customise_OneMuonAtHighEta(process):
    process.oneHSCPAtHighEtaFilter = cms.EDFilter("CandViewSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("abs(pdgId) == 13 && charge != 0 && status == 1 && 1.8 < abs(eta) < 2.5"),
        filter = cms.bool(True),
    )

    process.generation_step.replace(process.genParticles, process.genParticles+process.oneHSCPAtHighEtaFilter)
    process.simulation_step.replace(process.g4SimHits, process.oneHSCPAtHighEtaFilter+process.g4SimHits)

    return process

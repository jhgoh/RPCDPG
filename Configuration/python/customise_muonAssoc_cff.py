import FWCore.ParameterSet.Config as cms

def customise_muonAssoc(process):
    process.AODSIMoutput.outputCommands.extend([
        "keep SiPixelClusteredmNewDetSetVector_*_*_*",
        "keep Phase2TrackerCluster1DedmNewDetSetVector_*_*_*",
        "keep *_simMuonDTDigis_*_*",
        "keep *_simMuonCSCDigis_*_*",
        "keep *_simMuonRPCDigis_*_*",
        "keep *_simMuonGEMDigis_*_*",
        "keep *_mix_g4SimHits_*",
        "keep *_mix_g4SimHitsMuon*Hits_*",
        "keep *_mix_mergedTrackTruth_*",
        "keep *_simSiPixelDigis_Pixel_*",
        "keep *_simSiPixelDigis_Tracker_*",
        "keep *_g4SimHits_Muon*Hits_*",
        "keep *_g4SimHits__*",
        "keep *_cscSegments_*_*",
        "keep *_dt4DSegments_*_*",
        "keep *_dt1DRecHits_*_*",
        "keep *_globalMuons_*_*",
    ])

    return process

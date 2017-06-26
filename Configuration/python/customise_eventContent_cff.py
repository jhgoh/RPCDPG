import FWCore.ParameterSet.Config as cms

def customise_eventContent(process):
    process.setName_('reRECO')

    process.AODSIMoutput.outputCommands = [
        "drop *",
        "keep recoGenParticles_genParticles__*",
        "keep *_g4SimHits_MuonRPCHits_*",
        "keep *_simMuonRPCDigis_*_*",
        "keep *_simMuonRPCReDigis_*_*",

        "keep *_muons__*",

        'keep *_dt4DSegments_*_*', 'keep *_cscSegments_*_*', 'keep *_rpcRecHits_*_*',
        'keep *_gemRecHits_*_*', 'keep *_me0RecHits_*_*',

        'keep *_*_MergedTrackTruth_*',

        'keep recoTracks_generalTracks_*_*',
        'keep recoTrackExtras_generalTracks_*_*',
        'keep TrackingRecHitsOwned_generalTracks_*_*',

        'keep recoTracks_standAloneMuons_*_*',
        'keep recoTrackExtras_standAloneMuons_*_*',
        'keep TrackingRecHitsOwned_standAloneMuons_*_*',

        'keep recoTracks_globalMuons_*_*',
        'keep recoTrackExtras_globalMuons_*_*',
        'keep TrackingRecHitsOwned_globalMuons_*_*',
        'keep recoMuonTrackLinkss_globalMuons_*_*',

        'keep  *_offlinePrimaryVertices__*',
    ]

    return process

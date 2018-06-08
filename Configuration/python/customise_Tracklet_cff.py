import FWCore.ParameterSet.Config as cms

def customise_Tracklet(process):
    process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff") 
    process.reconstruction_step.insert(1, process.TTTracksFromTracklet)
    
    process.AODSIMoutput.outputCommands.append('keep *_TTTracks*_Level1TTTracks_*')

    return process

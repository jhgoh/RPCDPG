import FWCore.ParameterSet.Config as cms

def customise_EBeam14TeV(process):
    if not hasattr(process, 'generator'): return process
    process.generator.comEnergy = 14000

    return process

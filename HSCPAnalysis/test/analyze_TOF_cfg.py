import FWCore.ParameterSet.Config as cms
import sys, os

from Configuration.StandardSequences.Eras import eras

process = cms.Process('Analysis',eras.Phase2C1)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring(),
#    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)

process.load("RPCUpgrade.HSCPAnalysis.HSCPTreeMaker_cff")
process.HSCPTree.rpcDigis = "simMuonRPCReDigis"
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

#process.HSCPTree.signalPdgId = 13
process.HSCPTree.signalPdgId = 1000015

process.p = cms.Path(
    process.HSCPTree
)


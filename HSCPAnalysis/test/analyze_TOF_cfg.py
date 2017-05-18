import FWCore.ParameterSet.Config as cms
import sys, os

from Configuration.StandardSequences.Eras import eras

process = cms.Process('Analysis',eras.Phase2C1)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D12_cff')
process.load('Configuration.Geometry.GeometryExtended2023D12Reco_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring(),
#    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)

process.source.fileNames = [
"root://cmsxrootd.fnal.gov///store/mc/PhaseIIFall16DR82/HSCPppstau_M_1218_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/NoPU_90X_upgrade2023_realistic_v1-v2/90000/0C58D987-E5FD-E611-935C-02163E01A5AD.root",
"root://cmsxrootd.fnal.gov///store/mc/PhaseIIFall16DR82/HSCPppstau_M_1218_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/NoPU_90X_upgrade2023_realistic_v1-v2/90000/4EE8F8EF-E0FD-E611-B9DF-02163E014421.root",
"root://cmsxrootd.fnal.gov///store/mc/PhaseIIFall16DR82/HSCPppstau_M_1218_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/NoPU_90X_upgrade2023_realistic_v1-v2/90000/90D3BB3B-E6FD-E611-B189-02163E011CEE.root",
"root://cmsxrootd.fnal.gov///store/mc/PhaseIIFall16DR82/HSCPppstau_M_1218_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/NoPU_90X_upgrade2023_realistic_v1-v2/90000/9EE612A4-E7FD-E611-83F0-02163E0143D2.root",
"root://cmsxrootd.fnal.gov///store/mc/PhaseIIFall16DR82/HSCPppstau_M_1218_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/NoPU_90X_upgrade2023_realistic_v1-v2/90000/AC00BAF8-E1FD-E611-B58D-02163E01A710.root",
"root://cmsxrootd.fnal.gov///store/mc/PhaseIIFall16DR82/HSCPppstau_M_1218_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/NoPU_90X_upgrade2023_realistic_v1-v2/90000/D6372AD0-DEFD-E611-957F-02163E01A2F4.root",
"root://cmsxrootd.fnal.gov///store/mc/PhaseIIFall16DR82/HSCPppstau_M_1218_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/NoPU_90X_upgrade2023_realistic_v1-v2/90000/E4AD50CD-C3FD-E611-A9FE-02163E01A312.root",
]

process.load("RPCUpgrade.HSCPAnalysis.HSCPTreeMaker_cff")
process.HSCPTree.rpcDigis = "simMuonRPCDigis"
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(
    process.HSCPTree
)


import FWCore.ParameterSet.Config as cms
import sys, os

from Configuration.StandardSequences.Eras import eras

process = cms.Process('Analysis',eras.Phase2C1)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D1Reco_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring(),
#    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)


if len(sys.argv) > 2:
    module = sys.argv[2]
    basedir = "root://eoscms//eos/cms/store/user/jhgoh/RPCUpgrade/20160816_1/"
    process.source.fileNames.append("%s/%s/step3_000.root" % (basedir, module))
    for line in open("step2_%s.txt" % module):
        process.source.secondaryFileNames.append("root://eoscms//eos/cms/%s" % line)
else:
    module = "SingleMuPt100"
    basedir = "root://eoscms//eos/cms/store/user/jhgoh/RPCUpgrade/20160816_2/"
    process.source.fileNames.append("%s/out_reco.root" % basedir)
    process.source.secondaryFileNames.append("%s/out_digi.root" % basedir)

process.load("RPCUpgrade.HSCPAnalysis.tofAnalysis_cff")
process.tofAnalysis.rpcDigis = "simMuonRPCDigis"
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist_%s.root" % module),
)

process.p = cms.Path(
    process.tofAnalysis
)


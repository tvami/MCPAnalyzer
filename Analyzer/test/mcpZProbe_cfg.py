# Z->mumu MIP reference on RECO / RAW-RECO (e.g. /Muon0/Run2024*-ZMu-*/RAW-RECO)
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run3_2024_cff import Run3_2024

options = VarParsing('analysis')
options._tagOrder = []  # no "_numEventN" suffix on outputFile
options.maxEvents = -1
options.register('isData', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Run on data")
options.register('gtag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "Global Tag (default: MC or data GT)")
options.register('lumiMask', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "Golden JSON (data only)")
options.register('outputEvery', 1000, VarParsing.multiplicity.singleton, VarParsing.varType.int, "")
options.parseArguments()
if not options.gtag:
    options.gtag = '150X_dataRun3_v2' if options.isData else '150X_mcRun3_2024_realistic_v2'
if not options.outputFile or options.outputFile == 'output.root':
    options.outputFile = 'mcp_zprobe.root'

process = cms.Process("MCPZ", Run3_2024)
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')
process.load('RecoLocalTracker.SiPixelRecHits.SiPixelTemplateStoreESProducer_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = options.outputEvery

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.gtag, '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
if options.isData and options.lumiMask:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename=options.lumiMask).getVLuminosityBlockRange()

process.MCPAnalyzer = cms.EDAnalyzer("MCPZProbeAnalyzer")  # defaults from fillDescriptions; label matches MCPAnalyzer
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))
process.p = cms.Path(process.MCPAnalyzer)

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run3_2024_cff import Run3_2024

options = VarParsing('analysis')
options._tagOrder = []  # no "_numEventN" suffix on outputFile
options.maxEvents = -1
options.register('isData', False,
                 VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Run on data")
options.register('gtag', '',
                 VarParsing.multiplicity.singleton, VarParsing.varType.string, "Global Tag (default: MC or data GT)")
options.register('lumiMask', '',
                 VarParsing.multiplicity.singleton, VarParsing.varType.string, "Golden JSON (data only)")
options.register('saveTrigNames', True,
                 VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Store all HLT names/decisions")
options.register('pixelCPE', 'PixelCPETemplateReco',
                 VarParsing.multiplicity.singleton, VarParsing.varType.string, "Pixel CPE")
options.register('outputEvery', 100,
                 VarParsing.multiplicity.singleton, VarParsing.varType.int, "")
options.register('tpOnly', False,
                 VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Keep only Z tag-and-probe tracks")
options.register('requireDeDx', False,
                 VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Store only tracks with DeDxHitInfo")
options.register('requireTrig', False,
                 VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Keep only JetMET or Tau triggered events")
options.parseArguments()
if not options.gtag:
    options.gtag = '150X_dataRun3_v2' if options.isData else '150X_mcRun3_2024_realistic_v2'
if not options.outputFile or options.outputFile == 'output.root':
    options.outputFile = 'mcp_ntuple.root'

process = cms.Process("MCP", Run3_2024)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# Provides the PixelCPETemplateReco ESProducer (+ its SiPixelTemplateStore) used by probQ
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
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))

process.MCPAnalyzer = cms.EDAnalyzer(
    "MCPAnalyzer",
    slimmedMET = cms.InputTag("slimmedMETs"),
    slimmedPuppiMET = cms.InputTag("slimmedMETsPuppi"),
    isolatedTracks = cms.InputTag("isolatedTracks"),
    dedxHitInfo=cms.InputTag("isolatedTracks"),
    prunedGenParticles=cms.InputTag("prunedGenParticles"),
    primaryVertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    pixelCPE=cms.string(options.pixelCPE),
    mcpPdgId=cms.int32(10000200),
    triggerResults=cms.InputTag("TriggerResults", "", "HLT"),
    muons=cms.InputTag("slimmedMuons"),
    triggerObjects=cms.InputTag("slimmedPatTrigger"),
    saveTrigNames=cms.bool(options.saveTrigNames),
    tpOnly=cms.bool(options.tpOnly),
    requireTrig=cms.bool(options.requireTrig),
    requireDeDx=cms.bool(options.requireDeDx),
)

process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))

process.p = cms.Path(process.MCPAnalyzer)

import FWCore.ParameterSet.Config as cms

process = cms.Process('TestPuppi')
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.GlobalTag.globaltag = 'START53_V7G::All'

process.load('Dummy/Puppi/Puppi_cff')   

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.source = cms.Source("PoolSource",
  fileNames  = cms.untracked.vstring('file:/cms/home/brochero/IsolationPUPPITest/TaeTest/CMSSW_6_2_0_SLHC23_patch1/src/Dummy/Puppi/001B19B3-EBFA-E411-9E69-0025905B8606.root')
  #fileNames  = cms.untracked.vstring('file:/cms/home/tjkim/work/phase2/CMSSW_6_2_0_SLHC23_patch1/src/Dummy/Puppi/test/DYToMuMu-Ph2_003E69A0-89A7-E411-82FC-0025905B8606.root')
  #fileNames  = cms.untracked.vstring('file:/tmp/DYToMuMu-Ph2_003E69A0-89A7-E411-82FC-0025905B8606.root')
  #fileNames  = cms.untracked.vstring('root://cmsxrootd-site.fnal.gov//store/results/top/StoreResults/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/USER/Spring14dr_PU_S14_POSTLS170_V6AN1_miniAOD706p1_814812ec83fce2f620905d2bb30e9100-v2/00000/0012F41F-FA17-E411-A1FF-0025905A48B2.root')
)
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
)

# -- PUPPI
from pfPUPPISequence_cff import *
load_pfPUPPI_sequence(process, 'pfPUPPISequence', algo = 'PUPPI',
  src_puppi = 'pfAllHadronsAndPhotonsForPUPPI',
  cone_puppi_central = 0.5
)

# load user-defined muon PF-isolation values
muon_src, cone_size = 'muons', 0.4
from MuonPFIsolationSequence_cff import *

load_muonPFiso_sequence(process, 'MuonPFIsoSequencePUPPI', algo = 'R04PUPPI',
  src = muon_src,
  src_charged_hadron = 'pfPUPPIChargedHadrons',
  src_neutral_hadron = 'pfPUPPINeutralHadrons',
  src_photon         = 'pfPUPPIPhotons',
  coneR = cone_size
) 

process.TFileService = cms.Service("TFileService", fileName = cms.string("histo.root") )

process.pfIsoReader = cms.EDAnalyzer("PFIsoMuon",
                                     Muons = cms.InputTag('muons'),
                                     PrintMuons = cms.bool(True),
                                     IsoDepMuon = cms.VInputTag(cms.InputTag('muPFIsoDepositCHR04PUPPI'),
                                                                cms.InputTag('muPFIsoDepositPhR04PUPPI'),
                                                                cms.InputTag('muPFIsoDepositNHR04PUPPI')),
                                     IsoValMuon = cms.VInputTag(cms.InputTag('muPFIsoValueCHR04PUPPI'),
                                                                cms.InputTag('muPFIsoValuePhR04PUPPI'),
                                                                cms.InputTag('muPFIsoValueNHR04PUPPI')),
)

#process.puppiSequence = cms.Sequence(process.puppi)
#process.p = cms.Path(process.puppiSequence)
process.p = cms.Path(  process.pfPUPPISequence*
                       process.MuonPFIsoSequencePUPPI*
                       process.pfIsoReader
                    )

process.output = cms.OutputModule("PoolOutputModule", 
                                  outputCommands = cms.untracked.vstring('drop *',
                                                                         'keep *_muons_*_*',
                                                                         'keep *_generator_*_*',
                                                                         'keep *_g4SimHits_*_*',
                                                                         'keep *_genParticles_*_*',
                                                                         'keep *_addPileupInfo_*_*',
                                                                         'keep *_offlineSlimmedPrimaryVertices_*_*',
                                                                         #'drop *_*_Cleaned_*',
                                                                         #'drop *_puppi_*_*',
                                                                         'keep *_muPFIso*_*_*' ),
                                  fileName       = cms.untracked.string ("Output.root")                                                                                                                   
)
# schedule definition                                                                                                       
process.outpath  = cms.EndPath(process.output) 

import FWCore.ParameterSet.Config as cms

process = cms.Process("Iso")
# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.GlobalTag.globaltag = 'START53_V7G::All'

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Iso')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring( *(
        ##'file:/cms/home/brochero/MuonIsolation/CMSSW_6_2_0_SLHC23_patch2/src/PhysicsTools/PatAlgos/DYToMuMu-Ph2_003E69A0-89A7-E411-82FC-0025905B8606.root', 
        ##'file:/cms/home/brochero/MuonIsolation/CMSSW_6_2_0_SLHC23_patch2/src/PhysicsTools/PatAlgos/QCD-Ph2_00169ED7-6BA4-E411-BA59-0025905AA9CC.root', 
        ##'file:/cms/data/xrd/store/mc/GEM2019Upg14DR/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-RECO/Phase1age1kJan23_PU140BX25_PH1_1K_FB_V3-v1/', 
        #'root://xrootd.unl.edu//store/mc/TP2023SHCALDR/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-RECO/SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1/00000/003E69A0-89A7-E411-82FC-0025905B8606.root',
        'file:/cms/home/brochero/IsolationPUPPITest/TaeTest/CMSSW_6_2_0_SLHC23_patch1/src/Dummy/Puppi/001B19B3-EBFA-E411-9E69-0025905B8606.root',
    ))
)

process.load('Dummy/Puppi/Puppi_cff')

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



process.iso = cms.EDAnalyzer('MuonIsoAnalyzer',
                             vxtTag         = cms.InputTag("offlinePrimaryVertices"),
                             useIPxy        = cms.untracked.bool(False),
                             useIPz         = cms.untracked.bool(False),
                             useDYGenVtx    = cms.untracked.bool(False),
                             MaxSimRecodz   = cms.untracked.double(0.5), 
                             MaxSimRecodrho = cms.untracked.double(0.2), 
                             IsoDepMuon = cms.VInputTag(cms.InputTag('muPFIsoDepositCHR04PUPPI'),
                                                        cms.InputTag('muPFIsoDepositPhR04PUPPI'),
                                                        cms.InputTag('muPFIsoDepositNHR04PUPPI')),
                             IsoValMuon = cms.VInputTag(cms.InputTag('muPFIsoValueCHR04PUPPI'),
                                                        cms.InputTag('muPFIsoValuePhR04PUPPI'),
                                                        cms.InputTag('muPFIsoValueNHR04PUPPI')),
                             
                             )
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('TreeMuonIso.root')
                                   )


#process.Tracer = cms.Service("Tracer")
#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.demo*process.dump)
process.p = cms.Path(process.pfPUPPISequence*
                     process.MuonPFIsoSequencePUPPI*
                     process.iso
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
process.outpath  = cms.EndPath(process.output)

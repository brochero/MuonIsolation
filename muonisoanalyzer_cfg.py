import FWCore.ParameterSet.Config as cms
 
process = cms.Process("Iso")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Iso')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring( *(
        ##'file:/cms/home/brochero/MuonIsolation/CMSSW_6_2_0_SLHC23_patch2/src/PhysicsTools/PatAlgos/DYToMuMu-Ph2_003E69A0-89A7-E411-82FC-0025905B8606.root', 
        'file:/cms/home/brochero/MuonIsolation/CMSSW_6_2_0_SLHC23_patch2/src/PhysicsTools/PatAlgos/QCD-Ph2_00169ED7-6BA4-E411-BA59-0025905AA9CC.root', 
        #'root://xrootd.unl.edu//store/mc/TP2023SHCALDR/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-RECO/SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1/00000/003E69A0-89A7-E411-82FC-0025905B8606.root',
    ))
)

process.iso = cms.EDAnalyzer('MuonIsoAnalyzer',
   vxtTag      = cms.InputTag("offlinePrimaryVertices"),
   useIPxy     = cms.untracked.bool(True),
   useIPz      = cms.untracked.bool(True),
   useDYGenVtx = cms.untracked.bool(False)
)
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('TreeMuonIso.root')
                                   )


#process.Tracer = cms.Service("Tracer")
#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.demo*process.dump)
process.p = cms.Path(process.iso)

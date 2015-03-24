// -*- C++ -*-
//
// Package:    PatAlgos
// Class:      IsoAnalyzer
// 
/**\class IsoAnalyzer IsoAnalyzer.cc PhysicsTools/PatAlgos/plugins/IsoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Javier Brochero Cifuentes
//         Created:  Fri, 06 Feb 2015 13:06:17 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <algorithm> // max
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//////////////////////////////////////////////////////////

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "Validation/RecoMuon/plugins/ME0MuonTrackCollProducer.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////



//
// class declaration
//

class MuonIsoAnalyzer : public edm::EDAnalyzer {
public:
  explicit MuonIsoAnalyzer(const edm::ParameterSet&);
  ~MuonIsoAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  bool IsLooseMuon(const edm::Event& iEvent, reco::MuonCollection::const_iterator& i_muon_candidate);  
  bool IsTightMuon(const edm::Event& iEvent, reco::MuonCollection::const_iterator& i_muon_candidate);  
  std::vector<double> findSimVtx(const edm::Event& iEvent);
  
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  
  edm::InputTag vxtTag;
  bool useIPxy, useIPz;
  bool useDYGenVtx;

  TH1F *h_muon_tight_number;
  TH1F *h_muon_loose_number;

  // Tree
  TTree *MuonTree;

  // Vertex
  std::vector<float> *b_vertex_x;
  std::vector<float> *b_vertex_y;
  std::vector<float> *b_vertex_z;
  std::vector<bool>  *b_vertex_isGood;

  // Muon
  std::vector<float> *b_muon_px;
  std::vector<float> *b_muon_py;
  std::vector<float> *b_muon_pz;
  std::vector<float> *b_muon_E;

  std::vector<bool> *b_muon_tight;

  // Muon Isolation
  std::vector<float> *b_muon_ChHadrIso;
  std::vector<float> *b_muon_NeHadrIso;
  std::vector<float> *b_muon_PhotonIso;
  std::vector<float> *b_muon_PUpTIso;
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonIsoAnalyzer::MuonIsoAnalyzer(const edm::ParameterSet& iConfig):
  vxtTag     (iConfig.getParameter< edm::InputTag>("vxtTag")),
  useIPxy    (iConfig.getUntrackedParameter<bool >("useIPxy",  true)),  
  useIPz     (iConfig.getUntrackedParameter<bool >("useIPz",   true)), 
  useDYGenVtx(iConfig.getUntrackedParameter<bool >("useDYGenVtx", false))  
{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;

  // Some basic histos  
  h_muon_tight_number  = fs->make<TH1F>("h_muon_tight_number" , "Number of tight muons" , 100 , 0 , 100 );
  h_muon_loose_number  = fs->make<TH1F>("h_muon_loose_number" , "Number of loose muons" , 100 , 0 , 100 );
  
  // Tree
  MuonTree= fs->make<TTree>("MuonTree","Muon Isolation Studies");
  
  // Primary Vertex
  MuonTree->Branch("vertex_x",     "std::vector<float>", &b_vertex_x);
  MuonTree->Branch("vertex_y",     "std::vector<float>", &b_vertex_y);
  MuonTree->Branch("vertex_z",     "std::vector<float>", &b_vertex_z);
  MuonTree->Branch("vertex_isGood","std::vector<bool>",  &b_vertex_isGood);
  
  // Muon
  MuonTree->Branch("muon_px"   ,"std::vector<float>",&b_muon_px);
  MuonTree->Branch("muon_py"   ,"std::vector<float>",&b_muon_py);
  MuonTree->Branch("muon_pz"   ,"std::vector<float>",&b_muon_pz);
  MuonTree->Branch("muon_E"    ,"std::vector<float>",&b_muon_E);
  MuonTree->Branch("muon_tight","std::vector<bool>" ,&b_muon_tight);

  // Muon Isolation
  MuonTree->Branch("ChHadrIso","std::vector<float>",&b_muon_ChHadrIso);
  MuonTree->Branch("NeHadrIso","std::vector<float>",&b_muon_NeHadrIso);
  MuonTree->Branch("PhotonIso","std::vector<float>",&b_muon_PhotonIso);
  MuonTree->Branch("PUpTIso",  "std::vector<float>",&b_muon_PUpTIso);
    
}


MuonIsoAnalyzer::~MuonIsoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
MuonIsoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Vertex
   b_vertex_x      = new std::vector<float>; 
   b_vertex_y      = new std::vector<float>;
   b_vertex_z      = new std::vector<float>;   
   b_vertex_isGood = new std::vector<bool>;

   // Muon
   b_muon_px    = new std::vector<float>;
   b_muon_py    = new std::vector<float>;
   b_muon_pz    = new std::vector<float>;
   b_muon_E     = new std::vector<float>;
   b_muon_tight = new std::vector<bool>;

   // Muon Isolation
   b_muon_ChHadrIso = new std::vector<float>;
   b_muon_NeHadrIso = new std::vector<float>;
   b_muon_PhotonIso = new std::vector<float>;
   b_muon_PUpTIso   = new std::vector<float>;

   // Vertex  
   edm::Handle<reco::VertexCollection> vertex;
   iEvent.getByLabel(vxtTag, vertex);
   const reco::VertexCollection& vtxs = *(vertex.product());
   
   int firstGoodVertex = -999;
   
   // Loop over vertices
   if (vtxs.size() != 0) {
     for (size_t i=0; i<vtxs.size(); i++) {
       
       bool isGoodVertex  = false;
       
       if ( fabs(vtxs[i].z())        < 24 &&
	    vtxs[i].position().Rho() < 2  &&
	    vtxs[i].ndof()           > 4  &&
	    !(vtxs[i].isFake())              ) {
	 
	 isGoodVertex = true;
	 
	 if (firstGoodVertex < 0) firstGoodVertex = i;
       }
       
       b_vertex_x          -> push_back(vtxs[i].x());
       b_vertex_y          -> push_back(vtxs[i].y());
       b_vertex_z          -> push_back(vtxs[i].z());
       b_vertex_isGood     -> push_back(isGoodVertex);
     }// for(vertex)
   } // if(vertex)
   
   //std::cout << "First good Vertex: " << firstGoodVertex << std::endl;

   edm::Handle<reco::MuonCollection> n_muons;
   iEvent.getByLabel("muons", n_muons);
   
   //edm::LogInfo("Iso") << "number of Muons "<<n_muons->size();
   
   int n_muon_loose=0;
   int n_muon_tight=0;
   
   for(reco::MuonCollection::const_iterator i_muon = n_muons->begin();
       i_muon != n_muons->end(); 
       ++ i_muon){

     if (IsLooseMuon(iEvent, i_muon)) {
       
       // Number of loose muons
       n_muon_loose++;

       // Is a tight Muon?
       if (IsTightMuon(iEvent, i_muon)) {
	 // Number of tight muons
	 n_muon_tight++;
	 b_muon_tight->push_back(true);
       }
       else b_muon_tight->push_back(false);
       
       // Kinematic variables
       b_muon_px->push_back(i_muon->px());
       b_muon_py->push_back(i_muon->py());
       b_muon_pz->push_back(i_muon->pz());
       b_muon_E->push_back(i_muon->energy());
       
       // Isolation variables
       b_muon_PUpTIso->push_back(i_muon->pfIsolationR04().sumPUPt);
       b_muon_ChHadrIso->push_back(i_muon->pfIsolationR04().sumChargedHadronPt);
       b_muon_NeHadrIso->push_back(i_muon->pfIsolationR04().sumNeutralHadronEt);
       b_muon_PhotonIso->push_back(i_muon->pfIsolationR04().sumPhotonEt);
       
     }// if(loose_muon)
   } // for(i_muon)
   
   // Histos
   h_muon_loose_number->Fill(n_muon_loose);
   h_muon_tight_number->Fill(n_muon_tight);

   // Fill Tree
   MuonTree->Fill();
      
   // Vertex
   delete b_vertex_x;
   delete b_vertex_y;
   delete b_vertex_z;  
   delete b_vertex_isGood;

   // Muon
   delete b_muon_px;
   delete b_muon_py;
   delete b_muon_pz;
   delete b_muon_E;

   // Muon Isolation
   delete b_muon_PUpTIso;
   delete b_muon_ChHadrIso;
   delete b_muon_NeHadrIso;
   delete b_muon_PhotonIso;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}

//---------------------------------------------------------
//------------- Find Simulated Vertex ---------------------
//---------------------------------------------------------

std::vector<double> MuonIsoAnalyzer::findSimVtx(const edm::Event& iEvent){

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  
  std::vector<double> vtxCoord;
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);
  vtxCoord.push_back(0);

  if(genParticles.isValid()){

    for(reco::GenParticleCollection::const_iterator itg = genParticles->begin(); itg != genParticles->end(); ++itg ){

      int id = itg->pdgId();
      int status = itg->status();
      //int nDaughters = itg->numberOfDaughters();
      int nMothers = itg->numberOfMothers();
      //double phiGen = itg->phi();
      //double etaGen = itg->eta();
      //std::cout<<"id "<<id<<" "<<status<<" "<<nMothers<<" "<<phiGen<<" "<<etaGen<<std::endl;

      // Z = 23 ; gamma = 22
      // status = 3: identifies the "hard part" of the interaction, i.e. the partons that are used in the matrix element calculation, including immediate decays of resonances. "This includes the two incoming colliding particles and partons produced in hard interaction." 
      if((abs(id) == 23 || abs(id) == 22) && status == 3){

  	vtxCoord[0] = 1;

  	vtxCoord[4] = (double)(itg->vx()); 
  	vtxCoord[5] = (double)(itg->vy());
  	vtxCoord[6] = (double)(itg->vz());

      }

      // status = 1: particle not decayed or fragmented, represents the final state as given by the generator
      else if(abs(id) == 13 && status == 1 && nMothers == 0){

  	vtxCoord[0] = 2;

  	vtxCoord[1] = (double)(itg->vx()); 
  	vtxCoord[2] = (double)(itg->vy());
  	vtxCoord[3] = (double)(itg->vz());

      }

      // nu_e = 12... for WHAT?
      else if(abs(id) == 12 && status == 1 && nMothers == 0){

  	vtxCoord[0] = 3;

  	vtxCoord[1] = (double)(itg->vx()); 
  	vtxCoord[2] = (double)(itg->vy());
  	vtxCoord[3] = (double)(itg->vz());

      }

    }
    
  }

  return vtxCoord;

}

//---------------------------------------------------------
//------------ Loose Muon Selection -----------------------
//---------------------------------------------------------

bool MuonIsoAnalyzer::IsLooseMuon(const edm::Event& iEvent, reco::MuonCollection::const_iterator &i_muon_candidate)
{

  bool isPF  = i_muon_candidate->isPFMuon();
  bool isGLB = i_muon_candidate->isGlobalMuon();
  bool isTrk = i_muon_candidate->isTrackerMuon();

  return (isPF && (isGLB || isTrk) );

}


//---------------------------------------------------------
//------------ Tight Muon Selection -----------------------
//---------------------------------------------------------

bool MuonIsoAnalyzer::IsTightMuon(const edm::Event& iEvent, reco::MuonCollection::const_iterator &i_muon_candidate)
{
  bool result=false;
  
  if (i_muon_candidate->muonBestTrack().isNonnull() && 
      i_muon_candidate->innerTrack().isNonnull() && 
      i_muon_candidate->globalTrack().isNonnull()){

    std::vector<double> vtxCoord = findSimVtx(iEvent);
    GlobalPoint point(vtxCoord[1],vtxCoord[2],vtxCoord[3]);
    GlobalPoint pointDY(vtxCoord[4],vtxCoord[5],vtxCoord[6]);

    //double muonX = i_muon_candidate->vx();
    //double muonY = i_muon_candidate->vy();
    //double muonZ = i_muon_candidate->vz();
    double DYZ = pointDY.z();


    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByLabel(vxtTag,vertexHandle);
    const reco::VertexCollection* vertices = vertexHandle.product();
    
    double distInit = 24;
    int indexFinal = -1;
    
    for(int i = 0; i < (int)vertices->size(); i++){


      if(useDYGenVtx == true){
	//double vtxX = (*vertices)[i].x();
	//double vtxY = (*vertices)[i].y();
	double vtxZ = (*vertices)[i].z();
	
	double dist = fabs(DYZ - vtxZ);
	//std::cout<<"dist "<<dist<<std::endl;
	if(dist < distInit){
	
	  distInit = dist;
	  indexFinal = i; // Closed vertex to the Z/gamma object
	  
	}
      } // if(DYGenVertx)
      else if(useDYGenVtx == false){
	if ( fabs((*vertices)[i].z())        < 24 &&
	     (*vertices)[i].position().Rho() < 2  &&
	     (*vertices)[i].ndof()           > 4  &&
	     !((*vertices)[i].isFake())              ) {
	    
	    if (indexFinal < 0) indexFinal = i;
	  }
	  
	}// else
	
    }// for(vertices)

    bool ipxy = false;
    bool ipz  = false;

    if(useDYGenVtx == true){ 
      double ipxySim = 999;
      double ipzSim = 999;
      
      if(vtxCoord[0] == 2 || vtxCoord[0] == 3){
	ipxySim = fabs(i_muon_candidate->muonBestTrack()->dxy(math::XYZPoint(point.x(),point.y(),point.z())));
	ipzSim  = fabs(i_muon_candidate->muonBestTrack()-> dz(math::XYZPoint(point.x(),point.y(),point.z())));
      }
      
      else if(vtxCoord[0] == 1 ){
	ipxySim = fabs(i_muon_candidate->muonBestTrack()->dxy(math::XYZPoint(pointDY.x(),pointDY.y(),pointDY.z())));
	ipzSim  = fabs(i_muon_candidate->muonBestTrack()-> dz(math::XYZPoint(pointDY.x(),pointDY.y(),pointDY.z())));	
      }
      
      bool ipxySimBool = ipxySim < 0.2;
      bool ipzSimBool = ipzSim < 0.5;
      
      if(vertices->size() !=0 && useIPxy == true){
	if(vtxCoord[0] == 1)      ipxy = fabs(i_muon_candidate->muonBestTrack()->dxy((*vertices)[indexFinal].position())) < 0.2;
	else if(vtxCoord[0] == 3) ipxy = fabs(i_muon_candidate->muonBestTrack()->dxy((*vertices)[0].position())) < 0.2;
	else ipxy = ipxySimBool;
      }
      else if(vertices->size() == 0 && useIPxy == true) ipxy = false;
      else if(useIPxy == false) ipxy = true;
      
      if(vertices->size() !=0 && useIPz == true){
	if(vtxCoord[0] == 1)     ipz = fabs(i_muon_candidate->muonBestTrack()->dz((*vertices)[indexFinal].position())) < 0.5;
	else if(vtxCoord[0] > 3) ipz = fabs(i_muon_candidate->muonBestTrack()->dz((*vertices)[0].position())) < 0.5;
	else ipz = ipzSimBool;
      }
      else if(vertices->size() == 0 && useIPz == true) ipz = false;
      else if(useIPz == false) ipz = true;
   
    } //if(DYGenVertx)
    
    else if(useDYGenVtx == false){
      if(useIPz  == true && 
	 useIPxy == true){
	ipz  = fabs(i_muon_candidate->muonBestTrack()->dz ((*vertices)[indexFinal].position())) < 0.5;
	ipxy = fabs(i_muon_candidate->muonBestTrack()->dxy((*vertices)[indexFinal].position())) < 0.2;
      }
      else if(useIPz  == true && 
	      useIPxy == false){
	ipz  = fabs(i_muon_candidate->muonBestTrack()->dz ((*vertices)[indexFinal].position())) < 0.5;
        ipxy = true;
      }
      else if(useIPz  == false && 
	      useIPxy == true){
	ipz  = true;
	ipxy = fabs(i_muon_candidate->muonBestTrack()->dxy((*vertices)[indexFinal].position())) < 0.2;
      }
      else if (useIPz  == true &&
	       useIPxy == true){
	ipz  = true;
        ipxy = true;
      }
    }

    bool trkLayMeas = i_muon_candidate->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5; 
    bool isGlb = i_muon_candidate->isGlobalMuon(); 
    bool isPF = i_muon_candidate->isPFMuon(); 
    bool chi2 = i_muon_candidate->globalTrack()->normalizedChi2() < 10.; 
    bool validHits = i_muon_candidate->globalTrack()->hitPattern().numberOfValidMuonHits() > 0; 
    bool matchedSt = i_muon_candidate->numberOfMatchedStations() > 1;


    //bool validPxlHit = i_muon_candidate->innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
    bool validPxlHit = i_muon_candidate->innerTrack()->hitPattern().pixelLayersWithMeasurement(3,2) > 0;
    //bool validPxlHit = i_muon_candidate->innerTrack()->hitPattern().pixelLayersWithMeasurement(4,3) > 0;

    if(trkLayMeas && 
       isGlb && 
       isPF && 
       chi2 && 
       validHits && 
       matchedSt && 
       ipxy && 
       ipz && 
       validPxlHit) result = true;
    
  }

  return result;
}



// ------------ method called once each job just before starting event loop  ------------
void 
MuonIsoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonIsoAnalyzer::endJob() 
{
  MuonTree->Print();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MuonIsoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MuonIsoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MuonIsoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MuonIsoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonIsoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonIsoAnalyzer);

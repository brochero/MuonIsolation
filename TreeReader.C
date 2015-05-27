#ifndef __CINT__

#include<string>
#include<iostream>
#include<sstream>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<set>
#include<vector>
#include <algorithm>

#endif

// Root
#include "TDirectory.h"
#include "TROOT.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TTree.h"
#include "Math/Point3D.h"
//#include "TKey.h"
//#include "TPrint.h"
//#include <exception>
#include <sys/stat.h>


#ifndef __CINT__

TH1F *EfficiencyHisto(TH1F *histo1, TH1F *histo2, TString name, TString title);
float dz (TLorentzVector muon, ROOT::Math::XYZPoint muon_BestTrack, ROOT::Math::XYZPoint vertex);
float dxy(TLorentzVector muon, ROOT::Math::XYZPoint muon_BestTrack, ROOT::Math::XYZPoint vertex);

void display_usage()
{
  std::cout << "\033[1;37musage:\033[1;m skimfile cutindex [options]" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "    -i inputfile  Input file without .root" << std::endl;
  std::cout << "    -o name in the output file \"h_\"" << std::endl;
  std::cout << "    -d Input file directory. Default directory: InputTrees" << std::endl;
  std::cout << "    -h                 displays this help message and exits " << std::endl;
  std::cout << "" << std::endl;
}


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const TString currentDateTime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d at %X", &tstruct);

  return buf;
}

int main(int argc, const char* argv[]){

  gSystem->Load("libTree");

  const char * _output   = 0;
  const char * _input    = 0;
  const char * _dir      = "/home/brochero/MuonIsolation/InputTree/";
  const char * _tr       = 0;
  const char * _idiso    = 0;

  // Arguments used
  //std::set<int> usedargs;
  //Parsing input options
  if(argc == 1){
    display_usage();
    return -1;
  }

  else{
      //Argumet 1 must be a valid input fileName
      for (int i = 1; i < argc; i++){
	if( strcmp(argv[i],"-i") == 0 ){
	  _input = argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-d") == 0 ){
	  _dir = argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-o") == 0 ){
	  _output= argv[i+1];
	  i++;
	}
	if( strcmp(argv[i],"-h") == 0 ||
	    strcmp(argv[i],"--help") == 0 ){
	  display_usage();
	  return 0;
	}
      }
  }//else
  if( _input ==0 ){
    std::cerr << "\033[1;31mskimfile ERROR:\033[1;m The '-i' option is mandatory!"
	      << std::endl;
    display_usage();
    return -1;
  }
  
  // reassigning
  TString fname(_input);
  TString hname(_output);
  TString fdir(_dir);
  
  TChain theTree("iso/MuonTree"); 
  
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "-----------------------------------------" << std::endl;
  std::cout << "Signal: ";
  std::cout << fname + ".root" << std::endl;
  
  theTree.Add(fdir + fname + ".root");

  float weight, pTHat;
  float PUInter, PUTrueInter;
  std::vector<float> *vertex_x=0, *vertex_y=0,*vertex_z=0,*vertex_isGood=0;
  std::vector<float> *vertex_sim_x=0, *vertex_sim_y=0,*vertex_sim_z=0, *sim_vertex_isGood=0;
  std::vector<float> *muon_px=0, *muon_py=0, *muon_pz=0, *muon_E=0; 
  std::vector<float> *muon_BestTrack_vx=0, *muon_BestTrack_vy=0, *muon_BestTrack_vz=0; 
  std::vector<bool>  *muon_IsTight=0;
  std::vector<float> *ChHadrIso=0, *NeHadrIso=0, *PhotonIso=0, *PUpTIso=0;  
  std::vector<float> *ChHadrPUPPI=0, *NeHadrPUPPI=0, *PhotonPUPPI=0;  

  // Event Info
  theTree.SetBranchAddress("event_weight", &weight);
  theTree.SetBranchAddress("event_pTHat" , &pTHat);

  // PU information
  theTree.SetBranchAddress("pu_num_interactions",      &PUInter);
  theTree.SetBranchAddress("pu_true_num_interactions", &PUTrueInter);  

  // SIM Vertex
  theTree.SetBranchAddress("vertex_sim_x", &vertex_sim_x);
  theTree.SetBranchAddress("vertex_sim_y", &vertex_sim_y);
  theTree.SetBranchAddress("vertex_sim_z", &vertex_sim_z);
  theTree.SetBranchAddress("vertex_sim_isGood", &sim_vertex_isGood);

  // RECO Vertex
  theTree.SetBranchAddress("vertex_x",      &vertex_x);
  theTree.SetBranchAddress("vertex_y",      &vertex_y);
  theTree.SetBranchAddress("vertex_z",      &vertex_z);
  theTree.SetBranchAddress("vertex_isGood", &vertex_isGood);

  // Muon
  theTree.SetBranchAddress( "muon_px", &muon_px );
  theTree.SetBranchAddress( "muon_py", &muon_py );
  theTree.SetBranchAddress( "muon_pz", &muon_pz );
  theTree.SetBranchAddress( "muon_E",  &muon_E  );
  theTree.SetBranchAddress( "muon_tight",  &muon_IsTight  );

  // Muon Best Track
  theTree.SetBranchAddress( "muon_BestTrack_vx", &muon_BestTrack_vx );
  theTree.SetBranchAddress( "muon_BestTrack_vy", &muon_BestTrack_vy );
  theTree.SetBranchAddress( "muon_BestTrack_vz", &muon_BestTrack_vz );

  // Muon Isolation
  theTree.SetBranchAddress( "ChHadrIso", &ChHadrIso );
  theTree.SetBranchAddress( "NeHadrIso", &NeHadrIso );
  theTree.SetBranchAddress( "PhotonIso", &PhotonIso );
  theTree.SetBranchAddress( "PUpTIso"  , &PUpTIso   );

  // Muon Isolation
  theTree.SetBranchAddress( "ChHadrPUPPI", &ChHadrPUPPI );
  theTree.SetBranchAddress( "NeHadrPUPPI", &NeHadrPUPPI );
  theTree.SetBranchAddress( "PhotonPUPPI", &PhotonPUPPI );

  
  /*********************************
             Histograms
  **********************************/
  TH1::SetDefaultSumw2(kTRUE);

  TH1F *hevent_weight, *hevent_pTHat; 
  TH1F *hevent_muon_weight; 

  TH1F *hPU_Internumber;

  TH1F *hvertex_SIMnumber;
  TH1F *hvertex_number, *hvertex_goodnumber, *hvertex_goodnumber_index;
  TH2F *h2Dvertex_SIMRECO_x, *h2Dvertex_SIMRECO_y, *h2Dvertex_SIMRECO_z, *h2Dvertex_SIMRECO_rho;  
  TH1F *hvertex_SIMRECO_dz, *hvertex_SIMRECO_drho;
  TH2F *h2Dvertex_matchedSIMRECO;

  TH1F *hmuon_number;
  TH1F *hmuon_IsoCh, *hmuon_IsoNe, *hmuon_IsoPh, *hmuon_IsoPU;
  TH1F *hmuon_pT, *hmuon_eta,  *hmuon_vertex,  *hmuon_SIMvertex,  *hmuon_TruePU;
  TH1F *hmuon_RelIso, *hmuon_RelIsoNoPU, *hmuon_RelIsoCh, *hmuon_RelPUPPI, *hmuon_RelChPUPPI;

  TH1F *hmuon_PUPPICh, *hmuon_PUPPINe, *hmuon_PUPPIPh;

  TH1F *hevent_passmuonCh_weight; 
  TH1F *hpassmuonCh_number;
  TH1F *hpassmuonCh_pT, *hpassmuonCh_eta, *hpassmuonCh_vertex, *hpassmuonCh_SIMvertex, *hpassmuonCh_TruePU;

  TH1F *hevent_passmuon_weight; 
  TH1F *hpassmuon_number;
  TH1F *hpassmuon_pT, *hpassmuon_eta, *hpassmuon_vertex, *hpassmuon_SIMvertex, *hpassmuon_TruePU;

  TH1F *hevent_passmuonNoPU_weight; 
  TH1F *hpassmuonNoPU_number;
  TH1F *hpassmuonNoPU_pT, *hpassmuonNoPU_eta, *hpassmuonNoPU_vertex, *hpassmuonNoPU_SIMvertex, *hpassmuonNoPU_TruePU;

  TH1F *hevent_passmuonPUPPI_weight; 
  TH1F *hpassmuonPUPPI_number;
  TH1F *hpassmuonPUPPI_pT, *hpassmuonPUPPI_eta, *hpassmuonPUPPI_vertex, *hpassmuonPUPPI_SIMvertex, *hpassmuonPUPPI_TruePU;

  TH1F *hevent_passmuonChPUPPI_weight; 
  TH1F *hpassmuonChPUPPI_number;
  TH1F *hpassmuonChPUPPI_pT, *hpassmuonChPUPPI_eta, *hpassmuonChPUPPI_vertex, *hpassmuonChPUPPI_SIMvertex, *hpassmuonChPUPPI_TruePU;

  TH1F *hmuon_IsoChCut, *hmuon_IsoCut, *hmuon_IsoNoPUCut, *hmuon_PUPPICut, *hmuon_ChPUPPICut;
  TProfile *pfmuon_pT_RelIso, *pfmuon_eta_RelIso;

  hevent_weight = new TH1F("hevent_weight","Weights",200,0.0,2.0);
  hevent_pTHat  = new TH1F("hevent_pTHat" ,"pTHat",500,0.0,500.0);

  hPU_Internumber         = new TH1F("hPU_Internumber","Number of PU Interactions",100,0,200);

  hvertex_SIMnumber         = new TH1F("hvertex_SIMnumber","Number of SIM PV",50,0,100);

  hvertex_number            = new TH1F("hvertex_number","Number of PV",100,0,200);
  hvertex_goodnumber        = new TH1F("hvertex_goodnumber","Number of good PV",100,0,200);
  hvertex_goodnumber_index  = new TH1F("hvertex_goodnumber_index","Index of the good of the good PV",5,0,5);
  h2Dvertex_SIMRECO_x = new TH2F("h2Dvertex_SIMRECO_x","Vertex X position: SIM vs RECO",400,-0.02,0.02,400,-0.02,0.02);
  h2Dvertex_SIMRECO_y = new TH2F("h2Dvertex_SIMRECO_y","Vertex Y position: SIM vs RECO",400,-0.02,0.02,400,-0.02,0.02);
  h2Dvertex_SIMRECO_z = new TH2F("h2Dvertex_SIMRECO_z","Vertex Z position: SIM vs RECO",400,-20,20,400,-20,20);
  h2Dvertex_SIMRECO_rho = new TH2F("h2Dvertex_SIMRECO_rho","Vertex Rho position: SIM vs RECO",100,0.0,0.01,100,0.0,0.01);
  h2Dvertex_matchedSIMRECO = new TH2F("hvertex_matchedSIMRECO","Number of GEN and RECO vtx matched",20,0,20, 20,0,20);
  hvertex_SIMRECO_dz   = new TH1F("hvertex_SIMRECO_dz","dz(RECO-SIM)",100,-0.005,0.005);
  hvertex_SIMRECO_drho = new TH1F("hvertex_SIMRECO_drho","drho(RECO-SIM)",100,-0.005,0.005);

  hevent_muon_weight = new TH1F("hevent_muon_weight","Weights pass muons",200,0.0,2.0);
  hmuon_number = new TH1F("hmuon_number","Number of muons per event",20,0,20);
  hmuon_pT     = new TH1F("hmuon_pT","muon pT",8,20,100);
  hmuon_eta    = new TH1F("hmuon_eta","muon #eta",12,0,2.4);
  hmuon_vertex = new TH1F("hmuon_vertex","Vertices (1 entry per muon)",20,0,200);
  hmuon_SIMvertex = new TH1F("hmuon_SIMvertex","SIM Vertices (1 entry per muon)",10,0,100);
  hmuon_TruePU  = new TH1F("hmuon_TruePU","True PU (1 entry per muon)",20,0,200);

  hevent_passmuonCh_weight = new TH1F("hevent_passmuonCh_weight","Weights Pass Ch. Isolated muons",200,0.0,2.0);
  hpassmuonCh_number       = new TH1F("hpassmuonCh_number","Number of Ch. isolated muons per event",20,0,20);
  hpassmuonCh_pT           = new TH1F("hpassmuonCh_pT","muon pT with Ch. Iso cut",8,20,100); // bin = 50 
  hpassmuonCh_eta          = new TH1F("hpassmuonCh_eta","muon #eta with Ch. Iso cut",12,0,2.4); // bin = 50
  hpassmuonCh_vertex       = new TH1F("hpassmuonCh_vertex","Vertices (1 entry per muon Ch. isolated)",20,0,200); // bin = 50
  hpassmuonCh_SIMvertex    = new TH1F("hpassmuonCh_SIMvertex","SIM Vertices (1 entry per muon Ch. isolated)",10,0,100); // bin = 50
  hpassmuonCh_TruePU       = new TH1F("hpassmuonCh_TruePU","True PU (1 entry per muon Ch. isolated)",20,0,200); // bin = 50

  hevent_passmuon_weight = new TH1F("hevent_passmuon_weight","Weights Pass Isolated muons",200,0.0,2.0);
  hpassmuon_number       = new TH1F("hpassmuon_number","Number of isolated muons per event",20,0,20);
  hpassmuon_pT           = new TH1F("hpassmuon_pT","muon pT with Iso cut",8,20,100);
  hpassmuon_eta          = new TH1F("hpassmuon_eta","muon #eta with Iso cut",12,0,2.4);
  hpassmuon_vertex       = new TH1F("hpassmuon_vertex","Vertices (1 entry per muon isolated)",20,0,200);
  hpassmuon_SIMvertex    = new TH1F("hpassmuon_SIMvertex","SIM Vertices (1 entry per muon isolated)",10,0,100);
  hpassmuon_TruePU       = new TH1F("hpassmuon_TruePU","True PU (1 entry per muon isolated)",20,0,200);

  hevent_passmuonNoPU_weight = new TH1F("hevent_passmuonNoPU_weight","Weights Pass No PU Isolated muons",200,0.0,2.0);
  hpassmuonNoPU_number       = new TH1F("hpassmuonNoPU_number","Number of No PU isolated muons per event",20,0,20);
  hpassmuonNoPU_pT           = new TH1F("hpassmuonNoPU_pT","muon pT with No PU Iso cut",8,20,100);
  hpassmuonNoPU_eta          = new TH1F("hpassmuonNoPU_eta","muon #eta with No PU Iso cut",12,0,2.4);
  hpassmuonNoPU_vertex       = new TH1F("hpassmuonNoPU_vertex","Vertices (1 entry per muon No PU isolated)",20,0,200);
  hpassmuonNoPU_SIMvertex    = new TH1F("hpassmuonNoPU_SIMvertex","SIM Vertices (1 entry per muon No PU isolated)",10,0,100);
  hpassmuonNoPU_TruePU       = new TH1F("hpassmuonNoPU_TruePU","True PU (1 entry per muon No PU isolated)",20,0,200);

  hevent_passmuonPUPPI_weight = new TH1F("hevent_passmuonPUPPI_weight","Weights Pass PUPPI Isolated muons",200,0.0,2.0);
  hpassmuonPUPPI_number       = new TH1F("hpassmuonPUPPI_number","Number of PUPPI isolated muons per event",20,0,20);
  hpassmuonPUPPI_pT           = new TH1F("hpassmuonPUPPI_pT","muon pT with PUPPI Iso cut",8,20,100);
  hpassmuonPUPPI_eta          = new TH1F("hpassmuonPUPPI_eta","muon #eta with PUPPI Iso cut",12,0,2.4);
  hpassmuonPUPPI_vertex       = new TH1F("hpassmuonPUPPI_vertex","Vertices (1 entry per muon PUPPI isolated)",20,0,200);
  hpassmuonPUPPI_SIMvertex    = new TH1F("hpassmuonPUPPI_SIMvertex","SIM Vertices (1 entry per muon PUPPI isolated)",10,0,100);
  hpassmuonPUPPI_TruePU       = new TH1F("hpassmuonPUPPI_TruePU","True PU (1 entry per muon PUPPI isolated)",20,0,200);

  hevent_passmuonChPUPPI_weight = new TH1F("hevent_passmuonChPUPPI_weight","Weights Pass PUPPI Ch. Isolated muons",200,0.0,2.0);
  hpassmuonChPUPPI_number       = new TH1F("hpassmuonChPUPPI_number","Number of PUPPI Ch. isolated muons per event",20,0,20);
  hpassmuonChPUPPI_pT           = new TH1F("hpassmuonChPUPPI_pT","muon pT with PUPPI Ch. Iso cut",8,20,100);
  hpassmuonChPUPPI_eta          = new TH1F("hpassmuonChPUPPI_eta","muon #eta with PUPPI Ch. Iso cut",12,0,2.4);
  hpassmuonChPUPPI_vertex       = new TH1F("hpassmuonChPUPPI_vertex","Vertices (1 entry per muon PUPPI Ch. isolated)",20,0,200);
  hpassmuonChPUPPI_SIMvertex    = new TH1F("hpassmuonChPUPPI_SIMvertex","SIM Vertices (1 entry per muon PUPPI Ch. isolated)",10,0,100);
  hpassmuonChPUPPI_TruePU       = new TH1F("hpassmuonChPUPPI_TruePU","True PU (1 entry per muon PUPPI Ch. isolated)",20,0,200);
    

  hmuon_IsoCh = new TH1F("hmuon_IsoCh","muon Iso. Charged",20,0,10);
  hmuon_IsoNe = new TH1F("hmuon_IsoNe","muon Iso. Neutral",20,0,10);
  hmuon_IsoPh = new TH1F("hmuon_IsoPh","muon Iso. Photon" ,30,0,20);
  hmuon_IsoPU = new TH1F("hmuon_IsoPU","muon Iso. PU"     ,40,0,50);

  hmuon_PUPPICh = new TH1F("hmuon_PUPPICh","muon PUPPI Charged Iso.",20,0,10);
  hmuon_PUPPINe = new TH1F("hmuon_PUPPINe","muon PUPPI Neutral Iso.",20,0,10);
  hmuon_PUPPIPh = new TH1F("hmuon_PUPPIPh","muon PUPPI Photon Iso." ,30,0,20);

  hmuon_RelIsoNoPU = new TH1F("hmuon_RelIsoNoPU","muon Rel. Isolation w/o PU mitigation",20,0,10);
  hmuon_RelIso     = new TH1F("hmuon_RelIso","muon Relative Isolation",20,0,2);
  hmuon_RelIsoCh   = new TH1F("hmuon_RelIsoCh","muon Rel. Charged Iso.",20,0,2);
  hmuon_RelPUPPI   = new TH1F("hmuon_RelPUPPI","muon Rel. PUPPI Iso.",20,0,2);
  hmuon_RelChPUPPI = new TH1F("hmuon_RelChPUPPI","muon Rel. PUPPI charged Iso.",20,0,2);

  hmuon_IsoCut     = new TH1F("hmuon_IsoCut", "Muon Rel. Isolation Cut", 80, 0.0, 0.8);
  hmuon_IsoChCut   = new TH1F("hmuon_IsoChCut", "Muon Rel. Charged Hadron Isolation Cut", 50, 0.0, 0.5);
  hmuon_IsoNoPUCut = new TH1F("hmuon_IsoNoPUCut", "Muon Rel. Isolation no PU substraction Cut", 50, 0.0, 1.5);
  hmuon_PUPPICut   = new TH1F("hmuon_PUPPICut", "Muon Rel. PUPPI Isolation Cut", 50, 0.0, 1.5);
  hmuon_ChPUPPICut = new TH1F("hmuon_ChPUPPICut", "Muon Rel. PUPPI Charged Isolation Cut", 50, 0.0, 1.5);

  pfmuon_pT_RelIso  = new TProfile("pfmuon_pT_RelIso","muon pT Vs Relative Isolation",8,20,100,0,2);
  pfmuon_eta_RelIso = new TProfile("pfmuon_eta_RelIso","muon #eta Vs Relative Isolation",12,0,2.4,0,2);

  /***********************************************************************************
                             Relative Isolation Cuts
   ***********************************************************************************/
  float RelIsoCh_cut, RelIso_cut, RelIsoNoPU_cut, RelPUPPI_cut, RelChPUPPI_cut;

  // Phase I No Age
  if(fname.Contains("Phase1NoAged")) {
    std::cout << "Phase I No Age sample........................" << std::endl;
    // Bkg rejection: 60%
    // RelIsoCh_cut   = 0.263; 
    // RelIso_cut     = 0.450;
    // RelIsoNoPU_cut = 0.630;
    // RelPUPPI_cut   = 0.561;
    // RelChPUPPI_cut = 0.280;

    // Bkg rejection: 80%
    RelIsoCh_cut   = 0.113; 
    RelIso_cut     = 0.234;
    RelIsoNoPU_cut = 0.400;
    RelPUPPI_cut   = 0.340;
    RelChPUPPI_cut = 0.126;
  }

  // Phase I + Age
  else if(fname.Contains("Phase1age")) {
    std::cout << "Phase I + Age sample........................" << std::endl;

    // Bkg rejection: 60%
    // RelIsoCh_cut   = 0.155; 
    // RelIso_cut     = 0.300;
    // RelIsoNoPU_cut = 0.670;
    // RelPUPPI_cut   = 0.400;
    // RelChPUPPI_cut = 0.210;

    // Bkg rejection: 80%
    RelIsoCh_cut   = 0.046; 
    RelIso_cut     = 0.105;
    RelIsoNoPU_cut = 0.400;
    RelPUPPI_cut   = 0.191;
    RelChPUPPI_cut = 0.071;
  }

  // Phase II
  else if(fname.Contains("PH2_1K_FB_V6-v1")) {
    std::cout << "Phase II sample........................" << std::endl;
    // Bkg rejection: 60%
    // RelIsoCh_cut   = 0.221; 
    // RelIso_cut     = 0.430;
    // RelIsoNoPU_cut = 0.940;
    // RelPUPPI_cut   = 0.610; 
    // RelChPUPPI_cut = 0.265; 

    // Bkg rejection: 60%
    RelIsoCh_cut   = 0.078; 
    RelIso_cut     = 0.204;
    RelIsoNoPU_cut = 0.661;
    RelPUPPI_cut   = 0.364; 
    RelChPUPPI_cut = 0.115; 
  }

  else {
    std::cout << "!!!!!!!!!!!!!!!! Error !!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "Unrecognized Sample. Impossible to set a Isolation cut!" << std::endl;
    return(0);
  }
  /***********************************************************************************
   ***********************************************************************************/

  float n_loosemuons=0.0;
  float n_tightmuons=0.0;

  TStopwatch sw;
  sw.Start(kTRUE);
  
  std::cout << "--- Processing: " << theTree.GetEntries() << " events" << std::endl;
  for (Long64_t ievt=0; ievt<theTree.GetEntries();ievt++) {
     
    theTree.GetEntry( ievt ); 
 
    //////////////////////////////////////////////////////
    int step = theTree.GetEntries()/50;
    if (ievt%(step) == 0){ 
      float progress=(ievt)/(theTree.GetEntries()*1.0);
      int barWidth = 50;
      
      std::cout << "[";
      int pos = barWidth * progress;
    
      for (int i = 0; i < barWidth; ++i) {
      	if (i < pos) std::cout << "=";
      	else if (i == pos) std::cout << ">";
      	else std::cout << " ";
      }
      
      std::cout << "] " << int(progress * 100.0) << " %\r";
      std::cout.flush();
    }
    ////////////////////////////////////////////////////////

    hevent_weight->Fill(weight);
    hevent_pTHat ->Fill(pTHat);


    // Vertex variables
    float n_vertex = 0;
    float n_goodvertex = 0;
    float n_SIMvertex = 0;
    bool matched_vtx = false;
    int matched_vtx_reco_index = -999;
    int matched_vtx_sim_index  = -999;
    double vertex_dz_simreco;
    double vertex_drho_simreco;
    
    /********************************************
        Mathing between RECO and SIM vertex 
    ********************************************/
    ROOT::Math::XYZPoint vertex_reco;
    ROOT::Math::XYZPoint vertex_sim;

    bool firstGoodVtx=false;
    
    // SIM
    for(int i_vtx_sim=0; i_vtx_sim < vertex_sim_x->size(); i_vtx_sim++){
      
      vertex_sim.SetXYZ((*vertex_sim_x)[i_vtx_sim], (*vertex_sim_y)[i_vtx_sim], (*vertex_sim_z)[i_vtx_sim]);
      
      if((*sim_vertex_isGood)[i_vtx_sim]){ // Geometrical RECO cuts (z<24: Rho<2)
	
	n_SIMvertex++;
      
      }	//if(good SIM vtx)    
    }// for(SIM vtx)


    // RECO
    for(int i_vtx=0; i_vtx < vertex_x->size(); i_vtx++){
      n_vertex++;
      
      vertex_reco.SetXYZ((*vertex_x)[i_vtx], (*vertex_y)[i_vtx], (*vertex_z)[i_vtx]);
      
      if((*vertex_isGood)[i_vtx]){
	n_goodvertex++;

	if (firstGoodVtx == false){
	  hvertex_goodnumber_index->Fill(i_vtx);
	  firstGoodVtx=true;
	}

	// Compare ONLY with the first SIM vtx
	vertex_sim.SetXYZ((*vertex_sim_x)[0], (*vertex_sim_y)[0], (*vertex_sim_z)[0]); 
	
	vertex_dz_simreco   = vertex_sim.Z()   - vertex_reco.Z();	
	vertex_drho_simreco = vertex_sim.Rho() - vertex_reco.Rho();	

	if(i_vtx == 0 && // First RECO vtx
	   fabs(vertex_dz_simreco)   < 0.5 &&                                                                                                                                                                
           fabs(vertex_drho_simreco) < 0.2){ 
  
	  matched_vtx =true;
	  matched_vtx_sim_index  = 0;
	  matched_vtx_reco_index = i_vtx;

	  h2Dvertex_SIMRECO_x->Fill(vertex_sim.X() , vertex_reco.X());
	  h2Dvertex_SIMRECO_y->Fill(vertex_sim.Y() , vertex_reco.Y());
	  h2Dvertex_SIMRECO_z->Fill(vertex_sim.Z() , vertex_reco.Z());
	  
	  h2Dvertex_SIMRECO_rho->Fill(vertex_sim.Rho() , vertex_reco.Rho());
	  
	  hvertex_SIMRECO_dz->Fill(vertex_dz_simreco,weight);
	  hvertex_SIMRECO_drho->Fill(vertex_drho_simreco,weight);
	  
	} // if(RECO-SIM match)
	
	// Compare with all SIM vtx

	// if(matched_vtx == false){
	//   for(int i_vtx_sim=0; i_vtx_sim < vertex_sim_x->size(); i_vtx_sim++){
	    
	//     vertex_sim.SetXYZ((*vertex_sim_x)[i_vtx_sim], (*vertex_sim_y)[i_vtx_sim], (*vertex_sim_z)[i_vtx_sim]);
	    
	//     vertex_dz_simreco   = vertex_sim.Z()   - vertex_reco.Z();	
	//     vertex_drho_simreco = vertex_sim.Rho() - vertex_reco.Rho();	
	    
	//     if(fabs(vertex_dz_simreco)   < 0.5 &&
	//        fabs(vertex_drho_simreco) < 0.2){
	      
	//       matched_vtx =true;
	//       matched_vtx_sim_index  = i_vtx_sim;
	//       matched_vtx_reco_index = i_vtx;
	      
	//       h2Dvertex_SIMRECO_x->Fill(vertex_sim.X() , vertex_reco.X());
	//       h2Dvertex_SIMRECO_y->Fill(vertex_sim.Y() , vertex_reco.Y());
	//       h2Dvertex_SIMRECO_z->Fill(vertex_sim.Z() , vertex_reco.Z());
	      
	//       h2Dvertex_SIMRECO_rho->Fill(vertex_sim.Rho() , vertex_reco.Rho());
	      
	//       hvertex_SIMRECO_dz->Fill(vertex_dz_simreco,weight);
	//       hvertex_SIMRECO_drho->Fill(vertex_drho_simreco,weight);
	      
	//       break;
	//     }
	//   }//for(vtx_sim)
	// } //if(match SIMRECO)
      }// if(GooVertex)
    }// for(i_vtx)


    // RECO vertex matched to SIM vertex	
    if(matched_vtx){
      
      // PU Plot
      hPU_Internumber->Fill(PUInter);
      
      // Primary Vertex      
      hvertex_SIMnumber->Fill(n_SIMvertex, weight);
      hvertex_number->Fill(n_vertex, weight);    
      hvertex_goodnumber->Fill(n_goodvertex, weight);    
      
      vertex_reco.SetXYZ((*vertex_x)[matched_vtx_reco_index], (*vertex_y)[matched_vtx_reco_index], (*vertex_z)[matched_vtx_reco_index]);
      h2Dvertex_matchedSIMRECO->Fill(matched_vtx_reco_index,matched_vtx_sim_index);
      
      float n_muons = 0;

      float n_passmuonsCh      = 0;
      float n_passmuons        = 0;
      float n_passmuonsNoPU    = 0;
      float n_passmuonsPUPPI   = 0;
      float n_passmuonsChPUPPI = 0;
      
      /********************************************
          Muon Selection (including vtx cuts) 
      ********************************************/
      for(int i_muon=0; i_muon < muon_px->size(); i_muon++){

	TLorentzVector muon;
	muon.SetPxPyPzE( (*muon_px)[i_muon], (*muon_py)[i_muon], (*muon_pz)[i_muon], (*muon_E)[i_muon] );
	ROOT::Math::XYZPoint muon_BestTrack;
	muon_BestTrack.SetXYZ((*muon_BestTrack_vx)[i_muon], (*muon_BestTrack_vy)[i_muon], (*muon_BestTrack_vz)[i_muon]);

	float muon_dz  = dz (muon, muon_BestTrack, vertex_reco);
	float muon_dxy = dxy(muon, muon_BestTrack, vertex_reco);
	
	if(muon.Pt() > 20 &&
	   fabs(muon.Eta()) < 2.4) n_loosemuons+= weight;

	if((*muon_IsTight)[i_muon] && 
	   muon.Pt() > 20 &&
	   fabs(muon.Eta()) < 2.4 &&
	   muon_dz  < 0.5 &&
	   muon_dxy < 0.2
	   ){
	  
	  n_tightmuons+= weight;
	  n_muons++;
	  
	  hevent_muon_weight->Fill(weight);

	  hmuon_pT       ->Fill(muon.Pt(),       weight);
	  hmuon_eta      ->Fill(fabs(muon.Eta()),weight);
	  hmuon_vertex   ->Fill(n_goodvertex,    weight);
	  hmuon_SIMvertex->Fill(n_SIMvertex,     weight);
	  hmuon_TruePU   ->Fill(PUInter,     weight);
	  
	  hmuon_IsoCh->Fill( (*ChHadrIso)[i_muon],weight );
	  hmuon_IsoNe->Fill( (*NeHadrIso)[i_muon],weight );
	  hmuon_IsoPh->Fill( (*PhotonIso)[i_muon],weight );
	  hmuon_IsoPU->Fill( (*PUpTIso)[i_muon],  weight );

	  // PUPPI
	  hmuon_PUPPICh->Fill( (*ChHadrPUPPI)[i_muon],weight );
	  hmuon_PUPPINe->Fill( (*NeHadrPUPPI)[i_muon],weight );
	  hmuon_PUPPIPh->Fill( (*PhotonPUPPI)[i_muon],weight );

	  //----------------------------------------------------------------------------------	  
	  //----------------------------------------------------------------------------------	  
	  // Relative Charged Hadron Isolation (cone 0.4)
	  //----------------------------------------------------------------------------------	  
	  //----------------------------------------------------------------------------------	  
	  float RelIsoCh=0.0;
	  RelIsoCh=(*ChHadrIso)[i_muon]/muon.Pt();

	  hmuon_RelIsoCh->Fill(RelIsoCh,weight);
	  
	  // Cut for the Relative Charged Hadron Isolation
	  hmuon_IsoChCut->Fill(10,weight); // Overflows to count the total number of muons
	  double isoChcut_step = 0.5;
	  for(int step=0; step < 50; step++){
	    if(RelIsoCh < isoChcut_step) hmuon_IsoChCut->Fill(isoChcut_step - (0.005),weight);
	    else continue;
	    
	    isoChcut_step = isoChcut_step - 0.01; // binning
	  }

	  // Relative Charged Hadron Isolation Cut
	  if (RelIsoCh<RelIsoCh_cut){
	    n_passmuonsCh++;
	    hevent_passmuonCh_weight->Fill(weight);
	    
	    hpassmuonCh_pT       ->Fill(muon.Pt(),       weight);
	    hpassmuonCh_eta      ->Fill(fabs(muon.Eta()),weight);
	    hpassmuonCh_vertex   ->Fill(n_goodvertex,    weight);
	    hpassmuonCh_SIMvertex->Fill(n_SIMvertex,     weight);
	    hpassmuonCh_TruePU   ->Fill(PUInter,         weight);
	  }
	  	
	  //----------------------------------------------------------------------------------	  
	  //----------------------------------------------------------------------------------	  
	  // Relative Isolation (cone 0.4)
	  //----------------------------------------------------------------------------------	  
	  //----------------------------------------------------------------------------------	  
	  float RelIso=0.0;
	  RelIso = ( (*ChHadrIso)[i_muon] + 
		     std::max(0.0, (*NeHadrIso)[i_muon] + (*PhotonIso)[i_muon] - 0.5*(*PUpTIso)[i_muon]) );
	  RelIso = RelIso/muon.Pt();
	  
	  hmuon_RelIso     ->Fill(RelIso,weight);
	  pfmuon_pT_RelIso ->Fill(muon.Pt(),RelIso);
	  pfmuon_eta_RelIso->Fill(fabs(muon.Eta()),RelIso);

	  // Cut for the Relative Isolation
	  hmuon_IsoCut->Fill(10,weight); // Overflows to count the total number of muons
	  double isocut_step = 0.8;
	  for(int step=0; step < 80; step++){
	    if(RelIso < isocut_step) hmuon_IsoCut->Fill(isocut_step - (0.005),weight);
	    else continue;
	    
	    isocut_step = isocut_step - 0.01; // binning
	  }
	  
	  // Relative Isolation Cut
	  if (RelIso<RelIso_cut){
	    n_passmuons++;
	    hevent_passmuon_weight->Fill(weight);

	    hpassmuon_pT       ->Fill(muon.Pt(),       weight);
	    hpassmuon_eta      ->Fill(fabs(muon.Eta()),weight);
	    hpassmuon_vertex   ->Fill(n_goodvertex,    weight);
	    hpassmuon_SIMvertex->Fill(n_SIMvertex,     weight);
	    hpassmuon_TruePU   ->Fill(PUInter,         weight);
	  }

	  //----------------------------------------------------------------------------------	  
	  //----------------------------------------------------------------------------------	  
	  // Relative Isolation w/o PU substraction (cone 0.4)
	  //----------------------------------------------------------------------------------	  
	  //----------------------------------------------------------------------------------	  
	  float RelIsoNoPU=0.0;
	  RelIsoNoPU = (*ChHadrIso)[i_muon] + (*NeHadrIso)[i_muon] + (*PhotonIso)[i_muon] ;
	  RelIsoNoPU = RelIsoNoPU/muon.Pt();
	  
	  hmuon_RelIsoNoPU->Fill(RelIsoNoPU,weight);
		  
	  // Cut for the Relative Isolation w/o PU substraction
	  hmuon_IsoNoPUCut->Fill(10,weight); // Overflows to count the total number of muons
	  double isoNoPUcut_step = 1.5;
	  for(int step=0; step < 50; step++){
	    if(RelIsoNoPU < isoNoPUcut_step) hmuon_IsoNoPUCut->Fill(isoNoPUcut_step - (0.015),weight);
	    else continue;
	    isoNoPUcut_step = isoNoPUcut_step - 0.03; // binning
	  }
	  
	  // Relative Isolation No PU Cut
	  if (RelIsoNoPU<RelIsoNoPU_cut){
	    n_passmuonsNoPU++;
	    hevent_passmuonNoPU_weight->Fill(weight);

	    hpassmuonNoPU_pT       ->Fill(muon.Pt(),       weight);
	    hpassmuonNoPU_eta      ->Fill(fabs(muon.Eta()),weight);
	    hpassmuonNoPU_vertex   ->Fill(n_goodvertex,    weight);
	    hpassmuonNoPU_SIMvertex->Fill(n_SIMvertex,     weight);
	    hpassmuonNoPU_TruePU   ->Fill(PUInter,         weight);
	  }


	  //----------------------------------------------------------------------------------	  
	  //----------------------------------------------------------------------------------	  
	  // Relative PUPPI (cone 0.4)
	  //----------------------------------------------------------------------------------	  
	  //----------------------------------------------------------------------------------	  
	  float RelPUPPI=0.0;
	  RelPUPPI = (*ChHadrPUPPI)[i_muon] + (*NeHadrPUPPI)[i_muon] + (*PhotonPUPPI)[i_muon] ;
	  RelPUPPI = RelPUPPI/muon.Pt();
	  
	  hmuon_RelPUPPI->Fill(RelPUPPI,weight);
		  
	  // Cut for the Relative Isolation w/o PU substraction
	  hmuon_PUPPICut->Fill(10,weight); // Overflows to count the total number of muons
	  double PUPPIcut_step = 1.5;
	  for(int step=0; step < 50; step++){
	    if(RelPUPPI < PUPPIcut_step) hmuon_PUPPICut->Fill(PUPPIcut_step - (0.015),weight);
	    else continue;
	    PUPPIcut_step = PUPPIcut_step - 0.03; // binning
	  }
	  
	  // Relative PUPPI Isolation cut
	  if (RelPUPPI<RelPUPPI_cut){
	    n_passmuonsPUPPI++;
	    hevent_passmuonPUPPI_weight->Fill(weight);

	    hpassmuonPUPPI_pT       ->Fill(muon.Pt(),       weight);
	    hpassmuonPUPPI_eta      ->Fill(fabs(muon.Eta()),weight);
	    hpassmuonPUPPI_vertex   ->Fill(n_goodvertex,    weight);
	    hpassmuonPUPPI_SIMvertex->Fill(n_SIMvertex,     weight);
	    hpassmuonPUPPI_TruePU   ->Fill(PUInter,         weight);
	  }


	  //----------------------------------------------------------------------------------	  
	  //----------------------------------------------------------------------------------	  
	  // Relative Charged PUPPI (cone 0.4)
	  //----------------------------------------------------------------------------------	  
	  //----------------------------------------------------------------------------------	  
	  float RelChPUPPI=0.0;
	  RelChPUPPI = (*ChHadrPUPPI)[i_muon]/muon.Pt();
	  
	  hmuon_RelChPUPPI->Fill(RelChPUPPI,weight);
		  
	  // Cut for the Relative Isolation w/o PU substraction
	  hmuon_ChPUPPICut->Fill(10,weight); // Overflows to count the total number of muons
	  double ChPUPPIcut_step = 1.5;
	  for(int step=0; step < 50; step++){
	    if(RelChPUPPI < ChPUPPIcut_step) hmuon_ChPUPPICut->Fill(ChPUPPIcut_step - (0.015),weight);
	    else continue;
	    ChPUPPIcut_step = ChPUPPIcut_step - 0.03; // binning
	  }
	  
	  // Relative PUPPI Charged Isolation  Cut
	  if (RelChPUPPI<RelChPUPPI_cut){
	    n_passmuonsChPUPPI++;
	    hevent_passmuonChPUPPI_weight->Fill(weight);

	    hpassmuonChPUPPI_pT       ->Fill(muon.Pt(),       weight);
	    hpassmuonChPUPPI_eta      ->Fill(fabs(muon.Eta()),weight);
	    hpassmuonChPUPPI_vertex   ->Fill(n_goodvertex,    weight);
	    hpassmuonChPUPPI_SIMvertex->Fill(n_SIMvertex,     weight);
	    hpassmuonChPUPPI_TruePU   ->Fill(PUInter,         weight);
	  }

	}//if(tight)
      }//for(muons)
      
      // Number of muons per event
      hmuon_number    ->Fill(n_muons,    weight);

      // Number of PASS muons per event
      hpassmuonCh_number->Fill(n_passmuonsCh,weight);
      hpassmuon_number->Fill(n_passmuons,weight);
      hpassmuonNoPU_number->Fill(n_passmuonsNoPU,weight);
      hpassmuonPUPPI_number->Fill(n_passmuonsPUPPI,weight);
      hpassmuonChPUPPI_number->Fill(n_passmuonsChPUPPI,weight);
      
    } //if(matched_vtx)
    
  }//for(events)
  
  // Vertex
  delete vertex_sim_x;
  delete vertex_sim_y;
  delete vertex_sim_z;
  
  // Vertex
  delete vertex_x;
  delete vertex_y;
  delete vertex_z;
  delete vertex_isGood;
  
  // Muon
  delete muon_px;
  delete muon_py;
  delete muon_pz;
  delete muon_E;

  // Muon Best Track
  delete muon_BestTrack_vx;
  delete muon_BestTrack_vy;
  delete muon_BestTrack_vz;
  
  // Muon Isolation                                                                                     
  delete ChHadrIso;
  delete NeHadrIso;
  delete PhotonIso;
  delete PUpTIso;

  // Muon PUPPI Isolation
  delete ChHadrPUPPI;
  delete NeHadrPUPPI;
  delete PhotonPUPPI;
  
  /***********************************************
          Some Additional Histograms  
  **********************************************/
  // Normalize to 1 efficiency histogram
  if(hmuon_number->Integral() != 0.0){
    hmuon_IsoChCut  ->Scale(1.0/(hmuon_IsoChCut  ->GetBinContent(hmuon_IsoChCut->GetNbinsX() + 1)));
    hmuon_IsoCut    ->Scale(1.0/(hmuon_IsoCut    ->GetBinContent(hmuon_IsoCut->GetNbinsX() + 1)));
    hmuon_IsoNoPUCut->Scale(1.0/(hmuon_IsoNoPUCut->GetBinContent(hmuon_IsoNoPUCut->GetNbinsX() + 1)));
    hmuon_PUPPICut  ->Scale(1.0/(hmuon_PUPPICut  ->GetBinContent(hmuon_PUPPICut->GetNbinsX() + 1)));
    hmuon_ChPUPPICut->Scale(1.0/(hmuon_ChPUPPICut  ->GetBinContent(hmuon_ChPUPPICut->GetNbinsX() + 1)));
  }

  // Relative Charged Hadron Isolation (cut-> 60% bkg reje.) efficiency as a function of kinematic variables
  TH1F *hmuon_Cheff_pT;
  hmuon_Cheff_pT = EfficiencyHisto(hpassmuonCh_pT, hmuon_pT, "hmuon_Cheff_pT", "Ch. Isolation Efficiency p_{T}");
  TH1F *hmuon_Cheff_eta;  
  hmuon_Cheff_eta = EfficiencyHisto(hpassmuonCh_eta, hmuon_eta, "hmuon_Cheff_eta", "Ch. Isolation Efficiency #eta");
  TH1F *hmuon_Cheff_vertex;
  hmuon_Cheff_vertex = EfficiencyHisto(hpassmuonCh_vertex, hmuon_vertex, "hmuon_Cheff_vertex", "Ch. Isolation Efficiency vertex");
  TH1F *hmuon_Cheff_SIMvertex;
  hmuon_Cheff_SIMvertex = EfficiencyHisto(hpassmuonCh_SIMvertex, hmuon_SIMvertex, "hmuon_Cheff_SIMvertex", "Ch. Isolation Efficiency SIM vertex");
  TH1F *hmuon_Cheff_TruePU;
  hmuon_Cheff_TruePU = EfficiencyHisto(hpassmuonCh_TruePU, hmuon_TruePU, "hmuon_Cheff_TruePU", "Ch. Isolation Efficiency True PU");

  // Relative Isolation (cut-> 60% bkg reje.) efficiency as a function of kinematic variables
  TH1F *hmuon_eff_pT;
  hmuon_eff_pT = EfficiencyHisto(hpassmuon_pT, hmuon_pT, "hmuon_eff_pT", "Isolation Efficiency p_{T}");
  TH1F *hmuon_eff_eta;  
  hmuon_eff_eta = EfficiencyHisto(hpassmuon_eta, hmuon_eta, "hmuon_eff_eta", "Isolation Efficiency #eta");
  TH1F *hmuon_eff_vertex;
  hmuon_eff_vertex = EfficiencyHisto(hpassmuon_vertex, hmuon_vertex, "hmuon_eff_vertex", "Isolation Efficiency vertex");
  TH1F *hmuon_eff_SIMvertex;
  hmuon_eff_SIMvertex = EfficiencyHisto(hpassmuon_SIMvertex, hmuon_SIMvertex, "hmuon_eff_SIMvertex", "Isolation Efficiency SIM vertex");
  TH1F *hmuon_eff_TruePU;
  hmuon_eff_TruePU = EfficiencyHisto(hpassmuon_TruePU, hmuon_TruePU, "hmuon_eff_TruePU", "Isolation Efficiency True PU");

  // Relative Isolation no PU reje. (cut 60% bkg reje.) efficiency as a function of kinematic variables
  TH1F *hmuon_NoPUeff_pT;
  hmuon_NoPUeff_pT = EfficiencyHisto(hpassmuonNoPU_pT, hmuon_pT, "hmuon_NoPUeff_pT", "Isolation (No PU miti) Efficiency p_{T}");
  TH1F *hmuon_NoPUeff_eta;  
  hmuon_NoPUeff_eta = EfficiencyHisto(hpassmuonNoPU_eta, hmuon_eta, "hmuon_NoPUeff_eta", "Isolation (No PU miti) Efficiency #eta");
  TH1F *hmuon_NoPUeff_vertex;
  hmuon_NoPUeff_vertex = EfficiencyHisto(hpassmuonNoPU_vertex, hmuon_vertex, "hmuon_NoPUeff_vertex", "Isolation (No PU miti) Efficiency vertex");
  TH1F *hmuon_NoPUeff_SIMvertex;
  hmuon_NoPUeff_SIMvertex = EfficiencyHisto(hpassmuonNoPU_SIMvertex, hmuon_SIMvertex, "hmuon_NoPUeff_SIMvertex", "Isolation (No PU miti) Efficiency SIM vertex");
  TH1F *hmuon_NoPUeff_TruePU;
  hmuon_NoPUeff_TruePU = EfficiencyHisto(hpassmuonNoPU_TruePU, hmuon_TruePU, "hmuon_NoPUeff_TruePU", "Isolation (No PU miti) Efficiency True PU");

  // Relative PUPPI Isolation (cut 60% bkg reje.) efficiency as a function of kinematic variables
  TH1F *hmuon_PUPPIeff_pT;
  hmuon_PUPPIeff_pT = EfficiencyHisto(hpassmuonPUPPI_pT, hmuon_pT, "hmuon_PUPPIeff_pT", "PUPPI Isolation Efficiency p_{T}");
  TH1F *hmuon_PUPPIeff_eta;  
  hmuon_PUPPIeff_eta = EfficiencyHisto(hpassmuonPUPPI_eta, hmuon_eta, "hmuon_PUPPIeff_eta", "PUPPI Isolation Efficiency #eta");
  TH1F *hmuon_PUPPIeff_vertex;
  hmuon_PUPPIeff_vertex = EfficiencyHisto(hpassmuonPUPPI_vertex, hmuon_vertex, "hmuon_PUPPIeff_vertex", "PUPPI Isolation Efficiency vertex");
  TH1F *hmuon_PUPPIeff_SIMvertex;
  hmuon_PUPPIeff_SIMvertex = EfficiencyHisto(hpassmuonPUPPI_SIMvertex, hmuon_SIMvertex, "hmuon_PUPPIeff_SIMvertex", "PUPPI Isolation Efficiency SIM vertex");
  TH1F *hmuon_PUPPIeff_TruePU;
  hmuon_PUPPIeff_TruePU = EfficiencyHisto(hpassmuonPUPPI_TruePU, hmuon_TruePU, "hmuon_PUPPIeff_TruePU", "PUPPI Isolation Efficiency True PU");

  // Relative PUPPI Charged Isolation (cut 60% bkg reje.) efficiency as a function of kinematic variables
  TH1F *hmuon_ChPUPPIeff_pT;
  hmuon_ChPUPPIeff_pT = EfficiencyHisto(hpassmuonChPUPPI_pT, hmuon_pT, "hmuon_ChPUPPIeff_pT", "PUPPI Charged Isolation Efficiency p_{T}");
  TH1F *hmuon_ChPUPPIeff_eta;  
  hmuon_ChPUPPIeff_eta = EfficiencyHisto(hpassmuonChPUPPI_eta, hmuon_eta, "hmuon_ChPUPPIeff_eta", "PUPPI Charged Isolation Efficiency #eta");
  TH1F *hmuon_ChPUPPIeff_vertex;
  hmuon_ChPUPPIeff_vertex = EfficiencyHisto(hpassmuonChPUPPI_vertex, hmuon_vertex, "hmuon_ChPUPPIeff_vertex", "PUPPI Charged Isolation Efficiency vertex");
  TH1F *hmuon_ChPUPPIeff_SIMvertex;
  hmuon_ChPUPPIeff_SIMvertex = EfficiencyHisto(hpassmuonChPUPPI_SIMvertex, hmuon_SIMvertex, "hmuon_ChPUPPIeff_SIMvertex", "PUPPI Charged Isolation Efficiency SIM vertex");
  TH1F *hmuon_ChPUPPIeff_TruePU;
  hmuon_ChPUPPIeff_TruePU = EfficiencyHisto(hpassmuonChPUPPI_TruePU, hmuon_TruePU, "hmuon_ChPUPPIeff_TruePU", "PUPPI charged Isolation Efficiency True PU");


  /***********************************************
              Output and finish  
  **********************************************/

  // Get elapsed time
  sw.Stop();
  std::cout << "==================================================] 100% " << std::endl;
  std::cout << "--- End of event loop: "; sw.Print();
  
  // Number of Muons
  std::cout << std::endl;
  std::cout << "--- Number of Loose Muons: " << n_loosemuons << std::endl;
  std::cout << "--- Number of Tight Muons: " << n_tightmuons << std::endl;
  std::cout << "--- Muon Efficiency (Tight/Loose): " << n_tightmuons/n_loosemuons << std::endl;


  
  //Output Dir
  TString dirname="MuonResults";   
  // make a dir if it does not exist!!
  struct stat st;
  if(stat(dirname,&st) != 0) system("mkdir " + dirname);
  
  TString systname="";
  bool matchname=false;
  
  TString temp_fname=fname + ".root";
  
  
  // --- Write histograms
  
  TString outfname=dirname + "/hSF-" + hname + "_" + fname + ".root";
  TFile *target  = new TFile(outfname,"RECREATE" );  

  hevent_weight->Write();
  hevent_pTHat ->Write();
  
  hPU_Internumber->Write();

  hvertex_SIMnumber->Write();
  hvertex_number->Write();
  hvertex_goodnumber->Write();
  hvertex_goodnumber_index->Write();

  h2Dvertex_SIMRECO_x->Write();
  h2Dvertex_SIMRECO_y->Write();
  h2Dvertex_SIMRECO_z->Write();
  h2Dvertex_SIMRECO_rho->Write();
  h2Dvertex_matchedSIMRECO->Write();
  hvertex_SIMRECO_dz->Write();
  hvertex_SIMRECO_drho->Write();

  hmuon_number->Write();

  hevent_muon_weight->Write();
  hmuon_pT->Write();
  hmuon_eta->Write();
  hmuon_vertex->Write();
  hmuon_SIMvertex->Write();
  hmuon_TruePU->Write();

  hevent_passmuonCh_weight->Write();
  hpassmuonCh_number->Write();
  hpassmuonCh_pT->Write();
  hpassmuonCh_eta->Write();
  hpassmuonCh_vertex->Write();
  hpassmuonCh_SIMvertex->Write();  
  hpassmuonCh_TruePU->Write();

  hevent_passmuon_weight->Write();
  hpassmuon_number->Write();
  hpassmuon_pT->Write();
  hpassmuon_eta->Write();
  hpassmuon_vertex->Write();
  hpassmuon_SIMvertex->Write();
  hpassmuon_TruePU->Write();

  hevent_passmuonNoPU_weight->Write();
  hpassmuonNoPU_number->Write();
  hpassmuonNoPU_pT->Write();
  hpassmuonNoPU_eta->Write();
  hpassmuonNoPU_vertex->Write();
  hpassmuonNoPU_SIMvertex->Write();
  hpassmuonNoPU_TruePU->Write();

  hevent_passmuonPUPPI_weight->Write();
  hpassmuonPUPPI_number->Write();
  hpassmuonPUPPI_pT->Write();
  hpassmuonPUPPI_eta->Write();
  hpassmuonPUPPI_vertex->Write();
  hpassmuonPUPPI_SIMvertex->Write();
  hpassmuonPUPPI_TruePU->Write();

  hevent_passmuonChPUPPI_weight->Write();
  hpassmuonChPUPPI_number->Write();
  hpassmuonChPUPPI_pT->Write();
  hpassmuonChPUPPI_eta->Write();
  hpassmuonChPUPPI_vertex->Write();
  hpassmuonChPUPPI_SIMvertex->Write();
  hpassmuonChPUPPI_TruePU->Write();

  hmuon_IsoCh->Write();
  hmuon_IsoNe->Write();
  hmuon_IsoPh->Write();
  hmuon_IsoPU->Write();

  hmuon_PUPPICh->Write();
  hmuon_PUPPINe->Write();
  hmuon_PUPPIPh->Write();

  hmuon_RelIsoCh->Write();
  hmuon_RelIso->Write();
  hmuon_RelIsoNoPU->Write();
  hmuon_RelPUPPI->Write();
  hmuon_RelChPUPPI->Write();

  pfmuon_pT_RelIso->Write();
  pfmuon_eta_RelIso->Write();

  hmuon_IsoChCut->Write();
  hmuon_IsoCut->Write();
  hmuon_IsoNoPUCut->Write();
  hmuon_PUPPICut->Write();
  hmuon_ChPUPPICut->Write();

  hmuon_Cheff_pT->Write();
  hmuon_Cheff_eta->Write();
  hmuon_Cheff_vertex->Write();
  hmuon_Cheff_SIMvertex->Write();
  hmuon_Cheff_TruePU->Write();

  hmuon_eff_pT->Write();
  hmuon_eff_eta->Write();
  hmuon_eff_vertex->Write();
  hmuon_eff_SIMvertex->Write();
  hmuon_eff_TruePU->Write();

  hmuon_NoPUeff_pT->Write();
  hmuon_NoPUeff_eta->Write();
  hmuon_NoPUeff_vertex->Write();
  hmuon_NoPUeff_SIMvertex->Write();
  hmuon_NoPUeff_TruePU->Write();

  hmuon_PUPPIeff_pT->Write();
  hmuon_PUPPIeff_eta->Write();
  hmuon_PUPPIeff_vertex->Write();
  hmuon_PUPPIeff_SIMvertex->Write();
  hmuon_PUPPIeff_TruePU->Write();

  hmuon_ChPUPPIeff_pT->Write();
  hmuon_ChPUPPIeff_eta->Write();
  hmuon_ChPUPPIeff_vertex->Write();
  hmuon_ChPUPPIeff_SIMvertex->Write();
  hmuon_ChPUPPIeff_TruePU->Write();

  std::cout << "File saved as " << outfname << std::endl;

  target->Close();

}

  /***********************************************
             Efficiency Histograms
  **********************************************/
TH1F* EfficiencyHisto(TH1F *histo1, TH1F *histo2, TString name, TString title){
  
  std::vector<float>  v_eff, v_eff_error;
  float v_eff_bin[200];
  int i_bin_new = 0;
  
  float npass_hist        = 0.0;
  float npass_error_hist2 = 0.0;
  float n_hist            = 0.0;
  float n_error_hist2     = 0.0;
  float eff               = 0.0;
  float eff_error2        = 0.0;
  
  // Fill the beginning of the HISTO binning 
  v_eff_bin[i_bin_new]=(histo1->GetBinLowEdge(1)); // low-edge
  i_bin_new++;

  for(int i_bin=1; i_bin < (histo1->GetNbinsX() + 1) ; i_bin++){
    
    npass_hist += histo1->GetBinContent(i_bin); 
    n_hist     += histo2->GetBinContent(i_bin); 
    
    npass_error_hist2 += (histo1->GetBinError(i_bin))*(histo1->GetBinError(i_bin)); 
    n_error_hist2     += (histo2->GetBinError(i_bin))*(histo2->GetBinError(i_bin)); 
    
    if(n_hist != 0.0){
      eff        = npass_hist / n_hist;
      eff_error2 = (npass_hist / n_hist)*(npass_hist / n_hist)*
	( npass_error_hist2*(1.0/(npass_hist*npass_hist)) + (n_error_hist2*(1.0/(n_hist*n_hist))) );
    }
    else {
      eff        = 0.0;
      eff_error2 = 0.0;
    }
    
    //if (sqrt(eff_error2) < 0.2){ // max error to define the binning
    if (sqrt(eff_error2) < 1000000000000000000000.2){ // max error to define the binning
      v_eff.push_back(eff);
      v_eff_error.push_back(sqrt(eff_error2));
      v_eff_bin[i_bin_new]=(histo1->GetBinLowEdge(i_bin) + histo1->GetBinWidth(i_bin)); //upper-edge
      
      i_bin_new ++; // new number of bins

      npass_hist        = 0.0;
      n_hist            = 0.0;
      npass_error_hist2 = 0.0;
      n_error_hist2     = 0.0;
    } // if(error)

  } // for(i_bin)
  
  // Fill the last bin if it is still not filled.
  if (eff != 0.0){

    v_eff_bin[i_bin_new] = (histo1->GetBinLowEdge(histo1->GetNbinsX())  + histo1->GetBinWidth(histo1->GetNbinsX())); 
    v_eff.push_back(eff);
    v_eff_error.push_back(sqrt(eff_error2));
    i_bin_new ++; // new number of bins
  }
  i_bin_new-=2; // remove the last 2 entry (if the bin size is equal to the original histogram)
  
  
  TH1F *hist_eff;
  hist_eff = new TH1F(name, title, i_bin_new,v_eff_bin);
  // Fill the Eff. histo with the new binning
  for(int i_bin=1; i_bin<(v_eff.size() + 1); i_bin++){

    hist_eff->SetBinContent(i_bin,v_eff[i_bin-1]);
    hist_eff->SetBinError(i_bin,v_eff_error[i_bin-1]);

  }

  return hist_eff;
}

float dz (TLorentzVector muon, ROOT::Math::XYZPoint muon_BestTrack, ROOT::Math::XYZPoint vertex){
  float dz =  (muon_BestTrack.Z() - vertex.Z()) - 
    ( ( (muon_BestTrack.X() - vertex.X())*muon.Px() + 
	(muon_BestTrack.Y() - vertex.Y())*muon.Py() 
      )/muon.Pt()
    ) * muon.Pz() / muon.Pt();

  return fabs(dz);
}

float dxy(TLorentzVector muon, ROOT::Math::XYZPoint muon_BestTrack, ROOT::Math::XYZPoint vertex){
  float dxy =  -1.0 * 
    ( (muon_BestTrack.X() - vertex.X())*muon.Py() + 
      (muon_BestTrack.Y() - vertex.Y())*muon.Px() 
    )/muon.Pt() ;

  return fabs(dxy);
}

#endif


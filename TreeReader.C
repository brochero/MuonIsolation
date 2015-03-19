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
//#include "TKey.h"
//#include "TPrint.h"
//#include <exception>
#include <sys/stat.h>


#ifndef __CINT__

TH1F* EfficiencyHisto(TH1F *histo1, TH1F *histo2, TString name, TString title);

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

  
  std::vector<float> *vertex_x=0, *vertex_y=0,*vertex_z=0,*vertex_isGood=0;
  std::vector<float> *muon_px=0, *muon_py=0, *muon_pz=0, *muon_E=0; 
  std::vector<bool>  *muon_IsTight=0;
  std::vector<float> *ChHadrIso=0, *NeHadrIso=0, *PhotonIso=0, *PUpTIso=0;  

  // Vertex
  theTree.SetBranchAddress("b_vertex_x",      &vertex_x);
  theTree.SetBranchAddress("b_vertex_y",      &vertex_y);
  theTree.SetBranchAddress("b_vertex_z",      &vertex_z);
  theTree.SetBranchAddress("b_vertex_isGood", &vertex_isGood);

  // Muon
  theTree.SetBranchAddress( "muon_px", &muon_px );
  theTree.SetBranchAddress( "muon_py", &muon_py );
  theTree.SetBranchAddress( "muon_pz", &muon_pz );
  theTree.SetBranchAddress( "muon_E",  &muon_E  );
  theTree.SetBranchAddress( "muon_tight",  &muon_IsTight  );

  // Muon Isolation
  theTree.SetBranchAddress( "ChHadrIso", &ChHadrIso );
  theTree.SetBranchAddress( "NeHadrIso", &NeHadrIso );
  theTree.SetBranchAddress( "PhotonIso", &PhotonIso );
  theTree.SetBranchAddress( "PUpTIso"  , &PUpTIso   );

  
  /*********************************
             Histograms
  **********************************/
  
  TH1F *hvertex_number, *hvertex_goodnumber;
  TH1F * hmuon_number, *hpassmuon_number;
  TH1F *hmuon_IsoCh, *hmuon_IsoNe, *hmuon_IsoPh, *hmuon_IsoPU;
  TH1F *hmuon_pT, *hmuon_eta,  *hmuon_vertex, *hmuon_RelIso, *hmuon_RelIsoNoPU;
  TH1F *hpassmuon_pT, *hpassmuon_eta, *hpassmuon_vertex;
  TH1F *hmuon_IsoCut;
  TProfile *hmuon_pT_RelIso, *hmuon_eta_RelIso;

  hvertex_number     = new TH1F("hvertex_number","Number of PV",100,0,2000);
  hvertex_goodnumber = new TH1F("hvertex_goodnumber","Number of good PV",100,0,2000);

  hmuon_number = new TH1F("hmuon_number","Number of muons per event",20,0,20);
  hmuon_pT     = new TH1F("hmuon_pT","muon pT",100,0,200);
  hmuon_eta    = new TH1F("hmuon_eta","muon #eta",20,0,2.5);
  hmuon_vertex = new TH1F("hmuon_vertex","Vertices (1 entry per muon)",50,0,2000);

  hpassmuon_number = new TH1F("hpassmuon_number","Number of isolated muons per event",20,0,20);
  hpassmuon_pT     = new TH1F("hpassmuon_pT","muon pT with Iso cut",100,0,200);
  hpassmuon_eta    = new TH1F("hpassmuon_eta","muon #eta with Iso cut",20,0,2.5);
  hpassmuon_vertex = new TH1F("hpassmuon_vertex","Vertices (1 entry per muon isolated)",50,0,2000);
  
  hmuon_IsoCh = new TH1F("hmuon_IsoCh","muon Iso. Charged",40,0,10);
  hmuon_IsoNe = new TH1F("hmuon_IsoNe","muon Iso. Neutral",40,0,10);
  hmuon_IsoPh = new TH1F("hmuon_IsoPh","muon Iso. Photon" ,60,0,20);
  hmuon_IsoPU = new TH1F("hmuon_IsoPU","muon Iso. PU"     ,80,0,50);

  hmuon_RelIsoNoPU = new TH1F("hmuon_RelIsoNoPU","muon Rel. Isolation w/o PU mitigation",40,0,10);
  hmuon_RelIso     = new TH1F("hmuon_RelIso","muon Relative Isolation",40,0,2);

  hmuon_IsoCut  = new TH1F("hmuon_IsoCut", "Muon Rel. Isolation Cut", 50, 0.0, 0.5);

  hmuon_pT_RelIso  = new TProfile("hmuon_pT_RelIso","muon pT Vs Relative Isolation",40,0,200,0,2);
  hmuon_eta_RelIso = new TProfile("hmuon_eta_RelIso","muon #eta Vs Relative Isolation",20,0,2.5,0,2);


  TStopwatch sw;
  sw.Start(kTRUE);
  
  std::cout << "--- Processing: " << theTree.GetEntries() << " events" << std::endl;
  for (Long64_t ievt=0; ievt<theTree.GetEntries();ievt++) {
     
    theTree.GetEntry( ievt ); 
 
    //////////////////////////////////////////////////////
    if (ievt%100000 == 0){ 
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

    int n_vertex=0;
    int n_goodvertex=0;
    for(int i_vertex=0; i_vertex < vertex_x->size(); i_vertex++){

      TVector3 vertex;
      vertex.SetXYZ((*vertex_x)[i_vertex],(*vertex_y)[i_vertex],(*vertex_z)[i_vertex]);
      n_vertex++;

      if((*vertex_isGood)[i_vertex]) n_goodvertex++;
      
    }// for(i_vertex)

    hvertex_number->Fill(n_vertex);    
    hvertex_goodnumber->Fill(n_goodvertex);    

    int n_muons=0;
    int n_passmuons=0;

    for(int i_muon=0; i_muon < muon_px->size(); i_muon++){

      TLorentzVector muon;
      muon.SetPxPyPzE( (*muon_px)[i_muon], (*muon_py)[i_muon], (*muon_pz)[i_muon], (*muon_E)[i_muon] );
	
      if(muon_IsTight && 
	 muon.Pt() > 5 &&
	 fabs(muon.Eta()) < 2.4){

	n_muons++;

	hmuon_IsoCh->Fill( (*ChHadrIso)[i_muon] );
	hmuon_IsoNe->Fill( (*NeHadrIso)[i_muon] );
	hmuon_IsoPh->Fill( (*PhotonIso)[i_muon] );
	hmuon_IsoPU->Fill( (*PUpTIso)[i_muon] );
	
	float RelIso=0.0;
	RelIso = ( (*ChHadrIso)[i_muon] + 
		   std::max(0.0, (*NeHadrIso)[i_muon] + (*PhotonIso)[i_muon] - 0.5*(*PUpTIso)[i_muon]) );
	RelIso = RelIso/muon.Pt();
	
	float RelIsoNoPU=0.0;
	RelIsoNoPU = (*ChHadrIso)[i_muon] + (*NeHadrIso)[i_muon] + (*PhotonIso)[i_muon] ;
	RelIsoNoPU = RelIsoNoPU/muon.Pt();
	
	hmuon_pT->Fill(muon.Pt());
	hmuon_eta->Fill(fabs(muon.Eta()));
	hmuon_vertex->Fill(n_goodvertex);

	hmuon_RelIso->Fill(RelIso);
	hmuon_RelIsoNoPU->Fill(RelIsoNoPU);
	
	hmuon_pT_RelIso->Fill(muon.Pt(),RelIso);
	hmuon_eta_RelIso->Fill(fabs(muon.Eta()),RelIso);
	
	// Isolation Cut
	hmuon_IsoCut->Fill(10); // Overflows to count the total number of muons
	double isocut_step = 0.5;
	for(int step=0; step < 50; step++){
	  if(RelIso < isocut_step) hmuon_IsoCut->Fill(isocut_step - (0.005));
	  else continue;

	  isocut_step = isocut_step - 0.01; // binning
	}
	  
	// Isolation Cut
	if (RelIso<0.2){
	  n_passmuons++;
	  hpassmuon_pT ->Fill(muon.Pt());
	  hpassmuon_eta->Fill(fabs(muon.Eta()));
	  hpassmuon_vertex->Fill(n_goodvertex);
	}
	
      }//if(tight)
    }//for(muons)
    
    hmuon_number    ->Fill(n_muons);
    hpassmuon_number->Fill(n_passmuons);
    
  }//for(events)
  
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

  // Muon Isolation                                                                                     
  delete ChHadrIso;
  delete NeHadrIso;
  delete PhotonIso;
  delete PUpTIso;
  
  /***********************************************
          Some Additional Histograms  
  **********************************************/
  hmuon_IsoCut->Sumw2();
  if(hmuon_number->Integral() != 0.0) hmuon_IsoCut->Scale(1.0/(hmuon_IsoCut->GetBinContent(hmuon_IsoCut->GetNbinsX() + 1)));

  TH1F *hmuon_eff_pT;
  hmuon_eff_pT = EfficiencyHisto(hpassmuon_pT, hmuon_pT, "hmuon_eff_pT", "Isolation Efficiency p_{T}");
  TH1F *hmuon_eff_eta;  
  hmuon_eff_eta = EfficiencyHisto(hpassmuon_eta, hmuon_eta, "hmuon_eff_eta", "Isolation Efficiency #eta");
  TH1F *hmuon_eff_vertex;
  hmuon_eff_vertex = EfficiencyHisto(hpassmuon_vertex, hmuon_vertex, "hmuon_eff_vertex", "Isolation Efficiency vertex");

  /***********************************************
              Output and finish  
  **********************************************/

  // Get elapsed time
  sw.Stop();
  std::cout << "==================================================] 100% " << std::endl;
  std::cout << "--- End of event loop: "; sw.Print();
  
  
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

  hvertex_number->Write();
  hvertex_goodnumber->Write();

  hmuon_number->Write();
  hpassmuon_number->Write();

  hmuon_pT->Write();
  hmuon_eta->Write();
  hmuon_vertex->Write();
  hpassmuon_pT->Write();
  hpassmuon_eta->Write();
  hpassmuon_vertex->Write();

  hmuon_IsoCh->Write();
  hmuon_IsoNe->Write();
  hmuon_IsoPh->Write();
  hmuon_IsoPU->Write();

  hmuon_RelIsoNoPU->Write();
  hmuon_RelIso->Write();

  hmuon_pT_RelIso->Write();
  hmuon_eta_RelIso->Write();

  hmuon_IsoCut->Write();

  hmuon_eff_pT->Write();
  hmuon_eff_eta->Write();
  hmuon_eff_vertex->Write();


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
  
  float npass_hist       = 0.0;
  float n_hist           = 0.0;
  float npass_error_hist = 0.0;
  float n_error_hist     = 0.0;
  float eff              = 0.0;
  float eff_error        = 0.0;
  
  for(int i_bin=1; i_bin < (histo1->GetNbinsX() + 1) ; i_bin++){
    
    npass_hist += histo1->GetBinContent(i_bin); 
    n_hist     += histo2->GetBinContent(i_bin); 
    
    npass_error_hist += (histo1->GetBinError(i_bin))*(histo1->GetBinError(i_bin)); 
    n_error_hist     += (histo2->GetBinError(i_bin))*(histo2->GetBinError(i_bin)); 
    
    if(n_hist != 0.0){
      eff       = npass_hist / n_hist;
      eff_error = (npass_hist / n_hist)*(npass_hist / n_hist)*
	( npass_error_hist*(1.0/(npass_hist*npass_hist)) + (n_error_hist*(1.0/(n_hist*n_hist))) );
    }
    else {
      eff       = 0.0;
      eff_error = 0.0;
    }
    
    if (sqrt(eff_error) < 0.02){ // max error to define the binning
      v_eff.push_back(eff);
      v_eff_error.push_back(sqrt(eff_error));
      v_eff_bin[i_bin_new]=(histo1->GetBinLowEdge(i_bin));

      i_bin_new ++; // new number of bins
      npass_hist       = 0.0;
      n_hist           = 0.0;
      npass_error_hist = 0.0;
      n_error_hist     = 0.0;
    } // if(error)

  } // for(i_bin)
  
  // Fill the last bin if it is still not filled.
  if (eff != 0.0){
    v_eff.push_back(eff);
    v_eff_error.push_back(sqrt(eff_error));
    v_eff_bin[i_bin_new] = (histo1->GetBinLowEdge(histo1->GetNbinsX() + 1));
  }

  TH1F *hist_eff;
  hist_eff = new TH1F(name, title, i_bin_new,v_eff_bin);
  // Fill the Eff. histo with the new binning
  for(int i_bin=1; i_bin<(v_eff.size() + 1); i_bin++){

    hist_eff->SetBinContent(i_bin,v_eff[i_bin-1]);
    hist_eff->SetBinError(i_bin,v_eff_error[i_bin-1]);

  }

  return hist_eff;
}


#endif


#include <sys/stat.h>
#include <unistd.h>
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMatrix.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TLegend.h"
#include "THStack.h"
#include "TImage.h"
#include "TMarker.h"
#include "TLatex.h"
#include <vector>
#include <iostream>
#include <ostream>
#include <sstream>
#include <fstream>

void add_plot_label( char* label, double x, double y, double size = 0.05, int color = 1, int font = 62, int align = 22 ){

  TLatex *latex = new TLatex( x, y, label );
  latex->SetNDC();
  latex->SetTextSize(size);
  latex->SetTextColor(color);
  latex->SetTextFont(font);
  latex->SetTextAlign(align);
  latex->Draw();
}


int multiple_detector_fit()
{

  std::cout << "Beginning : ... " << std::endl;

  Int_t npoints = 1000;
  Double_t emin = 0.2; 
  Double_t emax = 3.0;

  bool use100m = true;
  bool use470m = true;
  bool use600m = true;


  std::vector<int> baselines;
  std::vector<double> scales;
  std::vector<std::string> names;
  std::vector<double> volume;

  if (use100m) baselines.push_back(100);
  if (use470m) baselines.push_back(470);
  if (use600m) baselines.push_back(600);


  double NULLVec[2][20];
  double NULLVec_beam[2][20];
  double NULLVec_pimin[2][20];
  double NULLVec_piplus[2][20];
  double NULLVec_K0[2][20];
  double NULLVec_Kmin[2][20];
  double NULLVec_Kplus[2][20];

  double OscVec[2][1001][7][20];
  double OscVec_beam[2][1001][7][20];
  double OscVec_pimin[2][1001][7][20];
  double OscVec_piplus[2][1001][7][20];
  double OscVec_K0[2][1001][7][20];
  double OscVec_Kmin[2][1001][7][20];
  double OscVec_Kplus[2][1001][7][20];


  for(int j = 0; j < 2; j++){
    for(int i = 0; i < 20; i++){
    
      NULLVec[j][i] = 0;
      NULLVec_beam[j][i] = 0;
      NULLVec_pimin[j][i] = 0;
      NULLVec_piplus[j][i] = 0;
      NULLVec_K0[j][i] = 0;
      NULLVec_Kmin[j][i] = 0;
      NULLVec_Kplus[j][i] = 0;

    }
  }

  for(int j = 0; j < 2; j++){
    for(int u = 0; u < 1000; u++){
      for(int s = 0; s < 7; s++){
	for(int i = 0; i < 20; i++){

	OscVec[j][u][s][i] = 0;
        OscVec_beam[j][u][s][i] = 0;
        OscVec_pimin[j][u][s][i] = 0;
        OscVec_piplus[j][u][s][i] = 0;
        OscVec_K0[j][u][s][i] = 0;
        OscVec_Kmin[j][u][s][i] = 0;
        OscVec_Kplus[j][u][s][i] = 0;

	}
      }
    }
  }


  int nbinsE = 0;
  int counter = 0;

  double bins[20] = {.200, .300, .400, .450, .500, .550, .600, .650, .700, .750, .800, .850, .900, .950, 1.000, 1.250, 1.500, 2.000, 2.500, 3.000};

  TH1D *ND_err = new TH1D("ND","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *ND_err_beam = new TH1D("ND_beam","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *ND_err_pimin = new TH1D("ND_pimin","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *ND_err_piplus = new TH1D("ND_piplus","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *ND_err_K0 = new TH1D("ND_K0","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *ND_err_Kmin = new TH1D("ND_Kmin","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *ND_err_Kplus = new TH1D("ND_Kplus","; Neutrino Energy [GeV];Events",19, bins);;

  TH1D *FD_err = new TH1D("FD","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *FD_err_beam = new TH1D("FD_beam","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *FD_err_pimin = new TH1D("FD_pimin","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *FD_err_piplus = new TH1D("FD_piplus","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *FD_err_K0 = new TH1D("FD_K0","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *FD_err_Kmin = new TH1D("FD_Kmin","; Neutrino Energy [GeV];Events",19, bins);;
  TH1D *FD_err_Kplus = new TH1D("FD_Kplus","; Neutrino Energy [GeV];Events",19, bins);;


  if (use100m){

    std::string temp_name = "../MatrixFiles/combined_ntuple_100m_nu_processed_SystSplit_numu.root";

    TFile temp_file(temp_name.c_str());
    TH1D *NULL_100;TH1D *NULL_100_beam;TH1D *NULL_100_pimin;TH1D *NULL_100_piplus;TH1D *NULL_100_K0;TH1D *NULL_100_Kmin;TH1D *NULL_100_Kplus;
    NULL_100 =        (TH1D*)(temp_file.Get("NumuCC"));
    //    NULL_100_beam =   (TH1D*)(temp_file.Get("NumuCC_beam"));
    // NULL_100_pimin =  (TH1D*)(temp_file.Get("NumuCC_pimin"));
    // NULL_100_piplus = (TH1D*)(temp_file.Get("NumuCC_piplus"));
    // NULL_100_K0 =     (TH1D*)(temp_file.Get("NumuCC_K0"));
    // NULL_100_Kmin =   (TH1D*)(temp_file.Get("NumuCC_Kmin"));
    //NULL_100_Kplus =  (TH1D*)(temp_file.Get("NumuCC_Kplus"));

    nbinsE = NULL_100->GetNbinsX();
    std::cout << nbinsE << std::endl;
      for(int i = 1; i <= nbinsE; i++){
	  NULLVec[counter][i-1] =        (NULL_100->GetBinContent(i));
	  //     NULLVec_beam[counter][i-1] =   (NULL_100_beam->GetBinContent(i));
	  // NULLVec_pimin[counter][i-1] =  (NULL_100_pimin->GetBinContent(i));
	  // NULLVec_piplus[counter][i-1] = (NULL_100_piplus->GetBinContent(i));
	  // NULLVec_K0[counter][i-1] =     (NULL_100_K0->GetBinContent(i));
	  // NULLVec_Kmin[counter][i-1] =   (NULL_100_Kmin->GetBinContent(i));
          //NULLVec_Kplus[counter][i-1] =  (NULL_100_Kplus->GetBinContent(i));


      }

    for(int u = 0; u < npoints; u++){
      for(int s = 0; s < 7; s++){
	TH1D *OSC_100;
	TString upoint = Form("%d",u);
	TString name = "Universe_";
	TString name2 = "_MultiSim_";
	TString mul = Form("%d",s);


	TString nom; nom = name+upoint+name2+mul;

	TString syst_beam; TString beam = "_beam"; syst_beam = name+upoint+name2+mul+beam;
	TString syst_pimin; TString pimin = "_pimin"; syst_pimin = name+upoint+name2+mul+pimin;
	TString syst_piplus; TString piplus = "_piplus"; syst_piplus = name+upoint+name2+mul+piplus;
	TString syst_K0; TString K0 = "_K0"; syst_K0 = name+upoint+name2+mul+K0;
	TString syst_Kmin; TString Kmin = "_Kmin"; syst_Kmin = name+upoint+name2+mul+Kmin;
	TString syst_Kplus; TString Kplus = "_Kplus"; syst_Kplus = name+upoint+name2+mul+Kplus;
       
	OSC_100 = (TH1D*)(temp_file.Get(nom));
	//        OSC_100_beam = (TH1D*)(temp_file.Get(syst_beam));
	// OSC_100_pimin = (TH1D*)(temp_file.Get(syst_pimin));
	// OSC_100_piplus = (TH1D*)(temp_file.Get(syst_piplus));
	// OSC_100_K0 = (TH1D*)(temp_file.Get(syst_K0));
        //OSC_100_Kmin = (TH1D*)(temp_file.Get(syst_Kmin));
        //OSC_100_Kplus = (TH1D*)(temp_file.Get(syst_Kplus));


          for(int i = 1; i <= nbinsE; i++){
              OscVec[counter][u][s][i-1] =        (OSC_100->GetBinContent(i));
	      //    OscVec_beam[counter][u][s][i-1] =   (OSC_100_beam->GetBinContent(i));
	      // OscVec_pimin[counter][u][s][i-1] =  (OSC_100_pimin->GetBinContent(i));
              //OscVec_piplus[counter][u][s][i-1] = (OSC_100_piplus->GetBinContent(i));
	      // OscVec_K0[counter][u][s][i-1] =     (OSC_100_K0->GetBinContent(i));
	      // OscVec_Kmin[counter][u][s][i-1] =   (OSC_100_Kmin->GetBinContent(i));
              //OscVec_Kplus[counter][u][s][i-1] =  (OSC_100_Kplus->GetBinContent(i));

	  }

	delete OSC_100;
	//	delete OSC_100_beam;
	// delete OSC_100_pimin;
	// delete OSC_100_piplus;
        //delete OSC_100_K0;
	// delete OSC_100_Kmin;
        //delete OSC_100_Kplus;

      }
    }
    
    delete NULL_100;
    //    delete NULL_100_beam;
    // delete NULL_100_pimin;
    //delete NULL_100_piplus;
    //delete NULL_100_K0;
    //delete NULL_100_Kmin;
    // delete NULL_100_Kplus;

    temp_file.Close();
    counter++;
  }

  if (use470m){

    std::string temp_name = "../MatrixFiles/combined_ntuple_600m_onaxis_nu_processed_SystSplit_numu.root";
    TFile temp_file(temp_name.c_str());
    TH1D *NULL_470;TH1D *NULL_470_beam;TH1D *NULL_470_pimin;TH1D *NULL_470_piplus;TH1D *NULL_470_K0;TH1D *NULL_470_Kmin;TH1D *NULL_470_Kplus;
    NULL_470 =        (TH1D*)(temp_file.Get("NumuCC"));
        NULL_470_beam =   (TH1D*)(temp_file.Get("NumuCC_beam"));
    NULL_470_pimin =  (TH1D*)(temp_file.Get("NumuCC_pimin"));
     NULL_470_piplus = (TH1D*)(temp_file.Get("NumuCC_piplus"));
    NULL_470_K0 =     (TH1D*)(temp_file.Get("NumuCC_K0"));
     NULL_470_Kmin =   (TH1D*)(temp_file.Get("NumuCC_Kmin"));
    NULL_470_Kplus =  (TH1D*)(temp_file.Get("NumuCC_Kplus"));

    nbinsE = NULL_470->GetNbinsX();
    std::cout << nbinsE << std::endl;
      for(int i = 1; i <= nbinsE; i++){
	  NULLVec[counter][i-1] =        (NULL_470->GetBinContent(i));
	      NULLVec_beam[counter][i-1] =   (NULL_470_beam->GetBinContent(i));
	   NULLVec_pimin[counter][i-1] =  (NULL_470_pimin->GetBinContent(i));
          NULLVec_piplus[counter][i-1] = (NULL_470_piplus->GetBinContent(i));
	   NULLVec_K0[counter][i-1] =     (NULL_470_K0->GetBinContent(i));
          NULLVec_Kmin[counter][i-1] =   (NULL_470_Kmin->GetBinContent(i));
          NULLVec_Kplus[counter][i-1] =  (NULL_470_Kplus->GetBinContent(i));


      }

    for(int u = 0; u < npoints; u++){
      for(int s = 0; s < 7; s++){
	TH1D *OSC_470;
	TString upoint = Form("%d",u);
	TString name = "Universe_";
	TString name2 = "_MultiSim_";
	TString mul = Form("%d",s);


	TString nom; nom = name+upoint+name2+mul;

	TString syst_beam; TString beam = "_beam"; syst_beam = name+upoint+name2+mul+beam;
	TString syst_pimin; TString pimin = "_pimin"; syst_pimin = name+upoint+name2+mul+pimin;
	TString syst_piplus; TString piplus = "_piplus"; syst_piplus = name+upoint+name2+mul+piplus;
	TString syst_K0; TString K0 = "_K0"; syst_K0 = name+upoint+name2+mul+K0;
	TString syst_Kmin; TString Kmin = "_Kmin"; syst_Kmin = name+upoint+name2+mul+Kmin;
	TString syst_Kplus; TString Kplus = "_Kplus"; syst_Kplus = name+upoint+name2+mul+Kplus;
       
	OSC_470 = (TH1D*)(temp_file.Get(nom));
	        OSC_470_beam = (TH1D*)(temp_file.Get(syst_beam));
	 OSC_470_pimin = (TH1D*)(temp_file.Get(syst_pimin));
        OSC_470_piplus = (TH1D*)(temp_file.Get(syst_piplus));
	 OSC_470_K0 = (TH1D*)(temp_file.Get(syst_K0));
        OSC_470_Kmin = (TH1D*)(temp_file.Get(syst_Kmin));
        OSC_470_Kplus = (TH1D*)(temp_file.Get(syst_Kplus));


          for(int i = 1; i <= nbinsE; i++){
              OscVec[counter][u][s][i-1] =        (OSC_470->GetBinContent(i));
	          OscVec_beam[counter][u][s][i-1] =   (OSC_470_beam->GetBinContent(i));
              OscVec_pimin[counter][u][s][i-1] =  (OSC_470_pimin->GetBinContent(i));
              OscVec_piplus[counter][u][s][i-1] = (OSC_470_piplus->GetBinContent(i));
              OscVec_K0[counter][u][s][i-1] =     (OSC_470_K0->GetBinContent(i));
              OscVec_Kmin[counter][u][s][i-1] =   (OSC_470_Kmin->GetBinContent(i));
              OscVec_Kplus[counter][u][s][i-1] =  (OSC_470_Kplus->GetBinContent(i));

	  }

	delete OSC_470;
		delete OSC_470_beam;
	 delete OSC_470_pimin;
        delete OSC_470_piplus;
        delete OSC_470_K0;
        delete OSC_470_Kmin;
	 delete OSC_470_Kplus;

      }
    }
    
    delete NULL_470;
        delete NULL_470_beam;
    delete NULL_470_pimin;
    delete NULL_470_piplus;
    delete NULL_470_K0;
    delete NULL_470_Kmin;
    delete NULL_470_Kplus;

    temp_file.Close();
    counter++;
  }

  int nL = 2;
  int mbins = (nbinsE*nL);
 
  TMatrix M6 (mbins,mbins);
  TMatrix M5 (mbins,mbins);
  TMatrix M4 (mbins,mbins);
  TMatrix M3 (mbins,mbins);
  TMatrix M2 (mbins,mbins);
  TMatrix M1 (mbins,mbins);
  TMatrix M0 (mbins,mbins);
  
  TMatrix C6 (mbins,mbins);
  TMatrix C5 (mbins,mbins);
  TMatrix C4 (mbins,mbins);
  TMatrix C3 (mbins,mbins);
  TMatrix C2 (mbins,mbins);
  TMatrix C1 (mbins,mbins);
  TMatrix C0 (mbins,mbins);

  int N = 0;

    TH1D *Fig6 = new TH1D("Fig6",";;",mbins,0,mbins);
    TH1D *Fig5 = new TH1D("Fig5",";;",mbins,0,mbins);
    TH1D *Fig4 = new TH1D("Fig4",";;",mbins,0,mbins);
    TH1D *Fig3 = new TH1D("Fig3",";;",mbins,0,mbins);
    TH1D *Fig2 = new TH1D("Fig2",";;",mbins,0,mbins);
    TH1D *Fig1 = new TH1D("Fig1",";;",mbins,0,mbins);
    TH1D *Fig0 = new TH1D("Fig0",";;",mbins,0,mbins);

  int Erri = 0, Errj = 0;

  std::cout << "Filling Error Matrix..." << std::endl;

  for(int Lrow = 0; Lrow < 2; Lrow++){
    for(int Erow = 0; Erow < nbinsE; Erow++){

      Errj = 0;

      for(int Lcol = 0; Lcol < 2; Lcol++){
        for(int Ecol = 0; Ecol < nbinsE; Ecol++){

          M6 (Erri,Errj) = 0;
          M5 (Erri,Errj) = 0;
          M4 (Erri,Errj) = 0;
          M3 (Erri,Errj) = 0;
          M2 (Erri,Errj) = 0;
          M1 (Erri,Errj) = 0;
          M0 (Erri,Errj) = 0;

	  N = 0;

	  for(int u = 0; u < npoints; u++){


	    
	    M6 (Erri,Errj) += (NULLVec[Lrow][Erow]-OscVec[Lrow][u][6][Erow])*(NULLVec[Lcol][Ecol]-OscVec[Lcol][u][6][Ecol]);
            M5 (Erri,Errj) += (NULLVec_Kplus[Lrow][Erow]-OscVec_Kplus[Lrow][u][5][Erow])*(NULLVec_Kplus[Lcol][Ecol]-OscVec_Kplus[Lcol][u][5][Ecol]);
            M4 (Erri,Errj) += (NULLVec_Kmin[Lrow][Erow]-OscVec_Kmin[Lrow][u][4][Erow])*(NULLVec_Kmin[Lcol][Ecol]-OscVec_Kmin[Lcol][u][4][Ecol]);
            M3 (Erri,Errj) += (NULLVec_K0[Lrow][Erow]-OscVec_K0[Lrow][u][3][Erow])*(NULLVec_K0[Lcol][Ecol]-OscVec_K0[Lcol][u][3][Ecol]);
            M2 (Erri,Errj) += (NULLVec_piplus[Lrow][Erow]-OscVec_piplus[Lrow][u][2][Erow])*(NULLVec_piplus[Lcol][Ecol]-OscVec_piplus[Lcol][u][2][Ecol]);
            M1 (Erri,Errj) += (NULLVec_pimin[Lrow][Erow]-OscVec_pimin[Lrow][u][1][Erow])*(NULLVec_pimin[Lcol][Ecol]-OscVec_pimin[Lcol][u][1][Ecol]);
            M0 (Erri,Errj) += (NULLVec_beam[Lrow][Erow]-OscVec_beam[Lrow][u][0][Erow])*(NULLVec_beam[Lcol][Ecol]-OscVec_beam[Lcol][u][0][Ecol]);
	    
	    N++;
	  }
	  
	  M6 (Erri,Errj) /= N;

          M5 (Erri,Errj) /= N;

          M4 (Erri,Errj) /= N;
	  
          M3 (Erri,Errj) /= N;
	  
          M2 (Erri,Errj) /= N;

          M1 (Erri,Errj) /= N;
	  
          M0 (Erri,Errj) /= N;
	  
	  
	  M6 (Erri,Errj) /= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];
	  M5 (Erri,Errj) /= NULLVec_Kplus[Lrow][Erow]*NULLVec_Kplus[Lcol][Ecol];
	
	  if( NULLVec_Kmin[Lrow][Erow]*NULLVec_Kmin[Lcol][Ecol] != 0){ M4 (Erri,Errj) /= NULLVec_Kmin[Lrow][Erow]*NULLVec_Kmin[Lcol][Ecol];}
	  else{ M4 (Erri,Errj) = 0; }

	  M3 (Erri,Errj) /= NULLVec_K0[Lrow][Erow]*NULLVec_K0[Lcol][Ecol];
	  M2 (Erri,Errj) /= NULLVec_piplus[Lrow][Erow]*NULLVec_piplus[Lcol][Ecol];
	  M1 (Erri,Errj) /= NULLVec_pimin[Lrow][Erow]*NULLVec_pimin[Lcol][Ecol];
	  M0 (Erri,Errj) /= NULLVec_beam[Lrow][Erow]*NULLVec_beam[Lcol][Ecol];
	  
	  if(Erri == Errj) Fig6->SetBinContent(Erri+1, sqrt(M6 (Erri,Errj)));
	  if(Erri == Errj) Fig5->SetBinContent(Erri+1, sqrt(M5 (Erri,Errj)));
	  if(Erri == Errj) Fig4->SetBinContent(Erri+1, sqrt(M4 (Erri,Errj)));
	  if(Erri == Errj) Fig3->SetBinContent(Erri+1, sqrt(M3 (Erri,Errj)));
	  if(Erri == Errj) Fig2->SetBinContent(Erri+1, sqrt(M2 (Erri,Errj)));
	  if(Erri == Errj) Fig1->SetBinContent(Erri+1, sqrt(M1 (Erri,Errj)));
	  if(Erri == Errj) Fig0->SetBinContent(Erri+1, sqrt(M0 (Erri,Errj)));
	  if(Erri == Errj && Lcol == 0){ 
	    std::cout << " ND bin " << Ecol+1 << std::endl;
	    ND_err->SetBinContent(Ecol+1, 100*sqrt(M6 (Erri,Errj)));
	    ND_err_Kplus->SetBinContent(Ecol+1, 100*sqrt(M5 (Erri,Errj)));
            ND_err_Kmin->SetBinContent(Ecol+1, 100*sqrt(M4 (Erri,Errj)));
            ND_err_K0->SetBinContent(Ecol+1, 100*sqrt(M3 (Erri,Errj)));
            ND_err_piplus->SetBinContent(Ecol+1, 100*sqrt(M2 (Erri,Errj)));
            ND_err_pimin->SetBinContent(Ecol+1, 100*sqrt(M1 (Erri,Errj)));
            ND_err_beam->SetBinContent(Ecol+1, 100*sqrt(M0 (Erri,Errj)));
	  }

          if(Erri == Errj && Lcol == 1){ 
	    std::cout << " FD bin " << Ecol+1 << std::endl;
	    FD_err->SetBinContent(Ecol+1, 100*sqrt(M6 (Erri,Errj)));
            FD_err_Kplus->SetBinContent(Ecol+1, 100*sqrt(M5 (Erri,Errj)));
            FD_err_Kmin->SetBinContent(Ecol+1, 100*sqrt(M4 (Erri,Errj)));
            FD_err_K0->SetBinContent(Ecol+1, 100*sqrt(M3 (Erri,Errj)));
            FD_err_piplus->SetBinContent(Ecol+1, 100*sqrt(M2 (Erri,Errj)));
            FD_err_pimin->SetBinContent(Ecol+1, 100*sqrt(M1 (Erri,Errj)));
            FD_err_beam->SetBinContent(Ecol+1, 100*sqrt(M0 (Erri,Errj)));
	  }

	  //	  std::cout << M6 (Erri,Errj) << "\t";

          Errj++;

	}}

      Erri++;

    }}

  for(int i = 0; i < Erri; i++){
    for(int j = 0; j < Errj; j++){

      C6 (i,j) = M6(i,j) / sqrt(M6 (i,i) * M6 (j,j));
      C5 (i,j) = M5(i,j) / sqrt(M5 (i,i) * M5 (j,j));
      C4 (i,j) = M4(i,j) / sqrt(M4 (i,i) * M4 (j,j));
      C3 (i,j) = M3(i,j) / sqrt(M3 (i,i) * M3 (j,j));
      C2 (i,j) = M2(i,j) / sqrt(M2 (i,i) * M2 (j,j));
      C1 (i,j) = M1(i,j) / sqrt(M1 (i,i) * M1 (j,j));
      C0 (i,j) = M0(i,j) / sqrt(M0 (i,i) * M0 (j,j));

    }
  }
  
  std::cout << "...Error Matrix Filled" << std::endl;


  TCanvas* c6 = new TCanvas("c6","",700,700);
  c6->SetLeftMargin(.1);
  c6->SetBottomMargin(.1);
  c6->SetTopMargin(.075);
  c6->SetRightMargin(.15);
  c6->cd();

  M6.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  //  TMatrixFBase->GetZaxis()->SetRangeUser(-0.05,0.4);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //  TMatrixFBase->GetZaxis()->SetTitle("Fractional Error Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kBlue);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);



  TLatex *ND = new TLatex(.15,.01,"LAr1-ND (200m) ");
  ND->SetNDC();
  ND->SetTextFont(62);
  ND->SetTextSize(0.04);
  ND->Draw();

  TLatex *MD = new TLatex(.5,.01,"T600 (600m, on axis)");
  MD->SetNDC();
  MD->SetTextFont(62);
  MD->SetTextSize(0.04);
  MD->Draw();

  TLatex *ND45 = new TLatex(.05,.15,"LAr1-ND (200m) ");
  ND45->SetNDC();
  ND45->SetTextAngle(90);
  ND45->SetTextFont(62);
  ND45->SetTextSize(0.04);
  ND45->Draw();

  TLatex *MD45 = new TLatex(.05,.54,"T600 (600m, on axis)");
  MD45->SetNDC();
  MD45->SetTextAngle(90);
  MD45->SetTextFont(62);
  MD45->SetTextSize(0.04);
  MD45->Draw();

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} Flux Fractional Error Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();

  c6->Print("total_matrix_absolute.pdf");


  /*  TCanvas* c8 = new TCanvas("c8","",700,700);
  c8->SetLeftMargin(.1);
  c8->SetBottomMargin(.1);
  c8->SetTopMargin(.1);
  c8->SetRightMargin(.2);
  c8->cd();
  Fig5->SetLineWidth(4);
  Fig5->Draw("");*/
  


  TCanvas* c61 = new TCanvas("c61","",700,700);
  c61->SetLeftMargin(.1);
  c61->SetBottomMargin(.1);
  c61->SetTopMargin(.075);
  c61->SetRightMargin(.15);
  c61->cd();

  C6.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-1,1);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //  TMatrixFBase->GetZaxis()->SetTitle("Combined Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);
  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kYellow);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);


  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} Flux Correlation Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();

  c61->Print("total_correlation_matrix_absolute.pdf");

 

  TCanvas* c5 = new TCanvas("c5","",700,700);
  c5->SetLeftMargin(.1);
  c5->SetBottomMargin(.1);
  c5->SetTopMargin(.075);
  c5->SetRightMargin(.15);
  c5->cd();

  M5.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
//   TMatrixFBase->GetZaxis()->SetRangeUser(-0.005,0.045);
  
TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //  TMatrixFBase->GetZaxis()->SetTitle("K^{+} Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kBlue);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} K#lower[-0.15]{+} Fractional Error Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();




  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();

  c5->Print("mult5_matrix_absolute.pdf");

  TCanvas* c51 = new TCanvas("c51","",700,700);
  c51->SetLeftMargin(.1);
  c51->SetBottomMargin(.1);
  c51->SetTopMargin(.075);
  c51->SetRightMargin(.15);
  c51->cd();

  C5.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-1,1);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //  TMatrixFBase->GetZaxis()->SetTitle("K#lower[-0.15]{+} Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kYellow);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} K#lower[-0.15]{+} Correlation Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();



  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();

  c51->Print("mult5_correlation_matrix_absolute.pdf");

  TCanvas* c4 = new TCanvas("c4","",700,700);
  c4->SetLeftMargin(.1);
  c4->SetBottomMargin(.1);
  c4->SetTopMargin(.075);
  c4->SetRightMargin(.15);
  c4->cd();

  M4.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
//   TMatrixFBase->GetZaxis()->SetRangeUser(-0.005,0.045);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //TMatrixFBase->GetZaxis()->SetTitle("K#lower[-0.15]{-} Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kBlue);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} K#lower[-0.15]{-} Fractional Error Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();



  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c4->Print("mult4_matrix_absolute.pdf");

  TCanvas* c41 = new TCanvas("c41","",700,700);
  c41->SetLeftMargin(.1);
  c41->SetBottomMargin(.1);
  c41->SetTopMargin(.075);
  c41->SetRightMargin(.15);
  c41->cd();

  C4.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-1,1);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //  TMatrixFBase->GetZaxis()->SetTitle("K#lower[-0.15]{-} Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kYellow);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} K#lower[-0.15]{-} Correlation Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c41->Print("mult4_correlation_matrix_absolute.pdf");



  TCanvas* c3 = new TCanvas("c3","",700,700);
  c3->SetLeftMargin(.1);
  c3->SetBottomMargin(.1);
  c3->SetTopMargin(.075);
  c3->SetRightMargin(.15);
  c3->cd();

  M3.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
//   TMatrixFBase->GetZaxis()->SetRangeUser(-0.005,0.045);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //TMatrixFBase->GetZaxis()->SetTitle("K#lower[-0.15]{0} Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kBlue);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} K#lower[-0.15]{0} Fractional Error Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();


  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c3->Print("mult3_matrix_absolute.pdf");

  TCanvas* c31 = new TCanvas("c31","",700,700);
  c31->SetLeftMargin(.1);
  c31->SetBottomMargin(.1);
  c31->SetTopMargin(.075);
  c31->SetRightMargin(.15);
  c31->cd();

  C3.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-1,1);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //TMatrixFBase->GetZaxis()->SetTitle("K#lower[-0.15]{0} Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kYellow);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} K#lower[-0.15]{0} Correlation Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();



  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c31->Print("mult3_correlation_matrix_absolute.pdf");


  TCanvas* c2 = new TCanvas("c2","",700,700);
  c2->SetLeftMargin(.1);
  c2->SetBottomMargin(.1);
  c2->SetTopMargin(.075);
  c2->SetRightMargin(.15);
  c2->cd();

  M2.Draw("COLZ");
  gStyle->SetPalette(56,0);
//   TMatrixFBase->GetZaxis()->SetRangeUser(-0.005,0.045);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //TMatrixFBase->GetZaxis()->SetTitle("#pi#lower[-0.15]{+} Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kBlue);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} #pi#lower[-0.15]{+} Fractional Error Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();



  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c2->Print("mult2_matrix_absolute.pdf");

  TCanvas* c21 = new TCanvas("c21","",700,700);
  c21->SetLeftMargin(.1);
  c21->SetBottomMargin(.1);
  c21->SetTopMargin(.075);
  c21->SetRightMargin(.15);
  c21->cd();

  C2.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-1,1);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //TMatrixFBase->GetZaxis()->SetTitle("#pi#lower[-0.15]{+} Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kYellow);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} #pi#lower[-0.15]{+} Correlation Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c21->Print("mult2_correlation_matrix_absolute.pdf");


  TCanvas* c1 = new TCanvas("c1","",700,700);
  c1->SetLeftMargin(.1);
  c1->SetBottomMargin(.1);
  c1->SetTopMargin(.075);
  c1->SetRightMargin(.15);
  c1->cd();

  M1.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
//   TMatrixFBase->GetZaxis()->SetRangeUser(-0.005,0.045);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //TMatrixFBase->GetZaxis()->SetTitle("#pi#lower[-0.15]{-} Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kBlue);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} #pi#lower[-0.15]{-} Fractional Error Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c1->Print("mult1_matrix_absolute.pdf");

  TCanvas* c11 = new TCanvas("c11","",700,700);
  c11->SetLeftMargin(.1);
  c11->SetBottomMargin(.1);
  c11->SetTopMargin(.075);
  c11->SetRightMargin(.15);
  c11->cd();

  C1.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-1,1);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //  TMatrixFBase->GetZaxis()->SetTitle("#pi#lower[-0.15]{-} Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kYellow);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} #pi#lower[-0.15]{-} Correlation Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c11->Print("mult1_correlation_matrix_absolute.pdf");


  TCanvas* c0 = new TCanvas("c0","",700,700);
  c0->SetLeftMargin(.1);
  c0->SetBottomMargin(.1);
  c0->SetTopMargin(.075);
  c0->SetRightMargin(.15);
  c0->cd();

  M0.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
//   TMatrixFBase->GetZaxis()->SetRangeUser(-0.005,0.045);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //  TMatrixFBase->GetZaxis()->SetTitle("Beam UniSim Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kBlue);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} Beam Fractional Error Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c0->Print("mult0_matrix_absolute.pdf");

  TCanvas* c01 = new TCanvas("c01","",700,700);
  c01->SetLeftMargin(.1);
  c01->SetBottomMargin(.1);
  c01->SetTopMargin(.075);
  c01->SetRightMargin(.15);
  c01->cd();

  C0.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-1,1);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  //  TMatrixFBase->GetZaxis()->SetTitle("Beam UniSim Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kYellow);
  split->DrawLineNDC(.1,.51,.849,.51);
  split->DrawLineNDC(.475,.101,.475,.930);
  add_plot_label("|            0.2 #minus 3.0 GeV            |            0.2 #minus 3.0 GeV            | ", 0.48,0.08,0.03);

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} Beam Correlation Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();

  c01->Print("mult0_correlation_matrix_absolute.pdf");
 
    TCanvas* c86 = new TCanvas("c86","",800,400);
    c86->SetLeftMargin(.1);
    c86->SetBottomMargin(.1);
    c86->SetTopMargin(.05);
    c86->SetRightMargin(.05);
    c86->SetLogy();
    c86->cd();
    
    Fig6->GetYaxis()->SetTitle("Fractional Error");
    Fig6->GetYaxis()->SetTitleFont(62);
    Fig6->GetXaxis()->SetTitleFont(62);
    Fig6->GetYaxis()->SetLabelFont(62);
    Fig6->GetXaxis()->SetLabelFont(62);
    Fig6->GetYaxis()->CenterTitle();
    Fig6->GetYaxis()->SetTitleSize(0.06);
    Fig6->GetYaxis()->SetTitleOffset(0.8);
    Fig6->GetXaxis()->SetLabelSize(0.06);
    Fig6->GetYaxis()->SetLabelSize(0.06);
    Fig6->GetXaxis()->SetTitleOffset(1.5);
    Fig6->SetStats(0);
    Fig6->SetMinimum(0.01);
    Fig6->SetMaximum(4);
    Fig6->SetMarkerStyle(8);
    Fig6->GetYaxis()->SetNdivisions(509);
    Fig6->GetXaxis()->SetNdivisions(509);
    Fig6->Draw("P");
    split->SetLineColor(1);
    split->SetLineWidth(2);
    split->DrawLine(19,0.01,19,4);
    
    TLatex *ND = new TLatex(.23,.85,"LAr1-ND (200m) ");
    ND->SetNDC();
    ND->SetTextFont(62);
    ND->SetTextSize(0.05);
    ND->Draw();
    
    TLatex *MD = new TLatex(.65,.85,"T600 (600m, on axis)");
    MD->SetNDC();
    MD->SetTextFont(62);
    MD->SetTextSize(0.05);
    MD->Draw();
    
    c86->Print("FractionalErrors_Total_absolute.pdf");
    
    TCanvas* c85 = new TCanvas("c85","",800,400);
    c85->SetLeftMargin(.1);
    c85->SetBottomMargin(.1);
    c85->SetTopMargin(.05);
    c85->SetRightMargin(.05);
    //    c85->SetLogy();
    c85->cd();
    
    Fig5->GetYaxis()->SetTitle("K#lower[-0.2]{+} Fractional Error");
    Fig5->GetYaxis()->SetTitleFont(62);
    Fig5->GetXaxis()->SetTitleFont(62);
    Fig5->GetYaxis()->SetLabelFont(62);
    Fig5->GetXaxis()->SetLabelFont(62);
    Fig5->GetYaxis()->CenterTitle();
    Fig5->GetYaxis()->SetTitleSize(0.06);
    Fig5->GetYaxis()->SetTitleOffset(0.8);
    Fig5->GetXaxis()->SetLabelSize(0.06);
    Fig5->GetYaxis()->SetLabelSize(0.06);
    Fig5->GetXaxis()->SetTitleOffset(1.5);
    Fig5->SetStats(0);
    Fig5->SetMinimum(0.01);
    Fig5->SetMaximum(4);
    Fig5->SetMarkerStyle(8);
    Fig5->GetYaxis()->SetNdivisions(509);
    Fig5->GetXaxis()->SetNdivisions(509);
    Fig5->Draw("P");
    split->SetLineColor(1);
    split->SetLineWidth(2);
    split->DrawLine(19,-0.01,19,0.21);
    ND->Draw();
    MD->Draw();
    
    //    c85->Print("FractionalErrors_Kplus_absolute.pdf");
    
    
    TCanvas* c84 = new TCanvas("c84","",800,400);
    c84->SetLeftMargin(.1);
    c84->SetBottomMargin(.1);
    c84->SetTopMargin(.05);
    c84->SetRightMargin(.05);
    c84->SetLogy();
    c84->cd();
    
    Fig4->GetYaxis()->SetTitle("K#lower[-0.2]{-} Fractional Error");
    Fig4->GetYaxis()->SetTitleFont(62);
    Fig4->GetXaxis()->SetTitleFont(62);
    Fig4->GetYaxis()->SetLabelFont(62);
    Fig4->GetXaxis()->SetLabelFont(62);
    Fig4->GetYaxis()->CenterTitle();
    Fig4->GetYaxis()->SetTitleSize(0.06);
    Fig4->GetYaxis()->SetTitleOffset(0.8);
    Fig4->GetXaxis()->SetLabelSize(0.06);
    Fig4->GetYaxis()->SetLabelSize(0.06);
    Fig4->GetXaxis()->SetTitleOffset(1.5);
    Fig4->SetStats(0);
    Fig4->SetMinimum(0.01);
    Fig4->SetMaximum(4);
    Fig4->SetMarkerStyle(8);
    Fig4->GetYaxis()->SetNdivisions(509);
    Fig4->GetXaxis()->SetNdivisions(509);
    Fig4->Draw("P");
    split->SetLineColor(1);
    split->SetLineWidth(2);
    split->DrawLine(19,0.01,19,4);
    ND->Draw();
    MD->Draw();
    
    c84->Print("FractionalErrors_Kmin_absolute.pdf");
    
    
    TCanvas* c83 = new TCanvas("c83","",800,400);
    c83->SetLeftMargin(.1);
    c83->SetBottomMargin(.1);
    c83->SetTopMargin(.05);
    c83->SetRightMargin(.05);
    c83->SetLogy();
    c83->cd();
    
    Fig3->GetYaxis()->SetTitle("K#lower[-0.2]{0} Fractional Error");
    Fig3->GetYaxis()->SetTitleFont(62);
    Fig3->GetXaxis()->SetTitleFont(62);
    Fig3->GetYaxis()->SetLabelFont(62);
    Fig3->GetXaxis()->SetLabelFont(62);
    Fig3->GetYaxis()->CenterTitle();
    Fig3->GetYaxis()->SetTitleSize(0.06);
    Fig3->GetYaxis()->SetTitleOffset(0.8);
    Fig3->GetXaxis()->SetLabelSize(0.06);
    Fig3->GetYaxis()->SetLabelSize(0.06);
    Fig3->GetXaxis()->SetTitleOffset(1.5);
    Fig3->SetStats(0);
    Fig3->SetMinimum(0.01);
    Fig3->SetMaximum(4);
    Fig3->SetMarkerStyle(8);
    Fig3->GetYaxis()->SetNdivisions(509);
    Fig3->GetXaxis()->SetNdivisions(509);
    Fig3->Draw("P");
    split->SetLineColor(1);
    split->SetLineWidth(2);
    split->DrawLine(19,0.01,19,4);
    ND->Draw();
    MD->Draw();
    
    c83->Print("FractionalErrors_K0_absolute.pdf");
    
    
    TCanvas* c82 = new TCanvas("c82","",800,400);
    c82->SetLeftMargin(.1);
    c82->SetBottomMargin(.1);
    c82->SetTopMargin(.05);
    c82->SetRightMargin(.05);
    c82->SetLogy();
    c82->cd();
    
    Fig2->GetYaxis()->SetTitle("#pi#lower[-0.2]{+} Fractional Error");
    Fig2->GetYaxis()->SetTitleFont(62);
    Fig2->GetXaxis()->SetTitleFont(62);
    Fig2->GetYaxis()->SetLabelFont(62);
    Fig2->GetXaxis()->SetLabelFont(62);
    Fig2->GetYaxis()->CenterTitle();
    Fig2->GetYaxis()->SetTitleSize(0.06);
    Fig2->GetYaxis()->SetTitleOffset(0.8);
    Fig2->GetXaxis()->SetLabelSize(0.06);
    Fig2->GetYaxis()->SetLabelSize(0.06);
    Fig2->GetXaxis()->SetTitleOffset(1.5);
    Fig2->SetStats(0);
    Fig2->SetMinimum(0.01);
    Fig2->SetMaximum(4);
    Fig2->SetMarkerStyle(8);
    Fig2->GetYaxis()->SetNdivisions(509);
    Fig2->GetXaxis()->SetNdivisions(509);
    Fig2->Draw("P");
    split->SetLineColor(1);
    split->SetLineWidth(2);
    split->DrawLine(19,0.01,19,4);
    ND->Draw();
    MD->Draw();
    
    c82->Print("FractionalErrors_piplus_absolute.pdf");
    
    
    TCanvas* c81 = new TCanvas("c81","",800,400);
    c81->SetLeftMargin(.1);
    c81->SetBottomMargin(.1);
    c81->SetTopMargin(.05);
    c81->SetRightMargin(.05);
    c81->SetLogy();
    c81->cd();
    
    Fig1->GetYaxis()->SetTitle("#pi#lower[-0.2]{-} Fractional Error");
    Fig1->GetYaxis()->SetTitleFont(62);
    Fig1->GetXaxis()->SetTitleFont(62);
    Fig1->GetYaxis()->SetLabelFont(62);
    Fig1->GetXaxis()->SetLabelFont(62);
    Fig1->GetYaxis()->CenterTitle();
    Fig1->GetYaxis()->SetTitleSize(0.06);
    Fig1->GetYaxis()->SetTitleOffset(0.8);
    Fig1->GetXaxis()->SetLabelSize(0.06);
    Fig1->GetYaxis()->SetLabelSize(0.06);
    Fig1->GetXaxis()->SetTitleOffset(1.5);
    Fig1->SetStats(0);
    Fig1->SetMinimum(0.01);
    Fig1->SetMaximum(4);
    Fig1->SetMarkerStyle(8);
    Fig1->GetYaxis()->SetNdivisions(509);
    Fig1->GetXaxis()->SetNdivisions(509);
    Fig1->Draw("P");
    split->SetLineColor(1);
    split->SetLineWidth(2);
    split->DrawLine(19,0.01,19,4);
    ND->Draw();
    MD->Draw();
    
    c81->Print("FractionalErrors_pimin_absolute.pdf");
    
    
    TCanvas* c80 = new TCanvas("c80","",800,400);
    c80->SetLeftMargin(.1);
    c80->SetBottomMargin(.1);
    c80->SetTopMargin(.05);
    c80->SetRightMargin(.05);
    c80->SetLogy();
    c80->cd();
    
    Fig0->GetYaxis()->SetTitle("Beam Fractional Error");
    Fig0->GetYaxis()->SetTitleFont(62);
    Fig0->GetXaxis()->SetTitleFont(62);
    Fig0->GetYaxis()->SetLabelFont(62);
    Fig0->GetXaxis()->SetLabelFont(62);
    Fig0->GetYaxis()->CenterTitle();
    Fig0->GetYaxis()->SetTitleSize(0.06);
    Fig0->GetYaxis()->SetTitleOffset(0.8);
    Fig0->GetXaxis()->SetLabelSize(0.06);
    Fig0->GetYaxis()->SetLabelSize(0.06);
    Fig0->GetXaxis()->SetTitleOffset(1.5);
    Fig0->SetStats(0);
    Fig0->SetMinimum(0.01);
    Fig0->SetMaximum(4);
    Fig0->SetMarkerStyle(8);
    Fig0->GetYaxis()->SetNdivisions(509);
    Fig0->GetXaxis()->SetNdivisions(509);
    Fig0->Draw("P");
    split->SetLineColor(1);
    split->SetLineWidth(2);
    split->DrawLine(19,0.01,19,4);
    ND->Draw();
    MD->Draw();
    
    c80->Print("FractionalErrors_beam_absolute.pdf");

    TCanvas* c9 = new TCanvas("c9","",800,400);
    c9->SetLeftMargin(.1);
    c9->SetBottomMargin(.15);
    c9->SetTopMargin(.05);
    c9->SetRightMargin(.05);
    c9->SetLogy();
    c9->cd();

    
    ND_err->GetYaxis()->SetTitle("Fractional Uncertainty [%]");
    ND_err->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    ND_err->GetYaxis()->SetTitleFont(62);
    ND_err->GetXaxis()->SetTitleFont(62);
    ND_err->GetYaxis()->SetLabelFont(62);
    ND_err->GetXaxis()->SetLabelFont(62);
    ND_err->GetYaxis()->CenterTitle();
    ND_err->GetYaxis()->SetTitleSize(0.06);
    ND_err->GetXaxis()->SetTitleSize(0.6);
    ND_err->GetXaxis()->SetLabelSize(0.06);
    ND_err->GetYaxis()->SetLabelSize(0.06);
    ND_err->GetYaxis()->SetTitleOffset(0.8);
    ND_err->GetXaxis()->SetTitleOffset(1.5);
    ND_err->SetStats(0);
    ND_err->SetMinimum(1.1);
    ND_err->SetMaximum(100);
    ND_err->GetYaxis()->SetNdivisions(509);
    ND_err->GetXaxis()->SetNdivisions(509);
    
    ND_err->SetLineWidth(4);
    ND_err->SetLineColor(kBlue-3);
    ND_err->Draw("h");
    FD_err->SetLineWidth(3);
    FD_err->SetLineStyle(2);
    FD_err->SetLineColor(kRed-3);
    FD_err->Draw("h same");
    
    TLegend* leg3=new TLegend(0.2,0.7,0.4,0.8);
    leg3->SetFillStyle(0);
    leg3->SetFillColor(0);
    leg3->SetBorderSize(0);
    leg3->SetTextFont(62);
    leg3->SetTextSize(0.03);
    leg3->AddEntry(ND_err,"LAr1-ND (200m)","L");
    leg3->AddEntry(FD_err,"T600 (600m, on axis)","L");
    leg3->Draw();
 

    TCanvas* c9 = new TCanvas("c9","",800,400);
    c9->SetLeftMargin(.1);
    c9->SetBottomMargin(.15);
    c9->SetTopMargin(.05);
    c9->SetRightMargin(.05);
    c9->SetLogy();
    c9->cd();

    ND_err_piplus->GetYaxis()->SetTitle("#pi#lower[-0.4]{-} Fractional Uncertainty [%]");
    ND_err_piplus->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    ND_err_piplus->GetYaxis()->SetTitleFont(62);
    ND_err_piplus->GetXaxis()->SetTitleFont(62);
    ND_err_piplus->GetYaxis()->SetLabelFont(62);
    ND_err_piplus->GetXaxis()->SetLabelFont(62);
    ND_err_piplus->GetYaxis()->CenterTitle();
    ND_err_piplus->GetYaxis()->SetTitleSize(0.06);
    ND_err_piplus->GetXaxis()->SetTitleSize(0.6);
    ND_err_piplus->GetXaxis()->SetLabelSize(0.06);
    ND_err_piplus->GetYaxis()->SetLabelSize(0.06);
    ND_err_piplus->GetYaxis()->SetTitleOffset(0.8);
    ND_err_piplus->GetXaxis()->SetTitleOffset(1.5);
    ND_err_piplus->SetStats(0);
    ND_err_piplus->SetMinimum(1.1);
    ND_err_piplus->SetMaximum(100);
    ND_err_piplus->GetYaxis()->SetNdivisions(509);
    ND_err_piplus->GetXaxis()->SetNdivisions(509);
    
    ND_err_piplus->SetLineWidth(4);
    ND_err_piplus->SetLineColor(kBlue-3);
    ND_err_piplus->Draw("h");
    FD_err_piplus->SetLineWidth(3);
    FD_err_piplus->SetLineStyle(2);
    FD_err_piplus->SetLineColor(kRed-3);
    FD_err_piplus->Draw("h same");
    leg3->Draw();

  
  cout<<"\nEnd of routine.\n";

  return 0;

}


#ifndef __CINT__
int main()
{
  multiple_detector_fit();
  return 0;
}
# endif

void Matrix_absolute(){
  multiple_detector_fit();
  return;
}

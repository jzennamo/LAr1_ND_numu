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

  bool use100m = false;
  bool use470m = true;
  bool use600m = false;


  std::vector<int> baselines;
  std::vector<double> scales;
  std::vector<std::string> names;
  std::vector<double> volume;

  if (use100m) baselines.push_back(100);
  if (use470m) baselines.push_back(470);
  if (use600m) baselines.push_back(600);
  int nL = baselines.size();

  double NULLVec[3][20];
  double OscVec[3][1001][7][20];

  for(int i = 0; i < 20; i++){
    NULLVec[0][i] = 0;
    NULLVec[1][i] = 0;
    NULLVec[2][i] = 0;

  }

  for(int u = 0; u < 1000; u++){
    for(int s = 0; s < 7; s++){
      for(int i = 0; i < 20; i++){

	OscVec[0][u][s][i] = 0;
        OscVec[1][u][s][i] = 0;
        OscVec[2][u][s][i] = 0;

      }
    }
  }


  int nbinsE = 0;

  int counter = 0;
  if (use100m){

    std::string temp_name = "../MatrixFiles/combined_ntuple_100m_nu_processed_numu.root";

    TFile temp_file(temp_name.c_str());
    TH1D *NULL_100;
    NULL_100 = (TH1D*)(temp_file.Get("NumuCC"));
    nbinsE = NULL_100->GetNbinsX();
    std::cout << nbinsE << std::endl;
      for(int i = 1; i <= nbinsE; i++){
	  NULLVec[counter][i-1] = (NULL_100->GetBinContent(i));
      }

    for(int u = 0; u < npoints; u++){
      for(int s = 0; s < 7; s++){
	TH1D *OSC_100;
	TString upoint = Form("%d",u);
	TString name = "Universe_";
	TString name2 = "_MultiSim_";
	TString mul = Form("%d",s);
	
	name += upoint;
	name += name2;
	name += mul;	
       
	OSC_100 = (TH1D*)(temp_file.Get(name));
          for(int i = 1; i <= nbinsE; i++){
              OscVec[counter][u][s][i-1] = (OSC_100->GetBinContent(i));
	      //	      if(OscVec[0][u][s][i-1] != OscVec[0][u][s][i-1]) std::cout << "erm" <<std::endl;

	  }

	delete OSC_100;
      }
    }
    counter++;
    delete NULL_100;
    temp_file.Close();
  }

  if (use470m){
    std::string temp_name = "../MatrixFiles/combined_ntuple_470m_nu_processed_numu.root";

    TFile temp_file(temp_name.c_str());
    TH1D *NULL_470;
    NULL_470 = (TH1D*)(temp_file.Get("NumuCC"));
    nbinsE = NULL_470->GetNbinsX();
    std::cout << nbinsE<< std::endl;
      for(int i = 1; i <= nbinsE; i++){
	  NULLVec[counter][i-1] = (NULL_470->GetBinContent(i));
      }

    for(int u = 0; u < npoints; u++){
      for(int s = 0; s < 7; s++){
	TH1D *OSC_470;
	TString upoint = Form("%d",u);//std::to_string(u);
	TString name = "Universe_";
	TString name2 = "_MultiSim_";
	TString mul = Form("%d",s);// = std::to_string(s);
	
	name += upoint;
	name += name2;
	name += mul;	
       
	OSC_470 = (TH1D*)(temp_file.Get(name));
          for(int i = 1; i <= nbinsE; i++){
	    OscVec[counter][u][s][i-1] = (OSC_470->GetBinContent(i));
	  }

	delete OSC_470;
      }
    }
    counter++;
    delete NULL_470;
    temp_file.Close();
  }

  if (use600m){
    std::string temp_name = "../MatrixFiles/combined_ntuple_600m_onaxis_nu_processed_numu.root";

    TFile temp_file(temp_name.c_str());
    TH1D *NULL_600;
    NULL_600 = (TH1D*)(temp_file.Get("NumuCC"));
    nbinsE = NULL_600->GetNbinsX();
    std::cout << nbinsE<< std::endl;
      for(int i = 1; i <= nbinsE; i++){
	  NULLVec[counter][i-1] = (NULL_600->GetBinContent(i));
      }

    for(int u = 0; u < npoints; u++){
      for(int s = 0; s < 7; s++){
	TH1D *OSC_600;
	TString upoint = Form("%d",u);//std::to_string(u);
	TString name = "Universe_";
	TString name2 = "_MultiSim_";
	TString mul = Form("%d",s);// = std::to_string(s);
	
	name += upoint;
	name += name2;
	name += mul;	
       
	OSC_600 = (TH1D*)(temp_file.Get(name));
          for(int i = 1; i <= nbinsE; i++){
	    OscVec[counter][u][s][i-1] = (OSC_600->GetBinContent(i));

	  }

	delete OSC_600;
      }
    }
    counter++;    
    delete NULL_600;
    temp_file.Close();
  }


  //  int nL = 3;
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

  for(int Lrow = 0; Lrow < nL; Lrow++){
    for(int Erow = 0; Erow < nbinsE; Erow++){

      Errj = 0;

      for(int Lcol = 0; Lcol < nL; Lcol++){
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
            M5 (Erri,Errj) += (NULLVec[Lrow][Erow]-OscVec[Lrow][u][5][Erow])*(NULLVec[Lcol][Ecol]-OscVec[Lcol][u][5][Ecol]);
            M4 (Erri,Errj) += (NULLVec[Lrow][Erow]-OscVec[Lrow][u][4][Erow])*(NULLVec[Lcol][Ecol]-OscVec[Lcol][u][4][Ecol]);
            M3 (Erri,Errj) += (NULLVec[Lrow][Erow]-OscVec[Lrow][u][3][Erow])*(NULLVec[Lcol][Ecol]-OscVec[Lcol][u][3][Ecol]);
            M2 (Erri,Errj) += (NULLVec[Lrow][Erow]-OscVec[Lrow][u][2][Erow])*(NULLVec[Lcol][Ecol]-OscVec[Lcol][u][2][Ecol]);
            M1 (Erri,Errj) += (NULLVec[Lrow][Erow]-OscVec[Lrow][u][1][Erow])*(NULLVec[Lcol][Ecol]-OscVec[Lcol][u][1][Ecol]);
            M0 (Erri,Errj) += (NULLVec[Lrow][Erow]-OscVec[Lrow][u][0][Erow])*(NULLVec[Lcol][Ecol]-OscVec[Lcol][u][0][Ecol]);

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
          M5 (Erri,Errj) /= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];
          M4 (Erri,Errj) /= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];
          M3 (Erri,Errj) /= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];
          M2 (Erri,Errj) /= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];
          M1 (Erri,Errj) /= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];
          M0 (Erri,Errj) /= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];
	  
	  if(Erri == Errj) Fig6->SetBinContent(Erri+1, sqrt(M6 (Erri,Errj)));
      if(Erri == Errj) Fig5->SetBinContent(Erri+1, sqrt(M5 (Erri,Errj)));
      if(Erri == Errj) Fig4->SetBinContent(Erri+1, sqrt(M4 (Erri,Errj)));
      if(Erri == Errj) Fig3->SetBinContent(Erri+1, sqrt(M3 (Erri,Errj)));
      if(Erri == Errj) Fig2->SetBinContent(Erri+1, sqrt(M2 (Erri,Errj)));
      if(Erri == Errj) Fig1->SetBinContent(Erri+1, sqrt(M1 (Erri,Errj)));
      if(Erri == Errj) Fig0->SetBinContent(Erri+1, sqrt(M0 (Erri,Errj)));


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
  split->DrawLineNDC(.1,.375,.849,.375);
  split->DrawLineNDC(.1,.65,.849,.65);
  split->DrawLineNDC(.6,.101,.6,.930);
  split->DrawLineNDC(.35,.101,.35,.930);

  add_plot_label("|    0.2 #minus 3.0 GeV    |    0.2 #minus 3.0 GeV    |    0.2 #minus 3.0 GeV    | ", 0.48,0.08,0.03);



  TLatex *ND = new TLatex(.12,.01,"LAr1-ND (100m) ");
  ND->SetNDC();
  ND->SetTextFont(62);
  ND->SetTextSize(0.0275);
  ND->Draw();

  TLatex *MD = new TLatex(.355,.01,"MicroBooNE (470m)");
  MD->SetNDC();
  MD->SetTextFont(62);
  MD->SetTextSize(0.0275);
  MD->Draw();

  TLatex *FD = new TLatex(.65,.01,"T600 (600m)");
  FD->SetNDC();
  FD->SetTextFont(62);
  FD->SetTextSize(0.0275);
  FD->Draw();

  TLatex *ND45 = new TLatex(.05,.1,"LAr1-ND (100m) ");
  ND45->SetNDC();
  ND45->SetTextAngle(90);
  ND45->SetTextFont(62);
  ND45->SetTextSize(0.03);
  ND45->Draw();

  TLatex *MD45 = new TLatex(.05,.37,"MicroBooNE (470m)");
  MD45->SetNDC();
  MD45->SetTextAngle(90);
  MD45->SetTextFont(62);
  MD45->SetTextSize(0.03);
  MD45->Draw();

  TLatex *FD45 = new TLatex(.05,.7,"T600 (600m)");
  FD45->SetNDC();
  FD45->SetTextAngle(90);
  FD45->SetTextFont(62);
  FD45->SetTextSize(0.03);
  FD45->Draw();



  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} Flux Fractional Error Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();

  c6->Print("total_3Det_matrix.pdf");



  TCanvas* c61 = new TCanvas("c61","",700,700);
  c61->SetLeftMargin(.1);
  c61->SetBottomMargin(.1);
  c61->SetTopMargin(.075);
  c61->SetRightMargin(.15);
  c61->cd();

  C6.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
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
  split->DrawLineNDC(.1,.375,.849,.375);
  split->DrawLineNDC(.1,.65,.849,.65);
  split->DrawLineNDC(.6,.101,.6,.930);
  split->DrawLineNDC(.35,.101,.35,.930);

  add_plot_label("|    0.2 #minus 3.0 GeV    |    0.2 #minus 3.0 GeV    |    0.2 #minus 3.0 GeV    | ", 0.48,0.08,0.03);
 
  ND->Draw();
  MD->Draw();
  FD->Draw();
  ND45->Draw();
  MD45->Draw();
  FD45->Draw();

  TLatex *Total = new TLatex(.2,.96,"#nu#lower[0.3]{#mu} Flux Correlation Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.045);
  Total->Draw();

  c61->Print("total_3Det_correlation_matrix.pdf");

 
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

void Matrix_3Det(){
  multiple_detector_fit();
  return;
}

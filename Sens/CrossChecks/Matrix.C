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
  double OscVec[2][1001][7][20];

  for(int i = 0; i < 20; i++){
    NULLVec[0][i] = 0;
    NULLVec[1][i] = 0;

  }

  for(int u = 0; u < 1000; u++){
    for(int s = 0; s < 7; s++){
      for(int i = 0; i < 20; i++){

	OscVec[0][u][s][i] = 0;
        OscVec[1][u][s][i] = 0;

      }
    }
  }


  int nbinsE = 0;


  if (use100m){

    std::string temp_name = "../MatrixFiles/combined_ntuple_100m_nu_processed_numu.root";

    TFile temp_file(temp_name.c_str());
    TH1D *NULL_100;
    NULL_100 = (TH1D*)(temp_file.Get("NumuCC"));
    nbinsE = NULL_100->GetNbinsX();
      for(int i = 1; i <= nbinsE; i++){
	  NULLVec[0][i-1] = (NULL_100->GetBinContent(i));
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
              OscVec[0][u][s][i-1] = (OSC_100->GetBinContent(i));

	  }

	delete OSC_100;
      }
    }
    
    delete NULL_100;
    temp_file.Close();
  }

  if (use470m){
    std::string temp_name = "../MatrixFiles/combined_ntuple_600m_onaxis_nu_processed_numu.root";

    TFile temp_file(temp_name.c_str());
    TH1D *NULL_470;
    NULL_470 = (TH1D*)(temp_file.Get("NumuCC"));
    nbinsE = NULL_470->GetNbinsX();
      for(int i = 1; i <= nbinsE; i++){
	  NULLVec[1][i-1] = (NULL_470->GetBinContent(i));
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
              OscVec[1][u][s][i-1] = (OSC_470->GetBinContent(i));

	  }

	delete OSC_470;
      }
    }
    
    delete NULL_470;
    temp_file.Close();
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

  TH1D *Fig5 = new TH1D("Fig5",";;",nbinsE,emin,emax);

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

	  if(Lcol == 0 && Lrow == 0 && (Erri == Errj)) Fig5->SetBinContent(Ecol, M6 (Erri,Errj));

	  std::cout << M6 (Erri,Errj) << "\t";

          Errj++;

	}}

      Erri++;

    }}
  std::endl;

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
  c6->SetTopMargin(.1);
  c6->SetRightMargin(.2);
  c6->cd();

  M6.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  //  TMatrixFBase->GetZaxis()->SetRangeUser(-0.05,0.4);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("Fractional Error Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);


  TLatex *ND = new TLatex(.2,.05,"LAr1-ND (100m) ");
  ND->SetNDC();
  ND->SetTextFont(62);
  ND->SetTextSize(0.05);
  ND->Draw();

  TLatex *MD = new TLatex(.57,.05,"T600 (on axis)");
  MD->SetNDC();
  MD->SetTextFont(62);
  MD->SetTextSize(0.05);
  MD->Draw();

  TLatex *ND45 = new TLatex(.075,.2,"LAr1-ND (100m) ");
  ND45->SetNDC();
  ND45->SetTextAngle(90);
  ND45->SetTextFont(62);
  ND45->SetTextSize(0.05);
  ND45->Draw();

  TLatex *MD45 = new TLatex(.075,.67,"T600 (on axis)");
  MD45->SetNDC();
  MD45->SetTextAngle(90);
  MD45->SetTextFont(62);
  MD45->SetTextSize(0.05);
  MD45->Draw();

  TLatex *Total = new TLatex(.2,.95,"Combined Covariance Matrix");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.05);


  c6->Print("total_matrix.pdf");

  TCanvas* c61 = new TCanvas("c61","",700,700);
  c61->SetLeftMargin(.1);
  c61->SetBottomMargin(.1);
  c61->SetTopMargin(.1);
  c61->SetRightMargin(.2);
  c61->cd();

  C6.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  //  TMatrixFBase->GetZaxis()->SetRangeUser(-0.6,1);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("Combined Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c61->Print("total_correlation_matrix.pdf");

  /*

  TCanvas* c5 = new TCanvas("c5","",700,700);
  c5->SetLeftMargin(.1);
  c5->SetBottomMargin(.1);
  c5->SetTopMargin(.1);
  c5->SetRightMargin(.2);
  c5->cd();

  M5.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.05,0.4);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("K^{+} Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  //  c5->Print("../Figs/mult5_matrix.pdf");

  TCanvas* c51 = new TCanvas("c51","",700,700);
  c51->SetLeftMargin(.1);
  c51->SetBottomMargin(.1);
  c51->SetTopMargin(.1);
  c51->SetRightMargin(.2);
  c51->cd();

  C5.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.6,1);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("K^{+} Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  //  c51->Print("../Figs/mult5_correlation_matrix.pdf");

  TCanvas* c4 = new TCanvas("c4","",700,700);
  c4->SetLeftMargin(.1);
  c4->SetBottomMargin(.1);
  c4->SetTopMargin(.1);
  c4->SetRightMargin(.2);
  c4->cd();

  M4.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.05,0.4);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("K^{-} Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  //  c4->Print("../Figs/mult4_matrix.pdf");

  TCanvas* c41 = new TCanvas("c41","",700,700);
  c41->SetLeftMargin(.1);
  c41->SetBottomMargin(.1);
  c41->SetTopMargin(.1);
  c41->SetRightMargin(.2);
  c41->cd();

  C4.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.6,1);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("K^{-} Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  //  c41->Print("../Figs/mult4_correlation_matrix.pdf");



  TCanvas* c3 = new TCanvas("c3","",700,700);
  c3->SetLeftMargin(.1);
  c3->SetBottomMargin(.1);
  c3->SetTopMargin(.1);
  c3->SetRightMargin(.2);
  c3->cd();

  M3.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.05,0.4);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("K^{0} Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  //  c3->Print("../Figs/mult3_matrix.pdf");

  TCanvas* c31 = new TCanvas("c31","",700,700);
  c31->SetLeftMargin(.1);
  c31->SetBottomMargin(.1);
  c31->SetTopMargin(.1);
  c31->SetRightMargin(.2);
  c31->cd();

  C3.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.6,1);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("K^{0} Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  //  c31->Print("../Figs/mult3_correlation_matrix.pdf");


  TCanvas* c2 = new TCanvas("c2","",700,700);
  c2->SetLeftMargin(.1);
  c2->SetBottomMargin(.1);
  c2->SetTopMargin(.1);
  c2->SetRightMargin(.2);
  c2->cd();

  M2.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.05,0.4);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("#pi^{+} Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  //  c2->Print("../Figs/mult2_matrix.pdf");

  TCanvas* c21 = new TCanvas("c21","",700,700);
  c21->SetLeftMargin(.1);
  c21->SetBottomMargin(.1);
  c21->SetTopMargin(.1);
  c21->SetRightMargin(.2);
  c21->cd();

  C2.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.6,1);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("#pi^{+} Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c21->Print("../Figs/mult2_correlation_matrix.pdf");


  TCanvas* c1 = new TCanvas("c1","",700,700);
  c1->SetLeftMargin(.1);
  c1->SetBottomMargin(.1);
  c1->SetTopMargin(.1);
  c1->SetRightMargin(.2);
  c1->cd();

  M1.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.05,0.4);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("#pi^{-} Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c1->Print("../Figs/mult1_matrix.pdf");

  TCanvas* c11 = new TCanvas("c11","",700,700);
  c11->SetLeftMargin(.1);
  c11->SetBottomMargin(.1);
  c11->SetTopMargin(.1);
  c11->SetRightMargin(.2);
  c11->cd();

  C1.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.6,1);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("#pi^{-} Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c11->Print("../Figs/mult1_correlation_matrix.pdf");


  TCanvas* c0 = new TCanvas("c0","",700,700);
  c0->SetLeftMargin(.1);
  c0->SetBottomMargin(.1);
  c0->SetTopMargin(.1);
  c0->SetRightMargin(.2);
  c0->cd();

  M0.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.05,0.4);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("Beam UniSim Covariance Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();


  c0->Print("../Figs/mult0_matrix.pdf");

  TCanvas* c01 = new TCanvas("c01","",700,700);
  c01->SetLeftMargin(.1);
  c01->SetBottomMargin(.1);
  c01->SetTopMargin(.1);
  c01->SetRightMargin(.2);
  c01->cd();

  C0.Draw("COLZ");
  gStyle->SetPalette(56,0);
  TMatrixFBase->SetContour(999);
  TMatrixFBase->GetZaxis()->SetRangeUser(-0.6,1);

  TMatrixFBase->GetZaxis()->SetTitleFont(62);
  TMatrixFBase->GetZaxis()->SetLabelFont(62);
  TMatrixFBase->GetZaxis()->SetTitleSize(0.045);
  TMatrixFBase->GetZaxis()->SetTitle("Beam UniSim Correlation Matrix");
  TMatrixFBase->GetZaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetXaxis()->SetTitle("");
  TMatrixFBase->GetXaxis()->SetLabelSize(0);
  TMatrixFBase->GetXaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetTitle("");
  TMatrixFBase->GetYaxis()->SetTitleOffset(1.5);
  TMatrixFBase->GetYaxis()->SetLabelSize(0);
  TMatrixFBase->SetStats(0);

  ND->Draw();
  MD->Draw();
  ND45->Draw();
  MD45->Draw();

  c01->Print("../Figs/mult0_correlation_matrix.pdf");
  */
 
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

void Matrix(){
  multiple_detector_fit();
  return;
}

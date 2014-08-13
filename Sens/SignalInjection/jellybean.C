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
#include "TRandom3.h"
#include "TGraph2D.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TROOT.h" //for gROOT
#include "TStyle.h" //for gStyle
#include "TH2D.h"
#include "TLegend.h"
#include "THStack.h"
#include "TImage.h"
#include "TMarker.h"
#include "TLatex.h"
#include <iostream>
#include <ostream>
#include <sstream>
#include <fstream>

std::vector< std::vector <float> >   NULLVec;				// Null hypothesis
std::vector< std::vector < std::vector <float> > >   OscVec;		// Oscillated hypothesis for each dm2 point
double NomVec[3][20];
double SystVec[3][1001][7][20];
TMatrixT <float> cov;

TGraph * fiveSigma;
TGraph * threeSigma;
TGraph * ninety;
TMarker * signalPt;
TMarker * bestFit;

TGraph2D * chisqSurf;
TH1D * cvHist;
TH1D * predHist;
TH1D * predScaledHist;

TH1D *smeared;
TH1D *unsmeared;

double predSin22th;
double predDm2;
bool drawSurf;

//plot properties
Int_t npoints = 500;			// Grid Pointsl
Int_t spoints = 1000;			// Universes

//grid boundaries
Double_t dm2min = 0.01;			//eV**2
Double_t dm2max = 100.;                 //eV**2
Double_t sin22thmin = 0.0001;
Double_t sin22thmax = 1.0;

int update_cov(bool shape_only, int nbinsE, int nL, int sigDm, int sigSin22th);
double chisq_calc(bool shape_only, int nL, int nbinsE, double signalSin22th, double signalDm2, bool smear, bool drawScale, bool drawSmear);
double getDMpt(int dm);
double getSin22THpt(int sint);

int multiple_detector_fit(){

  // -----------Declarations and modifiers--------------------

  std::string mode = "nu";	//beam mode to run in
  bool use100m = true;		//Include the detector at 100m?
  bool use470m = true;		//Include the detector at 470m?
  bool use700m = true;		//Include the detector at 700m?

  bool shape_only = true;

  double signal[2] = {.03,1.2};	//Injected signal, in form (sin22t,dm2)
  bool printJelly = false;
  std::string jellyTitle = "pngs/jellybean_dm2_10_sin22th_01.png";

  // Chisq Surface(testing)
  drawSurf = true;
  bool printSurf = false;
  std::string surfTitle = "chisqSurf_test.pdf";

  // Scale Factor Histogram (will only work with far detector and only really significant with shape-only)
  bool drawScale = false;
  bool printScale = false;
  std::string scaleTitle = "pngs/evd_scaleDiff.png";
  predSin22th = .05;
  predDm2 = 50.;

  // Smearing Histogram
  bool drawSmear = true;



  // ---------------------------------------------------------

  std::vector<std::string> baselines;
  if (use100m) baselines.push_back("100m");
  if (use470m) baselines.push_back("470m");
  if (use700m) baselines.push_back("600m_onaxis");
  Int_t nL = baselines.size();
  std::cout << " Number of baselines " << nL << std::endl;
  if (nL == 0){
    std::cout << "\nNo baselines selected!  What did you expect? Exiting.\n" << std::endl;
    return 2;
  }

  NULLVec.resize(nL);
  OscVec.resize(nL);

  for( int l = 0; l < nL; l++){
    OscVec[l].resize(npoints+1);
  }

  int counter = 0;
  int nbinsE = 0;

  std::cout << "Beginning : ... " << std::endl;

  for(int i = 0; i < nL; i++){
    std::string temp_name = "../MatrixFiles/combined_ntuple_" + baselines[i] + "_" + mode + "_processed_numu_Joseph_Smeared.root";
    TFile temp_file(temp_name.c_str());

    std::string temp_name_syst = "../MatrixFiles/combined_ntuple_" + baselines[i] + "_" + mode + "_processed_numu.root";
    TFile temp_file_syst(temp_name_syst.c_str());
    TH1D *NULL_BASE;

    NULL_BASE = (TH1D*)(temp_file.Get("NumuCC"));
    TH1D *Nom_base;

    Nom_base = (TH1D*)(temp_file_syst.Get("NumuCC"));
    nbinsE = NULL_BASE->GetNbinsX();

    for(int j = 1; j <= nbinsE; j++){
      NULLVec[i].push_back(1*(NULL_BASE->GetBinContent(j)));
      NomVec[i][j-1] = 1*(Nom_base->GetBinContent(j));
    }

    for(int dm = 0; dm <= npoints; dm++){
      TH1D *OSC_BASE;
      std::string dmpoint = std::to_string(dm);
      std::string name = "Osc_";
      name = name+dmpoint;
      OSC_BASE = (TH1D*)(temp_file.Get(name.c_str()));
      OSC_BASE->Rebin(1);
      for(int j = 1; j <= nbinsE; j++){
        OscVec[i][dm].push_back(1*(OSC_BASE->GetBinContent(j)));
      }
      delete OSC_BASE;
    }

    for(int u = 0; u < spoints; u++){
      for(int s = 0; s < 7; s++){
        TH1D *SYST_BASE;
    	TString upoint = Form("%d",u);
    	TString name = "Universe_";
    	TString name2 = "_MultiSim_";
    	TString mul = Form("%d",s);

    	name += upoint;
    	name += name2;
    	name += mul;

        SYST_BASE = (TH1D*)(temp_file_syst.Get(name));
        SYST_BASE->Rebin(1);
	for(int k = 1; k <= nbinsE; k++){
	  SystVec[i][u][s][k-1] = 1*(SYST_BASE->GetBinContent(k));
	}
	delete SYST_BASE;
      } 
    }
    delete NULL_BASE;
    temp_file.Close();
    counter++;
  }
  if(counter == 0){ std::cout << " OOPS! Didn't fill the vectors" << std::endl; return 2;}


  //Stuff for root, drawing options and such.
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(0000);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);

  TH1::AddDirectory(kFALSE);
  //This line is important.  But so stupid. 
  //If not included, when histos are read in they attach to the file they are from.
  //When the file goes out of scope, the histo dies and your code will seg fault.

// Build the covariance matrix
  update_cov(shape_only, nbinsE, nL, -1, -1);

// Calculate chisquared
  chisq_calc(shape_only,nL,nbinsE,signal[0],signal[1], true, drawScale, drawSmear);

// Update covariance matrix, and then do chisquared one more time with updated best fit points
  update_cov(shape_only, nbinsE, nL, bestFit->GetY(), bestFit->GetX());
  chisq_calc(shape_only,nL,nbinsE,signal[0],signal[1], false, false, false);

// Begin Drawing
  fiveSigma->SetMarkerSize(5);
  fiveSigma->SetMarkerColor(kBlue);
  fiveSigma->SetFillColor(kBlue);
  fiveSigma->SetMarkerStyle(7);
  threeSigma->SetMarkerSize(5);
  threeSigma->SetMarkerStyle(7);
  threeSigma->SetMarkerColor(kGreen + 1);
  threeSigma->SetFillColor(kGreen + 1);
  ninety->SetMarkerColor(kYellow);
  ninety->SetFillColor(kYellow);
  ninety->SetMarkerStyle(7);

  printf("\nDrawing best fits...\n");
	
  TCanvas* c3 = new TCanvas("c3","Sensitivity",700,700);
  c3->SetLeftMargin(.15);
  c3->SetBottomMargin(.15);
  c3->SetTopMargin(.05);
  c3->SetRightMargin(.05);
  c3->SetLogx();
  c3->SetLogy();

  TH2D* hr1=new TH2D("hr1","hr1",500,sin22thmin*10,sin22thmax,500,dm2min+0.0001,dm2max);
  hr1->Reset();
  hr1->SetFillColor(0);
  hr1->SetTitle(";sin^{2}2#theta#lower[0.4]{#mu#kern[-0.3]{#mu}};#Deltam_{41}^{2} [eV^{2}]");
  hr1->GetXaxis()->SetTitleOffset(1.2);
  hr1->GetYaxis()->SetTitleOffset(1.2);
  hr1->GetXaxis()->SetTitleFont(62);
  hr1->GetYaxis()->SetTitleFont(62);
  hr1->GetYaxis()->CenterTitle();
  hr1->GetXaxis()->CenterTitle();
  hr1->GetXaxis()->SetTitleSize(0.05);
  hr1->GetXaxis()->SetLabelSize(0.04);
  hr1->GetXaxis()->SetLabelOffset(0.001);
  hr1->GetYaxis()->SetTitleSize(0.05);
  hr1->GetYaxis()->SetLabelSize(0.04);
  hr1->SetStats(kFALSE);
  hr1->Draw();

  fiveSigma->Draw("p");
  threeSigma->Draw("psame");
  ninety->Draw("psame");
  bestFit->Draw("psame");
  signalPt->Draw("psame");
  c3->RedrawAxis();

  TLatex *tex_Detector = new TLatex(.2,.23,"#splitline{LAr1-ND (100m)}{and T600 (600m, on axis)}");
//  TLatex *tex_Detector = new TLatex(.2,.23,"MicroBooNE (470m)");
  tex_Detector->SetNDC();
  tex_Detector->SetTextFont(62);
  tex_Detector->SetTextSize(0.035);
  tex_Detector->Draw(); 
 
  TLatex *tex_pre = new TLatex(.18,.75,"PRELIMINARY");
  tex_pre->SetNDC();
  tex_pre->SetTextFont(62);
  tex_pre->SetTextColor(kRed-3);
  tex_pre->SetTextSize(0.03);
  tex_pre->Draw();

  //char str_signal [100];
  //sprintf (str_signal, "Injected Signal - #Delta m^{2}: %.1f, sin^{2}(2#theta): %.2f", signalPt->GetY(),signalPt->GetX());
  std::string str_signal;
  if(shape_only) str_signal = "Shape Only";
  else str_signal = "Shape + Rate";
  TLatex *tex_Signal = new TLatex(.18,.92,str_signal.c_str());
  tex_Signal->SetNDC();
  tex_Signal->SetTextFont(62);
  tex_Signal->SetTextSize(0.025);
  tex_Signal->Draw();

  TLatex *tex_mode = new TLatex(.18,.89,"#nu mode, CC Events");
  tex_mode->SetNDC();
  tex_mode->SetTextFont(62);
  tex_mode->SetTextSize(0.025);
  tex_mode->Draw();

  TLatex *tex_un = new TLatex(.18,.86,"Statistical Uncertainty Only");
  tex_un->SetNDC();
  tex_un->SetTextFont(62);
  tex_un->SetTextSize(0.025);
  tex_un->Draw();

  TLatex *tex_E = new TLatex(.18,.83,"Reconstructed Energy");
  tex_E->SetNDC();
  tex_E->SetTextFont(62);
  tex_E->SetTextSize(0.025);
  tex_E->Draw();

  TLatex *tex_eff = new TLatex(.18,.80,"80% #nu#lower[0.4]{#mu} Efficiency");
  tex_eff->SetNDC();
  tex_eff->SetTextFont(62);
  tex_eff->SetTextSize(0.025);
  tex_eff->Draw();

  TLegend* legt=new TLegend(0.2,0.3,0.6,0.45);
  legt->SetFillStyle(0);
  legt->SetFillColor(0);
  legt->SetBorderSize(0);
  legt->SetTextFont(62);
  legt->SetTextSize(0.03);
  legt->AddEntry(fiveSigma,"5 #sigma Confidence","f");
  legt->AddEntry(threeSigma,"3 #sigma Confidence","f");
  legt->AddEntry(ninety,"90% Confidence","f");
  legt->AddEntry(signalPt,"Signal","p");
  legt->Draw();

  if(printJelly){
    c3->Print(jellyTitle.c_str());
  }

// Print Chisq Surface!
  if(drawSurf){
    TCanvas* c4 = new TCanvas("c4","Sensitivity",700,700);
    c4->SetLeftMargin(.15);
    c4->SetBottomMargin(.15);
    c4->SetTopMargin(.05);
    c4->SetRightMargin(.05);
    c4->SetLogx();
    c4->SetLogy();

    chisqSurf->SetTitle("");
    chisqSurf->GetXaxis()->SetTitle("sin^{2}2#theta#lower[0.4]{#mu#kern[-0.3]{#mu}}");
    chisqSurf->GetXaxis()->CenterTitle(true);
    chisqSurf->GetXaxis()->SetLabelFont(62);
    chisqSurf->GetXaxis()->SetLabelOffset(0.003);
    chisqSurf->GetXaxis()->SetLabelSize(0.03);
    chisqSurf->GetXaxis()->SetTitleSize(0.05);
    chisqSurf->GetXaxis()->SetTitleFont(62);
    chisqSurf->GetXaxis()->SetTitleOffset(1.75);

    chisqSurf->GetYaxis()->SetTitle("#Deltam_{41}^{2} [eV^{2}]");
    chisqSurf->GetYaxis()->CenterTitle(true);
    chisqSurf->GetYaxis()->SetLabelFont(62);
    chisqSurf->GetYaxis()->SetLabelSize(0.03);
    chisqSurf->GetYaxis()->SetTitleSize(0.05);
    chisqSurf->GetYaxis()->SetTitleFont(62);
    chisqSurf->GetYaxis()->SetTitleOffset(2.0);

    chisqSurf->GetZaxis()->SetTitle("Confidence Level [#sigma]");
    chisqSurf->GetZaxis()->CenterTitle(true);
    chisqSurf->GetZaxis()->SetLabelFont(62);
    chisqSurf->GetZaxis()->SetLabelSize(0.03);
    chisqSurf->GetZaxis()->SetTitleSize(0.05);
    chisqSurf->GetZaxis()->SetTitleFont(62);
    chisqSurf->GetZaxis()->SetNdivisions(506);
    chisqSurf->GetZaxis()->SetTitleOffset(1.25);
    chisqSurf->GetZaxis()->SetRangeUser(0.0000001, 10.5);

    chisqSurf->SetMarkerColor(kBlue);
    chisqSurf->Draw("P");

    if(printSurf){
      c4->Print(surfTitle.c_str());
    }
  }

 

  // Draw Scale Plot
  if(drawScale){
    TCanvas* c5 = new TCanvas("c5","Scale Plot",700,700);
    c5->SetLeftMargin(.15);
    c5->SetBottomMargin(.15);
    c5->SetTopMargin(.05);
    c5->SetRightMargin(.05);
    c5->cd();

    cvHist->SetFillColor(0);
    predHist->SetTitle(";sin^{2}2#theta#lower[0.4]{#mu#kern[-0.3]{#mu}};#Deltam_{41}^{2} [eV^{2}]");
    predHist->GetXaxis()->SetTitleOffset(1.2);
    predHist->GetYaxis()->SetTitleOffset(1.2);
    predHist->GetXaxis()->SetTitleFont(62);
    predHist->GetYaxis()->SetTitleFont(62);
    predHist->GetYaxis()->CenterTitle();
    predHist->GetXaxis()->CenterTitle();
    predHist->GetXaxis()->SetTitleSize(0.05);
    predHist->GetXaxis()->SetLabelSize(0.04);
    predHist->GetXaxis()->SetLabelOffset(0.001);
    predHist->GetYaxis()->SetTitleSize(0.05);
    predHist->GetYaxis()->SetLabelSize(0.04);
    predHist->SetStats(kFALSE);

    TLegend* legS=new TLegend(0.6,0.3,0.85,0.45);
    legS->SetFillStyle(0);
    legS->SetFillColor(0);
    legS->SetBorderSize(0);
    legS->SetTextFont(62);
    legS->SetTextSize(0.03);
    legS->AddEntry(predHist,"Prediction Vector","h");
    legS->AddEntry(cvHist,"Signal Vector","h");
    legS->AddEntry(predScaledHist,"Scaled Prediction","p");

    cvHist->SetLineColor(kBlue-3);
    cvHist->SetFillColor(kBlue-10);
    predHist->SetLineColor(kRed-3);
    predHist->SetFillColor(kRed-10);
    predHist->SetMarkerStyle(20);
    predScaledHist->SetMarkerStyle(3);

    predHist->Draw("h");
    cvHist->Draw("hsame");
    predScaledHist->Draw("phsame");

    legS->Draw();

    char str_pred [100];
    sprintf (str_pred, "Prediction  - #Delta m^{2}: %.1f, sin^{2}(2#theta): %.2f", predDm2, predSin22th);
    TLatex *tex_pred = new TLatex(.54,.92,str_pred);
    tex_pred->SetNDC();
    tex_pred->SetTextFont(62);
    tex_pred->SetTextSize(0.025);
    tex_pred->Draw();
    tex_pre->SetX(.725575);
    tex_pre->Draw();
    tex_mode->SetX(.727011);
    tex_mode->Draw();
    tex_un->SetX(.640805);
    tex_un->Draw();
    tex_E->SetX(.691092);
    tex_E->Draw();
    tex_eff->SetX(.734195);
    tex_eff->Draw();

    TLatex *tex_Det = new TLatex(.2,.23,"T600 (600m)");
    tex_Det->SetNDC();
    tex_Det->SetTextFont(62);
    tex_Det->SetTextSize(0.035);
    tex_Det->Draw(); 

    if(printScale){
      c5->Print(scaleTitle.c_str());
    }
  }

  // Draw Smear Plot
  if(drawSmear){
    TCanvas* c6 = new TCanvas("c6","Smear Plot",700,700);
    c6->SetLeftMargin(.15);
    c6->SetBottomMargin(.15);
    c6->SetTopMargin(.05);
    c6->SetRightMargin(.05);
    c6->cd();

    /*  smeared->SetFillColor(0);
    unsmeared->SetTitle(";sin^{2}2#theta#lower[0.4]{#mu#kern[-0.3]{#mu}};#Deltam_{41}^{2} [eV^{2}]");
    unsmeared->GetXaxis()->SetTitleOffset(1.2);
    unsmeared->GetYaxis()->SetTitleOffset(1.2);
    unsmeared->GetXaxis()->SetTitleFont(62);
    unsmeared->GetYaxis()->SetTitleFont(62);
    unsmeared->GetYaxis()->CenterTitle();
    unsmeared->GetXaxis()->CenterTitle();
    unsmeared->GetXaxis()->SetTitleSize(0.05);
    unsmeared->GetXaxis()->SetLabelSize(0.04);
    unsmeared->GetXaxis()->SetLabelOffset(0.001);
    unsmeared->GetYaxis()->SetTitleSize(0.05);
    unsmeared->GetYaxis()->SetLabelSize(0.04);
    unsmeared->SetStats(kFALSE);*/

    TLegend* legS=new TLegend(0.6,0.3,0.85,0.45);
    legS->SetFillStyle(0);
    legS->SetFillColor(0);
    legS->SetBorderSize(0);
    legS->SetTextFont(62);
    legS->SetTextSize(0.03);
    legS->AddEntry(unsmeared,"Unsmeared","h");
    legS->AddEntry(smeared,"Statistically Smeared","h");

    smeared->SetLineColor(0);
    smeared->SetFillColor(kBlue-10);

    unsmeared->SetLineWidth(2);
    unsmeared->SetLineColor(kBlack);
    unsmeared->SetMarkerColor(kBlack);
    unsmeared->SetMarkerSize(0.1);
    unsmeared->SetMarkerStyle(2);

    smeared->Draw("h");
    unsmeared->Draw("ep same");

    legS->Draw();

  }

  return 0;
}

int update_cov(bool shape_only, int nbinsE, int nL, int sigDm, int sigSin22th){

//Build Matrices and Vectors
  TMatrixT <float> M(nbinsE*nL,nbinsE*nL);
  TMatrixT <float> M1(nbinsE*nL,nbinsE*nL);
  TMatrix N(nbinsE*nL,1);
  double NT = 0;

  M.Zero();
  M1.Zero();
  N.Zero();

  std::vector< float > CVvec;

// Fill Error Matrix
  int Erri = 0, Errj = 0;
  int n = 0;

  std::cout << "Filling Error Matrix..." << std::endl;

  //We are now filling the Error (Covariance) Matrix
  for(int Lrow = 0; Lrow < nL; Lrow++){
    for(int Erow = 0; Erow < nbinsE; Erow++){
      Errj = 0;
      for(int Lcol = 0; Lcol < nL; Lcol++){
        for(int Ecol = 0; Ecol < nbinsE; Ecol++){
	  for(int u = 0; u < spoints; u++){

	    //This is the Error Matrix Equation (See Dave's Thesis)
	    M (Erri,Errj) += (NomVec[Lrow][Erow]-SystVec[Lrow][u][6][Erow])*(NomVec[Lcol][Ecol]-SystVec[Lcol][u][6][Ecol]);
	    n++;
	  }
	  M (Erri,Errj) /= n;			// Average over universes
	  if(!shape_only){

	    // Build "Fractional" Error Matrix
	    M (Erri,Errj) /= NomVec[Lrow][Erow]*NomVec[Lcol][Ecol]; 
	    // Now  take the real statistics to fill the Error Matrix. Here is where we'll update with the best fit point.
	    if(sigDm == -1)  M(Erri, Errj) *= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];		
	    else M(Erri, Errj) *= (NULLVec[Lrow][Erow] - (OscVec[Lrow][sigDm][Erow] * getSin22THpt(sigSin22th)))*(NULLVec[Lcol][Ecol] - (OscVec[Lcol][sigDm][Ecol] * getSin22THpt(sigSin22th)));
	    // Add Stat. Errors
	    if(Erri == Errj){ 
	      M (Erri, Errj) += NULLVec[Lrow][Erow];
	    }
	  }
	  n = 0;
	  M1 (Erri, Errj) = M (Erri, Errj);
	  Errj++;
	}
      }

      // Fill vectors for chi^2 calaculations
      N(Erri,0) = (NULLVec[Lrow][Erow]);
      NT += (NULLVec[Lrow][Erow]);
      Erri++;
    }
  }

  // If we want to do a shape only analysis we need to factorize the Error Matrix
  // such that we remove any uncertainty associated with only the normalization (we will leave the mixed components)
  if(shape_only){

    TMatrix Mik(nbinsE*nL,1);
    TMatrix Mkj(1,nbinsE*nL);
    double Mkl = 0;
    Mik.Zero();
    Mkj.Zero();
    for(int i = 0; i < nbinsE*nL; i++){
      for(int j = 0; j < nbinsE*nL; j++){

	 Mik(i,0) += M(i,j);
	 Mkj(0,i) += M(j,i);
	 Mkl += M(i,j);
      }
    }

    //Subtract the Normalization component
    for(int i = 0; i < nbinsE*nL; i++){
      for(int j = 0; j < nbinsE*nL; j++){

	M(i, j) -= (N(i,0))*(N(j,0))*(Mkl)/(NT*NT);
      }
    }

    Errj = Erri = 0;
    //Now do what we did before to finalize the error matrix
    for(int Lrow = 0; Lrow < nL; Lrow++){
      for(int Erow = 0; Erow < nbinsE; Erow++){
	Errj = 0;
	for(int Lcol = 0; Lcol < nL; Lcol++){
	  for(int Ecol = 0; Ecol < nbinsE; Ecol++){

	    M (Erri,Errj) /= NomVec[Lrow][Erow]*NomVec[Lcol][Ecol];
	    if(sigDm == -1)  M(Erri, Errj) *= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];		
	    else M(Erri, Errj) *= (NULLVec[Lrow][Erow] - (OscVec[Lrow][sigDm][Erow] * getSin22THpt(sigSin22th)))*(NULLVec[Lcol][Ecol] - (OscVec[Lcol][sigDm][Ecol] * getSin22THpt(sigSin22th)));
	    if(Erri == Errj){ M (Erri, Errj) += NULLVec[Lrow][Erow];}
            Errj++;		
	  }
        }
        Erri++;
      }
    }
  }

  std::cout << "...Error Matrix Filled" << std::endl;

//Now create the inverse covariance matrix.
  std::cout << "Invert in...3...2...1..." << std::endl;
  cov.ResizeTo(nbinsE*nL,nbinsE*nL);
  cov = M.Invert();
  std::cout << "Inverted!" << std::endl;

  return 0;
}

double chisq_calc(bool shape_only, int nL, int nbinsE, double signalSin22th, double signalDm2, bool smear, bool drawScale, bool drawSmear){

  if(drawScale){
    double bins[20] = {.200, .300, .400, .450, .500, .550, .600, .650, .700, .750, .800, .850, .900, .950, 1.000, 1.250, 1.500, 2.000, 2.500, 3.000};
    cvHist = new TH1D("cvHist", "", nbinsE, bins);
    predHist = new TH1D("predHist","",nbinsE,bins);
    predScaledHist = new TH1D("predScaleHist","",nbinsE,bins);
  }
  if(drawSmear){
    double bins[20] = {.200, .300, .400, .450, .500, .550, .600, .650, .700, .750, .800, .850, .900, .950, 1.000, 1.250, 1.500, 2.000, 2.500, 3.000};
    smeared = new TH1D("smeared", "", nbinsE, bins);
    unsmeared = new TH1D("unsmeared","",nbinsE,bins);
  }

  int scaleDm2pt = 0, scaleSin22thpt = 0;

  // Fill Prediction Matrix
  std::vector< std::vector< TMatrixT <float> > >  Pred;
  Pred.resize(npoints+1);
  for(int i = 0; i < npoints+1; i++){
    Pred[i].resize(npoints+1);
    for(int j = 0; j < npoints+1; j++){
      Pred[i][j].ResizeTo(nbinsE*nL,1);
    }
  }

  double sin22th;
  int Predi = 0;
  std::cout << "Filling Prediction Matrices..." << std::endl;

  int mydm2pt = 0;
  int mysinpt = 0;
  double chiLow = -1;
  TMatrixT <float> CV(nbinsE*nL,1);
  CV.Zero();

  for(int dm = 0; dm <= npoints; dm++){
    // Here, we also try to find the closest estimate for our dm2 point within the set of npoints points
    if(getDMpt(dm) <= signalDm2)  mydm2pt = dm;
    if(drawScale){
      if(getDMpt(dm) <= predDm2) scaleDm2pt = dm;
    }
    for(int s = 0; s <= npoints; s++){
      Predi = 0; 
      for(int Lbin = 0; Lbin < nL; Lbin++){
        for(int Ebin = 0; Ebin < nbinsE; Ebin++){
	  sin22th = getSin22THpt(s);
	  if(sin22th <= signalSin22th && s > mysinpt)  mysinpt = s;
	  if(drawScale && sin22th <= predSin22th && s > scaleSin22thpt) scaleSin22thpt = s;
  	  Pred[dm][s] (Predi,0) = NULLVec[Lbin][Ebin] - (OscVec[Lbin][dm][Ebin])*(sin22th);
	  Predi++;
	}
      }
    }
  }

  signalDm2 = getDMpt(mydm2pt);				// Correct signal point for resolution
  signalSin22th = getSin22THpt(mysinpt);

  TMatrixT <float> CV_noSmear(nbinsE*nL,1);
  for(int Lbin = 0; Lbin < nL; Lbin ++){
    for(int Ebin = 0; Ebin < nbinsE; Ebin ++){
      CV_noSmear(Lbin*nbinsE + Ebin,0) = (NULLVec[Lbin][Ebin] - (OscVec[Lbin][mydm2pt][Ebin])*signalSin22th);
      if(smear){
	if(Lbin == 0 && Ebin == 0)std::cout << "... Smearing based on Statistics ..." << std::endl;
	TRandom3 * r3 = new TRandom3();
	CV(Lbin*nbinsE + Ebin,0) = r3->Gaus(CV_noSmear(Lbin*nbinsE + Ebin, 0),sqrt(CV_noSmear(Lbin*nbinsE + Ebin,0)));
	if(Lbin == 1){
	  smeared->SetBinContent(Ebin+1,CV(Lbin*nbinsE + Ebin,0));
	  unsmeared->SetBinContent(Ebin+1,CV_noSmear(Lbin*nbinsE + Ebin, 0));

	}
      }
      else{
	CV(Lbin*nbinsE + Ebin,0) = CV_noSmear(Lbin*nbinsE + Ebin, 0);
      }
    }
  }

  

  signalPt = new TMarker(signalSin22th,signalDm2,7);
  signalPt->SetMarkerColor(kRed);

  std::cout << "...Prediction and Signal Matrices Filled" << std::endl;

// Calculate Chi2
  fiveSigma = new TGraph();
  threeSigma = new TGraph();
  ninety = new TGraph();
  chisqSurf = new TGraph2D();
  Double_t chisq;

  int fiveSCounter = 0;
  int threeSCounter = 0;
  int surfCounter = 0;
  int ninetyCounter = 0;
  int ncounter = 0;
  double bestFit_sin22th, bestFit_dm2;

  for(int dm2 = 0; dm2 <= npoints; dm2++){
    for(int sint = 0; sint <= npoints; sint++){
      sin22th = getSin22THpt(sint);
      chisq = 0.0;
      double Scaling[nbinsE*nL];

      if(shape_only){
        //Scale all measurements to the ND NULL levels             
        for(int Edecbin = 0; Edecbin < nbinsE; Edecbin++){
	  Scaling[Edecbin]  = CV (Edecbin,0);
	  if(Pred[dm2][sint] (Edecbin,0) != 0){Scaling[Edecbin] /= Pred[dm2][sint] (Edecbin,0);}
	  if(Pred[dm2][sint] (Edecbin,0) == 0){ std::cout << "Scaling failed....Pred==0!" << std::endl; break;}
	}
        // Copy results from first baseline to the rest of these baselines
        int tick = 0;
        for(int Edecbin = nbinsE; Edecbin < nbinsE*nL; Edecbin++){
	  if(tick == nbinsE) tick = 0;
	  Scaling[Edecbin] = Scaling[tick];
	  tick++;
        }
      } //Shape only

      if(!shape_only){
        //We don't want to scale anything!
        for(int Edecbin = 0; Edecbin < nbinsE*nL; Edecbin++){
	  Scaling[Edecbin] = 1;
        }
      } // Shape and rate

      TMatrixT <float> PredScaled(nbinsE*nL,1);

      for(int i = 0; i < nbinsE*nL; i++){
        PredScaled(i,0) = Scaling[i]*(Pred[dm2][sint] (i,0));
      }

      if(drawScale && scaleDm2pt == dm2 && scaleSin22thpt == sint){
        for(int eBins = 0; eBins < nbinsE; eBins ++){
          cvHist->SetBinContent(eBins,(CV (((nL-1)*nbinsE) + eBins, 0)/(cvHist->GetXaxis()->GetBinWidth(eBins+1))));
          predHist->SetBinContent(eBins,(Pred[scaleDm2pt][scaleSin22thpt] (((nL-1)*nbinsE) + eBins, 0))/(cvHist->GetXaxis()->GetBinWidth(eBins+1)));
          predScaledHist->SetBinContent(eBins,(PredScaled (((nL-1)*nbinsE) + eBins,0))/(cvHist->GetXaxis()->GetBinWidth(eBins+1)));
        }
      }

      TMatrixT <float> Final(nbinsE*nL,1);
      TMatrixT <float> FinalT(1,nbinsE*nL);

      for(int Edecbin = 0; Edecbin < nbinsE*nL; Edecbin++){
        Final(Edecbin,0) = (CV(Edecbin,0)) - ((PredScaled (Edecbin,0))); 
        FinalT(0,Edecbin) = (CV(Edecbin,0)) - ((PredScaled (Edecbin,0)));
      } //Sum over all detectors and energies to determine the chi2 for that osc. point	      

      TMatrixT <float> middle(1, nbinsE*nL);
      middle.Mult(FinalT,cov);
      TMatrixT <float> Fin(1, 1);
      Fin.Mult(middle,Final);
      chisq = Fin(0,0);

  // Draw Jellybeans
      if(chisq <= 1.64){
	ninety->SetPoint(ninetyCounter,getSin22THpt(sint),getDMpt(dm2));
	ninetyCounter++;
      }
      if(chisq <= 9){
        threeSigma->SetPoint(threeSCounter,getSin22THpt(sint),getDMpt(dm2));
        threeSCounter++;
      }
      if(chisq <= 25) {
	  fiveSigma->SetPoint(fiveSCounter,getSin22THpt(sint),getDMpt(dm2));
	  fiveSCounter ++;
      }
     
  // Draw Surface
      if(drawSurf){
       if(chisq <= 100.) chisqSurf->SetPoint(surfCounter,getSin22THpt(sint),getDMpt(dm2),sqrt(chisq));
       else chisqSurf->SetPoint(surfCounter,getSin22THpt(sint),getDMpt(dm2),10.);
       surfCounter ++;
      }

  // Now, try to find best fit
      if(chiLow < 0) chiLow = chisq;
      if(chisq < chiLow){
        chiLow = chisq;
        bestFit_sin22th = getSin22THpt(sint);
        bestFit_dm2 = getDMpt(dm2);
      }
      ncounter ++;
    }
  }
  bestFit = new TMarker(bestFit_sin22th,bestFit_dm2, 3);
  return bestFit_sin22th;
}

double getDMpt(int dm){
  return pow(10.,(TMath::Log10(dm2min)+ (dm * (1./npoints))*TMath::Log10(dm2max/dm2min)));
}

double getSin22THpt(int sint){
  return pow(10., (TMath::Log10(sin22thmin)+ (sint * (1./npoints))*TMath::Log10(sin22thmax/sin22thmin)));
}

#ifndef __CINT__
int main()
{
  multiple_detector_fit();
  return 0;
}
# endif

void jellybean(){
  multiple_detector_fit();
  return;
}

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

std::vector< std::vector <float> >   NULLVec;						// Null hypothesis
std::vector< std::vector < std::vector <float> > >   OscVec;		// Oscillated hypothesis for each dm2 point
std::vector< std::vector< TMatrixT <float> > >  Pred;

double NomVec[3][20];
double SystVec[3][1001][7][20];
double chisqArray[251001];

TMatrixT <float> cov;
TMatrixT <float> CV;

// Here is our list of graphs and markers
TGraph2D * chisqSurf;
TGraph * fiveSigma;
TGraph * threeSigma;
TGraph * ninety;
TGraph * bestFits;
TMarker * signalPt;
TMarker * bestFit;

// Here is our list of histos
TH1D * cvHist;
TH1D * predHist;
TH1D * predScaledHist;
TH1D *smeared;
TH1D *unsmeared;
TH1D *bestfit;


bool use100m, use470m, use700m;
bool drawSurf, drawSmear, drawScale, smearBestFit, drawBestFitPlot;
bool shape_only;

double predSig[2];
int signalReso[2];
int bestFitReso[2];
int predSigReso[2];
double POT_weight[3];

//plot properties
Int_t npoints = 500;			// Grid Pointsl
Int_t spoints = 1000;			// Universes

//grid boundaries
Double_t dm2min = 0.01;			//eV**2
Double_t dm2max = 100.;         	//eV**2
Double_t sin22thmin = 0.0001;
Double_t sin22thmax = 1.0;

int update_cov(int nbinsE, int nL, bool useBestFit);
int build_vectors(int nL, int nbinsE, double signalSin22th, double signalDm2, bool smear);
double chisq_calc(int nL, int nbinsE);
int draw_jellybeans(double chiLow, bool exclude);
double getDMpt(int dm, bool inverse);
double getSin22THpt(int sint, bool inverse);

int multiple_detector_fit(){

  // -----------Declarations and modifiers--------------------

  std::string mode = "nu";	//beam mode to run in
  use100m = true;		//Include the detector at 100m?
  use470m = true;		//Include the detector at 470m?
  use700m = true;		//Include the detector at 700m?

  shape_only = true;
  bool smear = true;
  bool iterate = true;
  bool exclude = true;

  double signal[2] = {.015,2.5};	//Injected signal, in form (sin22t,dm2)

  // Chisq Surface(testing)
  drawSurf = false;

  // Scale Factor Histogram (will only work with far detector and only really significant with shape-only)
  drawScale = false;
  predSig[0] = .03;		// Prediction Vector sin22th
  predSig[1] = 1.2;		// prediction Vector dm2

  // Smearing Histogram
  drawSmear = false;
  smearBestFit = false;

  // Plot of Best Fit Points
  drawBestFitPlot = true;
  int bestFitIterations = 500;

  // ---------------------------------------------------------

  if(!smear) iterate = false;

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

  if( use100m &&  use470m && !use700m){POT_weight[0] = 0.23895; POT_weight[1] = 1.23895;}
  if( use100m && !use470m &&  use700m){POT_weight[0] = 1; 	POT_weight[1] = 1;}
  if( use100m &&  use470m &&  use700m){POT_weight[0] = 1;       POT_weight[1] = 2;       POT_weight[2] = 1;}
  if(!use100m &&  use470m && !use700m){POT_weight[0] = 1;}

  NULLVec.resize(nL);
  OscVec.resize(nL);
  for( int l = 0; l < nL; l++){
    OscVec[l].resize(npoints+1);
  }

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
      NULLVec[i].push_back(POT_weight[i]*(NULL_BASE->GetBinContent(j)));
      NomVec[i][j-1] = POT_weight[i]*(Nom_base->GetBinContent(j));
    }

    for(int dm = 0; dm <= npoints; dm++){
      TH1D *OSC_BASE;
      std::string dmpoint = std::to_string(dm);
      std::string name = "Osc_";
      name = name+dmpoint;
      OSC_BASE = (TH1D*)(temp_file.Get(name.c_str()));
      OSC_BASE->Rebin(1);
      for(int j = 1; j <= nbinsE; j++){
        OscVec[i][dm].push_back(POT_weight[i]*(OSC_BASE->GetBinContent(j)));
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
	  SystVec[i][u][s][k-1] = POT_weight[i]*(SYST_BASE->GetBinContent(k));
	}
	delete SYST_BASE;
      }
    }
    delete NULL_BASE;
    temp_file.Close();
  }


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


  if(!drawBestFitPlot){
    double chi_current, chiLow;
    // Build the prediction and signal vectors
    build_vectors(nL, nbinsE, signal[0], signal[1], smear);
    // Build the covariance matrix
    update_cov(nbinsE, nL, false);
    // Calculate chisquared for the first time.
    chiLow = chisq_calc(nL,nbinsE);
    std::cout << "... chisq = " << chiLow << "..." << std::endl;
    // Iteration
    if(iterate){
      for(int iter = 0; iter < 10; iter++){
        update_cov(nbinsE, nL, true);
        chi_current = chisq_calc(nL,nbinsE);
        std::cout << "... chisq = " << chi_current << "..." << std::endl;
        if(abs(chi_current - chiLow) < .002){
          chiLow = chi_current;
          break;
        }
        else if(iter == 9)	std::cout << "No convergence of best fit in 10 tries. Sorry." << std::endl;
        chiLow = chi_current;
      }
    }
    // Then, draw the jellybeans!
    draw_jellybeans(chiLow, exclude);
  }
  else{
    bestFits = new TGraph();
    for(int i = 0; i < bestFitIterations; i ++){
      build_vectors(nL, nbinsE, signal[0],signal[1],smear);
      update_cov(nbinsE,nL,false);
      chisq_calc(nL,nbinsE);
      update_cov(nbinsE, nL, true);
      chisq_calc(nL,nbinsE);
      bestFits->SetPoint(i,bestFit->GetX(),bestFit->GetY());
      std::cout << "............................." << (i+1) << "/" << bestFitIterations << "............................." << std::endl;
    }
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
    legS->AddEntry(predHist,"Prediction Vector","f");
    legS->AddEntry(cvHist,"Signal Vector","f");
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
  }

// Draw Smear Plot
  if(drawSmear){
    TCanvas* c6 = new TCanvas("c6","Smear Plot",700,700);
    c6->SetLeftMargin(.15);
    c6->SetBottomMargin(.15);
    c6->SetTopMargin(.05);
    c6->SetRightMargin(.05);
    c6->cd();

    unsmeared->SetTitle(";Reconstructed Energy [GeV];Events / Bin");
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
    unsmeared->SetStats(kFALSE);
    TGaxis::SetMaxDigits(3);

    TLegend* legS=new TLegend(0.475,0.82,0.8,0.94);
    legS->SetFillStyle(0);
    legS->SetFillColor(0);
    legS->SetBorderSize(0);
    legS->SetTextFont(62);
    legS->SetTextSize(0.02);
    legS->AddEntry(unsmeared,TString::Format("Unsmeared Signal : (%4.2f eV^{2}, %0.3f)",signal[1],signal[0]),"f");
    legS->AddEntry(smeared,TString::Format("Smeared Signal : (%4.2f eV^{2}, %0.3f)",signal[1],signal[0]),"fep");//"Smeared Signal","fep");
    if(smearBestFit)  legS->AddEntry(bestfit,TString::Format("#chi#lower[-0.3]{2} Best Fit : (%4.2f eV^{2}, %0.3f)",getDMpt(bestFitReso[1],false),getSin22THpt(bestFitReso[0],false)),"fep");
    //TString::Format("100*((%f)*pow(TMath::Sin((1.267)*(%f)*(.001*x/%f)),2))",sin2,m2,energy)

    unsmeared->SetLineColor(0);
    unsmeared->SetFillColor(kBlue-10);

    smeared->SetLineWidth(2);
    smeared->SetLineColor(kBlack);
    smeared->SetMarkerColor(kBlack);
    smeared->SetMarkerSize(0.5);
    smeared->SetMarkerStyle(kOpenCircle);

    if(smearBestFit){
      bestfit->SetLineWidth(2);
      bestfit->SetLineColor(kRed);
      bestfit->SetMarkerColor(kRed);
      bestfit->SetMarkerSize(0.5);
      bestfit->SetMarkerStyle(kOpenCircle);
    }

    unsmeared->Draw("h");
    smeared->Draw("ep same");
    if(smearBestFit) bestfit->Draw("p same");

    legS->Draw();
  }

// Draw Best Fit Plot
  if(drawBestFitPlot){
    bestFits->SetMarkerColor(kYellow);	bestFits->SetMarkerStyle(7);
    printf("\nDrawing best fits...\n");

    TCanvas* c9 = new TCanvas("c3","BestFits",700,700);
    c9->SetLeftMargin(.15);
    c9->SetBottomMargin(.15);
    c9->SetTopMargin(.05);
    c9->SetRightMargin(.05);
    c9->SetLogx();
    c9->SetLogy();

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

    bestFits->Draw("psame");
    signalPt->Draw("psame");
    c9->RedrawAxis();

    std::string det_str = "#splitline{";
    if(use100m && !use700m) det_str += "LAr1-ND (1.6 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}{";    
    else if(use100m && use700m) det_str +=  "LAr1-ND (6.6 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}{";    
    else det_str += "}{";                                          
    if(use470m && !use700m && use100m) det_str += "and MicroBooNE (8.2 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}";
    else if(use470m && use700m) det_str += "MicroBooNE (1.3 #times 10#lower[-0.5]{#scale[0.75]{21}} POT) and T600 (6.6 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}";
    else if(!use470m && use700m) det_str += "and T600 (6.6 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}";
    else if(use470m && !use700m && !use100m) det_str += "MicroBooNE (6.6 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}";
    else det_str += "}";

    TLatex *tex_Detector = new TLatex(.2,.23,det_str.c_str());
    tex_Detector->SetNDC();
    tex_Detector->SetTextFont(62);
    tex_Detector->SetTextSize(0.025);
    tex_Detector->Draw(); 

    TLatex *tex_pre = new TLatex(.18,.74,"PRELIMINARY");
    tex_pre->SetNDC();
    tex_pre->SetTextFont(62);
    tex_pre->SetTextColor(kRed-3);
    tex_pre->SetTextSize(0.03);
    tex_pre->Draw();

    TLatex *tex_mode = new TLatex(.18,.92,"#nu mode, CC Events");
    tex_mode->SetNDC();
    tex_mode->SetTextFont(62);
    tex_mode->SetTextSize(0.025);
    tex_mode->Draw();

    TLatex *tex_un = new TLatex(.18,.89,"Statistical and Flux Uncertainty");
    tex_un->SetNDC();
    tex_un->SetTextFont(62);
    tex_un->SetTextSize(0.025);
    tex_un->Draw();

    TLatex *tex_E = new TLatex(.18,.86,"Reconstructed Energy");
    tex_E->SetNDC();
    tex_E->SetTextFont(62);
    tex_E->SetTextSize(0.025);
    tex_E->Draw();

    TLatex *tex_eff = new TLatex(.18,.83,"80% #nu#lower[0.4]{#mu} Efficiency");
    tex_eff->SetNDC();
    tex_eff->SetTextFont(62);
    tex_eff->SetTextSize(0.025);
    tex_eff->Draw();

    std::string str_signal;
    if(shape_only) str_signal = "Shape Only Analysis";
    else str_signal = "Shape and Rate Analysis";
    TLatex *tex_Signal = new TLatex(.18,.785,str_signal.c_str());
    tex_Signal->SetNDC();
    tex_Signal->SetTextFont(62);
    tex_Signal->SetTextSize(0.025);
    tex_Signal->Draw();

    c9->Print("bestFitTest500.png");
  }

  return 0;
}

int update_cov(int nbinsE, int nL, bool useBestFit){

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
	    if(!useBestFit)  M(Erri, Errj) *= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];		
	    else M(Erri, Errj) *= (NULLVec[Lrow][Erow] - (OscVec[Lrow][bestFitReso[1]][Erow] * getSin22THpt(bestFitReso[0],false)))*(NULLVec[Lcol][Ecol] - (OscVec[Lcol][bestFitReso[1]][Ecol] * getSin22THpt(bestFitReso[0],false)));
	    // Add Stat. Errors
	    if(Erri == Errj){ 
	      M (Erri, Errj) += CV(Lrow * nbinsE + Ecol,0);
	    }
	  }
	  n = 0;
	  M1 (Erri, Errj) = M (Erri, Errj);
	  Errj++;
	}
      }

      // Fill vectors for chi^2 calaculations
      if(useBestFit){
        N(Erri,0) = NULLVec[Lrow][Erow] - (OscVec[Lrow][bestFitReso[1]][Erow] * getSin22THpt(bestFitReso[0],false));
        NT += NULLVec[Lrow][Erow] - (OscVec[Lrow][bestFitReso[1]][Erow] * getSin22THpt(bestFitReso[0],false));
      }
      else{
        N(Erri,0) = (NULLVec[Lrow][Erow]);
        NT += (NULLVec[Lrow][Erow]);
      }

      if(drawSmear){
	if(Lrow == 1){
	  bestfit->SetBinContent(Erow+1,(NULLVec[Lrow][Erow] - (OscVec[Lrow][bestFitReso[1]][Erow] * getSin22THpt(bestFitReso[0],false)))/(bestfit->GetXaxis()->GetBinWidth(Erow+1)));
	}
      }

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
	    if(!useBestFit)  M(Erri, Errj) *= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];		
	    else M(Erri, Errj) *= (NULLVec[Lrow][Erow] - (OscVec[Lrow][bestFitReso[1]][Erow] * getSin22THpt(bestFitReso[0],false)))*(NULLVec[Lcol][Ecol] - (OscVec[Lcol][bestFitReso[1]][Ecol] * getSin22THpt(bestFitReso[0],false)));
	    if(Erri == Errj){ 
	      M (Erri, Errj) += CV(Lrow * nbinsE + Ecol,0);
	    }
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

double chisq_calc(int nL, int nbinsE){

  if(drawScale){
    double bins[20] = {.200, .300, .400, .450, .500, .550, .600, .650, .700, .750, .800, .850, .900, .950, 1.000, 1.250, 1.500, 2.000, 2.500, 3.000};
    cvHist = new TH1D("cvHist", "", nbinsE, bins);
    predHist = new TH1D("predHist","",nbinsE,bins);
    predScaledHist = new TH1D("predScaleHist","",nbinsE,bins);
  }

// Now, calculate chisq
  chisqSurf = new TGraph2D();
  int chisqSurf_counter = 0;
  Double_t chisq;
  double chiLow = -1;
  double bestFit_sin22th, bestFit_dm2;

  for(int dm2 = 0; dm2 <= npoints; dm2++){
    for(int sint = 0; sint <= npoints; sint++){
      chisq = 0.;
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

      if(drawScale && dm2 == predSigReso[1] && sint == predSigReso[0]){
        for(int eBins = 0; eBins < nbinsE; eBins ++){
          cvHist->SetBinContent(eBins,(CV (((nL-1)*nbinsE) + eBins, 0)/(cvHist->GetXaxis()->GetBinWidth(eBins+1))));
          predHist->SetBinContent(eBins,(Pred[dm2][sint] (((nL-1)*nbinsE) + eBins, 0))/(cvHist->GetXaxis()->GetBinWidth(eBins+1)));
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

      chisqArray[dm2*(npoints+1)+sint] = chisq;

      if(chisq <= 100.) chisqSurf->SetPoint(chisqSurf_counter,getSin22THpt(sint,false),getDMpt(dm2,false),sqrt(chisq));
      else chisqSurf->SetPoint(chisqSurf_counter,getSin22THpt(sint,false),getDMpt(dm2,false),10.);
      chisqSurf_counter++;

      // Now, try to find best fit
      if(chiLow < 0) chiLow = chisq;
      if(chisq < chiLow){
        chiLow = chisq;
        bestFit_sin22th = getSin22THpt(sint,false);
	bestFitReso[0] = sint;
        bestFit_dm2 = getDMpt(dm2,false);
	bestFitReso[1] = dm2;
      }
    }
  }

  bestFit = new TMarker(bestFit_sin22th,bestFit_dm2, 3);

  std::cout << "Calculated Chisq Points" << std::endl;
  return chiLow;
}

int draw_jellybeans(double chiLow, bool exclude){

// All we need for this is the lowest chisq point, chisqSurf and bestFit
  fiveSigma = new TGraph();	threeSigma = new TGraph();	ninety = new TGraph();

  int fiveSCounter = 0;  int threeSCounter = 0;  int ninetyCounter = 0;
  int fiveExcludeN = 0, fiveExcludeS = 0, fiveExcludeW = 0, fiveExcludeE = 0;
  int threeExcludeN = 0, threeExcludeS = 0, threeExcludeW = 0, threeExcludeE = 0;

  double chisq;

  for(int dm2 = 0; dm2 < npoints; dm2 ++){
    for(int sint = 0; sint < npoints; sint ++){
      chisq = chisqArray[dm2 * (npoints + 1) + sint];

      if(chisq <= 4.61 + chiLow){
	ninety->SetPoint(ninetyCounter,getSin22THpt(sint,false),getDMpt(dm2,false));
	ninetyCounter++;
      }
      if(chisq <= 11.83 + chiLow){
        threeSigma->SetPoint(threeSCounter,getSin22THpt(sint, false),getDMpt(dm2, false));
        if(exclude){
          if(dm2 == 0)	threeExcludeS = 1;
	  if(dm2 == npoints-1) threeExcludeN = 1;
	  if(sint == 0)	threeExcludeW = 1;
	  if(sint == npoints-1) threeExcludeE = 1;
        }
        threeSCounter++;
      }
      if(chisq <= 28.23 + chiLow) {
	fiveSigma->SetPoint(fiveSCounter,getSin22THpt(sint, false),getDMpt(dm2, false));
	if(exclude){
	  if(dm2 == 0)	fiveExcludeS = 1;
	  if(dm2 == npoints-1) fiveExcludeN = 1;
	  if(sint == 0)	fiveExcludeW = 1;
	  if(sint == npoints-1) fiveExcludeE = 1;
	}
	fiveSCounter ++;
      }
    }
  }

// Begin Drawing
  fiveSigma->SetMarkerColor(kBlue);	fiveSigma->SetFillColor(kBlue);			fiveSigma->SetMarkerStyle(7);
  threeSigma->SetMarkerStyle(7);	threeSigma->SetMarkerColor(kGreen + 1);		threeSigma->SetFillColor(kGreen + 1);
  ninety->SetMarkerColor(kYellow);  	ninety->SetFillColor(kYellow);			ninety->SetMarkerStyle(7);

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

  if(!exclude || (exclude && fiveExcludeE + fiveExcludeN + fiveExcludeW + fiveExcludeS < 3))		fiveSigma->Draw("p");
  if(!exclude || (exclude && threeExcludeE + threeExcludeN + threeExcludeW + threeExcludeS < 3))	threeSigma->Draw("psame");
  ninety->Draw("psame");
  bestFit->Draw("psame");
  signalPt->Draw("psame");
  c3->RedrawAxis();

  std::string det_str = "#splitline{";
  if(use100m && !use700m) det_str += "LAr1-ND (1.6 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}{";    
  else if(use100m && use700m) det_str +=  "LAr1-ND (6.6 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}{";    
  else det_str += "}{";                                          
  if(use470m && !use700m && use100m) det_str += "and MicroBooNE (8.2 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}";
  else if(use470m && use700m) det_str += "MicroBooNE (1.3 #times 10#lower[-0.5]{#scale[0.75]{21}} POT) and T600 (6.6 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}";
  else if(!use470m && use700m) det_str += "and T600 (6.6 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}";
  else if(use470m && !use700m && !use100m) det_str += "MicroBooNE (6.6 #times 10#lower[-0.5]{#scale[0.75]{20}} POT)}";
  else det_str += "}";


  TLatex *tex_Detector = new TLatex(.2,.23,det_str.c_str());
  tex_Detector->SetNDC();
  tex_Detector->SetTextFont(62);
  tex_Detector->SetTextSize(0.025);
  tex_Detector->Draw(); 

  TLatex *tex_pre = new TLatex(.18,.74,"PRELIMINARY");
  tex_pre->SetNDC();
  tex_pre->SetTextFont(62);
  tex_pre->SetTextColor(kRed-3);
  tex_pre->SetTextSize(0.03);
  tex_pre->Draw();

  TLatex *tex_mode = new TLatex(.18,.92,"#nu mode, CC Events");
  tex_mode->SetNDC();
  tex_mode->SetTextFont(62);
  tex_mode->SetTextSize(0.025);
  tex_mode->Draw();

  TLatex *tex_un = new TLatex(.18,.89,"Statistical and Flux Uncertainty");
  tex_un->SetNDC();
  tex_un->SetTextFont(62);
  tex_un->SetTextSize(0.025);
  tex_un->Draw();

  TLatex *tex_E = new TLatex(.18,.86,"Reconstructed Energy");
  tex_E->SetNDC();
  tex_E->SetTextFont(62);
  tex_E->SetTextSize(0.025);
  tex_E->Draw();

  TLatex *tex_eff = new TLatex(.18,.83,"80% #nu#lower[0.4]{#mu} Efficiency");
  tex_eff->SetNDC();
  tex_eff->SetTextFont(62);
  tex_eff->SetTextSize(0.025);
  tex_eff->Draw();

  std::string str_signal;
  if(shape_only) str_signal = "Shape Only Analysis";
  else str_signal = "Shape and Rate Analysis";
  TLatex *tex_Signal = new TLatex(.18,.785,str_signal.c_str());
  tex_Signal->SetNDC();
  tex_Signal->SetTextFont(62);
  tex_Signal->SetTextSize(0.025);
  tex_Signal->Draw();

  TLegend* legt=new TLegend(0.2,0.3,0.6,0.45);
  legt->SetFillStyle(0);
  legt->SetFillColor(0);
  legt->SetBorderSize(0);
  legt->SetTextFont(62);
  legt->SetTextSize(0.03);
  if(!exclude || (exclude && fiveExcludeE + fiveExcludeN + fiveExcludeW + fiveExcludeS < 3))		legt->AddEntry(fiveSigma,"5#sigma Confidence","f");
  if(!exclude || (exclude && threeExcludeE + threeExcludeN + threeExcludeW + threeExcludeS < 3))	legt->AddEntry(threeSigma,"3#sigma Confidence","f");
  legt->AddEntry(ninety,"90% Confidence","f");
  legt->AddEntry(signalPt,"Signal","fp");
  legt->Draw();

  return 0; 
}

int build_vectors(int nL, int nbinsE, double signalSin22th, double signalDm2, bool smear){

  if(drawSmear){
    double bins[20] = {.200, .300, .400, .450, .500, .550, .600, .650, .700, .750, .800, .850, .900, .950, 1.000, 1.250, 1.500, 2.000, 2.500, 3.000};
    smeared = new TH1D("smeared", "", nbinsE, bins);
    unsmeared = new TH1D("unsmeared","",nbinsE,bins);
    bestfit = new TH1D("bestfit","",nbinsE,bins);
  }

  Pred.resize(npoints+1);
  for(int i = 0; i < npoints+1; i++){
    Pred[i].resize(npoints+1);
    for(int j = 0; j < npoints+1; j++){
      Pred[i][j].ResizeTo(nbinsE*nL,1);
    }
  }

  int Predi = 0;
  std::cout << "Filling Prediction Matrices..." << std::endl;

  signalReso[0] = 0;	signalReso[1] = 0;
  if(drawScale){
    predSigReso[0] = 0;	predSigReso[1] = 0;
  }
  for(int dm = 0; dm <= npoints; dm++){
    if(getDMpt(dm, false) <= signalDm2)	signalReso[1] = dm;
    if(drawScale && getDMpt(dm,false) <= predSig[1]) predSigReso[1] = dm;
    for(int s = 0; s <= npoints; s++){
      if(getSin22THpt(s,false) <= signalSin22th && s > signalReso[0]) signalReso[0] = s;
      if(drawScale && (getSin22THpt(s,false) <= predSig[0] && s > predSigReso[0]))  predSigReso[0] = s;	
      Predi = 0; 
      for(int Lbin = 0; Lbin < nL; Lbin++){
        for(int Ebin = 0; Ebin < nbinsE; Ebin++){
  	  Pred[dm][s] (Predi,0) = NULLVec[Lbin][Ebin] - (OscVec[Lbin][dm][Ebin])*(getSin22THpt(s,false));
	  Predi++;
	}
      }
    }
  }

  signalPt = new TMarker(getSin22THpt(signalReso[0],false),getDMpt(signalReso[1],false),7);
  signalPt->SetMarkerColor(kRed-3);

  CV.ResizeTo(nbinsE*nL,1);	CV.Zero();
  TMatrixT <float> CV_noSmear(nbinsE*nL,1);
  for(int Lbin = 0; Lbin < nL; Lbin ++){
    for(int Ebin = 0; Ebin < nbinsE; Ebin ++){
      CV_noSmear(Lbin*nbinsE + Ebin,0) = (NULLVec[Lbin][Ebin] - (OscVec[Lbin][signalReso[1]][Ebin])*signalPt->GetX());
      if(smear){
        if(Lbin == 0 && Ebin == 0)std::cout << "... Smearing based on Statistics ..." << std::endl;

	TRandom3 * r3 = new TRandom3(0);
	CV(Lbin*nbinsE + Ebin,0) = r3->Poisson(CV_noSmear(Lbin*nbinsE + Ebin, 0)); 
	  //r3->Gaus(CV_noSmear(Lbin*nbinsE + Ebin, 0),sqrt(CV_noSmear(Lbin*nbinsE + Ebin,0)));

        if(Lbin == 1 && drawSmear){
	  smeared->SetBinContent(Ebin+1,CV(Lbin*nbinsE + Ebin,0)/(smeared->GetXaxis()->GetBinWidth(Ebin+1)));
	  smeared->SetBinError(Ebin+1,sqrt(CV(Lbin*nbinsE + Ebin,0))/(smeared->GetXaxis()->GetBinWidth(Ebin+1)));
	  unsmeared->SetBinContent(Ebin+1,CV_noSmear(Lbin*nbinsE + Ebin, 0)/(smeared->GetXaxis()->GetBinWidth(Ebin+1)));
	}
      }
      else {
        CV(Lbin*nbinsE + Ebin,0) = CV_noSmear(Lbin*nbinsE + Ebin, 0);
      }
    }
  }
  std::cout << "...Prediction and Signal Matrices Filled" << std::endl;
  return 0;
}

double getDMpt(int dm, bool inverse){
  if(inverse)	return floor((TMath::Log10(dm) - TMath::Log10(dm2min))*npoints/(TMath::Log10(dm2max/dm2min)));
  else	return pow(10.,(TMath::Log10(dm2min)+ (dm * (1./npoints))*TMath::Log10(dm2max/dm2min)));
}

double getSin22THpt(int sint, bool inverse){
  if(inverse) return floor((TMath::Log10(sint) - TMath::Log10(sin22thmin))*npoints/(TMath::Log10(sin22thmax/sin22thmin)));
  else return pow(10., (TMath::Log10(sin22thmin)+ (sint * (1./npoints))*TMath::Log10(sin22thmax/sin22thmin)));
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

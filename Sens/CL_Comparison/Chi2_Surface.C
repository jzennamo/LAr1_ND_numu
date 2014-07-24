#include <sys/stat.h>
#include <unistd.h>
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TROOT.h" //for gROOT
#include "TStyle.h" //for gStyle
#include "TH2D.h"
#include "TH3D.h"
#include "TLegend.h"
#include "THStack.h"
#include "TImage.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TGraph2D.h"
#include <iostream>
#include <ostream>
#include <sstream>
#include <fstream>

int multiple_detector_fit()
{
/* ------------------------------ */

  // Fix DM2 point for slice plot (otherwise set to negative)
  int fixedDM = -1;

  bool shape_only = false;

  // Select one of the following modifiers
  bool statmatch = false;
  bool noflux_statmatch = false;
  bool full = true;

  // Select any (or all) detector setups
  // Note that multiple choices are only really significant for a fixed DM, since they'll all meld together in a 3d model
  bool use100 = true;
  bool use150 = false;
  bool use200 = false; 

/* ------------------------------ */


  TNtuple * thisNtuple[3];
  int counter = 0;


  // First, build the filename and import the file(s)
  std::string shapeName;
  std::string mod;

  if(shape_only)	shapeName = "ShapeOnly";
  else    shapeName = "Shape_and_Rate";

  if(statmatch)    mod = "_StatMatch.root";
  else if(noflux_statmatch)    mod = "_NoFlux_StatMatch.root";
  else if(full)    mod = ".root";

  if(use100){
    std::string name = "output/chi2_Surface_LAr1ND_100m_T600_";
    name = name + shapeName + mod;
    TFile *f = new TFile(name.c_str());	
    f->GetObject("chi2",thisNtuple[counter]);
    counter ++;
  }
  if(use150){
    std::string name = "output/chi2_Surface_LAr1ND_150m_T600_";
    name = name + shapeName + mod;
    TFile *f = new TFile(name.c_str());	
    f->GetObject("chi2",thisNtuple[counter]);
    counter ++;
  }
  if(use200){
    std::string name = "output/chi2_Surface_LAr1ND_200m_T600_";
    name = name + shapeName + mod;
    TFile *f = new TFile(name.c_str());	
    f->GetObject("chi2",thisNtuple[counter]);
    counter ++;
  }
 
  Float_t dm2_temp;
  Float_t sin22th_temp;
  Float_t chi2_temp;

  //extract the values from the ntuple
  int i_entry;

  // Initialize the surface plots
  TGraph2D *g[3];
  TGraph2D *ThreeSigma[3];
  TGraph2D *OneSigma[3];

  // Initialize the cross sectional plots
  TGraph *g2D[3];
  TGraph *ThreeSigma2D[3];
  TGraph *OneSigma2D[3];



for(int i = 0; i < counter; i++)				// Cycle through each detector setup
{
  thisNtuple[i]->SetBranchAddress("chisq", &chi2_temp);
  thisNtuple[i]->SetBranchAddress("dm2", &dm2_temp);
  thisNtuple[i]->SetBranchAddress("sin22th", &sin22th_temp);

  if(fixedDM > 0)
  {
    g2D[i] = new TGraph();
    ThreeSigma2D[i] = new TGraph();
    OneSigma2D[i] = new TGraph();
  }

  g[i] = new TGraph2D();
  ThreeSigma[i] = new TGraph2D();
  OneSigma[i] = new TGraph2D();
  
  i_entry = 0;

  std::cout<< "Entries : " << thisNtuple[i]->GetEntries() << std::endl;

  while (i_entry < thisNtuple[i]->GetEntries())			// Cycle through all entries
    {

      //select a point in the ntuple
      thisNtuple[i]->GetEntry(i_entry);

      //store them in the TGraph for plotting
      if(sqrt(chi2_temp) > 3)	ThreeSigma[i]->SetPoint(i_entry, sin22th_temp, dm2_temp, 3);
      else	ThreeSigma[i]->SetPoint(i_entry, sin22th_temp, dm2_temp, sqrt(chi2_temp));

      if(sqrt(chi2_temp) > 1)	OneSigma[i]->SetPoint(i_entry, sin22th_temp, dm2_temp, 1);
      else	OneSigma[i]->SetPoint(i_entry, sin22th_temp, dm2_temp, sqrt(chi2_temp));

      if(sqrt(chi2_temp) > 5)	g[i]->SetPoint(i_entry, sin22th_temp, dm2_temp, sqrt(25));
      else	g[i]->SetPoint(i_entry, sin22th_temp, dm2_temp, sqrt(chi2_temp));

      if(fixedDM > 0 && dm2_temp == fixedDM){	

	if(sqrt(chi2_temp) > 3)	ThreeSigma2D[i]->SetPoint(i_entry, sin22th_temp, 3);
        else	ThreeSigma2D[i]->SetPoint(i_entry, sin22th_temp, sqrt(chi2_temp));

        if(sqrt(chi2_temp) > 1)	OneSigma2D[i]->SetPoint(i_entry, sin22th_temp, 1);
        else	OneSigma2D[i]->SetPoint(i_entry, sin22th_temp, sqrt(chi2_temp));

        if(sqrt(chi2_temp) > 5)	g2D[i]->SetPoint(i_entry, sin22th_temp, 5);
        else	g2D[i]->SetPoint(i_entry, sin22th_temp, sqrt(chi2_temp));
      }

      i_entry++;
    }
}
//cout << "Counter = " << counter;
//cout << "...[0] points = " << g2D[0]->GetN();
//cout << "...[1] points = " << g2D[1]->GetN();
//cout << "...[2] points = " << g2D[2]->GetN();




  // Draw the surface plot
  if(fixedDM < 0 || counter == 1){
    std::cout << "Drawing ..." << std::endl; 
    gStyle->SetOptStat(0000);
    gStyle->SetOptFit(0000);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);
    
    TCanvas* c1 = new TCanvas("c1", "", 800, 800);     
    c1->SetLeftMargin(.15);
    c1->SetBottomMargin(.15);
    c1->SetTopMargin(.05);
    c1->SetRightMargin(.05);
    c1->SetLogy();
    c1->SetLogx();
    c1->cd();
    
    OneSigma[0]->SetTitle("");
    OneSigma[0]->GetXaxis()->SetTitle("sin^{2}2#theta#lower[0.4]{#mu#kern[-0.3]{#mu}}");
    OneSigma[0]->GetXaxis()->CenterTitle(true);
    OneSigma[0]->GetXaxis()->SetLabelFont(62);
    OneSigma[0]->GetXaxis()->SetLabelOffset(0.003);
    OneSigma[0]->GetXaxis()->SetLabelSize(0.03);
    OneSigma[0]->GetXaxis()->SetTitleSize(0.05);
    OneSigma[0]->GetXaxis()->SetTitleFont(62);
    OneSigma[0]->GetXaxis()->SetTitleOffset(1.75);
    
    OneSigma[0]->GetYaxis()->SetTitle("#Deltam_{41}^{2} [eV^{2}]");
    OneSigma[0]->GetYaxis()->CenterTitle(true);
    OneSigma[0]->GetYaxis()->SetLabelFont(62);
    OneSigma[0]->GetYaxis()->SetLabelSize(0.03);
    OneSigma[0]->GetYaxis()->SetTitleSize(0.05);
    OneSigma[0]->GetYaxis()->SetTitleFont(62);
    OneSigma[0]->GetYaxis()->SetTitleOffset(2.0);
    
    OneSigma[0]->GetZaxis()->SetTitle("Confidence Level [#sigma]");
    OneSigma[0]->GetZaxis()->CenterTitle(true);
    OneSigma[0]->GetZaxis()->SetLabelFont(62);
    OneSigma[0]->GetZaxis()->SetLabelSize(0.03);
    OneSigma[0]->GetZaxis()->SetTitleSize(0.05);
    OneSigma[0]->GetZaxis()->SetTitleFont(62);
    OneSigma[0]->GetZaxis()->SetNdivisions(506);
    OneSigma[0]->GetZaxis()->SetTitleOffset(1.25);
    OneSigma[0]->GetZaxis()->SetRangeUser(0.0000001, 5.5);

    for(int j = 0; j < counter; j ++)
    {
      OneSigma[j]->SetMarkerColor(kRed);
      if(j == 0) OneSigma[0]->Draw("P");
      else OneSigma[j]->Draw("Psame");
    
      ThreeSigma[j]->SetMarkerColor(kBlue);
      ThreeSigma[j]->Draw("Psame");

      g[j]->Draw("Psame");
    }
  }


if(fixedDM > 0){

  TCanvas *c2 = new TCanvas("c2","Fixed DM2",800,800);

  TMultiGraph *mg = new TMultiGraph();

  c2->SetLeftMargin(.15);
  c2->SetBottomMargin(.15);
  c2->SetTopMargin(.05);
  c2->SetRightMargin(.05);
  c2->SetLogx();
  c2->cd();

  OneSigma2D[0]->SetTitle("");
  OneSigma2D[0]->GetXaxis()->SetTitle("sin^{2}2#theta#lower[0.4]{#mu#kern[-0.3]{#mu}}");
  OneSigma2D[0]->GetXaxis()->CenterTitle(true);
  OneSigma2D[0]->GetXaxis()->SetLabelFont(62);
  OneSigma2D[0]->GetXaxis()->SetLabelOffset(0.003);
  OneSigma2D[0]->GetXaxis()->SetLabelSize(0.03);
  OneSigma2D[0]->GetXaxis()->SetTitleSize(0.05);
  OneSigma2D[0]->GetXaxis()->SetTitleFont(62);
  OneSigma2D[0]->GetXaxis()->SetTitleOffset(1.75);

  OneSigma2D[0]->GetYaxis()->SetTitle("Confidence Level [#sigma]");
  OneSigma2D[0]->GetYaxis()->CenterTitle(true);
  OneSigma2D[0]->GetYaxis()->SetLabelFont(62);
  OneSigma2D[0]->GetYaxis()->SetLabelSize(0.03);
  OneSigma2D[0]->GetYaxis()->SetTitleSize(0.05);
  OneSigma2D[0]->GetYaxis()->SetTitleFont(62);
  OneSigma2D[0]->GetYaxis()->SetNdivisions(506);
  OneSigma2D[0]->GetYaxis()->SetTitleOffset(1.25);
  OneSigma2D[0]->GetYaxis()->SetRangeUser(0.0000001, 5.5);

cout << "drawing...";


    OneSigma2D[0]->SetLineColor(kRed);
    mg->Add(OneSigma2D[0]);
    ThreeSigma2D[0]->SetLineColor(kBlue);
    mg->Add(ThreeSigma2D[0]);

    for(int j = 0; j < counter; j ++)
    {
      if(j == 1) g2D[j]->SetLineWidth(2);
      if(j == 2) g2D[j]->SetLineWidth(3);
      mg->Add(g2D[j]);
    }
    
    mg->Draw("AL");

  }

    return 0;
    
}

#ifndef __CINT__
int main()
{
  multiple_detector_fit();
  return 0;
}
# endif

void Chi2_Surface(){
  multiple_detector_fit();
  return;
}

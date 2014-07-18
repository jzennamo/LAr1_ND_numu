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


  //Choose input file you want
  TFile *f = new TFile("output/chi2_Surface_LAr1ND_100m_T600_ShapeOnly.root");

  //Read in the ntuple that you want to study
  TNtuple * thisNtuple;
  f->GetObject("chi2", thisNtuple);
  
  Float_t dm2_temp;
  Float_t sin22th_temp;
  Float_t chi2_temp;

  int npoints = 500;

  // configure the dm2, sin22th points:
  double dm2min(0.01), dm2max(100.0);
  double sin22thmin(0.001), sin22thmax(1.0);  

  //extract the values from the ntuple
  //The values are tied together as an entry such that each point has a dm2, sin22theta, and 
  thisNtuple->SetBranchAddress("chisq", &chi2_temp);
  thisNtuple->SetBranchAddress("dm2", &dm2_temp);
  thisNtuple->SetBranchAddress("sin22th", &sin22th_temp);

  int i_entry = 0;
  int i_dm2 = 0;
  int i_sin22th = 0;

  std::cout<< "Entries : " << thisNtuple->GetEntries() << std::endl;
  
  TGraph2D *g = new TGraph2D();

  //Cycle through each entry in the ntuple
  while (i_entry < thisNtuple->GetEntries())
    {

      //select a point in the ntuple
      thisNtuple->GetEntry(i_entry);

      //store them in the TGraph for plotting
      if(sqrt(chi2_temp) > 5){    //cut off at 5\sigma for presentation
      	g->SetPoint(i_entry, sin22th_temp, dm2_temp, sqrt(25));
      }
      else{
        g->SetPoint(i_entry, sin22th_temp, dm2_temp, sqrt(chi2_temp));
      }

      i_entry++;
    }
    
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
    
    g->SetTitle("");
    g->GetXaxis()->SetTitle("sin^{2}2#theta#lower[0.4]{#mu#kern[-0.3]{#mu}}");
    g->GetXaxis()->CenterTitle(true);
    g->GetXaxis()->SetLabelFont(62);
    g->GetXaxis()->SetLabelOffset(0.003);
    g->GetXaxis()->SetLabelSize(0.03);
    g->GetXaxis()->SetTitleSize(0.05);
    g->GetXaxis()->SetTitleFont(62);
    g->GetXaxis()->SetTitleOffset(1.75);
    
    g->GetYaxis()->SetTitle("#Deltam_{41}^{2} [eV^{2}]");
    g->GetYaxis()->CenterTitle(true);
    g->GetYaxis()->SetLabelFont(62);
    g->GetYaxis()->SetLabelSize(0.03);
    g->GetYaxis()->SetTitleSize(0.05);
    g->GetYaxis()->SetTitleFont(62);
    g->GetYaxis()->SetTitleOffset(2.0);
    
    g->GetZaxis()->SetTitle("Confidence Level [#sigma]");
    g->GetZaxis()->CenterTitle(true);
    g->GetZaxis()->SetLabelFont(62);
    g->GetZaxis()->SetLabelSize(0.03);
    g->GetZaxis()->SetTitleSize(0.05);
    g->GetZaxis()->SetTitleFont(62);
    g->GetZaxis()->SetNdivisions(506);
    g->GetZaxis()->SetTitleOffset(1.25);

    //Using "P" instead of "CONT" because it is so resource heavy it makes ROOT crash
    //Feel free to try it though! 
    g->Draw("P");

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

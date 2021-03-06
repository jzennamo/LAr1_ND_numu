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


  double NULLVec[100];
  double OscVec[1001][7][100];

  double MCVec[2][20];
  double DataVec[2][20];


  int nbinsE = 0;
  int Entries = 0;
  if (use470m){
    std::string temp_name = "../MatrixFiles/combined_ntuple_600m_onaxis_nu_processed_numu.root";

    TFile temp_file(temp_name.c_str());
    TH1D *NULL_470;
    TH1D *MC_470;
    NULL_470 = (TH1D*)(temp_file.Get("NumuCC"));
    MC_470 = (TH1D*)(temp_file.Get("NumuCC_MC"));
    NULL_470->Rebin(1);
    nbinsE = NULL_470->GetNbinsX();
    //    std::cout << "Entries in bin 1 : " << NULL_470->GetNEntries(1) << std::endl;
    // std::cout << "Events in bin 1 : " << NULL_470->GetBinContent(1) << std::endl;

    for(int i = 0; i <= nbinsE; i++){
      NULLVec[i] = (NULL_470->GetBinContent(i));
      MCVec[0][i] = (MC_470->GetBinContent(i));
      DataVec[0][i] = (NULL_470->GetBinContent(i));
    }

    for(int u = 0; u < npoints; u++){
      for(int s = 0; s < 7; s++){
	TH1D *OSC_470;
	std::string upoint = std::to_string(u);
	std::string name = "Universe_";
	std::string name2 = "_MultiSim_";
	std::string mul = std::to_string(s);
	
	name = name+upoint+name2+mul;	
       
	OSC_470 = (TH1D*)(temp_file.Get(name.c_str()));
	OSC_470->Rebin(1);
	for(int i = 0; i <= nbinsE; i++){
	  OscVec[u][s][i] = (OSC_470->GetBinContent(i));
	}
	delete OSC_470;
      }
    }

    delete MC_470;
    delete NULL_470;
    temp_file.Close();
  }

  if (use100m){

    std::string temp_name = "../MatrixFiles/combined_ntuple_200m_nu_processed_numu.root";

    TFile temp_file(temp_name.c_str());
    TH1D *NULL_100;
    TH1D *MC_100;
    NULL_100 = (TH1D*)(temp_file.Get("NumuCC"));
    MC_100 = (TH1D*)(temp_file.Get("NumuCC_MC"));
    NULL_100->Rebin(1);

    for(int i = 0; i <= nbinsE; i++){
      NULLVec[i] /= (NULL_100->GetBinContent(i));
      MCVec[1][i] = (MC_100->GetBinContent(i));
      DataVec[1][i] = (NULL_100->GetBinContent(i));
    }

    for(int u = 0; u < npoints; u++){
      for(int s = 0; s < 7; s++){
	TH1D *OSC_100;
	std::string upoint = std::to_string(u);
	std::string name = "Universe_";
	std::string name2 = "_MultiSim_";
	std::string mul = std::to_string(s);
	
	name = name+upoint+name2+mul;	
       
	OSC_100 = (TH1D*)(temp_file.Get(name.c_str()));
	OSC_100->Rebin(1);
	for(int i = 0; i <= nbinsE; i++){
	  OscVec[u][s][i] /= (OSC_100->GetBinContent(i));
	}
	delete OSC_100;
      }
    }
    delete MC_100;
    delete NULL_100;
    temp_file.Close();
  }

  double bins[20] = {.200, .300, .400, .450, .500, .550, .600, .650, .700, .750, .800, .850, .900, .950, 1.000, 1.250, 1.500, 2.000, 2.500, 3.000};

  TH1D *ratio = new TH1D("ratio",";;",19,bins);
  TH1D *RMS = new TH1D("RMS",";;",19,bins);
  TH1D *Err = new TH1D("Err",";;",19,bins);
  TH1D *Err_data = new TH1D("Err_data",";;",19,bins);
  
  std::vector <std::vector < TH1D*> >  ratio_Syst;
  ratio_Syst.resize(npoints+1);


  for(int i = 0; i <= npoints; i++){
    //    ratio_Syst[i].resize(7);
    for(int j = 0; j < 7; j++){
      ratio_Syst[i].push_back((TH1D*)ratio->Clone());
    }
   }


  TCanvas* c6 = new TCanvas("c6","",700,700);
  c6->SetLeftMargin(.15);
  c6->SetBottomMargin(.1);
  c6->SetTopMargin(.05);
  c6->SetRightMargin(.05);
  c6->cd();

  double Avg[nbinsE]; 
  double N[nbinsE];
  double FracErr[nbinsE];
  double FracErr_data[nbinsE];

  for(int E = 0; E <= nbinsE; E++) {

    ratio->SetBinContent(E,NULLVec[E]);
    Avg[E] = 0;
    N[E] = 0;
    Avg[E] += NULLVec[E];
    N[E] += 1;
    FracErr[E] = sqrt(pow(sqrt(MCVec[0][E])/MCVec[0][E],2)+pow(sqrt(MCVec[1][E])/MCVec[1][E],2));
    FracErr_data[E] = sqrt(pow(sqrt(DataVec[0][E])/DataVec[0][E],2)+pow(sqrt(DataVec[1][E])/DataVec[1][E],2));

  }
  
  
  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(3);
  ratio->SetLineStyle(1);
  ratio->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  ratio->GetXaxis()->SetLabelSize(0.04);
  ratio->GetXaxis()->SetTitleFont(62);
  ratio->GetXaxis()->SetLabelFont(62);
  ratio->GetYaxis()->SetLabelFont(62);
  ratio->GetXaxis()->SetTitleSize(0.04);
  ratio->GetYaxis()->SetTitleFont(62);
  ratio->GetYaxis()->SetTitleSize(0.04);
  ratio->GetXaxis()->CenterTitle();
  ratio->GetYaxis()->CenterTitle();
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetTitle("T600 (600m, on axis) / LAr1-ND (200m)");
  ratio->GetYaxis()->SetTitleOffset(1.7);
  ratio->GetYaxis()->SetLabelSize(0.04);
  ratio->SetStats(0);
  ratio->GetYaxis()->SetRangeUser(0.41,1); 
  ratio->GetYaxis()->SetNdivisions(506);
  ratio->Draw("h");

  for(int u = 0; u < 1000; u++){
    for(int E = 0; E <= nbinsE; E++){
 
      ratio_Syst[u][6]->SetBinContent(E,OscVec[u][6][E]);
      Avg[E] += OscVec[u][6][E];
      N[E] += 1;
      
    }

    ratio_Syst[u][6]->SetLineColor(kGray+1);
    ratio_Syst[u][6]->SetLineWidth(1);
    ratio_Syst[u][6]->SetLineStyle(3);
    ratio_Syst[u][6]->Draw("h same");
  
  }

  TLegend* legt1=new TLegend(0.2,0.8,0.6,0.85);
  legt1->SetFillStyle(0);
  legt1->SetFillColor(0);
  legt1->SetBorderSize(0);
  legt1->SetTextFont(62);
  legt1->SetTextSize(0.03);
  legt1->AddEntry(ratio,"Nominal #nu#lower[0.4]{#mu} Event Ratio","l");
  legt1->Draw();


  double rms = 0;
  int n = 0;

    std::cout << "Product of MultiSim Weights" << std::endl;    

  for(int E = 0; E <= nbinsE; E++){

      Avg[E] /= N[E];

    for(int u = 0; u < 1000; u++){

      rms += pow((ratio->GetBinContent(E) - ratio_Syst[u][6]->GetBinContent(E)),2);
      n++;
  }

    rms /= (n); 

    RMS->SetBinContent(E,100*sqrt(rms)/(ratio->GetBinContent(E)));

    std::cout << " RMS for bin " << ratio->GetXaxis()->GetBinCenter(E) << " is " << RMS->GetBinContent(E) << std::endl; 


    Err->SetBinContent(E,100*FracErr[E]);
    Err_data->SetBinContent(E,100*FracErr_data[E]);
    rms = 0;
    n = 0;}

  /*  TLatex *Total = new TLatex(.15,.95,"Product of MultiSim Weights");
  Total->SetNDC();
  Total->SetTextFont(62);
  Total->SetTextSize(0.05);
  Total->Draw();  */

  //ratio->GetYaxis()->SetRangeUser(0,0.15); 
   ratio->Draw("h same");
   c6->Print("total_ratio_200m.pdf");
  

  TCanvas* c10 = new TCanvas("c10","",700,700);
  c10->SetLeftMargin(.1);
  c10->SetBottomMargin(.1);
  c10->SetTopMargin(.05);
  c10->SetRightMargin(.05);
  c10->cd();

  RMS->SetLineColor(kBlack);
  RMS->SetLineWidth(3);
  RMS->SetLineStyle(1);
  RMS->SetLineWidth(3);
  RMS->SetLineStyle(1);
  RMS->GetXaxis()->SetLabelSize(0.04);
  RMS->GetXaxis()->SetTitleFont(62);
  RMS->GetXaxis()->SetTitleSize(0.04);
  RMS->GetYaxis()->SetTitleFont(62);
  RMS->GetYaxis()->SetTitleSize(0.04);
  RMS->GetXaxis()->SetLabelFont(62);
  RMS->GetYaxis()->SetLabelFont(62);
  RMS->GetXaxis()->CenterTitle();
  RMS->GetYaxis()->CenterTitle();
  RMS->GetXaxis()->SetTitleOffset(1.2);
  RMS->GetYaxis()->SetTitleOffset(1.2);
  RMS->GetYaxis()->SetLabelSize(0.04);
  RMS->SetStats(0);
  RMS->GetYaxis()->SetRangeUser(0.11,0.4);
  RMS->GetYaxis()->SetNdivisions(506);
  RMS->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  RMS->GetYaxis()->SetTitle("Percent Uncertainty [%]");
  RMS->GetYaxis()->SetRangeUser(0.1,6);
 RMS->Draw("h");
 Err->SetLineColor(kGray);
 Err->SetLineStyle(7);
 Err->SetLineWidth(3);
 Err->Draw("h same");
 Err_data->SetLineColor(kBlack);
 Err_data->SetLineStyle(7);
 Err_data->SetLineWidth(3);
 Err_data->Draw("h same");

 TLegend* legt=new TLegend(0.2,0.7,0.6,0.85);
 legt->SetFillStyle(0);
 legt->SetFillColor(0);
 legt->SetBorderSize(0);
 legt->SetTextFont(62);
 legt->SetTextSize(0.03);
 legt->AddEntry(RMS,"RMS for Ratio","l");
 legt->AddEntry(Err,"MC Stat Uncertainty","l");
 legt->AddEntry(Err_data,"Data Stat Uncertainty","l");
 legt->Draw();

  c10->Print("total_rms_200m.pdf");

 /*

  TCanvas* c0 = new TCanvas("c0","",700,700);
  c0->SetLeftMargin(.1);
  c0->SetBottomMargin(.1);
  c0->SetTopMargin(.15);
  c0->SetRightMargin(.1);
  c0->cd();

  for(int E = 0; E <= nbinsE; E++) {

    Avg[E] = 0;
    N[E] = 0;
    Avg[E] += NULLVec[E];
    N[E] += 1;
}
  
   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h");

  for(int u = 0; u < 1000; u++){
    for(int E = 0; E <= nbinsE; E++){
 
      ratio_Syst[u][0]->SetBinContent(E,OscVec[u][0][E]);
      Avg[E] += OscVec[u][0][E];
      N[E] += 1;
     
    }

    ratio_Syst[u][0]->SetLineColor(kGray+1);
    ratio_Syst[u][0]->SetLineWidth(1);
    ratio_Syst[u][0]->SetLineStyle(3);
    ratio_Syst[u][0]->Draw("h same");
   
  }

  rms = 0;
  n = 0;

  std::cout << "MultiSim 0 Weights" << std::endl;
  for(int E = 0; E <= nbinsE; E++){

      Avg[E] /= N[E];

    for(int u = 0; u < 1000; u++){

      rms += pow((ratio->GetBinContent(E) - ratio_Syst[u][0]->GetBinContent(E)),2);
      n++;
    }

    rms /= (n); 

    std::cout << " RMS for bin " << ratio->GetXaxis()->GetBinCenter(E) << " is " << 100*sqrt(rms)/ratio->GetBinContent(E) << std::endl; 
    RMS->SetBinContent(E,100*sqrt(rms)/ratio->GetBinContent(E));

    rms = 0;
    n = 0;}

  TLatex *piplus = new TLatex(.15,.95, "MultiSim 0 Weights");
  piplus->SetNDC();
  piplus->SetTextFont(62);
  piplus->SetTextSize(0.05);
  piplus->Draw();  

   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h same");

  c0->Print("../Figs/mult0_ratio.eps");

  TCanvas* c11 = new TCanvas("c11","",700,700);
  c11->SetLeftMargin(.1);
  c11->SetBottomMargin(.1);
  c11->SetTopMargin(.15);
  c11->SetRightMargin(.1);
  c11->cd();

  RMS->SetLineColor(kBlack);
  RMS->SetLineWidth(3);
  RMS->SetLineStyle(1);
  RMS->GetXaxis()->SetTitle("Smeared Neutrino Energy");
  RMS->GetXaxis()->SetLabelSize(0.03);
  RMS->GetXaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetTitle("RMS for MultiSim 0 Weights [%]");
  RMS->GetYaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetLabelSize(0.03);
  RMS->SetStats(0);
 RMS->Draw("h");

  c11->Print("../Figs/mult0_rms.eps");


  TCanvas* c1 = new TCanvas("c1","",700,700);
  c1->SetLeftMargin(.1);
  c1->SetBottomMargin(.1);
  c1->SetTopMargin(.15);
  c1->SetRightMargin(.1);
  c1->cd();

  for(int E = 0; E <= nbinsE; E++) {

    Avg[E] = 0;
    N[E] = 0;
    Avg[E] += NULLVec[E];
    N[E] += 1;
}
  
   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h");

  for(int u = 0; u < 1000; u++){
    for(int E = 0; E <= nbinsE; E++){
 
      ratio_Syst[u][1]->SetBinContent(E,OscVec[u][1][E]);
      Avg[E] += OscVec[u][1][E];
      N[E] += 1;
     
    }

    ratio_Syst[u][1]->SetLineColor(kGray+1);
    ratio_Syst[u][1]->SetLineWidth(1);
    ratio_Syst[u][1]->SetLineStyle(3);
    ratio_Syst[u][1]->Draw("h same");
   
  }

  rms = 0;
  n = 0;

  std::cout <<  "MultiSim 1 Weights" << std::endl;  
  for(int E = 0; E <= nbinsE; E++){

      Avg[E] /= N[E];

    for(int u = 0; u < 1000; u++){

      rms += pow((ratio->GetBinContent(E) - ratio_Syst[u][1]->GetBinContent(E)),2);
      n++;
  }

    rms /= (n); 


    std::cout << " RMS for bin " << ratio->GetXaxis()->GetBinCenter(E) << " is " << 100*sqrt(rms)/ratio->GetBinContent(E) << std::endl; 
    RMS->SetBinContent(E,100*sqrt(rms)/ratio->GetBinContent(E));
    rms = 0;
    n = 0;}

  TLatex *pimin = new TLatex(.15,.95, "MultiSim 1 Weights");
  pimin->SetNDC();
  pimin->SetTextFont(62);
  pimin->SetTextSize(0.05);
  pimin->Draw();  

   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h same");
  c1->Print("../Figs/mult1_ratio.eps");

  TCanvas* c12 = new TCanvas("c12","",700,700);
  c12->SetLeftMargin(.1);
  c12->SetBottomMargin(.1);
  c12->SetTopMargin(.15);
  c12->SetRightMargin(.1);
  c12->cd();

  RMS->SetLineColor(kBlack);
  RMS->SetLineWidth(3);
  RMS->SetLineStyle(1);
  RMS->GetXaxis()->SetTitle("Smeared Neutrino Energy");
  RMS->GetXaxis()->SetLabelSize(0.03);
  RMS->GetXaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetTitle("RMS for MultiSim 1 Weights [%]");
  RMS->GetYaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetLabelSize(0.03);
  RMS->SetStats(0);
 RMS->Draw("h");

  c11->Print("../Figs/mult1_rms.eps");


  TCanvas* c2 = new TCanvas("c2","",700,700);
  c2->SetLeftMargin(.1);
  c2->SetBottomMargin(.1);
  c2->SetTopMargin(.15);
  c2->SetRightMargin(.1);
  c2->cd();

  for(int E = 0; E <= nbinsE; E++) {

    Avg[E] = 0;
    N[E] = 0;
    Avg[E] += NULLVec[E];
    N[E] += 1;
}
  
   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h");

  for(int u = 0; u < 1000; u++){
    for(int E = 0; E <= nbinsE; E++){
 
      ratio_Syst[u][2]->SetBinContent(E,OscVec[u][2][E]);
      Avg[E] += OscVec[u][2][E];
      N[E] += 1;
     
    }

    ratio_Syst[u][2]->SetLineColor(kGray+1);
    ratio_Syst[u][2]->SetLineWidth(1);
    ratio_Syst[u][2]->SetLineStyle(3);
    ratio_Syst[u][2]->Draw("h same");
   
  }

  rms = 0;
  n = 0;
  std::cout << "MultiSim 2 Weights" << std::endl;    

  for(int E = 0; E <= nbinsE; E++){

      Avg[E] /= N[E];

    for(int u = 0; u < 1000; u++){

      rms += pow((ratio->GetBinContent(E) - ratio_Syst[u][2]->GetBinContent(E)),2);
      n++;
  }

    rms /= (n); 


    std::cout << " RMS for bin " << ratio->GetXaxis()->GetBinCenter(E) << " is " << 100*sqrt(rms)/ratio->GetBinContent(E) << std::endl; 
    RMS->SetBinContent(E,100*sqrt(rms)/ratio->GetBinContent(E));
    rms = 0;
    n = 0;}

  TLatex *Kplus = new TLatex(.15,.95,"MultiSim 2 Weights");
  Kplus->SetNDC();
  Kplus->SetTextFont(62);
  Kplus->SetTextSize(0.05);
  Kplus->Draw();  

   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h same");
  c2->Print("../Figs/mult2_ratio.eps");

  TCanvas* c13 = new TCanvas("c13","",700,700);
  c13->SetLeftMargin(.1);
  c13->SetBottomMargin(.1);
  c13->SetTopMargin(.15);
  c13->SetRightMargin(.1);
  c13->cd();

  RMS->SetLineColor(kBlack);
  RMS->SetLineWidth(3);
  RMS->SetLineStyle(1);
  RMS->GetXaxis()->SetTitle("Smeared Neutrino Energy");
  RMS->GetXaxis()->SetLabelSize(0.03);
  RMS->GetXaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetTitle("RMS for MultiSim 2 Weights [%]");
  RMS->GetYaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetLabelSize(0.03);
  RMS->SetStats(0);
 RMS->Draw("h");

  c13->Print("../Figs/mult2_rms.eps");


  TCanvas* c3 = new TCanvas("c3","",700,700);
  c3->SetLeftMargin(.1);
  c3->SetBottomMargin(.1);
  c3->SetTopMargin(.15);
  c3->SetRightMargin(.1);
  c3->cd();


  for(int E = 0; E <= nbinsE; E++) {

    Avg[E] = 0;
    N[E] = 0;
    Avg[E] += NULLVec[E];
    N[E] += 1;
}
  
   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h");



  for(int u = 0; u < 1000; u++){
    for(int E = 0; E <= nbinsE; E++){
 
      ratio_Syst[u][3]->SetBinContent(E,OscVec[u][3][E]);
      Avg[E] += OscVec[u][3][E];
      N[E] += 1;
     
    }

    ratio_Syst[u][3]->SetLineColor(kGray+1);
    ratio_Syst[u][3]->SetLineWidth(1);
    ratio_Syst[u][3]->SetLineStyle(3);
    ratio_Syst[u][3]->Draw("h same");
   
  }

  rms = 0;
  n = 0;

  std::cout << "MultiSim 3 Weights" << std::endl;    
  for(int E = 0; E <= nbinsE; E++){

      Avg[E] /= N[E];

    for(int u = 0; u < 1000; u++){

      rms += pow((ratio->GetBinContent(E) - ratio_Syst[u][3]->GetBinContent(E)),2);
      n++;
  }

    rms /= (n); 
    std::cout << " RMS for bin " << ratio->GetXaxis()->GetBinCenter(E) << " is " << 100*sqrt(rms)/ratio->GetBinContent(E) << std::endl; 
    RMS->SetBinContent(E,100*sqrt(rms)/ratio->GetBinContent(E));
    rms = 0;
    n = 0;}

  TLatex *Kmin = new TLatex(.15,.95,"MultiSim 3 Weights");
  Kmin->SetNDC();
  Kmin->SetTextFont(62);
  Kmin->SetTextSize(0.05);
  Kmin->Draw();  

   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h same");
  c3->Print("../Figs/mult3_ratio.eps");

  TCanvas* c14 = new TCanvas("c14","",700,700);
  c14->SetLeftMargin(.1);
  c14->SetBottomMargin(.1);
  c14->SetTopMargin(.15);
  c14->SetRightMargin(.1);
  c14->cd();

  RMS->SetLineColor(kBlack);
  RMS->SetLineWidth(3);
  RMS->SetLineStyle(1);
  RMS->GetXaxis()->SetTitle("Smeared Neutrino Energy");
  RMS->GetXaxis()->SetLabelSize(0.03);
  RMS->GetXaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetTitle("RMS for MultiSim 3 Weights [%]");
  RMS->GetYaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetLabelSize(0.03);
  RMS->SetStats(0);
 RMS->Draw("h");

  c14->Print("../Figs/mult3_rms.eps");

  TCanvas* c4 = new TCanvas("c4","",700,700);
  c4->SetLeftMargin(.1);
  c4->SetBottomMargin(.1);
  c4->SetTopMargin(.15);
  c4->SetRightMargin(.1);
  c4->cd();

  for(int E = 0; E <= nbinsE; E++) {

    Avg[E] = 0;
    N[E] = 0;
    Avg[E] += NULLVec[E];
    N[E] += 1;
}
  
   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h");

  for(int u = 0; u < 1000; u++){
    for(int E = 0; E <= nbinsE; E++){
 
      ratio_Syst[u][4]->SetBinContent(E,OscVec[u][4][E]);
      Avg[E] += OscVec[u][4][E];
      N[E] += 1;
     
    }

    ratio_Syst[u][4]->SetLineColor(kGray+1);
    ratio_Syst[u][4]->SetLineWidth(1);
    ratio_Syst[u][4]->SetLineStyle(3);
    ratio_Syst[u][4]->Draw("h same");
   
  }

  rms = 0;
  n = 0;

  std::cout << "MultiSim 4 Weights" << std::endl;      
  for(int E = 0; E <= nbinsE; E++){

      Avg[E] /= N[E];

    for(int u = 0; u < 1000; u++){

      rms += pow((ratio->GetBinContent(E) - ratio_Syst[u][4]->GetBinContent(E)),2);
      n++;
  }

    rms /= (n); 

    std::cout << " RMS for bin " << ratio->GetXaxis()->GetBinCenter(E) << " is " << 100*sqrt(rms)/ratio->GetBinContent(E) << std::endl; 
    RMS->SetBinContent(E,100*sqrt(rms)/ratio->GetBinContent(E));

    rms = 0;
    n = 0;}


  TLatex *KL = new TLatex(.15,.95, "MultiSim 4 Weights");
  KL->SetNDC();
  KL->SetTextFont(62);
  KL->SetTextSize(0.05);
  KL->Draw();  

   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h same");
  c4->Print("../Figs/mult4_ratio.eps");

  TCanvas* c15 = new TCanvas("c15","",700,700);
  c15->SetLeftMargin(.1);
  c15->SetBottomMargin(.1);
  c15->SetTopMargin(.15);
  c15->SetRightMargin(.1);
  c15->cd();

  RMS->SetLineColor(kBlack);
  RMS->SetLineWidth(3);
  RMS->SetLineStyle(1);
  RMS->GetXaxis()->SetTitle("Smeared Neutrino Energy");
  RMS->GetXaxis()->SetLabelSize(0.03);
  RMS->GetXaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetTitle("RMS for MultiSim 4 Weights [%]");
  RMS->GetYaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetLabelSize(0.03);
  RMS->SetStats(0);
 RMS->Draw("h");

  c15->Print("../Figs/mult4_rms.eps");

  TCanvas* c5 = new TCanvas("c5","",700,700);
  c5->SetLeftMargin(.1);
  c5->SetBottomMargin(.1);
  c5->SetTopMargin(.15);
  c5->SetRightMargin(.1);
  c5->cd();

  for(int E = 0; E <= nbinsE; E++) {

    Avg[E] = 0;
    N[E] = 0;
    Avg[E] += NULLVec[E];
    N[E] += 1;
}
  
   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h");

  for(int u = 0; u < 1000; u++){
    for(int E = 0; E <= nbinsE; E++){

      ratio_Syst[u][5]->SetBinContent(E,OscVec[u][5][E]);
      Avg[E] += OscVec[u][5][E];
      N[E] += 1;
     
    }

    ratio_Syst[u][5]->SetLineColor(kGray+1);
    ratio_Syst[u][5]->SetLineWidth(1);
    ratio_Syst[u][5]->SetLineStyle(3);
    ratio_Syst[u][5]->Draw("h same");
   
  }

  rms = 0;
  n = 0;


   std::cout << " Unisim Weights" << std::endl;    
  for(int E = 0; E <= nbinsE; E++){

      Avg[E] /= N[E];

    for(int u = 0; u < 1000; u++){

      rms += pow((ratio->GetBinContent(E) - ratio_Syst[u][5]->GetBinContent(E)),2);
      n++;
  }

    rms /= (n); 
 
    std::cout << " RMS for bin " << ratio->GetXaxis()->GetBinCenter(E) << " is " << 100*sqrt(rms)/ratio->GetBinContent(E) << std::endl; 
    RMS->SetBinContent(E,100*sqrt(rms)/ratio->GetBinContent(E));
    rms = 0;
    n = 0;}

  TLatex *uni = new TLatex(.15,.95,"MultiSim 5 Weights");
  uni->SetNDC();
  uni->SetTextFont(62);
  uni->SetTextSize(0.05);
  uni->Draw();  

   ratio->GetYaxis()->SetRangeUser(0,0.15); ratio->Draw("h same");
  c5->Print("../Figs/mult5_ratio.eps");

  TCanvas* c16 = new TCanvas("c16","",700,700);
  c16->SetLeftMargin(.1);
  c16->SetBottomMargin(.1);
  c16->SetTopMargin(.15);
  c16->SetRightMargin(.1);
  c16->cd();

  RMS->SetLineColor(kBlack);
  RMS->SetLineWidth(3);
  RMS->SetLineStyle(1);
  RMS->GetXaxis()->SetTitle("Smeared Neutrino Energy");
  RMS->GetXaxis()->SetLabelSize(0.03);
  RMS->GetXaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetTitle("RMS for MultiSim 5 Weights [%]");
  RMS->GetYaxis()->SetTitleOffset(1.5);
  RMS->GetYaxis()->SetLabelSize(0.03);
  RMS->SetStats(0);
  RMS->Draw("h");

  c16->Print("../Figs/mult5_rms.eps");
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

void Ratios(){
  multiple_detector_fit();
  return;
}

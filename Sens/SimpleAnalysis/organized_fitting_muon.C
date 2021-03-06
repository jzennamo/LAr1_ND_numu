//############################################################################
//Author: Georgia Karagiorgi (georgia@nevis.columbia.edu)
//Author: Corey Adams	(corey.adams@yale.edu)
//Author: Joseph Zennamo (jzennamo@uchicago.edu)
//Nov. 13, 2013
//
//This macro calculates nu_mu disapperance sensitivity for high-Dm2 oscillations
//by performing a chi2 calculation over reconstructed side-by-side nue+numu
//spectra. Up to three detector locations can be simulaneously fit, for
//constraining systematic uncertainties.
//
//Inputs:
//(1) ntuple of neutrinos with enutrue, enurec, baseline, weight, nutype, etc...
//   for each baseline
//(1) ntuple of fully oscillated (numu-->nue) neutrinos with enutrue, erec, baseline
//   weight, nutype, etc... for each baseline
//(2) full fractional systematics covariance matrix, including bin-to-bin and
//   sample-to-sample correlations within a certain baseline, as well as
//   baseline-to-baseline correlations
//
//	By default neutrino mode ntuples are normalized to 6.6e20 POT and antineutrino mode
//	is 10e20 POT.  You can effectively change this by adjusting the scaling of event rates
// 	at baselines
//
//	Original code used histograms for vectors, updated to use std::vectors for several reasons
//	Histograms are a pain to deal with in root, since root does not manage them as expected
//	std::vectors are leaner and meaner, should improve speed too.
//
//
//
//The chi2 is calculated over a 2D grid of Dm2 vs sin22th_mue, and a raster scan
//is used to determine the 90% (one-sided), 3sigma (two-sided) and 5sigma(two-sided)
// 1-dof sensitivity contours.
//############################################################################

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

//for filling the histograms:
//forward declaration of some functions:
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

	//Self explanatory, this is where the ntuples you want to use are located:
	char fileSource[200]="/uboone/app/users/jaz8600/work/lar1/SensitivityCode/Ntuple_Reprocessing/Files/";
	//Note: parts of this routine will write to that same folder, make sure its writable!

	//Set up the bins for the vector:
	//coming soon...

	bool verbose 	=true;	//verbose causes print outs of information that is good for checking
	bool debug 	=true;	//debug prints out info to make sure root is behaving...

	char *  label = "#nu mode, CC events";
    
	TString fileNameRoot = fileSource;
	fileNameRoot += "numu_dis_ND_MicroBooNE";
	double flatSystematicError = 0.15;		// Only used if nearDetStats = false.

	std::string mode = "nu";	//beam mode to run in
	bool use100m = true;		//Include the detector at 100m?
	bool use470m = false;		//Include the detector at 470m?
	bool use700m = true;		//Include the detector at 700m?

	bool shape_only =false;
	
	bool forceRemake = false;

	//Note: there is no infrastructure to handle more than 3 baselines.
	//So, try at your own risk!
	
	Double_t ubooneScale = 1.0;	//Scale the event rates (uniformly across all events), uboone
	Double_t LAr1NDScale = (1.0/3.0);	//Scale the event rates (uniformly across all events), near det
	Double_t LAr1FDScale = (1000.0/345.7);	//Scale the event rates (uniformly across all events), far det

	//Histogram definitions
	Double_t emin = 0.2;		//GeV
	Double_t emax = 3.0;		//GeV
	double bins[20] = {.200, .300, .400, .450, .500, .550, .600, .650, .700, .750, .800, .850, .900, .950, 1.000, 1.250, 1.500, 2.000, 2.500, 3.000};
	//	const Int_t nbinsE = 10;    //number of energy bins in each baseline, for each distribution (fosc, nue, or numu)
    
	// Options are etrue, eccqe, ecalo1, ecalo2;
	//	std::string energyType = "etrue";

	
	//How many points in the final graph do you want?  (symmetric in x and y)
	//More points means longer run time but better graphs
	//Note: most of the run time is in looping over ntuples, which only takes awhile
	//on the very first pass!  (subsequent runs are much faster)
	Int_t npoints = 500;
    
	//grid boundaries
	Double_t dm2min = 0.01;                       //eV**2
	Double_t dm2max = 100.;                       //eV**2
	Double_t sin22thmin = 0.0001;
	Double_t sin22thmax = 1.0;
	
	bool useNearDetStats = true; 			// Only matters if the covariance matrix vector is empty.
	std::vector<std::string>  cov_max_name;
	//cov_max_name.push_back(mode+"_covariance_matrix/output_fraccovmx_xsec.out");
    

    
//---------------End of section you should alter without care
//

	std::vector<int> baselines;
	std::vector<double> scales;
	std::vector<std::string> names;
	std::vector<double> volume;
  
  	//initialize all the baselines we'll use:
  	if (use100m) baselines.push_back(100);
  	if (use470m) baselines.push_back(470);
  	if (use700m) baselines.push_back(700);
  	Int_t nL = baselines.size();
	std::cout << " Number of baselines " << nL << std::endl;
  	if (nL == 0){
  		std::cout << "\nNo baselines selected!  What did you expect? Exiting.\n" << std::endl;
  		return 2;
  	}

  	//and they're corresponding scaling factors:
  	if (use100m) scales.push_back(LAr1NDScale);
  	if (use470m) scales.push_back(ubooneScale);
  	if (use700m) scales.push_back(LAr1FDScale);
  	//some detector names just for fun:
  	if (use100m) names.push_back("LAr1-ND");
  	if (use470m) names.push_back("MicroBooNE");
  	if (use700m) names.push_back("LAr1-FD");
	//fiducial volume FROM MC, also just for fun: (scale to get target volume)
  	if (use100m) volume.push_back(40);
  	if (use470m) volume.push_back(61.4);
  	if (use700m) volume.push_back(345.7);
  	

	double NULLVec = 0;
	std::vector <double>  OscVec;

	OscVec.resize(npoints+1);

	int nbinsE = 0;
	
	std::string temp_name = "../MatrixFiles/combined_ntuple_600m_onaxis_nu_processed_numu_Joseph_Smeared.root";
	
	TFile temp_file(temp_name.c_str());
	TH1D *NULL_700;
	NULL_700 = (TH1D*)(temp_file.Get("NumuCC"));
	nbinsE = NULL_700->GetNbinsX();
	for(int i = 1; i <= nbinsE; i++){
	  NULLVec += (NULL_700->GetBinContent(i));
	}
	
	for(int dm = 0; dm <= npoints; dm++){
            TH1D *OSC_700;
	    std::string dmpoint = std::to_string(dm);
	    std::string name = "Osc_";
            name = name+dmpoint;
            OSC_700 = (TH1D*)(temp_file.Get(name.c_str()));;
	    double m = 0;
            for(int i = 1; i <= nbinsE; i++){
              m += (OSC_700->GetBinContent(i));
	    }
	    OscVec[dm] = m;
	    delete OSC_700;
	}
	
	delete NULL_700;
	temp_file.Close();
	
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
	
	//Fitting variables
	Double_t chisq;
	Double_t deltachisq90=1.64, deltachisq3s=9.00, deltachisq5s=25.00, deltachisq99=2.58;
	
	//sensitivity contours
	
	TGraph *sens90; Double_t x90[npoints+1]; Double_t y90[npoints+1];
	TGraph *sens3s; Double_t x3s[npoints+1]; Double_t y3s[npoints+1];
	TGraph *sens5s; Double_t x5s[npoints+1]; Double_t y5s[npoints+1];
	

	//Diagonal Error matrix, add in Corey's stuff if we get the correlations

	double cov = 1/pow(0.04*NULLVec,2);
	std::cout << "...inverted..."  << std::endl;

	//Time to compute the chi2

	double DMpoint;
	double sin22th;

	double cvi, cvi_total;
	double predictioni, pred_total;
	for(int dm2 = 0; dm2 <= npoints; dm2++){

	  DMpoint = pow(10., (TMath::Log10(dm2min)+ (dm2 * (1./npoints))*TMath::Log10(dm2max/dm2min)) );
	  
	  for(int sint = 0; sint <= npoints; sint++){
	    
	    sin22th = pow(10., (TMath::Log10(sin22thmin)+ (sint * (1./npoints))*TMath::Log10(sin22thmax/sin22thmin)));
	   
	    
	    chisq = 0.0;
	      
	      cvi = predictioni = 0;
	      cvi = NULLVec;
	      predictioni = NULLVec-(sin22th)*OscVec[dm2];  		
	      chisq += (predictioni-cvi)*(predictioni-cvi)* (cov);
	      

	      if (chisq <= deltachisq5s){y5s[dm2] = DMpoint; x5s[dm2] = sin22th;}
	      if (chisq <= deltachisq3s){y3s[dm2] = DMpoint; x3s[dm2] = sin22th;}
	      if (chisq <= deltachisq90){y90[dm2] = DMpoint; x90[dm2] = sin22th;}	      
	      
	  }
	}
	
	//Plot Results:
	sens90 = new TGraph(npoints+1,x90,y90); sens90->SetLineColor(1); sens90->SetLineWidth(2);
	sens3s = new TGraph(npoints+1,x3s,y3s); sens3s->SetLineColor(9); sens3s->SetLineWidth(2);
	sens5s = new TGraph(npoints+1,x5s,y5s); sens5s->SetLineColor(9); sens5s->SetLineStyle(2); sens5s->SetLineWidth(1);
	
	//write the results to file:
	//======================================================
	printf("\nDrawing sensitivity curves...\n");
	
	
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
	hr1->SetTitle(";sin^{2}2#theta_{#mu#mu};#Deltam_{41}^{2} [eV^{2}]");
	hr1->GetXaxis()->SetTitleOffset(1.2);
	hr1->GetYaxis()->SetTitleOffset(1.2);
	hr1->GetXaxis()->SetTitleFont(62);
        hr1->GetYaxis()->SetTitleFont(62);
	hr1->GetYaxis()->CenterTitle();
	hr1->GetXaxis()->CenterTitle();
	hr1->GetXaxis()->SetTitleSize(0.05);
        hr1->GetXaxis()->SetLabelSize(0.04);
	hr1->GetXaxis()->SetLabelOffset(-0.003);
	hr1->GetYaxis()->SetTitleSize(0.05);
        hr1->GetYaxis()->SetLabelSize(0.04);
	hr1->SetStats(kFALSE);
	hr1->Draw();

	//	if(mode == "nu")        gROOT->ProcessLine(".x ./minisciboon_plot_numu.c(c3)");
	//if(mode == "nubar")     gROOT->ProcessLine(".x ./minisciboon_plot_numubar.c(c3)");
	
	//======================================================

       //	TLatex *tex_Detector = new TLatex(.2,.9,"2 #times LAr1-ND (200m), MicroBooNE (470m),");	
       // TLatex *tex_FDetector = new TLatex(.34,.86,"ICARUS (700m)");
       //TLatex *tex_Detector = new TLatex(.2,.9,"LAr1-ND (200m), T600 (700m)");
       //TLatex *tex_FDetector = new TLatex(.34,.86,"LAr1-FD (700m)");	  

	TLatex *tex_mode = new TLatex(.34,.97,"#nu mode, CC Events");
	tex_mode->SetNDC();
	tex_mode->SetTextFont(62);
	tex_mode->SetTextSize(0.04);
	tex_mode->Draw();

	TLatex *tex_Detector = new TLatex(.2,.9,"T600 (600m), 1 bin, 4% Uncert.");
	tex_Detector->SetNDC();
	tex_Detector->SetTextFont(62);
	tex_Detector->SetTextSize(0.035);
	tex_Detector->Draw();

       	TLegend* legt=new TLegend(0.2,0.2,0.6,0.45);
	legt->SetFillStyle(0);
	legt->SetFillColor(0);
	legt->SetBorderSize(0);
	legt->SetTextSize(0.04);

	/*
       	char name[200];
	bool added = false;
	for (int j=0; j<nL; j++) {
	  if (baselines[j] == 100){
	    sprintf(name + strlen(name), "LAr1-ND (%im)", baselines[j]);
	    added = true;
		}
	  else if (baselines[j] == 470){
	    
	    if (added) sprintf(name + strlen(name), ", MicroBooNE (%im)", baselines[j]);
	    else sprintf(name + strlen(name), "MicroBooNE (%im)", baselines[j]);
	    added = true;
	  }
	  else if (baselines[j] == 700){
	    if (added) sprintf(name + strlen(name), ", LAr1-FD (%im)", baselines[j]);
	    else sprintf(name + strlen(name), "LAr1-FD (%im)", baselines[j]);
	    added = true;
	  }
	}
	
	std::cout << "Name is " << name << std::endl;
	add_plot_label(name, 0.74, 0.85, 0.04, 1, 62);
	*/

	//	gROOT->ProcessLine(".x ./minisciboon_plot_numu.c(c3)");

	sens90->Draw("l same");
	legt->AddEntry(sens90,"90\% CL","l");
      	sens3s->Draw("l same");
	legt->AddEntry(sens3s,"3#sigma CL","l");
	sens5s->Draw("l same");
	legt->AddEntry(sens5s,"5#sigma CL","l");
	TLine *gdummy3 = new TLine();
	gdummy3->SetLineColor(kRed-3);
        gdummy3->SetLineStyle(2);
        gdummy3->SetLineWidth(2);
	legt->AddEntry(gdummy3,"MiniBooNE + SciBooNE 90% CL","l");
	legt->Draw();
	
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

void organized_fitting_muon(){
  multiple_detector_fit();
  return;
}


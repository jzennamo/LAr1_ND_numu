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
  	

	std::vector< std::vector <float> >   NULLVec;
	std::vector< std::vector <float> >   LEVec;
	std::vector< std::vector < std::vector <float> > >   OscVec;
	
	NULLVec.resize(nL);
	LEVec.resize(nL);
	OscVec.resize(nL);

	for( int l = 0; l < nL; l++){
	  OscVec[l].resize(npoints+1);
	}
	
	int counter = 0;
	int nbinsE = 0;

	std::cout << "Beginning : ... " << std::endl;

	if (use100m){
	  std::string temp_name ="../MatrixFiles/combined_ntuple_100m_nu_processed_numu_Joseph_Smeared.root";

	  TFile temp_file(temp_name.c_str());
	  TH1D *NULL_100;
	  NULL_100 = (TH1D*)(temp_file.Get("NumuCC"));      
	
	  nbinsE = NULL_100->GetNbinsX();

	  for(int i = 1; i <= nbinsE; i++){
	    NULLVec[counter].push_back(NULL_100->GetBinContent(i));
	  }

	  for(int dm = 0; dm <= npoints; dm++){
	    TH1D *OSC_100;
	    std::string dmpoint = std::to_string(dm);
	    std::string name = "Osc_";
	    name = name+dmpoint;
	    OSC_100 = (TH1D*)(temp_file.Get(name.c_str()));

	    for(int i = 1; i <= nbinsE; i++){	      
	      OscVec[counter][dm].push_back(OSC_100->GetBinContent(i));
	    }
	    delete OSC_100;
	  }

	  delete NULL_100;
	  temp_file.Close();
	  counter++;
	}

	if (use700m) {
	  std::string temp_name = "../MatrixFiles/combined_ntuple_600m_onaxis_nu_processed_numu_Joseph_Smeared.root";
        
	  TFile temp_file(temp_name.c_str());
          TH1D *NULL_700;
          NULL_700 = (TH1D*)(temp_file.Get("NumuCC"));

	  for(int i = 1; i <= nbinsE; i++){
            NULLVec[counter].push_back(NULL_700->GetBinContent(i));
          }

          for(int dm = 0; dm <= npoints; dm++){
            TH1D *OSC_700;
	    std::string dmpoint = std::to_string(dm);
	    std::string name = "Osc_";
            name = name+dmpoint;
            OSC_700 = (TH1D*)(temp_file.Get(name.c_str()));;
            for(int i = 1; i <= nbinsE; i++){
              OscVec[counter][dm].push_back(OSC_700->GetBinContent(i));
	    }
	    delete OSC_700;
          }

          delete NULL_700;
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
	
	//Fitting variables
	Double_t chisq;
	Double_t deltachisq90=1.64, deltachisq3s=9.00, deltachisq5s=25.00, deltachisq99=2.58;
	
	//output ntuple with chi2 values and sensitivity contour	
	TH1D* ND_null_lmls = new TH1D("ND_null_lmls", "", nbinsE, bins);
	TH1D* ND_null_hmls = new TH1D("ND_null_hmls", "", nbinsE, bins);
        TH1D* ND_null_lmhs = new TH1D("ND_null_lmhs", "", nbinsE, bins);
        TH1D* ND_null_hmhs = new TH1D("ND_null_hmhs", "", nbinsE, bins);

        TH1D* Ratio1 = new TH1D("Ratio1", "", nbinsE, bins);
        TH1D* Ratio2 = new TH1D("Ratio2", "", nbinsE, bins);
        TH1D* Ratio3 = new TH1D("Ratio3", "", nbinsE, bins);
        TH1D* Ratio4 = new TH1D("Ratio4", "", nbinsE, bins);

	TH1D* ND_lowm_lows = new TH1D("ND_lowm_lows", "", nbinsE, bins);
        TH1D* ND_highm_lows = new TH1D("ND_highm_lows", "", nbinsE, bins);
        TH1D* ND_lowm_highs = new TH1D("ND_lowm_highs", "", nbinsE, bins);
        TH1D* ND_highm_highs = new TH1D("ND_highm_highs", "", nbinsE, bins);

	TH1D* FD_null_lmls = new TH1D("FD_null_lmls", "", nbinsE, bins);
	TH1D* FD_null_hmls = new TH1D("FD_null_hmls", "", nbinsE, bins);
	TH1D* FD_null_lmhs = new TH1D("FD_null_lmhs", "", nbinsE, bins);
	TH1D* FD_null_hmhs = new TH1D("FD_null_hmhs", "", nbinsE, bins);

	TH1D* FD_lowm_lows = new TH1D("FD_lowm_lows", "", nbinsE, bins);
	TH1D* FD_highm_lows = new TH1D("FD_highm_lows", "", nbinsE, bins);
	TH1D* FD_lowm_highs = new TH1D("FD_lowm_highs", "", nbinsE, bins);
	TH1D* FD_highm_highs = new TH1D("FD_highm_highs", "", nbinsE, bins);

	
	
	//sensitivity contours
	
	TGraph *sens90; Double_t x90[npoints+1]; Double_t y90[npoints+1];
	TGraph *sens3s; Double_t x3s[npoints+1]; Double_t y3s[npoints+1];
	TGraph *sens5s; Double_t x5s[npoints+1]; Double_t y5s[npoints+1];
	

	//Diagonal Error matrix, add in Corey's stuff if we get the correlations
	TMatrix *M = new TMatrix(nbinsE*nL,nbinsE*nL);
	std::vector< float > CVvec;
	std::vector< std::vector< float > > Predvec;
	Predvec.resize(npoints+1);


	double high_DM2; 
	double low_DM2;
	double high_Sin2theta;
	double low_Sin2theta;
	std::cout << "Filling matrix..."  << std::endl;
	int Errbin = 0;
	for(int Lbin = 0; Lbin < nL; Lbin++){
	  for(int Ebin = 0; Ebin < nbinsE; Ebin++){
	    
	    //hhhmmmm...This is where I disagree with Davio
	    if(Lbin == 0){ (*M) (Errbin, Errbin) = NULLVec[Lbin][Ebin] + pow(flatSystematicError*NULLVec[Lbin][Ebin],2);}
	    else{ (*M) (Errbin, Errbin) = NULLVec[Lbin][Ebin] + pow((sqrt(NULLVec[0][Ebin])/NULLVec[0][Ebin])*NULLVec[Lbin][Ebin],2) + pow(0.04*NULLVec[Lbin][Ebin],2);}


	    CVvec.push_back(NULLVec[Lbin][Ebin]);	    

	    Errbin++;	    
	    
	  }//Energy bins
	}//Detectors

	std::cout << "...done with matrix..."  << std::endl;
	std::cout << "...filling prediction vector..."  << std::endl;

	for(int dm = 0; dm <= npoints; dm++){
	  for(int Lbin = 0; Lbin < nL; Lbin++){
	    for(int Ebin = 0; Ebin < nbinsE; Ebin++){
	      
	      Predvec[dm].push_back(OscVec[Lbin][dm][Ebin]);

	    }
	  }
	}

	std::cout << "...inverting matrix..."  << std::endl;

	//Now create the inverse covariance matrix.  Root can invert for us:
	//inverse cov matrix, M^{-1}, used in chi2 calculation
	TMatrix *cov = new TMatrix(nbinsE*nL,nbinsE*nL);
	cov = &(M->Invert());
	std::cout << "...inverted..."  << std::endl;

	//Time to compute the chi2

	double DMpoint;
	double sin22th;

	double cvi, cvi_total;
	double predictioni, pred_total;
	for(int dm2 = 0; dm2 <= npoints; dm2++){

	  DMpoint = pow(10., (TMath::Log10(dm2min)+ (dm2 * (1./npoints))*TMath::Log10(dm2max/dm2min)) );
	  
	  if(dm2 == 255) low_DM2 =  DMpoint;
	  if(dm2 ==  205) high_DM2 =  DMpoint;

		  
	  for(int sint = 0; sint <= npoints; sint++){
	    
	    sin22th = pow(10., (TMath::Log10(sin22thmin)+ (sint * (1./npoints))*TMath::Log10(sin22thmax/sin22thmin)));
	    
	   if(sint == 255 && dm2 == 255) low_Sin2theta = sin22th;
	   if(sint == 378 && dm2 ==  205) high_Sin2theta = sin22th;
	    
	    
	    chisq = 0.0;
	      
            double Scaling[nbinsE*nL];
            if(shape_only){
              //Scale all measurements to the ND NULL levels                           
              for(int Edecbin = 0; Edecbin < nbinsE; Edecbin++){

                Scaling[Edecbin]  = CVvec[Edecbin]-(sin22th)*Predvec[dm2][Edecbin];
                if((CVvec[Edecbin]) != 0){Scaling[Edecbin] /= CVvec[Edecbin];}
                if((CVvec[Edecbin]) == 0){ std::cout << "Scaling failed....CV==0!" << std::endl; break;}

              }
              // Copy results from first baseline to the rest of these baselines      
	      int tick = 0;
              for(int Edecbin = nbinsE; Edecbin < nbinsE*nL; Edecbin++){
		if(tick == nbinsE) tick = 0;
		Scaling[Edecbin] = Scaling[tick];
		tick++;
	      }
            }//Shape only
	    if(!shape_only){
              //We don't want to scale anything!
	      for(int Edecbin = 0; Edecbin < nbinsE*nL; Edecbin++){
		Scaling[Edecbin] = 1;
	      }

            }

	    for(int Edecbin = 0; Edecbin < nbinsE*nL; Edecbin++){
	      
	      cvi = predictioni = 0;
	      cvi = Scaling[Edecbin]*CVvec[Edecbin];
	      predictioni = CVvec[Edecbin]-(sin22th)*Predvec[dm2][Edecbin];  		
	      chisq += (predictioni-cvi)*(predictioni-cvi)* (*cov)(Edecbin,Edecbin);
	      
	      //To do...swich the scale factor. We measure the TEST and model the NULL, means I have to build another set of NULL_hists which are scaled to the prediction.

	      if(Edecbin < nbinsE && dm2 == 255 && sint == 255){    
		ND_lowm_lows->SetBinContent(Edecbin+1,predictioni/(Ratio1->GetXaxis()->GetBinWidth(Edecbin+1))); 
		ND_null_lmls->SetBinContent(Edecbin+1,cvi/(Ratio1->GetXaxis()->GetBinWidth(Edecbin+1)));
	      }// ND_null_lmls->SetBinError(Edecbin+1,1/sqrt((*cov)(Edecbin+1,Edecbin+1)));}

	      if(Edecbin < nbinsE && dm2 ==  205 && sint == 255){   
		ND_highm_lows->SetBinContent(Edecbin+1,predictioni/(Ratio1->GetXaxis()->GetBinWidth(Edecbin+1)));
		ND_null_hmls->SetBinContent(Edecbin+1,cvi/(Ratio1->GetXaxis()->GetBinWidth(Edecbin+1)));
	      }//ND_null_hmls->SetBinError(Edecbin+1,1/sqrt((*cov)(Edecbin+1,Edecbin+1)));}

	      if(Edecbin < nbinsE && dm2 == 206 && sint == 378){  
		ND_lowm_highs->SetBinContent(Edecbin+1,predictioni/(Ratio1->GetXaxis()->GetBinWidth(Edecbin+1)));
		ND_null_lmhs->SetBinContent(Edecbin+1,cvi/(Ratio1->GetXaxis()->GetBinWidth(Edecbin+1)));
	      }//ND_null_lmhs->SetBinError(Edecbin+1,1/sqrt((*cov)(Edecbin+1,Edecbin+1)));}

	      if(Edecbin < nbinsE && dm2 ==  205 && sint == 378){
		ND_highm_highs->SetBinContent(Edecbin+1,predictioni/(Ratio1->GetXaxis()->GetBinWidth(Edecbin+1)));
		ND_null_hmhs->SetBinContent(Edecbin+1,cvi/(Ratio1->GetXaxis()->GetBinWidth(Edecbin+1)));
	      }//ND_null_hmhs->SetBinError(Edecbin+1,1/sqrt((*cov)(Edecbin+1,Edecbin+1)));}

	      if(Edecbin >= nbinsE && Edecbin < 2*nbinsE && dm2 == 255 && sint == 255){    
		FD_lowm_lows->SetBinContent(Edecbin-nbinsE+1,predictioni/(Ratio1->GetXaxis()->GetBinWidth(Edecbin-nbinsE+1)));
		FD_null_lmls->SetBinContent(Edecbin-nbinsE+1,cvi/(Ratio1->GetXaxis()->GetBinWidth(Edecbin-nbinsE+1)));
	      }//FD_null_lmls->SetBinError(Edecbin-nbinsE+1+1,1/sqrt((*cov)(Edecbin+1,Edecbin+1)));}

	      if(Edecbin >= nbinsE && Edecbin < 2*nbinsE && dm2 ==  205 && sint == 255){ 
		FD_highm_lows->SetBinContent(Edecbin-nbinsE+1,predictioni/(Ratio1->GetXaxis()->GetBinWidth(Edecbin-nbinsE+1)));
		FD_null_hmls->SetBinContent(Edecbin-nbinsE+1,cvi/(Ratio1->GetXaxis()->GetBinWidth(Edecbin-nbinsE+1)));
	      }//FD_null_hmls->SetBinError(Edecbin-nbinsE+1+1,1/sqrt((*cov)(Edecbin+1,Edecbin+1)));}

	      if(Edecbin >= nbinsE && Edecbin < 2*nbinsE && dm2 == 255 && sint == 378){ 
		FD_lowm_highs->SetBinContent(Edecbin-nbinsE+1,predictioni/(Ratio1->GetXaxis()->GetBinWidth(Edecbin-nbinsE+1)));
		FD_null_lmhs->SetBinContent(Edecbin-nbinsE+1,cvi/(Ratio1->GetXaxis()->GetBinWidth(Edecbin-nbinsE+1)));
	      }//FD_null_lmhs->SetBinError(Edecbin-nbinsE+1+1,1/sqrt((*cov)(Edecbin+1,Edecbin+1)));}

	      if(Edecbin >= nbinsE && Edecbin < 2*nbinsE && dm2 ==  205 && sint == 378){  
		FD_highm_highs->SetBinContent(Edecbin-nbinsE+1,predictioni/(Ratio1->GetXaxis()->GetBinWidth(Edecbin-nbinsE+1)));
		FD_null_hmhs->SetBinContent(Edecbin-nbinsE+1,cvi/(Ratio1->GetXaxis()->GetBinWidth(Edecbin-nbinsE+1)));
	      }//FD_null_hmhs->SetBinError(Edecbin-nbinsE+1+1,1/sqrt((*cov)(Edecbin+1,Edecbin+1)));}

		


	    }//Sum over all detectors and energies to determine the chi2 for that osc. point	      
	      
	      //    	      chi2->Fill(sin22th,DMpoint,chisq);

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
	hr1->SetTitle(";sin^{2}2#theta#lower[0.4]{#mu#kern[-0.3]{#mu}};#Deltam_{41}^{2} [eV^{2}]");
	hr1->GetXaxis()->SetTitleOffset(1.2);
	hr1->GetYaxis()->SetTitleOffset(1.2);
	hr1->GetXaxis()->SetTitleFont(62);
        hr1->GetYaxis()->SetTitleFont(62);
	hr1->GetYaxis()->CenterTitle();
	hr1->GetXaxis()->CenterTitle();
	hr1->GetXaxis()->SetTitleSize(0.06);
        hr1->GetXaxis()->SetLabelSize(0.055);
	hr1->GetXaxis()->SetLabelOffset(-0.003);
	hr1->GetYaxis()->SetTitleSize(0.06);
        hr1->GetYaxis()->SetLabelSize(0.055);
	hr1->SetStats(kFALSE);
	hr1->Draw();

	//	if(mode == "nu")        gROOT->ProcessLine(".x ./minisciboon_plot_numu.c(c3)");
	//if(mode == "nubar")     gROOT->ProcessLine(".x ./minisciboon_plot_numubar.c(c3)");
	
	//======================================================

       //	TLatex *tex_Detector = new TLatex(.2,.9,"2 #times LAr1-ND (100m), MicroBooNE (470m),");	
       // TLatex *tex_FDetector = new TLatex(.34,.86,"ICARUS (700m)");
       //TLatex *tex_Detector = new TLatex(.2,.9,"LAr1-ND (100m), T600 (700m)");
       //TLatex *tex_FDetector = new TLatex(.34,.86,"LAr1-FD (700m)");	  

	TLatex *tex_mode = new TLatex(.34,.97,"#nu mode, CC Events");
	tex_mode->SetNDC();
	tex_mode->SetTextFont(62);
	tex_mode->SetTextSize(0.04);
	tex_mode->Draw();

	TLatex *tex_Detector = new TLatex(.2,.9,"#splitline{LAr1-ND (100m)}{and T600 (600m)}");
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
	
	//	c3->Print("Sensitivity.pdf");

       	double xS[2];
	double yS[2];
	xS[0] = 0;
	xS[1] = 3;
	yS[0] = 1;
	yS[1] = 1;

	TGraph* straight=new TGraph(2,xS,yS);
	straight->SetLineWidth(4.5);
	straight->SetLineColor(kBlack);
	straight->SetLineStyle(2);

	TCanvas* c4 = new TCanvas("c4","Near Detector - Small #Delta m^{2}, Small sin^{2}(2#theta)",700,700);
	c4->cd();

        TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
        pad1 -> SetBottomMargin(0);
        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.25);
        pad2 -> SetTopMargin(0);
	pad2 -> SetBottomMargin(0.35);

        pad1 -> Draw();
        pad2 -> Draw();	
        pad1 -> cd();
	//	TGaxis::SetMaxDigits(3);   
        ND_null_lmls->SetTitle(";Smeared Neutrino Energy [GeV];Events / Bin");
        ND_null_lmls->GetXaxis()->SetTitleOffset(1.2);
        ND_null_lmls->GetYaxis()->SetTitleOffset(1.2);
        ND_null_lmls->GetXaxis()->SetTitleFont(62);
        ND_null_lmls->GetYaxis()->SetTitleFont(62);
        ND_null_lmls->GetYaxis()->CenterTitle();
        ND_null_lmls->GetXaxis()->CenterTitle();
        ND_null_lmls->GetXaxis()->SetTitleSize(0.06);
        ND_null_lmls->GetXaxis()->SetLabelSize(0.055);
        ND_null_lmls->GetXaxis()->SetLabelOffset(-0.003);
        ND_null_lmls->GetYaxis()->SetTitleSize(0.06);
        ND_null_lmls->GetYaxis()->SetLabelSize(0.055);
        ND_null_lmls->SetStats(kFALSE);

	ND_null_lmls->SetLineWidth(2);
        ND_null_lmls->SetLineColor(kBlack);
        ND_null_lmls->SetMinimum(0.5);
	//	ND_null_lmls->Draw("h");

	ND_lowm_lows->SetMarkerColor(kRed);		
	ND_lowm_lows->SetMarkerStyle(22);
	ND_lowm_lows->SetMarkerSize(1.5);
	//	ND_lowm_lows->Draw("p1 same");


        TLatex *tex_pre = new TLatex(.6,.92,"PRELIMINARY");
	tex_pre->SetNDC();
        tex_pre->SetTextFont(62);
        tex_pre->SetTextColor(kRed-3);
        tex_pre->SetTextSize(0.05);
        tex_pre->Draw();
       
	TLatex *tex_ND = new TLatex(.6,.85,"LAr1-ND (100m)");
        tex_ND->SetNDC();
        tex_ND->SetTextFont(62);
        tex_ND->SetTextSize(0.05);
	tex_ND->Draw();

	TLatex *tex_NDPOT = new TLatex(.6,.77,"P.O.T. = 6.6 #times 10^{20}");
	tex_NDPOT->SetNDC();
	tex_NDPOT->SetTextFont(62);
	tex_NDPOT->SetTextSize(0.05);
	tex_NDPOT->Draw();

        TLatex *tex_stat = new TLatex(.6,.68,"Stat. Uncert. Only");
        tex_stat->SetNDC();
        tex_stat->SetTextFont(62);
        tex_stat->SetTextSize(0.05);
        tex_stat->Draw();

	char buff_pt_dm_lmls[100];
	sprintf(buff_pt_dm_lmls,"#Delta m_{41}^{2} = %4.2f eV^{2}", low_DM2);
	TLatex *tex_dm_lmls = new TLatex(.6,.6,buff_pt_dm_lmls);
	tex_dm_lmls->SetNDC();
	tex_dm_lmls->SetTextFont(62);
	tex_dm_lmls->SetTextSize(0.05);
	tex_dm_lmls->Draw();

        char buff_pt_sin_lmls[100];
        sprintf(buff_pt_sin_lmls,"sin^{2}(2#theta#lower[0.4]{#mu#kern[-0.3]{#mu}}) = %.2f", low_Sin2theta);
        TLatex *tex_sin_lmls = new TLatex(.6,.5,buff_pt_sin_lmls);
        tex_sin_lmls->SetNDC();
        tex_sin_lmls->SetTextFont(62);
        tex_sin_lmls->SetTextSize(0.05);
        tex_sin_lmls->Draw();
    
	TLegend* leg=new TLegend(0.62,0.2,0.82,0.45);
        leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextFont(62);
	leg->SetTextSize(0.04);
	leg->SetTextAlign(12);
        leg->AddEntry(ND_null_lmls,"Unoscillated","l");
	leg->AddEntry(ND_lowm_lows,"Oscillated","p");
        TMarker *gdummy = new TMarker();
        gdummy->SetMarkerColor(kBlue-3);
        gdummy->SetMarkerStyle(20);
        gdummy->SetMarkerSize(1);
        leg->AddEntry(gdummy,"#splitline{Ratio of}{Osc. to Unosc.}","p");
        leg->Draw();	

        pad2 -> cd();

	Ratio1->Divide(ND_lowm_lows,ND_null_lmls,1,1,"B");
        Ratio1->SetTitle(";Smeared Neutrino Energy [GeV];#splitline{Ratio}{Osc.-to-Unosc.} ");
        Ratio1->GetXaxis()->SetTitleOffset(0.95);
        Ratio1->GetYaxis()->SetTitleOffset(0.4);
        Ratio1->GetXaxis()->SetTitleFont(62);
	Ratio1->GetYaxis()->SetTitleFont(62);
	Ratio1->GetYaxis()->CenterTitle();
        Ratio1->GetXaxis()->CenterTitle();
        Ratio1->GetXaxis()->SetTitleSize(0.18);
	Ratio1->GetXaxis()->SetLabelSize(0.16);
	Ratio1->GetXaxis()->SetLabelOffset(0.008);
	Ratio1->GetYaxis()->SetTitleSize(0.1);
	Ratio1->GetYaxis()->SetTitleOffset(0.59);
        Ratio1->GetYaxis()->SetLabelSize(0.1);
	Ratio1->GetYaxis()->SetNdivisions(505);
        Ratio1->SetStats(kFALSE);
	Ratio1->SetMinimum(0.95);
	Ratio1->SetMaximum(1.05);
        Ratio1->SetMarkerColor(kBlue-3);
	Ratio1->SetLineColor(kBlue-3);
	Ratio1->SetMarkerStyle(20);
        Ratio1->SetMarkerSize(1);
	Ratio1->Draw("p");
	straight->Draw("l same");
	//	low_Sin2theta, high_DM2;

	
        TCanvas* c5 = new TCanvas("c5","Near Detector - High #Delta m^{2}, High sin^{2}(2#theta)",700,700);
        c5->cd();

        TPad *pad3 = new TPad("pad3","pad3",0,0.25,1,1);
        pad3 -> SetBottomMargin(0);
	TPad *pad4 = new TPad("pad4","pad4",0,0,1,0.25);
        pad4 -> SetTopMargin(0);
        pad4 -> SetBottomMargin(0.35);

        pad3 -> Draw();
	pad4 -> Draw();
        pad3 -> cd();

        ND_null_hmhs->SetTitle(";Smeared Neutrino Energy [GeV];Events / Bin");
        ND_null_hmhs->GetXaxis()->SetTitleOffset(1.2);
        ND_null_hmhs->GetYaxis()->SetTitleOffset(1.2);
        ND_null_hmhs->GetXaxis()->SetTitleFont(62);
        ND_null_hmhs->GetYaxis()->SetTitleFont(62);
        ND_null_hmhs->GetYaxis()->CenterTitle();
        ND_null_hmhs->GetXaxis()->CenterTitle();
        ND_null_hmhs->GetXaxis()->SetTitleSize(0.06);
        ND_null_hmhs->GetXaxis()->SetLabelSize(0.055);
        ND_null_hmhs->GetXaxis()->SetLabelOffset(-0.003);
        ND_null_hmhs->GetYaxis()->SetTitleSize(0.06);
        ND_null_hmhs->GetYaxis()->SetLabelSize(0.055);
	//	ND_null_hmhs->GetYaxis()->SetMaxDigits(3);
        ND_null_hmhs->SetStats(kFALSE);
        ND_null_hmhs->SetLineWidth(2);
        ND_null_hmhs->SetLineColor(kBlack);
        ND_null_hmhs->SetMinimum(0.5);
	ND_null_hmhs->GetYaxis()->SetRangeUser(10000,2000000);
	std::cout << "Total Event Count : " << ND_null_hmhs->Integral() << std::endl; 
        ND_null_hmhs->Draw("h");

        ND_highm_highs->SetMarkerColor(kRed);
	ND_highm_highs->SetMarkerStyle(22);
        ND_highm_highs->SetMarkerSize(1.5);
	ND_highm_highs->Draw("p1 same");

        tex_stat->Draw();
        tex_pre->Draw();
        tex_ND->Draw();
	tex_NDPOT->Draw();

        char buff_pt_dm_hmhs[100];
	sprintf(buff_pt_dm_hmhs,"#Delta m_{41}^{2} = %4.2f eV^{2}", high_DM2);
        TLatex *tex_dm_hmhs = new TLatex(.6,.60,buff_pt_dm_hmhs);
	tex_dm_hmhs->SetNDC();
        tex_dm_hmhs->SetTextFont(62);
        tex_dm_hmhs->SetTextSize(0.05);
        tex_dm_hmhs->Draw();

        char buff_pt_sin_hmhs[100];
        sprintf(buff_pt_sin_hmhs,"sin^{2}(2#theta#lower[0.4]{#mu#kern[-0.3]{#mu}}) = %.2f", high_Sin2theta);
        TLatex *tex_sin_hmhs = new TLatex(.6,.50,buff_pt_sin_hmhs);
        tex_sin_hmhs->SetNDC();
	tex_sin_hmhs->SetTextFont(62);
        tex_sin_hmhs->SetTextSize(0.05);
	tex_sin_hmhs->Draw();
        leg->Draw();


        pad4 -> cd();

        Ratio2->Divide(ND_highm_highs,ND_null_hmhs,1,1,"B");
        Ratio2->SetTitle(";Smeared Neutrino Energy [GeV];#splitline{Ratio}{Osc.-to-Unosc.} ");
        Ratio2->GetXaxis()->SetTitleOffset(0.95);
        Ratio2->GetYaxis()->SetTitleOffset(0.4);
        Ratio2->GetXaxis()->SetTitleFont(62);
        Ratio2->GetYaxis()->SetTitleFont(62);
        Ratio2->GetYaxis()->CenterTitle();
        Ratio2->GetXaxis()->CenterTitle();
        Ratio2->GetXaxis()->SetTitleSize(0.18);
        Ratio2->GetXaxis()->SetLabelSize(0.16);
        Ratio2->GetXaxis()->SetLabelOffset(0.008);
        Ratio2->GetYaxis()->SetTitleSize(0.1);
        Ratio2->GetYaxis()->SetTitleOffset(0.59);
        Ratio2->GetYaxis()->SetLabelSize(0.1);
	Ratio2->GetYaxis()->SetNdivisions(505);
	Ratio2->SetStats(kFALSE);
        Ratio2->SetMinimum(0.89);
        Ratio2->SetMaximum(1.11);
        Ratio2->SetMarkerColor(kBlue-3);
        Ratio2->SetLineColor(kBlue-3);
        Ratio2->SetMarkerStyle(20);
	Ratio2->SetMarkerSize(1);
        Ratio2->Draw("p");
        straight->Draw("l same");

	c5->Print("Evt_Dis_100m_44.pdf");


     

        TCanvas* c6 = new TCanvas("c6","Far Detector - Small #Delta m^{2}, Small sin^{2}(2#theta)",700,700);
        c6->cd();

        TPad *pad5 = new TPad("pad5","pad1",0,0.25,1,1);
        pad5 -> SetBottomMargin(0);
	TPad *pad6 = new TPad("pad6","pad2",0,0,1,0.25);
        pad6 -> SetTopMargin(0);
        pad6 -> SetBottomMargin(0.35);

        pad5 -> Draw();
	pad6 -> Draw();
        pad5 -> cd();


        FD_null_lmls->SetTitle(";Smeared Neutrino Energy [GeV];Events / Bin");
        FD_null_lmls->GetXaxis()->SetTitleOffset(1.2);
        FD_null_lmls->GetYaxis()->SetTitleOffset(1.2);
        FD_null_lmls->GetXaxis()->SetTitleFont(62);
        FD_null_lmls->GetYaxis()->SetTitleFont(62);
        FD_null_lmls->GetYaxis()->CenterTitle();
        FD_null_lmls->GetXaxis()->CenterTitle();
        FD_null_lmls->GetXaxis()->SetTitleSize(0.06);
        FD_null_lmls->GetXaxis()->SetLabelSize(0.055);
        FD_null_lmls->GetXaxis()->SetLabelOffset(-0.003);
        FD_null_lmls->GetYaxis()->SetTitleSize(0.06);
        FD_null_lmls->GetYaxis()->SetLabelSize(0.055);
        FD_null_lmls->SetStats(kFALSE);

        FD_null_lmls->SetLineWidth(2);
        FD_null_lmls->SetLineColor(kBlack);
        FD_null_lmls->SetMinimum(0.5);

	FD_null_lmls->Draw("h");

	FD_lowm_lows->SetMarkerColor(kRed);
	FD_lowm_lows->SetMarkerStyle(22);
	FD_lowm_lows->SetMarkerSize(1.5);
	FD_lowm_lows->Draw("p1 same");
        tex_stat->Draw();
        tex_pre->Draw();

        TLatex *tex_FD = new TLatex(.6,.85,"T600 (600m), on axis");
        tex_FD->SetNDC();
        tex_FD->SetTextFont(62);
	tex_FD->SetTextSize(0.047);
        tex_FD->Draw();

        TLatex *tex_FDPOT = new TLatex(.6,.77,"P.O.T. = 6.6 #times 10^{20}");
	tex_FDPOT->SetNDC();
        tex_FDPOT->SetTextFont(62);
        tex_FDPOT->SetTextSize(0.05);
        tex_FDPOT->Draw();


        tex_dm_lmls->Draw();
	tex_sin_lmls->Draw();
        leg->Draw();


        pad6 -> cd();

        Ratio3->Divide(FD_lowm_lows,FD_null_lmls,1,1,"B");
        Ratio3->SetTitle(";Smeared Neutrino Energy [GeV];#splitline{Ratio}{Osc.-to-Unosc.} ");
        Ratio3->GetXaxis()->SetTitleOffset(0.95);
        Ratio3->GetYaxis()->SetTitleOffset(0.4);
        Ratio3->GetXaxis()->SetTitleFont(62);
        Ratio3->GetYaxis()->SetTitleFont(62);
        Ratio3->GetYaxis()->CenterTitle();
        Ratio3->GetXaxis()->CenterTitle();
        Ratio3->GetXaxis()->SetTitleSize(0.18);
        Ratio3->GetXaxis()->SetLabelSize(0.16);
        Ratio3->GetXaxis()->SetLabelOffset(0.008);
        Ratio3->GetYaxis()->SetTitleSize(0.1);
        Ratio3->GetYaxis()->SetTitleOffset(0.59);
        Ratio3->GetYaxis()->SetLabelSize(0.1);
	Ratio3->GetYaxis()->SetNdivisions(505);
	Ratio3->SetStats(kFALSE);
        Ratio3->SetMinimum(0.95);
        Ratio3->SetMaximum(1.05);
        Ratio3->SetMarkerColor(kBlue-3);
        Ratio3->SetLineColor(kBlue-3);
        Ratio3->SetMarkerStyle(20);
	Ratio3->SetMarkerSize(1);
        Ratio3->Draw("p");
        straight->Draw("l same");

        TCanvas* c7 = new TCanvas("c7","Far Detector - High #Delta m^{2}, High sin^{2}(2#theta)",700,700);
        c7->cd();

        TPad *pad8 = new TPad("pad8","pad1",0,0.25,1,1);
        pad8 -> SetBottomMargin(0);
	TPad *pad9 = new TPad("pad9","pad2",0,0,1,0.25);
        pad9 -> SetTopMargin(0);
        pad9 -> SetBottomMargin(0.35);

        pad8 -> Draw();
	pad9 -> Draw();
        pad8 -> cd();


        FD_null_hmhs->SetTitle(";Smeared Neutrino Energy [GeV];Events / Bin");
        FD_null_hmhs->GetXaxis()->SetTitleOffset(1.2);
        FD_null_hmhs->GetYaxis()->SetTitleOffset(1.2);
        FD_null_hmhs->GetXaxis()->SetTitleFont(62);
        FD_null_hmhs->GetYaxis()->SetTitleFont(62);
        FD_null_hmhs->GetYaxis()->CenterTitle();
        FD_null_hmhs->GetXaxis()->CenterTitle();
        FD_null_hmhs->GetXaxis()->SetTitleSize(0.06);
        FD_null_hmhs->GetXaxis()->SetLabelSize(0.055);
        FD_null_hmhs->GetXaxis()->SetLabelOffset(-0.003);
        FD_null_hmhs->GetYaxis()->SetTitleSize(0.06);
        FD_null_hmhs->GetYaxis()->SetLabelSize(0.055);
        FD_null_hmhs->SetStats(kFALSE);

        FD_null_hmhs->SetLineWidth(2);
        FD_null_hmhs->SetLineColor(kBlack);
        FD_null_hmhs->SetMinimum(0.5);
	//        FD_null_hmhs->GetYaxis()->SetRangeUser(0.8*(FD_null_hmhs->GetBinContent(FD_null_hmhs->GetMinimumBin())),1.1*(FD_null_hmhs->GetBinContent(FD_null_hmhs->GetMaximumBin())));
	FD_null_hmhs->GetYaxis()->SetRangeUser(1000,400000);
	FD_null_hmhs->Draw("h");

        FD_highm_highs->SetMarkerColor(kRed);
        FD_highm_highs->SetMarkerStyle(22);
        FD_highm_highs->SetMarkerSize(1.5);
	FD_highm_highs->SetLineColor(kRed);
	FD_highm_highs->SetLineWidth(1.3);
	FD_highm_highs->Draw("p1 same");

	tex_stat->Draw();
        tex_pre->Draw();
        tex_FD->Draw();
        tex_FDPOT->Draw();


        tex_dm_hmhs->Draw();
	tex_sin_hmhs->Draw();

        leg->Draw();

        pad9 -> cd();

        Ratio4->Divide(FD_highm_highs,FD_null_hmhs,1,1,"B");
        Ratio4->SetTitle(";Smeared Neutrino Energy [GeV];#splitline{Ratio}{Osc.-to-Unosc.} ");
        Ratio4->GetXaxis()->SetTitleOffset(0.95);
        Ratio4->GetYaxis()->SetTitleOffset(0.4);
        Ratio4->GetXaxis()->SetTitleFont(62);
        Ratio4->GetYaxis()->SetTitleFont(62);
        Ratio4->GetYaxis()->CenterTitle();
        Ratio4->GetXaxis()->CenterTitle();
        Ratio4->GetXaxis()->SetTitleSize(0.18);
        Ratio4->GetXaxis()->SetLabelSize(0.16);
        Ratio4->GetXaxis()->SetLabelOffset(0.008);
        Ratio4->GetYaxis()->SetTitleSize(0.1);
        Ratio4->GetYaxis()->SetTitleOffset(0.59);
        Ratio4->GetYaxis()->SetLabelSize(0.1);
	Ratio4->GetYaxis()->SetNdivisions(505);
	Ratio4->SetStats(kFALSE);
        Ratio4->SetMinimum(0.9);
        Ratio4->SetMaximum(1.1);
        Ratio4->SetMarkerColor(kBlue-3);
        Ratio4->SetLineColor(kBlue-3);
        Ratio4->SetMarkerStyle(20);
	Ratio4->SetMarkerSize(1);
        Ratio4->Draw("p");
        straight->Draw("l same");

	c7->Print("Evt_Dis_600m_onaxis_44.pdf");
			
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


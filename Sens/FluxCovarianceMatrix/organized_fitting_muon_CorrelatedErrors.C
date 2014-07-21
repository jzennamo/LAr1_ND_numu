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
	char fileSource[200]="/uboone/app/users/jaz8600/work/lar1/SensitivityCode/Ntuple_Reprocessing/MatrixFiles/";
	//Note: parts of this routine will write to that same folder, make sure its writable!

	//Set up the bins for the vector:
	//coming soon...

	bool verbose 	=true;	//verbose causes print outs of information that is good for checking
	bool debug 	=true;	//debug prints out info to make sure root is behaving...

	char *  label = "#nu mode, CC events";
    
	TString fileNameRoot = fileSource;
	fileNameRoot += "numu_dis_ND_MicroBooNE";
	double flatSystematicError = 0.15;		// Only used if nearDetStats = false.
        double fluxSystematicError = 0.04;

	std::string mode = "nu";	//beam mode to run in
	bool use100m = true;		//Include the detector at 100m?
	bool use470m = false;		//Include the detector at 470m?
	bool use700m = true;		//Include the detector at 700m?

	bool shape_only = false;
	bool I_break_now = true;
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
	Int_t spoints = 1000;
     
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
	std::vector< std::vector <float> >   NCVec;
	std::vector< std::vector <float> >   LEVec;
	std::vector< std::vector < std::vector <float> > >   OscVec;

 //   std::vector< std::vector <float> >   NomVec;
  //  std::vector< std::vector < std::vector < std::vector < float > > > > SystVec;

	
	double NomVec[nL][20];
	double SystVec[nL][1001][7][20];

	NULLVec.resize(nL);
	NCVec.resize(nL);
	LEVec.resize(nL);
	OscVec.resize(nL);
//    NomVec.resize(nL);
//    SystVec.resize(nL);
    
	for( int l = 0; l < nL; l++){
	  OscVec[l].resize(npoints+1);
  //    SystVec[l].resize(1000);
  //      for(int u = 0; u < SystVec[0].size(); u++){
  //          SystVec[l][u].resize(7);
  //          }
  //     }
	}
	

    
    
	int counter = 0;
	int nbinsE = 0;

	std::cout << "Beginning : ... " << std::endl;

	if (use100m){
	  std::string temp_name = "../MatrixFiles/combined_ntuple_100m_nu_processed_numu_Joseph_Smeared.root";
	  TFile temp_file(temp_name.c_str());

	  std::string temp_name_syst = "../MatrixFiles/combined_ntuple_100m_nu_processed_numu.root";
	  TFile temp_file_syst(temp_name_syst.c_str());
	  TH1D *NULL_100;
        
	  NULL_100 = (TH1D*)(temp_file.Get("NumuCC"));
	  TH1D *Nom_100;

	  Nom_100 = (TH1D*)(temp_file_syst.Get("NumuCC"));
	  nbinsE = NULL_100->GetNbinsX();
	 
	  for(int i = 1; i <= nbinsE; i++){
	    NULLVec[counter].push_back(1*(NULL_100->GetBinContent(i)));
	    NomVec[counter][i-1] = 1*(Nom_100->GetBinContent(i));
	  }

	  for(int dm = 0; dm <= npoints; dm++){
	    TH1D *OSC_100;
	    std::string dmpoint = std::to_string(dm);
	    std::string name = "Osc_";
	    name = name+dmpoint;
	    OSC_100 = (TH1D*)(temp_file.Get(name.c_str()));
	    OSC_100->Rebin(1);
	    for(int i = 1; i <= nbinsE; i++){
	      OscVec[counter][dm].push_back(1*(OSC_100->GetBinContent(i)));
	    }
	    delete OSC_100;
	  }

	  for(int u = 0; u < spoints; u++){
	    for(int s = 0; s < 7; s++){
	      TH1D *SYST_100;
	      TString upoint = Form("%d",u);
	      TString name = "Universe_";
	      TString name2 = "_MultiSim_";
	      TString mul = Form("%d",s);
	      
	      name += upoint;
	      name += name2;
	      name += mul;

	      SYST_100 = (TH1D*)(temp_file_syst.Get(name));
	      SYST_100->Rebin(1);
	      for(int i = 1; i <= nbinsE; i++){
		SystVec[counter][u][s][i-1] = 1*(SYST_100->GetBinContent(i));
	      }
	      delete SYST_100;
	    } 
	  }


	  delete NULL_100;
	  temp_file.Close();
	  counter++;
	}

	if (use470m){
	  std::string temp_name = "../MatrixFiles/combined_ntuple_470m_nu_processed_numu_Joseph_Smeared.root";
	  TFile temp_file(temp_name.c_str());
	  
	  std::string temp_name_syst = "../MatrixFiles/combined_ntuple_470m_nu_processed_numu.root";
	  TFile temp_file_syst(temp_name_syst.c_str());
	  TH1D *NULL_470;
	  
	  NULL_470 = (TH1D*)(temp_file.Get("NumuCC"));
	  TH1D *Nom_470;
	  
	  Nom_470 = (TH1D*)(temp_file_syst.Get("NumuCC"));
	  nbinsE = NULL_470->GetNbinsX();
	  
	  for(int i = 1; i <= nbinsE; i++){
	    NULLVec[counter].push_back(NULL_470->GetBinContent(i));
	    NomVec[counter][i-1] = (Nom_470->GetBinContent(i));
	  }
	  
	  for(int dm = 0; dm <= npoints; dm++){
	    TH1D *OSC_470;
	    std::string dmpoint = std::to_string(dm);
	    std::string name = "Osc_";
	    name = name+dmpoint;
	    OSC_470 = (TH1D*)(temp_file.Get(name.c_str()));
	    OSC_470->Rebin(1);
	    for(int i = 1; i <= nbinsE; i++){
	      OscVec[counter][dm].push_back(OSC_470->GetBinContent(i));
	    }
	    delete OSC_470;
	  }

	  for(int u = 0; u < spoints; u++){
	    for(int s = 0; s < 7; s++){
	      TH1D *SYST_470;
	      TString upoint = Form("%d",u);
	      TString name = "Universe_";
	      TString name2 = "_MultiSim_";
	      TString mul = Form("%d",s);
	      
	      name += upoint;
	      name += name2;
	      name += mul;

	      SYST_470 = (TH1D*)(temp_file_syst.Get(name));
	      SYST_470->Rebin(1);
	      for(int i = 1; i <= nbinsE; i++){
		SystVec[counter][u][s][i-1] = (SYST_470->GetBinContent(i));
	      }
	      delete SYST_470;
	    } 
	  }


	  delete NULL_470;
	  temp_file.Close();
	  counter++;
	}


	//#####################################################################################
	//
	//   We will want to add MicroBooNE back in, but for now 2 detectors are enough
	//
	//###################################################################################


	if (use700m) {
	  std::string temp_name = "../MatrixFiles/combined_ntuple_600m_onaxis_nu_processed_numu_Joseph_Smeared.root";
	 TFile temp_file(temp_name.c_str());
	 std::string temp_name_syst = "../MatrixFiles/combined_ntuple_600m_onaxis_nu_processed_numu.root";
	  TFile temp_file_syst(temp_name_syst.c_str());
	  TH1D *NULL_600;
	  NULL_600 = (TH1D*)(temp_file.Get("NumuCC"));      

	  TH1D *Nom_600;
	  Nom_600 = (TH1D*)(temp_file_syst.Get("NumuCC"));
	  nbinsE = NULL_600->GetNbinsX();


	  for(int i = 1; i <= nbinsE; i++){
	    NULLVec[counter].push_back(NULL_600->GetBinContent(i));
	    NomVec[counter][i-1] = (Nom_600->GetBinContent(i));
	  }

	  for(int dm = 0; dm <= npoints; dm++){
	    TH1D *OSC_600;
	    std::string dmpoint = std::to_string(dm);
	    std::string name = "Osc_";
	    name = name+dmpoint;
	    OSC_600 = (TH1D*)(temp_file.Get(name.c_str()));
	    OSC_600->Rebin(1);
	    for(int i = 1; i <= nbinsE; i++){
	      OscVec[counter][dm].push_back(OSC_600->GetBinContent(i));
	    }
	    delete OSC_600;
	  }

	  for(int u = 0; u < spoints; u++){
	    for(int s = 0; s < 7; s++){
	      TH1D *SYST_600;
	      TString upoint = Form("%d",u);
	      TString name = "Universe_";
	      TString name2 = "_MultiSim_";
	      TString mul = Form("%d",s);
	      
	      name += upoint;
	      name += name2;
	      name += mul;

	      SYST_600 = (TH1D*)(temp_file_syst.Get(name));
	      SYST_600->Rebin(1);
	      for(int i = 1; i <= nbinsE; i++){
		SystVec[counter][u][s][i-1] = (SYST_600->GetBinContent(i));
	      }
	      delete SYST_600;
	    } 
	  }


	  delete NULL_600;
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
	

    //Build Matrices and Vectors
	TMatrixT <float> M(nbinsE*nL,nbinsE*nL);
	TMatrixT <float> M1(nbinsE*nL,nbinsE*nL);
	TMatrixT <float> CV(nbinsE*nL,1);
	TMatrix N(nbinsE*nL,1);
	double NT = 0;

	M.Zero();
	M1.Zero();
	CV.Zero();
	N.Zero();

	std::vector< std::vector< TMatrixT <float> > >  Pred;
	std::vector< float > CVvec;
	std::vector< std::vector< float > > Predvec;
	Predvec.resize(npoints+1);
	Pred.resize(npoints+1);

	for(int i = 0; i < npoints+1; i++){

	  Pred[i].resize(npoints+1);

	  for(int j = 0; j < npoints+1; j++){

	    Pred[i][j].ResizeTo(nbinsE*nL,1);

	  }
	}

	
	
	  
	  

	double high_DM2; 
	double low_DM2;
	double high_Sin2theta;
	double low_Sin2theta;

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
		
		M (Erri,Errj) /= n;
		

		// Asks if we are doing a "Shape+Rate" analysis 
		// If we are then we can just move forward
		if(!shape_only){

		  // Build "Fractional" Error Matrix
		  M (Erri,Errj) /= NomVec[Lrow][Erow]*NomVec[Lcol][Ecol]; 

		  // Now  take the real statistics to fill the Error Matrix
		  M(Erri, Errj) *= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];		

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

	    CV (Erri,0) = NULLVec[Lrow][Erow];
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
		  M(Erri, Errj) *= NULLVec[Lrow][Erow]*NULLVec[Lcol][Ecol];		
		  
		  if(Erri == Errj){ M (Erri, Errj) += NULLVec[Lrow][Erow];}


		  Errj++;		
	      }}
	      Erri++;
	    }}
			
	      }



       	TCanvas* c30 = new TCanvas("c30","",700,700);
	c30->SetLeftMargin(.15);
	c30->SetBottomMargin(.15);
	c30->SetTopMargin(.05);
	c30->SetRightMargin(.2);	
	M1.Draw("COLZ");
	/*

	TLatex *ND = new TLatex(.25,.075,"LAr1-ND");
	ND->SetNDC();
	ND->SetTextFont(62);
	ND->SetTextSize(0.05);
	ND->Draw();

	TLatex *MD = new TLatex(.5,.075,"T600");
	MD->SetNDC();
	MD->SetTextFont(62);
	MD->SetTextSize(0.05);
	MD->Draw();

	TLatex *ND45 = new TLatex(.1,.25,"LAr1-ND");
	ND45->SetNDC();
	ND45->SetTextAngle(90);
	ND45->SetTextFont(62);
	ND45->SetTextSize(0.05);
	ND45->Draw();

	TLatex *MD45 = new TLatex(.1,.65,"T600");
	MD45->SetNDC();
	MD45->SetTextAngle(90);
	MD45->SetTextFont(62);
	MD45->SetTextSize(0.05);
	MD45->Draw();
	
	*/
	std::cout << "...Error Matrix Filled" << std::endl;
       
	double sin22th;

	int Predi = 0;
	std::cout << "Filling Prediction Matrices..." << std::endl;
	for(int dm = 0; dm <= npoints; dm++){
	  for(int s = 0; s <= npoints; s++){
	    Predi = 0; 
	    for(int Lbin = 0; Lbin < nL; Lbin++){
	      for(int Ebin = 0; Ebin < nbinsE; Ebin++){
		
		sin22th = pow(10., (TMath::Log10(sin22thmin)+ (s * (1./npoints))*TMath::Log10(sin22thmax/sin22thmin)));

		Pred[dm][s] (Predi,0) = NULLVec[Lbin][Ebin] - (OscVec[Lbin][dm][Ebin])*(sin22th);
		
		Predi++;
		
	      }
	    }
	  }
	}

	std::cout << "...Prediction Matrices Filled" << std::endl;

	//Now create the inverse covariance matrix.  Root can invert for us:
	//inverse cov matrix, M^{-1}, used in chi2 calculation
	TMatrixT <float> cov(nbinsE*nL,nbinsE*nL);

	std::cout << "Invert in...3...2...1..." << std::endl;
	cov = M.Invert();
	std::cout << "Inverted!" << std::endl;

	//Time to compute the chi2
       
	double DMpoint;


	double cvi, cvi_total;
	double predictioni, pred_total;
	for(int dm2 = 0; dm2 <= npoints; dm2++){
	  DMpoint = pow(10., (TMath::Log10(dm2min)+ (dm2 * (1./npoints))*TMath::Log10(dm2max/dm2min)) );
	  
	  for(int sint = 0; sint <= npoints; sint++){
	    
	    sin22th = pow(10., (TMath::Log10(sin22thmin)+ (sint * (1./npoints))*TMath::Log10(sin22thmax/sin22thmin)));
	    
	    chisq = 0.0;
	      
	    double Scaling[nbinsE*nL];
	    
	    if(shape_only){
	      //Scale all measurements to the ND NULL levels             
	      for(int Edecbin = 0; Edecbin < nbinsE; Edecbin++){
		
		Scaling[Edecbin]  = Pred[dm2][sint] (Edecbin,0);
		if(CV(Edecbin,0) != 0){Scaling[Edecbin] /= CV(Edecbin,0);}
		if(CV(Edecbin,0) == 0){ std::cout << "Scaling failed....CV==0!" << std::endl; break;}
		
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
	    TMatrixT <float> CVscaled(nbinsE*nL,1);
	    

	    for(int i = 0; i < nbinsE*nL; i++){
	      CVscaled(i,0) = Scaling[i]*(CV(i,0));   
	    }

	    
	    TMatrixT <float> Final(nbinsE*nL,1);
            TMatrixT <float> FinalT(1,nbinsE*nL);



	    for(int Edecbin = 0; Edecbin < nbinsE*nL; Edecbin++){
	      
	      Final(Edecbin,0) = (CVscaled(Edecbin,0)) - ((Pred[dm2][sint] (Edecbin,0))); 
	      FinalT(0,Edecbin) = (CVscaled(Edecbin,0)) - ((Pred[dm2][sint] (Edecbin,0)));

	    }//Sum over all detectors and energies to determine the chi2 for that osc. point	      


	    TMatrixT <float> middle(1, nbinsE*nL);
	    
	    //middle = (&Final) * (&cov);
	    	
	    middle.Mult(FinalT,cov);
    
	    //	    TMatrixT(&Final, TMatrix::kMult, &cov);
	    // TMatrix *end = new TMatrixT(&middle, TMatrix::kMultTranspose, &Final);

	    //	    middle = (*Final)*(*cov);

            TMatrixT <float> Fin(1, 1);
	    
	    //end = (&middle)*(&(Final->Transpose()));
	  

  
	    Fin.Mult(middle,Final);

	    
	    chisq = Fin(0,0);
	    //	    std::cout << "chi2 : " << chisq << std::endl;
	    //    	      chi2->Fill(sin22th,DMpoint,chisq);

	      if (chisq <= deltachisq5s){
		y5s[dm2] = DMpoint; x5s[dm2] = sin22th;
		std::cout << "\Delta m^{2} : " << DMpoint << ", sin22theta : " << sin22th << std::endl;
		for(int i = 0; i < nbinsE; i++){
		  std::cout << i << ", SF = " << Scaling[i] << ";";
		}
		std::cout << " ####################### "<< std::endl;

	      }
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
	hr1->GetXaxis()->SetTitleSize(0.05);
        hr1->GetXaxis()->SetLabelSize(0.04);
	hr1->GetXaxis()->SetLabelOffset(0.001);
	hr1->GetYaxis()->SetTitleSize(0.05);
        hr1->GetYaxis()->SetLabelSize(0.04);
	hr1->SetStats(kFALSE);
	hr1->Draw();

	//	if(mode == "nu")        gROOT->ProcessLine(".x ./minisciboon_plot_numu.c(c3)");
	//if(mode == "nubar")     gROOT->ProcessLine(".x ./minisciboon_plot_numubar.c(c3)");
	
	//======================================================

       //	TLatex *tex_Detector = new TLatex(.2,.9,"2 #times LAr1-ND (100m), MicroBooNE (470m),");	
       // TLatex *tex_FDetector = new TLatex(.34,.86,"ICARUS (700m)");
       //TLatex *tex_Detector = new TLatex(.2,.9,"LAr1-ND (100m), T600 (700m)");
       //TLatex *tex_FDetector = new TLatex(.34,.86,"LAr1-FD (700m)");	  

	TLatex *tex_Detector = new TLatex(.2,.23,"#splitline{LAr1-ND (100m)}{and T600 (600m, on axis)}");
		//TLatex *tex_Detector = new TLatex(.2,.9, "T600 (600m)");
	tex_Detector->SetNDC();
	tex_Detector->SetTextFont(62);
        tex_Detector->SetTextSize(0.035);
	//	tex_Detector->SetTextSize(0.03);        
	tex_Detector->Draw(); 
	//tex_FDetector->SetNDC();
	//tex_FDetector->SetTextFont(62);
        //tex_FDetector->SetTextSize(0.03);
        //tex_FDetector->Draw();
   
	TLatex *tex_pre = new TLatex(.18,.78,"PRELIMINARY");
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

        TLatex *tex_un = new TLatex(.18,.89,"Statistical Uncertainty Only");
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


	/*        TLatex *tex_corr = new TLatex(.2,.94,"#rho = 0");
        tex_corr->SetNDC();
        tex_corr->SetTextFont(62);
        tex_corr->SetTextSize(0.04);
        tex_corr->Draw();*/

       	TLegend* legt=new TLegend(0.2,0.3,0.6,0.45);
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
	//sens3s->Draw("l same");
	legt->AddEntry(sens3s,"3#sigma CL","l");
       	//sens5s->Draw("l same");
	legt->AddEntry(sens5s,"5#sigma CL","l");
	TLine *gdummy3 = new TLine();
	gdummy3->SetLineColor(kRed-3);
        gdummy3->SetLineStyle(2);
        gdummy3->SetLineWidth(2);
	//	legt->AddEntry(gdummy3,"MiniBooNE + SciBooNE 90% CL","l");
	legt->Draw();
       
	/*    TLegend* leg3=new TLegend(0.2,0.2,0.4,0.35);
	      leg3->SetFillStyle(0);
	      leg3->SetFillColor(0);
	      leg3->SetBorderSize(0);
	      leg3->SetTextSize(0.03);
	      leg3->AddEntry("c3","MiniBooNE + SciBooNE 90% CL","L");
	      leg3->Draw();*/
	
	/*   TImage *img = TImage::Create();
	     img -> FromPad(c3);
	     if (useNearDetStats){
	     c3 -> Print(fileNameRoot+mode+"_nearDetStats.eps", "eps");
	     }
	     else{
	     c3 -> Print(fileNameRoot+mode+"_flatStats.eps", "eps");
	     }
	*/

        //c3 -> Print("ContourComparison/Sens_Matrix_LAr1ND_100m_T600_on_axis_ShapeOnly_NoFlux_StatMatch.C");

	c3 -> Print("Sens_Matrix_LAr1-ND_100m_T600_on_axis_Shape_and_Rate.pdf");
       	//c3 -> Print("Sens_Matrix_LAr1-ND_100m_T600_off_axis_ShapeOnly.pdf");


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

void organized_fitting_muon_CorrelatedErrors(){
  multiple_detector_fit();
  return;
}


#define NtupleReprocessing_cxx
#include "NtupleReprocessing.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <sys/stat.h>
#include <unistd.h>

void NtupleReprocessing::Loop(string signal, Int_t iDet, int ND_iDet, Long64_t max_entry)
{
//   In a ROOT session, you can do:
//      Root > .L NtupleReprocessing.C
//      Root > NtupleReprocessing t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  const static int npoints = 1000;
  double emin = 0.2, emax = 3.0; int ebins = 50;
  double bins[20] = {.200, .300, .400, .450, .500, .550, .600, .650, .700, .750, .800, .850, .900, .950, 1.000, 1.250, 1.500, 2.000, 2.500, 3.000};

  TH1D *numuCC = new TH1D("NumuCC","NumuCC;Reconstructed Neutrino Energy [GeV];Events",19,bins);    
      TH1D *numuCC_Syst[npoints+1][7];

      TH1D *pion_pT = new TH1D("pipT",";Pion pT [GeV];Events",50,0.1,1);

      TH1D *numuCC_MC = new TH1D("NumuCC_MC","NumuCC;Reconstructed Neutrino Energy [GeV];Events",19, bins);    
      TH1D *numuCC_Syst_MC[npoints+1][7];

      for(int i = 0; i <= npoints; i++){
	for(int j = 0; j <= 7; j++){
	numuCC_Syst[i][j] = (TH1D*)numuCC->Clone();
        numuCC_Syst_MC[i][j] = (TH1D*)numuCC_MC->Clone();	
	}
      }
      

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
  
   if(iDet != 0){ND_iDet = 1;}

   Double_t detect_dist = 0;   // 10000=ND, 47000=MicroBooNE, 70000=FD                                                                                                         
   if (iDet == 0 || iDet == 4) detect_dist = 10000.0;
   else if (iDet == 5) detect_dist = 15000.0;
   else if (iDet == 6) detect_dist = 20000.0;

   else if (iDet == 1) detect_dist = 47000.0;
   else if (iDet == 2) detect_dist = 70000.0;
   else if (iDet == 3) detect_dist = 60000.0;
   else if (iDet == 7) detect_dist = 60000.0;

   if (fChain == 0) return;
   
   b_iflux->GetEntry(0);
   
   bool isFid, isActive;
   int iChan = 0, iTop = 0;
   
   double nd_xmin;  double nd_xmax;
   double nd_ymin;  double nd_ymax;
   double nd_zmin;  double nd_zmax;
   
   if((iDet == 0 || iDet == 5 || iDet == 6) && ND_iDet == 1){
     
     nd_xmin = -200.0;
     nd_xmax =  200.0;
     nd_ymin = -242.0;
     nd_ymax =  158.0;
     nd_zmin =  0.0;
     nd_zmax =  365.0;
}
   
   if(iDet == 0 && ND_iDet == 4){
     nd_xmin = -55.0; nd_xmax =  315.0;
     nd_ymin = -172.0; nd_ymax =  172.0;
     nd_zmin =  0.0;   nd_zmax =  230.0;}

   
   if(iDet == 1){
     nd_xmin =  0.0;      nd_xmax =  256.0;
     nd_ymin = -116.5;      nd_ymax =  116.5;
     nd_zmin =  0.0;      nd_zmax =  1036.0;}
   
   if(iDet == 4){
     nd_xmin = -80.0; nd_xmax =  320.0;
     nd_ymin = -242.0; nd_ymax =  158.0;
     nd_zmin =  0.0;   nd_zmax =  2*365.0;
   }
   

   std::cout << "input file: " << InFile().Remove(InFile().Length()-5) << std::endl;

   double fluxweight = 1.0;
   double electron_cand_energy = 0.0;
   double electron_cand_angle = 0.0;
   double EnuReco = 0.;
   double enuccqe = 0.0;
   double enucalo1 = 0.0;
   double enucalo2 = 0.0;
   double LEccqe = 0.0;
   double LEcalo1 = 0.0;
   double LEcalo2 = 0.0;

   double potweight = utils.GetPOTNorm( iflux, iDet );



   if( max_entry == -1 ) max_entry = nentries;
   
   for (Long64_t jentry=0; jentry<nentries && jentry<max_entry;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if ( (jentry%100000) == 0){ std::cout << int(100*jentry/nentries) << "% done !"  << std::endl;}
   

      if(abs(inno) != 14) continue;
      
      TVector3 vtx(Vx, Vy, Vz);
      isFid = utils.IsFiducial( iDet, ND_iDet, vtx );
      isActive = utils.IsActive( iDet, ND_iDet, vtx );

      wgt = 0.0;
      
      if( isFid && isCC){
      
      fluxweight = utils.GetFluxWeight( energy, iflux, inno, ndecay );
      fluxweight *= potweight;
      wgt = fluxweight;
      
      }

      else{continue;}

     
     numuCC->Fill(enugen,wgt);
     numuCC_MC->Fill(enugen);
     pion_pT->Fill(sqrt(tpx*tpx+tpy*tpy),wgt);
     
     float syst = 0;

     for(int i = 0; i < npoints ; i++){
       for(int j = 0; j < 7; j++){
       
       syst = (MultiWeight->at(j)).at(i);
       
       numuCC_Syst[i][j]->Fill(enugen,syst*wgt);
       numuCC_Syst_MC[i][j]->Fill(enugen,syst);


       }
     }




   }//End of loop over events 

   TString outfile = InFile().Remove(InFile().Length()-5) + "_processed_" + signal + ".root";
   std::cout << "outfile file: " << outfile << std::endl;
   TFile *f = new TFile(outfile, "RECREATE");
   if (f->IsZombie()){
     std::cout << "Error, couldn't create output file." << std::endl;
     return;
   }

   numuCC->Write();
   numuCC_MC->Write();
   pion_pT->Write();
   for(int u = 0; u < npoints; u++){
     for(int s = 0; s < 7; s++){
     std::string uni = std::to_string(u);
     std::string name = "Universe_";
     std::string name2 = "_MultiSim_";
     std::string mul = std::to_string(s);
     std::string mc = "_MC";
     std::string start = ""; 
     
     start = name+uni+name2+mul+mc;
     name = name+uni+name2+mul;
     numuCC_Syst[u][s]->Write(name.c_str());
     numuCC_Syst_MC[u][s]->Write(start.c_str());

     }
   }

   f->Write();
   f->Close();

} //End of program

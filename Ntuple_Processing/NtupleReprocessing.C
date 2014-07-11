#define NtupleReprocessing_cxx
#include "NtupleReprocessing.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <sys/stat.h>
#include <unistd.h>


void NtupleReprocessing::Loop( string signal, Int_t iDet, int ND_iDet, Long64_t max_entry){

  // Signal options are:
  // numu - muon neutrino disappearance 

  //Energy Resolutions AND ND configurations 


  //ND_iDet corresponds to different configurations of the NearDetector
  //ND_iDet == 1 is just the argon, no ranger
  //ND_iDet == 2 is a non-magetized ranger tacked onto the end of the argon
  //ND_iDet == 3 is a magneized ranger tacked onto the end of the argon
  //ND_iDet == 4 is a totally magnetized argon volume 
  if(iDet != 0){ND_iDet = 1;}


  //Pion-Muon Separation:

  double muEff_pi = 1;
  double piEff_pi = 1;
 
  //

  double electronIDeff = 0.8;
  double muonIDeff     = 0.8;
  double photonMisID   = 0.06;
  double muonCCMisID   = 0.001;

  // NC photon vertex energy cuts
  double vtxEcut = 0.0; //0.025;   // GeV
  double convDistCut = 0.0; //5.0; // cm

  double egammaThreshold = 0.0; // GeV

  Double_t detect_dist = 0;   // 10000=ND, 47000=MicroBooNE, 70000=FD
  if (iDet == 0 || iDet == 4) detect_dist = 10000.0;
  else if (iDet == 1) detect_dist = 47000.0;
  else if (iDet == 2) detect_dist = 70000.0;
  else if (iDet == 3) detect_dist = 60000.0;
  else if (iDet == 5) detect_dist = 15000.0;
  else if (iDet == 6) detect_dist = 20000.0;
  else if (iDet == 7) detect_dist = 60000.0;



  //---------------------------------------------

  //   In a ROOT session, you can do:
  //      Root > .L NtupleReprocessing.C+ (+ to compile - this is necessary)
  //      Root > NtupleReprocessing t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //
  
  if (fChain == 0) return;

  //----------------------------------
  // Create the output file
  //----------------------------------
  std::cout << "input file: " << InFile().Remove(InFile().Length()-5) << std::endl;

  //------------------------------------------------------------------
  // These are new variables that get added at the reprocessing stage

  double efficiency = 0.0;
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

  std::vector<double> ConversionDist;

  //Get the POT weight.  Need to get iflux
  b_iflux->GetEntry(0); //iflux should now be filled

  double potweight = utils.GetPOTNorm( iflux, iDet );

  bool isFid, isActive;
  int iChan = 0, iTop = 0;

  TVector3 MuonExitPos(0,0,0), MuonExitMom(0,0,0), MuonInitialMom(0,0,0);
  TVector3 ChargedPionExitPos(0,0,0), ChargedPionExitMom(0,0,0), ChargedPionInitialMom(0,0,0);

  // initialize counters and histograms

  double muonCC = 0;
  double electronCC = 0;
  double NCpizero = 0;
  double foundBothPhotons = 0, foundOnePhoton = 0, foundNoPhotons = 0; 
  double multiPizeroEvents = 0;
  double NCsinglePhoton = 0;
  std::vector<double> singlePhotonE, singlePhotonTheta;

  double emin = 0.2, emax = 3.0;
  int ebins = 50;
  const int npoints = 500;
  Double_t dm2min = 0.01;                       //eV**2
  Double_t dm2max = 100.;                       //eV**2
  double mu_L = 0;
  double DMpoint;
  double mu_E = 0;
  double mu_Theta = 0;
  double smeared_mu_E = 0;

  double Ymax = 0, Ymin = 0, Xmin = 0, Xmax = 0, Zmin = 0, Zmax = 0, stopped = 0;
  double Ymax_F = 0, Ymin_F = 0, Xmin_F = 0, Xmax_F = 0, Zmin_F = 0, Zmax_F = 0;

  double Unmeas = 0, Cont = 0, MS = 0, LArB = 0, RangerB = 0;

  double nd_xmin;  double nd_xmax;
  double nd_ymin;  double nd_ymax;
  double nd_zmin;  double nd_zmax;

  if((iDet == 0 || iDet == 5 || iDet == 6)&& ND_iDet == 1){

    nd_xmin = -80.0; nd_xmax =  320.0;
    nd_ymin = -242.0; nd_ymax =  158.0;
    nd_zmin =  0.0;   nd_zmax =  365.0;}

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


  double Res = 0; 
  bool ranger = false;
  double bins[20] = {.200, .300, .400, .450, .500, .550, .600, .650, .700, .750, .800, .850, .900, .950, 1.000, 1.250, 1.500, 2.000, 2.500, 3.000};

  TH1D *MC_mup = new TH1D("MC_mup","; MC Truth Muon Energy [GeV];Events",60,0,1.5);
  TH1D *MC_mupZ = new TH1D("MC_mupZ","; MC Truth Muon Energy in z-Direction [GeV];Events",60,0,1.5);
  TH1D *MC_mupT = new TH1D("MC_mupT","; MC Truth Muon Energy in radial-Direction[GeV];Events",60,0,1.5);

  TH1D *MC_mur = new TH1D("MC_mur","; MC Truth Muon Length [cm];Events",1000,0,1000);
  TH1D *MC_murZ = new TH1D("MC_murZ","; MC Truth Muon Length in z-Direction [cm];Events",1000,0,1000);
  TH1D *MC_murT = new TH1D("MC_murT","; MC Truth Muon Length in radial-Direction [cm];Events",1000,0,1000);

  TH1D *MC_muT = new TH1D("MC_muT","; MC Truth Muon Angle (#theta);Events",300,0,180);
  TH1D *MC_mucosT = new TH1D("MC_mucosT","; MC Truth Muon cos(#theta);Events",250,0,6.2);


  TH1D *numuCC = new TH1D("NumuCC","NumuCC;Reconstructed Neutrino Energy [GeV];Events",19, bins);
  TH1D *numuCC_truth = new TH1D("NumuCC_True","NumuCC;True Neutrino Energy [GeV];Events",19, bins);
  TH1D *numuCC_cont = new TH1D("NumuCC_cont","NumuCC_cont;Reconstructed Neutrino Energy [GeV];Events",19, bins);
  TH1D *numuCC_cont_truth = new TH1D("NumuCC_cont_True","NumuCC_cont;Reconstructed Neutrino Energy [GeV];Events",19, bins);
  TH1D *numuCC_exit = new TH1D("NumuCC_exit","NumuCC_exit;Reconstructed Neutrino Energy [GeV];Events",19, bins);
  TH1D *numuCC_exit_truth = new TH1D("NumuCC_exit_True","NumuCC_exit;True Neutrino Energy [GeV];Events",19, bins);

  TH1D *numuNC = new TH1D("NumuNC","NumuCC;Reconstructed Neutrino Energy [GeV];Events",19, bins);
  TH1D *numuNC_truth = new TH1D("NumuNC_True","NumuCC;True Neutrino Energy [GeV];Events",19, bins);
  TH1D *numuNC_cont = new TH1D("NumuNC_cont","NumuCC_cont;Reconstructed Neutrino Energy [GeV];Events",19, bins);
  TH1D *numuNC_cont_truth = new TH1D("NumuNC_cont_True","NumuCC_cont;Reconstructed Neutrino Energy [GeV];Events",19, bins);
  TH1D *numuNC_exit = new TH1D("NumuNC_exit","NumuCC_exit;Reconstructed Neutrino Energy [GeV];Events",19, bins);
  TH1D *numuNC_exit_truth = new TH1D("NumuNC_exit_True","NumuCC_exit;True Neutrino Energy [GeV];Events",19, bins);

  TH1D *numuCCLepE = new TH1D("NumuCCLepE", "NumuCCLepE;Total Muon Energy [GeV];Events",19, bins);
  TH1D *numuCCLepE_cont = new TH1D("NumuCCLepE_cont", "NumuCCLepE_cont;Contained Muon Energy [GeV];Events",19, bins);
  TH1D *numuCCLepE_exit = new TH1D("NumuCCLepE_exit", "NumuCCLepE_exit;Exiting Muon Energy [GeV];Events",19, bins);

  TH1D *numuNCLepE = new TH1D("NumuNCLepE", "NumuCCLepE;Total Muon Energy [GeV];Events",19, bins);
  TH1D *numuNCLepE_cont = new TH1D("NumuNCLepE_cont", "NumuCCLepE_cont;Contained Muon Energy [GeV];Events",19, bins);
  TH1D *numuNCLepE_exit = new TH1D("NumuNCLepE_exit", "NumuCCLepE_exit;Exiting Muon Energy [GeV];Events",19, bins);

  TH1D *numuCCLepTheta = new TH1D("NumuCCLepTheta", "NumuCCLepE; Muon Angle (#theta);Events",20*ebins,0,180);
  TH1D *numuCCLepTheta_cont = new TH1D("NumuCCLepTheta_cont", "NumuCCLepE_cont;Contained Muon Angle (#theta);Events",20*ebins,0,180);
  TH1D *numuCCLepTheta_exit = new TH1D("NumuCCLepTheta_exit", "NumuCCLepE_exit;Exiting Muon Angle (#theta);Events",20*ebins,0,180);

  TH1D *numuNCLepTheta = new TH1D("NumuNCLepTheta", "NumuCCLepE; Muon Angle (#theta);Events",20*ebins,0,180);
  TH1D *numuNCLepTheta_cont = new TH1D("NumuNCLepTheta_cont", "NumuCCLepE_cont;Contained Muon Angle (#theta);Events",20*ebins,0,180);
  TH1D *numuNCLepTheta_exit = new TH1D("NumuNCLepTheta_exit", "NumuCCLepE_exit;Exiting Muon Angle (#theta);Events",20*ebins,0,180);

  TH1D *numuCCLepLength = new TH1D("NumuCCLepLength", "NumuCCLepE; Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuCCLepLength_cont = new TH1D("NumuCCLepLength_cont", "NumuCCLepE_cont;Contained Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuCCLepLength_exit = new TH1D("NumuCCLepLength_exit", "NumuCCLepE_exit;Exiting Muon Angle (#theta);Events",20*ebins,0,nd_zmax);

  TH1D *numuNCLepLength = new TH1D("NumuNCLepLength", "NumuCCLepE; Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuNCLepLength_cont = new TH1D("NumuNCLepLength_cont", "NumuCCLepE_cont;Contained Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuNCLepLength_exit = new TH1D("NumuNCLepLength_exit", "NumuCCLepE_exit;Exiting Muon Angle (#theta);Events",20*ebins,0,nd_zmax);

  TH1D *numuCCLepZLength = new TH1D("NumuCCLepZLength", "NumuCCLepE; Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuCCLepZLength_cont = new TH1D("NumuCCLepZLength_cont", "NumuCCLepE_cont;Contained Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuCCLepZLength_exit = new TH1D("NumuCCLepZLength_exit", "NumuCCLepE_exit;Exiting Muon Angle (#theta);Events",20*ebins,0,nd_zmax);

  TH1D *numuNCLepZLength = new TH1D("NumuNCLepZLength", "NumuCCLepE; Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuNCLepZLength_cont = new TH1D("NumuNCLepZLength_cont", "NumuCCLepE_cont;Contained Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuNCLepZLength_exit = new TH1D("NumuNCLepZLength_exit", "NumuCCLepE_exit;Exiting Muon Angle (#theta);Events",20*ebins,0,nd_zmax);

  TH1D *numuCCLepRLength = new TH1D("NumuCCLepRLength", "NumuCCLepE; Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuCCLepRLength_cont = new TH1D("NumuCCLepRLength_cont", "NumuCCLepE_cont;Contained Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuCCLepRLength_exit = new TH1D("NumuCCLepRLength_exit", "NumuCCLepE_exit;Exiting Muon Angle (#theta);Events",20*ebins,0,nd_zmax);

  TH1D *numuNCLepRLength = new TH1D("NumuNCLepRLength", "NumuCCLepE; Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuNCLepRLength_cont = new TH1D("NumuNCLepRLength_cont", "NumuCCLepE_cont;Contained Muon Angle (#theta);Events",20*ebins,0,nd_zmax);
  TH1D *numuNCLepRLength_exit = new TH1D("NumuNCLepRLength_exit", "NumuCCLepE_exit;Exiting Muon Angle (#theta);Events",20*ebins,0,nd_zmax);

  TH1D *nuE_all = new TH1D("NuEnergy_all",";True Neutrino Energy [GeV];Events",20,emin,emax);
  TH1D *nuE_accepted = new TH1D("NuEnergy_accepted",";True Neutrino Energy [GeV] Events",20,emin,emax);

  TH1D *muE_all = new TH1D("MuEnergy_all",";Inital Muon Energy [GeV];Events",20,0,emax);
  TH1D *muE_accepted = new TH1D("MuEnergy_accepted",";Initial Muon Energy [GeV];Events",20,0,emax);

  TH1D *muTheta_all = new TH1D("MuTheta_all",";Inital Muom Angle (#theta);Events",20,0,3.14);
  TH1D *muTheta_accepted = new TH1D("MuTheta_accepted",";Initial Muon Angle (#theta);Events",20,0,3.14);

  TH2D *muTE_all = new TH2D("MuThetaEnergy_all",";Inital Muom Angle (#theta);Events",120,0,3.14, 120,emin,emax);
  TH2D *muTE_accepted = new TH2D("MuThetaEnergy_accepted",";Initial Muon Angle (#theta);Events",120,0,3.14, 120,emin,emax);

  TH2D *muE_Reco_True = new TH2D("MuE_ReTr",";;Events",120,0,1, 120,0,1);
  TH2D *nuE_Reco_True = new TH2D("nuE_ReTr",";;Events",120,0,1, 120,0,1);

  TH1D *Res_cont = new TH1D("res_v_L_cont",";Internal Track Length [cm];#delta p / p",10,0,300);
  TH1D *Res_exit = new TH1D("res_v_L_exit",";Internal Track Length [cm];Events",10,0,300);
  TH1D *Res_all = new TH1D("res_v_L_all",";Internal Track Length [cm];Events",10,0,300);

  TH1D *Res_cont_N = new TH1D("res_v_L_cont_N",";Internal Track Length [cm];#delta p / p",10,0,300);
  TH1D *Res_exit_N = new TH1D("res_v_L_exit_N",";Internal Track Length [cm];Events",10,0,300);
  TH1D *Res_all_N = new TH1D("res_v_L_all_N",";Internal Track Length [cm];Events",10,0,300);

  TH1D *Pion_L = new TH1D("all_pi_L",";Internal Track Length [cm];Events",300,0,300);
  TH1D *Pion_cont_L = new TH1D("cont_pi_L",";Internal Track Length [cm];Events",300,0,300);

  TH1D *Muon_D = new TH1D("dis_vtx_muon",";Internal Track Length [cm];Events",300,0,300);
  TH1D *Pion_D = new TH1D("dis_vtx_pion",";Internal Track Length [cm];Events",300,0,300);
  TH1D *N_Exit_Pion = new TH1D("num_pion_exit",";Number of exiting pions;Events",7,-0.5,6.5);
  TH1D *N_Exit_Pion_CC = new TH1D("num_pion_exit_cc",";Number of exiting pions;Events",7,-0.5,6.5);
  TH1D *N_Cont_Pion_CC = new TH1D("num_pion_cont_cc",";Number of exiting pions;Events",7,-0.5,6.5);

  TH1D *muLx_all = new TH1D("MuLx_all",";Interaction Vertex Location [cm];Events",20,nd_xmin,nd_xmax);
  TH1D *muLx_accepted = new TH1D("MuLx_accepted",";Interaction Vertex Location [cm];Events",20,nd_xmin,nd_xmax);

  TH1D *muLy_all = new TH1D("MuLy_all",";Interaction Vertex Location [cm];Events",20,nd_ymin,nd_ymax);
  TH1D *muLy_accepted = new TH1D("MuLy_accepted",";Interaction Vertex Location [cm];Events",20,nd_ymin,nd_ymax);

  TH1D *muLz_all = new TH1D("MuLz_all",";Interaction Vertex Location [cm];Events",20,nd_zmin,nd_zmax);
  TH1D *muLz_accepted = new TH1D("MuLz_accepted",";Interaction Vertex Location [cm];Events",20,nd_zmin,nd_zmax);

  TH2D *nuEvmuE_Corr = new TH2D("nuEvmuE_Corr",";;Events", 20,emin,emax,20,emin,emax);
  TH2D *nuEvmuE_Corr_cont = new TH2D("nuEvmuE_Corr_cont",";;Events", 20,emin,emax,20,emin,emax);
  TH2D *nuEvmuE_Corr_exit = new TH2D("nuEvmuE_Corr_exit",";;Events", 20,emin,emax,20,emin,emax);

  TH2D *nuEvmuE_Corr_MC = new TH2D("nuEvmuE_Corr_MC",";;Events", 20,emin,emax,20,emin,emax);
  TH2D *nuEvmuE_Corr_MC_cont = new TH2D("nuEvmuE_Corr_MC_cont",";;Events", 20,emin,emax,20,emin,emax);
  TH2D *nuEvmuE_Corr_MC_exit = new TH2D("nuEvmuE_Corr_MC_exit",";;Events", 20,emin,emax,20,emin,emax);

  TH2D *nuParP_T_Z = new TH2D("nuParP_T_Z",";Transverse; Longitudinal;Events", 20*ebins,0,3,20*ebins,0,10);
  TH2D *nuParP_T_Z_cont = new TH2D("nuParP_T_Z_cont",";Transverse; Longitudinal;Events", 20*ebins,0,3,20*ebins,0,10);
  TH2D *nuParP_T_Z_exit = new TH2D("nuParP_T_Z_exit",";Transverse; Longitudinal;Events", 20*ebins,0,3,20*ebins,0,10);

  TH2D *muEvL_Corr = new TH2D("muEvL_Corr",";;Events", 90,0,1, 90, 0, 300);
  TH2D *muEvL_Corr_cont = new TH2D("muEvL_Corr_cont",";;Events", 90,0,1, 90, 0, 300);
  TH2D *muEvL_Corr_exit = new TH2D("muEvL_Corr_exit",";;Events", 90,0,1, 90, 100, 300);
  
  TH1D *nuE_all_cont = new TH1D("NuEnergy_all_cont",";True Neutrino Energy [GeV];Events",20,emin,emax);
  TH1D *nuE_accepted_cont = new TH1D("NuEnergy_accepted_cont",";True Neutrino Energy [GeV] Events",20,emin,emax);

  TH1D *muE_all_cont = new TH1D("MuEnergy_all_cont",";Inital Muon Energy [GeV];Events",20,0,emax);
  TH1D *muE_accepted_cont = new TH1D("MuEnergy_accepted_cont",";Initial Muon Energy [GeV];Events",20,0,emax);

  TH1D *muTheta_all_cont = new TH1D("MuTheta_all_cont",";Inital Muom Angle (#theta);Events",20,0,3.14);
  TH1D *muTheta_accepted_cont = new TH1D("MuTheta_accepted_cont",";Initial Muon Angle (#theta);Events",20,0,3.14);

  TH1D *nuE_all_exit = new TH1D("NuEnergy_all_exit",";True Neutrino Energy [GeV];Events",20,emin,emax);
  TH1D *nuE_accepted_exit  = new TH1D("NuEnergy_accepted_exit",";True Neutrino Energy [GeV] Events",20,emin,emax);

  TH1D *muE_all_exit  = new TH1D("MuEnergy_all_exit",";Inital Muon Energy [GeV];Events",20,0,emax);
  TH1D *muE_accepted_exit  = new TH1D("MuEnergy_accepted_exit",";Initial Muon Energy [GeV];Events",20,0,emax);

  TH1D *muTheta_all_exit  = new TH1D("MuTheta_all_exit",";Inital Muom Angle (#theta);Events",20,0,3.14);
  TH1D *muTheta_accepted_exit  = new TH1D("MuTheta_accepted_exit",";Initial Muon Angle (#theta);Events",20,0,3.14);

  //Right sign and wrong sign
  TH1D *RS_muE  = new TH1D("RS_muE",";Inital Muon Energy [GeV];Events",20,0,emax);
  TH1D *WS_muE  = new TH1D("WS_muE",";Inital Muon Energy [GeV];Events",20,0,emax);

  TH1D *RS_muE_cont  = new TH1D("RS_muE_cont",";Inital Muon Energy [GeV];Events",20,0,emax);
  TH1D *WS_muE_cont  = new TH1D("WS_muE_cont",";Inital Muon Energy [GeV];Events",20,0,emax);

  TH1D *RS_muE_exit  = new TH1D("RS_muE_exit",";Inital Muon Energy [GeV];Events",20,0,emax);
  TH1D *WS_muE_exit  = new TH1D("WS_muE_exit",";Inital Muon Energy [GeV];Events",20,0,emax);

  TH1D *RS_muL  = new TH1D("RS_muL",";Inital Muon Energy [GeV];Events",200, 0, 300);
  TH1D *WS_muL  = new TH1D("WS_muL",";Inital Muon Energy [GeV];Events",200, 0, 300);

  TH1D *RS_muL_cont  = new TH1D("RS_muL_cont",";Inital Muon Energy [GeV];Events",200, 0, 300);
  TH1D *WS_muL_cont  = new TH1D("WS_muL_cont",";Inital Muon Energy [GeV];Events",200, 0, 300);

  TH1D *RS_muL_exit  = new TH1D("RS_muL_exit",";Inital Muon Energy [GeV];Events",200, 0, 300);
  TH1D *WS_muL_exit  = new TH1D("WS_muL_exit",";Inital Muon Energy [GeV];Events",200, 0, 300);

  TH1D *numuCC_Osc[npoints+1];
  TH1D *numuCC_Osc_truth[npoints+1];
  TH1D *numuCC_Osc_cont[npoints+1];
  TH1D *numuCC_Osc_exit[npoints+1];

  TH1D *numuNC_Osc[npoints+1];
  TH1D *numuNC_Osc_cont[npoints+1];
  TH1D *numuNC_Osc_exit[npoints+1];
  TH1D *numuNC_Osc_truth[npoints+1];

  for(int i = 0; i <= npoints; i++){
    numuCC_Osc[i] = (TH1D*)numuCC->Clone();    
    numuCC_Osc_truth[i] = (TH1D*)numuCC_truth->Clone();
    numuCC_Osc_cont[i] = (TH1D*)numuCC_cont->Clone();
    numuCC_Osc_exit[i] = (TH1D*)numuCC_exit->Clone();
 
    numuNC_Osc[i] = (TH1D*)numuNC->Clone();
    numuNC_Osc_exit[i] = (TH1D*)numuNC_exit->Clone();
    numuNC_Osc_cont[i] = (TH1D*)numuNC_cont->Clone();
    numuNC_Osc_truth[i] = (TH1D*)numuNC_truth->Clone();

  }

  double cont_pion = 0, exit_pion = 0, NCevents = 0, CCevents = 0, NCselected = 0;

  double CC_wEpi = 0;
  double NC_wEpi = 0;

  //====================================================
  // Loop over entries in incoming ntuple
  //====================================================

  //Want to keep track of the L/E for the nu_mu oscillation studies
  int nu_LE_counter = 0;
  double nu_LE_wgt = 0;


  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  if( max_entry == -1 ) max_entry = nentries;

  for( Long64_t jentry=0; jentry<nentries && jentry<max_entry; jentry++ ){
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = GetEntry(jentry);   nbytes += nb;
    iTop = 0;
    iChan = 0;
    

    
    // Clear ntuple variables
    ConversionDist.clear();
    efficiency = 0.0;
    wgt = 0.0;
    electron_cand_energy = 0.0;
    electron_cand_angle = 0.0;
    enuccqe = 0.0;
    EnuReco = 0.;
    enucalo1 = 0.0;
    enucalo2 = 0.0;
    
    // Calculate the FS lepton's 4-vector and return neutrino path length
    nuleng = CalcLepton( detect_dist );

    // Get flux weight from FluxRW utilities
    
    fluxweight = utils.GetFluxWeight( energy, iflux, inno, ndecay );

    if(iDet == 3) {

      fluxweight *= pow(700,2)/pow(600,2);

      std::cout << "############################################################## Rescaling from 700m to 600m, USE ONLY IF RUNNING OVER 700m FILE!!!!! ###################" << std::endl;

    }


    // Lump in the Pot weight to flux weight so it isn't lost:
    fluxweight *= potweight;



    // Is the interaction in the fiducial volume?
    TVector3 vtx(Vx, Vy, Vz);
    isFid = utils.IsFiducial( iDet, ND_iDet, vtx );
    isActive = utils.IsActive( iDet, ND_iDet, vtx );

    if( !isActive ) continue;
    evt_wgt = true_nuleng = reco_energy = true_enugen = 0;
    
   
    //----------------------------------------------
    // Cycle through all pions in the event
    //----------------------------------------------
    /*    if (!isCC && ChargedPionMom->size() == 0) {
      continue;
    }
    
    if(signal == "numu" && isFid && !isCC && ChargedPionPos->size() > 0){
      
      nu_LE_wgt = 0; efficiency = muonIDeff; wgt = fluxweight*efficiency; nu_LE_wgt = (nuleng/100.)/(enugen*1000.);
   
      //      if(iDet == 0) wgt /= 3.0;

      NCevents += wgt;
      


      int pi_iter = -500, pi_exit = 0;
      ChargedPionExitPos *= 0.0;
      ChargedPionExitMom *= 0.0;
      ChargedPionInitialMom *= 0.0;

      int exit_pi_num = 0;

      bool Pion_contained = false;

      for(unsigned int pi = 0; pi < ChargedPionPos->size(); pi++){

        vector<gan::LorentzVectorLight> cPiPos = ChargedPionMom->at(pi);
	vector<gan::LorentzVectorLight> cPiMom = ChargedPionPos->at(pi);



        //Match the pion to the vertex

        unsigned int piPosIndex = 0;

        bool contained_pi = true;

        Pion_D->Fill(sqrt(pow((vtx.X()-cPiMom[0].X()),2)+
                          pow((vtx.Y()-cPiMom[0].Y()),2)+
                          pow((vtx.Z()-cPiMom[0].Z()),2)),wgt);

        if(sqrt(pow((vtx.X()-cPiMom[0].X()),2)+
		pow((vtx.Y()-cPiMom[0].Y()),2)+
                pow((vtx.Z()-cPiMom[0].Z()),2)) > 10) continue;

        while(contained_pi && piPosIndex < cPiMom.size()){

          TVector3 pipos(cPiMom[piPosIndex].X(),
                         cPiMom[piPosIndex].Y(),
                         cPiMom[piPosIndex].Z());
          contained_pi = utils.IsActive(iDet, ND_iDet,pipos,5);
          piPosIndex++;

        }

        if(!contained_pi) {
          exit_pi_num += 1; pi_iter = pi; pi_exit = piPosIndex;
	}//exiting pions

	if(contained_pi && exit_pi_num == 0){

          pi_iter = pi; pi_exit = piPosIndex;

	}//contained pions
      }//scan over pions and see if they are contained

      N_Exit_Pion->Fill(exit_pi_num,wgt);

      bool contained = false;
      if(pi_iter == -500) continue;
      if(exit_pi_num > 1){
        NC_wEpi += wgt;}
      //continue;}
      if(exit_pi_num == 1 ){
	exit_pion += wgt;
      }

      if(exit_pi_num == 0) { wgt *= piEff_pi; contained = true;}




      NCselected += wgt;

      ChargedPionExitPos = TVector3(ChargedPionPos->at(pi_iter)[pi_exit-1].X(),
                                    ChargedPionPos->at(pi_iter)[pi_exit-1].Y(),
                                    ChargedPionPos->at(pi_iter)[pi_exit-1].Z());
      ChargedPionExitMom = TVector3(ChargedPionMom->at(pi_iter)[pi_exit-1].X(),
                                    ChargedPionMom->at(pi_iter)[pi_exit-1].Y(),
                                    ChargedPionMom->at(pi_iter)[pi_exit-1].Z());
      ChargedPionInitialMom = TVector3(ChargedPionMom->at(pi_iter)[0].X(),
                                       ChargedPionMom->at(pi_iter)[0].Y(),
                                       ChargedPionMom->at(pi_iter)[0].Z());


      Double_t photon_energy = utils.TotalPhotonEnergy( iDet, ND_iDet, p1PhotonConversionPos, p1PhotonConversionMom,
                                                        p2PhotonConversionPos, p2PhotonConversionMom );


      //---------------------------------------------
      //After this point we don't know they are pions
      //---------------------------------------------

      //Unsmeared Lepton Energy
      mu_E = ChargedPionMom->at(pi_iter)[0].E();


      //internal muon length
      mu_L = 0;
      mu_L = sqrt(pow((vtx.X()-ChargedPionExitPos.X()),2)+
                  pow((vtx.Y()-ChargedPionExitPos.Y()),2)+
                  pow((vtx.Z()-ChargedPionExitPos.Z()),2));

      Pion_L->Fill(mu_L,wgt);

      //      std::cout << "Pion track length : " << mu_L << std::endl;

      //Initial Lepton production angle
      mu_Theta = utils.GetTheta(ChargedPionInitialMom.X(), ChargedPionInitialMom.Y(), ChargedPionInitialMom.Z());

      int HowMeasuresed = 0;

      //Smear the lepton energy and return the smearing for Syst studies
      ranger = utils.RangeOut(ChargedPionExitPos, ChargedPionExitMom, mu_L, ND_iDet);

      smeared_mu_E = utils.EnergyRes(contained, mu_E, mu_L, ND_iDet, ranger, MuonExitMom,
                                     MuonInitialMom, MuonExitPos, vtx, 0, HowMeasuresed);

      if(smeared_mu_E < 0) smeared_mu_E = 0;

      //Compute neutrino reco. energy
      if(smeared_mu_E > 0){
	EnuReco = utils.NuEnergyCalo_smeared( GeniePDG,
                                              GenieE,
                                              0.02, //proton threshold
                                              smeared_mu_E,
                                              mu_E,
                                              "NC",
                                              0.05)+photon_energy; //hadronic resolution
      }
      else{EnuReco = 0;}

      if(smeared_mu_E > 0){
	numuNC_truth->Fill(enugen, wgt);
        numuNC->Fill(EnuReco, wgt);
        numuNCLepE->Fill(smeared_mu_E, wgt);
	numuNCLepLength->Fill(mu_L, wgt);
	numuNCLepZLength->Fill(mu_L*(TMath::Cos(mu_Theta)), wgt);
	numuNCLepRLength->Fill(mu_L*(TMath::Sin(mu_Theta)), wgt);
	numuNCLepTheta->Fill(mu_Theta*(180/3.14), wgt);
	
        if(!contained){
	  numuNC_exit->Fill(EnuReco,wgt);
	  numuNCLepE_exit->Fill(smeared_mu_E, wgt);
	  numuNCLepLength_exit->Fill(mu_L, wgt);
	  numuNCLepZLength_exit->Fill(mu_L*(TMath::Cos(mu_Theta)), wgt);
	  numuNCLepRLength_exit->Fill(mu_L*(TMath::Sin(mu_Theta)), wgt);
	  numuNCLepTheta_exit->Fill(mu_Theta*(180/3.14), wgt);}

	else{
	  numuNC_cont->Fill(EnuReco,wgt);
	  numuNCLepE_cont->Fill(smeared_mu_E, wgt);
	  numuNCLepLength_cont->Fill(mu_L, wgt);
	  numuNCLepZLength_cont->Fill(mu_L*(TMath::Cos(mu_Theta)), wgt);
	  numuNCLepRLength_cont->Fill(mu_L*(TMath::Sin(mu_Theta)), wgt);
	  numuNCLepTheta_cont->Fill(mu_Theta*(180/3.14), wgt);
	}


        //Partial oscillation computation
	for(int dm2 = 0; dm2 <= npoints; dm2++){

          DMpoint = pow(10., (TMath::Log10(dm2min)+ (dm2 * (1./npoints))*TMath::Log10(dm2max/dm2min)) );
          numuNC_Osc[dm2]->Fill(EnuReco,wgt*pow(TMath::Sin((1.267)*(DMpoint)*nu_LE_wgt),2));

          if(!contained){numuNC_Osc_exit[dm2]->Fill(EnuReco,wgt*pow(TMath::Sin((1.267)*(DMpoint)*nu_LE_wgt),2));}
          else{numuNC_Osc_cont[dm2]->Fill(EnuReco,wgt*pow(TMath::Sin((1.267)*(DMpoint)*nu_LE_wgt),2));}


        }
      }// smeared E > 0!


    }//select NC event with some number of pions

    */
   
    //----------------------------------------------
    //
    // NOW MUONS!
    //
    //----------------------------------------------

    //----------------------------------------------
    // Track the muon till it exits the detector:
    //----------------------------------------------
      if (isCC && MuonMom->size() == 0) {
	  continue;
      }

    bool contained = true;
    unsigned int muonPosIndex = 0;
    MuonExitPos *= 0.0;
    MuonExitMom *= 0.0;
    MuonInitialMom *= 0.0;

    if(isCC && abs(inno) == 14){

      MC_mup->Fill(sqrt(pow(MuonMom->at(0).X(),2)+pow(MuonMom->at(0).Y(),2)+pow(MuonMom->at(0).Z(),2)),fluxweight);
      MC_mur->Fill(sqrt(pow(MuonMom->at(0).X(),2)+pow(MuonMom->at(0).Y(),2)+pow(MuonMom->at(0).Z(),2))/0.00212,fluxweight);
      MC_muT->Fill(double(180.0/3.14)*(utils.GetTheta(MuonMom->at(0).X(), MuonMom->at(0).Y(), MuonMom->at(0).Z())),fluxweight);
      MC_mucosT->Fill(fabs(TMath::Cos(utils.GetTheta(MuonMom->at(0).X(), MuonMom->at(0).Y(), MuonMom->at(0).Z()))),fluxweight);

      MC_mupZ->Fill(MuonMom->at(0).Z(),fluxweight);
      MC_murZ->Fill((MuonMom->at(0).Z())/0.00212,fluxweight);

      MC_mupT->Fill(sqrt(pow(MuonMom->at(0).X(),2)+pow(MuonMom->at(0).Y(),2)),fluxweight);
      MC_murT->Fill((sqrt(pow(MuonMom->at(0).X(),2)+pow(MuonMom->at(0).Y(),2)))/0.00212,fluxweight);




      while ( contained && muonPosIndex < MuonPos->size() ){
        TVector3 pos(MuonPos->at(muonPosIndex).X(), MuonPos->at(muonPosIndex).Y(), MuonPos->at(muonPosIndex).Z());
        contained = utils.IsActive(iDet, ND_iDet,pos,5);
        muonPosIndex++;
      }

      MuonExitPos = TVector3(MuonPos->at(muonPosIndex-1).X(), MuonPos->at(muonPosIndex-1).Y(), MuonPos->at(muonPosIndex-1).Z());
      MuonExitMom = TVector3(MuonMom->at(muonPosIndex-1).X(), MuonMom->at(muonPosIndex-1).Y(), MuonMom->at(muonPosIndex-1).Z());
      MuonInitialMom = TVector3(MuonMom->at(0).X(), MuonMom->at(0).Y(), MuonMom->at(0).Z());
    }
    

    unsigned int piPosIndex = 0;
     

    if(signal == "numu" && isFid && isCC && abs(inno) == 14 && MuonMom->at(0).E() > egammaThreshold){

      //weights

      nu_LE_wgt = 0; efficiency = muonIDeff; wgt = fluxweight*efficiency; nu_LE_wgt = (nuleng/100.)/(enugen*1000.);


      //      if(iDet == 0) wgt /= 3.0;
      if(contained) wgt *= muEff_pi;

      CCevents += wgt;

      Double_t photon_energy = utils.TotalPhotonEnergy( iDet, ND_iDet, p1PhotonConversionPos, p1PhotonConversionMom, 
							p2PhotonConversionPos, p2PhotonConversionMom );


      Muon_D->Fill(sqrt(pow((vtx.X()-MuonPos->at(0).X()),2)+
			pow((vtx.Y()-MuonPos->at(0).Y()),2)+
			pow((vtx.Z()-MuonPos->at(0).Z()),2)),wgt);

      int HowMeasuresed = 0;
      //Unsmeared Lepton Energy
      mu_E = MuonMom->at(0).E();
   
      //internal muon length
      mu_L = 0;
      mu_L = sqrt(pow((vtx.X()-MuonExitPos.X()),2)+pow((vtx.Y()-MuonExitPos.Y()),2)+pow((vtx.Z()-MuonExitPos.Z()),2));

      //Initial Lepton production angle 
      mu_Theta = utils.GetTheta(MuonMom->at(0).X(), MuonMom->at(0).Y(), MuonMom->at(0).Z());

      //Smear the lepton energy and return the smearing for Syst studies
      ranger = utils.RangeOut(MuonExitPos, MuonExitMom, mu_L, ND_iDet);
      smeared_mu_E = utils.EnergyRes(contained, mu_E, mu_L, ND_iDet, ranger, MuonExitMom, MuonInitialMom, MuonExitPos, vtx, 0, HowMeasuresed);

      if(smeared_mu_E < 0) smeared_mu_E = 0;

      if(smeared_mu_E > 0){
	EnuReco = utils.NuEnergyCalo_smeared( GeniePDG,
					      GenieE,
					      0.02,
					      smeared_mu_E,
					      mu_E,
					      "muon",
					      0.05)+photon_energy;
      }
      else{EnuReco = 0;}
      
      nuE_all->Fill(enugen, wgt);
      muE_all->Fill(mu_E, wgt);
      muTheta_all->Fill(mu_Theta, wgt);
      muTE_all->Fill(mu_Theta,mu_E,wgt);
      muLx_all->Fill(vtx.X(),wgt);
      muLy_all->Fill(vtx.Y(),wgt);
      muLz_all->Fill(vtx.Z(),wgt);
      nuEvmuE_Corr_MC->Fill(enugen,mu_E,wgt);

      if(contained){ nuEvmuE_Corr_MC_cont->Fill(enugen,mu_E,wgt);
	nuE_all_cont->Fill(enugen, wgt);
	muE_all_cont->Fill(mu_E, wgt);
	muTheta_all_cont->Fill(mu_Theta, wgt);
      }

      if(!contained){  nuEvmuE_Corr_MC_exit->Fill(enugen,mu_E,wgt);
	nuE_all_exit->Fill(enugen, wgt);
        muE_all_exit->Fill(mu_E, wgt);
	muTheta_all_exit->Fill(mu_Theta, wgt);}

      
      if(smeared_mu_E > 0){
	nuE_accepted->Fill(enugen, wgt);
	muE_accepted->Fill(mu_E, wgt);
	muTheta_accepted->Fill(mu_Theta, wgt);
	muTE_accepted->Fill(mu_Theta,mu_E,wgt);
	muLx_accepted->Fill(vtx.X(),wgt);
	muLy_accepted->Fill(vtx.Y(),wgt);
	muLz_accepted->Fill(vtx.Z(),wgt);
	if(contained){
	  nuE_accepted_cont->Fill(enugen, wgt);
	  muE_accepted_cont->Fill(mu_E, wgt);
	  muTheta_accepted_cont->Fill(mu_Theta, wgt);
	}
        if(!contained){
          nuE_accepted_exit->Fill(enugen, wgt);
          muE_accepted_exit->Fill(mu_E, wgt);
          muTheta_accepted_exit->Fill(mu_Theta, wgt);
	}

      }


      if(HowMeasuresed == 1) Cont += wgt;
      else if(HowMeasuresed == 2) MS += wgt;
      else if(HowMeasuresed == 3) LArB += wgt;
      else if(HowMeasuresed == 4) RangerB += wgt;
      else if(smeared_mu_E == 0) Unmeas += wgt;
      
      //Determine number of muons which exit where
      if(smeared_mu_E > 0){
	if(contained){ stopped += wgt;}
	else if(MuonExitPos.X() <= nd_xmin+5){Xmin += wgt;}
	else if(MuonExitPos.X() >= nd_xmax-5){Xmax += wgt;}
	else if(MuonExitPos.Y() <= nd_ymin+5){Ymin += wgt;}
	else if(MuonExitPos.Y() >= nd_ymax-5){Ymax += wgt;}
	else if(MuonExitPos.Z() <= nd_zmin+5){Zmin += wgt;}
	else if(MuonExitPos.Z() >= nd_zmax-5){Zmax += wgt;}
      }
      if(smeared_mu_E == 0){
	if(MuonExitPos.X() <= nd_xmin+5){Xmin_F += wgt;}
	else if(MuonExitPos.X() >= nd_xmax-5){Xmax_F += wgt;}
        else if(MuonExitPos.Y() <= nd_ymin+5){Ymin_F += wgt;}
	else if(MuonExitPos.Y() >= nd_ymax-5){Ymax_F += wgt;}
	else if(MuonExitPos.Z() <= nd_zmin+5){Zmin_F += wgt;}
        else if(MuonExitPos.Z() >= nd_zmax-5){Zmax_F += wgt;}
      }
      //Filling histograms
      //Wrong Sign
      if(smeared_mu_E > 0 && inno < 0){
	
	WS_muE->Fill(smeared_mu_E, wgt);
	WS_muL->Fill(mu_L, wgt);

	if(contained){
	  WS_muE_cont->Fill(smeared_mu_E, wgt);
	  WS_muL_cont->Fill(mu_L, wgt);
	}

	if(!contained){
	  WS_muE_exit->Fill(smeared_mu_E, wgt);
	  WS_muL_exit->Fill(mu_L, wgt);
	}
      }

      if(smeared_mu_E > 0 && inno > 0){

	RS_muE->Fill(smeared_mu_E, wgt);
	RS_muL->Fill(mu_L, wgt);

	if(contained){
	  RS_muE_cont->Fill(smeared_mu_E, wgt);
	  RS_muL_cont->Fill(mu_L, wgt);
	}

        if(!contained){
          RS_muE_exit->Fill(smeared_mu_E, wgt);
          RS_muL_exit->Fill(mu_L, wgt);
	}
      

	numuCC_truth->Fill(enugen, wgt);
	muE_Reco_True->Fill(smeared_mu_E, mu_E, wgt);
	nuE_Reco_True->Fill(EnuReco,enugen,wgt);
	
	numuCC->Fill(EnuReco,wgt);
	nuParP_T_Z->Fill(sqrt(pow(ParPx,2)+pow(ParPy,2)),ParPz,wgt);
	numuCCLepE->Fill(smeared_mu_E, wgt);
	numuCCLepTheta->Fill(mu_Theta*(180/3.14),wgt);
	numuCCLepLength->Fill(mu_L,wgt);
	numuCCLepZLength->Fill(mu_L*(TMath::Cos(mu_Theta)),wgt);
	numuCCLepRLength->Fill(mu_L*(TMath::Sin(mu_Theta)),wgt);


	
	muEvL_Corr->Fill(smeared_mu_E, mu_L, wgt);

	if(contained){ 
	  muEvL_Corr_cont->Fill(smeared_mu_E, mu_L, wgt);
	  nuParP_T_Z_cont->Fill(sqrt(pow(ParPx,2)+pow(ParPy,2)),ParPz,wgt);}
	if(!contained){
	  muEvL_Corr_exit->Fill(smeared_mu_E, mu_L, wgt);
	  nuParP_T_Z_exit->Fill(sqrt(pow(ParPx,2)+pow(ParPy,2)),ParPz,wgt);}


	nuEvmuE_Corr->Fill(enugen,mu_E,wgt);

	if(contained){ 
	  nuEvmuE_Corr_cont->Fill(enugen,mu_E,wgt);
	}
	
	if(!contained){ 
	  nuEvmuE_Corr_exit->Fill(enugen,mu_E,wgt);
	}


	double res2 = fabs(mu_E-smeared_mu_E)/(mu_E);
	if(smeared_mu_E != 0){ 
	  Res_all->Fill(mu_L, res2);
	  Res_all_N->Fill(mu_L);
	}

	if(contained){
	  numuCC_cont->Fill(EnuReco,wgt); 
	  numuCCLepE_cont->Fill(smeared_mu_E, wgt);
	  numuCCLepTheta_cont->Fill(mu_Theta,wgt);
	  numuCCLepLength_cont->Fill(mu_L,wgt);
	  numuCCLepZLength_cont->Fill(mu_L*(TMath::Cos(mu_Theta)),wgt);
	  numuCCLepRLength_cont->Fill(mu_L*(TMath::Sin(mu_Theta)),wgt);

	  Res_cont->Fill(mu_L,res2);
	  Res_cont_N->Fill(mu_L); 
	  
	  numuCC_cont_truth->Fill(enugen, wgt);
	  
	}
	
	if(!contained){
	  numuCC_exit->Fill(EnuReco,wgt);
	  numuCCLepE_exit->Fill(smeared_mu_E, wgt);
	  numuCCLepTheta_exit->Fill(mu_Theta,wgt);
          numuCCLepLength_exit->Fill(mu_L,wgt);
          numuCCLepZLength_exit->Fill(mu_L*(TMath::Cos(mu_Theta)),wgt);
          numuCCLepRLength_exit->Fill(mu_L*(TMath::Sin(mu_Theta)),wgt);


	  Res_exit->Fill(mu_L,res2);
	  Res_exit_N->Fill(mu_L);	  

	  numuCC_exit_truth->Fill(enugen, wgt);

	}  
      
	//Partial oscillation computation
	for(int dm2 = 0; dm2 <= npoints; dm2++){
	  
	  DMpoint = pow(10., (TMath::Log10(dm2min)+ (dm2 * (1./npoints))*TMath::Log10(dm2max/dm2min)) );
	  numuCC_Osc[dm2]->Fill(EnuReco,wgt*pow(TMath::Sin((1.267)*(DMpoint)*nu_LE_wgt),2));
	  numuCC_Osc_truth[dm2]->Fill(enugen,wgt*pow(TMath::Sin((1.267)*(DMpoint)*nu_LE_wgt),2));
	  if(!contained){
	    numuCC_Osc_exit[dm2]->Fill(EnuReco,wgt*pow(TMath::Sin((1.267)*(DMpoint)*nu_LE_wgt),2));
	  }
	  else{numuCC_Osc_cont[dm2]->Fill(EnuReco,wgt*pow(TMath::Sin((1.267)*(DMpoint)*nu_LE_wgt),2));}

	}
      }// smeared E > 0!
    }
    

    if ( (jentry%100000) == 0){ std::cout << int(100*jentry/nentries) << "% done !"  << std::endl;}

    
  } // end loop over entries 


  for( int L = 0; L <= Res_all->GetNbinsX(); L++){

    double bin = 0, bin_cont = 0, bin_exit = 0; 
    double binN = 0, binN_cont = 0, binN_exit = 0;

    bin = Res_all->GetBinContent(L);
    bin_cont = Res_cont->GetBinContent(L);
    bin_exit = Res_exit->GetBinContent(L);

    binN = Res_all_N->GetBinContent(L);
    binN_cont = Res_cont_N->GetBinContent(L);
    binN_exit = Res_exit_N->GetBinContent(L);

    if(binN_exit == 0 || bin_exit == 0){bin_exit = 0; binN_exit = 1;}

    Res_all->SetBinContent(L,(bin)/(binN));
    Res_cont->SetBinContent(L,(bin_cont)/(binN_cont));
    Res_exit->SetBinContent(L,(bin_exit)/(binN_exit));

    std::cout << "Average Exit Res : " << Res_exit->GetBinContent(L) << std::endl;

  }

  //Build an average L/E plot in bins of reconstructed energy so I can build osc. probs on the fly...Knock on wood

  TString outfile = InFile().Remove(InFile().Length()-5) + "_processed_CoreyBins_" + signal + "_Joseph_Smeared.root";
  std::cout << "outfile file: " << outfile << std::endl;
  TFile *f = new TFile(outfile, "RECREATE");
  if (f->IsZombie()){
    std::cout << "Error, couldn't create output file." << std::endl;
    return;
  }  
  //-- write histograms
  MC_mup->Write();
  MC_mur->Write();
  MC_muT->Write();
  MC_mucosT->Write();

  MC_mupZ->Write();
  MC_mupT->Write();

  MC_murZ->Write();
  MC_murT->Write();



  N_Exit_Pion->Write();
  N_Exit_Pion_CC->Write();
  N_Cont_Pion_CC->Write();
  Pion_cont_L->Write();
  Pion_L->Write();
  Muon_D->Write();
  Pion_D->Write();

  numuCC->Write();
  numuCC_cont->Write();
  numuCC_exit->Write();

  nuParP_T_Z->Write();
  nuParP_T_Z_cont->Write();
  nuParP_T_Z_exit->Write();

  numuCCLepE->Write();
  numuCCLepE_cont->Write();
  numuCCLepE_exit->Write();

  numuCCLepTheta->Write();
  numuCCLepTheta_cont->Write();
  numuCCLepTheta_exit->Write();

  nuE_all->Write();
  muE_all->Write();
  muTheta_all->Write();
  muTE_all->Write();

  nuE_accepted->Write();
  muE_accepted->Write();
  muTheta_accepted->Write();
  muTE_accepted->Write();

  numuNC_truth->Write();
  numuCC_truth->Write();
  numuCC_exit_truth->Write();

  muE_Reco_True->Write();
  nuE_Reco_True->Write();

  muLx_all->Write();
  muLx_accepted->Write();
  muLy_all->Write();
  muLy_accepted->Write();
  muLz_all->Write();
  muLz_accepted->Write();

  numuCC_cont_truth->Write();

  Res_all->Write();
  Res_cont->Write();
  Res_exit->Write();

  muEvL_Corr->Write();
  muEvL_Corr_cont->Write();
  muEvL_Corr_exit->Write();

  nuEvmuE_Corr->Write();
  nuEvmuE_Corr_cont->Write();
  nuEvmuE_Corr_exit->Write();

  nuEvmuE_Corr_MC->Write();
  nuEvmuE_Corr_MC_cont->Write();
  nuEvmuE_Corr_MC_exit->Write();

  RS_muE->Write();
  WS_muE->Write();
  RS_muE_cont->Write();
  WS_muE_cont->Write();
  RS_muE_exit->Write();
  WS_muE_exit->Write();

  RS_muL->Write();
  WS_muL->Write();
  RS_muL_cont->Write();
  WS_muL_cont->Write();
  RS_muL_exit->Write();
  WS_muL_exit->Write();

  nuE_all_cont->Write();
  nuE_accepted_cont->Write();
  muE_all_cont->Write();
  muE_accepted_cont->Write();
  muTheta_all_cont->Write();
  muTheta_accepted_cont->Write();
  nuE_all_exit->Write();
  nuE_accepted_exit->Write();
  muE_all_exit->Write();
  muE_accepted_exit->Write();
  muTheta_all_exit->Write();
  muTheta_accepted_exit->Write();

  numuNC->Write();
  numuNC_cont->Write();
  numuNC_exit->Write();


  numuCCLepTheta->Write();
  numuCCLepTheta_cont->Write();
  numuCCLepTheta_exit->Write();

  numuCCLepLength->Write();
  numuCCLepLength_cont->Write();
  numuCCLepLength_exit->Write();

  numuCCLepLength->Write();
  numuCCLepLength_cont->Write();
  numuCCLepLength_exit->Write();

  numuCCLepZLength->Write();
  numuCCLepZLength_cont->Write();
  numuCCLepZLength_exit->Write();

  numuCCLepRLength->Write();
  numuCCLepRLength_cont->Write();
  numuCCLepRLength_exit->Write();

  numuNCLepTheta->Write();
  numuNCLepTheta_cont->Write();
  numuNCLepTheta_exit->Write();

  numuNCLepLength->Write();
  numuNCLepLength_cont->Write();
  numuNCLepLength_exit->Write();

  numuNCLepLength->Write();
  numuNCLepLength_cont->Write();
  numuNCLepLength_exit->Write();

  numuNCLepZLength->Write();
  numuNCLepZLength_cont->Write();
  numuNCLepZLength_exit->Write();

  numuNCLepRLength->Write();
  numuNCLepRLength_cont->Write();
  numuNCLepRLength_exit->Write();

  for(int dm2 = 0; dm2 <= npoints; dm2++){
    std::string dmpoint = std::to_string(dm2);
    std::string name = "Osc_";
    std::string true_name = "Osc_True";
    std::string exit_name = "Osc_exit";
    std::string cont_name = "Osc_cont";
    std::string nc_name = "Osc_NC";
    std::string nc_exit_name = "Osc_NC_exit";
    std::string nc_cont_name = "Osc_NC_cont";

    name = name+dmpoint;
    true_name = true_name+dmpoint;
    exit_name = exit_name+dmpoint;
    cont_name = cont_name+dmpoint;
    nc_name = nc_name+dmpoint;
    nc_exit_name = nc_exit_name+dmpoint;
    nc_cont_name = nc_cont_name+dmpoint;

    numuCC_Osc[dm2]->Write(name.c_str());
    numuCC_Osc_truth[dm2]->Write(true_name.c_str());
    numuCC_Osc_exit[dm2]->Write(exit_name.c_str());
    numuCC_Osc_cont[dm2]->Write(cont_name.c_str());

    numuNC_Osc[dm2]->Write(nc_name.c_str());
    numuNC_Osc_exit[dm2]->Write(nc_exit_name.c_str());
    numuNC_Osc_cont[dm2]->Write(nc_cont_name.c_str());
  }

  f->Write();
  f->Close();

  std::cout << "==================================================" << std::endl;
  std::cout << "Rates for detector " << iDet << std::endl;
  std::cout << "Active volume = " << utils.GetActiveMass( iDet, ND_iDet ) << " tons" << std::endl;
  std::cout << "Fiducial volume = " << utils.GetFidMass( iDet, ND_iDet ) << " tons" << std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << "Analyzing " << signal << " events" << std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << "Of " << Xmin+Xmax+Ymin+Ymax+Zmin+Zmax+stopped << " events" << std::endl;
  std::cout << "Exited at X (%) " << (Xmin+Xmax)/(Xmin+Xmax+Ymin+Ymax+Zmin+Zmax+stopped) << std::endl; 
  std::cout << "Exited at Y (%) " << (Ymin+Ymax)/(Xmin+Xmax+Ymin+Ymax+Zmin+Zmax+stopped) <<std::endl;
  std::cout << "Exited at Zmin " << Zmin <<std::endl;
  std::cout << "Exited at Zmax " << Zmax <<std::endl;
  std::cout << "Exited at Zmin (%) " << Zmin/(Xmin+Xmax+Ymin+Ymax+Zmin+Zmax+stopped) <<std::endl;
  std::cout << "Exited at Zmax (%) " << Zmax/(Xmin+Xmax+Ymin+Ymax+Zmin+Zmax+stopped) <<std::endl;
  std::cout << "Contained muons (%) " << stopped/(Xmin+Xmax+Ymin+Ymax+Zmin+Zmax+stopped) <<std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << "Measured by : " <<std::endl;
  std::cout << "Contained " << Cont/(Unmeas+Cont+MS+LArB+RangerB) <<std::endl;
  std::cout << "MS " << MS/(Unmeas+Cont+MS+LArB+RangerB) <<std::endl;
  std::cout << "LAr B " << LArB/(Unmeas+Cont+MS+LArB+RangerB) <<std::endl;
  std::cout << "Ranger " << RangerB/(Unmeas+Cont+MS+LArB+RangerB) <<std::endl;
  std::cout << "Unmeasured " << Unmeas/(Unmeas+Cont+MS+LArB+RangerB) <<std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << "Of those that are Unmeasured :" << std::endl;
  std::cout << " Exit at X (%) " << (Xmin_F+Xmax_F)/(Xmin_F+Xmax_F+Ymin_F+Ymax_F+Zmin_F+Zmax_F) << std::endl;
  std::cout << " Exit at Y (%) " << (Ymin_F+Ymax_F)/(Xmin_F+Xmax_F+Ymin_F+Ymax_F+Zmin_F+Zmax_F)<< std::endl;
  std::cout << " Exit at Z max (%) " << (Zmax_F)/(Xmin_F+Xmax_F+Ymin_F+Ymax_F+Zmin_F+Zmax_F)<< std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << "Of all muons : " << std::endl;
  std::cout << "Contained " << stopped/(Xmin_F+Xmax_F+Ymin_F+Ymax_F+Zmin_F+Zmax_F+Xmin+Xmax+Ymin+Ymax+Zmin+Zmax+stopped) <<std::endl;
  std::cout << ">1 m, exit side " << (Xmin+Xmax+Ymin+Ymax)/(Xmin_F+Xmax_F+Ymin_F+Ymax_F+Zmin_F+Zmax_F+Xmin+Xmax+Ymin+Ymax+Zmin+Zmax+stopped) <<std::endl;
  std::cout << ">1 m, exit back " << (Zmax)/(Xmin_F+Xmax_F+Ymin_F+Ymax_F+Zmin_F+Zmax_F+Xmin+Xmax+Ymin+Ymax+Zmin+Zmax+stopped) <<std::endl;
  std::cout << "<1 m, exit side " << (Xmin_F+Xmax_F+Ymin_F+Ymax_F)/(Xmin_F+Xmax_F+Ymin_F+Ymax_F+Zmin_F+Zmax_F+Xmin+Xmax+Ymin+Ymax+Zmin+Zmax+stopped) <<std::endl;
  std::cout << "<1 m, exit back " << (Zmax_F)/(Xmin_F+Xmax_F+Ymin_F+Ymax_F+Zmin_F+Zmax_F+Xmin+Xmax+Ymin+Ymax+Zmin+Zmax+stopped) <<std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << "Number of total CC events : " << CCevents << std::endl;
  std::cout << "Number of CC events with exiting Pion and Muon : " << CC_wEpi << std::endl;
  std::cout << "==================================================" << std::endl;
  std::cout << "Number of CC events per spill: " << CCevents*(4.5e12/2.2e20) << std::endl;
  std::cout << "Output in " << outfile << std::endl;

}




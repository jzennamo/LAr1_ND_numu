//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 19 10:41:44 2014 by ROOT version 5.34/18
// from TTree EventsTot/Event info for ALL types
// found on file: /uboone/data/lar1/Ntuples/outputs_100m_nu/combined_ntuple_100m_nu.root
//////////////////////////////////////////////////////////

#ifndef NtupleReprocessing_h
#define NtupleReprocessing_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include "LorentzVectorLight.h"
#include <iostream>
#include <iomanip>
#include "Utils.cxx"

// Fixed size dimensions of array or collections stored in the TTree if any.

class NtupleReprocessing {
public :

   TFile          *infile;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           iflux;
   Int_t           ibkg;
   Int_t           nuchan;
   Int_t           inno;
   Double_t        enugen;
   Double_t        energy;
   Double_t        nuleng;
   Int_t           parid;
   Int_t           ndecay;
   Double_t        wgt;
   Int_t           NPi0;
   Int_t           NPi0FinalState;
   Int_t           NGamma;
   Char_t          FoundPhotons;
   vector<int>     *GeniePDG;
   vector<double>  *GenieE;
   vector<double>  *GeniePx;
   vector<double>  *GeniePy;
   vector<double>  *GeniePz;
   vector<string>  *GenieProc;
   Int_t           isCC;
   Int_t           mode;
   Double_t        Vx;
   Double_t        Vy;
   Double_t        Vz;
   Double_t        ParVx;
   Double_t        ParVy;
   Double_t        ParVz;
   Double_t        ParPx;
   Double_t        ParPy;
   Double_t        ParPz;
   Double_t        LepPx;
   Double_t        LepPy;
   Double_t        LepPz;
   Double_t        pdpx;
   Double_t        pdpy;
   Double_t        pdpz;
   Double_t        pppx;
   Double_t        pppy;
   Double_t        pppz;
   Double_t        tpx;
   Double_t        tpy;
   //  Double_t        tpx;
   Int_t           ptype;
   Int_t           tptype;
   gan::LorentzVectorLight *neutMom;
   Double_t        ThetaLep;
   Double_t        PhiLep;
   Double_t        ThetaLepSmeared;
   Double_t        PhiLepSmeared;
   Double_t        Elep;
   Double_t        ElepSmeared;
   vector<gan::LorentzVectorLight> *MuonPos;
   vector<gan::LorentzVectorLight> *MuonMom;
   vector<gan::LorentzVectorLight> *ElectronPos;
   vector<gan::LorentzVectorLight> *ElectronMom;
   vector<gan::LorentzVectorLight> *p1PhotonConversionPos;
   vector<gan::LorentzVectorLight> *p1PhotonConversionMom;
   vector<gan::LorentzVectorLight> *p2PhotonConversionPos;
   vector<gan::LorentzVectorLight> *p2PhotonConversionMom;
   vector<gan::LorentzVectorLight> *miscPhotonConversionPos;
   vector<gan::LorentzVectorLight> *miscPhotonConversionMom;
   vector<gan::LorentzVectorLight> *PionPos;
   vector<gan::LorentzVectorLight> *PionMom;
   vector<vector<gan::LorentzVectorLight> > *ChargedPionPos;
   vector<vector<gan::LorentzVectorLight> > *ChargedPionMom;
   vector< int>     *ChargePionSign;
   vector< vector<float> > *MultiWeight;

   // List of branches
   TBranch        *b_iflux;   //!
   TBranch        *b_ibkg;   //!
   TBranch        *b_nuchan;   //!
   TBranch        *b_inno;   //!
   TBranch        *b_enugen;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_nuleng;   //!
   TBranch        *b_parid;   //!
   TBranch        *b_ndecay;   //!
   TBranch        *b_wgt;   //!
   TBranch        *b_NPi0;   //!
   TBranch        *b_NPi0FinalState;   //!
   TBranch        *b_NGamma;   //!
   TBranch        *b_FoundAllPhotons;   //!
   TBranch        *b_GeniePDG;   //!
   TBranch        *b_GenieE;   //!
   TBranch        *b_GeniePx;   //!
   TBranch        *b_GeniePy;   //!
   TBranch        *b_GeniePz;   //!
   TBranch        *b_GenieProc;   //!
   TBranch        *b_isCC;   //!
   TBranch        *b_mode;   //!
   TBranch        *b_Vx;   //!
   TBranch        *b_Vy;   //!
   TBranch        *b_Vz;   //!
   TBranch        *b_ParVx;   //!
   TBranch        *b_ParVy;   //!
   TBranch        *b_ParVz;   //!
   TBranch        *b_ParPx;   //!
   TBranch        *b_ParPy;   //!
   TBranch        *b_ParPz;   //!
   TBranch        *b_LepPx;   //!
   TBranch        *b_LepPy;   //!
   TBranch        *b_LepPz;   //!
   TBranch        *b_pdpx;   //!
   TBranch        *b_pdpy;   //!
   TBranch        *b_pdpz;   //!
   TBranch        *b_pppx;   //!
   TBranch        *b_pppy;   //!
   TBranch        *b_pppz;   //!
   TBranch        *b_tpx;   //!
   TBranch        *b_tpy;   //!
   // TBranch        *b_tpx;   //!
   TBranch        *b_ptype;   //!
   TBranch        *b_tptype;   //!
   TBranch        *b_neutMom;   //!
   TBranch        *b_ThetaLep;   //!
   TBranch        *b_PhiLep;   //!
   TBranch        *b_ThetaLepSmeared;   //!
   TBranch        *b_PhiLepSmeared;   //!
   TBranch        *b_Elep;   //!
   TBranch        *b_ElepSmeared;   //!
   TBranch        *b_MuonPos;   //!
   TBranch        *b_MuonMom;   //!
   TBranch        *b_ElectronPos;   //!
   TBranch        *b_ElectronMom;   //!
   TBranch        *b_p1PhotonConversionPos;   //!
   TBranch        *b_p1PhotonConversionMom;   //!
   TBranch        *b_p2PhotonConversionPos;   //!
   TBranch        *b_p2PhotonConversionMom;   //!
   TBranch        *b_miscPhotonConversionPos;   //!
   TBranch        *b_miscPhotonConversionMom;   //!
   TBranch        *b_PionPos;   //!
   TBranch        *b_PionMom;   //!
   TBranch        *b_ChargedPionPos;   //!
   TBranch        *b_ChargedPionMom;   //!
   TBranch        *b_ChargePionSign;   //!
   TBranch        *b_MultiWeight;   //!

   NtupleReprocessing(TString file = "");
   //   NtupleReprocessing(TTree *tree=0);
   virtual ~NtupleReprocessing();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string, Int_t iDet = 1, int ND_iDet = 1,  Long64_t max_entry = -1);
   virtual Bool_t   Notify();
   virtual TString  InFile();
   virtual Double_t CalcLepton( Double_t detect_dist );
   virtual void     Show(Long64_t entry = -1);

   Utils utils;

};

#endif

#ifdef NtupleReprocessing_cxx
NtupleReprocessing::NtupleReprocessing(TString file) : fChain(0) {

  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.  
  gROOT -> ProcessLine(".L LorentzVectorLight.h+");

  TTree *tree;
  if (file.Length() == 0){
    std::cout << "Opening test file ntuple_example.root" << std::endl;
    infile = new TFile("ntuple_example.root");
    TDirectory * dir = (TDirectory*)infile->Get("ntuple_example.root:/nuana");
    dir->GetObject("EventsTot",tree);
  }
  else{
    std::cout << "Opening requested file " << file << std::endl;
    infile = new TFile( file );
    infile->GetObject("nuana/EventsTot",tree);
  }

  if (!infile || !infile->IsOpen()) {
    std::cout << "Could not find specified file!!" << std::endl;
    exit(1);
  }
  Init(tree);

}

NtupleReprocessing::~NtupleReprocessing()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t NtupleReprocessing::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NtupleReprocessing::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NtupleReprocessing::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   GeniePDG = 0;
   GenieE = 0;
   GeniePx = 0;
   GeniePy = 0;
   GeniePz = 0;
   GenieProc = 0;
   neutMom = 0;
   MuonPos = 0;
   MuonMom = 0;
   ElectronPos = 0;
   ElectronMom = 0;
   p1PhotonConversionPos = 0;
   p1PhotonConversionMom = 0;
   p2PhotonConversionPos = 0;
   p2PhotonConversionMom = 0;
   miscPhotonConversionPos = 0;
   miscPhotonConversionMom = 0;
   PionPos = 0;
   PionMom = 0;
   ChargedPionPos = 0;
   ChargedPionMom = 0;
   ChargePionSign = 0;
   MultiWeight = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("iflux", &iflux, &b_iflux);
   fChain->SetBranchAddress("ibkg", &ibkg, &b_ibkg);
   fChain->SetBranchAddress("nuchan", &nuchan, &b_nuchan);
   fChain->SetBranchAddress("inno", &inno, &b_inno);
   fChain->SetBranchAddress("enugen", &enugen, &b_enugen);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("nuleng", &nuleng, &b_nuleng);
   fChain->SetBranchAddress("parid", &parid, &b_parid);
   fChain->SetBranchAddress("ndecay", &ndecay, &b_ndecay);
   fChain->SetBranchAddress("wgt", &wgt, &b_wgt);
   fChain->SetBranchAddress("NPi0", &NPi0, &b_NPi0);
   fChain->SetBranchAddress("NPi0FinalState", &NPi0FinalState, &b_NPi0FinalState);
   fChain->SetBranchAddress("NGamma", &NGamma, &b_NGamma);
   fChain->SetBranchAddress("FoundPhotons", &FoundPhotons, &b_FoundAllPhotons);
   fChain->SetBranchAddress("GeniePDG", &GeniePDG, &b_GeniePDG);
   fChain->SetBranchAddress("GenieE", &GenieE, &b_GenieE);
   fChain->SetBranchAddress("GeniePx", &GeniePx, &b_GeniePx);
   fChain->SetBranchAddress("GeniePy", &GeniePy, &b_GeniePy);
   fChain->SetBranchAddress("GeniePz", &GeniePz, &b_GeniePz);
   fChain->SetBranchAddress("GenieProc", &GenieProc, &b_GenieProc);
   fChain->SetBranchAddress("isCC", &isCC, &b_isCC);
   fChain->SetBranchAddress("mode", &mode, &b_mode);
   fChain->SetBranchAddress("Vx", &Vx, &b_Vx);
   fChain->SetBranchAddress("Vy", &Vy, &b_Vy);
   fChain->SetBranchAddress("Vz", &Vz, &b_Vz);
   fChain->SetBranchAddress("ParVx", &ParVx, &b_ParVx);
   fChain->SetBranchAddress("ParVy", &ParVy, &b_ParVy);
   fChain->SetBranchAddress("ParVz", &ParVz, &b_ParVz);
   fChain->SetBranchAddress("ParPx", &ParPx, &b_ParPx);
   fChain->SetBranchAddress("ParPy", &ParPy, &b_ParPy);
   fChain->SetBranchAddress("ParPz", &ParPz, &b_ParPz);
   fChain->SetBranchAddress("LepPx", &LepPx, &b_LepPx);
   fChain->SetBranchAddress("LepPy", &LepPy, &b_LepPy);
   fChain->SetBranchAddress("LepPz", &LepPz, &b_LepPz);
   fChain->SetBranchAddress("pdpx", &pdpx, &b_pdpx);
   fChain->SetBranchAddress("pdpy", &pdpy, &b_pdpy);
   fChain->SetBranchAddress("pdpz", &pdpz, &b_pdpz);
   fChain->SetBranchAddress("pppx", &pppx, &b_pppx);
   fChain->SetBranchAddress("pppy", &pppy, &b_pppy);
   fChain->SetBranchAddress("pppz", &pppz, &b_pppz);
   fChain->SetBranchAddress("tpx", &tpx, &b_tpx);
   fChain->SetBranchAddress("tpy", &tpy, &b_tpy);
//    fChain->SetBranchAddress("tpx", &tpx, &b_tpx);
   fChain->SetBranchAddress("ptype", &ptype, &b_ptype);
   fChain->SetBranchAddress("tptype", &tptype, &b_tptype);
   fChain->SetBranchAddress("neutMom", &neutMom, &b_neutMom);
   fChain->SetBranchAddress("ThetaLep", &ThetaLep, &b_ThetaLep);
   fChain->SetBranchAddress("PhiLep", &PhiLep, &b_PhiLep);
   fChain->SetBranchAddress("ThetaLepSmeared", &ThetaLepSmeared, &b_ThetaLepSmeared);
   fChain->SetBranchAddress("PhiLepSmeared", &PhiLepSmeared, &b_PhiLepSmeared);
   fChain->SetBranchAddress("Elep", &Elep, &b_Elep);
   fChain->SetBranchAddress("ElepSmeared", &ElepSmeared, &b_ElepSmeared);
   fChain->SetBranchAddress("MuonPos", &MuonPos, &b_MuonPos);
   fChain->SetBranchAddress("MuonMom", &MuonMom, &b_MuonMom);
   fChain->SetBranchAddress("ElectronPos", &ElectronPos, &b_ElectronPos);
   fChain->SetBranchAddress("ElectronMom", &ElectronMom, &b_ElectronMom);
   fChain->SetBranchAddress("p1PhotonConversionPos", &p1PhotonConversionPos, &b_p1PhotonConversionPos);
   fChain->SetBranchAddress("p1PhotonConversionMom", &p1PhotonConversionMom, &b_p1PhotonConversionMom);
   fChain->SetBranchAddress("p2PhotonConversionPos", &p2PhotonConversionPos, &b_p2PhotonConversionPos);
   fChain->SetBranchAddress("p2PhotonConversionMom", &p2PhotonConversionMom, &b_p2PhotonConversionMom);
   fChain->SetBranchAddress("miscPhotonConversionPos", &miscPhotonConversionPos, &b_miscPhotonConversionPos);
   fChain->SetBranchAddress("miscPhotonConversionMom", &miscPhotonConversionMom, &b_miscPhotonConversionMom);
   fChain->SetBranchAddress("PionPos", &PionPos, &b_PionPos);
   fChain->SetBranchAddress("PionMom", &PionMom, &b_PionMom);
   fChain->SetBranchAddress("ChargedPionPos", &ChargedPionPos, &b_ChargedPionPos);
   fChain->SetBranchAddress("ChargedPionMom", &ChargedPionMom, &b_ChargedPionMom);
   fChain->SetBranchAddress("ChargePionSign", &ChargePionSign, &b_ChargePionSign);
   fChain->SetBranchAddress("MultiWeight", &MultiWeight, &b_MultiWeight);
   Notify();
}

Bool_t NtupleReprocessing::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NtupleReprocessing::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NtupleReprocessing::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

 TString NtupleReprocessing::InFile(){

   return infile->GetName();

 }

 Double_t NtupleReprocessing::CalcLepton( Double_t detect_dist ){

   LepPx = 0.0;
   LepPy = 0.0;
   LepPz = 0.0;


   for( unsigned long i = 0; i < GeniePDG->size(); i++ ){
     if( abs(GeniePDG->at(i)) == 11 || abs(GeniePDG->at(i)) == 12 ||
	 abs(GeniePDG->at(i)) == 13 || abs(GeniePDG->at(i)) == 14 ){

       LepPx = GeniePx->at(i);
       LepPy = GeniePy->at(i);
       LepPz = GeniePz->at(i);
     }
   }

   ThetaLep = utils.GetTheta( LepPx, LepPy, LepPz );
   PhiLep = utils.GetPhi( LepPx, LepPy );

   TVector3 length( ParVx - Vx, ParVy - Vy, detect_dist + Vz - ParVz);

   return length.Mag();

 }



#endif // #ifdef NtupleReprocessing_cxx

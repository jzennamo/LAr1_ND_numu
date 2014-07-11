

#include "Utils.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"

#include <iostream>
#include <iomanip>
#include <string>

Utils::Utils(){

  std::cout <<"Initializing Utils"<<endl;

  PotNormNubar  = 10e20;
  PotNormNu     = 6.6e20;

}

Utils::~Utils(){
}

Double_t Utils::GetPOTNorm(Int_t iflux, Int_t iDet){

  // Simple function, really, but put into utils to abstract it out of reprocessing.
  // 
  // ideally, this should involve a database lookup.  But I'm going to hardcode it.
  // Putting in values of 1 right now but that'll update with the final POT Numbers:
  // 
  

  Double_t POT_Sim = 1;
  Double_t POTnorm = 1;

  if (iflux == kNu)
  {
    POTnorm = PotNormNu;
    if (iDet == kND){//Near Det
      //POTnorm /= 3.0;
      POT_Sim = 7.879e19;
    }
    else if(iDet == kND_150){
      POT_Sim = 1.536e20;
    }
    else if(iDet == kND_200){
      POT_Sim = 2.591e20;
    }
    else if(iDet == kIC_offaxis){
      POT_Sim = 2.386e20;
    }
    else if (iDet == kND_long){
      POT_Sim = 2.5862e20;
    }
    else if (iDet == kuB){ //uboone
      POT_Sim = 6.053e21;
    }
    else if (iDet == kFD || iDet == kIC){ //LAr1-FD or ICARUS
      POT_Sim = 5.434e20;
    }
    else{
      POT_Sim = 1;
    }
  }
  else if (iflux == kNubar)
  {
    POTnorm = PotNormNubar;
    if (iDet == kND){//Near Det
      POT_Sim = 8.87e20;
    }
    else if (iDet == kuB){ //uboone
      POT_Sim = 3.70e22;
    }
    else if (iDet == kFD || iDet == kIC){ //LAr1-FD
      POT_Sim = 1.30e21;
    }
    else{ //unknown detector
      POT_Sim = 1;
    }
  }
  else if (iflux == kNu_Fosc)
  {
    POTnorm = PotNormNu;
    if (iDet == kND){//Near Det
      POT_Sim = 1.79e20;
    }
    else if (iDet == kuB){ //uboone
      POT_Sim = 5.37e21;
    }
    else if (iDet == kFD || iDet == kIC){ //LAr1-FD
      POT_Sim = 5.32e20;
    }
    else{ //unknown detector
      POT_Sim = 1;
    }
  }
  else if (iflux == kNubar_Fosc)
  {
    POTnorm = PotNormNubar;
    if (iDet == kND){//Near Det
      POT_Sim = 5.96e20;
    }
    else if (iDet == kuB){ //uboone
      POT_Sim = 2.88e22;
    }
    else if (iDet == kFD || iDet == kIC){ //LAr1-FD
      POT_Sim = 1.08e21;
    }
    else{ //unknown detector
      POT_Sim = 1;
    }
  }
  else
  {
    POT_Sim = -1; //unknown baseline
  }

  POTnorm /= POT_Sim;

  if (POTnorm > 1e-5) return POTnorm;
  else return 1.0;

}




//=======================================================================================
// Reweight the flux using histograms in FluxRW tools
//=======================================================================================
Double_t Utils::GetFluxWeight( Double_t energy, Int_t iflux, Int_t inno, Int_t ndecay ){

  Double_t wgt = 0;
  Int_t ntype = 0;
  Int_t ptype = 0;
  string ntypestr[4] = { "nue", "nuebar", "numu", "numubar" };
  string ptypestr = "";

  // normal samples
  if( iflux == kNu || iflux == kNubar ){  
    if      (inno == 12)  ntype = 1;
    else if (inno == -12) ntype = 2;
    else if (inno == 14)  ntype = 3;
    else if (inno == -14) ntype = 4;
  }
  // full osc samples
  else if( iflux == kNu_Fosc || iflux == kNubar_Fosc ){  
    if      (inno == 12)  ntype = 3;
    else if (inno == -12) ntype = 4;
    else if (inno == 14)  ntype = 3;
    else if (inno == -14) ntype = 4;
  }

  // Determine ptype using that list described in NtupleReprocessing.C
  if (ndecay > 0 && ndecay < 5){ ptype = 3; ptypestr = "K0L"; }
  else if (ndecay > 4 && ndecay < 11){ ptype = 4; ptypestr = "K+"; }
  else if (ndecay == 11 || ndecay == 12){ ptype = 1; ptypestr = "muon"; }
  else if (ndecay == 13 || ndecay == 14){ ptype = 2; ptypestr = "pion"; }
  
  if( ntype == 0 || ptype == 0) {
    std::cout << "Wrong type neutrino or parent is present. Returning weight = 0." << std::endl;
    std::cout << "ntype = " << ntype << ",  ptype = " << ptype << std::endl;
    return 0;
  }

  // neutrino mode
  if( iflux == kNu || iflux == kNu_Fosc ) 
    wgt = fluxrw_nu.GetWeight( energy, ntype, ptype );
  // antineutrino mode
  else if( iflux == kNubar || iflux == kNubar_Fosc )
    wgt = fluxrw_nubar.GetWeight( energy, ntype, ptype );

  //  std::cout << std::setprecision(3) << "neutrino = " << ntypestr[ntype-1] << " (" << energy << " GeV), parent = " 
  //	    << ptypestr << ", flux weight = " << wgt << std::endl;
  
  return wgt;

}


//=========================================================================================
// Polar angle Theta - angle between momentum vector and z-axis
//=========================================================================================
Double_t Utils::GetTheta( TVector3 mom, TVector3 ref ){

  return GetTheta( mom.X(), mom.Y(), mom.Z(), ref );
  
}

Double_t Utils::GetTheta( Double_t px, Double_t py, Double_t pz, TVector3 ref ){

  // a.b = |a| |b| cos(theta)
  Double_t cosTheta = (ref.X()*px + ref.Y()*py + ref.Z()*pz) / (sqrt(px*px + py*py + pz*pz) * ref.Mag());
  Double_t theta = TMath::ACos(cosTheta);
  return theta;

}

//=========================================================================================
// Polar angle Phi - angle in the x-y plane, x-axis = 0
//=========================================================================================
Double_t Utils::GetPhi( TVector3 mom, TVector3 ref ){

  mom.RotateUz(ref.Unit());
  return GetPhi( mom.X(), mom.Y() );
  
}

Double_t Utils::GetPhi( Double_t px, Double_t py ){

  return atan2( py, px );

}

//=========================================================================================
// CCQE neutrino energy
//=========================================================================================
Double_t Utils::NuEnergyCCQE( Double_t l_energy, Double_t l_p, Double_t l_theta, Double_t mass, Int_t mode ){
  
  Double_t M_n = 939.565;    // MeV/c2
  Double_t M_p = 938.272;    // MeV/c2
  Double_t bindingE = 30.0;  // MeV

  //  std::cout << "NuEnergyCCQE w/ binding E = " << bindingE << " and mass = " << mass << std::endl;
  // std::cout << "flux mode = " << mode << std::endl;
  // std::cout << "l_energy = " << l_energy << "; l_p = " << l_p << "; l_theta = " << l_theta << std::endl;
  
  Double_t nu_energy = -1000;
  
  if( mode == 0 || mode == 1 ) {
    Double_t nu_energy_num = pow(M_n,2) - pow(M_p - bindingE,2)
      - pow(mass,2) + 2.0*(M_p - bindingE)*l_energy;
    Double_t nu_energy_den = 2.0*(M_p - bindingE - l_energy + l_p*cos(l_theta));
    //    std::cout << " Nu E numerator   = " << nu_energy_num << std::endl;
    //  std::cout << " Nu E denominator = " << nu_energy_den << std::endl;
    if( nu_energy_den ) nu_energy = nu_energy_num / nu_energy_den;
 
  } else if ( mode == 2 || mode == 3 ) {
    Double_t nu_energy_num = pow(M_p,2) - pow(M_n - bindingE,2)
      - pow(mass,2) + 2.0*(M_n - bindingE)*l_energy;
    Double_t nu_energy_den = 2.0*(M_n - bindingE - l_energy + l_p*cos(l_theta));
    //    std::cout << " Nu E numerator   = " << nu_energy_num << std::endl;
    //  std::cout << " Nu E denominator = " << nu_energy_den << std::endl;
    if( nu_energy_den ) nu_energy = nu_energy_num / nu_energy_den;
  }
  //  std::cout << " Nu Energy = " << nu_energy << " MeV" << std::endl;

  if( nu_energy >= 0 && nu_energy < 10000 ) 
    return nu_energy;
  else
    return -1000;
}


//=========================================================================================
// Calorimetric energy reconstruction
//=========================================================================================
Double_t Utils::NuEnergyCalo(  std::vector<Int_t> *pdg, std::vector<Double_t> *energy,
			       Bool_t include_neutrons, Bool_t include_pizeros, Double_t prot_thresh ){

  //  std::cout << "Determining calorimetric energy" << std::endl;

  Double_t M_n = .939565;    // GeV/c2
  Double_t M_p = .938272;    // GeV/c2

  Double_t total_energy = 0;

  for( unsigned int i = 0; i < pdg->size(); ++i ){
   
    //std::cout << " part: " << pdg->at(i) << ", energy: " << energy->at(i) << std::endl;

    if( abs(pdg->at(i)) > 5000 ) continue;                           // ignore nuclear fragments
    if( abs(pdg->at(i)) == 12 || abs(pdg->at(i)) == 14 ) continue;   // ignore neutrinos

    if( !include_neutrons && abs(pdg->at(i)) == 2112 ) continue;                // neutrons
    if( !include_pizeros  && abs(pdg->at(i)) == 111  ) continue;                // pizeros
    if( abs(pdg->at(i)) == 2212 && energy->at(i)-M_p < prot_thresh ) continue;  // proton threshold

    total_energy += energy->at(i);

    if( abs(pdg->at(i)) == 2212 ) total_energy -= M_p;
    if( abs(pdg->at(i)) == 2112 ) total_energy -= M_n;
  }

  //  std::cout << " Total Energy = " << total_energy << " GeV" << std::endl;
  
  return total_energy;
  
}

//=========================================================================================
// Calorimetric energy reconstruction with smearing
//=========================================================================================
Double_t Utils::NuEnergyCalo_smeared(  std::vector<Int_t> *pdg, std::vector<Double_t> *energy,
				       Double_t prot_thresh, double smeared_lep_E, double true_E, 
				       std::string LEP, double HadronicRes ){


  Double_t M_n = .939565;    // GeV/c2
  Double_t M_p = .938272;    // GeV/c2
  Double_t total_energy = 0;
  TRandom *smear = new TRandom3(0);

  if(smeared_lep_E == 0){return total_energy;}

  for( unsigned int i = 0; i < pdg->size(); ++i ){
   

    if( abs(pdg->at(i)) > 5000 ) continue;                           // ignore nuclear fragments
    if( abs(pdg->at(i)) == 12 || abs(pdg->at(i)) == 14 ) continue;   // ignore neutrinos

    if( abs(pdg->at(i)) == 2112 ) continue;                // neutrons
    if( abs(pdg->at(i)) == 111  ) continue;                // pizeros
    if( abs(pdg->at(i)) == 2212 && energy->at(i)-M_p < prot_thresh ) continue;  // proton threshold
    

    if((LEP == "muon" && abs(pdg->at(i)) == 13) ||
       (LEP == "elec" && abs(pdg->at(i)) == 11)){ 
      
      total_energy += (smeared_lep_E); 
      continue;
   
    }
    if(LEP == "NC"){
    
      total_energy += (smeared_lep_E);
      total_energy -= (true_E);
    
    }

       
    total_energy += (smear->Gaus(energy->at(i), HadronicRes*energy->at(i)));

    if( abs(pdg->at(i)) == 2212 ) total_energy -= M_p;
    if( abs(pdg->at(i)) == 2112 ) total_energy -= M_n;

  }

  delete smear;
  return total_energy;
  
  
}


//=========================================================================================
// Get the visible energy near the interaction vertex
//=========================================================================================
Double_t Utils::VertexEnergy( std::vector<Int_t> *pdg, std::vector<Double_t> *energy, 
			      Double_t prot_thresh, Double_t pion_thresh ){

  std::cout << "Determining kinetic energy near vertex" << std::endl;

  Double_t M_p  = .938272;    // GeV/c2
  Double_t M_pi = .139570;    // GeV/c2

  Double_t total_energy = 0;
  
  for( unsigned int i = 0; i < pdg->size(); ++i ){
   
    std::cout << " part: " << pdg->at(i) << ", energy: " << energy->at(i) << std::endl;
    
    //if( abs(pdg->at(i)) > 5000 ) continue;                           // ignore nuclear fragments
    //if( abs(pdg->at(i)) == 12 || abs(pdg->at(i)) == 14 ) continue;   // ignore neutrinos
    //if( abs(pdg->at(i)) == 2112 ) continue;                          // neutrons
    //if( abs(pdg->at(i)) == 111  ) continue;                          // pizeros
   
    if( abs(pdg->at(i)) == 2212 && energy->at(i)-M_p > prot_thresh ){  // proton threshold
      total_energy += (energy->at(i) - M_p);
    }
    if( abs(pdg->at(i)) == 211  && energy->at(i)-M_pi > pion_thresh ){  // pion threshold
      total_energy += (energy->at(i) - M_pi);
    }
   
  }
  
  std::cout << " Total Kinetic Energy = " << total_energy << " GeV" << std::endl;
  
  return total_energy;
  
}


//=========================================================================================
// Get total photon energies in fiducial volume
//=========================================================================================
Double_t Utils::TotalPhotonEnergy( Int_t idet, Int_t ND_iDet, std::vector<gan::LorentzVectorLight> *p1pos,
				   std::vector<gan::LorentzVectorLight> *p1mom,
				   std::vector<gan::LorentzVectorLight> *p2pos,
				   std::vector<gan::LorentzVectorLight> *p2mom ){
  
  if( p1pos->size() != p1mom->size() || p2pos->size() != p2mom->size() ){
    std::cout << "photon vectors don't match!!" << std::endl;
    exit(1);
  }
  
  Double_t energy = 0.0;

  for( unsigned int i = 0; i < p1pos->size(); i++ ){

    TVector3 photon1Pos((*p1pos).at(i).X(), (*p1pos).at(i).Y(), (*p1pos).at(i).Z() );
    if( IsFiducial( idet, ND_iDet, photon1Pos ) )
      energy += p1mom->at(i).E();

  }
  for( unsigned int i = 0; i < p2pos->size(); i++ ){

    TVector3 photon2Pos((*p2pos).at(i).X(), (*p2pos).at(i).Y(), (*p2pos).at(i).Z() );
    if( IsFiducial( idet, ND_iDet, photon2Pos ) )
      energy += p2mom->at(i).E();

  }

  return energy;

} 

//=========================================================================================
// Check if point is in some fiducial volume definition
//==========================================================================================
void Utils::GetDetBoundary( Int_t idet, Int_t ND_iDet, Double_t &xmin, Double_t &xmax, 
			    Double_t &ymin, Double_t &ymax, Double_t &zmin, Double_t &zmax ){
 
  if(ND_iDet == 1 && idet == 0){
    
    nd_xmin = -200.0;
    nd_xmax =  200.0;
    nd_ymin = -242.0;
    nd_ymax =  158.0;
    nd_zmin =  0.0;
    nd_zmax =  365.0;

    fidCut_xmax = 10.0; fidCut_xmin = 10.0; fidCut_xAPA = 0.0;  
    fidCut_ymax = 10.0; fidCut_ymin = 10.0; 
    fidCut_zmax = 80.0; fidCut_zmin = 10.0; }

  else if(ND_iDet == 4 && idet == 0){
    

    /*    nd_xmin = -80.0; nd_xmax =  320.0;
    nd_ymin = -242.0; nd_ymax =  158.0;
    nd_zmin =  0.0;   nd_zmax =  365.0;
    */
    nd_xmin = -65.0; nd_xmax =  305.0;
    nd_ymin = -143; nd_ymax =  227.0;
    nd_zmin =  0.0;   nd_zmax =  300.0;
   
    fidCut_xmax = 10.0;  fidCut_xmin = 10.0; fidCut_xAPA = 0.0;
    fidCut_ymax = 10.0;  fidCut_ymin = 10.0;  
    fidCut_zmax = 80.00;  fidCut_zmin = 10.0; }

  else if(idet != 0){     
    nd_xmin = -200.0;
    nd_xmax =  200.0;
    nd_ymin = -242.0;
    nd_ymax =  158.0;
    nd_zmin =  0.0;
    nd_zmax =  365.0;

    fidCut_xmax = 10.0; fidCut_xmin = 10.0; fidCut_xAPA = 0.0;
    fidCut_ymax = 10.0; fidCut_ymin = 10.0;
    fidCut_zmax = 80.0; fidCut_zmin = 10.0;}

  ub_xmin =  0.0;  ub_xmax =  256.0;
  ub_ymin = -116.5;  ub_ymax =  116.5;
  ub_zmin =  0.0;  ub_zmax =  1036.0;
  
  fd_xmin =  -450.0;  fd_xmax =  450.0;
  fd_ymin =  -450.0;  fd_ymax =  450.0;
  fd_zmin =  0.0; fd_zmax =  900.0;

  ic_xmin =  -350.0; ic_xmax =  350.0;
  ic_ymin =  -158.0; ic_ymax =  158.0;
  ic_zmin =  0.0;    ic_zmax =  1795.0;

  nd_long_xmin = -80.0;  nd_long_xmax =  320.0;
  nd_long_ymin = -242.0; nd_long_ymax =  158.0;
  nd_long_zmin =  0.0;   nd_long_zmax =  2*365.0;

  // LAr1-ND
  if( idet == kND || idet == kND_150 || idet == kND_200  ){
    xmin = nd_xmin;
    xmax = nd_xmax;
    ymin = nd_ymin;
    ymax = nd_ymax;
    zmin = nd_zmin;
    zmax = nd_zmax;
  }
  // Uboone
  else if( idet == kuB ){
    xmin = ub_xmin;
    xmax = ub_xmax;
    ymin = ub_ymin;
    ymax = ub_ymax;
    zmin = ub_zmin;
    zmax = ub_zmax;
  }
  // LAr1-FD
  else if( idet == kFD ){
    xmin = fd_xmin;
    xmax = fd_xmax;
    ymin = fd_ymin;
    ymax = fd_ymax;
    zmin = fd_zmin;
    zmax = fd_zmax;
  }
  else if( idet == kIC || idet == kIC_offaxis){
    xmin = ic_xmin;
    xmax = ic_xmax;
    ymin = ic_ymin;
    ymax = ic_ymax;
    zmin = ic_zmin;
    zmax = ic_zmax;
  }
  else if( idet == kND_long ){
    xmin = nd_long_xmin;
    xmax = nd_long_xmax;
    ymin = nd_long_ymin;
    ymax = nd_long_ymax;
    zmin = nd_long_zmin;
    zmax = nd_long_zmax;
  }

  else{
    std::cout << "Don't recognize idet" << std::endl;
    exit(1);
  }
  return;

}

//=========================================================================================
// Check if point is in some fiducial volume definition
//==========================================================================================
Bool_t Utils::IsFiducial( Int_t idet, Int_t ND_iDet, TVector3 vtx){
  
  if(ND_iDet == 1 && idet == 0){
    
    nd_xmin = -200.0;
    nd_xmax =  200.0;
    nd_ymin = -242.0;
    nd_ymax =  158.0;
    nd_zmin =  0.0;
    nd_zmax =  365.0;   

    fidCut_xmax = 10.0; fidCut_xmin = 10.0; fidCut_xAPA = 0.0;  
    fidCut_ymax = 10.0; fidCut_ymin = 10.0; 
    fidCut_zmax = 80.0; fidCut_zmin = 10.0; }

  else if(ND_iDet == 4 && idet == 0){


    /*    nd_xmin = -80.0; nd_xmax =  320.0;
    nd_ymin = -242.0; nd_ymax =  158.0;
    nd_zmin =  0.0;   nd_zmax =  365.0;
    */
    nd_xmin = -65.0; nd_xmax =  305.0;
    nd_ymin = -143; nd_ymax =  227.0;
    nd_zmin =  0.0;   nd_zmax =  300.0;

    fidCut_xmax = 10.0;  fidCut_xmin = 10.0; fidCut_xAPA = 0.0;
    fidCut_ymax = 10.0;  fidCut_ymin = 10.0;  
    fidCut_zmax = 80.00;  fidCut_zmin = 10.0; }

  else if(idet != 0){
    nd_xmin = -200.0;
    nd_xmax =  200.0;
    nd_ymin = -242.0;
    nd_ymax =  158.0;
    nd_zmin =  0.0;
    nd_zmax =  365.0;

    fidCut_xmax = 10.0; fidCut_xmin = 10.0; fidCut_xAPA = 0.0;
    fidCut_ymax = 10.0; fidCut_ymin = 10.0;
    fidCut_zmax = 80.0; fidCut_zmin = 10.0;}

  ub_xmin =  0.0;  ub_xmax =  256.0;
  ub_ymin = -116.5;  ub_ymax =  116.5;
  ub_zmin =  0.0;  ub_zmax =  1036.0;
  
  fd_xmin =  -450.0;  fd_xmax =  450.0;
  fd_ymin =  -450.0;  fd_ymax =  450.0;
  fd_zmin =  0.0; fd_zmax =  900.0;

  ic_xmin =  -350.0; ic_xmax =  350.0;
  ic_ymin =  -158.0; ic_ymax =  158.0;
  ic_zmin =  0.0;    ic_zmax =  1795.0;

  nd_long_xmin = -80.0;  nd_long_xmax =  320.0;
  nd_long_ymin = -242.0; nd_long_ymax =  158.0;
  nd_long_zmin =  0.0;   nd_long_zmax =  2*365.0;

  Double_t xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), zmin(0.0), zmax(0.0);

  GetDetBoundary( idet, ND_iDet, xmin, xmax, ymin, ymax, zmin, zmax );

  // fidCut_x = fidCut;
  // fidCut_y = fidCut;
  // fidCut_z = fidCut;

  xmin += fidCut_xmin;
  xmax -= fidCut_xmax;
  ymin += fidCut_ymin;
  ymax -= fidCut_ymax;
  zmin += fidCut_zmin; 
  zmax -= fidCut_zmax;

  if( idet == kND || idet == kND_150 || idet == kND_200 ){  // Cut out volume around each side of the center of x:
    
    Double_t xcenter = (nd_xmax + nd_xmin)/2;  // gets the center.  If the cuts above on x aren't symmetric, this is wrong.
   
    if ( vtx.X() < (xcenter + fidCut_xAPA) && vtx.X() > (xcenter - fidCut_xAPA) ) return false; //cut around the center APA

    //    if( vtx.X() < (6 + fidCut_xmin) || vtx.X() > (245 - fidCut_xmax) ||
    //	vtx.Y() < (-161 + fidCut_ymin) || vtx.Y() > (77 - fidCut_ymin)) return false;

  }

  if( idet == kIC || idet == kIC_offaxis ){  // Cut out volume around each side of the center of x:
       
    // gets the center.  If the cuts above on x aren't symmetric, this is wrong.
    Double_t xcenter = (ic_xmax + ic_xmin)/2; 
    if ( vtx.X() < (xcenter + 50 + fidCut_xmax) && vtx.X() > (xcenter - 50 -fidCut_xmax) ) return false; //cut around the center APA
    Double_t zNCunit = 365;

    //    if( (vtx.Z() > (364 - fidCut_zmax)  && vtx.Z() < (366 + fidCut_zmin)) || 
    //	(vtx.Z() > (730 - fidCut_zmax)  && vtx.Z() < (732 + fidCut_zmin)) || 
    //	(vtx.Z() > (1096- fidCut_zmax)  && vtx.Z() < (1460+ fidCut_zmin))) return false;

  }

  if( vtx.X() > xmin && vtx.X() < xmax && vtx.Y() > ymin && vtx.Y() < ymax && vtx.Z() > zmin && vtx.Z() < zmax )
    return true;
  else
    return false;

}

//=========================================================================================
// Check if point is in some active volume definition
//==========================================================================================
Bool_t Utils::IsActive( Int_t idet, Int_t ND_iDet, TVector3 vtx, double cut ){
  
  if(ND_iDet == 1 && idet == 0){
    
    nd_xmin = -200.0;
    nd_xmax =  200.0;
    nd_ymin = -242.0;
    nd_ymax =  158.0;
    nd_zmin =  0.0;
    nd_zmax =  365.0;

    fidCut_xmax = 10.0; fidCut_xmin = 10.0; fidCut_xAPA = 0.0;  
    fidCut_ymax = 10.0; fidCut_ymin = 10.0; 
    fidCut_zmax = 80.0; fidCut_zmin = 10.0; }

  else if(ND_iDet == 4 && idet == 0){
    
    /*    nd_xmin = -80.0; nd_xmax =  320.0;
    nd_ymin = -242.0; nd_ymax =  158.0;
    nd_zmin =  0.0;   nd_zmax =  365.0;
    */
    nd_xmin = -65.0; nd_xmax =  305.0;
    nd_ymin = -143; nd_ymax =  227.0;
    nd_zmin =  0.0;   nd_zmax =  300.0;

    fidCut_xmax = 10.0;  fidCut_xmin = 10.0; fidCut_xAPA = 0.0;
    fidCut_ymax = 10.0;  fidCut_ymin = 10.0;  
    fidCut_zmax = 80.00;  fidCut_zmin = 10.0; }
    
  else if(idet != 0){
    nd_xmin = -200.0;
    nd_xmax =  200.0;
    nd_ymin = -242.0;
    nd_ymax =  158.0;
    nd_zmin =  0.0;
    nd_zmax =  365.0;
    
    fidCut_xmax = 10.0; fidCut_xmin = 10.0; fidCut_xAPA = 0.0;
    fidCut_ymax = 10.0; fidCut_ymin = 10.0;
    fidCut_zmax = 80.0; fidCut_zmin = 10.0;}
  
  ub_xmin =  0.0;  ub_xmax =  256.0;
  ub_ymin = -116.5;  ub_ymax =  116.5;
  ub_zmin =  0.0;  ub_zmax =  1036.0;
  
  fd_xmin =  -450.0;  fd_xmax =  450.0;
  fd_ymin =  -450.0;  fd_ymax =  450.0;
  fd_zmin =  0.0; fd_zmax =  900.0;

  ic_xmin =  -350.0; ic_xmax =  350.0;
  ic_ymin =  -158.0; ic_ymax =  158.0;
  ic_zmin =  0.0;    ic_zmax =  1795.0;

  nd_long_xmin = -80.0;  nd_long_xmax =  320.0;
  nd_long_ymin = -242.0; nd_long_ymax =  158.0;
  nd_long_zmin =  0.0;   nd_long_zmax =  2*365.0;

  Double_t xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), zmin(0.0), zmax(0.0);

  GetDetBoundary( idet, ND_iDet, xmin, xmax, ymin, ymax, zmin, zmax );

  if( idet == kND || idet == kND_150 || idet == kND_200){

    //    if(vtx.X() < (6 - cut) || vtx.X() > (245 + cut) ||
    //   vtx.Y() < (-161 - cut) || vtx.Y() > (77 + cut)) return false;
  }

  if( idet == kIC || idet == kIC_offaxis ){  // Cut out volume around each side of the center of x:   
    // gets the center.  If the cuts above on x aren't symmetric, this is wrong.
    Double_t xcenter = (ic_xmax + ic_xmin)/2; 
    //cut around the center APA
    if ( vtx.X() < (xcenter + 50 + cut) && vtx.X() > (xcenter - 50 - cut) ) return false; 

    //    if( (vtx.Z() > (364-cut)  && vtx.Z() < (366+cut)) ||
    //    (vtx.Z() > (730-cut)  && vtx.Z() < (732+cut)) ||
    //        (vtx.Z() > (1096-cut)  && vtx.Z() < (1460+cut))) return false;

  }

  if( vtx.X() > xmin+cut && vtx.X() < xmax-cut && vtx.Y() > ymin+cut && vtx.Y() < ymax-cut && vtx.Z() > zmin+cut && vtx.Z() < zmax-cut )
    return true;
  else
    return false;

}

//=========================================================================================
// Calculate fiducial mass
//=========================================================================================
Double_t Utils::GetFidMass( Int_t idet, Int_t ND_iDet ){
  
  if(ND_iDet == 1 && idet == 0){
    
    nd_xmin = -200.0;
    nd_xmax =  200.0;
    nd_ymin = -242.0;
    nd_ymax =  158.0;
    nd_zmin =  0.0;
    nd_zmax =  365.0;

    fidCut_xmax = 10.0; fidCut_xmin = 10.0; fidCut_xAPA = 0.0;  
    fidCut_ymax = 10.0; fidCut_ymin = 10.0; 
    fidCut_zmax = 80.0; fidCut_zmin = 10.0; }

  else if(ND_iDet == 4 && idet == 0){
    
    /*       nd_xmin = -80.0; nd_xmax =  320.0;
    nd_ymin = -242.0; nd_ymax =  158.0;
    nd_zmin =  0.0;   nd_zmax =  365.0;
    */

    nd_xmin = -65.0; nd_xmax =  305.0;
    nd_ymin = -143; nd_ymax =  227.0;
    nd_zmin =  0.0;   nd_zmax =  300.0;


    fidCut_xmax = 10.0;  fidCut_xmin = 10.0; fidCut_xAPA = 0.0;
    fidCut_ymax = 10.0;  fidCut_ymin = 10.0;  
    fidCut_zmax = 80.00;  fidCut_zmin = 10.0; }
    
  else if(idet != 0){
    nd_xmin = -140.0; nd_xmax =  400.0;
    nd_ymin = -172.0; nd_ymax =  172.0;
    nd_zmin =  0.0;   nd_zmax =  330.0;

    fidCut_xmax = 10.0; fidCut_xmin = 10.0; fidCut_xAPA = 0.0;
    fidCut_ymax = 10.0; fidCut_ymin = 10.0;
    fidCut_zmax = 80.0; fidCut_zmin = 10.0;}

  ub_xmin =  0.0;  ub_xmax =  256.0;
  ub_ymin = -116.5;  ub_ymax =  116.5;
  ub_zmin =  0.0;  ub_zmax =  1036.0;
  
  fd_xmin =  -450.0;  fd_xmax =  450.0;
  fd_ymin =  -450.0;  fd_ymax =  450.0;
  fd_zmin =  0.0; fd_zmax =  900.0;

  ic_xmin =  -350.0; ic_xmax =  350.0;
  ic_ymin =  -158.0; ic_ymax =  158.0;
  ic_zmin =  0.0;    ic_zmax =  1795.0;

  nd_long_xmin = -80.0;  nd_long_xmax =  320.0;
  nd_long_ymin = -242.0; nd_long_ymax =  158.0;
  nd_long_zmin =  0.0;   nd_long_zmax =  2*365.0;

  Double_t xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), zmin(0.0), zmax(0.0);

  GetDetBoundary( idet, ND_iDet, xmin, xmax, ymin, ymax, zmin, zmax );

  xmin += fidCut_xmin;
  xmax -= fidCut_xmax;
  ymin += fidCut_ymin;
  ymax -= fidCut_ymax;
  zmin += fidCut_zmin; 
  zmax -= fidCut_zmax;

  Double_t rho = 1.4;  // g/cm^3 = tons/m^3

  Double_t total_vol = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)/pow(100.0,3) * rho;

  if( idet == kND || idet == kND_150 || idet == kND_200 ){  // Cut out volume around each side of the center of x:
    
    Double_t xcenter = (nd_xmax + nd_xmin)/2;  // gets the center.  If the cuts above on x aren't symmetric, this is wrong.
    
    total_vol -= ((xcenter + fidCut_xAPA) - (xcenter - fidCut_xAPA))*(ymax-ymin)*(zmax-zmin)/pow(100.0,3) * rho;
    
  }

  if (idet == kIC || idet == kIC_offaxis) {total_vol = (xmax-xmin-100-34)*(ymax-ymin)*(zmax-zmin)/pow(100.0,3) * rho;}

  return total_vol;

}


//=========================================================================================
// Calculate fiducial mass
//=========================================================================================
Double_t Utils::GetActiveMass( Int_t idet, Int_t ND_iDet ){
  


  if(ND_iDet == 1 && idet == 0){
    
    nd_xmin = -200.0;
    nd_xmax =  200.0;
    nd_ymin = -242.0;
    nd_ymax =  158.0;
    nd_zmin =  0.0;
    nd_zmax =  365.0;


    fidCut_xmax = 10.0; fidCut_xmin = 10.0; fidCut_xAPA = 0.0;  
    fidCut_ymax = 10.0; fidCut_ymin = 10.0; 
    fidCut_zmax = 80.0; fidCut_zmin = 10.0; }

  else if(ND_iDet == 4 && idet == 0){
    
    /*     nd_xmin = -80.0; nd_xmax =  320.0;
    nd_ymin = -242.0; nd_ymax =  158.0;
    nd_zmin =  0.0;   nd_zmax =  365.0;
    */
    nd_xmin = -65.0; nd_xmax =  305.0;
    nd_ymin = -143; nd_ymax =  227.0;
    nd_zmin =  0.0;   nd_zmax =  300.0;



    fidCut_xmax = 10.0;  fidCut_xmin = 10.0; fidCut_xAPA = 0.0;
    fidCut_ymax = 10.0;  fidCut_ymin = 10.0;  
    fidCut_zmax = 80.00;  fidCut_zmin = 10.0; }
    
  else if(idet != 0){
    nd_xmin = -140.0; nd_xmax =  400.0;
    nd_ymin = -172.0; nd_ymax =  172.0;
    nd_zmin =  0.0;   nd_zmax =  330.0;

    fidCut_xmax = 10.0; fidCut_xmin = 10.0; fidCut_xAPA = 0.0;
    fidCut_ymax = 10.0; fidCut_ymin = 10.0;
    fidCut_zmax = 80.0; fidCut_zmin = 10.0;}

  ub_xmin =  0.0;  ub_xmax =  256.0;
  ub_ymin = -116.5;  ub_ymax =  116.5;
  ub_zmin =  0.0;  ub_zmax =  1036.0;
  
  fd_xmin =  -450.0;  fd_xmax =  450.0;
  fd_ymin =  -450.0;  fd_ymax =  450.0;
  fd_zmin =  0.0; fd_zmax =  900.0;

  ic_xmin =  -350.0; ic_xmax =  350.0;
  ic_ymin =  -158.0; ic_ymax =  158.0;
  ic_zmin =  0.0;    ic_zmax =  1795.0;

  nd_long_xmin = -80.0;  nd_long_xmax =  320.0;
  nd_long_ymin = -242.0; nd_long_ymax =  158.0;
  nd_long_zmin =  0.0;   nd_long_zmax =  2*365.0;

  Double_t xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), zmin(0.0), zmax(0.0);

  GetDetBoundary( idet, ND_iDet, xmin, xmax, ymin, ymax, zmin, zmax );

  Double_t rho = 1.4;  // g/cm^3 = tons/m^3

  double total_mass = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)/pow(100.0,3) * rho;
  if (idet == kIC || idet == kIC_offaxis) 
    total_mass = (xmax-xmin-100)*(ymax-ymin)*(zmax-zmin)/pow(100.0,3) * rho;
  return total_mass;

}

// 
//Find out if line intersects the plane
// 
Bool_t Utils::IntersectsPlane(const TVector3 & startPos, const TVector3 & startDir,
                              const TVector3 & planeCorner1,
                              const TVector3 & planeCorner2,
                              const TVector3 & planeCorner3,
                              const bool verbose) const{
                              // ) const{

  if (startDir.Mag() == 0) return false;
  std::cout << "Start Pos: (" << startPos.X() <<", "<< startPos.Y() <<", " << startPos.Z() << ")" << std::endl;
  std::cout << "Start Dir: (" << startDir.X()<<", "<< startDir.Y()<<", " << startDir.Z()<< ")" << std::endl;

  // a and b are defined to be vectors parallel to the plane:
  TVector3 a = planeCorner2 - planeCorner1;
  TVector3 b = planeCorner3 - planeCorner1;

  // Normal vector to the plane
  TVector3 n = a.Cross(b);

  //Check that the Normal vector isn't 0 - that the vectors aren't parallel
  if (n.Mag() == 0){
    if (verbose) std::cout << "Error, input vectors to define the plane are parallel." << std::endl;
    return false;
  }


  // Any point on the plane can be parametrized with the following equation:
  // p = planeCorner1 + u*a + v*b
  // Any point on the line can be parametrized with:
  // l = startPos + t*startDir
  // If this line intersects this plane, then for some values of (u,v,t) p = l
  // startPos + t*startDir = planeCorner1 + u*a + v*b
  // or
  // startPos - planeCorner1 = (-startDir, a, b) . (t, u, v) (Dot product)
  // (startDir, a, b) is a matrix with each column equal to the said vector.  It's 3x3 matrix

  // Invert the matrix, so we can solve for t, u, v directly.  Only singular if the line is parallel
  // to the plane.  In which case we can return false.
  // 
  // Check that the line is not parallel to the plane with dot product to normal vector:
  if (abs(n.Dot(startDir)/(n.Mag()*startDir.Mag())) < 1e-10){
    if (verbose) std::cout << "Direction is parallel to the plane, exiting." << std::endl;
    return false;
  }
  
  // Create the matrix:
  TMatrixD mat(3,3);
  // format is mat(row_n, col_n)
  mat(0,0) = -startDir.X();
  mat(0,1) = a.X();
  mat(0,2) = b.X();
  mat(1,0) = -startDir.Y();
  mat(1,1) = a.Y();
  mat(1,2) = b.Y();
  mat(2,0) = -startDir.Z();
  mat(2,1) = a.Z();
  mat(2,2) = b.Z();

  TMatrix inv = mat.Invert();

  // Create the vector to multiply.
  TMatrix t(3,1);
  t(0,0) = startPos.X() - planeCorner1.X();
  t(1,0) = startPos.Y() - planeCorner1.Y();
  t(2,0) = startPos.Z() - planeCorner1.Z();

  // Do the matrix multiplication:
  TMatrix result = inv*t;

  // t is  largely irrelevant in this function, but u and v are important.
  // if u and v are both within [0,1] then the line intersects the plane in the desired range

  if (verbose) std::cout << "\tt, u, v are "<< result(0,0) << ", " << result(1,0) << ", " << result(2,0) << std::endl;

  if (result(0,0) > 0){ // Intersects in forward direction
    if (result(1,0) >= 0 && result (1,0) <= 1)
    {
      if (result(2,0) >= 0 && result (2,0) <= 1)
      {
        return true;
      }
    }
    // then success! they intersect
  }
  return false;

}

//###########################################################################################
//Starting muon specific functions, but they can be generalized
//###########################################################################################

//###########################################################################################
// Used to determing how the muon behaves with respect to the Ranger
//###########################################################################################

//ND_iDet corresponds to different configurations of the NearDetector
//ND_iDet == 1 is just the argon, no ranger
//ND_iDet == 2 is a non-magetized ranger tacked onto the end of the argon
//ND_iDet == 3 is a magneized ranger tacked onto the end of the argon
//ND_iDet == 4 is a totally magnetized argon volume 

bool Utils::RangeOut(TVector3 MuXYZ, TVector3 MuP, Double_t track_L, int ND_iDet){
  
  bool useful = false;
  Double_t R_xmin,R_ymin,R_zmin,R_xmax,R_ymax,R_zmax;
  //Ranger dimensions
  //Right now it starts RIGHT at the back edge of the TPC 
  //To stop a 1GeV muon with Mitch's numbers (11/20) need ~6m ranger 
  //For now, same height and width as ND
  if(ND_iDet == 2 && MuP.Mag() < 1){

    R_xmin = -140.0;   R_xmax =  400.0; 
    R_ymin = -172.0;   R_ymax =  172.0;
    R_zmin = 330;   

    double mu_z = 0;

    if(((MuP.Mag()/.100)*59)+R_zmin < 930) R_zmax = ((MuP.Mag()/.100)*59)+R_zmin;
    else{R_zmax = 0; std::cout << " Fail ::::: " << MuP.Mag() << std::endl;}

    //Resolution due to this is 100MeV, from Mitch's talk (11/20), 100MeV per segment.
  }
  
  //Magnetized, we still want to be able to resolve 2GeV muons
  //need a 4m long detector to watch a 2GeV particle deflect 1m
  //because it needs to deflect a meter lets make the Ranger
  //1m taller and wider…? BUT this is cut for "acceptance" so it stays the same width

  else if(ND_iDet == 3){ 
    TVector3 B_field(0,1,0); //only in the Y!!!-direction

    R_xmin = -240.0;   R_xmax =  500.0; 
    R_ymin = -172.0;   R_ymax =  172.0;
    R_zmin = 330;    //R_zmax =  630.0;

    //This function gives the length necessary for a particle of 
    //transverse momentum p to deviate from straight by "Dev"
    Double_t Dev = .25;
    Double_t P = MuP.Perp(B_field)/(.3*B_field.Mag());
    Double_t muon_z = nd_zmax+100*TMath::Sin(TMath::ACos((P-Dev)/P))*P;

      //sqrt(((20*Dev*MuP.Perp(B_field))/(3*B_field.Mag()))-pow(Dev,2));//cm
    
    if(muon_z < (R_zmin+400)){ R_zmax = muon_z;} //std::cout << " Dev at range : " << muon_z << std::endl; }
    else{R_zmax = (R_zmin+400);}// std::cout << " Dev at back " << std::endl;}

    //Resolution due to this is TINY, 4cm for a 4m detector 
    //filled with air if the particle hit the back plane
    //NIM 24 (1963) 381-389
    //Let's say it is 4%
  }

  else if(ND_iDet == 4){ 
    TVector3 B_field(0,-1.5,0); //only in the Y!!!-direction


    nd_xmin = -65.0; nd_xmax =  305.0;
    nd_ymin = -143; nd_ymax =  227.0;
    nd_zmin =  0.0;   nd_zmax =  300.0;

    R_xmin = -65.0;   R_xmax =  305.0;  
    R_ymin = -172.0;   R_ymax =  172.0; 
    R_zmin = 300; 

    //This function gives the length necessary for a particle of 
    //transverse momentum p to deviate from straight by "Dev"
    Double_t Dev = .3;
    Double_t P = MuP.Perp(B_field)/(.3*B_field.Mag());
    Double_t muon_z = 100*TMath::Sin(TMath::ACos((P-Dev)/P))*P;

      //sqrt(((20*Dev*MuP.Perp(B_field))/(3*B_field.Mag()))-pow(Dev,2));//cm
    
    if(muon_z < 100){ R_zmax = R_zmin+muon_z;} //std::cout << " Dev at range : " << muon_z << std::endl; }
    else{R_zmax = (R_zmin+100);}

    //Resolution due to this is TINY, 4cm for a 4m detector 
    //filled with air if the particle hit the back plane
    //NIM 24 (1963) 381-389
    //Let's say it is 4%
  }
  else if(ND_iDet == 1){
    useful = true;
    return useful;
  }

  //at the back of the ranger
  TVector3    TopRight(R_xmax, R_ymin, R_zmax); 
  TVector3     TopLeft(R_xmax, R_ymax, R_zmax); 
  TVector3 BottomRight(R_xmin, R_ymin, R_zmax);
  
  useful = IntersectsPlane(MuXYZ, MuP, TopRight, TopLeft, BottomRight ,false);
  //  if(!useful){ cout << " ================ Missed getting measured with Muon System " << std::endl; }
  
  return useful;

}

//###########################################################################################
// Returns a smeared muon energy for use in the assignment of the reco nuetrino energy
//###########################################################################################

//
//ND_iDet corresponds to different configurations of the NearDetector
//ND_iDet == 1 is just the argon, no ranger
//ND_iDet == 2 is a non-magetized ranger tacked onto the end of the argon
//ND_iDet == 3 is a magneized ranger tacked onto the end of the argon
//ND_iDet == 4 is a totally magnetized argon volume 
Double_t Utils::EnergyRes(bool contained, Double_t  energy, Double_t track_L, int  ND_iDet, bool hit_back_of_ranger, TVector3  MuPExit, TVector3  MuPInitial, TVector3  MuXYZExit, TVector3  MuXYZInitial, double Syst, int & HowMeasuresed){
  
  Double_t smeared_energy;
  Double_t Res;
  TRandom *smear = new TRandom3(0);

  
  if(contained){    
    //Assumes resolution of 2% for stopped muons
    Res = 0.02;
    smeared_energy = smear->Gaus(energy, sqrt( pow(Res,2) + pow((Syst*Res),2))*energy);
    HowMeasuresed = 1;
    //    if(fabs(smeared_energy-energy)/energy > 0.1) std::cout << ":::::   Res: " << Res << ", Smear :" << fabs(smeared_energy-energy)/energy << ", Track L : " << track_L << std::endl;

  }

  else if(!contained && ND_iDet == 2 && hit_back_of_ranger){
    
    //Resolution is 100MeV per segment (per Mitch's talk - 11/20)
    Res = (.100)/(MuPExit.Mag()); 
    smeared_energy = smear->Gaus(energy, sqrt( pow(Res,2) + pow((Syst*Res),2))*energy); 

    if(Res > /*-0.100989*log(0.0004644*track_L)*/-0.102*log(0.000612*track_L) &&  track_L > 100){ 
      Res = /*-0.100989*log(0.0004644*track_L)*/-0.102*log(0.000612*track_L);
      smeared_energy = smear->Gaus(energy, sqrt( pow(Res,2) + pow((Syst*Res),2)));
    }
   else if(Res > /*-0.100989*log(0.0004644*track_L)*/-0.102*log(0.000612*track_L) && track_L < 100){
     smeared_energy = 0;}

  }
  
  else if(!contained && ND_iDet == 3 && hit_back_of_ranger){    
    //
    //.09m is the "uncertainty" on the posistion measurement
    //".29" is the number of "samplings" AKA every .29m we have a set of MRS's

    TVector3 B_field(0,1,0);

    if(MuPExit.Perp(B_field) < 0.05 && track_L > 100){
      Res = /*-0.100989*log(0.0004644*track_L)*/-0.102*log(0.000612*track_L);
      smeared_energy = smear->Gaus(energy, sqrt( pow(Res,2) + pow((Syst*Res),2)));
    }
    else if(MuPExit.Perp(B_field) < 0.05 && track_L < 100){
      smeared_energy = 0;
    }

    else{
      Double_t Dev = .25;
      Double_t P = MuPExit.Perp(B_field)/(.3*B_field.Mag());
      Double_t muon_z = 100*TMath::Sin(TMath::ACos((P-Dev)/P))*P;
      Res = sqrt(720/(((muon_z)/29)+4))*3*(pow(MuPExit.Perp(B_field),2)/(0.003*B_field.Mag()*pow(muon_z,2)));
      smeared_energy = smear->Gaus(energy, sqrt( pow(Res,2) + pow((Syst*Res),2)));
    }
    
  }

  else if(!contained && ND_iDet == 4){
    
    //This needs a proper handling based on the amount of time that the track spends in the 
    //argon, the track should have a very small resolution on        
    //.003m is the "uncertainty" on the posistion measurement, spacing between the wires
    //"L/.003" is the number of "samplings" AKA every .003m we have a wire
    //L is the length (in meters) of the track in the RANGER(!) for now I define this as 
    //distance nessecary for the track to deviate 1 meter from its inital trajectory
    
    TVector3 B_field(0,0.5,0);
    Double_t Dev = .05;//m
    Double_t P = MuPInitial.Perp(B_field)/(.3*B_field.Mag());
    Double_t muon_z = 100*TMath::Sin(TMath::ACos((P-Dev)/P))*P;//cm

    //(0.016)
    if(MuPInitial.Perp(B_field) > 0.02 && track_L > 30){
      Res = sqrt(pow(sqrt(720/(((muon_z)/.3)+4))*(.3/pow(muon_z,2)),2)+
		 pow(sqrt(track_L/14)*(0.016*sqrt(pow(.105,2)+pow(MuPInitial.Mag(),2)))/(track_L*pow(MuPInitial.Perp(B_field),2)),2))*(MuPInitial.Perp(B_field)/(0.003*B_field.Mag()));

      //   Res = sqrt(720/(((muon_z)/.3)+4))*(.3/pow(muon_z,2))*(MuPInitial.Perp(B_field)/(0.003*B_field.Mag()));

      smeared_energy = smear->Gaus(energy, sqrt( pow(Res,2) + pow((Syst*Res),2))*energy);
      HowMeasuresed = 3;
    } 
    else if(MuPInitial.Perp(B_field) < 0.02 && track_L > 100){
      Res = /*-0.100989*log(0.0004644*track_L)*/-0.102*log(0.000612*track_L);
      smeared_energy = smear->Gaus(energy, sqrt( pow(Res,2) + pow((Syst*Res),2))*energy);
      HowMeasuresed = 2;
    }
    else if(muon_z > fabs(MuXYZInitial.Z()-MuXYZExit.Z()) && track_L > 100){ 
      Res = /*-0.100989*log(0.0004644*track_L)*/-0.102*log(0.000612*track_L);
      smeared_energy = smear->Gaus(energy, sqrt( pow(Res,2) + pow((Syst*Res),2))*energy);
      HowMeasuresed = 2;
    }
    if(Res > /*-0.100989*log(0.0004644*track_L)*/-0.102*log(0.000612*track_L) && track_L < 100){smeared_energy = 0; Res = 0; HowMeasuresed = 0;}
    else if(muon_z > fabs(MuXYZInitial.Z()-MuXYZExit.Z()) && track_L < 100){smeared_energy = 0; Res = 0; HowMeasuresed = 0;}
    
   
    if(hit_back_of_ranger){//This means it IS measurable! 
      Double_t Res2 = 0;
      TVector3 B_field2(0,-1.5,0);
      Double_t Dev2 = .3;//m
      Double_t P2 = MuPExit.Perp(B_field2)/(.3*B_field2.Mag());
      Double_t muon_z2 = 100*TMath::Sin(TMath::ACos((P2-Dev2)/P2))*P2;//cm

      if(MuPExit.Perp(B_field2) > 0.1){
	
      Res2 = sqrt(pow(sqrt(720/(((muon_z2)/5)+4))*(1/pow(muon_z2,2)),2)+
		  pow(sqrt(track_L/1.76)*(0.016*sqrt(pow(.105,2)+pow(MuPExit.Mag(),2)))/(muon_z2*pow(MuPExit.Perp(B_field2),2)),2))*(MuPExit.Perp(B_field2)/(0.003*B_field2.Mag()));
      
	if(Res == 0){
	  smeared_energy = smear->Gaus(energy, sqrt( pow(Res2,2) + pow((Syst*Res2),2))*energy); 
	  HowMeasuresed = 4;
	
	}
	else{ //Resolutions are measured in both places!
	if(Res2 < Res){
	    smeared_energy = smear->Gaus(energy, sqrt( pow(Res2,2) + pow((Syst*Res2),2))*energy);
	    HowMeasuresed = 4;
	  }
	  else{smeared_energy = smear->Gaus(energy, sqrt( pow(Res,2) + pow((Syst*Res),2))*energy); HowMeasuresed = 3;}
	}	
	
	
      }//needs to be a pT = 100MeV muons(!)
    } //needs to hit back of ranger
    

    //    if(smeared_energy == 0){ std::cout << " Failed to be measured! " << std::endl;
    
  }
  
  else{ 
    //Assumes resolution function from ICARUS
    //NIM paper M.P.S for exiting muons
    //Wolfram this:
    //"fit[{50,.35},{100,.29},{150,.25},{200,.22},{250,.18}] logarithmic function"
    if(track_L > 100){
      Res = /*-0.100989*log(0.0004644*track_L)*/-0.102*log(0.000612*track_L);
      smeared_energy = smear->Gaus(energy, sqrt( pow(Res,2) + pow((Syst*Res),2))*energy);
      HowMeasuresed = 2;}
    else{smeared_energy = 0; HowMeasuresed = 0;}

    //    if(fabs(smeared_energy-energy)/energy < .01) std::cout << "Mp.S. :::::   Res: " << Res << ", Smear :" << fabs(smeared_energy-energy)/energy << ", Track L : " << track_L << ", True : " << energy << ", smeared E : " << smeared_energy  << std::endl;
    


    
  }

  delete smear;
  
  return smeared_energy; 
  
}//



//###########################################################################################





//  LocalWords:  xmax

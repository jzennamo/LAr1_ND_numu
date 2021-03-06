#ifndef UTILS_HH
#define UTILS_HH

#include "FluxRW/FluxRW.cxx"
#include "FluxRWnubar/FluxRWnubar.cxx"
#include "Rtypes.h"
#include "TVector3.h"
#include <vector>
#include "TDecompLU.h"
#include "TRandom3.h"

class TH1D;

class Utils {
  
 public:
  
  Utils();
  ~Utils();
  
  Double_t GetFluxWeight( Double_t energy, Int_t iflux, Int_t inno, Int_t ndecay );
  Double_t GetTheta( TVector3 mom, TVector3 ref = TVector3(0,0,1) );
  Double_t GetTheta( Double_t px, Double_t py, Double_t pz, TVector3 ref = TVector3(0,0,1) );
  Double_t GetPhi( TVector3 mom, TVector3 ref = TVector3(0,0,1) );
  Double_t GetPhi( Double_t px, Double_t py );
  Double_t NuEnergyCCQE( Double_t l_energy, Double_t l_p, Double_t l_theta, Double_t mass, Int_t mode );
  Double_t NuEnergyCalo( std::vector<Int_t> *pdg, std::vector<Double_t> *energy, 
			 Bool_t include_neutrons = false, Bool_t include_pizeros = false, Double_t prot_thresh = 0 );
  Double_t NuEnergyCalo_smeared(  std::vector<Int_t> *pdg, std::vector<Double_t> *energy,
				  Double_t prot_thresh, double smeared_lep_E, double true_E, std::string LEP, double HadronicRes );
  Double_t VertexEnergy( std::vector<Int_t> *pdg, std::vector<Double_t> *energy, 
			 Double_t prot_thresh = 0.0, Double_t pion_thresh = 0.0 );
  Double_t TotalPhotonEnergy( Int_t idet, Int_t ND_iDet, std::vector<gan::LorentzVectorLight> *p1pos,
			      std::vector<gan::LorentzVectorLight> *p1mom,
			      std::vector<gan::LorentzVectorLight> *p2pos,
			      std::vector<gan::LorentzVectorLight> *p2mom );
  Bool_t   IsFiducial( Int_t idet, Int_t ND_iDet, TVector3 vtx);
  Bool_t   IsActive( Int_t idet, Int_t ND_iDet, TVector3 vtx, double cut = 0);
  Double_t GetFidMass( Int_t idet, Int_t ND_iDet);
  Double_t GetActiveMass( Int_t idet, Int_t ND_iDet);
  void     GetDetBoundary( Int_t idet, Int_t ND_iDet, Double_t &xmin, Double_t &xmax, 
			     Double_t &ymin, Double_t &ymax, Double_t &zmin, Double_t &zmax ); 
  Double_t GetPOTNorm( Int_t iflux, Int_t iDet );
  //  Bool_t MuonExit();
  //This function returns true if the line along direction startDir that starts at startPos
  //intersects the parallelogram defined by the 3 points given.
  //NB a parallelogram needs 4 points.  The last point is defined as p1 + a +b, where a and b
  //are the vectors from p1 to p2, p1 to p3.  
  Bool_t   IntersectsPlane(const TVector3 & startPos, const TVector3 & startDir,
                           const TVector3 & planeCorner1,
                           const TVector3 & planeCorner2,
                           const TVector3 & planeCorner3,
                           const bool verbose = false) const;

 bool RangeOut(TVector3 MuXYZ, TVector3 MuP, Double_t track_L, int ND_iDet);
 Double_t EnergyRes(bool contained, Double_t  energy, Double_t  track_L, int  ND_iDet, bool hit_back_of_ranger, 
		    TVector3  MuPExit,   TVector3  MuPInitial, 
		    TVector3  MuXYZExit, TVector3  MuXYZInitial, double Syst, int & HowMeasuresed);

 private:

  FluxRW      fluxrw_nu;
  FluxRWnubar fluxrw_nubar;

  double PotNormNubar;
  double PotNormNu;

  const static Int_t kNu         = 0;
  const static Int_t kNu_Fosc 	 = 1;
  const static Int_t kNubar 	 = 2;
  const static Int_t kNubar_Fosc = 3; 

  const static Int_t kND = 0;
  const static Int_t kuB = 1;
  const static Int_t kFD = 2; 
  const static Int_t kIC = 3; 
  const static Int_t kND_long = 4;
  const static Int_t kND_150 = 5;
  const static Int_t kND_200 = 6;
  const static Int_t kIC_offaxis = 7;


  //  double fidCut_x, fidCut_y, fidCut_z;
  double fidCut_xmax, fidCut_xmin, fidCut_xAPA, fidCut_ymax, fidCut_ymin, fidCut_zmax, fidCut_zmin;

  double nd_xmin, nd_xmax, nd_ymin, nd_ymax, nd_zmin, nd_zmax;
  double ub_xmin, ub_xmax, ub_ymin, ub_ymax, ub_zmin, ub_zmax;
  double fd_xmin, fd_xmax, fd_ymin, fd_ymax, fd_zmin, fd_zmax;
  double ic_xmin, ic_xmax,ic_ymin, ic_ymax, ic_zmin, ic_zmax;
  double nd_long_xmin, nd_long_xmax,  nd_long_ymin, nd_long_ymax, nd_long_zmin, nd_long_zmax;

};
#endif

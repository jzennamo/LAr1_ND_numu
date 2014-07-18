void Chi2_Surface(){

  //  float m2 = 0.3, float sin2 = 0.01;

  TFile f("output/chi2_Surface_LAr1ND_100m_T600_Shape_and_Rate.root", "READ");

  TNtuple * thisNtuple;

  f.GetObject("chi2", thisNtuple);

  Float_t dm2_temp;
  Float_t sin22th_temp;
  Float_t chi2_temp;

  int npoints = 500;
  // configure the dm2, sin22th points:
  double dm2min(0.01), dm2max(100.0);
  double sin22thmin(0.0001), sin22thmax(1.0);  

  thisNtuple->SetBranchAddress("chisq", &chi2_temp);
  thisNtuple->SetBranchAddress("dm2", &dm2_temp);
  thisNtuple->SetBranchAddress("sin22th", &sin22th_temp);

  int i_entry = 0;
  int i_dm2 = 0;
  int i_sin22th = 0;

  TH2D * chi2map = new TH2D("c","",  npoints, dm2min,dm2max, npoints, sin22thmin, sin22thmax);

  std::cout<< "Entries : " << thisNtuple->GetEntries() << std::endl;

  while (i_entry < thisNtuple->GetEntries())
    {

      thisNtuple->GetEntry(i_entry);

      //      double DM = pow(10.,(TMath::Log10(dm2min)+(dm2_temp*1./npoints)*TMath::Log10(dm2max/dm2min)));
      //      double Theta = pow(10.,(TMath::Log10(sin22thmin)+(sin22th_temp*1./npoints)*TMath::Log10(sin22thmax/sin22thmin)));

      chi2map->SetBinContent(i_dm2+1, i_sin22th+1, chi2_temp);

      std::cout << "dm2 : " << i_dm2 << ", sin2 : " << i_sin22th << ", chi2 : " << (chi2_temp - chi2map->GetBinContent(i_dm2+1, i_sin22th+1)) << std::endl;

      i_sin22th++;

      if (i_sin22th % (npoints+1) == 0){
	i_sin22th = 0;
	i_dm2++;
      }
      i_entry++;
    }
  
  
  /*  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(0000);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetMarkerStyle(21);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetPalette(51,0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleYSize(0.06);*/
  //TGaxis::SetMaxDigits(3);
  
  TCanvas* c1 = new TCanvas("c1", "", 800, 800);     
  c1->SetLeftMargin(.15);
  c1->SetBottomMargin(.15);
  c1->SetTopMargin(.05);
  c1->SetRightMargin(.05);
  c1->SetLogy();
  c1->SetLogx();
  c1->cd();



  chi2map->Draw("CONTZ");






}

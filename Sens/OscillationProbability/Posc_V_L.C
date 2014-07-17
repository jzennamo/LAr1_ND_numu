void Posc_V_L( float m2 = 0.3, float sin2 = 0.01, float energy = 0.7 ){

  gStyle->SetOptStat(0000);
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
  gStyle->SetTitleYSize(0.06);
  //TGaxis::SetMaxDigits(3);
  
  TCanvas* c1 = new TCanvas("c1","",800,500);     
  c1->SetLeftMargin(.15);
  c1->SetBottomMargin(.15);
  c1->SetTopMargin(.05);
  c1->SetRightMargin(.05);
  c1->cd();

  TH1 *hr = new TH2F("hr","",100,0,900,100,0,10);
  hr->SetStats(0);
  hr->GetXaxis()->SetTitle("Length of Neutrino Flight [m]");
  hr->GetXaxis()->SetLabelFont(42);
  hr->GetXaxis()->SetLabelSize(0.06);
  hr->GetXaxis()->SetTitleSize(0.06);
  hr->GetXaxis()->SetTitleOffset(1.1);
  hr->GetXaxis()->CenterTitle();
  hr->GetYaxis()->SetTitle("P(#nu_{#mu} #rightarrow #nu_{s}) [%]");
  hr->GetYaxis()->SetNdivisions(506);
  hr->GetXaxis()->SetNdivisions(506);
  hr->GetYaxis()->SetLabelFont(42);
  hr->GetYaxis()->SetLabelOffset(0.01);
  hr->GetYaxis()->SetLabelSize(0.05);
  hr->GetYaxis()->SetTitleSize(0.07);
  hr->GetYaxis()->SetTitleOffset(0.8);
  hr->GetYaxis()->CenterTitle();
  hr->GetYaxis()->SetRangeUser(0,200*sin2);
  hr->Draw();


  TF1 *Posc = new TF1("Posc", TString::Format("100*((%f)*pow(TMath::Sin((1.267)*(%f)*(.001*x/%f)),2))",sin2,m2,energy), 0, 900);
  Posc->Draw("same");

  TLatex *  dM = new TLatex(0.43, 0.82,TString::Format("#Delta m_{41}^{2} = %2.2f eV^{2}",m2));
  dM->SetNDC();
  dM->SetTextFont(42);
  dM->SetTextSize(0.05);
  TLatex *  sins = new TLatex(0.43, 0.74,TString::Format("sin^{2}(2 #theta_{#mu#mu}) = %1.3f",sin2));
  sins->SetNDC();
  sins->SetTextFont(42);
  sins->SetTextSize(0.05);
  dM->Draw();
  sins->Draw();
  TLatex *  beam = new TLatex(0.43, 0.9,TString::Format("Neutrino Energy: %4.0f MeV",1000*energy));
  beam->SetNDC();
  beam->SetTextFont(42);
  beam->SetTextSize(0.05);
  beam->Draw();

  }

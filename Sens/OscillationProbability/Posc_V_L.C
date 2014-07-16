void Posc_V_L(){

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
  
  
  TCanvas* c1 = new TCanvas("c1","",300,300,600,700);     
  c1->SetLeftMargin(.15);
  c1->SetBottomMargin(.15);
  c1->SetTopMargin(.05);
  c1->SetRightMargin(.05);
  c1->cd();

						 
  TH1 *hr = new TH2F("hr","",2,0,900,2,0,1);
  hr->SetStats(0);
  hr->GetXaxis()->SetTitle("Length of Neutrino Flight [m]");
  hr->GetXaxis()->SetLabelFont(42);
  hr->GetXaxis()->SetLabelSize(0.05);
  hr->GetXaxis()->SetTitleSize(0.06);
  hr->GetXaxis()->SetTitleOffset(1.1);
  hr->GetYaxis()->SetTitle("P(#nu_{#mu} #rightarrow #nu_{s}) [%]");
  hr->GetYaxis()->SetNdivisions(506);
  hr->GetXaxis()->SetNdivisions(506);
  hr->GetYaxis()->SetLabelFont(42);
  hr->GetYaxis()->SetLabelOffset(0.014);
  hr->GetYaxis()->SetLabelSize(0.05);
  hr->GetYaxis()->SetTitleSize(0.07);
  hr->Draw("");


  TF1 *Posc = new TF1("Posc", "100*((0.01)*pow(TMath::Sin((1.267)*(0.3)*(.001*x/0.8)),2))", 0, 900);
  Posc->Draw("same");

  TLatex *  dM = new TLatex(0.43, 0.82,"#Delta m_{41}^{2} = 0.3 eV^{2}");
  dM->SetNDC();
  dM->SetTextFont(42);
  dM->SetTextSize(0.05);
  TLatex *  sins = new TLatex(0.43, 0.74,"sin^{2}(2 #theta_{#mu#mu}) = 0.01");
  sins->SetNDC();
  sins->SetTextFont(42);
  sins->SetTextSize(0.05);
  dM->Draw();
  sins->Draw();
  TLatex *  beam = new TLatex(0.43, 0.9,"Booster Neutrino Beam");
  beam->SetNDC();
  beam->SetTextFont(42);
  beam->SetTextSize(0.05);
  beam->Draw();


  }



void Total_Rate(){

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
  c1->SetTopMargin(.1);
  c1->SetRightMargin(.1);
  c1->cd();

  TH1 *hr = new TH2F("hr","",100,0,2192,100,1000,25000000);
  hr->SetStats(0);
  hr->GetXaxis()->SetTitle("Days Since MicroBooNE Starts Data Taking");
  hr->GetXaxis()->SetLabelFont(62);
  hr->GetXaxis()->SetTitleFont(62);
  hr->GetXaxis()->SetLabelSize(0.05);
  hr->GetXaxis()->SetTitleSize(0.05);
  hr->GetXaxis()->SetTitleOffset(1.3);
  hr->GetXaxis()->CenterTitle();
  hr->GetYaxis()->SetTitle("Accumulated #nu#lower[0.3]{#mu} Events");
  hr->GetYaxis()->SetNdivisions(506);
  hr->GetXaxis()->SetNdivisions(512);
  hr->GetYaxis()->SetLabelFont(62);
  hr->GetYaxis()->SetTitleFont(62);
  hr->GetYaxis()->SetLabelOffset(0.01);
  hr->GetYaxis()->SetLabelSize(0.05);
  hr->GetYaxis()->SetTitleSize(0.05);
  hr->GetYaxis()->SetTitleOffset(1.1);
  hr->GetYaxis()->CenterTitle();
  hr->Draw();
  

  TF1 *uB = new TF1("uB", "68.2*(x)", 0, 2192);

  TF1 *ND = new TF1("ND", "1566.5*(x-1095.74)", 0, 2192);

  TF1 *IC = new TF1("IC", "353.6*(x-1095.74)", 0, 2192);


  IC->SetLineColor(kGreen+3);
  IC->Draw("same");


  ND->SetLineColor(kRed-3);
  ND->Draw("same");


  uB->SetLineColor(kBlue-3);
  uB->Draw("same");  


  TLegend *leg = new TLegend(0.2,0.7,0.5,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(62);
  leg->AddEntry(uB,"MicroBooNE","l");
  leg->AddEntry(ND,"LAr1-ND","l");
  leg->AddEntry(IC,"ICARUS T600","l");
  leg->Draw();


  }

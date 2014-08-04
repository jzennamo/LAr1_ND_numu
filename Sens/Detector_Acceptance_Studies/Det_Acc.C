void Acceptance(double Num, double N_err, double Denom, double D_err, double & out, double & err){

  if(Num == 0){ out = 0;}
  else if(Denom == 0){  out = 0;}
  else{
    out = Num/Denom;

    err = out*sqrt(pow((N_err/Num),2)+
		   pow((D_err/Denom),2));
  }
  return;
}

void Det_Acc(){

  double emin = 0.2, emax = 3.0;  int ebins = 50;
  int _rebin = 25;

  //LAr1 ND
  
  double nd_xmin = -200.0;
  double nd_xmax =  200.0;
  double nd_ymin = -242.0;
  double nd_ymax =  158.0;
  double nd_zmin =  0.0;
  double nd_zmax =  365.0;
  

  //MicroBooNE
 
  /*
  double nd_xmin =  -128.0;
  double nd_xmax =  128.0;
  double nd_ymin = -116.5;
  double nd_ymax =  116.5;
  double nd_zmin =  0.0;
  double nd_zmax =  1036.0;   
  */

  std::string temp_name = "../MatrixFiles/combined_ntuple_100m_nu_processed_numu_Joseph_Smeared.root"; 

  //TLatex *tex_ND = new TLatex(.4,.93,"MicroBooNE");
  TLatex *tex_ND = new TLatex(.375,.93,"LAr1-ND (100m)");
  tex_ND->SetNDC();
  tex_ND->SetTextFont(42);
  tex_ND->SetTextSize(0.05);

  
  TFile *temp_file = TFile::Open(temp_name.c_str());

  nuE_all = (TH1D*)(temp_file->Get("NuEnergy_all"));
  nuE_accepted = (TH1D*)(temp_file->Get("NuEnergy_accepted"));

  muE_all = (TH1D*)(temp_file->Get("MuEnergy_all"));
  muE_accepted = (TH1D*)(temp_file->Get("MuEnergy_accepted"));

  muTheta_all = (TH1D*)(temp_file->Get("MuTheta_all"));
  muTheta_accepted = (TH1D*)(temp_file->Get("MuTheta_accepted"));

  muTE_all = (TH2D*)(temp_file->Get("MuThetaEnergy_all"));
  muTE_accepted = (TH2D*)(temp_file->Get("MuThetaEnergy_accepted"));

  muLx_all = (TH1D*)(temp_file->Get("MuLx_all"));
  muLx_accepted = (TH1D*)(temp_file->Get("MuLx_accepted"));

  muLy_all = (TH1D*)(temp_file->Get("MuLy_all"));
  muLy_accepted = (TH1D*)(temp_file->Get("MuLy_accepted"));

  muLz_all = (TH1D*)(temp_file->Get("MuLz_all"));
  muLz_accepted = (TH1D*)(temp_file->Get("MuLz_accepted"));

  muE_TR = (TH2D*)(temp_file->Get("MuE_ReTr"));
  muE_TR->GetYaxis()->SetTitleSize(0.05);
  muE_TR->GetXaxis()->SetTitleSize(0.05);
  muE_TR->GetZaxis()->SetTitleSize(0.05);

  nuE_TR = (TH2D*)(temp_file->Get("nuE_ReTr"));
  nuE_TR->GetYaxis()->SetTitleSize(0.05);
  nuE_TR->GetXaxis()->SetTitleSize(0.05);
  nuE_TR->GetZaxis()->SetTitleSize(0.05);


  TH1D *NuE_Eff = new TH1D("nuE_eff",";True Neutrino Energy [GeV]; Acceptance",20,emin,emax);
  NuE_Eff->GetYaxis()->SetTitleSize(0.06);
  NuE_Eff->GetXaxis()->SetTitleSize(0.06);
  TH1D *muE_Eff = new TH1D("muE_eff",";Initial Muon Energy [GeV]; Acceptance",20,emin,emax);
  muE_Eff->GetYaxis()->SetTitleSize(0.06);
  muE_Eff->GetXaxis()->SetTitleSize(0.06);
  TH1D *muT_Eff = new TH1D("muT_eff",";Initial Muon Angle  (#theta); Acceptance",20,0,3.14);
  muT_Eff->GetYaxis()->SetTitleSize(0.06);
  muT_Eff->GetXaxis()->SetTitleSize(0.06);
  TH1D *muLx_Eff = new TH1D("muLx_eff",";Vertex X Postion [cm]; Acceptance",20,nd_xmin,nd_xmax);
  muLx_Eff->GetYaxis()->SetTitleSize(0.06);
  muLx_Eff->GetXaxis()->SetTitleSize(0.06);
  TH1D *muLy_Eff = new TH1D("muLy_eff",";Vertex Y Postion [cm]; Acceptance",20,nd_ymin,nd_ymax);
  muLy_Eff->GetYaxis()->SetTitleSize(0.06);
  muLy_Eff->GetXaxis()->SetTitleSize(0.06);
  TH1D *muLz_Eff = new TH1D("muLz_eff",";Vertex Z Postion [cm]; Acceptance",20,nd_zmin,nd_zmax);
  muLz_Eff->GetYaxis()->SetTitleSize(0.06);
  muLz_Eff->GetXaxis()->SetTitleSize(0.06);



  for( int i = 0; i <= muLz_Eff->GetNbinsX(); i++){
    double bin1 = 0, err1 = 0;
    Acceptance(muLz_accepted->GetBinContent(i), muLz_accepted->GetBinError(i),
               muLz_all->GetBinContent(i), muLz_all->GetBinError(i),
	       bin1, err1);
    muLz_Eff->SetBinContent(i,bin1); muLz_Eff->SetBinError(i,err1);
  }


  for( int i = 0; i <= muLy_Eff->GetNbinsX(); i++){
    double bin1 = 0, err1 = 0;
    Acceptance(muLy_accepted->GetBinContent(i), muLy_accepted->GetBinError(i),
               muLy_all->GetBinContent(i), muLy_all->GetBinError(i),
	       bin1, err1);
    muLy_Eff->SetBinContent(i,bin1); muLy_Eff->SetBinError(i,err1);
  }


  for( int i = 0; i <= muLx_Eff->GetNbinsX(); i++){
    double bin1 = 0, err1 = 0;
    Acceptance(muLx_accepted->GetBinContent(i), muLx_accepted->GetBinError(i),
               muLx_all->GetBinContent(i), muLx_all->GetBinError(i),
               bin1, err1);
    muLx_Eff->SetBinContent(i,bin1); muLx_Eff->SetBinError(i,err1);
  }

  for( int i = 0; i <= NuE_Eff->GetNbinsX(); i++){
    
    double bin1 = 0, err1 = 0;
    Acceptance(nuE_accepted->GetBinContent(i), nuE_accepted->GetBinError(i),
               nuE_all->GetBinContent(i), nuE_all->GetBinError(i),
               bin1, err1);
    NuE_Eff->SetBinContent(i,bin1); NuE_Eff->SetBinError(i,err1);
    
  }

  for( int i = 0; i <= muE_Eff->GetNbinsX(); i++){

    double bin1 = 0, err1 = 0;

    Acceptance(muE_accepted->GetBinContent(i), muE_accepted->GetBinError(i),
	       muE_all->GetBinContent(i), muE_all->GetBinError(i),
               bin1, err1);

    muE_Eff->SetBinContent(i,bin1); muE_Eff->SetBinError(i,err1);

  }

  for( int i = 0; i <= muT_Eff->GetNbinsX(); i++){

    double bin1 = 0, err1 = 0;

    Acceptance(muTheta_accepted->GetBinContent(i), muTheta_accepted->GetBinError(i),
               muTheta_all->GetBinContent(i), muTheta_all->GetBinError(i),
	       bin1, err1);

    muT_Eff->SetBinContent(i,bin1); muT_Eff->SetBinError(i,err1);

  }


  

  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(0000);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  //  gStyle->SetPalette(1);
  gStyle->SetPalette(51,0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  //  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleYSize(0.06);
  //  TGaxis::SetMaxDigits(3);
  
  /*  TCanvas* c4 = new TCanvas("c4","",700,700);
  c4->cd();
  NuE_Eff->SetMarkerStyle(21);
  NuE_Eff->SetMarkerSize(1.2);
  NuE_Eff->SetMaximum(1.09);
  NuE_Eff->SetMinimum(0.51);
  NuE_Eff->Draw("ep");
  tex_ND->Draw();
  c4->Print("NuE_Acc_MicroBooNE.pdf");
  //  c4->Print("NuE_eff.gif");

  TCanvas* c5 = new TCanvas("c5","",700,700);
  c5->cd();
  muE_Eff->SetMarkerStyle(21);
  muE_Eff->SetMarkerSize(1.2);
  muE_Eff->SetMaximum(1.09);
  muE_Eff->SetMinimum(0.51);
  muE_Eff->Draw("ep");
  tex_ND->Draw();
  c5->Print("MuE_Acc_MicroBooNE.pdf");
  //  c5->Print("MuE_eff.gif");

  TCanvas* c6 = new TCanvas("c6","",700,700);     
  c6->cd();
  muT_Eff->SetMarkerStyle(21);
  muT_Eff->SetMarkerSize(1.2);
  muT_Eff->SetMaximum(1.09);
  muT_Eff->SetMinimum(0.51);
  muT_Eff->Draw("ep");
  tex_ND->Draw();
  c6->Print("MuTheta_Acc_MicroBooNE.pdf");
  //  c6->Print("MuTheta_eff.gif");
  */
  TCanvas* c7 = new TCanvas("c7","",700,700);
  c7->cd();
  muLx_Eff->SetMarkerStyle(21);
  muLx_Eff->SetMarkerSize(1.2);
  muLx_Eff->SetMaximum(1.09);
  muLx_Eff->SetMinimum(0.0);
  muLx_Eff->Draw("ep");
  tex_ND->Draw();
  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kGray);
  split->DrawLine(nd_xmin+10,0,nd_xmin+10,1.09);
  split->DrawLine(nd_xmax-10,0,nd_xmax-10,1.09);
  muLx_Eff->Draw("ep same");
  c7->Print("MuLx_Acc_LAr1ND.pdf");
  //  c7->Print("MuLx_eff.gif");

  TCanvas* c8 = new TCanvas("c8","",700,700);
  c8->cd();
  muLy_Eff->SetMarkerStyle(21);
  muLy_Eff->SetMarkerSize(1.2);
  muLy_Eff->SetMaximum(1.09);
  muLy_Eff->SetMinimum(0.0);
  muLy_Eff->Draw("ep");
  tex_ND->Draw();
  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kGray);
  split->DrawLine(nd_ymin+10,0,nd_ymin+10,1.09);
  split->DrawLine(nd_ymax-10,0,nd_ymax-10,1.09); 
  muLy_Eff->Draw("ep same");
 c8->Print("MuLy_Acc_LAr1ND.pdf");
  //  c8->Print("MuLy_eff.gif");


  TCanvas* c9 = new TCanvas("c9","",700,700);
  c9->cd();
  muLz_Eff->SetMarkerStyle(21);
  muLz_Eff->SetMarkerSize(1.2);
  muLz_Eff->SetMaximum(1.09);
  muLz_Eff->SetMinimum(0.0);
  muLz_Eff->Draw("ep");
  tex_ND->Draw();
  TLine *split = new TLine();
  split->SetLineStyle(2);
  split->SetLineWidth(5);
  split->SetLineColor(kGray);
  split->DrawLine(nd_zmin+10,0,nd_zmin+10,1.09);
  split->DrawLine(nd_zmax-80,0,nd_zmax-80,1.09);
  muLz_Eff->Draw("ep same");
  c9->Print("MuLz_Acc_LAr1ND.pdf");
  //  c9->Print("MuLz_eff.gif");
  
}

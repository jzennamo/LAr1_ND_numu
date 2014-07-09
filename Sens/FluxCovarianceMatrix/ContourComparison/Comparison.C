
{
//=========Macro generated from canvas: c3/Sensitivity
//=========  (Mon Jan  6 15:48:30 2014) by ROOT version5.34/05
   TCanvas *c3 = new TCanvas("c3", "Sensitivity",616,22,700,700);
   gStyle->SetOptStat(0);
   c3->Range(-3.5625,-2.744868,0.1875,2.24973);
   c3->SetFillColor(0);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   c3->SetLogx();
   c3->SetLogy();
   c3->SetLeftMargin(0.15);
   c3->SetRightMargin(0.05);
   c3->SetTopMargin(0.05);
   c3->SetBottomMargin(0.15);
   c3->SetFrameBorderMode(0);
   c3->SetFrameBorderMode(0);
   
   TH2D *hr1__1 = new TH2D("hr1__1","",500,0.001,1,500,0.0101,100);
   hr1__1->SetDirectory(0);
   hr1__1->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   hr1__1->SetLineColor(ci);
   hr1__1->GetXaxis()->SetTitle("sin^{2}2#theta_{#mu#mu}");
   hr1__1->GetXaxis()->CenterTitle(true);
   hr1__1->GetXaxis()->SetLabelFont(42);
   hr1__1->GetXaxis()->SetLabelOffset(-0.003);
   hr1__1->GetXaxis()->SetLabelSize(0.04);
   hr1__1->GetXaxis()->SetTitleSize(0.05);
   hr1__1->GetXaxis()->SetTitleOffset(1.2);
   hr1__1->GetYaxis()->SetTitle("#Deltam_{41}^{2} [eV^{2}]");
   hr1__1->GetYaxis()->CenterTitle(true);
   hr1__1->GetYaxis()->SetLabelFont(42);
   hr1__1->GetYaxis()->SetLabelSize(0.04);
   hr1__1->GetYaxis()->SetTitleSize(0.05);
   hr1__1->GetYaxis()->SetTitleOffset(1.2);
   hr1__1->GetZaxis()->SetLabelFont(42);
   hr1__1->GetZaxis()->SetLabelSize(0.035);
   hr1__1->GetZaxis()->SetTitleSize(0.035);
   hr1__1->GetZaxis()->SetTitleFont(42);
   hr1__1->Draw("");

   tex = new TLatex(0.3,0.9,"#nu mode, CC Events");
   tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   
   gROOT->ProcessLine(".x ./Sens_Matrix_LAr1ND_100m_T600_on_axis_Shape_and_Rate.C(c3)");
   gROOT->ProcessLine(".x ./Sens_Matrix_LAr1ND_150m_T600_on_axis_Shape_and_Rate.C(c3)");
   gROOT->ProcessLine(".x ./Sens_Matrix_LAr1ND_200m_T600_on_axis_Shape_and_Rate.C(c3)");

   TLegend *leg = new TLegend(0.2,0.175,0.5,0.375,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.025);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLine *gdummy3 = new TLine();
   gdummy3->SetLineColor(kBlack);
   gdummy3->SetLineStyle(1);
   gdummy3->SetLineWidth(2);
   TLegendEntry *entry=leg->AddEntry(gdummy3,"LAr1-ND (100m), T600 (on axis): 90% CL","l");
   TLine *gdummy3 = new TLine();
   gdummy3->SetLineColor(kRed-3);
   gdummy3->SetLineStyle(1);
   gdummy3->SetLineWidth(2);
   TLegendEntry *entry=leg->AddEntry(gdummy3,"LAr1-ND (150m), T600 (on axis): 90% CL","l");
   TLine *gdummy3 = new TLine();
   gdummy3->SetLineColor(kBlue-3);
   gdummy3->SetLineStyle(1);
   gdummy3->SetLineWidth(2);
   TLegendEntry *entry=leg->AddEntry(gdummy3,"LAr1-ND (200m), T600 (on axis): 90% CL","l");

   leg->Draw();
   c3->Modified();
   c3->cd();
   c3->SetSelected(c3);

}

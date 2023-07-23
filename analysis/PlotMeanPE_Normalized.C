#include <TGraphErrors.h>
#include "TMath.h"


void PlotMeanPE_Normalized(){

  TCanvas *c1 = new TCanvas("c1", "c1", 100,100,700,500);
  //c1->SetFillColor(42);
  TCanvas *c2 = new TCanvas("c2", "c2", 150, 100, 500, 500);
  TCanvas *c3 = new TCanvas("c3", "c3", 200, 100, 500, 500);
  

  Double_t pos[5] = {-18.0,-9.0, 0.0, 9.0, 18.0};
  Double_t posErr[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  Double_t mean[5] = {13.09, 13.5, 13.84, 14.19, 14.97};
  Double_t meanErr[5] = {0.0317, 0.0325, 0.0333, 0.0338, 0.0356};
  Double_t meanNM[5] = {10.03, 10.65, 11.25, 12.01, 12.82};
  Double_t meanErrNM[5] = {0.0246, 0.0257, 0.0270, 0.0288, 0.0308};

  Int_t index = 4;

  Double_t normPoint = mean[index];
  Double_t pointErr = meanErr[index];
  Double_t meanNorm[5], meanNormErr[5], meanNormNM[5], meanNormErrNM[5];

  for(Int_t i = 0; i < 5; i++)
  {
	meanNorm[i] = mean[i]/normPoint;
	meanNormErr[i] = meanNorm[i]*TMath::Sqrt( (meanErr[i]*meanErr[i])/(mean[i]*mean[i]) + (pointErr*pointErr)/(normPoint*normPoint) );

	meanNormNM[i] = meanNM[i]/normPoint;
	meanNormErrNM[i] = meanNormNM[i]*TMath::Sqrt( (meanErrNM[i]*meanErrNM[i])/(meanNM[i]*meanNM[i]) + (pointErr*pointErr)/(normPoint*normPoint) );

  }

  c1->cd();
  c1->SetGrid();

  // With mirror
  TGraphErrors *plot = new TGraphErrors(5, pos, meanNorm, posErr, meanNormErr);
  //plot->SetTitle("Normalized Mean Photon Number vs. Position Along Detector");
  plot->SetTitle(" ");
  plot->SetMinimum(0.0); // sets minimum y value of graph display
  plot->SetMaximum(1.2); // sets maximum y value of graph display
  plot->SetMarkerStyle(20);
  plot->SetMarkerColor(4);
  plot->SetMarkerSize(1.4);
  plot->GetXaxis()->SetTitle("Position along bar (cm)");
  plot->GetYaxis()->SetTitle("Relative PE quantity produced");

  plot->Draw("ALP");

  // Without mirror (NM = No Mirror)
  TGraphErrors *plot_NM = new TGraphErrors(5, pos, meanNormNM, posErr, meanNormErrNM);
  plot_NM->SetMarkerStyle(20);
  plot_NM->SetMarkerColor(2);
  plot_NM->SetMarkerSize(1.4);
  plot_NM->Draw("LP");



  c2->cd();
  c2->SetGrid();
  plot->Draw("ALP");



  c3->cd();
  c3->SetGrid();  
  
  normPoint = meanNM[index];
  pointErr = meanErrNM[index];

  for(Int_t i = 0; i < 5; i++)
  {
	meanNormNM[i] = meanNM[i]/normPoint;
	meanNormErrNM[i] = meanNormNM[i]*TMath::Sqrt( (meanErrNM[i]*meanErrNM[i])/(meanNM[i]*meanNM[i]) + (pointErr*pointErr)/(normPoint*normPoint) );
  }

  plot_NM = new TGraphErrors(5, pos, meanNormNM, posErr, meanNormErrNM);
  plot_NM->SetTitle(" ");
  plot_NM->SetMinimum(0.0); // sets minimum y value of graph display
  plot_NM->SetMaximum(1.2); // sets maximum y value of graph display
  plot_NM->SetMarkerStyle(20);
  plot_NM->SetMarkerColor(2);
  plot_NM->SetMarkerSize(1.4);
  plot_NM->GetXaxis()->SetTitle("Position along bar (cm)");
  plot_NM->GetYaxis()->SetTitle("Relative PE quantity produced");
  plot_NM->Draw("ALP");

}


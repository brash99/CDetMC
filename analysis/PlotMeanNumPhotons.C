// Should comment what this macro does

#include <TGraphErrors.h>


void PlotMeanNumPhotons(){

  TCanvas *grc = new TCanvas("grc", "grc", 100,100,700,500);
  grc->SetGrid();
  grc->SetFillColor(42);


  // With mirror
  TGraphErrors *plot = new TGraphErrors("meanNumPhotVsDetPos");
  plot->SetTitle("Mean Num. Photons vs. Position Along Detector");
  plot->SetMinimum(0.0); // sets minimum y value of graph display
  plot->SetMaximum(16); // sets maximum y value of graph display
  plot->SetMarkerStyle(20);
  plot->SetMarkerColor(4);
  plot->SetMarkerSize(1.4);
  plot->Draw("ALP");

  // Without mirror (NM = No Mirror)
  TGraphErrors *plot_NM = new TGraphErrors("meanNumPhotVsDetPos_NM");
  plot_NM->SetMarkerStyle(20);
  plot_NM->SetMarkerColor(2);
  plot_NM->SetMarkerSize(1.4);
  plot_NM->Draw("LP");
}



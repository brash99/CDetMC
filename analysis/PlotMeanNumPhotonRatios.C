

void PlotMeanNumPhotonRatios(){

  TCanvas *grc = new TCanvas("grc", "grc", 100,100,700,500);
  grc->SetGrid();
  //grc->SetFillColor(42);
  //grc->SetFillColor(0);


  TGraph *plot_R = new TGraph("meanNumPhotonRatios");
  //plot_R->SetTitle("Ratio of Mean Num. Photons (With Mirror / Without Mirror)");
  plot_R->SetTitle(" ");
  plot_R->SetMinimum(1.0); // sets minimum y value of graph display
  plot_R->SetMaximum(1.4); // sets maximum y value of graph display
  plot_R->SetMarkerStyle(20);
  plot_R->SetMarkerColor(4);
  //plot_R->SetMarkerColor(kViolet);
  plot_R->SetMarkerSize(1.4);
  plot_R->GetXaxis()->SetTitle("Position along bar (cm)");
  plot_R->GetYaxis()->SetTitle("Ratio of mean num. PE, With/Without Mirror");

  plot_R->Draw("ALP");

  //grc->Update();
  //grc->GetFrame()->SetFillColor(42);
  //grc->GetFrame()->SetBorderSize(12);
  //grc->Modified();

}

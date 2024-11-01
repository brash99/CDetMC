#include <iostream>
#include <TF1.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TLinearFitter.h>

using namespace std;
using RNode = ROOT::RDF::RNode;

std::vector<RNode> v;
TList* myGeometryData2;
int global_run_number;

#include "AnalyseSignalsRDataFrameNoKE.h"

//RNode AnalyseSignalsRDataFrameNoKE(int run_number = 4000) {
void AnalyseSignalsRDataFrameNoKEHeaders(int run_number = 9999) {

	global_run_number = run_number;
	std::cout << run_number << std::endl;
	//TString fileName;
	//fileName.Form("data/AnaBarMC_%d.root",run_number);
	auto fileName = "data/AnaBarMC_"+std::to_string(run_number)+".root";
	auto treeName = "T";

        TFile* f = new TFile((TString)fileName,"READ");
        TTree* t = (TTree*)f->Get(treeName);

	ROOT::RDataFrame d(treeName,fileName);

        myGeometryData2 = (TList*)t->GetUserInfo()->FindObject("myGeometryData");

        //myGeometryData->Print();

	//auto entries = d.Count();
	//cout << *entries << " entries in Tree with no filter" << endl;

	auto fdf = d.Define("trigger", "getTrigger(Detector_Nhits, &Detector_id[0])")
       			.Define("fMass", "getMass(Prim_pdg)")
       			.Define("fMomentum","getMomentum(Prim_E,fMass)")
       			.Define("fPx", "getPx(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fPy", "getPy(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fPz", "getPz(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fNewTheta", "getNewTheta(fMomentum,fPy)")
       			.Define("fNewPhi", "getNewPhi(fMomentum,fPx,fPz)")
       			.Define("trigger2", "getTrigger2(trigger,fNewTheta)")
       			.Define("trigger3", "getTrigger3(trigger,fNewTheta)")
       			.Define("fingerXVec","getFingerXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)")
       			.Define("fingerYVec","getFingerYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)")
       			.Define("fingerZVec","getFingerZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)")
       			.Define("fingerTVec","getFingerTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)")
       			.Define("anaBarXVec","getAnaBarXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)")
       			.Define("anaBarYVec","getAnaBarYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)")
       			.Define("anaBarZVec","getAnaBarZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)")
       			.Define("anaBarTVec","getAnaBarTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)")
       			.Define("fingerID","getFingerID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("fingerPDG","getFingerPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("anaBarID","getAnaBarID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("anaBarPDG","getAnaBarPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
                        .Define("anaBarPMTID","getAnaBarPMTID(trigger,&PMT_Nphotons[0])") \
                        .Define("fingerPMTID","getFingerPMTID(trigger,&PMT_Nphotons[0])") \
                        .Define("fingerPMTNPhotons","getFingerPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarPMTNPhotons","getAnaBarPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarXPMT","getAnaBarXPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarZPMT","getAnaBarZPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarPXPMT","getAnaBarPXPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarPZPMT","getAnaBarPZPMT(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarPMTTimeOverThreshold","getAnaBarPMTTimeOverThreshold(trigger,&PMT_Nphotons[0],&PMT_TimeOverThreshold[0])")
       			.Define("anaBarPMTTime","getAnaBarPMTTime(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarPMTTimeTop","getAnaBarPMTTimeTop(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarPMTTimeBottom","getAnaBarPMTTimeBottom(trigger,&PMT_Nphotons[0],&PMT_Time[0])")
       			.Define("anaBarNPhotonsTotal","getAnaBarNPhotonsTotal(trigger,&PMT_Nphotons[0])")
       			.Define("imult","getAnaBarMult(trigger,&PMT_Nphotons[0])")
       			.Define("fingerEd","getFingerEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEd","getAnaBarEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEdTotal","getAnaBarEdTotal(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])");

	//auto entries2 = fdf.Count();
	//cout << *entries2 << " entries in Expanded Dataframe with no filter" << endl;

	auto triggers = fdf.Filter("trigger==true").Count();
	cout << *triggers << " entries passed Main trigger" << endl;

        auto fdft = fdf.Filter("trigger==true");
        v.push_back(fdft);

	//return fdft;
}

TCanvas* plotC1(){

  //RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

  auto hFingerX = v[0].Histo1D("fingerXVec");
  auto hFingerY = v[0].Histo1D("fingerYVec");
  auto hFingerZ = v[0].Histo1D("fingerZVec");
  auto hFingerT = v[0].Histo1D("fingerTVec");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 100,100,500,270);
  c1->Divide(2,2, 0.01, 0.01, 0);

  c1->cd(1);
  hFingerX->Draw();
  c1->cd(2);
  hFingerY->Draw();
  c1->cd(3);
  hFingerZ->Draw();
  c1->cd(4);
  hFingerT->Draw();

  c1->DrawClone();
  c1->Print("plots/c1.pdf");

  return c1;

}

TCanvas* plotC2(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hPrimE = v[0].Histo1D("Prim_E");
	auto hPrimTh = v[0].Histo1D("fNewTheta");
	auto hPrimPh = v[0].Histo1D("fNewPhi");
	auto hPrimPdg = v[0].Histo1D("Prim_pdg");

	TCanvas *c2 = new TCanvas("c2","c2",800,800);
	c2->Divide(2,2,0.01,0.01,0);

	c2->cd(1);
	hPrimE->Draw();
	c2->cd(2);
	hPrimTh->Draw();
	c2->cd(3);
	hPrimPh->Draw();
	c2->cd(4);
	hPrimPdg->Draw();

	c2->DrawClone();
	c2->Print("plots/c2RA.pdf");

	return c2;

}

void plotDetector(ROOT::RDF::RResultPtr<TH2D> hist) {

    double opacity=0.2;
    double x1 = 55.0;
    double y1 = -61.6;
    double x2 = -45.2;
    double y2 = -8.74;
    TBox *rect1 = new TBox(x1, y1, x2, y2);
    rect1->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect1);
    x1 = 55.0;
    y1 = -8.74;
    x2 = -45.2;
    y2 = 44.14;
    TBox *rect2 = new TBox(x1, y1, x2, y2);
    rect2->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect2);
    x1 = 62.5;
    y1 = -114.54;
    x2 = -37.7;
    y2 = -61.62;
    TBox *rect3 = new TBox(x1, y1, x2, y2);
    rect3->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect3);
    x1 = 70.0;
    y1 = -167.46;
    x2 = -30.2;
    y2 = -114.54;
    TBox *rect4 = new TBox(x1, y1, x2, y2);
    rect4->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect4);
    x1 = 62.5;
    y1 = 44.14;
    x2 = -37.7;
    y2 = 97.50;
    TBox *rect5 = new TBox(x1, y1, x2, y2);
    rect5->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect5);
    x1 = 70.0;
    y1 = 97.50;
    x2 = -30.2;
    y2 = 150.0;
    TBox *rect6 = new TBox(x1, y1, x2, y2);
    rect6->SetFillColorAlpha(kRed, opacity);
    hist->GetListOfFunctions()->Add(rect6);

}

void plotSinglePoints(ROOT::RDF::RResultPtr<TH2D> hist) {

    double x1 = 60.0;
    double y1 = 115.0;
    double x2 = 70.0;
    double y2 = 125.0;
    double opacity = 0.9;
    TBox *rect7 = new TBox(x1, y1, x2, y2);
    rect7->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect7);
    x1 = 60.0;
    y1 = -115.0;
    x2 = 70.0;
    y2 = -125.0;
    opacity = 0.9;
    TBox *rect8 = new TBox(x1, y1, x2, y2);
    rect8->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect8);
    x1 = -35.0;
    y1 = -5.0;
    x2 = -45.0;
    y2 = 5.0;
    opacity = 0.9;
    TBox *rect9 = new TBox(x1, y1, x2, y2);
    rect9->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect9);
    x1 = -5.0;
    y1 = 55.0;
    x2 = 5.0;
    y2 = 65.0;
    opacity = 0.9;
    TBox *rect10 = new TBox(x1, y1, x2, y2);
    rect10->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect10);
    x1 = -5.0;
    y1 = -55.0;
    x2 = 5.0;
    y2 = -65.0;
    opacity = 0.9;
    TBox *rect11 = new TBox(x1, y1, x2, y2);
    rect11->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect11);
    x1 = 25.0;
    y1 = -85.0;
    x2 = 35.0;
    y2 = -95.0;
    opacity = 0.9;
    TBox *rect12 = new TBox(x1, y1, x2, y2);
    rect12->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect12);
    x1 = 25.0;
    y1 = 85.0;
    x2 = 35.0;
    y2 = 95.0;
    opacity = 0.9;
    TBox *rect13 = new TBox(x1, y1, x2, y2);
    rect13->SetFillColorAlpha(kGreen, opacity);
    hist->GetListOfFunctions()->Add(rect13);

}


TCanvas* plotC33(){

    auto hAnaBarPMTTime_vs_ID = v[0].Histo2D({"h1", "AnaBar Time vs ID", 100, 0.0, 2500.0,100,0.0,20.0},"anaBarPMTID","anaBarPMTTime");
    auto hAnaBarXZPMT = v[0].Histo2D({"h1","AnaBar XvsZPMT",100,-100.0,100.0,100,-200,200},"anaBarXPMT","anaBarZPMT");
    auto hAnaBarPMTTime = v[0].Histo1D({"h1","AnaBar Time",100,0.0,20.0},"anaBarPMTTime");
    auto hAnaBarPMTTimeTop = v[0].Histo1D({"h1","AnaBar Time Top",100,0.0,20.0},"anaBarPMTTimeTop");
    auto hAnaBarPMTTimeBottom = v[0].Histo1D({"h1","AnaBar Time Bottom",100,0.0,20.0},"anaBarPMTTimeBottom");
    auto hAnaBarPMTTime_vs_Nphoton = v[0].Histo2D({"h1", "AnaBar Time vs Nphotons", 100, 0.0, 300.0,100,0.0,20.0},"anaBarPMTNPhotons","anaBarPMTTime");

    TCanvas* c33 = new TCanvas("c33","c33",800,800);
    c33->Divide(3,2,0.01,0.01,0);

    c33->cd(1);
    hAnaBarPMTTime_vs_ID->Draw("COLZ");
    c33->cd(2);
    hAnaBarPMTTime_vs_Nphoton->Draw("COLZ");
    c33->cd(3);
    hAnaBarPMTTime->Draw();
    hAnaBarPMTTimeTop->Draw("SAME");
    hAnaBarPMTTimeBottom->Draw("SAME");
    c33->cd(4);
    hAnaBarPMTTimeTop->Draw();
    c33->cd(5);
    hAnaBarPMTTimeBottom->Draw();
    c33->cd(6);
    hAnaBarXZPMT->Draw("COLZ");

    plotDetector(hAnaBarXZPMT);
    plotSinglePoints(hAnaBarXZPMT);

    c33->DrawClone();
    c33->Print("plots/c33.pdf");

    return c33;

}

TCanvas* plotC34(){

    auto hAnaBarXZPMT = v[0].Histo2D({"h1","AnaBar XvsZPMT",100,-100.0,100.0,100,-200,200},"anaBarXPMT","anaBarZPMT");

    TCanvas* c34 = new TCanvas("c34","c34",800,800);
    c34->Divide(1,1,0.01,0.01,0);

    c34->cd(1);
    hAnaBarXZPMT->Draw("COLZ");

    plotDetector(hAnaBarXZPMT);
    plotSinglePoints(hAnaBarXZPMT);

    c34->DrawClone();
    c34->Print("plots/c34.pdf");

    return c34;

}

TCanvas* plotC333() {

    auto hAnaBarPMTTimeOverThreshold_vs_ID = v[0].Histo2D({"h1", "AnaBar Time Over Threshold vs ID", 100, 0.0, 2500.0,100,0.0,50.0},"anaBarPMTID","anaBarPMTTimeOverThreshold");
    auto hAnaBarPMTTimeOverThreshold = v[0].Histo1D({"h1", "AnaBar Time Over Threshold", 100, 0.0, 50.0},"anaBarPMTTimeOverThreshold");
    auto hAnaBarPMTTimeTop = v[0].Histo1D({"h1", "AnaBar Time Top", 100, 0.0, 13.0},"anaBarPMTTimeTop");
    auto hAnaBarPMTTimeBottom = v[0].Histo1D({"h1", "AnaBar Time Bottom", 100, 0.0, 13.0},"anaBarPMTTimeBottom");
    auto hAnaBarPMTTimeOverThreshold_vs_Nphoton = v[0].Histo2D({"h1", "AnaBar Time Over Threshold vs Nphotons", 100, 0.0, 300.0,100,0.0,50.0},"anaBarPMTNPhotons","anaBarPMTTimeOverThreshold");
    auto hAnaBarXZPMT = v[0].Histo2D({"h1","AnaBar XvsZPMT",100,-100.0,100.0,100,-200,200},"anaBarXPMT","anaBarZPMT");

    TCanvas* c333 = new TCanvas("c33","c33",800,800);
    c333->Divide(2,2,0.01,0.01,0);

    c333->cd(1);
    hAnaBarPMTTimeOverThreshold_vs_ID->Draw("COLZ");
    c333->cd(2);
    hAnaBarPMTTimeOverThreshold_vs_Nphoton->Draw("COLZ");
    c333->cd(3);
    hAnaBarPMTTimeOverThreshold->Draw();

    c333->DrawClone();
    c333->Print("plots/c333.pdf");

    return c333;

}

TCanvas* plotC3(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hDetectorNhits = v[0].Histo1D("Detector_Nhits");
	auto hDetectorPdg = v[0].Histo1D("anaBarPDG");
	auto hDetectorID = v[0].Histo1D("anaBarID");
	auto hFingerPdg = v[0].Histo1D("fingerPDG");
	auto hFingerID = v[0].Histo1D("fingerID");
	auto hPMTID = v[0].Histo1D("PMT_id");
	auto hAnaBarPMTID = v[0].Histo1D("anaBarPMTID");
	auto hFingerPMTID = v[0].Histo1D("fingerPMTID");
	auto hAnaBarPMTTime = v[0].Histo1D("anaBarPMTTime");

	TCanvas *c3 = new TCanvas("c3","c3",800,800);
	c3->Divide(3,3,0.01,0.01,0);

	c3->cd(1);
	hDetectorNhits->Draw();
	c3->cd(2);
	hFingerPdg->Draw();
	c3->cd(3);
	hDetectorPdg->Draw();
	c3->cd(4);
	hFingerID->Draw();
	c3->cd(5);
	hDetectorID->Draw();
	c3->cd(6);
	hPMTID->Draw();
	c3->cd(7);
	hFingerPMTID->Draw();
	c3->cd(8);
	hAnaBarPMTID->Draw();
	c3->cd(9);
	hAnaBarPMTTime->Draw();

	c3->DrawClone();
	c3->Print("plots/c3RA.pdf");

	return c3;

}


TCanvas* plotC4(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hFingerEd = v[0].Histo1D("fingerEd");
	auto hFingerPMTNphot = v[0].Histo1D("fingerPMTNPhotons");
	auto hAnaBarPMTNphot = v[0].Histo1D("anaBarPMTNPhotons");
	auto hAnaBarEd = v[0].Histo1D("anaBarEd");

	TCanvas *c4 = new TCanvas("c4","c4",800,800);

	c4->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.01,0.51,0.50,0.99);
	pad1->Draw();
	pad1->cd();
	hFingerEd->GetXaxis()->SetRangeUser(1.0,10);
	hFingerEd->Draw();

	c4->cd();
	TPad *pad2 = new TPad("pad2","pad2",0.51,0.51,0.99,0.99);
	pad2->Draw();
	pad2->cd();
	hFingerPMTNphot->GetXaxis()->SetRangeUser(-10,250);
	hFingerPMTNphot->Draw();

	c4->cd();
	TPad *pad3 = new TPad("pad3","pad3",0.01,0.01,0.50,0.50);
	//pad3->SetLogy();
	pad3->Draw();
	pad3->cd();
	hAnaBarEd->GetXaxis()->SetRangeUser(1.0,10);
	hAnaBarEd->Draw();

	c4->cd();
	TPad *pad4 = new TPad("pad4","pad4",0.51,0.01,0.99,0.50);
	//pad4->SetLogy();
	pad4->Draw();
	pad4->cd();
	hAnaBarPMTNphot->GetXaxis()->SetRangeUser(-20,180);
	hAnaBarPMTNphot->Draw();

	c4->DrawClone();
	c4->Print("plots/c4RA.pdf");

	return c4;
}


TCanvas* plotC5(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hAnaBarX = v[0].Histo1D("anaBarXVec");
	auto hAnaBarY = v[0].Histo1D("anaBarYVec");
	auto hAnaBarZ = v[0].Histo1D("anaBarZVec");
	auto hAnaBarT = v[0].Histo1D("anaBarTVec");

	TCanvas *c5 = new TCanvas("c5","c5",800,800);
	c5->Divide(2,2,0.01,0.01,0);

	c5->cd(1);
	hAnaBarX->Draw();
	c5->cd(2);
	hAnaBarY->Draw();
	c5->cd(3);
	hAnaBarZ->Draw();
	c5->cd(4);
	hAnaBarT->Draw();

	c5->DrawClone();
	c5->Print("plots/c5RA.pdf");

	return c5;

}

TCanvas* plotC6(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hE1vsE2 = v[0].Histo2D({"h2", "E1 vs E2", 100, 0.01, 10.0, 100, 0.01, 30.0},"fingerEd","anaBarEdTotal");

	TCanvas *c6 = new TCanvas("c6","c6",800,800);
	c6->Divide(1,1,0.01,0.01,0);

	c6->cd(1);
	hE1vsE2->Draw("COLZ");

	c6->DrawClone();
	c6->Print("plots/c6RA.pdf");

	return c6;

}

TCanvas* plotC7(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hFinger_Edep_vs_Nphot = v[0].Filter("trigger2").Histo2D({"h3", "Finger Edep vs Nphot", 100, 0.01, 500.0, 100, 0.01, 10.0},"fingerPMTNPhotons","fingerEd");
	auto hAnaBar_Edep_vs_Nphot = v[0].Filter("trigger2").Histo2D({"h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0},"anaBarEdTotal","anaBarNPhotonsTotal");
	auto hNphot0_vs_Nphot1 = v[0].Filter("trigger2").Histo2D({"h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 500.0},"anaBarNPhotonsTotal","fingerPMTNPhotons");

	TCanvas *c7 = new TCanvas("c7","c7",800,800);
	c7->Divide(2,2,0.01,0.01,0);

	c7->cd(1);
	hFinger_Edep_vs_Nphot->Draw("COLZ");
	c7->cd(2);
	hAnaBar_Edep_vs_Nphot->Draw("COLZ");
	c7->cd(3);
	hNphot0_vs_Nphot1->Draw("COLZ");
	c7->cd(4);
	TProfile *prof = hAnaBar_Edep_vs_Nphot->ProfileX();
	prof->Fit("pol1");

	c7->DrawClone();
	c7->Print("plots/c7RA.pdf");

	return c7;

}

TCanvas* plotC8(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hFinger_Edep_vs_NphotCut = v[0].Filter("trigger3").Histo2D({"h3", "Finger Edep vs Nphot", 100, 0.01, 500.0, 100, 0.01, 10.0},"fingerPMTNPhotons","fingerEd");
	auto hAnaBar_Edep_vs_NphotCut = v[0].Filter("trigger3").Histo2D({"h4", "AnaBar Edep vs NphotTotal", 100, 0.01, 30.0, 100, 0.01, 500.0},"anaBarEdTotal","anaBarNPhotonsTotal");
	auto hNphot0_vs_Nphot1Cut = v[0].Filter("trigger3").Histo2D({"h5", "AnaBar NphotTotal vs Finger Nphot", 100, 0.01, 500.0, 100, 0.01, 500.0},"anaBarNPhotonsTotal","fingerPMTNPhotons");

	TCanvas *c8 = new TCanvas("c8","c8",800,800);
	c8->Divide(2,2,0.01,0.01,0);

	c8->cd(1);
	hFinger_Edep_vs_NphotCut->Draw("COLZ");
	c8->cd(2);
	hAnaBar_Edep_vs_NphotCut->Draw("COLZ");
	c8->cd(3);
	hNphot0_vs_Nphot1Cut->Draw("COLZ");
	c8->cd(4);
	TProfile *profCut = hAnaBar_Edep_vs_NphotCut->ProfileX();
	profCut->Fit("pol1");

	c8->DrawClone();
	c8->Print("plots/c8RA.pdf");

	return c8;

}

TCanvas* plotC11(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hAnaBarMult = v[0].Histo1D("imult");

	TCanvas *c11 = new TCanvas("c11", "c11", 800,800);
	c11->Divide(1,1, 0.01, 0.01, 0);

	c11->cd(1);
	hAnaBarMult->Draw();

	c11->DrawClone();
	c11->Print("plots/c11RA.pdf");

	return c11;

}


TCanvas* plotC12(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hPrimPx = v[0].Histo1D("fPx");
	auto hPrimPy = v[0].Histo1D("fPy");
	auto hPrimPz = v[0].Histo1D("fPz");

	TCanvas *c12 = new TCanvas("c12","c12",800,800);
	c12->Divide(2,2,0.01,0.01,0);

	c12->cd(1);
	hPrimPx->Draw();
	c12->cd(2);
	hPrimPy->Draw();
	c12->cd(3);
	hPrimPz->Draw();

	c12->DrawClone();
	c12->Print("plots/c12.pdf");

	return c12;

}

TCanvas* plotC13() {


	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hPx_vs_x = v[0].Filter("trigger2").Histo2D({"h33", "G4SBS Px vs x", 100, -800.0, 800.0, 100, -800.0, 800.0},"anaBarXVec","fPx");
	auto hPz_vs_z = v[0].Filter("trigger2").Histo2D({"h34", "G4SBS Pz vs z", 100, -2400.0, 2400.0, 100, -2400.0, 2400.0},"anaBarZVec","fPz");
	auto hz_vs_x = v[0].Filter("trigger2").Histo2D({"h35", "CDet z vs x", 100, -800.0, 800.0, 100, -2400.0, 2400.0},"anaBarXVec","anaBarZVec");
	auto hPrimXZ = v[0].Histo2D({"h99", "G4SBS z vs x", 100, -80.0, 80.0, 100, -240.0, 240.0},"Prim_X","Prim_Z");



	TCanvas *c13 = new TCanvas("c13","c13",800,800);
	c13->Divide(2,2,0.01,0.01,0);

	c13->cd(1);
	hPx_vs_x->Draw("COLZ");
	c13->cd(2);
	hPz_vs_z->Draw("COLZ");
	c13->cd(3);
	hz_vs_x->Draw("COLZ");
	c13->cd(4);
	hPrimXZ->Draw("COLZ");
    plotDetector(hPrimXZ);

	c13->DrawClone();
	c13->Print("plots/c13.pdf");

	return c13;

}


TCanvas* plotC14(){

	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hPrimX = v[0].Histo1D("Prim_X");
	auto hPrimY = v[0].Histo1D("Prim_Y");
	auto hPrimZ = v[0].Histo1D("Prim_Z");
	auto hPrimXZ = v[0].Histo2D({"h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0},"Prim_X","Prim_Z");
	
	TCanvas *c14 = new TCanvas("c14","c14",800,800);
	c14->Divide(2,2,0.01,0.01,0);

	c14->cd(1);
	hPrimZ->Draw();
	c14->cd(2);
	hPrimY->Draw();
	c14->cd(3);
	hPrimZ->Draw();
	c14->cd(4);
	hPrimXZ->Draw("COLZ");
        plotDetector(hPrimXZ);

	c14->DrawClone();
	c14->Print("plots/c14.pdf");

	return c14;

}

TCanvas* plotC15() {


	//RNode fdft = AnalyseSignalsRDataFrameNoKE(global_run_number);

	auto hPx_vs_x = v[0].Filter("trigger2").Histo2D({"h33", "Px vs x", 100, -80.0, 80.0, 100, -1.0, 1.0},"anaBarXPMT","anaBarPXPMT");
	auto hPz_vs_z = v[0].Filter("trigger2").Histo2D({"h34", "Pz vs z", 100, -240.0, 240.0, 100, -1.0, 1.0},"anaBarZPMT","anaBarPZPMT");
	auto hz_vs_x = v[0].Filter("trigger2").Histo2D({"h35", "z vs x", 100, -80.0, 80.0, 100, -240.0, 240.0},"anaBarXPMT","anaBarZPMT");
	auto hPrimXZ = v[0].Histo2D({"h99", "z vs z", 100, -80.0, 80.0, 100, -240.0, 240.0},"Prim_X","Prim_Z");



	TCanvas *c15 = new TCanvas("c15","c15",800,800);
	c15->Divide(2,2,0.01,0.01,0);

	c15->cd(1);
	hPx_vs_x->Draw("COLZ");
	c15->cd(2);
	hPz_vs_z->Draw("COLZ");
	c15->cd(3);
	hz_vs_x->Draw("COLZ");
	c15->cd(4);
	hPrimXZ->Draw("COLZ");
        plotDetector(hPrimXZ);

	c15->DrawClone();
	c15->Print("plots/c15.pdf");

	return c15;

}

TCanvas* plotC16() {


	auto hx_vs_x = v[0].Filter("trigger2").Histo2D({"h33", "X vs XPMT", 100, -80.0, 80.0, 100, -80.0, 80.0},"anaBarXPMT","Prim_X");
	auto hz_vs_z = v[0].Filter("trigger2").Histo2D({"h34", "Z vs ZPMT", 100, -240.0, 240.0, 100, -240.0, 240.0},"anaBarZPMT","Prim_Z");
	auto hPx_vs_Px = v[0].Filter("trigger2").Histo2D({"h35", "Px vs PxPMT", 100, -1.0, 1.0, 100, -800.0, 800.0},"anaBarPXPMT","fPx");
	auto hPz_vs_Pz = v[0].Filter("trigger2").Histo2D({"h99", "Pz vs PzPMT", 100, -0.6, 0.6, 100, -2400.0, 2400.0},"anaBarPZPMT","fPz");

	TCanvas *c16 = new TCanvas("c16","c16",800,800);
	c16->Divide(2,2,0.01,0.01,0);

	c16->cd(1);
	hx_vs_x->Draw("COLZ");
	c16->cd(2);
	hz_vs_z->Draw("COLZ");
	c16->cd(3);
	hPx_vs_Px->Draw("COLZ");
	c16->cd(4);
	hPz_vs_Pz->Draw("COLZ");

	c16->DrawClone();
	c16->Print("plots/c16.pdf");

	return c16;

}

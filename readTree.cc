#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iomanip>

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TApplication.h"
#include "TGraphAsymmErrors.h"

using namespace std;

int _multiWP[3] = {1,0,0};//Medium - ele, Loose - mu

void showHist(THStack *, string, string, string, double);
double calculateZbi(double signal, double bkg, double unc);

int SRID (int njets, int nbjets) {
    int index = 0;
    
    if(njets == 3 && nbjets == 0)
        index = 0;

    if(njets == 3 && nbjets == 1)
        index = 1;

    if(njets == 3 && nbjets > 1)
        index = 2;   

    if(njets > 3 && nbjets == 0)
        index = 3;

    if(njets > 3 && nbjets == 1)
        index = 4;

    if(njets > 3 && nbjets > 1)
        index = 5;
    
    /*
     if(njets > 2 && nbjets > 1)
        index = 0;
    */
    return index;

}

void readTree()
{
    
    const int nLeptonsMax = 10;
    int fontToUse = 42;
    gStyle->SetOptFit(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetNdivisions(505,"XY");
    
    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
    
    gStyle->SetLabelFont(fontToUse,"XYZ");
    gStyle->SetLabelSize(0.05,"XYZ");
    gStyle->SetLabelOffset(0.001,"X");
    
    gStyle->SetTitleFont(fontToUse,"");
    gStyle->SetTitleFont(fontToUse,"XYZ");
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetTitleSize(0.06,"XY");
    //gStyle->SetTitleXOffset(1.0);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYOffset(1.2);
    
    gStyle->SetErrorX(0.);
    
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadBottomMargin(0.12);
    //gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadLeftMargin(0.15);
    
    gStyle->SetStatFont(fontToUse);
    gStyle->SetStatColor(10);
    gStyle->SetStatFontSize(0.08);
    gStyle->SetTitleFont(fontToUse);
    gStyle->SetTitleFontSize(0.08);
    
    gStyle->SetMarkerSize(1.);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerColor(gStyle->GetHistLineColor());
    
    gStyle->SetPalette(1,0);
    
    gStyle->SetFuncColor(kRed);

    // List of branches
    TBranch        *b__eventNb;   //!
    TBranch        *b__runNb;   //!
    TBranch        *b__lumiBlock;   //!
    TBranch        *b__weight;   //!
    TBranch        *b__genqpt;   //!
    TBranch        *b__nLeptons;   //!
    TBranch        *b__lPt;   //!
    TBranch        *b__lEta;   //!
    TBranch        *b__lPhi;   //!
    TBranch        *b__lE;   //!

    TBranch        *b__nEle;   //!
    TBranch        *b__nMu;   //!
    TBranch        *b__nTau;   //!
    TBranch        *b__flavors;   //!
    TBranch        *b__charges;   //!
    TBranch        *b__indeces;   //!

    TBranch        *b__isolation;   //!
    TBranch        *b__miniisolation;   //!
    TBranch        *b__multiisolation;   //!

    TBranch        *b__ptrel;   //!
    TBranch        *b__ptratio;   //!

    TBranch        *b__origin;   //!
    TBranch        *b__originReduced;   //!

    TBranch        *b__PVchi2;   //!
    TBranch        *b__PVerr;   //!
    TBranch        *b__ipPV;   //!
    TBranch        *b__ipPVerr;   //!
    TBranch        *b__ipZPV;   //!
    TBranch        *b__ipZPVerr;   //!
    TBranch        *b__ipPVmc;   //!
    TBranch        *b__3dIP;   //!
    TBranch        *b__3dIPerr;   //!
    TBranch        *b__3dIPsig;   //!

    TBranch        *b__mt;   //!

    TBranch        *b__isloose;   //!
    TBranch        *b__istight;   //!
    TBranch        *b__chargeConst; 
    TBranch        *b__istightID;   //!

    TBranch        *b__met;   //!
    TBranch        *b__met_phi;   //!
    TBranch        *b_HT;   //!
    TBranch        *b__genmet; 
    TBranch        *b__genmet_phi;
    TBranch        *b__n_bJets;   //!
    TBranch        *b__n_Jets;   //!
    TBranch        *b__bTagged;   //!
    TBranch        *b__jetEta;   //!
    TBranch        *b__jetPhi;   //!
    TBranch        *b__jetPt;   //!
    TBranch        *b__jetE; 
    TBranch        *b__jetFlavour;
    TBranch        *b__csv;   //!
    TBranch        *b__jetDeltaR;   //!

    TBranch        *b__trigDiMuIso;
    TBranch        *b__trigMu8Ele23Iso;
    TBranch        *b__trigMu23Ele12Iso;
    TBranch        *b__trigEle23Ele12Iso;

    
    // Declaration of leaf types
    ULong64_t       _eventNb;
    ULong64_t       _runNb;
    ULong64_t       _lumiBlock;
    Double_t        _weight;
    Double_t        _genqpt;
    Int_t           _nLeptons;
    Double_t        _lPt[10];
    Double_t        _lEta[10];
    Double_t        _lPhi[10];
    Double_t        _lE[10];

    Int_t           _nEle;
    Int_t           _nMu;
    Int_t           _nTau;
    Int_t           _flavors[10];
    Double_t        _charges[10];
    Int_t           _indeces[10];

    Double_t        _isolation[10];
    Double_t        _miniisolation[10][2];
    Bool_t          _multiisolation[10][3];
    
    Double_t        _ptrel[10];
    Double_t        _ptratio[10];

    Int_t           _origin[10];
    Int_t           _originReduced[10];

    Double_t        _PVchi2;
    Double_t        _PVerr[3];
    Double_t        _ipPV[10];
    Double_t        _ipPVerr[10];
    Double_t        _ipZPV[10];
    Double_t        _ipZPVerr[10];
    Double_t        _ipPVmc[10];
    Double_t        _3dIP[10];
    Double_t        _3dIPerr[10];
    Double_t        _3dIPsig[10];

    Double_t        _mt[10];

    Bool_t          _isloose[10];
    Bool_t          _istight[10];
    Bool_t          _istightID[10];
    Bool_t          _chargeConst[10];

    Double_t        _met;
    Double_t        _met_phi;
    Double_t        HT;
    Double_t        _genmet;
    Double_t        _genmet_phi;

    Int_t           _n_bJets;
    Int_t           _n_Jets;
    Bool_t          _bTagged[20];
    Double_t        _jetEta[20];
    Double_t        _jetPhi[20];
    Double_t        _jetPt[20];
    Double_t        _jetE[20];
    Int_t           _jetFlavour[20];
    Double_t        _csv[20];
    Double_t        _jetDeltaR[20][10];
    Int_t           _n_PV;
    Double_t        _n_MCTruth_PV;

    int _hitsNumber[nLeptonsMax];

    bool _vtxFitConversion[nLeptonsMax];

    double _mvaValue[nLeptonsMax];

    bool _fromHardProcessFinalState[nLeptonsMax];
    bool _isPromptFinalState[nLeptonsMax];

    Bool_t          _trigDiMuIso;
    Bool_t          _trigMu8Ele23Iso;
    Bool_t          _trigMu23Ele12Iso;
    Bool_t          _trigEle23Ele12Iso;

    
    TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
    
    const int nSamples = 38;
    
    
    TString fileList[nSamples] = {"ttbar_2l_13.root", 
        "dy550_13.root", "dy550_1_13.root", "dy550_2_13.root", "dy550_3_13.root", "dy550_4_13.root", 
        "dy50_13.root", "dy50_1_13.root", "dy50_2_13.root", "dy50_3_13.root", "dy50_4_13.root", 
        "ttbar_2l_13.root", "singleTop.root", "singleTopBar.root",
        "dy550_13.root", "dy550_1_13.root", "dy550_2_13.root", "dy550_3_13.root", "dy550_4_13.root", 
        "dy50_13.root", "dy50_1_13.root", "dy50_2_13.root", "dy50_3_13.root", "dy50_4_13.root",   
        "wjets_13.root", "wjets_1_13.root", "wjets_2_13.root", "wjets_3_13.root", "wjets_4_13.root", "wjets_5_13.root", "wjets_6_13.root", "wjets_7_13.root",
        "WZ13.root", "ttZ13.root", "ttH13.root", "tZq13.root",
        "ttW13.root",
        "data_combine.root"
    };
    

    double xSections[nSamples] = {
        831.76 * 0.326 * 0.326, 
        71310, 224.2, 37.2, 3.581, 1.124,
        6025.2, 139.4, 42.75, 5.497, 2.21,
        831.76 * 0.326 * 0.326, 831.76*0.326*(1-0.326), 831.76*0.326*(1-0.326), 
        71310, 224.2, 37.2, 3.581, 1.124,
        6025.2, 139.4, 42.75, 5.497, 2.21,
        61526, 1347, 360, 48.9, 12.8, 5.26, 1.33, 0.03089,
        4.42965, 0.2529, 0.2, 0.0758,
        0.2043, 
        1.
    };

    //Color_t colsStack[nSamples] = {kRed, kOrange+3, kOrange+3, kOrange+3, kOrange+3, kOrange+3, kYellow, kYellow, kYellow, kYellow, kYellow, kBlue, kGreen, kGreen+3};
    Color_t colsStack[nSamples] = {kMagenta-5, kMagenta-5, kMagenta-5, kMagenta-5, kMagenta-5, kMagenta-5, kMagenta-5, kMagenta-5, kMagenta-5, kMagenta-5, kMagenta-5, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, kBlue, kGreen, kGreen+3, kMagenta, kRed, 
        kBlack
    };
    Color_t colsStackJet[3] = {kRed, kOrange, kBlue};
    

    TFile *hfile[nSamples];
    TTree* inputTreeFake[nSamples];
    
    //TString names[nSamples] = {"t#bar{t}","dy_M-5_50","dy_M-10_50_1","dy_M-10_50_2","dy_M-10_50_3","dy_M-10_50_4","dy_M-50","dy_M-50_1","dy_M-50_2","dy_M-50_3","dy_M-50_4","WZ","t#bar{t}Z","t#bar{t}W"};
    TString names[nSamples] = {"charge mis ID", "DY","dy_M-550_1","dy_M-550_2","dy_M-550_3","dy_M-550_4","dy_M-50","dy_M-50_1","dy_M-50_2","dy_M-50_3","dy_M-50_4","fakes", "ttbar_1l_Top", "ttbar_1l_TopBar", "DY","dy_M-550_1","dy_M-550_2","dy_M-550_3","dy_M-550_4","dy_M-50","dy_M-50_1","dy_M-50_2","dy_M-50_3","dy_M-50_4","WJets", "WJets_1_13", "WJets_2_13", "WJets_3_13", "WJets_4_13", "WJets_5_13", "WJets_6_13", "WJets_7_13","WZ","t#bar{t}Z", "t#bar{t}H", "tZq","t#bar{t}W", 
    "Data(2.11/fb)"
};

    const int numberBinspt = 12;
    double ptb23[numberBinspt] = {10., 20., 30.,40.,50., 60., 70., 80., 90., 100., 200., 300.}; 

    TH1D* h_leadpt[nSamples];
    TH1D* h_2ndleadpt[nSamples];
    TH1D* h_trailpt[nSamples];
    TH1D* h_sumpt[nSamples];
    TH1D* h_mll[nSamples];
    
    TH1D* h_njets[7][nSamples];
    
    TH1D* h_SR[nSamples];
    
    TH1D* h_dilep[nSamples];

    TH1D* h_ele_pt[nSamples];
    TH1D* h_ele_ptc[nSamples];

    TH1D* h_mu_pt[nSamples];
    TH1D* h_mu_ptc[nSamples];

    TGraphAsymmErrors *gh_ele_ptc[nSamples];
    TGraphAsymmErrors *gh_mu_ptc[nSamples];

    double weights[nSamples];
    
    for (int sample=0; sample!=nSamples; ++sample)  {

        weights[sample] = 0;

        TString name;
        name = Form("lead_leptons_%d",sample);
        h_leadpt[sample] = new TH1D(name, "Leading lepton p_{T} "+names[sample]+";Leading lepton p_{T} [GeV]; events / 10 GeV", 20, 0, 200);
        h_leadpt[sample]->SetLineColor(colsStack[sample]);
        h_leadpt[sample]->SetFillColor(colsStack[sample]);
        h_leadpt[sample]->SetFillStyle(1001);
        h_leadpt[sample]->SetMarkerColor(colsStack[sample]);
        //h_leadpt[sample]->SetMarkerStyle(20+sample*5);
        h_leadpt[sample]->Sumw2();
        
        name = Form("2ndlead_leptons_%d",sample);
        h_2ndleadpt[sample] = new TH1D(name, "2nd Leading lepton p_{T} "+names[sample]+";2nd Leading lepton p_{T} [GeV]; events / 10 GeV", 10, 0, 100);
        h_2ndleadpt[sample]->SetLineColor(colsStack[sample]);
        h_2ndleadpt[sample]->SetFillColor(colsStack[sample]);
        h_2ndleadpt[sample]->SetFillStyle(1001);
        h_2ndleadpt[sample]->SetMarkerColor(colsStack[sample]);
        //h_leadpt[sample]->SetMarkerStyle(20+sample*5);
        h_2ndleadpt[sample]->Sumw2();
        
        name = Form("trail_leptons_%d",sample);
        h_trailpt[sample] = new TH1D(name, "Trailing lepton p_{T} "+names[sample]+";Trailing lepton p_{T} [GeV]; events / 10 GeV", 10, 0, 100);
        h_trailpt[sample]->SetLineColor(colsStack[sample]);
        h_trailpt[sample]->SetFillColor(colsStack[sample]);
        h_trailpt[sample]->SetFillStyle(1001);
        h_trailpt[sample]->SetMarkerColor(colsStack[sample]);
        //h_trailpt[sample]->SetMarkerStyle(20+sample*5);
        h_trailpt[sample]->Sumw2();
        
        name = Form("sumpt_leptons_%d",sample);
        h_sumpt[sample] = new TH1D(name, "Sum lepton p_{T} "+names[sample]+";Sum lepton p_{T} [GeV]; events / 10 GeV", 40, 0, 400);
        h_sumpt[sample]->SetLineColor(colsStack[sample]);
        h_sumpt[sample]->SetFillColor(colsStack[sample]);
        h_sumpt[sample]->SetFillStyle(1001);
        h_sumpt[sample]->SetMarkerColor(colsStack[sample]);
        //h_trailpt[sample]->SetMarkerStyle(20+sample*5);
        h_sumpt[sample]->Sumw2();
        
        name = Form("mll_%d",sample);
        h_mll[sample] = new TH1D(name, "Invariant mass of 2 lepton "+names[sample]+";Invariant mass of 2 lepton [GeV]; events / 10 GeV", 25, 70, 120);
        h_mll[sample]->SetLineColor(colsStack[sample]);
        h_mll[sample]->SetFillColor(colsStack[sample]);
        h_mll[sample]->SetFillStyle(1001);
        h_mll[sample]->SetMarkerColor(colsStack[sample]);
        h_mll[sample]->Sumw2();
        
        for(int i = 0; i < 7; i++){
            name = Form("h_njets_%d_%d",sample, i);
            h_njets[i][sample] = new TH1D(name, "N_{jets} "+names[sample]+";N_{jets}; events / 1", 10, 0, 10);
            h_njets[i][sample]->SetLineColor(colsStack[sample]);
            h_njets[i][sample]->SetFillColor(colsStack[sample]);
            h_njets[i][sample]->SetFillStyle(1001);
            h_njets[i][sample]->SetMarkerColor(colsStack[sample]);
            //h_njets[sample]->SetMarkerStyle(20+sample*5);
            h_njets[i][sample]->Sumw2();
        }

        name = Form("h_SR_%d",sample);
        h_SR[sample] = new TH1D(name, "SR "+names[sample]+";SR; events / 1", 6, -0.5, 5.5);
        h_SR[sample]->SetLineColor(colsStack[sample]);
        h_SR[sample]->SetFillColor(colsStack[sample]);
        h_SR[sample]->SetFillStyle(1001);
        h_SR[sample]->SetMarkerColor(colsStack[sample]);

        h_SR[sample]->GetXaxis()->SetBinLabel(1, "3 jets, 0 bjets");
        h_SR[sample]->GetXaxis()->SetBinLabel(2, "3 jets, 1 bjets");
        h_SR[sample]->GetXaxis()->SetBinLabel(3, "3 jets, #geq 2 bjets");
        h_SR[sample]->GetXaxis()->SetBinLabel(4, "#geq 3 jets, 0 bjets");
        h_SR[sample]->GetXaxis()->SetBinLabel(5, "#geq 3 jets, 1 bjets");
        h_SR[sample]->GetXaxis()->SetBinLabel(6, "#geq 3 jets, #geq 2 bjets");
        //h_SR[sample]->GetXaxis()->SetLabelSize(0.04);

        
        name = Form("h_each_dilepton_%d",sample);
        h_dilep[sample] = new TH1D(name, "Di-lepton events"+names[sample]+";dilep; events / 1", 7, -0.5, 6.5);
        h_dilep[sample]->SetLineColor(colsStack[sample]);
        h_dilep[sample]->SetFillColor(colsStack[sample]);
        h_dilep[sample]->SetFillStyle(1001);
        h_dilep[sample]->SetMarkerColor(colsStack[sample]);
        h_dilep[sample]->Sumw2();

        h_dilep[sample]->GetXaxis()->SetBinLabel(1, "Total");
        h_dilep[sample]->GetXaxis()->SetBinLabel(2, "#mu^{+}#mu^{+}");
        h_dilep[sample]->GetXaxis()->SetBinLabel(3, "e^{+}#mu^{+}");
        h_dilep[sample]->GetXaxis()->SetBinLabel(4, "e^{+}e^{+}");
        h_dilep[sample]->GetXaxis()->SetBinLabel(5, "#mu^{-}#mu^{-}");
        h_dilep[sample]->GetXaxis()->SetBinLabel(6, "e^{-}#mu^{-}");
        h_dilep[sample]->GetXaxis()->SetBinLabel(7, "e^{-}e^{-}");
        //h_dilep[sample]->GetXaxis()->SetLabelSize(1);
        //h_dilep[sample]->GetXaxis()->SetLabelOffset(0.1);
        

        name = Form("electron_pt_%d", sample);
        h_ele_pt[sample] = new TH1D(name, "Electron p_{T} "+names[sample]+";Electron p_{T} (GeV); events", numberBinspt - 1, ptb23);
        h_ele_pt[sample]->SetLineColor(colsStack[sample]);
        h_ele_pt[sample]->SetFillColor(colsStack[sample]);
        h_ele_pt[sample]->SetMarkerColor(colsStack[sample]);
        h_ele_pt[sample]->Sumw2();

        name = Form("electron_pt_cut_%d", sample);
        h_ele_ptc[sample] = new TH1D(name, "Electron p_{T} cut "+names[sample]+";Electron p_{T} (GeV); events", numberBinspt - 1, ptb23);
        h_ele_ptc[sample]->SetLineColor(colsStack[sample]);
        h_ele_ptc[sample]->SetFillColor(colsStack[sample]);
        h_ele_ptc[sample]->SetMarkerColor(colsStack[sample]);
        h_ele_ptc[sample]->Sumw2();

        name = Form("muon_pt_%d", sample);
        h_mu_pt[sample] = new TH1D(name, "Muon p_{T} "+names[sample]+";Muon p_{T} (GeV); events", numberBinspt - 1, ptb23);
        h_mu_pt[sample]->SetLineColor(colsStack[sample]);
        h_mu_pt[sample]->SetFillColor(colsStack[sample]);
        h_mu_pt[sample]->SetMarkerColor(colsStack[sample]);
        h_mu_pt[sample]->Sumw2();

        name = Form("muon_pt_cut_%d", sample);
        h_mu_ptc[sample] = new TH1D(name, "Muon p_{T} cut"+names[sample]+";Muon p_{T} (GeV); events", numberBinspt - 1, ptb23);
        h_mu_ptc[sample]->SetLineColor(colsStack[sample]);
        h_mu_ptc[sample]->SetFillColor(colsStack[sample]);
        h_mu_ptc[sample]->SetMarkerColor(colsStack[sample]);
        h_mu_ptc[sample]->Sumw2();

    }
    
    THStack* st_leadpt = new THStack("st_leadpt","Leading lepton p_{T}");
    THStack* st_2ndleadpt = new THStack("st_2ndleadpt","2nd Leading lepton p_{T}");
    THStack* st_trailpt = new THStack("st_trailpt","Trailing lepton p_{T}");
    THStack* st_sumpt = new THStack("st_sumpt","Sum lepton p_{T}");
    THStack* st_mll = new THStack("st_mll","Invariant mass of 2 leptons");
    THStack* st_njets = new THStack("st_njets","N_{jets}");
    THStack* st_SR = new THStack("st_SR","SR");
    THStack* st_dilep = new THStack("st_dilep","Di-lepton events");
    
    for (int i=nSamples-3; i!=-1; --i) {
        st_leadpt->Add(h_leadpt[i]);
        //st_2ndleadpt->Add(h_2ndleadpt[i]);
        st_trailpt->Add(h_trailpt[i]);
        //st_sumpt->Add(h_sumpt[i]);
        st_mll->Add(h_mll[i]);
        st_njets->Add(h_njets[4][i]);
        st_SR->Add(h_SR[i]);
        st_dilep->Add(h_dilep[i]);
    }

    st_leadpt->Add(h_leadpt[nSamples-2]);
    st_trailpt->Add(h_leadpt[nSamples-2]);
    st_mll->Add(h_mll[nSamples-2]);
    st_njets->Add(h_njets[4][nSamples-2]);
    st_SR->Add(h_SR[nSamples-2]);
    st_dilep->Add(h_dilep[nSamples-2]);
    
    const int nVars  = 20;
    TH1D* distribs[nVars][nSamples];
    TH1D* distribs_n_1[nVars][nSamples-1];
    THStack* distribsST[nVars];
    THStack* distribsST_n_1[nVars];
    TString varN[nVars] = {"p_{T}^{leading} [GeV]","p_{T}^{trailing} [GeV]",
        "|#eta|^{max}","M_{ll}(OSSF)","M_{T} [GeV]","M_{T} [GeV]","R_{mini}",
        "E_{T}^{miss} [GeV]","H_{T} [GeV]","N_{jets}","N_{jets40}","N_{b jets}","M_{T2}","M_{T2}^{true}","#DeltaM(OSSF,Z)", "JET p_{T}^{lead}", "csv_{l}", "JET p_{T}^{sublead}", "csv_{subl}",
        "NPV"};
    double varMin[nVars] = {40,20,
        0,0,0,0,0,
        50,150,3,0,0,0,0,0, 30, 0., 30, 0.,
        0.};
    
    double varMax[nVars] = {200,100,
        2.4,200,200,200,10,
        450,750,10,10,6,200,200,100, 400, 1., 400, 1.,
        50.
    };
    
    int nBins[nVars] = {20,10,
        12,50,20,40,50,
        20,30,7,10,6,40,40,20, 37, 100, 37, 100,
        50
    };
    
    for (int i=0; i!=nVars;++i) {
        TString name = Form("varST_%d",i);
        distribsST[i] = new THStack(name,varN[i] + ";" + varN[i]);
        //distribsST[i]->GetXaxis()->SetTitle(varN[i]);
        for (int j=nSamples-1; j!=-1; --j) {
            name = Form("var_%d_%d",i,j);
            distribs[i][j] = new TH1D(name,name+";"+varN[i]+";events",nBins[i],varMin[i],varMax[i]);
            distribs[i][j]->SetLineColor(colsStack[j]);
            if (j < nSamples - 1)
                distribs[i][j]->SetFillColor(colsStack[j]);
            distribs[i][j]->SetMarkerColor(colsStack[j]);
            //distribs[i][j]->SetMarkerStyle(1);
            distribs[i][j]->SetMarkerStyle(20);
            distribs[i][j]->SetMarkerSize(0.5);
            distribs[i][j]->SetLineWidth(1);
            distribs[i][j]->Sumw2();
            if (j < nSamples - 2)
                distribsST[i]->Add(distribs[i][j]);
        }

        distribsST[i]->Add(distribs[i][nSamples-2]);
        
    }

    for (int i=0; i!=nVars;++i) {
        TString name = Form("varST_n_1_%d",i);
        distribsST_n_1[i] = new THStack(name,varN[i]);
        for (int j=nSamples-2; j!=-1; --j) {
            name = Form("var_n_1_%d_%d",i,j);
            distribs_n_1[i][j] = new TH1D(name,name+";"+varN[i],nBins[i],varMin[i],varMax[i]);
            distribs_n_1[i][j]->SetTitle(varN[i]);
            distribs_n_1[i][j]->SetLineColor(colsStack[j]);
            if (j < nSamples - 2)
                distribs_n_1[i][j]->SetFillColor(colsStack[j]);
            distribs_n_1[i][j]->SetMarkerColor(colsStack[j]);
            //distribs[i][j]->SetMarkerStyle(1);
            distribs_n_1[i][j]->SetMarkerStyle(20);
            distribs_n_1[i][j]->SetMarkerSize(0.5);
            distribs_n_1[i][j]->SetLineWidth(1);
            distribs_n_1[i][j]->Sumw2();
            if (j < nSamples - 2)
                distribsST_n_1[i]->Add(distribs_n_1[i][j]);
        }
        
    }

    //0 - vloose, 1 - vloose FOIDEmu, 2 - vlooseFOIDISOEMU, 3 - tight
    double valuesMVA[3][4];
    valuesMVA[0][0] = -0.16;
    valuesMVA[1][0] = -0.65;
    valuesMVA[2][0] = -0.74;

    valuesMVA[0][1] = -0.70;
    valuesMVA[1][1] = -0.83; 
    valuesMVA[2][1] = -0.92; 

    valuesMVA[0][2] = -0.155;
    valuesMVA[1][2] = -0.56;
    valuesMVA[2][2] = -0.76; 

    valuesMVA[0][3] = 0.87;
    valuesMVA[1][3] = 0.60;
    valuesMVA[2][3] = 0.17;

    //double stiso[2] = {0.0893, 0.121};
    double stiso[2] = {0.0354, 0.0646};
    
    TFile *file_dataMC = TFile::Open("/Users/ikhvastu/Desktop/CERN/pileupCorr/dataMC.root","READ");
    TH1D *h_dataMC = (TH1D*)file_dataMC->Get("h3");    

    ofstream myEvents;
    myEvents.open("myEvents.txt");
    
    for (int sam = 0; sam != nSamples; ++sam) {

        //if(sam != 24 && sam != 20) continue;
        //if(sam == 1) continue;
        //if(sam > 24) 
        //if(sam > 13 && sam < 24) continue;
        TString addString;
        if(sam < 37)
            addString = "addTrueNPV/";
        else
            addString = "";

        hfile[sam] = new TFile("/Users/ikhvastu/Desktop/CERN/MCsamples/" + addString + fileList[sam],"read");
        /*
        else
            hfile[sam] = new TFile("data/"+fileList[sam],"read");
            */
        
        hfile[sam]->cd("FakeElectrons");
        inputTreeFake[sam] = static_cast<TTree*>(hfile[sam]->Get("FakeElectrons/fakeTree"));
        
        inputTreeFake[sam]->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
        inputTreeFake[sam]->SetBranchAddress("_runNb", &_runNb, &b__runNb);
        inputTreeFake[sam]->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
        inputTreeFake[sam]->SetBranchAddress("_weight", &_weight, &b__weight);
        if(sam < nSamples-1){
            inputTreeFake[sam]->SetBranchAddress("_genqpt", &_genqpt, &b__genqpt);
            inputTreeFake[sam]->SetBranchAddress("_originReduced", _originReduced, &b__originReduced);
            inputTreeFake[sam]->SetBranchAddress("_isPromptFinalState", &_isPromptFinalState);
            inputTreeFake[sam]->SetBranchAddress("_fromHardProcessFinalState", &_fromHardProcessFinalState);
            inputTreeFake[sam]->SetBranchAddress("_n_MCTruth_PV", &_n_MCTruth_PV);

        }
        inputTreeFake[sam]->SetBranchAddress("_nLeptons", &_nLeptons, &b__nLeptons);
        inputTreeFake[sam]->SetBranchAddress("_lPt", _lPt, &b__lPt);
        inputTreeFake[sam]->SetBranchAddress("_lEta", _lEta, &b__lEta);
        inputTreeFake[sam]->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
        inputTreeFake[sam]->SetBranchAddress("_lE", _lE, &b__lE);

        //inputTreeFake[sam]->SetBranchAddress("_nEle", &_nEle, &b__nEle);
        //inputTreeFake[sam]->SetBranchAddress("_nMu", &_nMu, &b__nMu);
        //inputTreeFake[sam]->SetBranchAddress("_nTau", &_nTau, &b__nTau);
        inputTreeFake[sam]->SetBranchAddress("_flavors", _flavors, &b__flavors);
        inputTreeFake[sam]->SetBranchAddress("_charges", _charges, &b__charges);
        //inputTreeFake[sam]->SetBranchAddress("_indeces", _indeces, &b__indeces);

        inputTreeFake[sam]->SetBranchAddress("_isolation", _isolation, &b__isolation);
        inputTreeFake[sam]->SetBranchAddress("_miniisolation", _miniisolation, &b__miniisolation);
        inputTreeFake[sam]->SetBranchAddress("_multiisolation", _multiisolation, &b__multiisolation);

        inputTreeFake[sam]->SetBranchAddress("_ptrel", _ptrel, &b__ptrel);
        inputTreeFake[sam]->SetBranchAddress("_ptratio", _ptratio, &b__ptratio);

        //inputTreeFake[sam]->SetBranchAddress("_origin", _origin, &b__origin);
        

        inputTreeFake[sam]->SetBranchAddress("_PVchi2", &_PVchi2, &b__PVchi2);
        inputTreeFake[sam]->SetBranchAddress("_PVerr", _PVerr, &b__PVerr);
        inputTreeFake[sam]->SetBranchAddress("_ipPV", _ipPV, &b__ipPV);
        inputTreeFake[sam]->SetBranchAddress("_ipPVerr", _ipPVerr, &b__ipPVerr);
        inputTreeFake[sam]->SetBranchAddress("_ipZPV", _ipZPV, &b__ipZPV);
        inputTreeFake[sam]->SetBranchAddress("_ipZPVerr", _ipZPVerr, &b__ipZPVerr);
        inputTreeFake[sam]->SetBranchAddress("_ipPVmc", _ipPVmc, &b__ipPVmc);
        inputTreeFake[sam]->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
        inputTreeFake[sam]->SetBranchAddress("_3dIPerr", _3dIPerr, &b__3dIPerr);
        inputTreeFake[sam]->SetBranchAddress("_3dIPsig", _3dIPsig, &b__3dIPsig);

        inputTreeFake[sam]->SetBranchAddress("_mt", _mt, &b__mt);

        //inputTreeFake[sam]->SetBranchAddress("_isloose", _isloose, &b__isloose);
        inputTreeFake[sam]->SetBranchAddress("_istight", _istight, &b__istight);
        inputTreeFake[sam]->SetBranchAddress("_chargeConst", _chargeConst, &b__chargeConst);
        
        //inputTreeFake[sam]->SetBranchAddress("_istightID", _istightID, &b__istightID);

        inputTreeFake[sam]->SetBranchAddress("_met", &_met, &b__met);
        inputTreeFake[sam]->SetBranchAddress("_met_phi", &_met_phi, &b__met_phi);

        //inputTreeFake[sam]->SetBranchAddress("HT", &HT, &b_HT);
        inputTreeFake[sam]->SetBranchAddress("_genmet", &_genmet, &b__genmet);
        inputTreeFake[sam]->SetBranchAddress("_genmet_phi", &_genmet_phi, &b__genmet_phi);

        inputTreeFake[sam]->SetBranchAddress("_n_bJets", &_n_bJets, &b__n_bJets);
        inputTreeFake[sam]->SetBranchAddress("_n_Jets", &_n_Jets, &b__n_Jets);
        //inputTreeFake[sam]->SetBranchAddress("_bTagged", _bTagged, &b__bTagged);
        inputTreeFake[sam]->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
        inputTreeFake[sam]->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
        inputTreeFake[sam]->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
        inputTreeFake[sam]->SetBranchAddress("_jetE", _jetE, &b__jetE);
        //inputTreeFake[sam]->SetBranchAddress("_jetFlavour", _jetFlavour, &b__jetFlavour);
        inputTreeFake[sam]->SetBranchAddress("_csv", _csv, &b__csv);
        inputTreeFake[sam]->SetBranchAddress("_jetDeltaR", _jetDeltaR, &b__jetDeltaR);

        inputTreeFake[sam]->SetBranchAddress("_hitsNumber", &_hitsNumber);

        inputTreeFake[sam]->SetBranchAddress("_vtxFitConversion", &_vtxFitConversion);

        inputTreeFake[sam]->SetBranchAddress("_mvaValue", &_mvaValue);

        inputTreeFake[sam]->SetBranchAddress("_trigDiMuIso", &_trigDiMuIso , &b__trigDiMuIso ); // Di-muon isolation trigger
        inputTreeFake[sam]->SetBranchAddress("_trigMu8Ele23Iso", &_trigMu8Ele23Iso , &b__trigMu8Ele23Iso ); // Di-muon isolation trigger
        inputTreeFake[sam]->SetBranchAddress("_trigMu23Ele12Iso", &_trigMu23Ele12Iso , &b__trigMu23Ele12Iso ); // Di-muon isolation trigger
        inputTreeFake[sam]->SetBranchAddress("_trigEle23Ele12Iso", &_trigEle23Ele12Iso , &b__trigEle23Ele12Iso ); // Di-muon isolation trigger

        //inputTreeFake[sam]->SetBranchAddress("_trigEmulator", &_trigEmulator);
        //inputTreeFake[sam]->SetBranchAddress("_isotrigEmulator", &_isotrigEmulator);

        inputTreeFake[sam]->SetBranchAddress("_n_PV", &_n_PV);
        
        
        inputTreeFake[sam]->SetBranchAddress("_mvaValue", &_mvaValue);

        
        _hCounter->Read("hCounter");
        Double_t scale = xSections[sam]*2110/(_hCounter->GetBinContent(1));
        //Double_t scale = xSections[sam]*5000/10000;
        
        Long64_t nEntries = inputTreeFake[sam]->GetEntries();
        std::cout<<"Entries in "<<fileList[sam]<<" "<<nEntries<<std::endl;
        std::cout<<xSections[sam]<<" "<<_hCounter->GetBinContent(1)<<" "<<scale<<std::endl;

        int allEvents = 0;
        int negEvents = 0;
        int counter = 0;

        int nLoc = 0;
        int nLocLoose = 0;
        int nLocEle = 0;
        int nLocMu = 0;
        int nLoc20 = 0;

        int posCharge = 0;
        int negCharge = 0;

        int nJLoc = 0;
        int nBLoc = 0;
        int nBLooseLoc = 0;
        int nBMedLoc = 0;
        double HTLoc = 0;
        int leptIndLoose[6];
        int leptInd[6];
        int leptInd20[6];
    
        TLorentzVector l0p4, l1p4;

        //double weight = 1;

        for (Long64_t it=0; it!=nEntries; ++it) {

            inputTreeFake[sam]->GetEntry(it);

            if(it%100000 == 0)
                cout<<'.'<<flush;

            //if(it == 10000) break;
            
            if(sam == nSamples-1){
                scale = 1.;
                _weight = 1.;
            }

            //if(_eventNb != 99516) continue;
            //if(_eventNb == 99516) cout << "good event" << endl;
            //std::cout << "New event" << std::endl;
 
            if (_nLeptons < 2) continue;
            
            nLoc = 0;
            nLocLoose = 0;
            nLocEle = 0;
            nLocMu = 0;
            nLoc20 = 0;

            posCharge = 0;
            negCharge = 0;
            if(sam == 1 && _genqpt > 100)
                continue;
            if(sam == 6 && _genqpt > 100)
                continue;
            if(sam == 14 && _genqpt > 100)
                continue;
            if(sam == 19 && _genqpt > 100)
                continue;
            if(sam == 24 && _genqpt > 100)
                continue;

            for (int i=0; i!=_nLeptons; ++i) {
                if (_flavors[i] > 1) continue;
                
                if(_3dIPsig[i] > 4) continue;

                if(_miniisolation[i][0] > 0.4) continue;

                leptIndLoose[nLocLoose] = i;
                nLocLoose++;

                //if(!_istight[i]) continue;

                //if(_originReduced[i] != 0) continue;
                //if(_lPt[i] < 20) continue;
                //bool passedMVA = true;
                /*
                if (sam == 27){

                    passedMVA = false;
                    if (TMath::Abs(_lEta[i]) < 0.8 ) {
                        passedMVA = _mvaValue[i]> valuesMVA[0][3];
                    } else if (TMath::Abs(_lEta[i]) < 1.479 ) {
                        passedMVA = _mvaValue[i]> valuesMVA[1][3];
                    } else {
                        passedMVA = _mvaValue[i]> valuesMVA[2][3];
                    }
                }
                */
              
                //if(_flavors[i] == 0 && ((_miniisolation[i][0] < 0.12) && ((_ptratio[i] > 0.80) || (_ptrel[i] > 7.2)))){
                if(_flavors[i] == 0 && _isolation[i] < stiso[TMath::Abs(_lEta[i]) > 1.479]){
                    /*
                    if(sam < 25){
                        if(_hitsNumber[i] != 0) continue;
                        if(_vtxFitConversion[i]) continue;
                        if(!_chargeConst[i]) continue;
                    }
                    */
                    if(_hitsNumber[i] != 0) continue;
                    if(_vtxFitConversion[i]) continue;
                    if(!_chargeConst[i]) continue;

                    bool passedMVA = false;
                    if (TMath::Abs(_lEta[i]) < 0.8 ) {
                        passedMVA = _mvaValue[i]> valuesMVA[0][3];
                    } else if (TMath::Abs(_lEta[i]) < 1.479 ) {
                        passedMVA = _mvaValue[i]> valuesMVA[1][3];
                    } else {
                        passedMVA = _mvaValue[i]> valuesMVA[2][3];
                    }

                    if(!passedMVA) continue;

                    leptInd[nLoc] = i;
                    nLoc++;
                    nLocEle++;
                    if(_lPt[i] > 20){
                        leptInd20[nLoc20] = i;
                        nLoc20++;
                        if(_charges[i] == 1)
                            posCharge++;
                        if(_charges[i] == -1)
                            negCharge++;
                        
                    }
                    
                    
                }
          
                
                //if(_flavors[i] == 1 && ((_miniisolation[i][0] < 0.16) && ((_ptratio[i] > 0.76) || (_ptrel[i] > 7.2)))){
                if(_flavors[i] == 1 && _isolation[i] < 0.1){
                    leptInd[nLoc] = i;
                    nLoc++;
                    nLocMu++;
                    if(_lPt[i] > 20){
                        leptInd20[nLoc20] = i;
                        nLoc20++;
                        if(_charges[i] == 1)
                            posCharge++;
                        if(_charges[i] == -1)
                            negCharge++;
                        
                    }
                    
                   
                }
                
                
            }

             
            if (nLoc20 < 2) continue;
            if (posCharge < 2 && negCharge < 2) continue;
            //if (sam < 25)
            if((!_trigDiMuIso) && (!_trigMu8Ele23Iso) && (!_trigMu23Ele12Iso) && (!_trigEle23Ele12Iso)) continue; 

            
            
            double maxPt = -9999;
            int maxPtInd = -1;

            for(int i = 0; i < nLoc20; i++){
                if (_lPt[leptInd20[i]] > maxPt) {
                    maxPt = _lPt[leptInd20[i]];
                    maxPtInd = leptInd20[i];
                }
            }

            

            double max2ndPt = -9999;
            int max2ndPtInd = -1;
            for(int i = 0; i < nLoc20; i++){
                if (leptInd20[i] == maxPtInd) continue;
                if (_lPt[leptInd20[i]] > max2ndPt && _charges[leptInd20[i]] == _charges[maxPtInd]) {
                    max2ndPt = _lPt[leptInd20[i]];
                    max2ndPtInd = leptInd20[i];
                }
            }

            //std::cout << "Max and second Max Pt: " << maxPtInd << " " << max2ndPtInd << std::endl;

            if(max2ndPtInd == -1){
                maxPt = -9999;
                int maxPtIndBad = maxPtInd;
                for(int i = 0; i < nLoc20; i++){
                    if (leptInd20[i] == maxPtIndBad) continue;
                    if (_lPt[leptInd20[i]] > maxPt) {
                        maxPt = _lPt[leptInd20[i]];
                        maxPtInd = leptInd20[i];
                    }
                }

                for(int i = 0; i < nLoc20; i++){
                    if (leptInd20[i] == maxPtInd) continue;
                    if (_lPt[leptInd20[i]] > max2ndPt && _charges[leptInd20[i]] == _charges[maxPtInd]) {
                        max2ndPt = _lPt[leptInd20[i]];
                        max2ndPtInd = leptInd20[i];
                    }
                }
                

            }

            if(maxPt < 40) continue;
            //if(max2ndPt < 40) continue;


            bool fakeVetoDecision = ((_originReduced[maxPtInd] == 0 || _isPromptFinalState[maxPtInd] || _fromHardProcessFinalState[maxPtInd]) && ( _originReduced[max2ndPtInd] == 0 || _isPromptFinalState[max2ndPtInd] || _fromHardProcessFinalState[max2ndPtInd]));
            bool chargeMisIDVetoDecision = ((_originReduced[maxPtInd] != 0 && !_isPromptFinalState[maxPtInd] && !_fromHardProcessFinalState[maxPtInd]) || ( _originReduced[max2ndPtInd] != 0 && !_isPromptFinalState[max2ndPtInd] && !_fromHardProcessFinalState[max2ndPtInd]));

            //bool fakeVetoDecision = ((_originReduced[maxPtInd] == 0) && (_originReduced[max2ndPtInd] == 0));
            //bool chargeMisIDVetoDecision = ((_originReduced[maxPtInd] != 0) || (_originReduced[max2ndPtInd] != 0));

            //std::cout << "Two highest pt lepton:" << _lPt[maxPtInd] << " " << _charges[maxPtInd] << " " << _lPt[max2ndPtInd] << " " << _charges[max2ndPtInd] << std::endl; 
            //if(_charges[maxPtInd] != _charges[max2ndPtInd]) continue;
            if(sam == 13 && fakeVetoDecision)
                continue;
            if(sam == 12 && fakeVetoDecision)
                continue;
            if(sam == 11 && fakeVetoDecision)
                continue;
            if(sam < 11 && chargeMisIDVetoDecision)
                continue;

            if((sam > 13 && sam < 24) && fakeVetoDecision)
                continue;
            if((sam > 23 && sam < 32) && fakeVetoDecision)
                continue;
            //if((sam > 11 && sam < 18) && (( _originReduced[maxPtInd] == 0) && ( _originReduced[max2ndPtInd] == 0)))
            //    continue;
            
            
            //if(sam == 0 && (( _originReduced[maxPtInd] != 0 && _flavors[maxPtInd] == 0) || ( _originReduced[max2ndPtInd] != 0 && _flavors[max2ndPtInd] == 0)))
            //    continue;
            //if(sam == 1 && (( _originReduced[maxPtInd] == 0 && _flavors[maxPtInd] == 0) && ( _originReduced[max2ndPtInd] == 0 && _flavors[max2ndPtInd] == 0)))
            //    continue;
            
            //if((sam > 1 && sam < 12) && (( _originReduced[maxPtInd] != 0 && _flavors[maxPtInd] == 0) || ( _originReduced[max2ndPtInd] != 0 && _flavors[max2ndPtInd] == 0)))
            //    continue;
            //if((sam > 11 && sam < 19) && (( _originReduced[maxPtInd] == 0 && _flavors[maxPtInd] == 0) && ( _originReduced[max2ndPtInd] == 0 && _flavors[max2ndPtInd] == 0)))
            //    continue;
                
            

            int maxInd[2] = {maxPtInd, max2ndPtInd};
            
            nJLoc = 0;
            nBLoc = 0;
            nBLooseLoc = 0;
            nBMedLoc = 0;
            HTLoc = 0;
            int jetBMedLocInd[5];
            int jetBLooseLocInd[5];
            double csv_max = -9999.;
            
            for (int i=0; i!=_n_Jets; ++i) {
                //std::cout << _csv[i] << std::endl; 
                bool clean = true;

                TLorentzVector jet;
                jet.SetPtEtaPhiE(_jetPt[i], _jetEta[i], _jetPhi[i], 0.);

                //std::cout << _jetDeltaRloose[i] << std::endl;
                
                for (int j=0; j!=nLocLoose; ++j) {
                    
                    TLorentzVector lep;
                    lep.SetPtEtaPhiE(_lPt[leptIndLoose[j]], _lEta[leptIndLoose[j]], _lPhi[leptIndLoose[j]], _lE[leptIndLoose[j]]);
                   
                    //std::cout << "Pt of jets: " << _jetPt[i] << " ; Delta R: " << _jetDeltaR[i][leptInd[j]] << std::endl;
                    //if (_jetDeltaR[i][leptInd[j]] < 0.4) {
                    if(jet.DeltaR(lep) < 0.4){
                            //removeLep = true;
                            //_jetPt[i]-=_lPt[leptInd[j]];
                        clean = false;
                        break;
                    } else clean = true;
                    
                }

                
                if (clean && _jetPt[i] > 30 && _jetEta[i] < 2.5) {
                    nJLoc++;
                    HTLoc+=_jetPt[i];
                    //std::cout << _bTagged[i] << std::endl;
                    //if (_bTagged[i])
                    if(_csv[i] > csv_max)
                        csv_max = _csv[i];
                    if(_csv[i] > 0.605){
                        jetBLooseLocInd[nBLoc] = i;
                        nBLoc++;
                    }
                    if(_csv[i] > 0.89){
                        jetBMedLocInd[nBMedLoc] = i;
                        nBMedLoc++;
                    }
                }
            }

            h_njets[0][sam]->Fill(TMath::Min(double(nJLoc),varMax[9]-1),scale*_weight);



            //std::cout << nJLoc << " " << nBLoc << std::endl;

            //if (nJLoc < 3) continue;
            //if (nBLoc < 1) continue;
            //if (HTLoc < 155) continue;

            double deltaMZ = 999999.;
            double minll = 99999.;
            int index3 = -1;
            
            for (int l0 = 0; l0<nLoc; ++l0) {
                if(leptInd[l0] == maxPtInd || leptInd[l0] == max2ndPtInd) continue;    
                l0p4.SetPtEtaPhiE(_lPt[leptInd[l0]],_lEta[leptInd[l0]],_lPhi[leptInd[l0]],_lE[leptInd[l0]]);

                for(int l1 = 0; l1 < 2; ++l1){
                    if (_charges[leptInd[l0]] != _charges[maxInd[l1]]) {
                        l1p4.SetPtEtaPhiE(_lPt[maxInd[l1]],_lEta[maxInd[l1]],_lPhi[maxInd[l1]],_lE[maxInd[l1]]);
                        l1p4+=l0p4;
                        double mdiL = l1p4.M();
                        if (_flavors[leptInd[l0]] == _flavors[maxInd[l1]] ) {
                            if (mdiL < minll) minll = mdiL;
                            if (fabs(mdiL - 91) < deltaMZ) {
                                deltaMZ = fabs(mdiL - 91);
                                index3 = leptInd[l0];
                            }
                        }
                    }
                }
            }
            
            
            //------ Filling histos --------------------------

            double ele_mll = -1.;

            l0p4.SetPtEtaPhiE(_lPt[maxInd[0]],_lEta[maxInd[0]],_lPhi[maxInd[0]],_lE[maxInd[0]]);
            l1p4.SetPtEtaPhiE(_lPt[maxInd[1]],_lEta[maxInd[1]],_lPhi[maxInd[1]],_lE[maxInd[1]]);

            ele_mll = (l0p4+l1p4).M();
            
            if(ele_mll > 120.)
                ele_mll = 119.9;
            if(ele_mll < 70.)
                ele_mll = 70.1;

            // cuts applied here
            //if (nJLoc < 3) continue;
            //if(nBLoc < 1) continue;
            //if (HTLoc < 150) continue;
            //if (deltaMZ < 15) continue;

            //if (_met < 50) continue;

            double dataMC_SF;
            if(sam < nSamples-1)
                dataMC_SF = h_dataMC->GetBinContent(h_dataMC->GetXaxis()->FindBin(_n_MCTruth_PV));
            else
                dataMC_SF = 1;

            h_SR[sam]->Fill(SRID(nJLoc,nBLoc),scale*_weight);
            distribs[11][sam]->Fill(TMath::Min(double(nBLoc),varMax[11]-1),scale*_weight);

            h_dilep[sam]->Fill(0.,scale*_weight*dataMC_SF);
            if(_charges[maxPtInd] == 1){
                
                if(_flavors[maxPtInd] == 1 &&  _flavors[max2ndPtInd] == 1)
                    h_dilep[sam]->Fill(1.,scale*_weight);
                if(_flavors[maxPtInd] == 1 &&  _flavors[max2ndPtInd] == 0)
                    h_dilep[sam]->Fill(2.,scale*_weight);
                if(_flavors[maxPtInd] == 0 &&  _flavors[max2ndPtInd] == 1)
                    h_dilep[sam]->Fill(2.,scale*_weight);
                if(_flavors[maxPtInd] == 0 &&  _flavors[max2ndPtInd] == 0)
                    h_dilep[sam]->Fill(3.,scale*_weight);
            }
            
            if(_charges[maxPtInd] == -1){
                
                if(_flavors[maxPtInd] == 1 &&  _flavors[max2ndPtInd] == 1)
                    h_dilep[sam]->Fill(4.,scale*_weight);
                if(_flavors[maxPtInd] == 1 &&  _flavors[max2ndPtInd] == 0)
                    h_dilep[sam]->Fill(5.,scale*_weight);
                if(_flavors[maxPtInd] == 0 &&  _flavors[max2ndPtInd] == 1)
                    h_dilep[sam]->Fill(5.,scale*_weight);
                if(_flavors[maxPtInd] == 0 &&  _flavors[max2ndPtInd] == 0)
                    h_dilep[sam]->Fill(6.,scale*_weight);
            }

            
            //if(nBLoc < 2) continue;

            h_mll[sam]->Fill(ele_mll, scale*_weight);
            distribs[4][sam]->Fill(TMath::Min(TMath::Min(_mt[maxPtInd], _mt[max2ndPtInd]),varMax[4]-0.1),scale*_weight);

            //if(deltaMZ > 15)
            distribs[9][sam]->Fill(TMath::Min(double(nJLoc),varMax[9]-1),scale*_weight);
            distribs[8][sam]->Fill(TMath::Min(HTLoc,varMax[8]-1),scale*_weight);
            //
            //if(deltaMZ > 15)
            

            h_njets[1][sam]->Fill(TMath::Min(double(nJLoc),varMax[9]-1),scale*_weight);
            
            if(sam < nSamples - 1){
                //distribs_n_1[15][sam]->Fill(TMath::Min(_jetPt[jetBMedLocInd[0]],varMax[15]-0.1),scale*_weight);
                distribs_n_1[16][sam]->Fill(TMath::Min(csv_max,varMax[16]-0.0001),scale*_weight);
            }

            //if (nBMedLoc < 1) continue;
            //if(nBLoc < 1) continue;
            if(nBLoc > 0){
                TLorentzVector bjet;
                //bjet.SetPtEtaPhiE(_jetPt[jetBMedLocInd[0]],_jetEta[jetBMedLocInd[0]],_jetPhi[jetBMedLocInd[0]],_jetE[jetBMedLocInd[0]]);
                bjet.SetPtEtaPhiE(_jetPt[jetBLooseLocInd[0]],_jetEta[jetBLooseLocInd[0]],_jetPhi[jetBLooseLocInd[0]],_jetE[jetBLooseLocInd[0]]);
                double deltaR_bjetLept_min = bjet.DeltaR(l0p4);
                double deltaR_bjetLept_max = bjet.DeltaR(l1p4);
                if(deltaR_bjetLept_max < deltaR_bjetLept_min)
                    deltaR_bjetLept_min = deltaR_bjetLept_max;
                distribs[6][sam]->Fill(TMath::Min(deltaR_bjetLept_min,varMax[6]-1),scale*_weight);
            }

            if(sam < nSamples - 1){
                distribs_n_1[15][sam]->Fill(TMath::Min(_jetPt[jetBLooseLocInd[0]],varMax[15]-0.1),scale*_weight);
                //distribs_n_1[16][sam]->Fill(TMath::Min(csv_max,varMax[16]-0.0001),scale*_weight);
            }
            //if(_jetPt[jetBLooseLocInd[0]] < 50) continue;

            h_njets[2][sam]->Fill(TMath::Min(double(nJLoc),varMax[9]-1),scale*_weight);

            //if (HTLoc < 150) continue;

            h_njets[3][sam]->Fill(TMath::Min(double(nJLoc),varMax[9]-1),scale*_weight);

            //

            h_njets[4][sam]->Fill(TMath::Min(double(nJLoc),varMax[9]-1),scale*_weight);

            distribs[0][sam]->Fill(TMath::Min(_lPt[maxInd[0]],varMax[0]-1),scale*_weight);
            distribs[1][sam]->Fill(TMath::Min(_lPt[maxInd[1]],varMax[0]-1),scale*_weight);
            
            

            distribs[7][sam]->Fill(TMath::Min(_met,varMax[7]-0.1),scale*_weight);

            
            //
            //if (nBLoc < 2 && nBMedLoc < 1) continue;
            h_njets[5][sam]->Fill(TMath::Min(double(nJLoc),varMax[9]-1),scale*_weight);
              
            
            //if(_jetPt[jetBLooseLocInd[0]] < )

            //if(index3 != -1 && _mt[index3] < 50) continue;
            double min_mt = TMath::Min(_mt[maxPtInd], _mt[max2ndPtInd]);
            //if(min_mt < 50) continue;
            //
            h_njets[6][sam]->Fill(TMath::Min(double(nJLoc),varMax[9]-1),scale*_weight);
            //if (nBLoc < 2 && nBMedLoc < 1) continue;

            if(sam < nSamples - 1 && nBLoc > 1){
                distribs_n_1[17][sam]->Fill(TMath::Min(_jetPt[jetBLooseLocInd[1]],varMax[15]-0.1),scale*_weight);
                distribs_n_1[18][sam]->Fill(TMath::Min(_csv[jetBLooseLocInd[1]],varMax[16]-0.0001),scale*_weight);
            }

            
            
                

            //double dataMC_SF = 1.;

            distribs[19][sam]->Fill(TMath::Min(double(_n_PV),varMax[19]-0.1),scale*_weight*dataMC_SF);
            
            allEvents++;
            if(_weight < 0)
                negEvents++;
            
            //std::cout << _runNb << " " << _lumiBlock << " " << _eventNb << "\n";
                
        }
        std::cout << "Scale and weight: " << scale << " " << _weight << std::endl; 
        std::cout << "Total nEv: " << (allEvents - 2 *negEvents) * scale * TMath::Abs(_weight) << " ; negEvents: " << negEvents * scale * TMath::Abs(_weight) << std::endl;
        cout<<endl;

        weights[sam] = scale*TMath::Abs(_weight);
        
        //gh_ele_ptc[sam] = new TGraphAsymmErrors(h_ele_ptc[sam], h_ele_pt[sam]);
        //gh_mu_ptc[sam] = new TGraphAsymmErrors(h_mu_ptc[sam], h_mu_pt[sam]);

    }
    cout<<endl;

    std::cout<<"Done"<<std::endl;
    
    TLegend* mtleg = new TLegend(0.5,0.88,0.85,0.55);
    mtleg->SetFillColor(0);
    mtleg->SetFillStyle(0);
    mtleg->SetBorderSize(0);
    for (int i=0; i!=nSamples-1; ++i) {
        //if(i == 1) continue;
        if(i > 0 && i < 11)
            continue;
        if(i > 11 && i < 32)
            continue;
        if(i > 36)
            continue;
        mtleg->AddEntry(h_leadpt[i],names[i],"f");
    }

    //mtleg->AddEntry(h_leadpt[nSamples-1],names[nSamples-1],"l");


/*
    TCanvas* allPlots = new TCanvas("allPlots","allPlots",1600,1200);
    allPlots->Divide(5,4,0.001,0.001);
    for (int i=0; i!=nVars;++i) {
        allPlots->cd(i+1);
        //distribsST[i]->SetMaximum(distribsST[i]->GetMaximum() * 10);
        distribsST[i]->Draw("hist");
        if ((i == 15) || (i == 16))
            allPlots->cd(i+1)->SetLogy();
        mtleg->Draw("same");
    }
    */

    double SFnumber = 1.6;


    TCanvas* plot = new TCanvas("ptLep","ptLep",1600,1200);
    plot->Divide(3,2);

    plot->cd(1);
    //plot->cd(1)->SetLogy();
    //distribsST[9]->SetMaximum(2000);
    showHist(distribsST[0],"p_{leading}","","",SFnumber);
    mtleg->Draw("same");
    //distribs[9][nSamples-1]->Add(distribs[9][nSamples-2]);
    //distribs[9][nSamples-1]->Add(distribs[9][nSamples-3]);
    
    distribs[0][nSamples-1]->SetLineWidth(2.);
    distribs[0][nSamples-1]->SetFillColor(0);
    distribs[0][nSamples-1]->SetLineColor(1);
    distribs[0][nSamples-1]->SetLineStyle(1);
    distribs[0][nSamples-1]->Draw("Psame");
    

    plot->cd(2);
    //plot->cd(2)->SetLogy();
    //distribsST[11]->SetMaximum(2000);
    showHist(distribsST[1],"p_{trailing}","","",SFnumber);
    mtleg->Draw("same");
    //distribs[11][nSamples-1]->Add(distribs[11][nSamples-2]);
    //distribs[11][nSamples-1]->Add(distribs[11][nSamples-3]);
    
    distribs[1][nSamples-1]->SetLineWidth(2.);
    distribs[1][nSamples-1]->SetFillColor(0);
    distribs[1][nSamples-1]->SetLineColor(1);
    distribs[1][nSamples-1]->SetLineStyle(1);
    distribs[1][nSamples-1]->Draw("Psame");
    

    plot->cd(3);
    //plot->cd(1)->SetLogy();
    //distribsST[9]->SetMaximum(2000);
    showHist(distribsST[9],"N_{jets}","","",SFnumber);
    mtleg->Draw("same");
    //distribs[9][nSamples-1]->Add(distribs[9][nSamples-2]);
    //distribs[9][nSamples-1]->Add(distribs[9][nSamples-3]);
    
    distribs[9][nSamples-1]->SetLineWidth(2.);
    distribs[9][nSamples-1]->SetFillColor(0);
    distribs[9][nSamples-1]->SetLineColor(1);
    distribs[9][nSamples-1]->SetLineStyle(1);
    distribs[9][nSamples-1]->Draw("Psame");
    

    plot->cd(4);
    //plot->cd(2)->SetLogy();
    //distribsST[11]->SetMaximum(2000);
    showHist(distribsST[11],"N_{bjets}","","",SFnumber);
    mtleg->Draw("same");
    //distribs[11][nSamples-1]->Add(distribs[11][nSamples-2]);
    //distribs[11][nSamples-1]->Add(distribs[11][nSamples-3]);
    
    distribs[11][nSamples-1]->SetLineWidth(2.);
    distribs[11][nSamples-1]->SetFillColor(0);
    distribs[11][nSamples-1]->SetLineColor(1);
    distribs[11][nSamples-1]->SetLineStyle(1);
    distribs[11][nSamples-1]->Draw("Psame");
    

    //plot->cd(5);
    //plot->cd(3)->SetLogy();
    //st_mll->SetMaximum(2000);
    //showHist(st_mll,"","","Events / 1 ",1.);
    /*
    h_mll[nSamples-1]->SetLineWidth(2.);
    h_mll[nSamples-1]->SetFillColor(0);
    h_mll[nSamples-1]->SetLineColor(1);
    h_mll[nSamples-1]->SetLineStyle(1);
    h_mll[nSamples-1]->Draw("Psame");
    */

    plot->cd(5);
    //plot->cd(4)->SetLogy();
    //distribsST[7]->SetMaximum(2000);
    showHist(distribsST[7],"#slash{E}_{T}","GeV","Events / 50 GeV",SFnumber);
    mtleg->Draw("same");
    //distribs[7][nSamples-1]->Add(distribs[7][nSamples-2]);
    //distribs[7][nSamples-1]->Add(distribs[7][nSamples-3]);
    
    distribs[7][nSamples-1]->SetLineWidth(2.);
    distribs[7][nSamples-1]->SetFillColor(0);
    distribs[7][nSamples-1]->SetLineColor(1);
    distribs[7][nSamples-1]->SetLineStyle(1);
    distribs[7][nSamples-1]->Draw("Psame");
    

      
    //st_dilep->SetMaximum(50);
    
    
    /*
    h_dilep[nSamples-1]->SetLineWidth(2.);
    h_dilep[nSamples-1]->SetFillColor(0);
    h_dilep[nSamples-1]->SetLineColor(1);
    h_dilep[nSamples-1]->SetLineStyle(1);
    h_dilep[nSamples-1]->Draw("Psame");
    */

    //plot->cd(8);  
    //showHist(distribsST[4],"{M}_{T}","GeV","Events / 50 GeV",10);
    /*
    distribs[4][nSamples-1]->SetLineWidth(2.);
    distribs[4][nSamples-1]->SetFillColor(0);
    distribs[4][nSamples-1]->SetLineColor(1);
    distribs[4][nSamples-1]->SetLineStyle(1);
    distribs[4][nSamples-1]->Draw("Psame");
    */

    //plot->cd(9);  
    //showHist(distribsST[6],"#Delta R","","Events",10);

    plot->cd(6);
    showHist(distribsST[8],"H_{T}","","Events",SFnumber);
    distribs[8][nSamples-1]->SetLineWidth(2.);
    distribs[8][nSamples-1]->SetFillColor(0);
    distribs[8][nSamples-1]->SetLineColor(1);
    distribs[8][nSamples-1]->SetLineStyle(1);
    distribs[8][nSamples-1]->Draw("Psame");
    mtleg->Draw("same");

    
    TCanvas* plotplot = new TCanvas("ptLepplot","ptLepplot",600,400);
    plotplot->Divide(2,1);
    plotplot->cd(1);
    showHist(st_dilep,"","","Events / 1 ",SFnumber);
    st_dilep->GetXaxis()->SetLabelSize(0.1);
    st_dilep->GetXaxis()->SetLabelOffset(0.01);
    h_dilep[nSamples-1]->SetLineWidth(2.);
    h_dilep[nSamples-1]->SetFillColor(0);
    h_dilep[nSamples-1]->SetLineColor(1);
    h_dilep[nSamples-1]->SetLineStyle(1);
    h_dilep[nSamples-1]->Draw("Psame");
    mtleg->Draw("same");

    plotplot->cd(2);
    showHist(distribsST[19],"norm","","Events",SFnumber);
    distribs[19][nSamples-1]->SetLineWidth(2.);
    distribs[19][nSamples-1]->SetFillColor(0);
    distribs[19][nSamples-1]->SetLineColor(1);
    distribs[19][nSamples-1]->SetLineStyle(1);
    distribs[19][nSamples-1]->DrawNormalized("Psame");
    mtleg->Draw("same");


    //plotplot->cd(2);
    TCanvas* plotplot1 = new TCanvas("ptLepplot1","ptLepplot1",600,400);
    showHist(st_SR,"SR","","Events",10);
    st_SR->GetXaxis()->SetLabelSize(0.04);
    st_SR->GetXaxis()->SetLabelOffset(0.01);
    mtleg->Draw("same");
    //plot->cd(12);
    //mtleg->Draw();
    /*
    distribs[6][nSamples-1]->SetLineWidth(2.);
    distribs[6][nSamples-1]->SetFillColor(0);
    distribs[6][nSamples-1]->SetLineColor(1);
    distribs[6][nSamples-1]->SetLineStyle(1);
    distribs[6][nSamples-1]->Draw("Psame");
    */
    TCanvas* plot1 = new TCanvas("ptLep1","ptLep1",1600,1200);
    plot1->Divide(2,2);
    plot1->cd(1);  
    showHist(distribsST_n_1[15],"norm","","Events",10);
    distribs_n_1[15][nSamples-2]->SetLineWidth(2.);
    distribs_n_1[15][nSamples-2]->SetFillColor(0);
    distribs_n_1[15][nSamples-2]->SetLineColor(1);
    distribs_n_1[15][nSamples-2]->SetLineStyle(1);
    distribs_n_1[15][nSamples-2]->DrawNormalized("Psame");
    
    plot1->cd(2);  
    showHist(distribsST_n_1[16],"norm","","Events",10);
    distribs_n_1[16][nSamples-2]->SetLineWidth(2.);
    distribs_n_1[16][nSamples-2]->SetFillColor(0);
    distribs_n_1[16][nSamples-2]->SetLineColor(1);
    distribs_n_1[16][nSamples-2]->SetLineStyle(1);
    distribs_n_1[16][nSamples-2]->DrawNormalized("Psame");
    
    plot1->cd(3);  
    showHist(distribsST_n_1[17],"norm","","Events",10);
    distribs_n_1[17][nSamples-2]->SetLineWidth(2.);
    distribs_n_1[17][nSamples-2]->SetFillColor(0);
    distribs_n_1[17][nSamples-2]->SetLineColor(1);
    distribs_n_1[17][nSamples-2]->SetLineStyle(1);
    distribs_n_1[17][nSamples-2]->DrawNormalized("Psame");
    
    plot1->cd(4);  
    showHist(distribsST_n_1[18],"norm","","Events",10);
    distribs_n_1[18][nSamples-2]->SetLineWidth(2.);
    distribs_n_1[18][nSamples-2]->SetFillColor(0);
    distribs_n_1[18][nSamples-2]->SetLineColor(1);
    distribs_n_1[18][nSamples-2]->SetLineStyle(1);
    distribs_n_1[18][nSamples-2]->DrawNormalized("Psame");
    

    
    //plot->cd(2)->SetLogy();
    //showHist(st_trailpt,"Trailing lepton p_{T}","GeV","Events / 10 GeV",10);
    

    /*
    plot->cd(9);
    plot->cd(9)->SetLogy();
    showHist(distribsST[4],"M_{T} [GeV]","","",1.2);
    distribs[4][nSamples-1]->SetLineWidth(2.);
    distribs[4][nSamples-1]->SetFillColor(0);
    distribs[4][nSamples-1]->SetLineColor(1);
    distribs[4][nSamples-1]->SetLineStyle(2);
    distribs[4][nSamples-1]->Draw("histsame");

    plot->cd(10);
    plot->cd(10)->SetLogy();
    showHist(distribsST[5],"M_{T} [GeV]","","",1.2);
    distribs[5][nSamples-1]->SetLineWidth(2.);
    distribs[5][nSamples-1]->SetFillColor(0);
    distribs[5][nSamples-1]->SetLineColor(1);
    distribs[5][nSamples-1]->SetLineStyle(2);
    distribs[5][nSamples-1]->Draw("histsame");
    */


    /*
    plot->cd(7);
    plot->cd(7)->SetLogy();
    showHist(st_sumpt,"Sum lepton p_{T}","GeV","Events / 10 GeV",10);

    plot->cd(7);
    plot->cd(7)->SetLogy();
    showHist(distribsST[15],"M_{T2ll}","GeV","Events / 10 GeV",10);

    plot->cd(8);
    plot->cd(8)->SetLogy();
    showHist(distribsST[16],"M_{T2blbl}","GeV","Events / 10 GeV",10);
    */


    //plot->cd(9);

    /*
    TCanvas * ceff = new TCanvas("ceff","ceff",900,900);
    ceff->Divide(2,1);

    ceff->cd(1);
    gh_ele_ptc[0]->GetXaxis()->SetTitle("electron pt, [GeV]");
    gh_ele_ptc[0]->GetYaxis()->SetTitle("eff");
    gh_ele_ptc[0]->SetTitle("ele eff charge non consistent prompt");
    gh_ele_ptc[0]->Draw();
    
    ceff->cd(2);
    gh_mu_ptc[0]->GetXaxis()->SetTitle("muon pt, [GeV]");
    gh_mu_ptc[0]->GetYaxis()->SetTitle("eff");
    gh_mu_ptc[0]->SetTitle("mu eff prompt");
    gh_mu_ptc[0]->Draw();
    */
    
    
    std::cout<<"Done 2"<<std::endl;

    
//____________________________________________________________________________________________________

    
    int ind[8] = {0, 11, 32, 33, 34, 35, 36, 37};

    double sum_int[7][7] = {0.};
    double sum_err[7][7] = {0.};

    for(int i = 0; i < 7; i++){
        for(int j = ind[i]; j < ind[i+1]; j++){
        
            Double_t error;
            double integ = h_njets[0][j]->IntegralAndError(1,10,error);

            sum_int[0][i] += integ;
            sum_err[0][i] += error*error;

            integ = h_njets[1][j]->IntegralAndError(1,10,error);

            sum_int[1][i] += integ;
            sum_err[1][i] += error*error;

            integ = h_njets[2][j]->IntegralAndError(1,10,error);

            sum_int[2][i] += integ;
            sum_err[2][i] += error*error;

            integ = h_njets[3][j]->IntegralAndError(1,10,error);

            sum_int[3][i] += integ;
            sum_err[3][i] += error*error;
        
            integ = h_njets[4][j]->IntegralAndError(1,10,error);

            sum_int[4][i] += integ;
            sum_err[4][i] += error*error;

            integ = h_njets[5][j]->IntegralAndError(1,10,error);

            sum_int[5][i] += integ;
            sum_err[5][i] += error*error;

            integ = h_njets[6][j]->IntegralAndError(1,10,error);

            sum_int[6][i] += integ;
            sum_err[6][i] += error*error;
        
        
        }
    }

    ofstream tableNF;
    tableNF.open("tableNF.txt");

    tableNF << setprecision(3) << "\n";

    tableNF<< "2 SS leptons";
    for(int i = 0; i < 7; i++)
        tableNF << " & $" << sum_int[0][i] << "\\pm" << TMath::Sqrt(sum_err[0][i]) << "$ ";
    double backgr = 0.;
    for(int i = 0; i < 6; i++)
        backgr += sum_int[0][i];
    tableNF<< "& " << calculateZbi(sum_int[0][6], backgr, 0.3) << "\\\\ \\hline \n";

    tableNF << "3 loose jets";
    for(int i = 0; i < 7; i++)
        tableNF << " & $" << sum_int[1][i] << "\\pm" << TMath::Sqrt(sum_err[1][i]) << "$";
    backgr = 0.;
    for(int i = 0; i < 6; i++)
        backgr += sum_int[1][i];
    tableNF<< "& " << calculateZbi(sum_int[1][6], backgr, 0.3) << "\\\\ \\hline \n";

    tableNF << "1 medium b-jets";
    for(int i = 0; i < 7; i++)
        tableNF << "& $" << sum_int[2][i] << "\\pm" << TMath::Sqrt(sum_err[2][i]) << "$ ";
    backgr = 0.;
    for(int i = 0; i < 6; i++)
        backgr += sum_int[2][i];
    tableNF<< "& " << calculateZbi(sum_int[2][6], backgr, 0.3) << "\\\\ \\hline \n";


    tableNF << "$H_{\\rm T} >$ 150 GeV ";
    for(int i = 0; i < 7; i++)
        tableNF << "& $" << sum_int[3][i] << "\\pm" << TMath::Sqrt(sum_err[3][i]) << "$ ";
    backgr = 0.;
    for(int i = 0; i < 6; i++)
        backgr += sum_int[3][i];
    tableNF<< "& " << calculateZbi(sum_int[3][6], backgr, 0.3) << "\\\\ \\hline \n";


    tableNF << "$|m_{ll} - M_{Z}| > $ 15 GeV  ";
    for(int i = 0; i < 7; i++)
        tableNF << "& $" << sum_int[4][i] << "\\pm" << TMath::Sqrt(sum_err[4][i]) << "$ ";
    backgr = 0.;
    for(int i = 0; i < 6; i++)
        backgr += sum_int[4][i];
    tableNF<< "& " << calculateZbi(sum_int[4][6], backgr, 0.3) << "\\\\ \\hline \n";

    
    tableNF << " 2 medium b-jets  ";
    for(int i = 0; i < 7; i++)
        tableNF << "& $" << sum_int[5][i] << "\\pm" << TMath::Sqrt(sum_err[5][i]) << "$ ";
    backgr = 0.;
    for(int i = 0; i < 6; i++)
        backgr += sum_int[5][i];
    tableNF<< "& " << calculateZbi(sum_int[5][6], backgr, 0.3) << "\\\\ \\hline \n";

    tableNF << " MET $> 50$ GeV  ";
    for(int i = 0; i < 7; i++)
        tableNF << "& $" << sum_int[6][i] << "\\pm" << TMath::Sqrt(sum_err[6][i]) << "$ ";
    backgr = 0.;
    for(int i = 0; i < 6; i++)
        backgr += sum_int[6][i];
    tableNF<< "& " << calculateZbi(sum_int[6][6], backgr, 0.3) << "\\\\ \\hline \n";

    //____________________________________________________________________________________________________

    const int SRNumber = 6;

    double yields[2][SRNumber] = {0.};//fakes, rare, s1, s2 
    double yieldsErrs[2][SRNumber] = {0.};//fakes, rare, s1, s2
    int yInd[3] = {0,36,37};

    for (int i=0; i!=2; ++i) {
        for (int k=0; k!=SRNumber; ++k) {
            yields[i][k] = 0;
            yieldsErrs[i][k] = 0;
        }
    }

    for (int i=0; i!=2; ++i) {
        for (int j=yInd[i]; j!=yInd[i+1]; ++j) {
            for (int k=0; k!=SRNumber; ++k) {
                //if (h_SR[j]->GetBinContent(indSR[k]+1) >=0) {
                    yields[i][k]+=h_SR[j]->GetBinContent(k+1);
                    if (h_SR[j]->GetBinContent(k+1) == 0)
                        yieldsErrs[i][k]+=weights[j]*weights[j];
                    else
                        yieldsErrs[i][k]+=(h_SR[j]->GetBinError(k+1))*(h_SR[j]->GetBinError(k+1));
                //}
            }
        }
    }

    TH1F* h_SR_yield[2];
    for (int i =0; i!=2; ++i) {
        TString name = Form("h_SR_yield_%d",i);
        h_SR_yield[i] = new TH1F(name,name,400,0,400);
        h_SR_yield[i]->Sumw2();
        for (int k=0; k!=SRNumber; ++k) {
            h_SR_yield[i]->SetBinContent(k+1,yields[i][k]);
            h_SR_yield[i]->SetBinError(k+1,sqrt(yieldsErrs[i][k]));
        }
    }

 //_____________________________________________________________________________________________________    
    // DATACARD creation
    //gSystem->Exec("rm datacard.txt"); // delete previous tex file
    ofstream datacard;
    datacard.open("datacard.txt"); // create a new datacard file
    datacard << fixed << showpoint << setprecision(4);
    datacard << "Date: May, 26, 2015" << endl;
    datacard << "Description: SUSY search, first probe" << endl;
    datacard << "lumi 10 pb^-1" << endl;
    const int Nback = 1;
    datacard << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
    datacard << "imax " <<  SRNumber << " number of channels" << endl;
    datacard << "jmax " <<  Nback << " number of backgrounds" << endl;
    datacard << "kmax " <<  (Nback + 1) * SRNumber + 2 << " number of nuisance parameters" << endl;
    datacard << "shapes * * FAKE" << endl;
    //datacard << "kmax " << 3 << " number of nuisance parameters" << endl;
    datacard << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
    int iBin = 1;
    datacard << "Observation ";
    for (int i = 0; i < SRNumber; i++) {
        //datacard << 0 << " ";
        if(h_SR_yield[0]->GetBinContent(iBin) <= 0)
            datacard << 0.0001 << " ";
        else
            datacard << h_SR_yield[0]->GetBinContent(iBin) << " ";
        iBin++;
    }

    datacard << endl;
    datacard << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
    datacard << "bin ";
    for (int i = 0; i < SRNumber ; i++)
        datacard << i + 1 << " " << i + 1 << " " ;
    datacard << endl;
    
    datacard << "process " ;
    for (int i = 0; i < SRNumber ; i++)
        datacard << "S B1 ";
    datacard << endl;
    
    datacard << "process ";
    for (int i = 0; i < SRNumber ; i++)
        datacard << "0 1 ";
    datacard << endl;
    
    datacard << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
    datacard << "rate ";
    iBin = 1;
    for (int i = 0; i < SRNumber; i++) {
        
        if ( h_SR_yield[1]->GetBinContent(iBin) <= 0.)
            datacard << 0.0001 << " ";
        else
            datacard << h_SR_yield[1]->GetBinContent(iBin) << " ";
        
        
        if (h_SR_yield[0]->GetBinContent(iBin) <= 0.)
            datacard << 0.0001 << " ";
        else
            datacard << h_SR_yield[0]->GetBinContent(iBin) << " ";
        
        iBin++;
    }
    datacard << endl;
    datacard << "------------------------------------------------------------------------------------------------------------------------------------" << endl;
    
    iBin = 1;
    
    for (int i = 0; i < SRNumber; i++) {
                    
        datacard << (Nback + 1) * i+1 << " lnN ";
                    
        for (int num = 0; num < (Nback + 1) * i; num++)
        datacard << "1.00 ";
                    
        if (h_SR_yield[1]->GetBinContent(iBin) <= 0)
            datacard << "1.0";
        else
            datacard << 1 +  TMath::Abs(h_SR_yield[1]->GetBinError(iBin)/h_SR_yield[1]->GetBinContent(iBin));
                    
        for (int num = (Nback + 1) * i + 1; num < (Nback + 1) * SRNumber; num++)
            datacard << " 1.00";
        datacard << endl;
                    
        //_________________________________________________________________________
                    
        datacard << (Nback + 1) * i + 2 << " lnN ";
                    
        for (int num = 0; num < (Nback + 1) * i + 1; num++)
             datacard << "1.00 ";
                    
        if ( h_SR_yield[0]->GetBinContent(iBin) <= 0)
             datacard << "1.0";
        else
             datacard << 1 +  TMath::Abs(1/sqrt(10 * h_SR_yield[0]->GetBinContent(iBin))) ;
            //datacard << 1 +  TMath::Abs(h_SR_yield[0]->GetBinError(iBin)/h_SR_yield[0]->GetBinContent(iBin));
                    
        for (int num = (Nback + 1) * i + 2; num < (Nback + 1) * SRNumber; num++)
            datacard << " 1.00";
        datacard << endl;
                    
        iBin++;
    }
    
    //_________________________________________________________________________
    
    datacard << (Nback + 1) * SRNumber + 1 << " lnN "; // without split in MT2
    
    for (int num = 0; num < SRNumber ; num++)
        datacard << "1.00 1.50 ";
    
    datacard << endl;
    
    datacard << (Nback + 1) * SRNumber + 2 << " lnN "; // without split in MT2
    
    for (int num = 0; num < SRNumber; num++)
        datacard << "1.20 1.00 ";
    
    datacard << endl;
    
    datacard.close();
    
    std::cout << "datacard DONE" << std::endl;

    
    //h_SR_yield[0]  - TH1F, which has 30 bins, - non prompt backgrund
    //h_SR_yield[1] - prompt background
    //h_SR_yield[signal] - SUSY signal

    
    //return 0;

}


void showHist(THStack *stack, string title, string titleX, string titleY, double num){ 

    if(title == "norm"){
        /*
        TList *histos = stack->GetHists();
        TH1F *sum = new TH1F("sum","sum of histograms",100,-4,4);
        TIter next(histos);
        TH1F *hist;
        while ((hist =(TH1F*)next())) {
            cout << "Adding " << hist->GetName() << endl;
            sum->Add(hist);
        }
        */
        ((TH1D*)stack->GetStack()->Last())->DrawNormalized("hist");
        stack->SetMaximum(stack->GetMaximum() * num);
        //histo_stack->Draw("hist");
        //return;
    }
    else{
        //stack->GetXaxis()->SetTitle(titleX.c_str());
        stack->SetMaximum(stack->GetMaximum() * num);
        stack->Draw("hist");
    }
}

int main(int argc, char *argv[]){

    TApplication *rootapp = new TApplication("example", &argc, argv);

    readTree();

    rootapp->Run();

    return 0;
}

double calculateZbi(double signal, double bkg, double unc){

 double n_on = signal+bkg;
 double mu_b_hat=bkg;
 double sigma_b=unc*bkg;
 double tau = mu_b_hat/(sigma_b*sigma_b);
 double n_off = tau*mu_b_hat;
 double P_Bi = TMath::BetaIncomplete(1./(1.+tau),n_on,n_off+1);
 double Z_Bi = sqrt(2)*TMath::ErfInverse(1 - 2*P_Bi);
 return Z_Bi;
 //  std::cout<<"The calculated Zbi for a signal of "<<signal<<" events and background of "<<bkg<<" events with a systematic uncertainty of "<<unc*100<<"% is "<<Z_Bi<<std::endl;

}


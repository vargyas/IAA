#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TLegend.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLine.h>
#include <TMatrixTSym.h>
#include <TFitResult.h>

#include "mplot.h"
#include "mcorr.h"
#include "loadFilip.h"
#include "loadFilipIAAfinal.h"
#include "miaa.h"
#include "JFiete.h"



// -------------------------------------
// Filip's prelimiary IAA(pta)
// obtained with plotdigitizer...
// from fig. IAA.png in talks/introfigs
// -------------------------------------
TGraphAsymmErrors * filipIAApta(bool stat)
{
    TGraphAsymmErrors * fiaa = nullptr;

    double xval[] = { 0.85, 1.5, 2.5, 3.5, 5., 7 };
    double yval[] = {1.53677, 1.41644, 1.32706, 1.12765, 1.2136, 1.1861};
    double xerrminus[]     = { 0.15, 0.5, 0.5, 0.5, 1.0, 1.0 }; // symm. bins, same as plus
    double ystaterrminus[] = { 0.28535, 0.10314, 0.06188, 0.03782, 0.03782, 0.04469 }; // symm. bins, same as plus
    double ysysterrplus[]  = { 0.56383, 0.13065, 0.11001, 0.08939, 0.09627, 0.09282 };
    double ysysterrminus[] = { 0.15127, 0.07563, 0.06533, 0.05844, 0.05157, 0.04126 };

    if(stat) // pp is not background removed
    {
        fiaa = new TGraphAsymmErrors(6, xval, yval, xerrminus, xerrminus, ystaterrminus, ystaterrminus);
        fiaa->SetName("gFilipPrelimIAAstat");
    }
    if(!stat)
    {
        fiaa = new TGraphAsymmErrors(6, xval, yval, xerrminus, xerrminus, ysysterrminus, ysysterrplus );
        fiaa->SetName("gFilipPrelimIAAsyst");
        fiaa->SetFillColorAlpha(41,0.2) ;
    }
    fiaa->SetTitle("Filip's preliminary I_{AA}(pta)");

    return fiaa;
}
// -------------------------------------
// published near-side I_AA, pTt 8-15GeV, C:0-5%
// http://arxiv.org/pdf/1110.0121.pdf
// -------------------------------------
TGraphAsymmErrors * publishedIAA(int icase)
{
    TGraphAsymmErrors * g = nullptr; // return graph

    double xval[] = {3.5, 5.0, 7.0, 9.0};
    double xerrminus[] = {0.5, 1.0, 1.0, 1.0};
    double xerrplus[] = {0.5, 1.0, 1.0, 1.0};

    double yval1[] = {1.4, 1.21, 1.12, 1.25};
    double yerrminus1[] = {0.1118033988749895, 0.09848857801796104, 0.08944271909999159, 0.12041594578792295};
    double yerrplus1[] = {0.1118033988749895, 0.09848857801796104, 0.08944271909999159, 0.12041594578792295};
    int numpoints = 4;
    TGraphAsymmErrors * p8088_d1x1y1 = new TGraphAsymmErrors(numpoints, xval, yval1, xerrminus, xerrplus, yerrminus1, yerrplus1);
    p8088_d1x1y1->SetName("/HepData/8088/d1x1y1");
    p8088_d1x1y1->SetTitle("PRL flat backg.");

    double yval2[] = {1.29, 1.19, 1.12, 1.25};
    double yerrminus2[] ={0.10295630140987, 0.08544003745317531, 0.08944271909999159, 0.12041594578792295};
    double yerrplus2[] = {0.10295630140987, 0.08544003745317531, 0.08944271909999159, 0.12041594578792295};
    TGraphAsymmErrors * p8088_d1x1y2 = new TGraphAsymmErrors(numpoints, xval, yval2, xerrminus, xerrplus, yerrminus2, yerrplus2);
    p8088_d1x1y2->SetName("/HepData/8088/d1x1y2");
    p8088_d1x1y2->SetTitle("PRL v2 backg.");

    double yval3[] = {1.25, 1.22, 1.09, 1.3};
    double yerrminus3[] = {0.1140175425099138, 0.09848857801796104, 0.09433981132056604, 0.12727922061357855};
    double yerrplus3[] = {0.1140175425099138, 0.09848857801796104, 0.09433981132056604, 0.12727922061357855};
    TGraphAsymmErrors * p8088_d1x1y3 = new TGraphAsymmErrors(numpoints, xval, yval3, xerrminus, xerrplus, yerrminus3, yerrplus3);
    p8088_d1x1y3->SetName("/HepData/8088/d1x1y3");
    p8088_d1x1y3->SetTitle("PRL #eta-gap backg.");

    if(icase == 0)
        g = new TGraphAsymmErrors( *p8088_d1x1y1 ); // 0: flat back.
    if(icase == 1)
        g = new TGraphAsymmErrors( *p8088_d1x1y2 ); // 1: v2 backg.
    if(icase == 2)
        g = new TGraphAsymmErrors( *p8088_d1x1y3 ); // 2: etaGap backg.

    return g;
}
// -------------------------------------
// Published $I_{\rm AA}$ from pi0-hadron
// correlation: arxiv:1608.07201
// -------------------------------------
TGraphAsymmErrors * publishedIAA_pi0()
{
    // TODO handle syst errors
    //TFile * pubFile = TFile::Open("../data/published/HEPData-ins1483164-1-root.root");
    //TH1F * h = (TH1F*) pubFile->Get("Table 3/Hist1D_y1");
    //return h;
    double p9206_d3x1y1_xval[] = { 0.75, 1.5, 2.5, 3.5, 5.0, 7.0, 9.0 };
    double p9206_d3x1y1_xerrminus[] = { 0.25, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0 };
    double p9206_d3x1y1_xerrplus[] = { 0.25, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0 };
    double p9206_d3x1y1_yval[] = { 1.805, 1.716, 1.668, 1.382, 1.181, 1.044, 1.196 };
    double p9206_d3x1y1_yerrminus[] = { 0.2830282671395209, 0.21306337085477645, 0.24039342753078752, 0.20972362766269326, 0.14494136745594752, 0.14458561477546789, 0.16617159805454118 };
    double p9206_d3x1y1_yerrplus[] = { 0.2830282671395209, 0.21306337085477645, 0.24039342753078752, 0.20972362766269326, 0.14494136745594752, 0.14458561477546789, 0.16617159805454118 };
    double p9206_d3x1y1_ystatminus[] = { 0.173, 0.114, 0.117, 0.12, 0.068, 0.076, 0.117 };
    double p9206_d3x1y1_ystatplus[] = { 0.173, 0.114, 0.117, 0.12, 0.068, 0.076, 0.117 };
    int p9206_d3x1y1_numpoints = 7;
    TGraphAsymmErrors * p9206_d3x1y1 = new TGraphAsymmErrors(p9206_d3x1y1_numpoints, p9206_d3x1y1_xval, p9206_d3x1y1_yval, p9206_d3x1y1_xerrminus, p9206_d3x1y1_xerrplus, p9206_d3x1y1_yerrminus, p9206_d3x1y1_yerrplus);
    p9206_d3x1y1->SetName("/HepData/9206/d3x1y1");
    p9206_d3x1y1->SetTitle("/HepData/9206/d3x1y1");
    //p9206_d3x1y1.Draw("AP");
    return p9206_d3x1y1;
}

double GetRelErr(TH1D * h, int ib)
{
    if( h->GetBinContent(ib)!=0 ) return h->GetBinError(ib)/h->GetBinContent(ib);
    else return 0;
}

void ChangeIAAErrors(int ic, int iptt, int ipta, MCorr * mpp, MCorr * maa, TH1D * hIAA)
{
    MTools mt;
    TH1D * _hPP_raw = (TH1D*)mt.Flip( mpp->hDEtaRaw2D[0][iptt][ipta] );
    TH1D * _hPP_mix = (TH1D*)mt.Flip( mpp->hDEtaMix2D[0][iptt][ipta] );
    TH1D * _hAA_raw = (TH1D*)mt.Flip( maa->hDEtaRaw2D[ic][iptt][ipta] );
    TH1D * _hAA_mix = (TH1D*)mt.Flip( maa->hDEtaMix2D[ic][iptt][ipta] );

    _hPP_raw->Rebin(2); _hPP_raw->Scale(1./2.);
    _hPP_mix->Rebin(2); _hPP_mix->Scale(1./2.);
    _hAA_raw->Rebin(2); _hAA_raw->Scale(1./2.);
    _hAA_mix->Rebin(2); _hAA_mix->Scale(1./2.);

    double errPP_raw=0;
    double errPP_mix=0;
    double errAA_raw=0;
    double errAA_mix=0;
    double calc=0;

    int nbinsx = hIAA->GetNbinsX();
    for(int ib=1; ib<=nbinsx; ib++)
    {
        // calculate the error by hand
        errPP_raw = GetRelErr(_hPP_raw,ib);
        errPP_mix = GetRelErr(_hPP_mix,ib);

        errAA_raw = GetRelErr(_hAA_raw,ib);
        errAA_mix = GetRelErr(_hAA_mix,ib);

        calc = TMath::Sqrt(errPP_raw*errPP_raw +errPP_mix*errPP_mix+ errAA_raw*errAA_raw+errAA_mix*errAA_mix);
        //calc = errPP_raw+errPP_mix+errAA_raw+errAA_mix;

        hIAA->SetBinContent(ib, hIAA->GetBinContent(ib));
        hIAA->SetBinError(ib, calc);
    }
    delete _hPP_raw;    delete _hPP_mix;    delete _hAA_raw;    delete _hAA_mix;
}

// -------------------------------------
// merge centrality dependent
// histo of AA with pp
// -------------------------------------
TH1D * MergeCentAAPP(TH1D * hAA, TH1D * hPP)
{
    // copy original binning
    const int nbins = hAA->GetNbinsX();
    double * xbins = new double[nbins+2];
    for(int ib=1;ib<=nbins+1;ib++) xbins[ib-1] = hAA->GetBinLowEdge(ib);

    // increase x proportionally to original bins
    // to make room for PP result
    xbins[nbins+1] = xbins[nbins]+(xbins[nbins]-xbins[nbins-1]);

    TH1D * hmer = new TH1D(Form("%s_merged",hAA->GetName()),"",nbins+1,xbins);

    // copy AA and PP and set the label text
    for(int ib=1; ib<=nbins;ib++)
    {
        hmer->SetBinContent(ib, hAA->GetBinContent(ib) );
        hmer->SetBinError(ib, hAA->GetBinError(ib) );
        hmer->GetXaxis()->SetBinLabel(ib, Form("%.0f-%.0f",xbins[ib-1],xbins[ib]));
    }
    hmer->SetBinContent(nbins+1, hPP->GetBinContent(1) );
    hmer->SetBinError(nbins+1, hPP->GetBinError(1) );
    hmer->GetXaxis()->SetBinLabel(nbins+1,"p-p");
    return hmer;
}

void DrawMiniYieldPlot(std::vector<double> yields, std::vector<TString> legList)
{
    Int_t colors[] = { 601, 634, 417, 433, 635, 619, 603, 636, 719, 435};

    static int imyp=0;
    const int N = yields.size();
    double xmax = *max_element(std::begin(yields), std::end(yields))*1.2;

    TH1D * h = new TH1D(Form("h%d",imyp++),"",N,0,1);
    h->GetXaxis()->SetRangeUser(0, xmax);
    for(Int_t ib=1; ib<=N; ib++)
        h->GetYaxis()->SetBinLabel(ib, legList.at(ib-1));
    h->GetXaxis()->SetTitle("yield");
    h->Draw();

    TH1D * hn[N];


    for(Int_t ileg=0; ileg<N; ileg++)
    {
        double xrat = 1/double(N);
        hn[ileg] = new TH1D(Form("h%d",imyp++),"",1,ileg*xrat, (ileg+1)*xrat);
        hn[ileg]->SetBinContent(1, yields.at(ileg));
        //hn[ileg]->SetMarkerColor(colors[ileg]);
        //hn[ileg]->SetLineColor(colors[ileg]);
        //hn[ileg]->SetLineWidth(2);

        hn[ileg]->SetFillColor(colors[ileg]);
        hn[ileg]->DrawCopy("hbar same");
    }
}

// -------------------------------------
//
//        M A I N     M A C R O
//
// -------------------------------------

//int main()
void processDraw()
{
    // -------------------------------------
    // processing correlation from JCORRAN
    // ROOT file with mcorr.h class
    // -------------------------------------

    int imcorr=0;
    MCorr * mc[2][6];

    AliJBin * Cent[2], * Vtx[2], * Eta[2], * Phi[2], * PTt[2], * PTa[2];
    int NumCent[2], NumVtx[2], NumPtt[2], NumPta[2], NumEta[2], NumPhi[2];

    // -------------------------------------
    // SETUP
    // canonical values and constants
    const double minpt = 3; // GeV
    const double minptFilip = 4;

    const double etaCut = 1.0;
    const double phiCut = 0.2*TMath::Pi();//*TMath::Pi();
    TString SSetup[2][10];
    const TString dataDir = "../data/jcorran/";
    int kTC = 0; // 0=TPCOnly, 1=RAA, 2=GlobalSDD

    // initalizing  pp
    // all setups are in: JDiHadronCorr_legotrain_allTrigg-CF_pp-1235_20170302-1857-2760GeV_LHC11a_p4_AOD113_withSDD.root

    // ALL: H0_T0, H0_T1, Hm_T0, Hm_T1, Hp_T1, Hm_T1
    //TString inNamePP[] = {dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_pp-1270_20170330-1749-2760GeV_LHC11a_p4_AOD113_withSDD.root",  // TPCOnly
    //                      dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_pp-1312_20170424-1753-2760GeV_LHC11a_p4_AOD113_withSDD.root"}; // GlobalSDD
    // WITHOUT SDD
    TString inNamePP[] = {dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_pp-1269_20170330-1749-2760GeV_LHC11a_p4_AOD113_noSDD.root", // TPCOnly
                          dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_pp-1413_20170830-1640-2760GeV_LHC11a_p4_AOD113_noSDD.root", // RAA
                          dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_pp-1311_20170424-1753-2760GeV_LHC11a_p4_AOD113_noSDD.root"}; // GlobalSDD
    int isetPP = 0;

    mc[0][isetPP] = new MCorr( imcorr++, kPP, inNamePP[kTC], "LHC11a", "113", "0", "noSDD", "TPCOnly_H0_T0", phiCut, etaCut, -10, minpt, 15., 0, true, false );
    mc[0][isetPP]->Initialize();
    SSetup[0][isetPP++] = {"all"};

    /*
    mc[0][isetPP] = new MCorr( imcorr++, kPP, inNamePP[kTC], "LHC11a", "113", "0", "withSDD", "TPCOnly_Hp_T0", phiCut, etaCut, -10, minpt, 15., 1, true, false );
    mc[0][isetPP]->Initialize();
    SSetup[0][isetPP++] = {"unlike-sign"};

    mc[0][isetPP] = new MCorr( imcorr++, kPP, inNamePP[kTC], "LHC11a", "113", "0", "withSDD", "TPCOnly_Hm_T0", phiCut, etaCut, -10, minpt, 15., 1, true, false );
    mc[0][isetPP]->Initialize();
    SSetup[0][isetPP++] = {"like-sign"};
    */

    // initializing aa
    TString inNameAA00[] = { dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3528_20170404-0250_runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-4035_20170830-2051-runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3606_20170424-1958_runlist_3-LHC10h_AOD86_MgFpMgFm.root"}; // H0_T0
    TString inNameAA01[] = { dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3529_20170404-0250_runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-4034_20170830-2052-runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3607_20170424-1958_runlist_3-LHC10h_AOD86_MgFpMgFm.root"}; // H0_T1
    TString inNameAAm0[] = { dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3530_20170404-0251_runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-4033_20170830-2052-runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3608_20170424-1959_runlist_3-LHC10h_AOD86_MgFpMgFm.root"}; // Hm_T0
    TString inNameAAm1[] = { dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3531_20170404-0251_runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-4036_20170830-2053-runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3609_20170424-2013_runlist_3-LHC10h_AOD86_MgFpMgFm.root"}; // Hm_T1
    TString inNameAAp0[] = { dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3532_20170404-0252_runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-4037_20170830-2053-runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3610_20170424-2013_runlist_3-LHC10h_AOD86_MgFpMgFm.root"}; // Hp_T0
    TString inNameAAp1[] = { dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3533_20170404-0252_runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-4032_20170830-2052-runlist_3-LHC10h_AOD86_MgFpMgFm.root",
                             dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_PbPb-3611_20170424-2013_runlist_3-LHC10h_AOD86_MgFpMgFm.root"}; // Hp_T1

    int isetAA=0;

    mc[1][isetAA] = new MCorr( imcorr++, kPbPb, inNameAA00[kTC], "LHC10h", "086", "3", "", "TPCOnly_H0_T0", phiCut, etaCut, -8, minpt, 15., 0, true, true);
    mc[1][isetAA]->Initialize();
    SSetup[1][isetAA++] = {"all"};
    /*
    mc[1][isetAA] = new MCorr( imcorr++, kPbPb, inNameAAm0[kTC], "LHC10h", "086", "3", "", "TPCOnly_Hm_T0", phiCut, etaCut, -8, minpt, 15., 1, true, false );
    mc[1][isetAA]->Initialize();
    SSetup[1][isetAA++] = {"unlike-sign"};

    mc[1][isetAA] = new MCorr( imcorr++, kPbPb, inNameAAp0[kTC], "LHC10h", "086", "3", "", "TPCOnly_Hp_T0", phiCut, etaCut, -8, minpt, 15., 1, true, false );
    mc[1][isetAA]->Initialize();
    SSetup[1][isetAA++] = {"like-sign"};
    */
    // -------------------------------------

    const int NSetup[2] = {isetPP, isetAA};
    //const int NSetup[2] = {1, 1};


    // fit, cent, setup, trigger, assoc
    int kF=3;// kC=5, kS=3, kT=5, kA=7;

    for(int it=0; it<2; it++)
    {
        Cent[it] = mc[it][0]->fhst->GetBin("Cent");   NumCent[it] = Cent[it]->Size();
        Vtx[it]  = mc[it][0]->fhst->GetBin("Vtx");    NumVtx[it]  = Vtx[it]->Size();
        PTt[it]  = mc[it][0]->fhst->GetBin("PTt");    NumPtt[it]  = PTt[it]->Size();
        PTa[it]  = mc[it][0]->fhst->GetBin("PTa");    NumPta[it]  = PTa[it]->Size();
        Eta[it]  = mc[it][0]->fhst->GetBin("EtaGap"); NumEta[it]  = Eta[it]->Size();
        Phi[it]  = mc[it][0]->fhst->GetBin("PhiGap"); NumPhi[it]  = Phi[it]->Size();
    }
    const int minPtt = PTt[0]->GetBin(minpt);
    const int minPttFilip = PTt[0]->GetBin(minptFilip);

    std::cout << "\nloading input files done...\n" << endl;



    // -------------------------------------
    // processing correlation
    // -------------------------------------
    const double fitrange_eta = 1.6;

    for( int itype=0; itype<2; itype++ ){
        for( int isetup=0; isetup<NSetup[itype]; isetup++ )
        {


            mc[itype][isetup]->LoadDEtaDPhi();
            mc[itype][isetup]->ProjectDEtaDPhi();
            mc[itype][isetup]->FitDEtaHistos("IEMRNSQ", fitrange_eta); //IEMRNSQ or WLRNSQ


            mc[itype][isetup]->DoQA();


            mc[itype][isetup]->DrawDEtaFlipAll(1);
            ////mc[itype][isetup]->DrawDEtaFlipAll(2);

            //mc[itype][isetup]->SaveOutput();

            //if(isetup==0) mc[itype][isetup]->DrawDEtaAll(1);

            //mc[itype][isetup]->DrawDEtaAll(2);

            //mc[itype][isetup]->Draw_Fit_INT_YieldPta();
            //
            //mc[itype][isetup]->DrawDEtaWingAll();
            //mc[itype][isetup]->DrawWingCorr();


            //mc[itype][isetup]->DrawMixed();
            //mc[itype][isetup]->DrawMixedVtx();

            //mc[itype][isetup]->DrawRaw2D1D();
            //mc[itype][isetup]->hDEtaRealFlip[0][1][1]->Print();
            //mc[itype][isetup]->DrawDEta(1,0,1,1); // 1-1D, 2-2D
            //mc[itype][isetup]->DrawDEta(1,0,2,1); // 1-1D, 2-2D
            //mc[itype][isetup]->DrawDEta(1,0,3,1); // 1-1D, 2-2D
            //mc[itype][isetup]->DrawDEta(2); // 1-1D, 2-2D
            //mc[itype][isetup]->DrawDEta();

            ////mc[itype][isetup]->DrawFitQA(1); //1-1D, 2-2D
            ////mc[itype][isetup]->DrawFitQA(2); //1-1D, 2-2D

            //mc[itype][isetup]->Draw2DHistos();
            //mc[itype][isetup]->DrawRawMixed();
            //delete mc[itype][isetup];

        }
    }




    // -------------------------------------
    //  PLOTTING PART
    // -------------------------------------
    int iplot = 1000;
    std::vector<TH1*>    hList;
    std::vector<TString> legList;

    // -------------------------------------
    // Draw inclusive pt comparison of
    // different setups (or cuts)
    // -------------------------------------

    MTools * mt = new MTools();

	double dEta = Eta[0]->At(NumEta[0]) - Eta[0]->At(0);

    /*
    for(int it=0; it<2; it++)
	{
		for(int ic=0; ic<NumCent[it]; ic++)
		{

			MPlot * mptt = new MPlot(iplot++, "p_{T}", "1/N_{eve} 1/(2#pip_{T}|#Delta#eta|) dN/dp_{T} [ (GeV/c)^{-2} ]", true);
			mptt->AddInfo( mc[it][0]->BuildInfo() );
			mt->DivideWithX( mc[it][0]->fhChargedPt[ic] ); // correcting with 1/p_T

			mc[it][0]->fhChargedPt[ic]->Scale(1./2./TMath::Pi()/dEta/mc[it][0]->GetNoEvents(ic), "width");

            hList  = { mc[it][0]->fhChargedPt[ic] };
            legList ={ SSetup[0][0] };
			mptt->addHList(hList, legList);

			mptt->SetLimitsXY(0.6, 30, 1E-10, 1E1);
			mptt->SetRatioLimitsXY(0.6, 100., 0.8, 2);
			mptt->SetLog(true, true);
			mptt->Draw();
			mptt->Save( Form("../figs/QA/incl_pt_%s_C%d","pp", ic));
		}
	}
    */



	// -------------------------------------
    // create IAA from the MCorr objects
    // -------------------------------------
	MIaa * miaa[3];

	for(int is=0; is<NSetup[0]; is++)
	{
		miaa[is] = new MIaa(mc[0][is], mc[1][is]);
		miaa[is]->CreateIAAdeta();
		miaa[is]->CreateIAApta();
	}

    // -------------------------------------
    // load Filip's IAA preliminary histos
    // with loadFilipIAAfinal.h
    // -------------------------------------
    FilipIAA * iaafilip = new FilipIAA();

    // -------------------------------------
    // load Filip's correlation
    // histos with loadFilip.h
    // -------------------------------------
    FilipHistos * hfilip[2];
    hfilip[0] = new FilipHistos("PP");
    hfilip[1] = new FilipHistos("AA");




    //void DrawFitBackg()
    /*
	{
        for(int iptt=minPtt; iptt<NumPtt[1]; iptt++){
            for(int ipta=0; ipta<NumPta[1]; ipta++){
                if( PTt[1]->At(iptt) <= PTa[1]->At(ipta) )
                    continue;
                hList.clear(); legList.clear();
                TH1D * h = nullptr;
                for(int ifit=0; ifit<2; ifit++) {
					h = (TH1D*) MergeCentAAPP(mc[1][0]->GetNotUniqueFitResultCent("hBackg",1,ifit,iptt,ipta), mc[0][0]->GetNotUniqueFitResultCent("hBackg",1,ifit,iptt,ipta));
                    hList.push_back( h );
                    legList.push_back( mc[0][0]->mfit_eta_1d[ifit][0][iptt][0]->GetName() );
                }
                MPlot * mbp_c = new MPlot(++iplot, "centrality [%]", "background (fit)", false);
                mbp_c->addHList(hList, legList, "PE");
                //myp_c->SetLimitsXY(0, 90, 0, 5);
                //mbp_c->AddInfo( mc[0][0]->BuildInfo() );
                mbp_c->AddInfo( mc[0][0]->BuildPTtTitle(iptt));
                mbp_c->AddInfo( mc[0][0]->BuildPTaTitle(ipta));
                mbp_c->CopyBinLabelsX(h);
                mbp_c->Draw();
                mbp_c->Save(Form("../figs/Fit/Backg_centALL_T0%dA0%d", iptt, ipta));
            }
        }
    }

    //void DrawFitWidthCent()
    {
        for(int iptt=minPtt; iptt<NumPtt[1]; iptt++){
			for(int ipta=0; ipta<NumPta[1]; ipta++){

				//FilipIAA * iaafilip = new FilipIAA();

                if( PTt[1]->At(iptt) <= PTa[1]->At(ipta) )
                    continue;

                hList.clear(); legList.clear();
				TH1D * h1 = nullptr;
				TH1D * h2 = nullptr;
                for(int ifit=0; ifit<2; ifit++) {
                    // only for Gaus, GenGaus, DoubleGaus
                    //if(ifit==0 || ifit==1 || ifit==4)
                    //{
						h1 = (TH1D*) MergeCentAAPP(mc[1][0]->GetNotUniqueFitResultCent("hWidth",1,ifit,iptt,ipta), mc[0][0]->GetNotUniqueFitResultCent("hWidth",1,ifit,iptt,ipta));
                        //h2 = (TH1D*) MergeCentAAPP(mc[1][0]->GetNotUniqueFitResultCent("hWidth",2,ifit,iptt,ipta), mc[0][0]->GetNotUniqueFitResultCent("hWidth",2,ifit,iptt,ipta));
						hList.push_back(h1);
                        //hList.push_back(h2);

                        legList.push_back( mc[0][0]->mfit_eta_1d[ifit][0][iptt][0]->GetName() );
                        //legList.push_back( mc[0][0]->mfit_eta_1d[ifit][0][iptt][0]->GetName()+"2D" );
					//}
                }
                MPlot * mwp_c = new MPlot(++iplot, "centrality [%]", "#sigma (fit)", false);
                mwp_c->addHList(hList, legList, "PE");
                //mwp_c->SetLimitsY(0, 1);
                mwp_c->AddInfo( mc[0][0]->BuildPTtTitle(iptt));
                mwp_c->AddInfo( mc[0][0]->BuildPTaTitle(ipta));
                mwp_c->CopyBinLabelsX(h1);
                mwp_c->Draw();
                mwp_c->Save(Form("../figs/Fit/Width_centALL_T0%dA0%d",iptt, ipta));
            }
        }
    }
    for(int iptt=minPtt; iptt<NumPtt[1]; iptt++){
        for(int ipta=0; ipta<NumPta[1]; ipta++){

            if( PTt[1]->At(iptt) <= PTa[1]->At(ipta) )
                continue;

            hList.clear(); legList.clear();
            TH1D * h1 = nullptr;
            int ifit=1;

            h1 = (TH1D*) MergeCentAAPP(mc[1][0]->GetNotUniqueFitResultCent("hExpo",1,ifit,iptt,ipta), mc[0][0]->GetNotUniqueFitResultCent("hExpo",1,ifit,iptt,ipta));
            hList.push_back(h1);
            legList.push_back( mc[0][0]->mfit_eta_1d[ifit][0][iptt][0]->GetName() );

            MPlot * mwp_c = new MPlot(++iplot, "centrality [%]", "exponent (fit)", false);
            mwp_c->addHList(hList, legList, "PE");
            mwp_c->SetLimitsY(0.5, 3);
            mwp_c->AddInfo( mc[0][0]->BuildPTtTitle(iptt));
            mwp_c->AddInfo( mc[0][0]->BuildPTaTitle(ipta));
            mwp_c->CopyBinLabelsX(h1);
            mwp_c->Draw();
            mwp_c->Save(Form("../figs/Fit/Expo_centALL_T0%dA0%d",iptt, ipta));
        }
    }
    */
    // -------------------------------------
	// PLOTTING Filip Delta Eta ALL
    // -------------------------------------
    const int ifit = 1;
    for(int itype=0; itype<2; itype++)
    {
        for(int ic=0; ic<NumCent[itype]; ic++){ //
			for(int iptt=minPttFilip; iptt<NumPtt[itype]; iptt++){
				for(int ipta=2; ipta<NumPta[itype]; ipta++){
			//for(int iptt=minPtt; iptt<NumPtt[1]; iptt++){
			//	for(int ipta=0; ipta<NumPta[1]; ipta++){
					if( PTt[itype]->At(iptt) <= PTa[itype]->At(ipta) )
                        continue;



					hList.clear(); legList.clear();

                    //if(itype==1)mc[itype][0]->hDEtaSig1D[0][ic][iptt][ipta]->Scale(0.5);
                    hList = { (TH1D*)mc[itype][0]->hDEtaSig1D[ifit][ic][iptt][ipta], hfilip[itype]->hDEta[ic][iptt-1][ipta-2] };

                    //hList ={ (TH1D*)mc[itype][0]->hDEtaRealFlip[ic][iptt][ipta], (TH1D*)mc[itype][0]->hDEtaReal2DFlip[ic][iptt][ipta], hfilip[itype]->hDEta[ic][iptt-1][ipta-2] };
                    legList ={ "AN(2016)", "AN(2012)" };
                    MPlot * mcorrfilip = new MPlot(++iplot, "|#Delta#eta|", "1/N_{trigg.}dN/d#Delta#eta",true);

                    mcorrfilip->addHList(hList, legList);
                    if(itype>0) mcorrfilip->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
                    mcorrfilip->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[itype]->At(iptt), PTt[itype]->At(iptt+1) ) );
                    mcorrfilip->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTa[itype]->At(ipta), PTa[itype]->At(ipta+1) ) );
                    mcorrfilip->SetRatioLimits(0.8, 1.2);
					mcorrfilip->SetGridXY();
                    mcorrfilip->Draw();

					//std::vector<double> yields = { mc[itype][0]->mfit_eta_1d[0][ic][iptt][ipta]->GetYield(), mc[itype][0]->mfit_eta_2d[0][ic][iptt][ipta]->GetYield() };
					//mcorrfilip->fCanvas->cd();
					//TPad *subpad = new TPad("subpad","",0.55,0.75,0.9,0.95);
					//subpad->Draw();
					//subpad->cd();
					//DrawMiniYieldPlot(yields, legList);

                    //mcorrfilip->SetRatioLimits(0., 2.);
					mcorrfilip->Save(Form("../figs/Corr/DEta_T%d_C0%dT0%dA0%d", itype, ic,iptt,ipta));
                }
            }
        }
    }



	std::cout << "starting  PLOTTING Filip IAA(dEta)\n";


    // -------------------------------------
    // PLOTTING Filip IAA(dEta)
    // -------------------------------------
    int drawfit = 1;
    for(int ic=0; ic<NumCent[1]; ic++){
		for(int iptt=minPttFilip; iptt<NumPtt[1]; iptt++){
			for(int ipta=2; ipta<NumPta[1]; ipta++) {
                if( PTt[1]->At(iptt) <= PTa[1]->At(ipta) )
                    continue;

				MPlot * miaafilip = new MPlot(++iplot, "|#Delta#eta|", "I_{AA}", false);

				int ifit=drawfit;
                //miaa[0]->hIAA_deta_1d[ifit][ic][iptt][ipta]->Scale(0.5);
                hList = { miaa[0]->hIAA_deta_1d[ifit][ic][iptt][ipta] };
                legList = { "AN(2016)"};
                miaafilip->addHList(hList, legList, "pe", "p");
                miaafilip->AddInfo( mc[1][0]->BuildInfo() );
                miaafilip->AddInfo( mc[1][0]->BuildCentTitle(ic) );
                miaafilip->AddInfo( mc[1][0]->BuildPTtTitle(iptt) );
                miaafilip->AddInfo( mc[1][0]->BuildPTaTitle(ipta)  );
				miaafilip->SetLimitsXY(0, 0.28, 0,2.5);
				//miaafilip->SetRatioLimitsXY(0, 0.28,0.9,1.1);
				//miaafilip->SetLimitsXY(0, 0.28, 0.5, 2.5);
				//miaafilip->SetRatioLimitsXY(0,0.28,0.9,1.15);
				miaafilip->SetGridXY();
                miaafilip->Draw();
				miaafilip->DrawThisTGraphAsymmErrors( iaafilip->gIAAsyst[ic][iptt-1][ipta-2], "E2 same", 24, 42 );
				miaafilip->DrawThisTGraphAsymmErrors( iaafilip->gIAAstat[ic][iptt-1][ipta-2], "PE same", 24, 1 );
				miaafilip->AddThisEntry(iaafilip->gIAAstat[ic][iptt-1][ipta-2], "AN(2012)");
				miaafilip->Save(Form("../figs/IAA/IAA_DEta_C0%dT0%dA0%d", ic,iptt,ipta));
            }
        }
    }
    cout << "\n Filip IAA(dEta) done...\n";

    // -------------------------------------
    // PLOTTING published IAA(pta)
    // -------------------------------------
    TGraphAsymmErrors * gPubIAA_c[3];  // published IAA (all 3 backg. subtraction)
	for(int it=0; it<3; it++) { gPubIAA_c[it] = publishedIAA(it); }
	TGraphAsymmErrors * gPubIAA_pi0 = publishedIAA_pi0();

    for(int ifit=0; ifit<2; ifit++)
    {
        for(int ic=0; ic<NumCent[1]; ic++)
        {
            int iptt = mc[0][0]->fPTt->GetBin(8.);
            // draw IAA from integral --------------------
            MPlot * mpiaa = new MPlot(iplot++, "p_{T, assoc} [GeV]", "I_{AA}", false);
			mpiaa->SetPalette(0); // set to open symbols

            hList = { miaa[0]->hIAA_eta_1d[ifit][0][iptt], miaa[0]->hIAA_eta_INT_1d[ifit][0][iptt] };
            legList = { "AN(2016) bc", "AN(2016) fit" };

			mpiaa->addHList(hList, legList, "PE");
			mpiaa->Draw();

			mpiaa->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
			mpiaa->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1)) );

			mpiaa->DrawThisTGraphAsymmErrors( gPubIAA_c[0], "PE same", 24, 1 );
			mpiaa->DrawThisTGraphAsymmErrors( gPubIAA_c[1], "PE same", 25, 1 );
			mpiaa->DrawThisTGraphAsymmErrors( gPubIAA_c[2], "PE same", 26, 1 );
			mpiaa->DrawThisTGraphAsymmErrors( gPubIAA_pi0, "PE same", 24, 417 );

			mpiaa->AddThisEntry(gPubIAA_c[0], "PRL data (flat)");
			mpiaa->AddThisEntry(gPubIAA_c[1], "PRL data (v_{2})");
			mpiaa->AddThisEntry(gPubIAA_c[2], "PRL data (#eta-gap)");
			mpiaa->AddThisEntry(gPubIAA_pi0, "#pi^{0} paper");

			//mpiaa->DrawThisTH1(miaa[0]->hIAA_eta_1d[ifit][ic][iptt], "PE same", 24, 601, 1);
			//mpiaa->DrawThisTH1(miaa[0]->hIAA_eta_INT_1d[ifit][ic][iptt], "PE same", 25, 634, 1);

			mpiaa->SetLimitsXY(0,10,0.5,3.5);
			mpiaa->Save(Form("../figs/IAA/IAA_C0%dF0%d",ic,ifit));
        }
    }
    cout << "\n IAA(pta) done...\n";




    int iptt_pub = PTt[0]->GetBin(8.); // always plot IAA and ICP for 8-15GeV bin
    for(int ifit=0; ifit<kF; ifit++)
    {
        MPlot * miaa_f_pta = new MPlot(iplot++, "p_{T, assoc} [GeV]", "I_{AA}", false);
		miaa_f_pta->SetPalette(0); // set to open symbols
        // draw IAA from with and compare with Filip's preliminary --------------------

		//hList   = { miaa[0]->hIAA_eta_1d[ifit][0][iptt_pub], miaa[0]->hIAA_eta_INT_1d[ifit][0][iptt_pub]  };
		//legList = { mc[0][0]->mfit_eta_1d[ifit][0][iptt_pub][2]->GetName(), "bin counting" };

        hList   = { miaa[0]->hIAA_eta_1d[ifit][0][iptt_pub], miaa[0]->hIAA_eta_INT_1d[ifit][0][iptt_pub] };
        legList = { "AN(2016) bc", "AN(2016) fit" };

        miaa_f_pta->addHList(hList, legList, "PE");
        miaa_f_pta->Draw();
		miaa_f_pta->DrawThisTGraphAsymmErrors( filipIAApta(true), "PE same", 26, 1 );
		miaa_f_pta->DrawThisTGraphAsymmErrors( filipIAApta(false), "5 same", 25, 41 );
		//miaa_f_pta->fPad->cd();
		//filipIAApta(0)->Draw("PE same");
		//filipIAApta(1)->Draw("2 same");
		miaa_f_pta->AddThisEntry(filipIAApta(true), "AN(2012) preliminary", "pe");
		//miaa_f_pta->AddThisEntry((TObject*)0, mc[0][0]->mfit_eta_2d[ifit][0][iptt_pub][0]->GetName(), "");
        miaa_f_pta->SetLimitsXY(0,8,0.5,3.5);
		//miaa_f_pta->DrawThisTH1(miaa[0]->hIAA_eta_1d[ifit][0][iptt_pub], "PE same", 24, 601, 1);
		//miaa_f_pta->DrawThisTH1(miaa[0]->hIAA_eta_INT_1d[ifit][0][iptt_pub], "PE same", 25, 634, 1);
        miaa_f_pta->Save(Form("../figs/IAA/IAA_FilipPrel_FIT%d", ifit));
    }
    cout << "\n Filip IAA(pta) done...\n";



















/*
        TH1D * hIAAFilip[5][5][5]; // histogram which

        for(int ic=0; ic<NumCent[1]; ic++){
            for(int iptt=minPtt; iptt<NumPtt[1]; iptt++){ 
                for(int ipta=2; ipta<NumPta[1]; ipta++) { // start from 2-3
                    if( PTt[1]->At(iptt) <= PTa[1]->At(ipta) )
                        continue;

                    //hIAAdEta[0][ic][iptt][ipta]->Rebin(2);hIAAdEta[0][ic][iptt][ipta]->Scale(1./2.);
                    MPlot * miaafilip = new MPlot(++iplot, "|#Delta#eta|", "I_{AA}", true);
                    //hList = { hfilip[1]->hIAA[ic][iptt-1][ipta-2], hIAAdEta[0][ic][iptt][ipta], hIAAdEta_ppfit[0][ic][iptt][ipta], hIAAdEta_ppFilip[0][ic][iptt][ipta] };
                    //legList = {"AN(2012)", "AN(2016)", "AN(2016): pp from fit", "PbPb:AN(2016), pp:AN(2012)"};



                    //hIAAdEta_newbinning[0][ic][iptt][ipta]->Scale(5.);
                    //hIAAdEta_ppFilip_newbinning[0][ic][iptt][ipta]->Scale(5.);
                    //hIAAdEta_ppfit_newbinning[0][ic][iptt][ipta]->Scale(5.);

                    // Self-consistency check, reproduce Filip's I_AA (dEta) from correlation histos and compare to his I_AA
                    cout << ic << "\t" << iptt << "\t" << ipta << endl;
                    hIAAFilip[ic][iptt][ipta] = (TH1D*) hfilip[1]->hDEta[ic][iptt-1][ipta-2]->Clone();
                    MFit * fitAA = new MFit(0,0,hIAAFilip[ic][iptt][ipta], -1.6, 1.6, false);
                    hIAAFilip[ic][iptt][ipta]->Fit( fitAA->ffit, "EMRNSQ");
                    mt->subtractConstTH1(hIAAFilip[ic][iptt][ipta], fitAA->ffit->GetParameter(0));

                    TH1D * htmpPP = (TH1D*) hfilip[0]->hDEta[0][iptt-1][ipta-2]->Clone();
                    MFit * fitPP = new MFit(0,0,htmpPP, -1.6, 1.6, false);
                    htmpPP->Fit( fitPP->ffit, "EMRNSQ");
                    mt->subtractConstTH1(htmpPP, fitPP->ffit->GetParameter(0));

                    hIAAFilip[ic][iptt][ipta]->Divide( htmpPP );
                    // -----------------------------------------------

                    hList = { hfilip[1]->hIAA[ic][iptt-1][ipta-2], hIAAFilip[ic][iptt][ipta], hIAAdEta[0][ic][iptt][ipta] };
                    legList = {"AN(2012)", "AN(2012) repr.", "AN(2016)"};

                    //hList = { hfilip[1]->hIAA[ic][iptt-1][ipta-2], hIAAdEta_newbinning[0][ic][iptt][ipta], hIAAdEta_ppFilip_newbinning[0][ic][iptt][ipta], hIAAdEta_ppfit_newbinning[0][ic][iptt][ipta] };
                    //legList = {"AN(2012)", "AN(2016)", "AA:(2016) PP:(2012)", "AN(2016) pp Fit"};
                    miaafilip->addHList(hList, legList, "pe", "p");
                    miaafilip->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
                    miaafilip->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1) ) );
                    miaafilip->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTa[1]->At(ipta), PTa[1]->At(ipta+1) ) );
                    miaafilip->SetLimitsXY(0, 0.28, 0.5, 2.5);
                    miaafilip->SetRatioLimitsXY(0, 0.28, 0.5, 1.5);



                    miaafilip->Draw();
                    miaafilip->DrawThisTGraphAsymmErrors( iaafilip->gIAAsyst[ic][iptt-1][ipta-2], "E2 same", 24, 42 );
                    miaafilip->DrawThisTGraphAsymmErrors( iaafilip->gIAAstat[ic][iptt-1][ipta-2], "PE same", 24, 1 );

                    miaafilip->Save(Form("../figs/IAA/IAA_DEta_Filip_C0%dT0%dA0%d", ic, iptt, ipta));

                    // now plot the same with bars and no ratio
                    MPlot * miaafilip2 = new MPlot(++iplot, "|#Delta#eta|", "I_{AA}",false);
                    //hList = {hIAAdEta[0][ic][iptt][ipta], hIAAdEta_ppfit[0][ic][iptt][ipta], hIAAdEta_ppFilip[0][ic][iptt][ipta] };
                    //legList = {"AN(2016)", "AN(2016): pp from fit", "PbPb:AN(2016), pp:AN(2012)"};
                    hList = {hIAAdEta[0][ic][iptt][ipta]};
                    legList = {"AN(2016)"};
                    miaafilip2->addHList(hList, legList, "pe", "p");
                    miaafilip2->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
                    miaafilip2->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1) ) );
                    miaafilip2->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTa[1]->At(ipta), PTa[1]->At(ipta+1) ) );
                    miaafilip2->SetLimitsXY(0, 0.28, 0.5, 2.5);
                    miaafilip2->Draw();
                    miaafilip2->DrawThisTGraphAsymmErrors( iaafilip->gIAAsyst[ic][iptt-1][ipta-2], "E2 same", 24, 42 );
                    miaafilip2->DrawThisTGraphAsymmErrors( iaafilip->gIAAstat[ic][iptt-1][ipta-2], "PE same", 24, 1 );
                    miaafilip2->DrawThisTH1(hIAAdEta_newbinning[0][ic][iptt][ipta],"PE same",24,2);
                    miaafilip2->AddThisEntry(iaafilip->gIAAstat[ic][iptt-1][ipta-2], "AN(2012)");
                    miaafilip2->Save(Form("../figs/IAA/IAA_DEta_Filip2_C0%dT0%dA0%d", ic,iptt,ipta));
                }
            }
        }





        return;










        // compare DEta with Filip
        std::cout << "Compare DEta with Filip \n";

        newnbins = hfilip[1]->hDEta[1][1][1]->GetNbinsX();
        for(int ibins=1; ibins<=newnbins; ibins++){ newxbins[ibins-1] = hfilip[1]->hDEta[1][1][1]->GetBinLowEdge(ibins); }

        for(int itype=0; itype<2; itype++){
            for(int iptt=minPtt; iptt<NumPtt[itype]; iptt++){
                for(int ipta=2; ipta<NumPta[itype]; ipta++)
                { // start from 2-3
                    if( PTt[itype]->At(iptt) <= PTa[itype]->At(ipta) )
                        continue;
                    for(int ic=0; ic<NumCent[itype]; ic++)
                    {
                        htmp = (TH1D*)mc[itype][0]->hDEtaRealFlip[ic][iptt][ipta]->Clone(Form("htmpC%dT%dA%d",ic,iptt,ipta));
                        MFit * mfit_ours = new MFit(0,0,htmp, 0,1.4, false);
                        htmp->Fit(mfit_ours->ffit, "IEMRNSQ");
                        mt->subtractConstTH1(htmp, mfit_ours->ffit->GetParameter(0));

                        TH1D * htmp_new = (TH1D*)mt->NewHistoWithUniqueBins(htmp, newnbins-1, newxbins);

                        TString name = hfilip[itype]->hDEta[ic][iptt-1][ipta-2]->GetName();
                        TH1D * htmp_filip = (TH1D*) hfilip[itype]->hDEta[ic][iptt-1][ipta-2]->Clone(name+"_resub");
                        MFit * mfit_filip = new MFit(0,0,htmp_filip, 0, 1.4, false);
                        htmp_filip->Fit(mfit_filip->ffit, "IEMRNSQ");
                        mt->subtractConstTH1( htmp_filip, mfit_filip->ffit->GetParameter(0) );

                        if(itype==0)
                            // normalize pp to each other
                            mt->shiftToThisTail(htmp_filip, htmp_new, 1.2, 1.6);

                        //htmp_filip->Scale(1., "width");
                        //double scaleCorr = mc[itype][0]->fFilipCorr_eta[ic][iptt][ipta] ;
                        //htmp->Scale( scaleCorr );
                        //cout << "scaling " << Types[itype] << "\t" << scaleCorr << endl;
                        //htmp_new->Scale( scaleCorr );
                        //htmp_new->Scale(1./2.);

                        MPlot * mdetafilip = new MPlot(++iplot, "|#Delta#eta|", "1/N_{trigg.}dN/d|#Delta#eta|", true);
                        //hList = { hfilip[itype]->hDEta[ic][iptt-1][ipta-2], htmp_new };
                        hList = { htmp_filip, htmp_new };
                        legList = {"AN(2012)",  "AN(2016)"};

                        mdetafilip->addHList(hList, legList, "PE", "p");
                        mdetafilip->SetLimitsX(0, 1.4);
                        if(itype>0) mdetafilip->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
                        mdetafilip->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1) ) );
                        mdetafilip->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTa[1]->At(ipta), PTa[1]->At(ipta+1) ) );
                        mdetafilip->SetRatioLimitsXY(0, 1.4, 0.0, 2.0);
                        mdetafilip->Draw();
                        mdetafilip->DrawThisTF1(mfit_filip->ffit, "same");
                        mdetafilip->Save(Form("../figs/Corr/DEta_Filip_%s_C0%dT0%dA0%d", Types[itype].Data(),ic,iptt,ipta));

                        delete mfit_filip;
                        delete mfit_ours;
                    }
                }
            }
        }
        for(int iptt=minPtt; iptt<NumPtt[1]; iptt++){
            for(int ipta=2; ipta<NumPta[1]; ipta++) { // start from 2-3
                if( PTt[1]->At(iptt) <= PTa[1]->At(ipta) )
                    continue;

                TH1D * htmp_pp = (TH1D*)mc[0][0]->hDEtaRealFlip[0][iptt][ipta]->Clone();
                mt->subtractConstTH1( htmp_pp, 2 * mc[0][0]->mfit_eta_ggc[0][iptt][ipta]->ffit->GetParameter(0) );

                for(int ic=0; ic<NumCent[1]; ic++)
                {
                    TH1D * htmp_AA = (TH1D*)mc[1][0]->hDEtaRealFlip[ic][iptt][ipta]->Clone();
                    mt->subtractConstTH1( htmp_AA, 2 * mc[1][0]->mfit_eta_ggc[ic][iptt][ipta]->ffit->GetParameter(0) );

                    MPlot * mdetaAApp = new MPlot(++iplot, "|#Delta#eta|", "1/N_{trigg.}dN/d|#Delta#eta|", true);
                    hList = {hfilip[1]->hDEta[ic][iptt-1][ipta-2], hfilip[0]->hDEta[0][iptt-1][ipta-2] };
                    legList={"AA","pp"};
                    mdetaAApp->addHList(hList, legList, "PE", "p");
                    mdetaAApp->SetLimitsX(0, 1.0);
                    mdetaAApp->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
                    mdetaAApp->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1) ) );
                    mdetaAApp->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTa[1]->At(ipta), PTa[1]->At(ipta+1) ) );
                    mdetaAApp->SetRatioLimitsXY(0, 1.0, 0.0, 2.0);
                    mdetaAApp->Draw();
                    mdetaAApp->Save(Form("../figs/Corr/DEta_Filip_AA-PP_C0%dT0%dA0%d", ic,iptt,ipta));

                    MPlot * mdetaAApp2 = new MPlot(++iplot, "|#Delta#eta|", "1/N_{trigg.}dN/d|#Delta#eta|", true);
                    hList = {htmp_AA, htmp_pp};
                    mdetaAApp2->addHList(hList, legList, "PE", "p");
                    mdetaAApp2->SetLimitsX(0, 1.0);
                    mdetaAApp2->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
                    mdetaAApp2->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1) ) );
                    mdetaAApp2->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTa[1]->At(ipta), PTa[1]->At(ipta+1) ) );
                    mdetaAApp2->SetRatioLimitsXY(0, 1.0, 0.0, 2.0);
                    mdetaAApp2->Draw();
                    mdetaAApp2->Save(Form("../figs/Corr/DEta_OURS_AA-PP_C0%dT0%dA0%d", ic,iptt,ipta));
                }
            }
        }

    }




return;




    if(bPlotFitResults)
    {

        // -------------------------------------
        // Plot width(pta) from fit
        // -------------------------------------
        TH1D * hWidthCent[2][kT][kA], * hExpoCent[kT][kA];
        
        double * centbins = new double[10];
        for(int i=0;i<=kC;i++) { centbins[i]=Cent[1]->At(i);  }
        centbins[kC+1]=100;

        for(int iptt=minPtt; iptt<NumPtt[1]; iptt++)
        {
            for(int ipta=0; ipta<NumPta[1]; ipta++)
            {
                if( PTt[1]->At(iptt) < PTa[1]->At(ipta) )
                    continue;

                for(int ifit=0; ifit<2; ifit++) {
                    hWidthCent[ifit][iptt][ipta] = new TH1D( Form("hWidthCent_F0%dT0%dA0%d",ifit,iptt,ipta),"",kC+1,centbins);
                }
                for(int ic=0; ic<NumCent[1]; ic++) {
                    hWidthCent[0][iptt][ipta]->SetBinContent(ic+1, mc[1][0]->mfit_eta_gc[ic][iptt][ipta]->GetWidth() );
                    hWidthCent[0][iptt][ipta]->SetBinError(ic+1, 0 );
                    hWidthCent[1][iptt][ipta]->SetBinContent(ic+1, mc[1][0]->mfit_eta_ggc[ic][iptt][ipta]->GetWidth() );
                    hWidthCent[1][iptt][ipta]->SetBinError(ic+1, 0 );
                }
                hWidthCent[0][iptt][ipta]->SetBinContent(NumCent[1]+1, mc[0][0]->mfit_eta_gc[0][iptt][ipta]->GetWidth()); // pp data
                hWidthCent[0][iptt][ipta]->SetBinError(NumCent[1]+1, 0 ); // pp data
                hWidthCent[1][iptt][ipta]->SetBinContent(NumCent[1]+1, mc[0][0]->mfit_eta_ggc[0][iptt][ipta]->GetWidth()); // pp data
                hWidthCent[1][iptt][ipta]->SetBinError(NumCent[1]+1, 0); // pp data

                hWidthCent[0][iptt][ipta]->Print();
                hWidthCent[1][iptt][ipta]->Print();
            }
        }
        for(int iptt=minPtt; iptt<NumPtt[1]; iptt++) {
            for(int ipta=0; ipta<NumPta[1]; ipta++) {
                if( PTt[1]->At(iptt) < PTa[1]->At(ipta) )
                    continue;

                MPlot * mwidth1 = new MPlot(++iplot, "|#Delta#eta|", "1/N_{trigg.}dN/d|#Delta#eta|", false);
                hList = {hWidthCent[0][iptt][ipta], hWidthCent[1][iptt][ipta]};
                legList = {"Gaus", "Gen.Gaus"};
                mwidth1->addHList(hList, legList, "PE", "p");
                mwidth1->SetLimitsY(0, 1.0);
                mwidth1->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1) ) );
                mwidth1->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTa[1]->At(ipta), PTa[1]->At(ipta+1) ) );
                mwidth1->Draw();
                mwidth1->Save(Form("../figs/Fit/Width_Cent_T0%dA0%d",iptt,ipta));
            }
        }

        cout << "end the fucking width_cent\n";

    }



    */






    /*
    cout << "Yield (dEta): ours vs. Filip" << endl;
    for(int itype=0; itype<2; itype++){
        for(int isetup=0; isetup<NSetup[itype]; isetup++){
            for(int ic=0; ic<NumCent[itype]; ic++){

                int ic_ours = 0;
                if(itype==1) ic_ours=ic+1;
                if(ic_ours==4) 
                    break;

                for(int iptt=0; iptt<NumPttNew; iptt++){
                // ---------
                    for(int ipta=0; ipta<NumPtaNew; ipta++){
                        if( PTtBo[iptt] <= PTaBo[ipta] )
                            continue;

                        MPlot * metaFilip = new MPlot(iplot++, "|#Delta#eta|", "1/N_{trigg.}dN/d|#Delta#eta|", true);
                        hList   = { hfilip[itype]->hDEta[ic][iptt][ipta], hDEtaRawFlip[itype][isetup][ic_ours][iptt][ipta], hDEtaRealFlip[itype][isetup][ic_ours][iptt][ipta] };
                        legList = {"Filip", "Our raw", "Our real"};
                        metaFilip->addHList(hList, legList, "l", "P", "l"); // fasz: why is this not working??? 
                        metaFilip->Draw();
                        if(itype==0) metaFilip->AddInfo( "pp" );
                        if(itype==1) metaFilip->AddInfo( Form("Cent: %.0f-%.0f %%",CentBo[ic], CentBo[ic+1])); 
                        metaFilip->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTtBo[iptt], PTtBo[iptt+1]) );
                        metaFilip->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTaBo[ipta], PTaBo[ipta+1]) );
                        metaFilip->SetRatioLimitsXY(0,1.6,-2,3);
                        metaFilip->Save(Form("../figs/PubYield/PubYied_Filip_%s_S0%d_C0%dT0%dA0%d", Types[itype].Data(), isetup, ic_ours, iptt, ipta));
                    }
                }
            }
        }
    }

    // draw ratio plot here

    for(int ic=1; ic<5; ic++){
        for(int ipta=0; ipta<4; ipta++){
            mep[ic][ipta] = new MPlot(iplot++, "|#Delta#eta|", "yield", true, 0.5);
            hList = { hDEtaFlip[0][0][1][ipta], hDEtaFlip[1][ic][1][ipta] };
            legList = {"pp", "PbPb" };
            mep[ic][ipta]->addHList(hList, legList, "PE", "p");
            mep[ic][ipta]->SetLimitsX(0, 0.3);
            mep[ic][ipta]->SetRatioLimitsXY(0,0.3, 0,2);
            mep[ic][ipta]->SetRatioLabel("I_{AA}");

            mep[ic][ipta]->Draw();
            mep[ic][ipta]->Save(Form("../figs/PubYield/IAA_DEta_C0%dA0%d", ic, ipta));
        }
    }


    cout << "I_AA (dEta): ours vs. Filip" << endl;
    // draw separate plot here
    // our+Filip:
    // ptt = 0: 4-6, 1:6-8, 2:8-15 Filip
    // pta: 0:2-3, 1:3-4, 2:4-6, 3:6-8
    for(int iptt=0; iptt<NumPttNew; iptt++){
        for(int ipta=0; ipta<NumPtaNew-1; ipta++){ // omit last because of Filip
            if( PTtBo[iptt] <= PTaBo[ipta] ) 
                continue;
            for(int ic=1; ic<5; ic++){
                for(int isetup=0; isetup<NSetup[1]; isetup++)
                {
                    hDEtaRealFlip[1][isetup][ic][iptt][ipta]->Divide( hDEtaRealFlip[0][isetup][0][iptt][ipta] );
                    hDEtaRawFlip[1][isetup][ic][iptt][ipta]->Divide( hDEtaRawFlip[0][isetup][0][iptt][ipta] );
                }
                MPlot * miaa = new MPlot(iplot++, "|#Delta#eta|", "I_{AA}", true);
                hList = { hfilip[1]->hIAA[ic-1][iptt][ipta], hDEtaRawFlip[1][1][ic][iptt][ipta], hDEtaRealFlip[1][1][ic][iptt][ipta] };
                legList = {"Filip", "RAW: #phi<0.5 |v_{z}|<8 cm", "REAL: #phi<0.5 |v_{z}|<8 cm"};

                miaa->addHList(hList, legList, "PE", "p");
                miaa->SetLimitsXY(0, 0.28, 0.0, 2.5);
                miaa->SetAppearance( hfilip[1]->hIAA[ic-1][iptt][ipta], 2, 2, 24 );
                miaa->Draw();
                
                miaa->SetRatioLimitsXY(0, 0.28, 0.5, 1.5);
                miaa->AddInfo( Form("Cent: %.0f-%.0f %%",CentBo[ic-1], CentBo[ic])); 
                miaa->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTtBo[iptt], PTtBo[iptt+1]) );
                miaa->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTaBo[ipta], PTaBo[ipta+1]) );
                miaa->Save(Form("../figs/PubYield/IAA_DEta_Filip_C0%dT0%dA0%d", ic, iptt, ipta));
            }
        }
    }

    cout << "I_AA (dEta): reproduced from Filip's histos" << endl;

    TH1D * IAA_Filip[kT][kA];
    for(int iptt=0; iptt<NumPttNew; iptt++){
        for(int ipta=0; ipta<NumPtaNew-1; ipta++){ // omit last because of Filip
            if( PTtBo[iptt] <= PTaBo[ipta] ) 
                continue;
            IAA_Filip[iptt][ipta] = (TH1D*) hfilip[1]->hDEta[0][iptt][ipta]->Clone();
            IAA_Filip[iptt][ipta]->Divide( (TH1D*)hfilip[0]->hDEta[0][iptt][ipta]);
            
            
            MPlot * miaa_filip = new MPlot( iplot++, "|#Delta#eta|", "I_{AA}", true);
            hList = { hfilip[1]->hIAA[0][iptt][ipta], IAA_Filip[iptt][ipta] };
            legList={ "Filip's original IAA", "IAA from Filip's histos" };
            miaa_filip->addHList(hList, legList, "PE", "p");
                miaa_filip->SetLimitsXY(0, 0.28, 0.0, 2.5);
            miaa_filip->Draw();
                miaa_filip->SetRatioLimitsXY(0, 0.28, 0.7, 1.3);
            miaa_filip->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTtBo[iptt], PTtBo[iptt+1]) );
            miaa_filip->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTaBo[ipta], PTaBo[ipta+1]) );
            miaa_filip->Save(Form("../figs/PubYield/IAA_DEta_FilipVSFilip_C0%dT0%dA0%d", 0, iptt, ipta));
        }
    }
    */



        // -------------------------------------
        // Compare original + 2D method's projection
        // -------------------------------------

    /*
        int etaMinBin, etaMaxBin, phiMinBin, phiMaxBin;

        TString histname;
        for(int itype=0; itype<2; itype++){
            for(int iptt=minPtt; iptt<NumPtt[itype]; iptt++){
                for(int ipta=0; ipta<NumPta[itype]; ipta++) { 
                    if( PTt[itype]->At(iptt) <= PTa[itype]->At(ipta) )
                        continue;
                    for(int ic=0; ic<NumCent[itype]; ic++){
                        etaMinBin = mc[itype][0]->hDEtaDPhiReal[ic][iptt][ipta]->GetXaxis()->FindBin(-etaCut);
                        etaMaxBin = mc[itype][0]->hDEtaDPhiReal[ic][iptt][ipta]->GetXaxis()->FindBin(etaCut);
                        phiMinBin = mc[itype][0]->hDEtaDPhiReal[ic][iptt][ipta]->GetYaxis()->FindBin(-phiCut);
                        phiMaxBin = mc[itype][0]->hDEtaDPhiReal[ic][iptt][ipta]->GetYaxis()->FindBin(phiCut);
                        
                        histname = mc[itype][0]->hDEtaDPhiReal[ic][iptt][ipta]->GetName();
                        
                        he = (TH1D*) mc[itype][0]->hDEtaDPhiReal[ic][iptt][ipta]->ProjectionX( Form("%s_eta",histname.Data()), phiMinBin, phiMaxBin );
                        hp = (TH1D*) mc[itype][0]->hDEtaDPhiReal[ic][iptt][ipta]->ProjectionY( Form("%s_phi",histname.Data()), etaMinBin, etaMaxBin );
                        //he->Print();
                        //hp->Print();
                        //he->Scale(1./(phiMaxBin-phiMinBin));
                        //hp->Scale(2.*TMath::Pi()/(etaMaxBin-etaMinBin), "width");

                        //he->Scale(0.02);
                        //hp->Scale(27);

                        //he->Scale( 0.3 );
                        //hp = (TH1D*) mt->shrinkHist( hp, TMath::Pi() );
                        //hp->Scale( 0.6, "width" );
                        
                        MPlot * meta_c = new MPlot(iplot++, "#Delta#eta", "1/N_{trigg.}dN/d#Delta#eta", true);
                        hList   = { he, mc[itype][0]->hDEtaReal[ic][iptt][ipta] };
                        legList = { "2D correction", "1D correction" }; 
                        meta_c->addHList(hList, legList);
                        meta_c->AddInfo( mc[itype][0]->BuildInfo() );
                        if(itype==1) meta_c->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
                        meta_c->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1) ) );
                        meta_c->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTa[1]->At(ipta), PTa[1]->At(ipta+1) ) );
                        meta_c->Draw();
                        meta_c->Save(Form("../figs/Corr/2D1D_eta_%s_C%dT%dA%d",Types[itype].Data(),ic,iptt,ipta));

                        MPlot * mphi_c = new MPlot(iplot++, "#Delta#phi", "1/N_{trigg.}dN/d#Delta#phi", true);
                        hList   = { hp, mc[itype][0]->hDPhiReal[0][ic][iptt][ipta] };
                        legList = { "2D correction", "1D correction" }; 

                        mphi_c->AddInfo(  mc[itype][0]->BuildInfo() );
                        if(itype==1) mphi_c->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
                        mphi_c->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1) ) );
                        mphi_c->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTa[1]->At(ipta), PTa[1]->At(ipta+1) ) );

                        mphi_c->addHList(hList, legList);
                        mphi_c->SetLimitsX(-0.5,0.5);

                        mphi_c->Draw();
                        mphi_c->Save(Form("../figs/Corr/2D1D_phi_%s_C%dT%dA%d",Types[itype].Data(),ic,iptt,ipta));
                    }
                }
            }
        }
    */
}


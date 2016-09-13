#include <TH1.h>
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

#include "../inc/mplot.h"
#include "../inc/mcorr.h"
#include "../inc/loadFilip.h"
#include "../inc/loadFilipIAAfinal.h"
#include "../inc/miaa.h"
#include "../inc/JFiete.h"

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

    if(stat)
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
// gets yield of histogram after subtracting 
// constant from fit
// -------------------------------------
//double getYield(TH1 * h1, TF1 * fit)
//{
//    MTools * mt = new MTools();
//    TH1D * h2 = (TH1D*)h1->Clone();
//    mt->subtractConstTH1(h2, fit->GetParameter(0) );
//    double yield = h2->Integral(h2->FindBin(1e-6), h2->FindBin(1.1));
//    delete h2;
//    delete mt;
//    return yield;
//}

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
    MCorr * mc[2][16]; // [type: pp/AA][setup: below]

    AliJBin * Cent[2], * Vtx[2], * Eta[2], * Phi[2], * PTt[2], * PTa[2];
    int NumCent[2], NumVtx[2], NumPtt[2], NumPta[2], NumEta[2], NumPhi[2];

    // -------------------------------------
    // SETUP
    const double minpt = 4; // GeV
    const TString Types[] = {"pp", "PbPb"};
    //enum kTypes = {kPP, kPBPB};
    const double etaCut = 1.0;
    const double phiCut = 0.2; // TMath::Pi()*0.2

    const bool bPlotFitResults  = true;
    bool bPlotFilip             = true;

    const TString dataDir = "/Users/vargyas/Work/ALICE/IAA/data/jcorran/"; 

    // initalizing  pp
//    const TString inNamePP = dataDir+"newJCORRAN/JDiHadronCorr_legotrain_allTrigg-CF_pp-903_20160428-2137-2760GeV_LHC11a_p4_AOD113_withSDD.root";
    //const TString inNamePP = dataDir+"modularJCORRAN/JDiHadronCorr_legotrain_allTrigg-CF_pp-1116_20160818-1832-2760GeV_LHC11a_p4_AOD113_withSDD.root";
    //const TString inNamePP = dataDir+"modularJCORRANvtx/JDiHadronCorr_legotrain_allTrigg-CF_pp-1139_20160905-1807-2760GeV_LHC11a_p4_AOD113_noSDD.root";
    const TString inNamePP = dataDir+"modularJCORRANvtx/JDiHadronCorr_legotrain_allTrigg-CF_pp-1146_1146_20160908-1821-2760GeV_LHC11a_p4_AOD113_withSDD.root"; // dEta < 1.6
    mc[0][0] = new MCorr( kPP, inNamePP, "LHC11a", "113", "0", "withSDD", "TPCOnly", phiCut, etaCut, -10, minpt, 15, 0, true );
    mc[0][0]->Initialize();

    // initalizing  PbPb 
//    const TString inNameAA = dataDir+"newJCORRAN/JDiHadronCorr_legotrain_allTrigg-CF_PbPb-2594_20160513-1811_LHC10h_AOD160_runlist_3_rebinned4Filip.root";
    //const TString inNameAA = dataDir+"newJCORRAN/JDiHadronCorr_legotrain_allTrigg-CF_PbPb-2594_20160513-1811_LHC10h_AOD160_runlist_2.root";
    //const TString inNameAA = dataDir+"modularJCORRAN/JDiHadronCorr_legotrain_allTrigg-CF_PbPb-2910_20160818-2237_runlist_3-LHC10h_AOD86_MgFpMgFm_rebinned4Filip.root";
    // const TString inNameAA = dataDir+"modularJCORRAN/JDiHadronCorr_legotrain_allTrigg-CF_PbPb-2871_20160804-1855_runlist_3-LHC10h_AOD86_MgFpMgFm.root";
    //const TString inNameAA = dataDir+"noGeomAccCorr/JDiHadronCorr_legotrain_allTrigg-CF_PbPb-2009_20151123-1916_runlist_3-LHC10h_AOD86_MgFpMgFm.root";

    const TString inNameAA = dataDir+"modularJCORRANvtx/JDiHadronCorr_legotrain_allTrigg-CF_PbPb-2948_20160908-1830_runlist_3-LHC10h_AOD86_MgFpMgFm.root";

    //const TString inNameAA = dataDir+"modularJCORRANvtx/JDiHadronCorr_legotrain_allTrigg-CF_PbPb-2950_20160908-1834-LHC11h_AOD115_fullTPCFlow.root";
    //const TString inNameAA = dataDir+"modularJCORRANvtx/JDiHadronCorr_legotrain_allTrigg-CF_PbPb-2956_20160912-1725_runlist_1-LHC11h_AOD115_FemtoMinusPlusFemtoPlus.root";
    mc[1][0] = new MCorr( kPbPb, inNameAA, "LHC10h", "86", "3", "", "TPCOnly", phiCut, etaCut, -8, minpt, 15., 0, true );
    mc[1][0]->Initialize();
    // -------------------------------------


    const int NSetup[] = {1, 1}; // {pp, AA}
    // fit, cent, setup, trigger, assoc
    int kF=5, kC=5, kS=3, kT=5, kA=7;

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

    std::cout << "\nloading input files done...\n" << endl;

    if( inNameAA.Contains("rebinned4Filip") ) kC=4;
    //else bPlotFilip = false;


    //const TString fitnames[] = {"Gauss fit", "Gen.Gauss fit", "Kaplan fit", "Cauchy fit", "QGauss fit"};
    //const TString subnames[] = {"(flat)", "(v_{2})", "(#eta-gap)"};

    // -------------------------------------
    // processing correlation
    // -------------------------------------
    const double fitrange_eta = 1.4;
    const double fitrange_phi = 0.4;

    for( int itype=0; itype<2; itype++ ){
        for( int isetup=0; isetup<1; isetup++ )
        {
            mc[itype][isetup]->DoQA();

            mc[itype][isetup]->LoadDEtaDPhi();
            mc[itype][isetup]->ProjectDEtaDPhi(phiCut);
            //mc[itype][isetup]->DrawWingCorr();


            mc[itype][isetup]->LoadDEtaHistos();
//          mc[itype][isetup]->LoadDPhiHistos();



    //      mc[itype][isetup]->DrawCorr2D1D(true,kDEta);
//          mc[itype][isetup]->DrawCorr2D1D(true,kDPhi);
//          mc[itype][isetup]->Draw2DHistos();

            mc[itype][isetup]->FitDEtaHistos("EMRNSQ", fitrange_eta); //IEMRNSQ or WLRNSQ
//          mc[itype][isetup]->FitDPhiHistos("WLRNSQ", fitrange_phi);

            mc[itype][isetup]->DrawDEta(true);
            mc[itype][isetup]->DrawRaw2D1D();
//          mc[itype][isetup]->DrawRawMix();
//          mc[itype][isetup]->DrawDPhi(true);
//          mc[itype][isetup]->DrawFitQA(true);
        }
    }  


    // -------------------------------------
    // create IAA from the MCorr objects
    // -------------------------------------
    MIaa * miaa = new MIaa(mc[0][0], mc[1][0]);
    miaa->CreateIAAdeta();
    miaa->CreateIAApta();

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


    MTools * mt = new MTools();




    // -------------------------------------
    //  PLOTTING PART
    // -------------------------------------
    int iplot = 1000;
    std::vector<TH1*>    hList;
    std::vector<TString> legList;

    // -------------------------------------
    // PLOTTING Filip Delta Eta
    // -------------------------------------
    TH1D * htmp_deta = nullptr;
    for(int itype=0; itype<2; itype++)
    {
        for(int ic=0; ic<NumCent[itype]; ic++){
            for(int iptt=minPtt; iptt<NumPtt[1]; iptt++){
                for(int ipta=2; ipta<NumPta[1]; ipta++){
                    if( PTt[1]->At(iptt) <= PTa[1]->At(ipta) )
                        continue;

                    htmp_deta = mt->RebinHistoToOther((TH1D*)mc[itype][0]->hDEtaSig[2][ic][iptt][ipta], (TH1D*)hfilip[itype]->hDEta[ic][iptt-1][ipta-2]);
                    hList={ (TH1D*)hfilip[itype]->hDEta[ic][iptt-1][ipta-2], htmp_deta };
                    legList={ "AN(2012)", "AN(2016)" };

                    MPlot * mcorrfilip = new MPlot(++iplot, "|#Delta#eta|", "1/N_{trigg.}dN/d#Delta#eta",true);

                    mcorrfilip->addHList(hList, legList);
                    mcorrfilip->SetLimitsX(0, 1.6);
                    if(itype>0) mcorrfilip->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
                    mcorrfilip->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1) ) );
                    mcorrfilip->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTa[1]->At(ipta), PTa[1]->At(ipta+1) ) );
                    mcorrfilip->Draw();
                    mcorrfilip->SetRatioLimits(0., 2.);
                    mcorrfilip->Save(Form("figs/Corr/Filip_DEta_T%d_C0%dT0%dA0%d", itype, ic,iptt,ipta));
                }
            }
        }
    }
    cout << "\n Filip DeltaEta done...\n\n";

    // -------------------------------------
    // PLOTTING Filip IAA(dEta)
    // -------------------------------------
    TH1D * htmp_iaa_fit[5];
    for(int ic=0; ic<NumCent[1]; ic++){
        for(int iptt=minPtt; iptt<NumPtt[1]; iptt++){
            for(int ipta=2; ipta<NumPta[1]; ipta++){
                if( PTt[1]->At(iptt) <= PTa[1]->At(ipta) )
                    continue;

                MPlot * miaafilip = new MPlot(++iplot, "|#Delta#eta|", "I_{AA}",false);
                hList.clear(); legList.clear();

                for(int ifit=0; ifit<5; ifit++)
                {
                    htmp_iaa_fit[ifit] = (TH1D*)mt->RebinHistoToOther( miaa->hIAA_deta[ifit][ic][iptt][ipta], hfilip[1]->hIAA[0][0][0] );
                    if(ifit==4)
                    {
                        hList.push_back( htmp_iaa_fit[ifit] );
                        legList.push_back( mc[0][0]->mfit_eta[ifit][0][iptt][ipta]->GetName() );
                    }
                }


                miaafilip->addHList(hList, legList, "pe", "p");
                miaafilip->AddInfo( mc[1][0]->BuildInfo() );
                miaafilip->AddInfo( mc[1][0]->BuildCentTitle(ic) );
                miaafilip->AddInfo( mc[1][0]->BuildPTtTitle(iptt) );
                miaafilip->AddInfo( mc[1][0]->BuildPTaTitle(ipta)  );
                miaafilip->SetLimitsXY(0, 0.28, 0.5, 2.5);
                miaafilip->Draw();
                //miaafilip->DrawThisTGraphAsymmErrors( iaafilip->gIAAsyst[ic][iptt-1][ipta-2], "E2 same", 24, 42 );
                miaafilip->DrawThisTGraphAsymmErrors( iaafilip->gIAAstat[ic][iptt-1][ipta-2], "PE same", 24, 1 );
                miaafilip->AddThisEntry(iaafilip->gIAAstat[ic][iptt-1][ipta-2], "AN(2012)");
                miaafilip->Save(Form("figs/IAA/IAA_DEta_Filip2_C0%dT0%dA0%d", ic,iptt,ipta));
            }
        }
    }
    cout << "\n Filip IAA(dEta) done...\n\n";

    // -------------------------------------
    // PLOTTING published IAA(pta)
    // -------------------------------------
    TGraphAsymmErrors * gPubIAA_c[3];  // published IAA (all 3 backg. subtraction)
    for(int it=0; it<3; it++){ gPubIAA_c[it] = publishedIAA(it);}

    for(int ic=0; ic<NumCent[1]; ic++)
    {
        int iptt = mc[0][0]->fPTt->GetBin(8.);
        // draw IAA from integral --------------------
        MPlot * mpiaa = new MPlot(iplot++, "p_{T, assoc} [GeV]", "I_{AA}", false);
        mpiaa->SetPalette(1); // set to full symbols
        hList.clear(); legList.clear();
        cout << ic << "\t" << iptt<<"\t" << endl;
        for(int ifit=0; ifit<5; ifit++)
        {
            hList.push_back( (TH1D*)miaa->hIAA_eta_INT[ifit][ic][iptt]);
            //hList.push_back( (TH1D*)miaa->hIAA_eta[ifit][ic][iptt] );
            legList.push_back( mc[0][0]->mfit_eta[ifit][0][iptt][2]->GetName() );
        }
        mpiaa->addHList(hList, legList, "PE");
        mpiaa->Draw();

        mpiaa->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
        mpiaa->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1)) );

        mpiaa->AddThisEntry(gPubIAA_c[0], "PRL data (flat)");
        mpiaa->AddThisEntry(gPubIAA_c[1], "PRL data (v_{2})");
        mpiaa->AddThisEntry(gPubIAA_c[2], "PRL data (#eta-gap)");
        mpiaa->DrawThisTGraphAsymmErrors( gPubIAA_c[0], "PE same", 24, 1 );
        mpiaa->DrawThisTGraphAsymmErrors( gPubIAA_c[1], "PE same", 25, 1 );
        mpiaa->DrawThisTGraphAsymmErrors( gPubIAA_c[2], "PE same", 26, 1 );

        mpiaa->SetLimitsXY(0,10,0.5,3.5);
        mpiaa->Save(Form("figs/IAA/IAA_INT_C0%d", ic));
    }
    cout << "\n Filip IAA(pta) done...\n\n";

    return;
/*








    // -------------------------------------
    // create IAA(pta) histos from yields 
    // extracted by MCorr for various fits 
    // and back. subtraction
    // -------------------------------------
    TH1D * hIAA_eta[NSetup[1]][kF][kC][kT];
    TH1D * hIAA_eta_INT[NSetup[1]][kF][kC][kT];
    TH1D * hIAA_phi[NSetup[1]][kF][kS][kC][kT];
    TH1D * hIAA_phi_INT[NSetup[1]][kF][kS][kC][kT];

    TH1D * hICP[2][kF][kS][kT];
    std::cout << "Processing IAA(pta)\n";

    for(int ifit=0; ifit<kF; ifit++)
    {
        for(int isetup=0; isetup<NSetup[1]; isetup++)
        {
            for(int ic=0; ic<NumCent[1]; ic++)
            {
                for(int iptt=minPtt; iptt<NumPtt[1]; iptt++)
                {
                    hIAA_eta[isetup][ifit][ic][iptt] = (TH1D*)mc[1][isetup]->hYield_eta[ifit][ic][iptt]->Clone();
                    hIAA_eta[isetup][ifit][ic][iptt]->Divide( (TH1D*)mc[0][isetup]->hYield_eta[ifit][0][iptt] );

                    hIAA_eta_INT[isetup][ifit][ic][iptt] = (TH1D*)mc[1][isetup]->hYield_eta_Int[ifit][ic][iptt]->Clone();
                    hIAA_eta_INT[isetup][ifit][ic][iptt]->Divide( (TH1D*)mc[0][isetup]->hYield_eta_Int[ifit][0][iptt] );

                    for(int isub=0; isub<kS; isub++)
                    {
                        hIAA_phi[isetup][ifit][isub][ic][iptt] = (TH1D*)mc[1][isetup]->hYield_phi[ifit][isub][ic][iptt]->Clone();
                        hIAA_phi[isetup][ifit][isub][ic][iptt]->Divide( (TH1D*)mc[0][isetup]->hYield_phi[ifit][isub][0][iptt] );

                        hIAA_phi_INT[isetup][ifit][isub][ic][iptt] = (TH1D*)mc[1][isetup]->hYield_phi_Int[ifit][isub][ic][iptt]->Clone();
                        hIAA_phi_INT[isetup][ifit][isub][ic][iptt]->Divide( (TH1D*)mc[0][isetup]->hYield_phi_Int[ifit][isub][0][iptt] );
                    }
                }
            }
        }
    }

    // -------------------------------------
    // create IAA(DEta) histos from the
    // correlation histos themselves
    // -------------------------------------
    std::cout << "Processing IAA(deta)\n";
    TH1D * hIAAdEta[NSetup[1]][kC][kT][kA];
    TH1D * hIAAdEta_newbinning[NSetup[1]][kC][kT][kA];
    TH1D * hIAAdEta_ppfit[NSetup[1]][kC][kT][kA];
    TH1D * hIAAdEta_ppfit_newbinning[NSetup[1]][kC][kT][kA];
    TH1D * hIAAdEta_ppFilip[NSetup[1]][kC][kT][kA];
    TH1D * hIAAdEta_ppFilip_newbinning[NSetup[1]][kC][kT][kA];
    int newnbins = hfilip[1]->hIAA[1][1][1]->GetNbinsX();
    double newxbins[999];
    for(int ibins=1; ibins<=newnbins; ibins++){ newxbins[ibins-1] = hfilip[1]->hIAA[1][1][1]->GetBinLowEdge(ibins); }

    for(int isetup=0; isetup<NSetup[1]; isetup++){
        for(int ic=0; ic<NumCent[1]; ic++){
            for(int iptt=minPtt; iptt<NumPtt[1]; iptt++){
                for(int ipta=0; ipta<NumPta[1]; ipta++){
                    if( PTt[1]->At(iptt) < PTa[1]->At(ipta) )
                        continue;
                    TString name = Form("hIAAdEta_S%dC%dT%dA%d",isetup,ic,iptt,ipta);
                    
                    //clone AA for IAA
                    TH1D * htmpAA = (TH1D*) mc[1][isetup]->hDEtaRealFlip[ic][iptt][ipta]->Clone(name);
                    MFit * fitAA = new MFit(0,0,htmpAA, -1.6, 1.6, false); // larger range for better backg. estimation
                    htmpAA->Fit( fitAA->ffit, "IEMRNSQ");
                    mt->subtractConstTH1(htmpAA, fitAA->ffit->GetParameter(0));

                    hIAAdEta[isetup][ic][iptt][ipta]         = (TH1D*) htmpAA->Clone(name);
                    hIAAdEta_ppfit[isetup][ic][iptt][ipta]   = (TH1D*) htmpAA->Clone(name);
                    hIAAdEta_ppFilip[isetup][ic][iptt][ipta] = (TH1D*) htmpAA->Clone(name);

                    //use pp histo for IAA
                    TH1D * htmpPP = (TH1D*) mc[0][isetup]->hDEtaRealFlip[0][iptt][ipta]->Clone(name+"tmp");
                    MFit * fitPP = new MFit(0,0,htmpPP, -1.6, 1.6, false);
                    htmpPP->Fit( fitPP->ffit, "IEMRNSQ");
                    mt->subtractConstTH1( htmpPP, fitPP->ffit->GetParameter(0));

                    hIAAdEta[isetup][ic][iptt][ipta]->Divide(htmpPP);
                    hIAAdEta_newbinning[isetup][ic][iptt][ipta] = (TH1D*)mt->NewHistoWithUniqueBins(hIAAdEta[isetup][ic][iptt][ipta], newnbins-1, newxbins);

                    //use pp fit for IAA
                    htmpPP->Reset();
                    htmpPP->Add( mc[0][isetup]->mfit_eta_ggc[0][iptt][ipta]->ffit );
                    mt->subtractConstTH1( htmpPP, mc[0][isetup]->mfit_eta_ggc[0][iptt][ipta]->ffit->GetParameter(0) );
                    htmpPP->Scale(2);
                    hIAAdEta_ppfit[isetup][ic][iptt][ipta]->Divide(htmpPP);
                    hIAAdEta_ppfit_newbinning[isetup][ic][iptt][ipta] = (TH1D*)mt->NewHistoWithUniqueBins(hIAAdEta_ppfit[isetup][ic][iptt][ipta], newnbins-1, newxbins);

                    //use pp from Filip for IAA
                    htmpPP->Reset();
                    if( ipta>=2 && iptt>=1 && PTt[1]->At(iptt) > PTa[1]->At(ipta) ) {
                        htmpPP = (TH1D*) hfilip[0]->hDEta[0][iptt-1][ipta-2]->Clone(name+"filip");
                        //rebin to ours:
                        hIAAdEta_ppFilip_newbinning[isetup][ic][iptt][ipta] = (TH1D*)mt->NewHistoWithUniqueBins(hIAAdEta_ppFilip[isetup][ic][iptt][ipta], newnbins-1, newxbins);
                        TH1D * htmpPP_new =  (TH1D*)mt->NewHistoWithUniqueBins(htmpPP, newnbins-1,newxbins);
                        hIAAdEta_ppFilip_newbinning[isetup][ic][iptt][ipta]->Divide(htmpPP_new);
                        hIAAdEta_ppFilip[isetup][ic][iptt][ipta]->Divide(htmpPP);
                    }
                    delete htmpAA;
                    delete htmpPP;
                    delete fitAA;
                    delete fitPP;
//                    cout << "NBINS: " << ourhistos[0][isetup]->hDEtaReal[ic][iptt][ipta]->GetNbinsX() << "\t" << ourhistos[1][isetup]->hDEtaReal[ic][iptt][ipta]->GetNbinsX() << endl;

                }
            }
        }
    }

    std::cout << "\n\n loading input and histos done...\n starting plotting\n";


    // -------------------------------------
    //  PLOTTING PART
    // -------------------------------------
    int iplot = 1000;
    std::vector<TH1*>    hList;  
    std::vector<TString> legList;
    TH1D * he, * hp, * hep; // eta and phi histos plotted together (+ratio)
    TH1D * htmp;


    // -------------------------------------
    // plot yield (pT, assoc) in published 
    // ptt bin 8-15 GeV
    // -------------------------------------
    
    TGraphAsymmErrors * gPubIAA_c[3];  // published IAA (all 3 backg. subtraction)
    for(int it=0; it<3; it++){ gPubIAA_c[it] = publishedIAA(it);} 

    for(int ifit=0; ifit<kF; ifit++)
    {
        int iptt = PTt[0]->GetBin(8.); // always plot IAA and ICP for 8-15GeV bin

        for(int ic=0; ic<NumCent[1]; ic++)
        {
            // draw IAA from integral --------------------
            MPlot * miaa = new MPlot(iplot++, "p_{T, assoc} [GeV]", "I_{AA}", false);
            miaa->SetPalette(1); // set to full symbols
            hList   = { hIAA_eta_INT[0][ifit][ic][iptt] };
            legList = { Form("#eta corr. |#Delta#phi|<%.1f",mc[0][0]->fPhiCut) };
            miaa->addHList(hList, legList, "PE");

            hIAA_phi_INT[0][ifit][0][ic][iptt]->SetMarkerStyle(20);
            hIAA_phi_INT[0][ifit][1][ic][iptt]->SetMarkerStyle(21);
            hIAA_phi_INT[0][ifit][2][ic][iptt]->SetMarkerStyle(22);

            miaa->Draw();
            miaa->AddInfo( fitnames[ifit].Data() );
            miaa->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
            miaa->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1)) );
            for(int isub=0; isub<kS; isub++){
                miaa->AddThisEntry( hIAA_phi_INT[0][ifit][isub][ic][iptt], ( Form("#phi corr. %s |#Delta#eta|<%.1f",subnames[isub].Data(), mc[0][0]->fEtaCut ) ) );
            }
            miaa->AddThisEntry(gPubIAA_c[0], "PRL data (flat)");
            miaa->AddThisEntry(gPubIAA_c[1], "PRL data (v_{2})");
            miaa->AddThisEntry(gPubIAA_c[2], "PRL data (#eta-gap)");
            miaa->DrawThisTGraphAsymmErrors( gPubIAA_c[0], "PE same", 24, 1 );
            miaa->DrawThisTGraphAsymmErrors( gPubIAA_c[1], "PE same", 25, 1 );
            miaa->DrawThisTGraphAsymmErrors( gPubIAA_c[2], "PE same", 26, 1 );
            miaa->DrawThisTH1(hIAA_phi_INT[0][ifit][0][ic][iptt], "PE same", 24, 2);
            miaa->DrawThisTH1(hIAA_phi_INT[0][ifit][1][ic][iptt], "PE same", 25, 3);
            miaa->DrawThisTH1(hIAA_phi_INT[0][ifit][2][ic][iptt], "PE same", 26, 4);
            hIAA_phi_INT[0][ifit][0][ic][iptt]->SetMarkerColor(kBlue+1);
            hIAA_phi_INT[0][ifit][1][ic][iptt]->SetMarkerColor(kRed+1);
            hIAA_phi_INT[0][ifit][2][ic][iptt]->SetMarkerColor(kGreen+1);
            hIAA_phi_INT[0][ifit][0][ic][iptt]->SetLineColor(kBlue+1);
            hIAA_phi_INT[0][ifit][1][ic][iptt]->SetLineColor(kRed+1);
            hIAA_phi_INT[0][ifit][2][ic][iptt]->SetLineColor(kGreen+1);

            hIAA_eta_INT[0][ifit][ic][iptt]->SetMarkerColor(kBlack);  hIAA_eta[0][ifit][ic][iptt]->SetLineColor(kBlack);
            hIAA_eta_INT[0][ifit][ic][iptt]->SetMarkerSize(1.);

            miaa->SetLimitsXY(0,10,0.5,3.5);
            miaa->Save(Form("figs/IAA/IAA_INT_FIT%d_C0%d", ifit, ic));

            // draw IAA from fit --------------------
            hList   = { hIAA_eta[0][ifit][ic][iptt] };
            legList = { Form("#eta corr. |#Delta#phi|<%.1f",mc[0][0]->fPhiCut) };
            miaa->resetHList();
            miaa->addHList(hList, legList, "PE");

            hIAA_phi[0][ifit][0][ic][iptt]->SetMarkerStyle(20);
            hIAA_phi[0][ifit][1][ic][iptt]->SetMarkerStyle(21);
            hIAA_phi[0][ifit][2][ic][iptt]->SetMarkerStyle(22);

            miaa->Draw();
            for(int isub=0; isub<kS; isub++){
                miaa->AddThisEntry( hIAA_phi[0][ifit][isub][ic][iptt], Form("#phi corr. %s \t |#Delta#eta|<%.1f",subnames[isub].Data(), mc[0][0]->fEtaCut ) );
            }
            miaa->AddThisEntry(gPubIAA_c[0], "PRL data (flat)");
            miaa->AddThisEntry(gPubIAA_c[1], "PRL data (v_{2})");
            miaa->AddThisEntry(gPubIAA_c[2], "PRL data (#eta-gap)");

            miaa->DrawThisTGraphAsymmErrors( gPubIAA_c[0], "PE same", 24, 1 );
            miaa->DrawThisTGraphAsymmErrors( gPubIAA_c[1], "PE same", 25, 1 );
            miaa->DrawThisTGraphAsymmErrors( gPubIAA_c[2], "PE same", 26, 1 );
            miaa->DrawThisTH1(hIAA_phi[0][ifit][0][ic][iptt], "PE same", 24, 2);
            miaa->DrawThisTH1(hIAA_phi[0][ifit][1][ic][iptt], "PE same", 25, 3);
            miaa->DrawThisTH1(hIAA_phi[0][ifit][2][ic][iptt], "PE same", 26, 4);
            hIAA_phi[0][ifit][0][ic][iptt]->SetMarkerColor(kBlue+1);
            hIAA_phi[0][ifit][1][ic][iptt]->SetMarkerColor(kRed+1);
            hIAA_phi[0][ifit][2][ic][iptt]->SetMarkerColor(kGreen+1);
            hIAA_phi[0][ifit][0][ic][iptt]->SetLineColor(kBlue+1);
            hIAA_phi[0][ifit][1][ic][iptt]->SetLineColor(kRed+1);
            hIAA_phi[0][ifit][2][ic][iptt]->SetLineColor(kGreen+1);

            hIAA_eta[0][ifit][ic][iptt]->SetMarkerColor(kBlack);  hIAA_eta[0][ifit][ic][iptt]->SetLineColor(kBlack);
            hIAA_eta[0][ifit][ic][iptt]->SetMarkerSize(1.);

            miaa->SetLimitsXY(0,10,0.5,3.5);
            miaa->Save(Form("figs/IAA/IAA_FIT%d_C0%d", ifit, ic));

        }

        // draw ICP ---------------------
        hICP[0][ifit][0][iptt] = (TH1D*) hIAA_eta[0][ifit][0][iptt]->Clone(Form("hICP_eta_F%dT%d",ifit,iptt));
        hICP[0][ifit][0][iptt]->Divide( hIAA_eta[0][ifit][NumCent[1]-1][iptt] );
        for(int isub=0; isub<kS; isub++)
        {
            hICP[1][ifit][isub][iptt] = (TH1D*) hIAA_phi[0][ifit][isub][0][iptt]->Clone(Form("hICP_phi_F%dS%dT%d",ifit,isub,iptt));
            hICP[1][ifit][isub][iptt]->Divide( hIAA_phi[0][ifit][isub][NumCent[1]-1][iptt] );
        }
        TString icp_label = Form( "I_{CP} (%.0f-%.0f %%/%.0f-%.0f %%)", Cent[1]->At(0), Cent[1]->At(1), Cent[1]->At(NumCent[1]-1), Cent[1]->At(NumCent[1]) );
        MPlot * micp = new MPlot(iplot++, "p_{T, assoc} [GeV]", icp_label, false);
        micp->SetPalette(1); // set to full symbols
        hList    = { hICP[0][ifit][0][iptt], hICP[1][ifit][0][iptt], hICP[1][ifit][1][iptt], hICP[1][ifit][2][iptt] };
        legList  = { "#eta "+subnames[0],"#phi "+subnames[0], "#phi "+subnames[1], "#phi "+subnames[2] };
        micp->addHList(hList, legList, "PE");
        micp->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1)) );
        micp->AddInfo( fitnames[ifit] );
        micp->Draw();
        micp->SetLimitsXY(1,10,0.5,3);
        micp->Save(Form("figs/IAA/ICP_FIT%d", ifit));
    }

    if(bPlotFilip){

        // plot IAA (delta eta)
        // ___________________________
        // Filip's histos
        // ptt: 0:4-6, 1:6-8, 2:8-15
        // pta: 0:2-3, 1:3-4, 2:4-6, 3:6-8

        // our histos:
        // ptt: 0:3-4, 1:4-6, 2:6-8, 3:8-15
        // pta: 0:0-1, 1:1-2, 2:2-3, 3:3-4, 4:4-6, 5:6-8, 6:8-10

        std::cout << "Compare I_AA with Filip's\n";

        
        int iptt_pub = PTt[0]->GetBin(8.); // always plot IAA and ICP for 8-15GeV bin
        for(int ifit=0; ifit<kF; ifit++)
        {
            MPlot * miaa_f_pta = new MPlot(iplot++, "p_{T, assoc} [GeV]", "I_{AA}", false);
            miaa_f_pta->SetPalette(1); // set to full symbols
            // draw IAA from with and compare with Filip's preliminary --------------------
            hList   = { hIAA_eta[0][ifit][0][iptt_pub], hIAA_eta_INT[0][ifit][0][iptt_pub]  };
            legList = { Form("#eta corr. |#Delta#phi|<%.1f (fit)",mc[0][0]->fPhiCut), Form("#eta corr. |#Delta#phi|<%.1f (integral)",mc[0][0]->fPhiCut)  };
            miaa_f_pta->addHList(hList, legList, "PE");
            miaa_f_pta->Draw();
            miaa_f_pta->AddThisEntry(gFilipIAAstat, "AN(2012) preliminary (R<0.2)");
            miaa_f_pta->DrawThisTGraphAsymmErrors( gFilipIAAstat, "PE same", 26, 1 );
            miaa_f_pta->DrawThisTGraphAsymmErrors( gFilipIAAsyst, "5 same", 25, 1 );
            miaa_f_pta->SetLimitsXY(0,8,0.5,3.5);
            miaa_f_pta->Save(Form("figs/IAA/IAA_FilipPrel_FIT%d", ifit));
        }

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

                    miaafilip->Save(Form("figs/IAA/IAA_DEta_Filip_C0%dT0%dA0%d", ic, iptt, ipta));

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
                    miaafilip2->Save(Form("figs/IAA/IAA_DEta_Filip2_C0%dT0%dA0%d", ic,iptt,ipta));
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
                        mdetafilip->Save(Form("figs/Corr/DEta_Filip_%s_C0%dT0%dA0%d", Types[itype].Data(),ic,iptt,ipta));

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
                    mdetaAApp->Save(Form("figs/Corr/DEta_Filip_AA-PP_C0%dT0%dA0%d", ic,iptt,ipta));

                    MPlot * mdetaAApp2 = new MPlot(++iplot, "|#Delta#eta|", "1/N_{trigg.}dN/d|#Delta#eta|", true);
                    hList = {htmp_AA, htmp_pp};
                    mdetaAApp2->addHList(hList, legList, "PE", "p");
                    mdetaAApp2->SetLimitsX(0, 1.0);
                    mdetaAApp2->AddInfo( Form("Cent: %.0f-%.0f %%",Cent[1]->At(ic), Cent[1]->At(ic+1)));
                    mdetaAApp2->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", PTt[1]->At(iptt), PTt[1]->At(iptt+1) ) );
                    mdetaAApp2->AddInfo( Form("p_{Ta}#in %.0f-%.0f GeV", PTa[1]->At(ipta), PTa[1]->At(ipta+1) ) );
                    mdetaAApp2->SetRatioLimitsXY(0, 1.0, 0.0, 2.0);
                    mdetaAApp2->Draw();
                    mdetaAApp2->Save(Form("figs/Corr/DEta_OURS_AA-PP_C0%dT0%dA0%d", ic,iptt,ipta));
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
                mwidth1->Save(Form("figs/Fit/Width_Cent_T0%dA0%d",iptt,ipta));
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
                        metaFilip->Save(Form("figs/PubYield/PubYied_Filip_%s_S0%d_C0%dT0%dA0%d", Types[itype].Data(), isetup, ic_ours, iptt, ipta));
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
            mep[ic][ipta]->Save(Form("figs/PubYield/IAA_DEta_C0%dA0%d", ic, ipta));
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
                miaa->Save(Form("figs/PubYield/IAA_DEta_Filip_C0%dT0%dA0%d", ic, iptt, ipta));
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
            miaa_filip->Save(Form("figs/PubYield/IAA_DEta_FilipVSFilip_C0%dT0%dA0%d", 0, iptt, ipta));
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
                        meta_c->Save(Form("figs/Corr/2D1D_eta_%s_C%dT%dA%d",Types[itype].Data(),ic,iptt,ipta));

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
                        mphi_c->Save(Form("figs/Corr/2D1D_phi_%s_C%dT%dA%d",Types[itype].Data(),ic,iptt,ipta));
                    }
                }
            }
        }
    */
}


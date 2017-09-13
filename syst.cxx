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
#include <map>

#include "mplot.h"
#include "mcorr.h"
#include "loadFilip.h"
#include "loadFilipIAAfinal.h"
#include "miaa.h"
#include "JFiete.h"


// -------------------------------------
//
//        M A I N     M A C R O
//
// -------------------------------------

void ProcessMC(MCorr * mc, double fitrange, TString fitopt)
{

	if(mc->OutPutExists(fitopt, fitrange)) {
		std::cout << "output generated already, skipping \n";
		return;
    } else {
		mc->Initialize();
		mc->LoadDEtaDPhi();
        mc->ProjectDEtaDPhi();
		mc->FitDEtaHistos(fitopt, fitrange); //IEMRNSQ or WLRNSQ
		mc->SaveOutput();
		delete mc;
	}
}

void syst()
{
    // -------------------------------------
    // processing correlation from JCORRAN 
    // ROOT file with mcorr.h class
    // -------------------------------------

    std::vector<bool>   vertexcorr   = {true, false};
    std::vector<double> vertexcut_pp = {-5, -10};
	std::vector<double> vertexcut_aa = {-8, -6, -4};
    std::vector<double> phiproj   = {0.15, 0.2, 0.25, 0.4};
    std::vector<double> fitrange  = {1.0, 1.2, 1.4, 1.6};
	std::vector<int> resonance = {0, 1};
	std::vector<TString> fitopt = {"IEMRNSQ", "WLRNSQ"};

    int imcorr=0;

    // -------------------------------------
    // SETUP OF DEFAULT VALUES
    const double minpt = 3; // GeV
    const double etaCut = 1.0;
    const double phiCut = 0.2 * TMath::Pi();

    //const int iTC = 0;

    const TString dataDir = "~/cernbox/Work/ALICE/IAA/data/jcorran/";

    const TString trackCut[] = {"TPCOnly", "RAA", "GlobalSDD"};
    TString inNamePP[] = {dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_pp-1269_20170330-1749-2760GeV_LHC11a_p4_AOD113_noSDD.root", // TPCOnly
                          dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_pp-1413_20170830-1640-2760GeV_LHC11a_p4_AOD113_noSDD.root", // RAA
                          dataDir+"JDiHadronCorr_legotrain_allTrigg-CF_pp-1311_20170424-1753-2760GeV_LHC11a_p4_AOD113_noSDD.root"}; // GlobalSDD


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


	MCorr * mc_pp, * mc_aa;

	bool ivcorr=true;

    for(double iphi : {0.2}) {
        for(double ifitrange : {1.6}) {
            for(int ir : {0}) {
                for(TString fopt : {"IEMRNSQ"})
                {
                    for(int iTC: {0}) {
                        for(double iv : {-10})
                        {
                            mc_pp = new MCorr( imcorr++, 0, inNamePP[iTC], "LHC11a", "113", "0", "noSDD", Form("%s_H0_T1",trackCut[iTC].Data()), iphi*TMath::Pi(), 1.0, iv, minpt, 15., ir, ivcorr, true);
                            ProcessMC( mc_pp, ifitrange, fopt );

                        }
                        for(double iv : {-8})
                        {
                            mc_aa = new MCorr( imcorr++, 1, inNameAA01[iTC], "LHC10h", "86", "3", "", Form("%s_H0_T1",trackCut[iTC].Data()), iphi*TMath::Pi(), 1.0, iv, minpt, 15., ir, ivcorr, true );
                            ProcessMC( mc_aa, ifitrange, fopt );


                        }
                    }

				}
			}
		}
	}

	return;

    // -------------------------------------
    // create and save IAA from the MCorr objects
    // -------------------------------------
    MIaa * miaa;
	for(int i=0; i<4; ++i)
	{
		miaa = new MIaa(mc_pp, mc_aa);
		miaa->CreateIAAdeta();
		miaa->CreateIAApta();
		miaa->SaveIAA();

	}
}

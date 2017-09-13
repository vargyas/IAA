#ifndef MCORR_H
#define MCORR_H

/* ***************************************************************************
 * This macro to extract I_{AA} vs. delta eta and p_{T,assoc} from ALICE data.
 * It also compares our results to published
 * http://hepdata.cedar.ac.uk/view/ins930312
 * **************************************************************************/

#include <TString.h>
#include <TFitResultPtr.h>
#include <THashList.h>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <TVirtualFitter.h>

#include <sys/stat.h>
#include <unistd.h>
#include <string>

#include "AliJHistManagerROOT6.cxx"
#include "mfit.h"
#include "mtools.h"
#include "mplot.h"


enum kType { kPP, kPbPb };

class MCorr
{
    private:
        // info strings related to I/O:
        int fMCorrCount;        // class index
        int fType;              // 0=pp, 1=PbPb
        int fResonance;         // 0: no resonance cut (e.g. merge with hResonanceCut), 1: resonance cut

        bool fLoadedDEtaDPhi, fProjectedDEtaDPhi, fFittedDEta, fFittedDPhi;
        bool fVertexCorr, fWingCorr;

		double fFitRange;

        double fsumAssocBinsForMixAbove;
        double fsumTriggBinsForMixAbove;

        double fNTrigg[6][15]; // C, T
        double fNTrigg_noVCut[6][15]; // C, T (no vertex cut)
        double fNEve[6]; // C

        TString fTypeName;    // pp, PbPb
        TString fInFileName;  // name of JCORRAN input file (assuming dir: ../data/jcorran)
        TString fOutFileName; // output ROOT file name (generated in WriteOutput() )
        TString fTrackCut;    // name of track cut
        TString fFitOpt;

        TFile * fInFile, * fOutFile;

        std::vector<TString> fFitNames;

        MTools * mt;

    public:
        int iplot=0;
        std::vector<TH1*>    hList;
        std::vector<TString> legList;

        // info string related to input file
        TString fPeriod, fAOD, fRL, fComment;
        // with array dimensions
        const static int kT = 5, // PTt
              kA = 7,            // PTa
              kC = 5,            // Cent
              kV = 11,           // Vertex
              kF = 3;            // Fit functions (Gauss+Gen.Gauss+SmallRangeGauss)

        // various cuts (-1 is all):
        double fPhiCut, fEtaCut, fVertexCut;

        // hist manager and histos from input file
        AliJHistManager * fhst;
        AliJTH1D fhDEtaNearRaw, fhDEtaNearMix, fhTriggPt, fhZVert;
        AliJTH1D fhChargedPt, fhiCentr, fhCentr, fhIetaTrigg, fhIetaAssoc, fhIphiTrigg, fhIphiAssoc;
        AliJTProfile fhTrackingEff;
        AliJTH2D fhDphiDetaPta, fhResonanceCut;
        //AliJTH1D fhDEtaNear2D, fhDEtaNearFar2D, fhDPhiNear2D, fhDPhiFar2D;
        AliJBin *fCent, *fVtx, *fPTt, *fPTa, *fEta, *fPhi;

        // AliJBin indices:
        int fNumPhi, fNumEta, fNumPta, fNumPtt, fNumVtx, fNumCent;
        int fMinPTtBin, fMaxPTtBin;
        int fVertexSkip, fPhiSkip, fEtaSkip;
        double fMinPTt, fMaxPTt;

        // correction coming from different mixed event correction:
        // Filip: range/hMix->Integral(), so average is around 1
        // Ours:  1./hMix->GetBinContent() at (0,0), so hMix is efficiency correction,
        // and its maximum should be <= 1
        // double fFilipCorr_eta[kC][kT][kA];

        // HISTOGRAM DEFINITIONS:

        // ntrigg
        TH1D * hTriggVSum[kC][kT];

        // fit results
        THashList * FitResultsEta;

		// histos for DEta
        //TH1D * hDEtaRaw[kC][kT][kA],
        TH1D * hDEtaRaw2D[kC][kT][kA],
        // * hDEtaMix[kC][kT][kA],
             * hDEtaMix2D[kC][kT][kA],
             * hDEtaReal[kC][kT][kA],
             * hDEtaReal2D[kC][kT][kA],
             * hDEtaRealFlip[kC][kT][kA],
             * hDEtaRealFlip2D[kC][kT][kA],
             * hDEtaRealFlipErr[kF][kC][kT][kA],
             * hDEtaRealFlipErr2D[kF][kC][kT][kA],
             * hDEtaSig1D[kF][kC][kT][kA],
             * hDEtaSig2D[kF][kC][kT][kA],
             * hDEtaReal2DFar[kC][kT][kA];

		TH2D * hDEtaDPhiRaw[kC][kT][kA],
			 * hDEtaDPhiMix[kC][kT][kA],
             * hDEtaDPhiReal[kC][kT][kA];

        //TH1D * hDEtaWingCorr[kC][kT][kA];

        TH1D * hChargedPtPub[kC]; // publised inclusive p_t

        // fit functions
		MFit * mfit_eta_1d[kF][kC][kT][kA],
			 * mfit_eta_2d[kF][kC][kT][kA];

        // DEFAULT CONSTRUCTOR
        MCorr() :
            fMCorrCount(0),
            fType(-1),
            fResonance(0),
            fLoadedDEtaDPhi(false),
            fProjectedDEtaDPhi(false),
            fFittedDEta(false),
            fFittedDPhi(false),
            fVertexCorr(false),
            fWingCorr(false),
            fFitRange(0),
            fsumTriggBinsForMixAbove(6.), // same as Filip
            fsumAssocBinsForMixAbove(4.), // same as Filip
            fPhiCut(0.0),
            fEtaCut(0.0),
            fVertexCut(0.0),
            fMinPTt(-1),
            fMaxPTt(-1),
            fInFileName(""),
            fPeriod(""),
            fAOD(""),
            fRL(""),
            fComment(""),
            fTrackCut(""),
            fTypeName("NOTYPE")
        {
            // default const.
			mt = new MTools();
        }

        // CONSTRUCTOR
        MCorr(int imcorr, int itype, TString inname, TString per, TString aod, TString rl, TString comment,
                TString tcut, double phicut, double etacut, double vertexcut, double minptt, double maxptt, int ires, bool vcorr, bool wcorr) :
            fMCorrCount(imcorr),
            fType(itype),
            fResonance(ires),
            fLoadedDEtaDPhi(false),
            fProjectedDEtaDPhi(false),
            fFittedDEta(false),
            fFittedDPhi(false),
            fVertexCorr(vcorr),
            fWingCorr(wcorr),
            fFitRange(0),
            fsumTriggBinsForMixAbove(6.), // same as Filip
            fsumAssocBinsForMixAbove(4.), // same as Filip
            fPhiCut(phicut),
            fEtaCut(etacut),
            fVertexCut(vertexcut),
            fMinPTt(minptt),
            fMaxPTt(maxptt),
            fInFileName(inname),
            fPeriod(per),
            fAOD(aod),
            fRL(rl),
            fComment(comment),
            fTrackCut(tcut)
        {
            switch(fType)
            {
                case kPP:   fTypeName = "pp";     break;
                case kPbPb: fTypeName = "PbPb";   break;
                default:    fTypeName = "NOTYPE"; break;
            }

            mt = new MTools();
            //TH1::SetDefaultSumw2(false);

        }

        // COPY CONSTRUCTOR
        // TODO: fix order
        MCorr(const MCorr& obj) :
            fMCorrCount(obj.fMCorrCount),
            fType(obj.fType),
            fInFileName(obj.fInFileName),
            fPeriod(obj.fPeriod),
            fAOD(obj.fAOD),
            fRL(obj.fRL),
            fComment(obj.fComment),
            fTrackCut(obj.fTrackCut),
            fPhiCut(obj.fPhiCut),
            fEtaCut(obj.fEtaCut),
            fVertexCut(obj.fVertexCut),
            fVertexCorr(obj.fVertexCorr),
            fWingCorr(obj.fWingCorr),
			fFitRange(obj.fFitRange),
            fMinPTt(obj.fMinPTt),
            fMaxPTt(obj.fMaxPTt),
            fResonance(obj.fResonance),
            fLoadedDEtaDPhi(obj.fLoadedDEtaDPhi),
            fProjectedDEtaDPhi(obj.fProjectedDEtaDPhi),
            fFittedDEta(obj.fFittedDEta),
            fFittedDPhi(obj.fFittedDPhi),
            fsumTriggBinsForMixAbove(obj.fsumTriggBinsForMixAbove),
            fsumAssocBinsForMixAbove(obj.fsumAssocBinsForMixAbove)
        {
            // copy
        }

        void Initialize()
        {
            LoadInputFile();

            LoadNtriggers();

            if(fType==kPP) {fsumTriggBinsForMixAbove=4; fsumAssocBinsForMixAbove=3;}
            std::cout <<"\n"<< Form("MCorr %d initialization done... \n",fMCorrCount);
        }

        void LoadInputFile()
        {
            std::cout << "LoadInputFile()\n";

            fInFile = TFile::Open( fInFileName );
            fInFile->cd();

            fhst = new AliJHistManager( "hst" );
			fhst->LoadConfig( Form("JDiHadronIAA_%s/AliJHistos", fTrackCut.Data()) );

            fhTriggPt     = fhst->GetTH1D("hTriggPtBin"); //!        // 3D: C, V, T
            fhZVert       = fhst->GetTH1D("hZVert");      //!        // 1D: C
            fhDphiDetaPta = fhst->GetTH2D("hDphiDetaPta"); //!      // 4D: D, C, T, A
			if(fResonance==0) fhResonanceCut = fhst->GetTH2D("hResonanceCut"); //!      // 4D: D, C, T, A

			fhChargedPt = fhst->GetTH1D("hChargedPt");   //! // 1D: C
            fhTrackingEff = fhst->GetTProfile("hTrackingEff");   //! // 1D: C
            fhiCentr    = fhst->GetTH1D("hiCentr");      //! // 0D
            fhCentr     = fhst->GetTH1D("hCentr");       //! // 0D
            fhIetaTrigg = fhst->GetTH1D("hIetaTrigg");   //! // 2D: C, T
            fhIetaAssoc = fhst->GetTH1D("hIetaAssoc");   //! // 2D: C, T
            fhIphiTrigg = fhst->GetTH1D("fhIphiTrigg");   //! // 2D: C, T
            fhIphiAssoc = fhst->GetTH1D("fhIphiAssoc");   //! // 2D: C, T

            fCent = fhst->GetBin("Cent");   fNumCent = fCent->Size();
            fVtx  = fhst->GetBin("Vtx");    fNumVtx  = fVtx->Size();
            fPTt  = fhst->GetBin("PTt");    fNumPtt  = fPTt->Size();
            fPTa  = fhst->GetBin("PTa");    fNumPta  = fPTa->Size()-1;
            fEta  = fhst->GetBin("EtaGap"); fNumEta  = fEta->Size();
            fPhi  = fhst->GetBin("PhiGap"); fNumPhi  = fPhi->Size();

            if( fPhiCut==-1 )    fPhiCut = fPhi->At( fNumPhi );
            if( fEtaCut==-1 )    fEtaCut = fEta->At( fNumEta );
            if( fVertexCut==-1 ) fVertexCut = fVtx->At( fNumVtx );

            fVertexSkip = fVtx->GetBin( fVertexCut );
			fPhiSkip    = fPhi->GetBin( fPhiCut-0.01 )+1; // this is bin, but loop will be up to bin upper edge
            fEtaSkip    = fEta->GetBin( fEtaCut )+1;
            fMinPTtBin  = fPTt->GetBin( fMinPTt );
            fMaxPTtBin  = fPTt->GetBin( fMaxPTt );

            std::cout << "LoadInputFile() done...\n";
		}


        // DESTRUCTOR
        ~MCorr()
        {
            std::cout << "deleting MCorr instances... \n";

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
                        // eta related histos
						if(fProjectedDEtaDPhi){
                            //delete hDEtaRaw[ic][iptt][ipta];
                            //delete hDEtaMix[ic][iptt][ipta];
                            delete hDEtaReal[ic][iptt][ipta];
                            delete hDEtaRealFlip[ic][iptt][ipta];

                            delete hDEtaRaw2D[ic][iptt][ipta];
                            delete hDEtaMix2D[ic][iptt][ipta];
                            delete hDEtaReal2D[ic][iptt][ipta];
                            delete hDEtaRealFlip2D[ic][iptt][ipta];

                            delete hDEtaReal2DFar[ic][iptt][ipta];


                        }

                        if(fFittedDEta)
                        {
                            for(int ifit=0; ifit<kF; ifit++)
							{
                                //delete mfit_eta_1d[ifit][ic][iptt][ipta];
                                //delete mfit_eta_2d[ifit][ic][iptt][ipta];
                                hDEtaSig1D[ifit][ic][iptt][ipta];
                                hDEtaSig2D[ifit][ic][iptt][ipta];
							}
                        }

                        if(fLoadedDEtaDPhi)
                        {
                            delete hDEtaDPhiRaw[ic][iptt][ipta];
                            delete hDEtaDPhiMix[ic][iptt][ipta];
                            delete hDEtaDPhiReal[ic][iptt][ipta];
                        }
                    }
                }
            }

            delete fhst;
            FitResultsEta->Clear("C");
			delete FitResultsEta;
			delete mt;
			fInFile->Close();
        }

        void LoadNtriggers(){
            std::cout << "LoadNtriggers()\n";

            // print ntrigger stat to separate file
            TString fOutTextFileName = fInFileName;
            fOutTextFileName.ReplaceAll("root", "txt");

            printf("saving file to: %s\n", fOutTextFileName.Data());
            ofstream outF;
            outF.open( fOutTextFileName );

            //for(int iv=0; iv<(fNumVtx/2); iv++)
            //    outF << "  &" << fhst->GetBin("Vtx")->BuildTitle(iv);
            //for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
            //    outF << std::endl << fhst->GetBin("PTt")->BuildTitle(iptt) << " & ";
            //    for(int iv=0; iv<(fNumVtx); iv++){
            //        outF << fhTriggPt[0][iv][iptt]->GetEntries() << "\t &";
            //    }
            //}
            outF << std::endl << std::endl;

            for(int ic=0; ic<fNumCent; ic++)
            {
                fhiCentr[0]->Print();
                //std::cout << endl << ic << std::endl;
                fNEve[ic]=fhiCentr[0]->GetBinContent(ic+1);
                std::cout << "No. of events in ic " << ic << " = " << fNEve[ic] << std::endl;
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                {
                    fNTrigg[ic][iptt]=0;
                    fNTrigg_noVCut[ic][iptt]=0;
                    hTriggVSum[ic][iptt]=(TH1D*)fhTriggPt[ic][0][iptt]->Clone(); hTriggVSum[ic][iptt]->Reset();
                    for(int iv=0; iv<fNumVtx; iv++){
						fNTrigg_noVCut[ic][iptt] += fhTriggPt[ic][iv][iptt]->Integral(); // count ntrigg without vertex cut (phi)
                        if( iv>=fVertexSkip && iv<(fNumVtx-fVertexSkip) )
                        {
                            hTriggVSum[ic][iptt]->Add( fhTriggPt[ic][iv][iptt] );
							fNTrigg[ic][iptt] += fhTriggPt[ic][iv][iptt]->Integral(); // count ntrigg with vertex cut (eta)
                        }
                    }
                }
            }

            double nevents = 0;
            for(int ic=0; ic<fNumCent; ic++) {
                outF <<  "  & " << fhst->GetBin("Cent")->BuildTitle(ic);
                nevents += fNEve[ic];
            }
            outF << std::endl;
            for(int ic=0; ic<fNumCent; ic++) outF << " & " << fNEve[ic];
            outF << std::endl << std::endl << nevents << std::endl << std::endl;
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                outF << std::endl << fhst->GetBin("PTt")->BuildTitle(iptt) << " & ";
                for(int ic=0; ic<fNumCent; ic++)
                    outF << " &" << fNTrigg[ic][iptt];
            }
            outF.close();
            std::cout << "LoadNtriggers() done...\n";

        }

        void DoQA()
        {
            if(fType==kPbPb) DrawCentr();

            DrawIetaIphi();

            DrawInclPt();
        }

        void DrawCentr()
        {
            double noeve = 0;
            for(int ic=0;ic<fNumCent;ic++) noeve += fNEve[ic];

            //plot centrality "centrality[%]","counts",
			MPlot * mcentr = new MPlot(iplot++, "centrality[%]", "counts", false);
            hList = { fhCentr[0] };
            legList={ Form("%.2f M",fhCentr[0]->Integral()/1000000.) };

            mcentr->addHList(hList, legList, "l");
            mcentr->AddInfo( BuildInfo() );
            mcentr->SetLimitsXY(0,100,noeve/1000.,noeve/10.);
            mcentr->SetLog(false, true);
            mcentr->Draw();
            mcentr->Save( "../figs/QA/centrality" );

            //print event stat. table
            for(int ic=0; ic<fNumCent; ic++)
            {
                std::cout << fCent->BuildTitle(ic) << "\t" << fNEve[ic] << std::endl;
            }
        }

        // TODO: phi is missing
        void DrawIetaIphi()
        {
            TString xtit = "incl. ?";
            TString ytit = "1/N_{eve.} dN/d ?";
            double integral[2][5][6];
            TString outN = "outpicture";

            for(int iep=0; iep<2; iep++) // 0: eta, 1: phi
            {
                for(int ic=0; ic<fNumCent; ic++)
                {
                    if(iep==0) xtit = "incl. #eta"; ytit = "dN/d#eta (scaled)";
                    if(iep==1) xtit = "incl. #phi"; ytit = "dN/d#phi (scaled)";
                    MPlot * mincl = new MPlot(iplot++, xtit, ytit, false);
                    hList = {}; legList = {};
                    if(iep==0)
                    {
                        for( int ipt=0; ipt<3; ipt++ )
                        {
                            integral[0][ic][ipt] = fhIetaAssoc[ic][ipt]->Integral();
                            fhIetaAssoc[ic][ipt]->Scale(1./integral[0][ic][ipt]);
                            hList.push_back(fhIetaAssoc[ic][ipt]);
                            legList.push_back(fPTa->BuildTitle(ipt));
                        } for( int ipt=0; ipt<3; ipt++ ) {
                            integral[0][ic][ipt+3] = fhIetaTrigg[ic][ipt]->Integral();
                            fhIetaTrigg[ic][ipt]->Scale(1./integral[0][ic][ipt+3]);
                            hList.push_back(fhIetaTrigg[ic][ipt]);
                            legList.push_back(fPTt->BuildTitle(ipt));
                        }
                        std::cout << "iep " << iep << " done...\n";
                    }
                    else if(iep==1)
                    {
                        for( int ipt=0; ipt<3; ipt++ )
                        {
                            fhIphiAssoc[ic][ipt]->Print();
                            integral[1][ic][ipt] = fhIphiAssoc[ic][ipt]->Integral();
                            fhIphiAssoc[ic][ipt]->Scale(1./integral[1][ic][ipt]);
                            hList.push_back(fhIphiAssoc[ic][ipt]);
                            legList.push_back(fPTa->BuildTitle(ipt));
                        } for( int ipt=0; ipt<3; ipt++ ) {
                            fhIphiTrigg[ic][ipt]->Print();
                            integral[1][ic][ipt+3] = fhIphiTrigg[ic][ipt]->Integral();
                            fhIphiTrigg[ic][ipt]->Scale(1./integral[1][ic][ipt+3]);
                            hList.push_back(fhIphiTrigg[ic][ipt]);
                            legList.push_back(fPTt->BuildTitle(ipt));
                        }
                        std::cout << "iep " << iep << " done...\n";

                    }
                    mincl->addHList(hList, legList, "f");
                    mincl->AddInfo( BuildInfo() );
                    if(fType!=kPP) mincl->AddInfo( BuildCentTitle(ic) );
                    mincl->Draw();
                    iep==0?mincl->SetLimitsY(0.011, 0.015):mincl->SetLimitsY(0., 0.01);
                    iep==0 ? outN = Form("../figs/QA/ieta_%s_C%d",fTypeName.Data(), ic) : outN = Form("../figs/QA/iphi_%s_C%d",fTypeName.Data(), ic);
                    mincl->Save( outN );
                    if(iep==0)
                    {
                        for( int ipt=0; ipt<3; ipt++ )
                        {
                            fhIetaAssoc[ic][ipt]->Scale(integral[0][ic][ipt]);
                            fhIetaTrigg[ic][ipt]->Scale(integral[0][ic][ipt+3]);
                        }
                    }
                    if(iep==1)
                    {
                        for( int ipt=0; ipt<3; ipt++ )
                        {
                            fhIphiAssoc[ic][ipt]->Scale(integral[1][ic][ipt]);
                            fhIphiTrigg[ic][ipt]->Scale(integral[1][ic][ipt+3]);
                        }
                    }
                }
            }
        }

        void DrawInclPt()
        {
            GetPublishedInclPt();
            std::cout << "DrawInclPt()\n";

            const double dEta = fEta->At(fNumEta) - fEta->At(0);
            std::cout << "DETA=" <<dEta <<std::endl;

            TString title = "";
            for(int ic=0; ic<fNumCent; ic++)
            {
                MPlot * mptt = new MPlot(iplot++, "p_{T}", "1/N_{eve} 1/(2#pip_{T}|#Delta#eta|) dN/dp_{T} [ (GeV/c)^{-2} ]", true);
                mptt->AddInfo( BuildInfo() );
                mptt->SetRatioFit(0, true, 0.6, 30);
                //fhChargedPt[ic]->Rebin(2); fhChargedPt[ic]->Scale(1./2.);
                fhChargedPt[ic]->Scale(1./2./TMath::Pi()/dEta/fNEve[ic], "width");
                mt->DivideWithX( fhChargedPt[ic] ); // correcting with 1/p_T

                //fhChargedPt[ic]->Multiply(fhTrackingEff[ic]);
                hList  = { hChargedPtPub[ic], fhChargedPt[ic] };
                if(fType==kPP) {
                    title = "arxiv:1307.1093";
                    fhChargedPt[ic]->Scale(62.2); // published data is Ed^3sigma/dp, so scaling with sigma
                    fhChargedPt[ic]->Scale(0.88); // trigger efficiency for TPCOnly
                }
                else if(fType==kPbPb) {
                    title = "arxiv:1208.2711";
                    mptt->AddInfo( BuildCentTitle(ic) );
                    //fhChargedPt[ic]->Scale(0.9); // trigger efficiency for TPCOnly
                }
                legList ={ title, "this analysis" };
                mptt->addHList(hList, legList);

				mptt->SetLimitsXY(0.6, 30, 1E-9, 1E4);
                mptt->SetRatioLimitsXY(0.6, 30., 0.8, 1.2);
                mptt->SetLog(true, true);
                mptt->Draw();
                mptt->Save( Form("../figs/QA/incl_pt_%s_C%d",fTypeName.Data(), ic));
            }
        }

        void GetPublishedInclPt()
        {
            TFile * infile;
            switch(fType)
            {
                case kPbPb:
                    infile = TFile::Open("../data/published/ALICEPbPbReference2760.root");
                    for(int ic=0; ic<fNumCent; ic++)
                    {
                        if(fCent->At(ic)==0 && fCent->At(ic+1)==5) { hChargedPtPub[ic] = (TH1D*)infile->Get("h8210_d1x1y1"); }
                        else if(fCent->At(ic)==5 && fCent->At(ic+1)==10)  { hChargedPtPub[ic] = (TH1D*)infile->Get("h8210_d2x1y1"); }
                        else if(fCent->At(ic)==0 && fCent->At(ic+1)==10)  { hChargedPtPub[ic] = (TH1D*)infile->Get("h8210_d10x1y1"); }
                        else if(fCent->At(ic)==10 && fCent->At(ic+1)==20) { hChargedPtPub[ic] = (TH1D*)infile->Get("h8210_d3x1y1"); }
                        else if(fCent->At(ic)==20 ) { hChargedPtPub[ic] = (TH1D*)infile->Get("h8210_d15x1y1"); } // upper limit missing (40)
                        else if(fCent->At(ic)==40 ) { hChargedPtPub[ic] = (TH1D*)infile->Get("h8210_d13x1y1"); } // upper limit missing (60)
                        else if(fCent->At(ic)==60)  { hChargedPtPub[ic] = (TH1D*)infile->Get("h8210_d12x1y1"); } // upper limit missing (80/90)
                    }
                    break;
                case kPP:
                    infile = TFile::Open("../data/published/ALICEppReference2760.root");
                    hChargedPtPub[0] = (TH1D*) infile->Get("h8462_d1x1y2");
                    //hChargedPtPub[0]->Print();
                    break;
            }
            std::cout << "GetPublishedInclPt() done...\n";
        }
		double GetNoEvents(int ic) {return fNEve[ic]; }

		TString GetNotUniqueCTA(int ic, int iptt, int ipta)
		{
			return Form("%s_C%dT%dA%d",fTypeName.Data(),ic,iptt,ipta);
		}
        TString GetCA(int ic, int ipta)
        {
            return Form("%s_C%dA%d_count%d",fTypeName.Data(),ic,ipta,fMCorrCount);
        }
        TString GetCT(int ic, int iptt)
        {
            return Form("%s_C%dT%d_count%d",fTypeName.Data(),ic,iptt,fMCorrCount);
        }

        TString GetCTA(int ic, int iptt, int ipta)
        {
            return Form("%s_C%dT%dA%d_count%d",fTypeName.Data(),ic,iptt,ipta,fMCorrCount);
        }

        TString GetCVTA(int ic, int iv, int iptt, int ipta)
        {
            //return Form("%s_C%dV%dT%dA%d",fTypeName.Data(),ic,iv,iptt,ipta);
            return Form("%s_C%dV%dT%dA%d_count%d",fTypeName.Data(),ic,iv,iptt,ipta,fMCorrCount);
        }


        void LoadDEtaDPhi()
        {
            if(fType==0) LoadDEtaDPhi_CombinePtAbove();
            if(fType==1) LoadDEtaDPhi_CombinePtAbove();

        }

        void LoadDEtaDPhi_PP()
        {
            std::cout << "MCorr::LoadDEtaDPhi_PP()...\n";

            int ic = 0;
            double deta = fEta->At(fNumEta);
            TString cta = "";
            TString newname = "";
            TH2D * htmp_raw = nullptr;
            TH2D * htmp[2][10][10][10]; // DVTA
            TH2D * htmp_mix[10]; // for pp only vertex dependence, rest is summed
            TH2D * htmp_mixall = nullptr;

            for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
            {
                htmp_mix[iv] = (TH2D*) fhDphiDetaPta[1][0][iv][0][0]->Clone(Form("hmix_pp_tmp_V%d",iv));
                htmp_mix[iv]->RebinX(2);
                htmp_mix[iv]->Reset();
            }
            htmp_mixall = (TH2D*) fhDphiDetaPta[1][0][fVertexSkip][0][0]->Clone("hmixall_pp_tmp");
            htmp_mixall->RebinX(2);
            htmp_mixall->Reset();

            // load raw and mixed into htmp
            // then create empty histos for raw/mixed/real
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt) < fPTa->At(ipta))
                        continue; // PTa upper border should be smaller than PTt lower

                    for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                    {

                        for(int itype=0; itype<2; itype++)
                        {
                            newname = Form("%s_tmp",fhDphiDetaPta[itype][ic][iv][iptt][ipta]->GetName());
                            htmp[itype][iv][iptt][ipta] = (TH2D*) fhDphiDetaPta[itype][ic][iv][iptt][ipta]->Clone(newname);
                            // restore original without resonance cut
                            if(fResonance==0)
                                htmp[itype][iv][iptt][ipta]->Add(fhResonanceCut[itype][ic][iv][iptt][ipta]);
                            htmp[itype][iv][iptt][ipta]->RebinX(2);
                            htmp[itype][iv][iptt][ipta]->Scale(1./2.);
                        }
                        htmp_mix[iv]->Add( htmp[1][iv][iptt][ipta] );

                    }
                    cta = GetCTA(ic,iptt,ipta);
                    hDEtaDPhiRaw[ic][iptt][ipta]  = (TH2D*) htmp[0][fVertexSkip][iptt][ipta]->Clone(Form("hDEtaDPhiRaw_%s",cta.Data()));
                    hDEtaDPhiMix[ic][iptt][ipta]  = (TH2D*) htmp[1][fVertexSkip][iptt][ipta]->Clone(Form("hDEtaDPhiMix_%s",cta.Data()));
                    hDEtaDPhiReal[ic][iptt][ipta] = (TH2D*) htmp[0][fVertexSkip][iptt][ipta]->Clone(Form("hDEtaDPhiReal_%s",cta.Data()));
                    hDEtaDPhiRaw[ic][iptt][ipta]->Reset();
                    hDEtaDPhiMix[ic][iptt][ipta]->Reset();
                    hDEtaDPhiReal[ic][iptt][ipta]->Reset();

                }
            }


            for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
            {
                htmp_mixall->Add(htmp_mix[iv]);
                htmp_mix[iv]->Scale(1./mt->GetMixedNorm2D(htmp_mix[iv], deta));
            }
            htmp_mixall->Scale(1./mt->GetMixedNorm2D(htmp_mixall, deta));

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt) < fPTa->At(ipta))
                        continue; // PTa upper border should be smaller than PTt lower

                    for(int ic=0; ic<fNumCent; ic++)
                    {
                        for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                        {
                            // Add raw correlation
                            hDEtaDPhiRaw[ic][iptt][ipta]->Add( htmp[0][iv][iptt][ipta] );
                            hDEtaDPhiMix[ic][iptt][ipta]->Add( htmp[1][iv][iptt][ipta] );
                        }
                        // Correct with mixed in vertex bins if requested
                        if(fVertexCorr)
                        {
                            for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                            {
                                htmp_raw = (TH2D*) htmp[0][iv][iptt][ipta]->Clone(Form("%s_tmp_%d",htmp[0][iv][iptt][ipta]->GetName(),fMCorrCount));
                                htmp_raw->Divide(htmp_mix[iv]);
                                hDEtaDPhiReal[ic][iptt][ipta]->Add( htmp_raw );
                            }
                        }
                        else
                        {
                            hDEtaDPhiReal[ic][iptt][ipta]->Add( hDEtaDPhiRaw[ic][iptt][ipta] );
                            hDEtaDPhiReal[ic][iptt][ipta]->Divide( htmp_mixall );
                        }
                        hDEtaDPhiReal[ic][iptt][ipta]->Scale( 1./fNTrigg[ic][iptt], "width" );
                    }
                }
            }
            fLoadedDEtaDPhi=true;
        }





        void LoadDEtaDPhi_AA()
        {
            std::cout << "MCorr::LoadDEtaDPhi_PP()...\n";

            double deta = fEta->At(fNumEta);
            TString cta = "";
            TString newname = "";
            TH2D * htmp_raw = nullptr;
            TH2D * htmp[2][5][10][10][10]; // DCVTA
            TH2D * htmp_mix[5][10]; // for PbPb vertex and centrality dependence, rest is summed
            TH2D * htmp_mixall[5];;

            for(int ic=0; ic<fNumCent; ic++)
            {
                for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                {
                    htmp_mix[ic][iv] = (TH2D*) fhDphiDetaPta[1][ic][iv][0][0]->Clone(Form("hmix_T%d_tmp_C%dV%d",fType,ic,iv));
                    htmp_mix[ic][iv]->RebinX(2);
                    htmp_mix[ic][iv]->Reset();
                }

                htmp_mixall[ic] = (TH2D*) fhDphiDetaPta[1][ic][fVertexSkip][0][0]->Clone(Form("hmixall_T%d_tmp_C%d",fType,ic));
                htmp_mixall[ic]->RebinX(2);
                htmp_mixall[ic]->Reset();
            }

            // load raw and mixed into htmp
            // then create empty histos for raw/mixed/real
            for(int ic=0; ic<fNumCent; ic++)
            {
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                {
                    for(int ipta=0; ipta<fNumPta; ipta++)
                    {
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                        {
                            for(int itype=0; itype<2; itype++)
                            {
                                newname = Form("%s_tmp",fhDphiDetaPta[itype][ic][iv][iptt][ipta]->GetName());
                                htmp[itype][ic][iv][iptt][ipta] = (TH2D*) fhDphiDetaPta[itype][ic][iv][iptt][ipta]->Clone(newname);
                                // restore original without resonance cut
                                if(fResonance==0)
                                    htmp[itype][ic][iv][iptt][ipta]->Add(fhResonanceCut[itype][ic][iv][iptt][ipta]);
                                htmp[itype][ic][iv][iptt][ipta]->RebinX(2);
                                htmp[itype][ic][iv][iptt][ipta]->Scale(1./2.);
                            }
                            htmp_mix[ic][iv]->Add( htmp[1][ic][iv][iptt][ipta] );

                        }
                        cta = GetCTA(ic,iptt,ipta);
                        hDEtaDPhiRaw[ic][iptt][ipta]  = (TH2D*) htmp[0][ic][fVertexSkip][iptt][ipta]->Clone(Form("hDEtaDPhiRaw_%s",cta.Data()));
                        hDEtaDPhiMix[ic][iptt][ipta]  = (TH2D*) htmp[1][ic][fVertexSkip][iptt][ipta]->Clone(Form("hDEtaDPhiMix_%s",cta.Data()));
                        hDEtaDPhiReal[ic][iptt][ipta] = (TH2D*) htmp[0][ic][fVertexSkip][iptt][ipta]->Clone(Form("hDEtaDPhiReal_%s",cta.Data()));
                        hDEtaDPhiRaw[ic][iptt][ipta]->Reset();
                        hDEtaDPhiMix[ic][iptt][ipta]->Reset();
                        hDEtaDPhiReal[ic][iptt][ipta]->Reset();

                    }
                }
            }
            for(int ic=0; ic<fNumCent; ic++)
            {
                for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                {
                    htmp_mixall[ic]->Add(htmp_mix[ic][iv]);
                    htmp_mix[ic][iv]->Scale(1./mt->GetMixedNorm2D(htmp_mix[ic][iv], deta));
                }

                htmp_mixall[ic]->Scale(1./mt->GetMixedNorm2D(htmp_mixall[ic], deta));

                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                {
                    for(int ipta=0; ipta<fNumPta; ipta++)
                    {
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        for(int ic=0; ic<fNumCent; ic++)
                        {
                            for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                            {
                                // Add raw correlation
                                hDEtaDPhiRaw[ic][iptt][ipta]->Add( htmp[0][ic][iv][iptt][ipta] );
                                hDEtaDPhiMix[ic][iptt][ipta]->Add( htmp[1][ic][iv][iptt][ipta] );
                            }
                            // Correct with mixed in vertex bins if requested
                            if(fVertexCorr)
                            {
                                for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                                {
                                    htmp_raw = (TH2D*) htmp[0][ic][iv][iptt][ipta]->Clone(Form("%s_tmp_%d",htmp[0][ic][iv][iptt][ipta]->GetName(),fMCorrCount));
                                    htmp_raw->Divide(htmp_mix[ic][iv]);
                                    hDEtaDPhiReal[ic][iptt][ipta]->Add( htmp_raw );
                                }
                            }
                            else
                            {
                                hDEtaDPhiReal[ic][iptt][ipta]->Add( hDEtaDPhiRaw[ic][iptt][ipta] );
                                hDEtaDPhiReal[ic][iptt][ipta]->Divide( htmp_mixall[ic] );
                            }
                            hDEtaDPhiReal[ic][iptt][ipta]->Scale( 1./fNTrigg[ic][iptt], "width" );
                        }
                    }
                }
            }
            fLoadedDEtaDPhi=true;
        }









        void LoadDEtaDPhi_CombinePtAbove()
        {
            std::cout << "MCorr::LoadDEtaDPhi()...\n";

            double scaleMix = 1.;
            TString cta = "";
            TString newname = "";
            TH2D * htmp[2][5][10][10][10];
            TH2D * htmp_raw = nullptr;
			TH2D * htmp_mix[5][10][10][10];
            //TH2D * hmix_pp[10];

            double deta = fEta->At(fNumEta);

            //TH2D * htmp_mix_pp[10]; // for pp only vertex dependence, rest is summed

            const int iptt_mixed = fPTt->GetBin( fsumTriggBinsForMixAbove );
            const int ipta_mixed = fPTa->GetBin( fsumAssocBinsForMixAbove );

			// load raw and mixed into htmp
			// then create empty histos for raw/mixed/real
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt) < fPTa->At(ipta))
                        continue; // PTa upper border should be smaller than PTt lower

                    for(int ic=0; ic<fNumCent; ic++)
                    {
                        for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                        {
							for(int itype=0; itype<2; itype++)
                            {
                                newname = Form("%s_T%d_tmp",fhDphiDetaPta[itype][ic][iv][iptt][ipta]->GetName(), fType);
                                htmp[itype][ic][iv][iptt][ipta] = (TH2D*) fhDphiDetaPta[itype][ic][iv][iptt][ipta]->Clone(newname);
								// restore original without resonance cut
                                if(fResonance==0)
                                    htmp[itype][ic][iv][iptt][ipta]->Add(fhResonanceCut[itype][ic][iv][iptt][ipta]);
                                htmp[itype][ic][iv][iptt][ipta]->RebinX(2);
                                htmp[itype][ic][iv][iptt][ipta]->Scale(1./2.);
                            }
                        }
						cta = GetCTA(ic,iptt,ipta);
						hDEtaDPhiRaw[ic][iptt][ipta]  = (TH2D*) htmp[0][ic][fVertexSkip][iptt][ipta]->Clone(Form("hDEtaDPhiRaw_%s",cta.Data()));
						hDEtaDPhiMix[ic][iptt][ipta]  = (TH2D*) htmp[1][ic][fVertexSkip][iptt][ipta]->Clone(Form("hDEtaDPhiMix_%s",cta.Data()));
						hDEtaDPhiReal[ic][iptt][ipta] = (TH2D*) htmp[0][ic][fVertexSkip][iptt][ipta]->Clone(Form("hDEtaDPhiReal_%s",cta.Data()));
						hDEtaDPhiRaw[ic][iptt][ipta]->Reset();
						hDEtaDPhiMix[ic][iptt][ipta]->Reset();
						hDEtaDPhiReal[ic][iptt][ipta]->Reset();
                    }
                }
            }
            //for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)


            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt) < fPTa->At(ipta))
                        continue; // PTa upper border should be smaller than PTt lower

                    for(int ic=0; ic<fNumCent; ic++)
                    {
                        for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                        {
                            // Add raw correlation
                            hDEtaDPhiRaw[ic][iptt][ipta]->Add( htmp[0][ic][iv][iptt][ipta] );
                            hDEtaDPhiMix[ic][iptt][ipta]->Add( htmp[1][ic][iv][iptt][ipta] ); // this should be deleted once add large pt combination is enabled

                            // Add mixed (combine large pT bin)

							if(fPTa->At(ipta)+0.1 < fsumAssocBinsForMixAbove && fPTt->At(iptt) + 0.1 < fsumTriggBinsForMixAbove )
                            {
								htmp_mix[ic][iv][iptt][ipta] = (TH2D*) htmp[1][ic][iv][iptt][ipta]->Clone(Form("%s_tmp",htmp[1][ic][iv][iptt][ipta]->GetName()));
                                hDEtaDPhiMix[ic][iptt][ipta]->Add(htmp[1][ic][iv][iptt][ipta]);
                            } else {
								htmp_mix[ic][iv][iptt][ipta] = (TH2D*) htmp[1][ic][iv][iptt][ipta]->Clone(Form("%s_tmp",htmp[1][ic][iv][iptt][ipta]->GetName()));
								htmp_mix[ic][iv][iptt][ipta]->Reset();
                                for(int it=fMinPTtBin; it<fNumPtt; it++ ) {
                                    for(int ia=0; ia<fNumPta; ia++)
                                    {
                                        if(fPTt->At(it) < fPTa->At(ia)) continue;
                                        if(fPTa->At(ia)+0.1 < fsumAssocBinsForMixAbove || fPTt->At(it) + 0.1 < fsumTriggBinsForMixAbove ) continue;
										htmp_mix[ic][iv][iptt][ipta]->Add( htmp[1][ic][iv][it][ia] );
                                        hDEtaDPhiMix[ic][iptt][ipta]->Add( htmp[1][ic][iv][it][ia] );
                                    }
                                }
                            }
                        }
                        // Correct with mixed in vertex bins if requested
                        if(fVertexCorr)
                        {
                            for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                            {
                                htmp_raw = (TH2D*) htmp[0][ic][iv][iptt][ipta]->Clone(Form("%s_T%d_tmp",htmp[0][ic][iv][iptt][ipta]->GetName(),fType));
                                //htmp_mix[ic][iv][iptt][ipta]->Scale( 1./mt->GetMixedNorm2D(htmp_mix[ic][iv][iptt][ipta]) );
                                //htmp_raw->Divide(htmp_mix[ic][iv][iptt][ipta]);
                                scaleMix = 1./mt->GetMixedNorm2D(htmp[1][ic][iv][iptt][ipta], deta);
                                htmp[1][ic][iv][iptt][ipta]->Scale(scaleMix);
                                htmp_raw->Divide(htmp[1][ic][iv][iptt][ipta]);
                                htmp[1][ic][iv][iptt][ipta]->Scale(1./scaleMix);
                                hDEtaDPhiReal[ic][iptt][ipta]->Add( htmp_raw );
                            }
                        }
                        else
                        {
                            hDEtaDPhiReal[ic][iptt][ipta]->Add( hDEtaDPhiRaw[ic][iptt][ipta] );
                            scaleMix = 1./( mt->GetMixedNorm2D(hDEtaDPhiMix[ic][iptt][ipta], deta) );
                            hDEtaDPhiMix[ic][iptt][ipta]->Scale(scaleMix);
                            hDEtaDPhiReal[ic][iptt][ipta]->Divide( hDEtaDPhiMix[ic][iptt][ipta] );
                            hDEtaDPhiMix[ic][iptt][ipta]->Scale(1./scaleMix);
                        }
                        hDEtaDPhiReal[ic][iptt][ipta]->Scale( 1./fNTrigg[ic][iptt], "width" );
                    }
                }
            }
            fLoadedDEtaDPhi=true;
        }


		void ProjectDEtaDPhi()
        {
            std::cout << "MCorr::ProjectDEtaDPhi()...\n";
			TString cta = "";
            Int_t phi_firstbin = -1;
            Int_t phi_lastbin = -1;
            Int_t phi_firstmixedbin = -1;
            Int_t phi_lastmixedbin = -1;
            Double_t firstmixed = 0.6*TMath::Pi()+0.001;
            Double_t lastmixed  = 1.4*TMath::Pi()-0.001;
            Int_t phi_firstbin_far = -1;
            Int_t phi_lastbin_far = -1;
            double deta = fEta->At(fNumEta);

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++) {
                for(int ipta=0; ipta<fNumPta; ipta++) {
                    if(fPTt->At(iptt) < fPTa->At(ipta))
                        continue; // PTa upper border should be smaller than PTt lower
                    for(int ic=0; ic<fNumCent; ic++) {
                        cta = GetCTA(ic,iptt,ipta);

						phi_firstbin = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->FindBin(-fPhiCut+0.001);
						phi_lastbin  = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->FindBin(fPhiCut-0.001);
                        phi_firstmixedbin = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->FindBin(firstmixed);
                        phi_lastmixedbin  = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->FindBin(lastmixed);

						// ------ 2D correlation

						hDEtaRaw2D[ic][iptt][ipta]  = (TH1D*) hDEtaDPhiRaw[ic][iptt][ipta]->ProjectionX( Form("hDEtaRaw2D_%s",cta.Data()), phi_firstbin, phi_lastbin,"e" );
                        hDEtaRaw2D[ic][iptt][ipta]->Scale(1./double(phi_lastbin-phi_firstbin+1));

                        // preserve correct 2D normalisation
                        Double_t scaleMix = 1./mt->GetMixedNorm2D(hDEtaDPhiMix[ic][iptt][ipta], deta);
                        hDEtaDPhiMix[ic][iptt][ipta]->Scale(scaleMix);

                        hDEtaMix2D[ic][iptt][ipta]  = (TH1D*) hDEtaDPhiMix[ic][iptt][ipta]->ProjectionX( Form("hDEtaMix2D_%s",cta.Data()), phi_firstmixedbin, phi_lastmixedbin,"e" );
                        hDEtaMix2D[ic][iptt][ipta]->Scale(1./double(phi_lastmixedbin-phi_firstmixedbin+1));
                        hDEtaDPhiMix[ic][iptt][ipta]->Scale(1./scaleMix); // restore original

                        hDEtaReal2D[ic][iptt][ipta] = (TH1D*) hDEtaDPhiReal[ic][iptt][ipta]->ProjectionX( Form("hDEtaReal2D_%s",cta.Data()), phi_firstbin, phi_lastbin,"e" );
                        hDEtaReal2D[ic][iptt][ipta]->Scale(1./2./TMath::Pi());
                        //hDEtaReal2D[ic][iptt][ipta]->Scale(1./2.);

                        // ------ prepare wing correction
						phi_firstbin_far = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->FindBin( (1.-0.4)*TMath::Pi() );
                        phi_lastbin_far  = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->FindBin( (1.+0.4)*TMath::Pi() );

                        hDEtaReal2DFar[ic][iptt][ipta] = (TH1D*) hDEtaDPhiReal[ic][iptt][ipta]->ProjectionX( Form("hDEtaReal2DFar_%s",cta.Data()), phi_firstbin_far, phi_lastbin_far,"e" );
                        hDEtaReal2DFar[ic][iptt][ipta]->Scale(1./mt->GetMixedNorm1D(hDEtaReal2DFar[ic][iptt][ipta], 8));

                        // set error manually to zero so wing correction doesn't introduce error
                        for(int ib=1; ib<=hDEtaReal2DFar[ic][iptt][ipta]->GetNbinsX(); ib++)
                            hDEtaReal2DFar[ic][iptt][ipta]->SetBinError(ib,0);

                        // wing correction only for the 2D mixed event correction method
                        if(fWingCorr) {
                            hDEtaReal2D[ic][iptt][ipta]->Divide( hDEtaReal2DFar[ic][iptt][ipta] );
                        }
                        hDEtaRealFlip2D[ic][iptt][ipta] = (TH1D*)mt->Flip( hDEtaReal2D[ic][iptt][ipta] );

						// ------ 1D correlation

						hDEtaReal[ic][iptt][ipta] = (TH1D*) hDEtaRaw2D[ic][iptt][ipta]->Clone(Form("hDEtaReal2Dmixed1D_%s",cta.Data()));
						hDEtaReal[ic][iptt][ipta]->Divide(hDEtaMix2D[ic][iptt][ipta]);

                        hDEtaReal[ic][iptt][ipta]->Scale(1./fNTrigg[ic][iptt], "width");
                        hDEtaReal[ic][iptt][ipta]->Scale(1., "width");
                        //if(fType==0) hDEtaReal[ic][iptt][ipta]->Scale(4.0);
                        hDEtaRealFlip[ic][iptt][ipta] = (TH1D*) mt->Flip(hDEtaReal[ic][iptt][ipta]);
                    }
                }
            }
            fProjectedDEtaDPhi = true;
        }

        void DrawWingCorr()
        {
            if(!fProjectedDEtaDPhi) return;

            TString fname;
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        MPlot * mwcorr = new MPlot(iplot++, "#Delta#eta", "wing correction", false);
						hDEtaReal2DFar[ic][iptt][ipta]->Rebin(2); hDEtaReal2DFar[ic][iptt][ipta]->Scale(1./2.);
                        hList = { hDEtaReal2DFar[ic][iptt][ipta] };
                        legList={""};
                        mwcorr->addHList(hList, legList, "pe");
                        mwcorr->SetLimitsX(-1.6, 1.6);
                        //mwcorr->SetRatioLimitsXY(-1.6, 1.6, 0, 2.5);
                        mwcorr->AddInfo( BuildInfo() );
                        if(fType != kPP) mwcorr->AddInfo( BuildCentTitle(ic) );
                        mwcorr->AddInfo( BuildPTtTitle(iptt) );
                        mwcorr->AddInfo( BuildPTaTitle(ipta) );
						mwcorr->SetLimitsX(-1.6,1.6);
                        mwcorr->Draw();
                        fname = Form("../figs/Corr/wingcorr_%s_%s", fTypeName.Data(), GetCTA(ic,iptt,ipta).Data());
                        mwcorr->Save( fname );
                    }
                }
            }
        }

        void DrawDEtaWingAll()
        {
            if(!fProjectedDEtaDPhi)
                return;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
                        DrawDEtaWing(ic,iptt,ipta);
                    }
                }
            }
        }

        void DrawDEtaWing(int ic, int iptt, int ipta)
        {
            // needs wing correction to be turned off, otherwise it corrects twice
            if(fWingCorr) {
                std::cerr << "Turn off wing correction to use DrawDEtaWing()\n";
                return;
            }
            TH1D * htmp_on = (TH1D*)hDEtaReal2D[ic][iptt][ipta]->Clone(Form("hwing_on_%s",hDEtaReal2D[ic][iptt][ipta]->GetName()));
            TH1D * htmp_off = (TH1D*)hDEtaReal2D[ic][iptt][ipta]->Clone(Form("hwing_off_%s",hDEtaReal2D[ic][iptt][ipta]->GetName()));
            htmp_on->Divide(hDEtaReal2DFar[ic][iptt][ipta]);

            htmp_on->Rebin(2); htmp_on->Scale(1./2.);
            htmp_off->Rebin(2); htmp_off->Scale(1./2.);

            hList   = { htmp_on, htmp_off };

            legList = {  "wing on", "wing off" };
            MPlot * mw = new MPlot(iplot++, "#Delta#eta", "1/N_{trigg.}dN/d#Delta#eta", true);

            mw->addHList(hList, legList, "ple");
            mw->SetLimitsX(-1.6, 1.6);

            mw->AddInfo( BuildInfo() );
            if(fType != kPP) mw->AddInfo( BuildCentTitle(ic) );
            mw->AddInfo( BuildPTtTitle(iptt) );
            mw->AddInfo( BuildPTaTitle(ipta) );
            mw->Draw();

            mw->Save( Form("../figs/Corr/wingcorr_corr_%s_%s", fTypeName.Data(), GetCTA(ic,iptt,ipta).Data()) );

        }

        // -------------------------------------
        // Plot 2D histograms (without MPlot, as
        // TH2::Draw is not implemented yet)
        // -------------------------------------
        void Draw2DHistos()
        {
            TCanvas * c[100];
            int iplot2d = 1000;
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<1/*fNumCent*/; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        c[iplot2d] = new TCanvas(Form("canvasfor2dnear%d",iplot2d),"c",800,700);
                        c[iplot2d]->SetTheta(23);
                        c[iplot2d]->SetPhi(-60);
                        hDEtaDPhiMix[ic][iptt][ipta]->SetAxisRange(-1.5, 1.5,"x");
                        hDEtaDPhiMix[ic][iptt][ipta]->SetAxisRange(-0.45, 1.45,"y");
                        //hDEtaDPhiReal[ic][iptt][ipta]->Draw("surf2");
                        hDEtaDPhiMix[ic][iptt][ipta]->Draw("colz");

                        c[iplot2d]->SaveAs(Form("../figs/Corr/2D_Near_%s_C%dT%dA%d.pdf",fTypeName.Data(),ic,iptt,ipta));
                        ++iplot2d;

                        //c[iplot2d] = new TCanvas(Form("canvasfor2dfar%d",iplot2d),"c",800,700);
                        //hDEtaDPhiReal[ic][iptt][ipta]->SetAxisRange(0.45, 1.4,"y");
                        //hDEtaDPhiReal[ic][iptt][ipta]->Draw("surf5");
                        //c[iplot2d]->SaveAs(Form("../figs/Corr/2D_Far_%s_C%dT%dA%d.pdf",fTypeName.Data(),ic,iptt,ipta));
                        //++iplot2d;

                    }
                }
            }
        }
        void DrawRawMixed()
        {
			if(!fProjectedDEtaDPhi)
                return;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        TH1D * htemp = (TH1D*) hDEtaMix2D[ic][iptt][ipta]->Clone(Form("hmixtemp_%s",hDEtaMix2D[ic][iptt][ipta]->GetName()));
                        //mt->scaleToThisTail(htemp, hDEtaRaw2D[ic][iptt][ipta], 1.,1.6);
                        //hList = {  hDEtaRaw2D[ic][iptt][ipta], htemp };
                        //legList = {"raw", "mixed"};

                        hList = { htemp };
                        legList = {"mixed"};

                        MPlot * mrawm = new MPlot(iplot++, "#Delta#eta", "counts", false);
                        mrawm->addHList(hList, legList, "pe");
                        mrawm->AddInfo( BuildInfo() );

                        if(fType != kPP) mrawm->AddInfo( BuildCentTitle(ic) );
                        mrawm->AddInfo( BuildPTtTitle(iptt) );
                        mrawm->AddInfo( BuildPTaTitle(ipta) );
                        double y1 = hDEtaMix2D[ic][iptt][ipta]->GetBinContent(hDEtaMix2D[ic][iptt][ipta]->FindBin(-0.001));
                        double y2 = hDEtaMix2D[ic][iptt][ipta]->GetBinContent(hDEtaMix2D[ic][iptt][ipta]->FindBin(0.001));
                        mrawm->AddInfo( Form("Mixed at 0 = %.1f, %.1f", y1, y2) );
                        mrawm->Draw();
                        TString fname = Form("../figs/Corr/eta_rawmix_%s_%s", fTypeName.Data(), GetCTA(ic,iptt,ipta).Data());
                        mrawm->Save( fname );
                    }
                }
            }
        }

        // check mixed event symmetry for opposite-side vertex bins
        void DrawMixedVtx()
        {
            TH1D * htmp, * h1, * h2;
            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip)/2; iv++) {
                            h1 = (TH1D*) fhDphiDetaPta[1][ic][iv][iptt][ipta]->ProjectionX(Form("hmvert_%s",GetCVTA(ic,iv,iptt,ipta).Data()));
                            h2 = (TH1D*) fhDphiDetaPta[1][ic][fNumVtx-iv-1][iptt][ipta]->ProjectionX(Form("hmvert_%s",GetCVTA(ic,fNumVtx-iv-1,iptt,ipta).Data()));
                            htmp = (TH1D*)h1->Clone(Form("hallmixedvert%s",GetCVTA(ic,fNumVtx-iv-1,iptt,ipta).Data()));
                            htmp->Add(h2);
                            htmp->Scale(1./2.);

                            hList = {  h1, h2 };
                            legList = {  fVtx->BuildTitle(iv), fVtx->BuildTitle(fNumVtx-iv-1) };

                            MPlot * meta = new MPlot(iplot++, "#Delta#eta", "mixed event", true);

                            meta->addHList(hList, legList, "pe");
                            meta->SetLimitsX(-1.6, 1.6);

                            meta->AddInfo( BuildInfo() );
                            if(fType != kPP) meta->AddInfo( BuildCentTitle(ic) );
                            meta->AddInfo( BuildPTtTitle(iptt) );
                            meta->AddInfo( BuildPTaTitle(ipta) );
                            meta->Draw();

                            meta->Save( Form("../figs/Test/mixed_vert/mix_%s", GetCVTA(ic,iv,iptt,ipta).Data()) );
                        }
                    }
                }
            }
        }




        // check mixed event's ptt and pta dependence:
        // conclusion: matter is Pb-Pb, can be merged in pp
        void DrawMixed()
        {
            if(!fProjectedDEtaDPhi)
                return;

            //new TCanvas(); hDEtaDPhiMix[0][2][2]->Draw("surf5");

            TH1D * hall[5][10];
            TH1D * hmix[5][10][10];

            int minpttbin=fPTt->GetBin(6.);
            int minptabin=fPTa->GetBin(4.);
            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){

                    hall[ic][iptt] = (TH1D*)hDEtaMix2D[ic][iptt][0]->Clone(Form("hallC%dT%d",ic,iptt));
                    hall[ic][iptt]->Reset();

                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
                        hmix[ic][iptt][ipta] = (TH1D*)hDEtaMix2D[ic][iptt][ipta]->Clone(Form("%s_tmp",hDEtaMix2D[ic][iptt][ipta]->GetName()));
                        hall[ic][iptt]->Add((TH1D*)hDEtaMix2D[ic][iptt][ipta]);
                        hmix[ic][iptt][ipta]->Scale(1./mt->GetMixedNorm1D(hmix[ic][iptt][ipta]));

                    }
                    hall[ic][iptt]->Scale(1./mt->GetMixedNorm1D(hall[ic][iptt]));

                    hList   = { hall[ic][iptt] };
                    legList = { "all" };

                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
                        hList.push_back(hmix[ic][iptt][ipta]);
                        legList.push_back(fPTa->BuildTitle(ipta));
                    }


                    MPlot * meta = new MPlot(iplot++, "#Delta#eta", "mixed event", true);

                    meta->addHList(hList, legList, "pe");
                    meta->SetLimitsX(-1.6, 1.6);

                    meta->AddInfo( BuildInfo() );
                    if(fType != kPP) meta->AddInfo( BuildCentTitle(ic) );
                    meta->AddInfo( BuildPTtTitle(iptt) );
                    meta->Draw();

                    meta->Save( Form("../figs/Test/mixed_merge/mix_%s", GetCT(ic,iptt).Data()) );
                }


                for(int ipta=0; ipta<fNumPta; ipta++){
                    hall[ic][ipta] = (TH1D*)hDEtaMix2D[ic][0][0]->Clone(Form("hallC%dA%d",ic,ipta));
                    hall[ic][ipta]->Reset();
                    for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
                        hmix[ic][iptt][ipta] = (TH1D*)hDEtaMix2D[ic][iptt][ipta]->Clone(Form("%s_tmp",hDEtaMix2D[ic][iptt][ipta]->GetName()));
                        hall[ic][ipta]->Add((TH1D*)hDEtaMix2D[ic][iptt][ipta]);
                        hmix[ic][iptt][ipta]->Scale(1./mt->GetMixedNorm1D(hmix[ic][iptt][ipta]));

                    }
                    hall[ic][ipta]->Scale(1./mt->GetMixedNorm1D(hall[ic][ipta]));

                    hList   = { hall[ic][ipta] };
                    legList = { "all" };

                    for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
                        hList.push_back(hmix[ic][iptt][ipta]);
                        legList.push_back(fPTt->BuildTitle(iptt));
                    }


                    MPlot * meta2 = new MPlot(iplot++, "#Delta#eta", "mixed event", true);

                    meta2->addHList(hList, legList, "pe");
                    meta2->SetLimitsX(-1.6, 1.6);

                    meta2->AddInfo( BuildInfo() );
                    if(fType != kPP) meta2->AddInfo( BuildCentTitle(ic) );
                    meta2->AddInfo( BuildPTaTitle(ipta) );
                    meta2->Draw();

                    meta2->Save( Form("../figs/Test/mixed_merge/mix_%s", GetCA(ic,ipta).Data()) );

                }
            }
        }


        void DrawDEtaFlip()
        {
			if(!fProjectedDEtaDPhi)
                return;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        hList   = { hDEtaReal[ic][iptt][ipta] };
                        legList = {  Form("|#Delta#phi|<%.1f", fPhiCut) };
                        MPlot * meta = new MPlot(iplot++, "#Delta#eta", "1/N_{trigg.}dN/d#Delta#eta", false);

                        meta->addHList(hList, legList, "pe");
                        meta->SetLimitsX(-1.6, 1.6);

                        meta->AddInfo( BuildInfo() );
                        if(fType != kPP) meta->AddInfo( BuildCentTitle(ic) );
                        meta->AddInfo( BuildPTtTitle(iptt) );
                        meta->AddInfo( BuildPTaTitle(ipta) );
                        meta->Draw();

                        meta->Save( Form("../figs/Corr/eta_sym_%s", GetCTA(ic,iptt,ipta).Data()) );
                    }
                }
            }
        }

        void DrawCorr2D1D(bool save, int which)
        {
			if(!fProjectedDEtaDPhi)
                return;

            TString xtit, ytit, fname;
            double xlimit = 0;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        switch( which )
                        {
                            case kDEta:
                                //htmp = mt->RebinHistoToOther( (TH1D*)hDEtaSig[ic][iptt][ipta], (TH1D*)hDEtaSig2D[ic][iptt][ipta] );
                                xtit =  "#Delta#eta"; ytit = "1/N_{trigg.}dN/d#Delta#eta";
                                fname = Form("../figs/Corr/eta_2D1D_%s_%s", fTypeName.Data(), GetCTA(ic,iptt,ipta).Data());

								hList   = {  hDEtaSig1D[0][ic][iptt][ipta], hDEtaSig2D[0][ic][iptt][ipta] };
                                legList = {  "1D", "2D" };
                                xlimit = 1.4;
                                break;
                                /*
                            case kDPhi:
                                htmp =  mt->RebinHistoToOther( (TH1D*)hDPhiReal[0][ic][iptt][ipta], (TH1D*)hDPhiReal2D[0][ic][iptt][ipta] );
                                xtit =  "#Delta#phi/#pi"; ytit =  "1/N_{trigg.}dN/d#Delta#phi/#pi";
                                fname = Form("figs/Corr/phi_2D1D_%s_%s", fTypeName.Data(), GetCTA(ic,iptt,ipta).Data());

                                mfit_sub0 = new MFit(0,1,htmp, -0.45, 0.45, false);
                                htmp->Fit(mfit_sub0->ffit, opt);
                                mt->subtractConstTH1(htmp, mfit_sub0->ffit->GetParameter(0));

                                htmp_2D = (TH1D*) hDPhiReal2D[0][ic][iptt][ipta]->Clone();
                                htmp_2D->Fit(mfit_sub0->ffit, opt);
                                mt->subtractConstTH1(htmp_2D, mfit_sub0->ffit->GetParameter(0));


                                hList   = {  htmp, htmp_2D };
                                legList = {  "1D histo", "2D histo (2D corr.)" };
                                xlimit = 0.45;
                                break;
                                */
                        }

                        MPlot * mcorr = new MPlot(iplot++, xtit, ytit, true);
                        mcorr->addHList(hList, legList, "pe");
                        mcorr->SetLimitsX(0, xlimit);
                        mcorr->SetRatioLimitsXY(0, xlimit, 0, 2.5);
                        mcorr->AddInfo( BuildInfo() );
                        if(fType != kPP) mcorr->AddInfo( BuildCentTitle(ic) );
                        mcorr->AddInfo( BuildPTtTitle(iptt) );
                        mcorr->AddInfo( BuildPTaTitle(ipta) );
                        mcorr->Draw();
                        if(save) mcorr->Save( fname );
                    }
                }
            }
        }


		bool checkChi2NDF(TF1 * f, double limit)
		{
            return ((f->GetChisquare())/((double)f->GetNDF()) <= limit);
		}

        void FitDEtaHistos( TString opt="EMRNS", double fitmax = 1.6 )
		{
			fFitOpt = opt;
			std::cout << "MCorr::FitDEtaHistos\n";
			fFitRange = fitmax;
            TVirtualFitter::SetDefaultFitter("Minuit");

			FitResultsEta = new THashList( 50000 );
			FitResultsEta->SetOwner(kTRUE);

            CreateFitContainers( FitResultsEta );

			TFitResultPtr r1[kF][kC][kT][kA]; // 1D
			TFitResultPtr r2[kF][kC][kT][kA]; // 2D

			double fitmin = 0;
            int fitcount = 0;

            for(int ic=0; ic<fNumCent; ic++)
            {
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                {
                    for(int ipta=0; ipta<fNumPta; ipta++)
                    {
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        int ifit = 0;

                        cout << ic << "\t" << iptt << "\t" << ipta << endl;
                        cout << hDEtaRealFlip[ic][iptt][ipta]->Integral() << endl;
						// ---------------------- Gauss + Constant

                        mfit_eta_1d[ifit][ic][iptt][ipta] = new MFit(fitcount++, kOneGenGaussConst,kDEta,hDEtaRealFlip[ic][iptt][ipta],0,fitmax,true);
						r1[ifit][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, opt);
                        //hDEtaRealFlipErr[ifit][ic][iptt][ipta] = (TH1D*)hDEtaRealFlip[ic][iptt][ipta]->Clone(Form("%sF%d_err",hDEtaRealFlip[ic][iptt][ipta]->GetName(),ifit));
                        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(mfit_eta_1d[ifit][ic][iptt][ipta]->hDEtaErr);

                        mfit_eta_2d[ifit][ic][iptt][ipta]  = new MFit(fitcount++, kOneGenGaussConst,kDEta,hDEtaRealFlip2D[ic][iptt][ipta], 0, fitmax, true);
                        r2[ifit][ic][iptt][ipta] = hDEtaRealFlip2D[ic][iptt][ipta]->Fit( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, opt);
                        //hDEtaRealFlipErr2D[ifit][ic][iptt][ipta] = (TH1D*)hDEtaRealFlip2D[ic][iptt][ipta]->Clone(Form("%sF%d_err",hDEtaRealFlip2D[ic][iptt][ipta]->GetName(),ifit));
                        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(mfit_eta_2d[ifit][ic][iptt][ipta]->hDEtaErr);
                        ++ifit;


						// ---------------------- Generalized Gauss + Constant

                        mfit_eta_1d[ifit][ic][iptt][ipta] = new MFit(fitcount++, kOneGenGaussConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax, false);
						if( checkChi2NDF(mfit_eta_1d[0][ic][iptt][ipta]->ffit, 3.) ) // init from Gaussian fit if converged
							mfit_eta_1d[ifit][ic][iptt][ipta]->ffit->SetParameters( mfit_eta_1d[0][ic][iptt][ipta]->ffit->GetParameters() );
						r1[ifit][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, opt );
                        //hDEtaRealFlipErr[ifit][ic][iptt][ipta] = (TH1D*)hDEtaRealFlip[ic][iptt][ipta]->Clone(Form("%sF%d_err",hDEtaRealFlip[ic][iptt][ipta]->GetName(),ifit));
                        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(mfit_eta_1d[ifit][ic][iptt][ipta]->hDEtaErr);

                        mfit_eta_2d[ifit][ic][iptt][ipta] = new MFit(fitcount++, kOneGenGaussConst,kDEta,hDEtaRealFlip2D[ic][iptt][ipta], 0, fitmax, false);
                        if( checkChi2NDF(mfit_eta_2d[0][ic][iptt][ipta]->ffit, 3.) ) // init from Gaussian fit if converged
                            mfit_eta_2d[ifit][ic][iptt][ipta]->ffit->SetParameters( mfit_eta_2d[0][ic][iptt][ipta]->ffit->GetParameters() );
                        r2[ifit][ic][iptt][ipta] = hDEtaRealFlip2D[ic][iptt][ipta]->Fit( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, opt );
                        //hDEtaRealFlipErr2D[ifit][ic][iptt][ipta] = (TH1D*)hDEtaRealFlip2D[ic][iptt][ipta]->Clone(Form("%sF%d_err",hDEtaRealFlip2D[ic][iptt][ipta]->GetName(),ifit));
                        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(mfit_eta_2d[ifit][ic][iptt][ipta]->hDEtaErr);
                        ++ifit;


                        // ----------------------  Small range Gauss + constant

                        fitmin = 2.*mfit_eta_1d[0][ic][iptt][ipta]->GetWidth();

                        mfit_eta_1d[ifit][ic][iptt][ipta]  = new MFit(fitcount++, kOneSmallGaussConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], fitmin, fitmax, true);
                        //if( checkChi2NDF(mfit_eta_1d[0][ic][iptt][ipta]->ffit, 3.) ) // init from Gaussian fit if converged
                        //    mfit_eta_1d[ifit][ic][iptt][ipta]->ffit->SetParameters( mfit_eta_1d[0][ic][iptt][ipta]->ffit->GetParameters() );
                        r1[ifit][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, opt);
                        //hDEtaRealFlipErr[ifit][ic][iptt][ipta] = (TH1D*)hDEtaRealFlip[ic][iptt][ipta]->Clone(Form("%sF%d_err",hDEtaRealFlip[ic][iptt][ipta]->GetName(),ifit));
                        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(mfit_eta_1d[ifit][ic][iptt][ipta]->hDEtaErr);

                        mfit_eta_2d[ifit][ic][iptt][ipta]  = new MFit(fitcount++, kOneSmallGaussConst,kDEta,hDEtaRealFlip2D[ic][iptt][ipta], fitmin, fitmax, true);
                        //if( checkChi2NDF(mfit_eta_2d[0][ic][iptt][ipta]->ffit, 3.) ) // init from Gaussian fit if converged
                        //    mfit_eta_2d[ifit][ic][iptt][ipta]->ffit->SetParameters( mfit_eta_2d[0][ic][iptt][ipta]->ffit->GetParameters() );
                        r2[ifit][ic][iptt][ipta] = hDEtaRealFlip2D[ic][iptt][ipta]->Fit( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, opt);
                        //hDEtaRealFlipErr2D[ifit][ic][iptt][ipta] = (TH1D*)hDEtaRealFlip2D[ic][iptt][ipta]->Clone(Form("%sF%d_err",hDEtaRealFlip2D[ic][iptt][ipta]->GetName(),ifit));
                        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(mfit_eta_2d[ifit][ic][iptt][ipta]->hDEtaErr);
                        ++ifit;

                        // Kaplan + Constant
                        //mfit_eta_1d[ifit][ic][iptt][ipta] = new MFit(kKaplanConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax);
                        //r1[2][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, opt ); //EMRNSQ
                        //mfit_eta_2d[ifit][ic][iptt][ipta] = new MFit(kKaplanConst,kDEta,hDEtaReal2DFlip[ic][iptt][ipta], 0, fitmax);
                        //r2[2][ic][iptt][ipta] = hDEtaReal2DFlip[ic][iptt][ipta]->Fit( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, opt ); //EMRNSQ
                        //++ifit;

                        // Cauchy + Constant
                        //mfit_eta_1d[ifit][ic][iptt][ipta]  = new MFit(kCauchyConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax);
                        //r1[3][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, "RNSQ" ); // EMRNSQ
                        //mfit_eta_2d[ifit][ic][iptt][ipta]  = new MFit(kCauchyConst,kDEta,hDEtaReal2DFlip[ic][iptt][ipta], 0, fitmax);
                        //r2[3][ic][iptt][ipta] = hDEtaReal2DFlip[ic][iptt][ipta]->Fit( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, "RNSQ" ); // EMRNSQ
                        //++ifit;

                        // Two Gauss + Constant
                        //mfit_eta_1d[ifit][ic][iptt][ipta] = new MFit(kTwoGenGaussConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax, true);
                        //r1[ifit][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, opt ); //EMRNSQ
                        //mfit_eta_2d[ifit][ic][iptt][ipta] = new MFit(kTwoGenGaussConst,kDEta,hDEtaReal2DFlip[ic][iptt][ipta], 0, fitmax, true);
                        //r2[ifit][ic][iptt][ipta] = hDEtaReal2DFlip[ic][iptt][ipta]->Fit( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, opt ); //EMRNSQ
                        //++ifit;

                        for(int i=0; i<ifit; ++i)
                        {
							if(ic==0 && iptt==fMinPTtBin && ipta==0) {fFitNames.push_back(mfit_eta_1d[i][ic][iptt][ipta]->GetName());}
							hDEtaSig1D[i][ic][iptt][ipta] = mt->subtractConstTH1( hDEtaRealFlip[ic][iptt][ipta], mfit_eta_1d[i][ic][iptt][ipta]->ffit );
                            hDEtaSig2D[i][ic][iptt][ipta] = mt->subtractConstTH1( hDEtaRealFlip2D[ic][iptt][ipta], mfit_eta_2d[i][ic][iptt][ipta]->ffit );
						}
                    }
                }
            }

			FillAllFitHistos(1, FitResultsEta, mfit_eta_1d, r1);
			FillAllFitHistos(2, FitResultsEta, mfit_eta_2d, r2);

            fFittedDEta=true;

        }

        void CreateFitContainers(THashList * results)
        {
            std::cout << "MCorr::CreateFitContainers()...\n";

            double ptaBO[fNumPta+1];
            for(int i=0;i<(fNumPta+1);i++) {
                ptaBO[i] = fPTa->At(i);
            }
            double centBO[fNumCent+1];
            for(int i=0;i<(fNumCent+1);i++) {
                centBO[i] = fCent->At(i);
            }
            TString name;

            for(int id=1; id<3; id++)
            {
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                {
                    for(int ifit = 0;ifit<kF; ifit++)
                    {
                        for(int ic=0; ic<fNumCent; ic++)
                        {
							name = GetNotUniqueFitResultKeyCent(id,ifit,ic,iptt);

                            results->Add( new TH1D(Form("hYield_%s",name.Data()), "", fNumPta, ptaBO) );
                            results->Add( new TH1D(Form("hYield_INT_%s",name.Data()),"", fNumPta, ptaBO) );
                            results->Add( new TH1D(Form("hWidth_%s",name.Data()),"", fNumPta, ptaBO) );
                            results->Add( new TH1D(Form("hBackg_%s",name.Data()),"", fNumPta, ptaBO) );
                            results->Add( new TH1D(Form("hChi2NDF_%s",name.Data()),"", fNumPta, ptaBO) );
                            results->Add( new TH1D(Form("hExpo_%s",name.Data()),"", fNumPta, ptaBO) );
                        }
                        for(int ipta=0; ipta<fNumPta; ipta++)
                        {
                            if(fPTt->At(iptt) < fPTa->At(ipta))
                                continue; // PTa upper border should be smaller than PTt lower

							name = GetNotUniqueFitResultKeyPta(id,ifit,iptt,ipta);
                            results->Add( new TH1D(Form("hYield_%s",name.Data()), "", fNumCent, centBO) );
                            results->Add( new TH1D(Form("hYield_INT_%s",name.Data()),"", fNumCent, centBO) );
                            results->Add( new TH1D(Form("hWidth_%s",name.Data()),"", fNumCent, centBO) );
                            results->Add( new TH1D(Form("hBackg_%s",name.Data()),"", fNumCent, centBO) );
                            results->Add( new TH1D(Form("hChi2NDF_%s",name.Data()),"", fNumCent, centBO) );
                            results->Add( new TH1D(Form("hExpo_%s",name.Data()),"", fNumCent, centBO) );
                        }
                    }
                }
            }
        }

        void FillAllFitHistos(int id, THashList * results, MFit * feta[kF][kC][kT][kA], TFitResultPtr r[5][kC][kT][kA])
        {
			std::cout << "MCorr::FillAllFitHistos()...\n";

            TString name = "";
            double val = 0;
            double valerr = 0;

            // fill pta dependent histos
            for(int ifit = 0;ifit<kF; ifit++)
            {
                for(int ic=0; ic<fNumCent; ic++)
                {
                    for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                    {
						name = GetNotUniqueFitResultKeyCent(id,ifit,ic,iptt);
                        TH1D * hyield = (TH1D*) results->FindObject( Form("hYield_%s",name.Data()) );
                        TH1D * hyield_int = (TH1D*) results->FindObject( Form("hYield_INT_%s",name.Data()) );
                        TH1D * hwidth = (TH1D*) results->FindObject( Form("hWidth_%s",name.Data()) );
                        TH1D * hbackg = (TH1D*) results->FindObject( Form("hBackg_%s",name.Data()) );
                        TH1D * hchi2ndf = (TH1D*) results->FindObject( Form("hChi2NDF_%s",name.Data()) );
                        TH1D * hexpo = (TH1D*) results->FindObject( Form("hExpo_%s",name.Data()) );
                        for(int ipta=0; ipta<fNumPta; ipta++)
                        {
                            if(fPTt->At(iptt) < fPTa->At(ipta))
                                continue; // PTa upper border should be smaller than PTt lower

                            if(r[ifit][ic][iptt][ipta])
                            {
                                hyield->SetBinContent(ipta+1, feta[ifit][ic][iptt][ipta]->GetYield());
                                hyield->SetBinError(ipta+1, feta[ifit][ic][iptt][ipta]->GetYieldError(r[ifit][ic][iptt][ipta]) );

                                int int_binmin = hDEtaSig1D[ifit][ic][iptt][ipta]->FindBin(0.0);
                                int int_binmax = hDEtaSig1D[ifit][ic][iptt][ipta]->FindBin(fFitRange);
                                if(id==1) val = hDEtaSig1D[ifit][ic][iptt][ipta]->IntegralAndError(int_binmin, int_binmax, valerr, "width" );
                                if(id==2) val = hDEtaSig2D[ifit][ic][iptt][ipta]->IntegralAndError(int_binmin, int_binmax, valerr, "width" );
                                hyield_int->SetBinContent(ipta+1, val);
                                hyield_int->SetBinError(ipta+1, valerr);
                                hyield_int->Scale(2.);

                                hwidth->SetBinContent(ipta+1, feta[ifit][ic][iptt][ipta]->GetWidth() );
                                hwidth->SetBinError(ipta+1, feta[ifit][ic][iptt][ipta]->GetWidthError(r[ifit][ic][iptt][ipta]) );

                                hbackg->SetBinContent(ipta+1, feta[ifit][ic][iptt][ipta]->ffit->GetParameter(0) );
                                hbackg->SetBinError(ipta+1, feta[ifit][ic][iptt][ipta]->ffit->GetParError(0) );

                                val = feta[ifit][ic][iptt][ipta]->ffit->GetChisquare();
                                val /= (double)feta[ifit][ic][iptt][ipta]->ffit->GetNDF();
                                hchi2ndf->SetBinContent(ipta+1,  val);
                                hchi2ndf->SetBinError(ipta+1, 0);

                                if(ifit==1) {

                                    val = feta[ifit][ic][iptt][ipta]->GetExpo();
                                    valerr = feta[ifit][ic][iptt][ipta]->GetExpoError();
                                    hexpo->SetBinContent(ipta+1, val);
                                    hexpo->SetBinError(ipta+1, valerr);
                                }
                            }
                        }
                    }
                }
            }
            // fill centrality dependent histos

            for(int ifit = 0;ifit<kF; ifit++)
            {
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                {
                    for(int ipta=0; ipta<fNumPta; ipta++)
                    {
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
						name = GetNotUniqueFitResultKeyPta(id,ifit,iptt,ipta);
                        TH1D * hyield = (TH1D*) results->FindObject( Form("hYield_%s",name.Data()) );
                        TH1D * hyield_int = (TH1D*) results->FindObject( Form("hYield_INT_%s",name.Data()) );
                        TH1D * hwidth = (TH1D*) results->FindObject( Form("hWidth_%s",name.Data()) );
                        TH1D * hbackg = (TH1D*) results->FindObject( Form("hBackg_%s",name.Data()) );
                        TH1D * hchi2ndf = (TH1D*) results->FindObject( Form("hChi2NDF_%s",name.Data()) );
                        TH1D * hexpo = (TH1D*) results->FindObject( Form("hExpo_%s",name.Data()) );

                        for(int ic=0; ic<fNumCent; ic++)
                        {
                            if(r[ifit][ic][iptt][ipta])
                            {
                                hyield->SetBinContent(ic+1,  feta[ifit][ic][iptt][ipta]->GetYield() );
                                hyield->SetBinError(ic+1, feta[ifit][ic][iptt][ipta]->GetYieldError(r[ifit][ic][iptt][ipta]));

                                int int_binmin = hDEtaSig1D[ifit][ic][iptt][ipta]->FindBin(0.0);
                                int int_binmax = hDEtaSig1D[ifit][ic][iptt][ipta]->FindBin(0.5);
                                if(id==1) val = hDEtaSig1D[ifit][ic][iptt][ipta]->IntegralAndError(int_binmin, int_binmax, valerr, "width" );
                                if(id==2) val = hDEtaSig2D[ifit][ic][iptt][ipta]->IntegralAndError(int_binmin, int_binmax, valerr, "width" );
                                hyield_int->SetBinContent(ic+1, val);
                                hyield_int->SetBinError(ic+1, valerr);

                                hwidth->SetBinContent(ic+1, feta[ifit][ic][iptt][ipta]->GetWidth() );
                                hwidth->SetBinError(ic+1, feta[ifit][ic][iptt][ipta]->GetWidthError(r[ifit][ic][iptt][ipta]));

                                hbackg->SetBinContent(ic+1, feta[ifit][ic][iptt][ipta]->ffit->GetParameter(0) );
                                hbackg->SetBinError(ic+1, feta[ifit][ic][iptt][ipta]->ffit->GetParError(0) );

                                val = feta[ifit][ic][iptt][ipta]->ffit->GetChisquare();
                                val /= (double)feta[ifit][ic][iptt][ipta]->ffit->GetNDF();
                                hchi2ndf->SetBinContent(ic+1,  val);
                                hchi2ndf->SetBinError(ic+1, 0);

                                if(ifit==1) {
                                    val = feta[ifit][ic][iptt][ipta]->GetExpo();
                                    valerr = feta[ifit][ic][iptt][ipta]->GetExpoError();
                                    hexpo->SetBinContent(ic+1, val);
                                    hexpo->SetBinError(ic+1, valerr);
                                }
                            }

                        }
                    }
                }
            }
        }

        void WriteOutput(){
            // create output
            // fPhiCut replace "." -> print decimals
            TString outdir = "../data/processed";
            TString phi = Form("%.1f", fPhiCut);
            phi.ReplaceAll(".", "dot");
            if(fType == kPbPb) {
                fOutFileName = Form("%s/corr_%s_%sAOD%sRL%s_%s_PhiCut%sVSkip%.0f.root", outdir.Data(), fTypeName.Data(), fPeriod.Data(), fAOD.Data(), fRL.Data(), fTrackCut.Data(), phi.Data(), fVertexCut);
            }
            if(fType == kPP) {
                fOutFileName = Form("%s/corr_%s_%sAOD%s_%s_PhiCut%s.root", outdir.Data(), fTypeName.Data(), fPeriod.Data(), fAOD.Data(), fTrackCut.Data(), phi.Data());
            }

            fOutFile = new TFile( fOutFileName, "RECREATE" );
            fOutFile->cd();
            printf( "Writing output to: %s\n", fOutFileName.Data());

            // write ntrigg
            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                    hTriggVSum[ic][iptt]->Write(Form("hTrigg_C0%dT0%d", ic, iptt));
                }
            }
            fOutFile->Write();
            fOutFile->Close();
        }

        // Basic plotting will take place here
        // more elaborate plots from other macro
        void DrawDEtaFlipAll(int id)
        {
			if(!fProjectedDEtaDPhi)
                return;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
                        DrawDEtaFlip(id,ic,iptt,ipta);
                    }
                }
            }
        }

        void DrawDEtaFlip(int id, int ic, int iptt, int ipta)
        {
            if(id==1) hList   = { hDEtaRealFlip[ic][iptt][ipta] };
            if(id==2) hList   = { hDEtaRealFlip2D[ic][iptt][ipta] };
            legList = {  Form("%dD |#Delta#phi|<%.1f", id, fPhiCut) };
            MPlot * meta = new MPlot(iplot++, "|#Delta#eta|", "1/N_{trigg.}dN/d|#Delta#eta|", false);

            meta->addHList(hList, legList, "pe");
            meta->SetLimitsX(0, 1.6);

            meta->AddInfo( BuildInfo() );
            if(fType != kPP) meta->AddInfo( BuildCentTitle(ic) );
            meta->AddInfo( BuildPTtTitle(iptt) );
            meta->AddInfo( BuildPTaTitle(ipta) );
            meta->Draw();

            if(fFittedDEta)
            {
                meta->fPad->cd();
                for(int ifit=0; ifit<3; ifit++)
                {
                    if(id==1)
                    {
                        mfit_eta_1d[ifit][ic][iptt][ipta]->Draw();
						meta->AddThisEntry( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, mfit_eta_1d[ifit][ic][iptt][ipta]->GetName(), "l");
						meta->AddThisEntry((TObject*)0, Chi2NDF(mfit_eta_1d[ifit][ic][iptt][ipta]->ffit),"");

                        //if(ifit!=2){
                        //	hDEtaRealFlipErr[ifit][ic][iptt][ipta]->SetFillColorAlpha(mfit_eta_1d[ifit][ic][iptt][ipta]->ffit->GetLineColor(), 0.5);
                        //    hDEtaRealFlipErr[ifit][ic][iptt][ipta]->Draw("e3 same");
                        //}
                    }
                    if(id==2)
                    {
                        mfit_eta_2d[ifit][ic][iptt][ipta]->Draw();
                        meta->AddThisEntry( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, mfit_eta_2d[ifit][ic][iptt][ipta]->GetName(), "l");
                        meta->AddThisEntry((TObject*)0, Chi2NDF(mfit_eta_2d[ifit][ic][iptt][ipta]->ffit),"");

                        //meta->DrawThisTH1((TH1D*)hDEtaRealFlip2D[ic][iptt][ipta], "PE same", 24, 601);
                        //if(ifit!=2) {
                        //    hDEtaRealFlipErr2D[ifit][ic][iptt][ipta]->SetFillColorAlpha(mfit_eta_2d[ifit][ic][iptt][ipta]->ffit->GetLineColor(), 0.5);
                        //if(ifit==0) hDEtaRealFlipErr2D[ifit][ic][iptt][ipta]->SetFillStyle(3004);
                        //if(ifit==1) hDEtaRealFlipErr2D[ifit][ic][iptt][ipta]->SetFillStyle(3007);
                        //hDEtaRealFlipErr2D[ifit][ic][iptt][ipta]->Draw("e3 same");
                        //}
                    }
                }
            }
            meta->Save( Form("../figs/Corr/etaflip_%dD_%s", id, GetCTA(ic,iptt,ipta).Data()) );
        }


        void DrawDEtaAll(int id)
        {
            if(!fProjectedDEtaDPhi)
                return;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
                        DrawDEta(id,ic,iptt,ipta);
                    }
                }
            }
        }

        void DrawDEta(int id, int ic, int iptt, int ipta)
        {
            if(id==1) hList   = { hDEtaReal[ic][iptt][ipta] };
            if(id==2) hList   = { hDEtaReal2D[ic][iptt][ipta] };
            legList = {  Form("%dD |#Delta#phi|<%.1f", id, fPhiCut) };
            MPlot * meta = new MPlot(iplot++, "#Delta#eta", "1/N_{trigg.}dN/d#Delta#eta", false);

            meta->addHList(hList, legList, "pe");
            meta->SetLimitsX(-1.6, 1.6);

            meta->AddInfo( BuildInfo() );
            if(fType != kPP) meta->AddInfo( BuildCentTitle(ic) );
            meta->AddInfo( BuildPTtTitle(iptt) );
            meta->AddInfo( BuildPTaTitle(ipta) );
            meta->Draw();

            meta->Save( Form("../figs/Corr/eta_%dD_%s", id, GetCTA(ic,iptt,ipta).Data()) );
        }



		void DrawFitQA(int id=1)
        {
            DrawFitWidth(id);
            DrawFitYield(id);
            DrawFitBackg(id);
            DrawFitQuality(id);
            DrawFitExpo(id);

            //DrawDataFitRatios(id);
        }
		void DrawDataFitRatios()
        {
            if(!fFittedDEta)
                return;

			double fmin = mfit_eta_1d[0][0][fMinPTtBin][0]->fitmin;
			double fmax = mfit_eta_1d[0][0][fMinPTtBin][0]->fitmax;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
                        TH1D * htmp_ggc_ratio = (TH1D*) hDEtaRealFlip[ic][iptt][ipta]->Clone("htmp_ggc_ratio");
                        TH1D * htmp_gc_ratio = (TH1D*) hDEtaRealFlip[ic][iptt][ipta]->Clone("htmp_gc_ratio");
                        // since cloning doesn't work for this constructor, copying things manually...
                        MFit * mfit_gc  = new MFit(-1, 0,0,htmp_gc_ratio,  fmin, fmax, true);
                        MFit * mfit_ggc = new MFit(-2, 0,0,htmp_ggc_ratio, fmin, fmax, false);
						mfit_gc->ffit->SetParameters(mfit_eta_1d[0][ic][iptt][ipta]->ffit->GetParameters());
						mfit_ggc->ffit->SetParameters(mfit_eta_1d[1][ic][iptt][ipta]->ffit->GetParameters());

                        htmp_gc_ratio->Divide((TF1*) mfit_gc->ffit);
                        htmp_ggc_ratio->Divide((TF1*) mfit_ggc->ffit);

                        MPlot * meta_r = new MPlot(iplot++, "#Delta#eta", "data/fit", false);

                        hList   = { htmp_gc_ratio, htmp_ggc_ratio };
                        legList = { "Gauss", "Generalized Gauss" };
                        meta_r->addHList(hList,legList);
                        meta_r->AddInfo( BuildInfo() );
                        if(fType==kPbPb) meta_r->AddInfo( BuildCentTitle(ic) );
                        meta_r->AddInfo( BuildPTtTitle(iptt));
                        meta_r->AddInfo( BuildPTaTitle(ipta));
                        meta_r->SetLimitsXY(0, 1.2, 0.8, 1.2);
                        meta_r->Draw();
                        meta_r->Save( Form("../figs/Corr/eta_ratio_%s", GetCTA(ic,iptt,ipta).Data()) );
                    }
                }
            }
        }




        // draws yield
		void DrawFitYield(int id=1)
        {
			DrawFitYieldPta(id);
			if(fType!=kPP) DrawFitYieldCent(id);
        }
		void DrawFitYieldPta(int id=1)
        {
            TString name = "";
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ic=0; ic<fNumCent; ic++)
                {
                    hList.clear(); legList.clear();
                    for(int ifit=0; ifit<2; ifit++)
                    {
                        name = GetNotUniqueFitResultKeyCent(id,ifit,ic,iptt);
                        hList.push_back( (TH1D*) FitResultsEta->FindObject( Form("hYield_%s",name.Data()) ) );
						legList.push_back( mfit_eta_1d[ifit][0][iptt][0]->GetName() );
                    }
                    MPlot * myp_p = new MPlot(++iplot, "p_{T, assoc} [GeV]", "yield", false);
                    myp_p->addHList(hList, legList, "PE");
                    myp_p->AddInfo( BuildInfo() );
                    if(fType==kPbPb) myp_p->AddInfo( BuildCentTitle(ic) );
                    myp_p->AddInfo( BuildPTtTitle(iptt));
                    //myp_p->SetLimitsXY(0, 10, 0, 5);
                    myp_p->Draw();
					myp_p->Save(Form("../figs/Fit/Yield_pta_%s_%dD_T0%dC0%d",fTypeName.Data(), id, iptt, ic));
                }
            }
        }
		void Draw_Fit_INT_YieldPta(int ifit=0)
		{
			TString name = "";
			for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
			{
				for(int ic=0; ic<fNumCent; ic++)
				{
					//hList.clear(); legList.clear();

						name = GetNotUniqueFitResultKeyCent(1,ifit,ic,iptt);
						hList = {(TH1D*) FitResultsEta->FindObject( Form("hYield_%s",name.Data()) ),  (TH1D*) FitResultsEta->FindObject( Form("hYield_INT_%s",name.Data()) ) };
						legList = { mfit_eta_1d[ifit][0][iptt][0]->GetName(), "bin counting" };

					MPlot * myp_p = new MPlot(++iplot, "p_{T, assoc} [GeV]", "yield", false);
					myp_p->addHList(hList, legList, "PE");
					myp_p->AddInfo( BuildInfo() );
					if(fType==kPbPb) myp_p->AddInfo( BuildCentTitle(ic) );
					myp_p->AddInfo( BuildPTtTitle(iptt));
					//myp_p->SetLimitsXY(0, 10, 0, 5);
					myp_p->Draw();
					myp_p->Save(Form("../figs/Fit/Yield_pta_%s_T0%dC0%d",fTypeName.Data(), iptt, ic));
				}
			}
		}
		void DrawFitYieldCent(int id=1)
        {
            TString name = "";
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt)<fPTa->At(ipta))
                        continue;

                    hList.clear(); legList.clear();
                    for(int ifit=0; ifit<2; ifit++) {
						name = GetNotUniqueFitResultKeyPta(id,ifit,iptt,ipta);
                        TH1D * h = (TH1D*) FitResultsEta->FindObject( Form("hYield_%s",name.Data()) );
                        hList.push_back( h );
						legList.push_back( mfit_eta_1d[ifit][0][iptt][0]->GetName() );
                    }
                    MPlot * myp_c = new MPlot(++iplot, "centrality [%]", "yield", false);
                    myp_c->addHList(hList, legList, "PE");
                    //myp_c->SetLimitsXY(0, 90, 0, 5);
                    myp_c->AddInfo( BuildInfo() );
                    myp_c->AddInfo( BuildPTtTitle(iptt));
                    myp_c->AddInfo( BuildPTaTitle(ipta));
                    myp_c->Draw();
					myp_c->Save(Form("../figs/Fit/Yield_cent_%s_%dD_T0%dA0%d",fTypeName.Data(), id, iptt, ipta));
                }
            }
        }

        // draws exponent of the gen.gaussian
        void DrawFitExpo(int id=1)
        {
            if(fType!=kPP) DrawFitExpoCent(id);
            DrawFitExpoPta(id);
        }
        void DrawFitExpoPta(int id=1)
        {
            TString name = "";
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ic=0; ic<fNumCent; ic++)
                {
                    hList.clear(); legList.clear();

                    name = GetNotUniqueFitResultKeyCent(id,1,ic,iptt);
                    hList.push_back( (TH1D*) FitResultsEta->FindObject( Form("hExpo_%s",name.Data()) ) );
                    legList.push_back( mfit_eta_1d[1][0][iptt][0]->GetName() );

                    MPlot * mep_p = new MPlot(++iplot, "p_{T, assoc} [GeV]", "exponent (fit)", false);
                    mep_p->addHList(hList, legList, "PE");
                    mep_p->SetLimitsXY(0, 10, 0, 3);
                    mep_p->AddInfo( BuildInfo() );
                    if(fType==kPbPb) mep_p->AddInfo( BuildCentTitle(ic) );
                    mep_p->AddInfo( BuildPTtTitle(iptt));
                    mep_p->Draw();
                    mep_p->Save(Form("../figs/Fit/Expo_pta_%s_%dD_T0%dC0%d",fTypeName.Data(), id, iptt, ic));
                }
            }
        }
        void DrawFitExpoCent(int id=1)
        {
            TString name = "";
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt)<fPTa->At(ipta))
                        continue;

                    hList.clear(); legList.clear();
                    name = GetNotUniqueFitResultKeyPta(id,1,iptt,ipta);
                    TH1D * h = (TH1D*) FitResultsEta->FindObject( Form("hExpo_%s",name.Data()) );
                    hList.push_back( h );
                    legList.push_back( mfit_eta_1d[1][0][iptt][0]->GetName() );

                    MPlot * mwp_c = new MPlot(++iplot, "centrality [%]", "exponent (fit)", false);
                    mwp_c->addHList(hList, legList, "PE");
                    mwp_c->SetLimitsXY(0, 90, 0, 3);
                    mwp_c->AddInfo( BuildInfo() );
                    mwp_c->AddInfo( BuildPTtTitle(iptt));
                    mwp_c->AddInfo( BuildPTaTitle(ipta));
                    mwp_c->Draw();
                    mwp_c->Save(Form("../figs/Fit/Expo_cent_%s_%dD_T0%dA0%d",fTypeName.Data(), id, iptt, ipta));
                }
            }
        }
        // draws width
		void DrawFitWidth(int id=1)
        {
            if(fType!=kPP) DrawFitWidthCent(id);
			DrawFitWidthPta(id);
        }

		void DrawFitWidthPta(int id=1)
        {
            TString name = "";
            TString legend = "";
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ic=0; ic<fNumCent; ic++)
                {
                    hList.clear(); legList.clear();
                    for(int ifit=0; ifit<2; ifit++)
                    {
						name = GetNotUniqueFitResultKeyCent(id,ifit,ic,iptt);
                        hList.push_back( (TH1D*) FitResultsEta->FindObject( Form("hWidth_%s",name.Data()) ) );
                        legList.push_back( mfit_eta_1d[ifit][0][iptt][0]->GetName() );
                    }
                    MPlot * mwp_p = new MPlot(++iplot, "p_{T, assoc} [GeV]", "#sigma (fit)", false);
                    mwp_p->addHList(hList, legList, "PE");
                    //myp_p->SetLimitsXY(0, 10, 0, 5);
                    mwp_p->AddInfo( BuildInfo() );
                    if(fType==kPbPb) mwp_p->AddInfo( BuildCentTitle(ic) );
                    mwp_p->AddInfo( BuildPTtTitle(iptt));
                    mwp_p->Draw();
					mwp_p->Save(Form("../figs/Fit/Width_pta_%s_%dD_T0%dC0%d",fTypeName.Data(), id, iptt, ic));
                }
            }
        }
		void DrawFitWidthCent(int id=1)
        {
            TString name = "";
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt)<fPTa->At(ipta))
                        continue;

                    hList.clear(); legList.clear();
                    for(int ifit=0; ifit<2; ifit++) {
						name = GetNotUniqueFitResultKeyPta(id,ifit,iptt,ipta);
                        TH1D * h = (TH1D*) FitResultsEta->FindObject( Form("hWidth_%s",name.Data()) );
                        hList.push_back( h );
						legList.push_back( mfit_eta_1d[ifit][0][iptt][0]->GetName() );
                    }
                    MPlot * mwp_c = new MPlot(++iplot, "centrality [%]", "#sigma (fit)", false);
                    mwp_c->addHList(hList, legList, "PE");
                    //myp_c->SetLimitsXY(0, 90, 0, 5);
                    mwp_c->AddInfo( BuildInfo() );
                    mwp_c->AddInfo( BuildPTtTitle(iptt));
                    mwp_c->AddInfo( BuildPTaTitle(ipta));
                    mwp_c->Draw();
					mwp_c->Save(Form("../figs/Fit/Width_cent_%s_%dD_T0%dA0%d",fTypeName.Data(), id, iptt, ipta));
                }
            }
        }

        // drawing constant background of the fit
		void DrawFitBackg(int id=1)
        {
			DrawFitBackgPta(id);
            if(fType!=kPP) DrawFitBackgCent(id);
        }
		void DrawFitBackgPta(int id=1)
        {
            TString name = "";
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ic=0; ic<fNumCent; ic++)
                {
                    hList.clear(); legList.clear();
                    for(int ifit=0; ifit<2; ifit++)
                    {
						name = GetNotUniqueFitResultKeyCent(id,ifit,ic,iptt);
                        hList.push_back( (TH1D*) FitResultsEta->FindObject( Form("hBackg_%s",name.Data()) ) );
						legList.push_back( mfit_eta_1d[ifit][0][iptt][0]->GetName() );
                    }
                    MPlot * mbp_p = new MPlot(++iplot, "p_{T, assoc} [GeV]", "background (fit)", false);
                    mbp_p->addHList(hList, legList, "PE");
                    mbp_p->SetLimitsX(0, fPTt->At(iptt));
                    mbp_p->AddInfo( BuildInfo() );
                    if(fType==kPbPb) mbp_p->AddInfo( BuildCentTitle(ic) );
                    mbp_p->AddInfo( BuildPTtTitle(iptt));
                    mbp_p->Draw();
					mbp_p->Save(Form("../figs/Fit/Backg_pta_%s_%dD_T0%dC0%d",fTypeName.Data(), id, iptt, ic));
                }
            }
        }
		void DrawFitBackgCent(int id=1)
        {
            TString name = "";
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt)<fPTa->At(ipta))
                        continue;

                    hList.clear(); legList.clear();
                    for(int ifit=0; ifit<2; ifit++) {
						name = GetNotUniqueFitResultKeyPta(id,ifit,iptt,ipta);
						TH1D * h = (TH1D*) FitResultsEta->FindObject( Form("hBackg_%s",name.Data()) );
                        hList.push_back( h );
						legList.push_back( mfit_eta_1d[ifit][0][iptt][0]->GetName() );
                    }
                    MPlot * mbp_c = new MPlot(++iplot, "centrality [%]", "background (fit)", false);
                    mbp_c->addHList(hList, legList, "PE");
                    //myp_c->SetLimitsXY(0, 90, 0, 5);
                    mbp_c->AddInfo( BuildInfo() );
                    mbp_c->AddInfo( BuildPTtTitle(iptt));
                    mbp_c->AddInfo( BuildPTaTitle(ipta));
                    mbp_c->Draw();
					mbp_c->Save(Form("../figs/Fit/Backg_cent_%s_%dD_T0%dA0%d",fTypeName.Data(),id, iptt, ipta));
                }
            }
        }
        TH1D * GetFitResultCent(TString what, int id,int ifit,int iptt,int ipta)
        {
            // what == hBackg, hYield, hWidth
            TString name = GetFitResultKeyPta(id,ifit,iptt,ipta);
            return (TH1D*) FitResultsEta->FindObject( Form("%s_%s",what.Data(), name.Data()) );
        }
		TH1D * GetNotUniqueFitResultCent(TString what, int id,int ifit,int iptt,int ipta)
		{
			// what == hBackg, hYield, hWidth
			TString name = GetNotUniqueFitResultKeyPta(id,ifit,iptt,ipta);
			return (TH1D*) FitResultsEta->FindObject( Form("%s_%s",what.Data(), name.Data()) );
		}

        // chi2/ndf plots
        void DrawFitQuality(int id=1)
        {
            TString name = "";
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ic=0; ic<fNumCent; ic++)
                {
                    hList.clear(); legList.clear();
                    for(int ifit=0; ifit<2; ifit++)
                    {
                        name = GetNotUniqueFitResultKeyCent(id,ifit,ic,iptt);
                        hList.push_back( (TH1D*) FitResultsEta->FindObject( Form("hChi2NDF_%s",name.Data()) ) );
						legList.push_back( mfit_eta_1d[ifit][0][iptt][0]->GetName() );
                    }
                    MPlot * mcp_p = new MPlot(++iplot, "p_{T, assoc} [GeV]", "#chi^{2}/NDF (fit)", false);
                    mcp_p->addHList(hList, legList, "PE");
                    mcp_p->SetLimitsXY(0, 10, 0, 5);
                    mcp_p->AddInfo( BuildInfo() );
                    if(fType==kPbPb) mcp_p->AddInfo( BuildCentTitle(ic) );
                    mcp_p->AddInfo( BuildPTtTitle(iptt));
                    mcp_p->Draw();
                    mcp_p->Save(Form("../figs/Fit/Chi2NDF_pta_%s_%dD_T0%dC0%d",fTypeName.Data(),id, iptt, ic));
                }
            }
        }


        TString GetFitResultKeyPta(int id, int ifit, int iptt, int ipta) {
            return Form("eta_cent_%s_%dD_F%d_T%dA%d_count%d",fTypeName.Data(),id,ifit,iptt,ipta,fMCorrCount);
        }
		TString GetNotUniqueFitResultKeyPta(int id, int ifit, int iptt, int ipta) {
			return Form("eta_cent_%s_%dD_F%d_T%dA%d",fTypeName.Data(),id,ifit,iptt,ipta);
		}
        TString GetFitResultKeyCent(int id, int ifit, int ic, int iptt) {
            return Form("eta_pta_%s_%dD_F%d_C%dT%d_count%d",fTypeName.Data(),id,ifit,ic,iptt,fMCorrCount);
        }
		TString GetNotUniqueFitResultKeyCent(int id, int ifit, int ic, int iptt) {
			return Form("eta_pta_%s_%dD_F%d_C%dT%d",fTypeName.Data(),id,ifit,ic,iptt);
		}
        TString BuildCentTitle(int ic)
        {
            return Form("Cent: %.0f-%.0f %%",fCent->At(ic), fCent->At(ic+1));
        }
        TString BuildPTtTitle(int iptt)
        {
            return Form("p_{Tt}#in %.0f-%.0f GeV", fPTt->At(iptt), fPTt->At(iptt+1));
        }
        TString BuildPTaTitle(int ipta)
        {
            if(ipta==0) return Form("p_{Ta}#in %.1f-%.0f GeV", fPTa->At(ipta), fPTa->At(ipta+1));
            else return Form("p_{Ta}#in %.0f-%.0f GeV", fPTa->At(ipta), fPTa->At(ipta+1));
        }
        TString BuildInfo()
        {
            return Form("%s AOD%s %s",fPeriod.Data(),fAOD.Data(),fTypeName.Data());
        }
        TString Chi2NDF(TF1 * f)
        {
            double chi2 = f->GetChisquare();
            int ndf     = f->GetNDF();
            return Form("%.2f/%d=%.2f",chi2, ndf,chi2/double(ndf));
        }

		TString GenerateOutputName()
		{
			return Form("%s_%s_R%d_v%.0f_p%.0f_fit%.0f_%s",fTypeName.Data(),fTrackCut.Data(), fResonance, TMath::Abs(fVertexCut),100*fPhiCut, fFitRange*10, fFitOpt.Data());
		}

		void SaveOutput()
        {
			TString outfilename = Form("../data/syst/%s.root",GenerateOutputName().Data());
			std::cout << "MCorr::SaveOutput() to " << outfilename << std::endl;

			TFile * foutfile = new TFile(outfilename,"RECREATE","title");
			TString name = "";
            // correlation 1D-2D and wing correction
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt)<fPTa->At(ipta))
                        continue;
                    for(int ic=0;ic<fNumCent; ic++)
                    {
						name = GetNotUniqueCTA(ic,iptt,ipta);
						hDEtaRealFlip[ic][iptt][ipta]->Write(Form("hReal_1D_%s",name.Data()));
                        hDEtaRealFlip2D[ic][iptt][ipta]->Write(Form("hReal_2D_%s",name.Data()));
						hDEtaReal2DFar[ic][iptt][ipta]->Write(Form("hReal2DFar_%s",name.Data()));

						for(int ifit=0; ifit<kF; ifit++)
						{
							hDEtaSig1D[ifit][ic][iptt][ipta]->Write(Form("hSig_1D_F%d_%s",ifit, name.Data()));
							hDEtaSig2D[ifit][ic][iptt][ipta]->Write(Form("hSig_2D_F%d_%s",ifit, name.Data()));
                        }

                    }
                }
            }

            // fit results
			FitResultsEta->Write();

			foutfile->Write();
            foutfile->Close();

        }
		// check if given setup has been run before by checking output directory for generated name. if true, one can
		// skip processing again. Arguments are needed as these are traditionally set after loading, during fitting.
		bool OutPutExists(TString fitopt, double fitrange){
			fFitOpt = fitopt;
			fFitRange = fitrange;
			TString outfilename = Form("../data/syst/%s.root",GenerateOutputName().Data());
			const std::string& name = std::string(outfilename.Data());
			struct stat buffer;
			return (stat (name.c_str(), &buffer) == 0);
		}

};

#endif /* MCORR_H */

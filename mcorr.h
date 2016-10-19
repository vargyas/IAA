#ifndef MCORR_H
#define MCORR_H

/*
 * ***************************************************************************
 * This macro compares our results to published
 * http://hepdata.cedar.ac.uk/view/ins930312
 * correction in vertex bins can be chosen
 * ***************************************************************************
 */

#include <TString.h>
#include <TFitResultPtr.h>
#include <THashList.h>
#include<iostream>
#include<fstream>

#include "AliJHistManagerROOT6.cxx"
#include "mfit.h"
#include "mtools.h"
#include "mplot.h"
//#include "test_pythonic_range.cxx"


// Overloaded helpers for fit histogram filling
/*
void FillFitHistos( TH1D * hY, int ib, MFit * mfit, TFitResultPtr r )
{
    hY->SetBinContent( ib, mfit->GetYield() );
    hY->SetBinError( ib, mfit->GetYieldError(r) );
}
void FillFitHistos( TH1D * hY, TH1D * hW, int ib, MFit * mfit, TFitResultPtr r )
{
    hY->SetBinContent( ib, mfit->GetYield() );
    hY->SetBinError( ib, mfit->GetYieldError(r) );
    hW->SetBinContent( ib, mfit->GetWidth() );
    hW->SetBinError( ib, mfit->GetWidthError(r) );
}
void FillFitHistos( TH1D * hY, TH1D * hW, TH1D * hE, int ib, MFit * mfit, TFitResultPtr r )
{
    hY->SetBinContent( ib, mfit->GetYield() );
    hY->SetBinError( ib, mfit->GetYieldError(r) );
    hW->SetBinContent( ib, mfit->GetWidth() );
    hW->SetBinError( ib, mfit->GetWidthError(r) );
    hE->SetBinContent( ib, mfit->GetExpo() );
    hE->SetBinError( ib, mfit->GetExpoError() );
}

void FillIntHistos( TH1D * hY, TH1D * hCorr, int ib, MFit * mfit)
{
    double val	  = 0;
    double valerr = 0;
    TH1D * h = (TH1D*) hCorr->Clone();

    // subtract constant before integrating
    TF1 * fUE = mfit->GetUE();
    h->Add(fUE, -1.0);
    delete fUE;

    int int_binmin = 0, int_binmax = 1;
    switch(mfit->GetHistIndex())
    {
    case 0:
        int_binmin = hCorr->FindBin(-0.6);
        int_binmax = hCorr->FindBin( 0.6);
    case 1:
        int_binmin = hCorr->FindBin(-0.3);
        int_binmax = hCorr->FindBin( 0.3);
    }

    val = h->IntegralAndError(int_binmin, int_binmax, valerr, "width" );
    hY->SetBinContent(ib, val);
    hY->SetBinError(ib, valerr);
    delete h;
}
*/
enum kType { kPP, kPbPb };


class MCorr
{
    private:
        // info strings related to I/O:
        int fType;			  // 0=pp, 1=PbPb
        TString fTypeName;    // pp, PbPb
        TString fInFileName;  // name of JCORRAN input file (assuming dir: ../data/jcorran)
        TString fOutFileName; // output ROOT file name (generated in WriteOutput() )
        TFile * fInFile, * fOutFile;      //
        TString fTrackCut;    // name of track cut

        double fsumAssocBinsForMixAbove;
        double fsumTriggBinsForMixAbove;

        double fNTrigg[6][15]; // C, T
        double fNTrigg_noVCut[6][15]; // C, T (no vertex cut in phi)
        double fNEve[6]; // C

        bool fLoadedDEta, fLoadedDPhi, fLoadedDEtaDPhi, fProjectedDEtaDPhi, fFittedDEta, fFittedDPhi;
        bool fVertexCorr;

        int fWhichMixed; // choose mixed event correction: 0=hDEtaNearM, 1=3DMixed

        std::vector<TString> fFits;

    public:
        int iplot=0;
        std::vector<TH1*>    hList;
        std::vector<TString> legList;

        TString fPeriod, fAOD, fRL, fComment; // string related to input file
        // with array dimension
        const static int kT = 5, // PTt
              kA = 10,    // PTa
              kC = 5,     // Cent
              kV = 30,    // Vertex
              kF = 5,     // Fit functions
              kS = 3;     // Subtracion method

        // various cuts (-1 is all):
        double fPhiCut, fEtaCut, fVertexCut;
        // hist manager and histos from input file
        AliJHistManager * fhst;
        AliJTH1D fhDEtaNearRaw, fhDEtaNearMix, fhTriggPt, fhZVert, fhDphiAssoc;
        AliJTH1D fhChargedPt, fhiCentr, fhCentr, fhIetaTrigg, fhIetaAssoc, fhIphiTrigg, fhIphiAssoc;
        AliJTH2D fhDphiDetaPta;
        //AliJTH1D fhDEtaNear2D, fhDEtaNearFar2D, fhDPhiNear2D, fhDPhiFar2D;
        AliJBin *fCent, *fVtx, *fPTt, *fPTa, *fEta, *fPhi;

        // AliJBin indices:
        int fNumPhi, fNumEta, fNumPta, fNumPtt, fNumVtx, fNumCent;
        int fMinPTtBin, fMaxPTtBin;
        double fMinPTt, fMaxPTt;
        int fVertexSkip, fPhiSkip, fEtaSkip;


        // correction coming from different mixed event correction:
        // Filip: range/hMix->Integral(), so average is around 1
        // Ours:  1./hMix->GetBinContent() at (0,0), so hMix is efficiency correction,
        // and its maximum should be <= 1
        double fFilipCorr_eta[kC][kT][kA];

        // HISTOGRAM DEFINITIONS:

        // ntrigg
        TH1D * hTriggVSum[kC][kT];

        // fit results
        THashList * FitResultsEta;

        // histos for DEta
        TH1D * hDEtaRaw[kC][kT][kA],
             * hDEtaMix[kC][kT][kA],
             * hDEtaReal[kC][kT][kA],
             * hDEtaRealFlip[kC][kT][kA],
             * hDEtaSig[kF][kC][kT][kA],
             * hDEtaSig2D[kF][kC][kT][kA],
             * hDEtaRawVtx[kC][kV][kT][kA],
             * hDEtaMixVtx[kC][kV][kT][kA],
             * hDEtaRealVtx[kC][kV][kT][kA],
             * hDEtaMixVtxTmp[kC][kV][kT][kA],
             * hDEtaRaw2D[kC][kT][kA],
             * hDEtaMix2D[kC][kT][kA],
             * hDEtaReal2D[kC][kT][kA],
             * hDEtaReal2DFlip[kC][kT][kA],
             * hDEtaReal2Dmixed1D[kC][kT][kA],
             * hDEtaRaw2DRcut[kC][kT][kA],
             * hDEtaMix2DRcut[kC][kT][kA],
             * hDEtaReal2DRcut[kC][kT][kA],
             * hDEtaRaw2DFar[kC][kT][kA],
             * hDEtaMix2DFar[kC][kT][kA],
             * hDEtaReal2DFar[kC][kT][kA];

        // histos for DPhi
        // [Eta<1, Eta>1, subtracted], C, T, A
        TH1D * hDPhiRaw[2][kC][kT][kA],
             * hDPhiMix[2][kC][kT][kA],
             * hDPhiReal[3][kC][kT][kA],
             * hDPhiRaw2D[2][kC][kT][kA],
             * hDPhiMix2D[2][kC][kT][kA],
             * hDPhiReal2D[3][kC][kT][kA],
             * hDPhiReal2Dmixed1D[3][kC][kT][kA],
             * hDPhiRaw2DRcut[2][kC][kT][kA],
             * hDPhiMix2DRcut[2][kC][kT][kA],
             * hDPhiReal2DRcut[3][kC][kT][kA];

        TH2D * hDEtaDPhiRaw[kC][kT][kA],
             * hDEtaDPhiMix[kC][kT][kA],
             * hDEtaDPhiReal[kC][kT][kA];

        TH1D * hDEtaWingCorr[kC][kT][kA];

        TH1D * hChargedPtPub[kC]; // publised inclusive p_t

        // fit functions
        MFit * mfit_eta_1d[kF][kC][kT][kA],
             * mfit_eta_2d[kF][kC][kT][kA],
             * mfit_phi_gc[kC][kT][kA][kS],
             * mfit_phi_ggc[kC][kT][kA][kS],
             * mfit_phi_kc[kC][kT][kA][kS],
             * mfit_phi_cc[kC][kT][kA][kS],
             * mfit_phi_qgc[kC][kT][kA][kS];

        // DEFAULT CONSTRUCTOR
        MCorr() :

            fType(-1),
            fInFileName(""),
            fPeriod(""),
            fAOD(""),
            fRL(""),
            fComment(""),
            fTrackCut(""),
            fPhiCut(0.0),
            fEtaCut(0.0),
            fVertexCut(0.0),
            fVertexCorr(false),
            fMinPTt(-1),
            fMaxPTt(-1),
            fWhichMixed(0),
            fLoadedDEta(false),
            fLoadedDPhi(false),
            fLoadedDEtaDPhi(false),
            fProjectedDEtaDPhi(false),
            fFittedDEta(false),
            fFittedDPhi(false)
        {
            // default const.
        }

        // CONSTRUCTOR
        MCorr(int itype, TString inname, TString per, TString aod, TString rl, TString comment,
                TString tcut, double phicut, double etacut, double vertexcut, double minptt, double maxptt, int imix, bool vcorr) :

            fType(itype),
            fInFileName(inname),
            fPeriod(per),
            fAOD(aod),
            fRL(rl),
            fComment(comment),
            fTrackCut(tcut),
            fPhiCut(phicut),
            fEtaCut(etacut),
            fVertexCut(vertexcut),
            fVertexCorr(vcorr),
            fMinPTt(minptt),
            fMaxPTt(maxptt),
            fWhichMixed(imix),
            fLoadedDEta(false),
            fLoadedDPhi(false),
            fLoadedDEtaDPhi(false),
            fProjectedDEtaDPhi(false),
            fFittedDEta(false),
            fFittedDPhi(false),
            fsumTriggBinsForMixAbove(6.), // same as Filip
            fsumAssocBinsForMixAbove(4.) // same as Filip
        {
            switch(fType)
            {
                case kPP:	fTypeName = "pp";	  break;
                case kPbPb: fTypeName = "PbPb";	  break;
                default:	fTypeName = "NOTYPE"; break;
            }

        }

        // COPY CONSTRUCTOR
        MCorr(const MCorr& obj) :
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
            fMinPTt(obj.fMinPTt),
            fMaxPTt(obj.fMaxPTt),
            fWhichMixed(obj.fWhichMixed),
            fLoadedDEta(obj.fLoadedDEta),
            fLoadedDPhi(obj.fLoadedDPhi),
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

            //CreateOutHistos();

            std::cout <<"\nMCorr initialization done... \n";
        }

        void LoadInputFile()
        {
            std::cout << "LoadInputFile()\n";

            TString jdir = "../data/jcorran/";
            fInFile = TFile::Open( fInFileName );
            fInFile->cd();

            fhst = new AliJHistManager( "hst" );

            if( -1 != fhst->LoadConfig( Form("JDiHadronCorr_%s/AliJHistos", fTrackCut.Data()) ) ) {
                fhst->LoadConfig( Form("JDiHadronCorr_%s/AliJHistos", fTrackCut.Data()) );
            }
            else if( -1 != fhst->LoadConfig( Form("JDiHadronIAA_%s/AliJHistos", fTrackCut.Data()) ) ) {
                fhst->LoadConfig( Form("JDiHadronIAA_%s/AliJHistos", fTrackCut.Data()) );
            }

            fhDEtaNearRaw = fhst->GetTH1D("hDEtaNear");   //!        // 5D: C, V, E, T, A
            if( fWhichMixed==0 ) fhDEtaNearMix = fhst->GetTH1D("hDEtaNearM");  //!            // 5D: C, V, E, T, A
            if( fWhichMixed==1 ) fhDEtaNearMix = fhst->GetTH1D("hDEtaNearMixAcceptance");  //!        // 4D: C, V, T, A
            fhTriggPt     = fhst->GetTH1D("hTriggPtBin"); //!        // 3D: C, V, T
            fhZVert       = fhst->GetTH1D("hZVert");      //!        // 1D: C
            fhDphiAssoc   = fhst->GetTH1D("hDphiAssoc");  //!        // 5D: D, C, E, T, A
            fhDphiDetaPta = fhst->GetTH2D("hDphiDetaPta"); //!      // 4D: D, C, T, A

            fhChargedPt = fhst->GetTH1D("hChargedPt"); 	 //! // 1D: C
            fhiCentr	= fhst->GetTH1D("hiCentr");	     //! // 0D
            fhCentr	    = fhst->GetTH1D("hCentr");		 //! // 0D
            fhIetaTrigg = fhst->GetTH1D("hIetaTrigg");   //! // 2D: C, T
            fhIetaAssoc = fhst->GetTH1D("hIetaAssoc");   //! // 2D: C, T
            fhIphiTrigg = fhst->GetTH1D("fhIphiTrigg");   //! // 2D: C, T
            fhIphiAssoc = fhst->GetTH1D("fhIphiAssoc");   //! // 2D: C, T

            fCent = fhst->GetBin("Cent");   fNumCent = fCent->Size();
            fVtx  = fhst->GetBin("Vtx");    fNumVtx  = fVtx->Size();
            fPTt  = fhst->GetBin("PTt");    fNumPtt  = fPTt->Size();
            fPTa  = fhst->GetBin("PTa");    fNumPta  = fPTa->Size();
            fEta  = fhst->GetBin("EtaGap"); fNumEta  = fEta->Size();
            fPhi  = fhst->GetBin("PhiGap"); fNumPhi  = fPhi->Size();

            if( fPhiCut==-1 )    fPhiCut = fPhi->At( fNumPhi );
            if( fEtaCut==-1 )    fEtaCut = fEta->At( fNumEta );
            if( fVertexCut==-1 ) fVertexCut = fVtx->At( fNumVtx );

            fVertexSkip = fVtx->GetBin( fVertexCut );
            fPhiSkip    = fPhi->GetBin( fPhiCut );
            fEtaSkip    = fEta->GetBin( fEtaCut );
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
                        if(fLoadedDEta){
                            delete hDEtaRaw[ic][iptt][ipta];
                            delete hDEtaMix[ic][iptt][ipta];
                            delete hDEtaReal[ic][iptt][ipta];
                            delete hDEtaRealFlip[ic][iptt][ipta];
                            for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++) {
                                delete hDEtaRawVtx[ic][iv][iptt][ipta];
                                delete hDEtaMixVtx[ic][iv][iptt][ipta];
                                delete hDEtaRealVtx[ic][iv][iptt][ipta];
                            }
                        }
                        if(fProjectedDEtaDPhi) {
                            delete hDEtaRaw2DFar[ic][iptt][ipta];
                            delete hDEtaMix2DFar[ic][iptt][ipta];
                            delete hDEtaReal2DFar[ic][iptt][ipta];
                            delete hDEtaRaw2DFar[ic][iptt][ipta];
                            delete hDEtaMix2DFar[ic][iptt][ipta];
                            delete hDEtaReal2DFar[ic][iptt][ipta];
                        }
                        if(fFittedDEta)
                        {

                            //delete mfit_eta_gc[ic][iptt][ipta];
                            //delete mfit_eta_ggc[ic][iptt][ipta];
                            //delete mfit_eta_kc[ic][iptt][ipta];
                            //delete mfit_eta_cc[ic][iptt][ipta];
                            //delete mfit_eta_qgc[ic][iptt][ipta];
                        }
                        if(fLoadedDEtaDPhi)
                        {
                            delete hDEtaDPhiRaw[ic][iptt][ipta];
                            delete hDEtaDPhiMix[ic][iptt][ipta];
                            delete hDEtaDPhiReal[ic][iptt][ipta];
                        }
                        // phi related histos
                        if(fLoadedDPhi){
                            for(int ie=0; ie<2; ie++)
                            {
                                delete hDPhiRaw[ie][ic][iptt][ipta];
                                delete hDPhiMix[ie][ic][iptt][ipta];
                                delete hDPhiReal[ie][ic][iptt][ipta];
                            }
                            if(fFittedDPhi)
                            {
                                for(int isub=0; isub<kS; isub++)
                                {
                                    delete mfit_phi_gc[ic][iptt][ipta][isub];
                                    delete mfit_phi_ggc[ic][iptt][ipta][isub];
                                    delete mfit_phi_kc[ic][iptt][ipta][isub];
                                    delete mfit_phi_cc[ic][iptt][ipta][isub];
                                    delete mfit_phi_qgc[ic][iptt][ipta][isub];
                                }
                            }
                        }
                    }
                }
            }
            FitResultsEta->Clear("C");

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
				std::cout << endl << ic << std::endl;
                fNEve[ic]=fhiCentr[0]->GetBinContent(ic+1);
				std::cout << fNEve[ic] << std::endl;
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                {
                    fNTrigg[ic][iptt]=0;
                    fNTrigg_noVCut[ic][iptt]=0;
                    hTriggVSum[ic][iptt]=(TH1D*)fhTriggPt[ic][0][iptt]->Clone(); hTriggVSum[ic][iptt]->Reset();
                    for(int iv=0; iv<fNumVtx; iv++){
                        fNTrigg_noVCut[ic][iptt] += fhTriggPt[ic][iv][iptt]->GetEntries(); // count ntrigg without vertex cut (phi)
                        if( iv>=fVertexSkip && iv<(fNumVtx-fVertexSkip) )
                        {
                            hTriggVSum[ic][iptt]->Add( fhTriggPt[ic][iv][iptt] );
                            fNTrigg[ic][iptt] += fhTriggPt[ic][iv][iptt]->GetEntries(); // count ntrigg with vertex cut (eta)
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
            for(int ic=0; ic<fNumCent; ic++)
            {
                MPlot * mincl = new MPlot(iplot++, "incl. #eta", "1/N_{eve.} dN/d#eta", false);
                hList = {}; legList = {};
                for( int ipt=0; ipt<3; ipt++ )
                {
                    hList.push_back(fhIetaAssoc[ic][ipt]);
                    legList.push_back(fPTa->BuildTitle(ipt));
                } for( int ipt=0; ipt<3; ipt++ ) {
                    hList.push_back(fhIetaTrigg[ic][ipt]);
                    legList.push_back(fPTt->BuildTitle(ipt));
                }
                mincl->addHList(hList, legList, "f");
                mincl->AddInfo( BuildInfo() );
                mincl->Draw();
				mincl->Save( Form("../figs/QA/ieta_%s_C%d",fTypeName.Data(), ic) );
            }
        }

        void DrawInclPt()
        {
            GetPublishedInclPt();
            std::cout << "DrawInclPt()\n";
            MTools mt;

            const double dEta = fEta->At(fNumEta) - fEta->At(0);
            TString title = "";
            for(int ic=0; ic<fNumCent; ic++)
            {
                MPlot * mptt = new MPlot(iplot++, "p_{T}", "1/N_{eve} 1/(2#pip_{T}|#Delta#eta|) dN/dp_{T} [ (GeV/c)^{-2} ]", true);
                mptt->AddInfo( BuildInfo() );
                mt.DivideWithX( fhChargedPt[ic] ); // correcting with 1/p_T
                fhChargedPt[ic]->Scale(1./2./TMath::Pi()/dEta/fNEve[ic], "width");
                hList  = { hChargedPtPub[ic], fhChargedPt[ic] };
                if(fType==kPP) {
                    title = "arxiv:1307.1093";
                    fhChargedPt[ic]->Scale(62.2); // published data is Ed^3sigma/dp, so scaling with sigma
                    //fhChargedPt[ic]->Scale(0.7); // trigger efficiency for TPCOnly
                }
                else if(fType==kPbPb) {
                    title = "arxiv:1208.2711";
                    mptt->AddInfo( BuildCentTitle(ic) );
                    //fhChargedPt[ic]->Scale(0.9); // trigger efficiency for TPCOnly
                }
                legList ={ title, "this analysis" };
                mptt->addHList(hList, legList);

                mptt->SetLimitsXY(0.6, 30, 1E-6, 1E4);
                mptt->SetRatioLimitsXY(0.6, 30., 0, 2);
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
                    hChargedPtPub[0]->Print();
                    break;
            }
            std::cout << "GetPublishedInclPt() done...\n";
        }

        void LoadDPhiHistos()
        {
            std::cout << "MCorr::LoadDPhiHistos()...\n";

            const int ipta_mixed = fPTa->GetBin(3.);
			//const double dPhi = 1.;
            double scaleMix = 1;
			//double dEtaGapWidth = fEta->At(1)-fEta->At(0);
            TString cta;
			//TH1D * htmp = nullptr;
			//TH1 * htmp_mix = nullptr;
            MTools * mt = new MTools();

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        cta = Form("C%.0fT%.0fA%.0f",fCent->At(ic),fPTt->At(iptt),fPTa->At(ipta));
                        for(int etagap=0; etagap<2; etagap++)
                        {
                            hDPhiRaw[etagap][ic][iptt][ipta] = (TH1D*)fhDphiAssoc[0][ic][0][iptt][ipta]->Clone(Form("hDPhiRaw_E0%d%s",etagap,cta.Data()));
                            hDPhiMix[etagap][ic][iptt][ipta] = (TH1D*)fhDphiAssoc[1][ic][0][iptt][ipta]->Clone(Form("hDPhiMix_E0%d%s",etagap,cta.Data()));
                            hDPhiReal[etagap][ic][iptt][ipta] = (TH1D*)hDPhiRaw[etagap][ic][iptt][ipta]->Clone(Form("hDPhiReal_E0%d%s",etagap,cta.Data()));
                            hDPhiRaw[etagap][ic][iptt][ipta]->Reset();
                            hDPhiMix[etagap][ic][iptt][ipta]->Reset();
                            hDPhiReal[etagap][ic][iptt][ipta]->Reset();
                        }
                        for(int ie=0; ie<fNumEta; ie++)
                        {
                            if(fEta->At(ie+1) <= fEtaCut)
                            {
                                hDPhiRaw[0][ic][iptt][ipta]->Add(fhDphiAssoc[0][ic][ie][iptt][ipta]);
                                hDPhiMix[0][ic][iptt][ipta]->Add(fhDphiAssoc[1][ic][ie][iptt][ipta]);
                            }
                            if(fEta->At(ie) >= fEtaCut )
                            {
                                hDPhiRaw[1][ic][iptt][ipta]->Add(fhDphiAssoc[0][ic][ie][iptt][ipta]);
                                hDPhiMix[1][ic][iptt][ipta]->Add(fhDphiAssoc[1][ic][ie][iptt][ipta]);
                            }
                        }
                        // correcting with mixed event after the summation of etagap bins
                        for(int ie=0; ie<2; ie++)
                        {
                            //scaleMix = 2./hDPhiMix[ie][ic][iptt][ipta]->Integral();
                            //int binAtZero = hDPhiMix[0][ic][iptt][ipta]->FindBin(0.);
                            //scaleMix = hDPhiMix[0][ic][iptt][ipta]->GetBinContent( binAtZero );
                            scaleMix = 1./( GetMixedNorm1D( hDPhiMix[ie][ic][iptt][ipta]) );

                            hDPhiMix[ie][ic][iptt][ipta]->Scale( scaleMix );
                            hDPhiReal[ie][ic][iptt][ipta]->Add( hDPhiRaw[ie][ic][iptt][ipta] );
                            if(ipta > ipta_mixed) { hDPhiReal[ie][ic][iptt][ipta]->Divide( hDPhiMix[ie][ic][iptt][ipta_mixed] ); }
                            else				  { hDPhiReal[ie][ic][iptt][ipta]->Divide( hDPhiMix[ie][ic][iptt][ipta] ); }

                        }


                        // amekkora range-ben veszem, azzal el kell osztani

//                        std::cout << "DEBUG_PHIBIN" << "\t"<< fEta->At(fEtaSkip) << "\t" << fEta->At(fNumEta) << std::endl;
//                        if( fEtaCut!=fNumEta ) {
//                            hDPhiRaw[1][ic][iptt][ipta]->Scale( 2*(fEta->At(fEtaSkip)/(fEta->At(fNumEta-1)-fEta->At(fEtaSkip)) ));
//                            hDPhiReal[1][ic][iptt][ipta]->Scale( 2*(fEta->At(fEtaSkip)/(fEta->At(fNumEta-1)-fEta->At(fEtaSkip)) ));
//                        }
                    }
                    for(int ipta=0; ipta<fNumPta; ipta++)
                    {
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        for(int etagap=0; etagap<2; etagap++)
                        {
                            hDPhiRaw[etagap][ic][iptt][ipta]->Scale( 1./fNTrigg_noVCut[ic][iptt], "width" );
                            hDPhiReal[etagap][ic][iptt][ipta]->Scale( 1./fNTrigg_noVCut[ic][iptt], "width" );
                            //hDPhiReal[etagap][ic][iptt][ipta]->Rebin(2); hDPhiReal[etagap][ic][iptt][ipta]->Scale(1./2.);
                        }

                        // slightly shift large etagap histo to signal's tail
                        // mt->shiftToThisTail(hDPhiReal[1][ic][iptt][ipta], hDPhiReal[0][ic][iptt][ipta], 0.5, 1.4);
                        // save etagap background as the 3rd index
                        hDPhiReal[2][ic][iptt][ipta] = (TH1D*) hDPhiReal[0][ic][iptt][ipta]->Clone(Form("hDPhiReal_E02%s",cta.Data()));
                        hDPhiReal[2][ic][iptt][ipta]->Add(hDPhiReal[1][ic][iptt][ipta], -1.);
                    }
                }
            }
            delete mt;
            fLoadedDPhi=true;
        }

        TString GetCTA(int ic, int iptt, int ipta)
        {
            return Form("C%.0fT%.0fA%.0f",fCent->At(ic),fPTt->At(iptt),fPTa->At(ipta));
        }

        TString GetCVTA(int ic, int iv, int iptt, int ipta)
        {
            return Form("C%.0fV%dT%.0fA%.0f",fCent->At(ic),(iv),fPTt->At(iptt),fPTa->At(ipta));
        }

        double GetMixedNorm1D(TH1D * hmix)
        // Average over mixed event's 2 bins around (0)
        // (instead of using value at (0) )
        {
            double scaleMix = 0;
            double eps = 0.001;
            Int_t scaleMixBins[2];
            scaleMixBins[0] = hmix->FindBin(eps);
            scaleMixBins[1] = hmix->FindBin(-1.*eps);

            for(Int_t ib=0; ib<2; ib++)
                scaleMix += hmix->GetBinContent( scaleMixBins[ib] );

            return scaleMix/2.;
        }

        double GetMixedNorm2D(TH2D * hmix)
        // Average over mixed event's 4 bins around (0,0)
        // (instead of using value at (0,0) )
        {
            double scaleMix = 0;
            double eps = 0.001;
            Int_t scaleMixBins[4];
            scaleMixBins[0] = hmix->FindBin(eps, eps);
            scaleMixBins[1] = hmix->FindBin(eps, -1.*eps);
            scaleMixBins[2] = hmix->FindBin(-1.*eps, eps);
            scaleMixBins[3] = hmix->FindBin(-1.*eps, -1.*eps);

            for(Int_t ib=0; ib<4; ib++)
                scaleMix += hmix->GetBinContent( scaleMixBins[ib] );

            return scaleMix/4.;
        }

        void LoadDEtaDPhi()
        {
            std::cout << "MCorr::LoadDEtaDPhi()...\n";

            double scaleMix = 1.;
            TString cta = "";
            TH2D * htmp_raw = nullptr;
            TH2D * htmp_mix = nullptr;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt) < fPTa->At(ipta))
                        continue; // PTa upper border should be smaller than PTt lower

                    for(int ic=0; ic<fNumCent; ic++){
                        cta = GetCTA(ic,iptt,ipta);
                        hDEtaDPhiRaw[ic][iptt][ipta]  = (TH2D*) fhDphiDetaPta[0][ic][0][iptt][ipta]->Clone(Form("hDEtaDPhiRaw_%s",cta.Data()));
                        hDEtaDPhiMix[ic][iptt][ipta]  = (TH2D*) fhDphiDetaPta[1][ic][0][iptt][ipta]->Clone(Form("hDEtaDPhiMix_%s",cta.Data()));
                        hDEtaDPhiReal[ic][iptt][ipta] = (TH2D*) fhDphiDetaPta[0][ic][0][iptt][ipta]->Clone(Form("hDEtaDPhiReal_%s",cta.Data()));
                        hDEtaDPhiRaw[ic][iptt][ipta]->Reset();
                        hDEtaDPhiMix[ic][iptt][ipta]->Reset();
                        hDEtaDPhiReal[ic][iptt][ipta]->Reset();


                        for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                        {
                            // Add raw correlation
                            hDEtaDPhiRaw[ic][iptt][ipta]->Add( fhDphiDetaPta[0][ic][iv][iptt][ipta] );

                            // Add mixed (combine large pT bin)
                            if(fPTa->At(ipta)+0.1 < fsumAssocBinsForMixAbove || fPTt->At(iptt) + 0.1 < fsumTriggBinsForMixAbove )
                            {
                                htmp_mix = (TH2D*) fhDphiDetaPta[1][ic][iv][iptt][ipta]->Clone(Form("%s_tmp",fhDphiDetaPta[1][ic][iv][iptt][ipta]->GetName()));
                                hDEtaDPhiMix[ic][iptt][ipta]->Add(fhDphiDetaPta[1][ic][iv][iptt][ipta]);
                            } else {
                                htmp_mix = (TH2D*) fhDphiDetaPta[1][ic][iv][iptt][ipta]->Clone(Form("%s_tmp",fhDphiDetaPta[1][ic][iv][iptt][ipta]->GetName()));
                                htmp_mix->Reset();
                                for(int it=fMinPTtBin; it<fNumPtt; it++ ) {
                                    for(int ia=0; ia<fNumPta; ia++)
                                    {
                                        if(fPTt->At(it) < fPTa->At(ia)) continue;
                                        if(fPTa->At(ia)+0.1 < fsumAssocBinsForMixAbove || fPTt->At(it) + 0.1 < fsumTriggBinsForMixAbove ) continue;
                                        htmp_mix->Add( fhDphiDetaPta[1][ic][iv][it][ia] );
                                        hDEtaDPhiMix[ic][iptt][ipta]->Add( fhDphiDetaPta[1][ic][iv][it][ia] );
                                    }
                                }
                            }

                            // Correct with mixed in vertex bins if requested
                            if(fVertexCorr)
                            {
                                htmp_raw = (TH2D*) fhDphiDetaPta[0][ic][iv][iptt][ipta]->Clone(Form("%s_tmp",fhDphiDetaPta[0][ic][iv][iptt][ipta]->GetName()));
                                htmp_mix->Scale( 1./GetMixedNorm2D(htmp_mix) );
                                htmp_raw->Divide(htmp_mix);
                                hDEtaDPhiReal[ic][iptt][ipta]->Add( htmp_raw );
                                //hDEtaDPhiRaw[ic][iptt][ipta]->Add( fhDphiDetaPta[0][ic][iv][iptt][ipta] );
                            }
                        }
                        if(!fVertexCorr)
                        {
                            hDEtaDPhiReal[ic][iptt][ipta]->Add( hDEtaDPhiRaw[ic][iptt][ipta] );
                            // normalize mixed event around (0,0) average
                            scaleMix = 1./( GetMixedNorm2D(hDEtaDPhiMix[ic][iptt][ipta]) );
                            hDEtaDPhiMix[ic][iptt][ipta]->Scale(scaleMix);
                            hDEtaDPhiReal[ic][iptt][ipta]->Divide( hDEtaDPhiMix[ic][iptt][ipta] );
                        }
                        hDEtaDPhiReal[ic][iptt][ipta]->Scale( 1./fNTrigg[ic][iptt] );
                        //hDEtaDPhiRaw[ic][iptt][ipta]->Scale( 1., "width" );

                    }
                }
            }
            fLoadedDEtaDPhi=true;
        }


        void ProjectDEtaDPhi(Double_t Rcut)
        {
            std::cout << "MCorr::ProjectDEtaDPhi()...\n";
            TString cta;
            Int_t phi_firstbin, phi_lastbin;
            Int_t phi_firstbin_far, phi_lastbin_far;

            Int_t eta_firstbin, eta_midbin, eta_lastbin;
			//const double dEta = fEta->At(fNumEta) - fEta->At(0);
            MTools mt;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ipta=0; ipta<fNumPta; ipta++){
                    if(fPTt->At(iptt) < fPTa->At(ipta))
                        continue; // PTa upper border should be smaller than PTt lower
                    for(int ic=0; ic<fNumCent; ic++)
                    {
                        cta = GetCTA(ic,iptt,ipta);

						phi_firstbin = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->FindBin(-fPhiCut );
                        phi_lastbin  = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->FindBin(fPhiCut );
                        eta_firstbin = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetXaxis()->FindBin(-fEtaCut );
                        eta_midbin   = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetXaxis()->FindBin(fEtaCut );
                        eta_lastbin  = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetXaxis()->FindBin(fEta->At(fNumEta));

                        //cout << "\t" << fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->GetNbins() << "\t" << fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->GetBinWidth(2) << endl;
						//cout << "\t" << fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->GetBinLowEdge(phi_firstbin) << " - " << fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->GetBinUpEdge(phi_lastbin) << "\t" << fPhiCut << endl;

                        hDEtaRaw2D[ic][iptt][ipta]  = (TH1D*) hDEtaDPhiRaw[ic][iptt][ipta]->ProjectionX( Form("hDEtaRaw2D_%s",cta.Data()), phi_firstbin, phi_lastbin,"e" );
                        hDEtaMix2D[ic][iptt][ipta]  = (TH1D*) hDEtaDPhiMix[ic][iptt][ipta]->ProjectionX( Form("hDEtaMix2D_%s",cta.Data()), phi_firstbin, phi_lastbin,"e" );
                        hDEtaReal2D[ic][iptt][ipta] = (TH1D*) hDEtaDPhiReal[ic][iptt][ipta]->ProjectionX( Form("hDEtaReal2D_%s",cta.Data()), phi_firstbin, phi_lastbin,"e" );
                        //hDEtaReal2D[ic][iptt][ipta]->Scale(1./double(phi_lastbin-phi_firstbin));
                        //hDEtaReal2D[ic][iptt][ipta]->Scale(2.*fPhiCut);
						hDEtaReal2D[ic][iptt][ipta]->RebinX(2);
                        hDEtaReal2D[ic][iptt][ipta]->Scale(1., "width");
                        hDEtaReal2DFlip[ic][iptt][ipta] = (TH1D*)mt.Flip( hDEtaReal2D[ic][iptt][ipta] );


                        //hDEtaRaw2DRcut[ic][iptt][ipta]  = (TH1D*) mt->DoProjectionCircle(hDEtaDPhiRaw[ic][iptt][ipta], true, Form("hDEtaRaw2DCirc_%s",cta.Data()), phi_firstbin, phi_lastbin, Rcut, "e" );
                        //hDEtaMix2DRcut[ic][iptt][ipta]  = (TH1D*) mt->DoProjectionCircle(hDEtaDPhiMix[ic][iptt][ipta], true, Form("hDEtaMix2DCirc_%s",cta.Data()), phi_firstbin, phi_lastbin, Rcut, "e" );
                        //hDEtaReal2DRcut[ic][iptt][ipta] = (TH1D*) mt->DoProjectionCircle(hDEtaDPhiReal[ic][iptt][ipta], true, Form("hDEtaReal2DCirc_%s",cta.Data()), phi_firstbin, phi_lastbin, Rcut, "e" );
                        //hDEtaReal2DRcut[ic][iptt][ipta]->Scale(1./double(phi_lastbin-phi_firstbin));
                        //hDEtaReal2DRcut[ic][iptt][ipta]->Scale(2.*fPhiCut);


                        phi_firstbin_far = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->FindBin( 1.-fPhiCut );
                        phi_lastbin_far  = fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->FindBin( 1.+fPhiCut );
                        hDEtaRaw2DFar[ic][iptt][ipta]  = (TH1D*) hDEtaDPhiRaw[ic][iptt][ipta]->ProjectionX( Form("hDEtaRaw2DFar_%s",cta.Data()), phi_firstbin_far, phi_lastbin_far,"e" );
                        hDEtaMix2DFar[ic][iptt][ipta]  = (TH1D*) hDEtaDPhiMix[ic][iptt][ipta]->ProjectionX( Form("hDEtaMix2DFar_%s",cta.Data()), phi_firstbin_far, phi_lastbin_far,"e" );
                        hDEtaReal2DFar[ic][iptt][ipta] = (TH1D*) hDEtaDPhiReal[ic][iptt][ipta]->ProjectionX( Form("hDEtaReal2DFar_%s",cta.Data()), phi_firstbin_far, phi_lastbin_far,"e" );
                        //hDEtaMix2DFar[ic][iptt][ipta]->Scale(1./GetMixedNorm1D( hDEtaMix2DFar[ic][iptt][ipta]));
                        //hDEtaReal2DFar[ic][iptt][ipta]->Divide( hDEtaMix2DFar[ic][iptt][ipta] );
                        hDEtaReal2DFar[ic][iptt][ipta]->Scale(1./GetMixedNorm1D( hDEtaReal2DFar[ic][iptt][ipta] ));

/*

                        hDEtaReal2Dmixed1D[ic][iptt][ipta] = (TH1D*) hDEtaRaw2D[ic][iptt][ipta]->Clone();
                        hDEtaReal2Dmixed1D[ic][iptt][ipta]->Divide(hDEtaMix2D[ic][iptt][ipta]);
                        hDEtaReal2Dmixed1D[ic][iptt][ipta]->Scale(1./fNTrigg[ic][iptt]);
                        hDEtaReal2Dmixed1D[ic][iptt][ipta]->Scale(2.*fPhiCut);


                        hDPhiRaw2D[0][ic][iptt][ipta]  = (TH1D*) hDEtaDPhiRaw[ic][iptt][ipta]->ProjectionY( Form("hDPhiRaw2D_E00%s",cta.Data()), eta_firstbin, eta_midbin,"e" );
                        hDPhiMix2D[0][ic][iptt][ipta]  = (TH1D*) hDEtaDPhiMix[ic][iptt][ipta]->ProjectionY( Form("hDPhiMix2D_E00%s",cta.Data()), eta_firstbin, eta_midbin,"e" );
                        hDPhiReal2D[0][ic][iptt][ipta] = (TH1D*) hDEtaDPhiReal[ic][iptt][ipta]->ProjectionY( Form("hDPhiReal2D_E00%s",cta.Data()), eta_firstbin, eta_midbin,"e" );
                        hDPhiReal2D[0][ic][iptt][ipta]->Scale(1./double(eta_midbin-eta_firstbin));
                        hDPhiReal2D[0][ic][iptt][ipta]->Scale(2.*fEtaCut);

                        hDPhiRaw2D[1][ic][iptt][ipta]  = (TH1D*) hDEtaDPhiRaw[ic][iptt][ipta]->ProjectionY( Form("hDPhiRaw2D_E01%s",cta.Data()), eta_midbin, eta_lastbin,"e" );
                        hDPhiMix2D[1][ic][iptt][ipta]  = (TH1D*) hDEtaDPhiMix[ic][iptt][ipta]->ProjectionY( Form("hDPhiMix2D_E01%s",cta.Data()), eta_midbin, eta_lastbin,"e" );
                        hDPhiReal2D[1][ic][iptt][ipta] = (TH1D*) hDEtaDPhiReal[ic][iptt][ipta]->ProjectionY( Form("hDPhiReal2D_E01%s",cta.Data()), eta_midbin, eta_lastbin,"e" );
                        hDPhiReal2D[1][ic][iptt][ipta]->Scale(1./double(eta_lastbin-eta_midbin));
                        hDPhiReal2D[1][ic][iptt][ipta]->Scale(2.*fEtaCut);

                        hDPhiReal2D[2][ic][iptt][ipta] = (TH1D*)hDPhiReal2D[0][ic][iptt][ipta]->Clone(Form("hDPhiReal2D_E02%s", cta.Data()));
                        hDPhiReal2D[2][ic][iptt][ipta]->Add(hDPhiReal2D[1][ic][iptt][ipta], -1.);

                        hDPhiRaw2DRcut[0][ic][iptt][ipta]  = (TH1D*) mt->DoProjectionCircle(hDEtaDPhiRaw[ic][iptt][ipta], false, Form("hDPhiRaw2DCirc_%s",cta.Data()), eta_firstbin, eta_midbin, Rcut, "e" );
                        hDPhiMix2DRcut[0][ic][iptt][ipta]  = (TH1D*) mt->DoProjectionCircle(hDEtaDPhiMix[ic][iptt][ipta], false, Form("hDPhiMix2DCirc_%s",cta.Data()), eta_firstbin, eta_midbin, Rcut, "e" );
                        hDPhiReal2DRcut[0][ic][iptt][ipta] = (TH1D*) mt->DoProjectionCircle(hDEtaDPhiReal[ic][iptt][ipta], false, Form("hDPhiReal2DCirc_%s",cta.Data()), eta_firstbin, eta_midbin, Rcut, "e" );
                        hDPhiReal2DRcut[0][ic][iptt][ipta]->Scale(1./double(eta_midbin-eta_firstbin));
                        hDPhiReal2DRcut[0][ic][iptt][ipta]->Scale(2.*fEtaCut);

                        hDPhiRaw2DRcut[1][ic][iptt][ipta]  = (TH1D*) mt->DoProjectionCircle(hDEtaDPhiRaw[ic][iptt][ipta], false, Form("hDPhiRaw2DCirc_%s",cta.Data()), eta_midbin, eta_lastbin, Rcut, "e" );
                        hDPhiMix2DRcut[1][ic][iptt][ipta]  = (TH1D*) mt->DoProjectionCircle(hDEtaDPhiMix[ic][iptt][ipta], false, Form("hDPhiMix2DCirc_%s",cta.Data()), eta_midbin, eta_lastbin, Rcut, "e" );
                        hDPhiReal2DRcut[1][ic][iptt][ipta] = (TH1D*) mt->DoProjectionCircle(hDEtaDPhiReal[ic][iptt][ipta], false, Form("hDPhiReal2DCirc_%s",cta.Data()), eta_midbin, eta_lastbin, Rcut, "e" );
                        hDPhiReal2DRcut[1][ic][iptt][ipta]->Scale(1./double(eta_lastbin-eta_midbin));
                        hDPhiReal2DRcut[1][ic][iptt][ipta]->Scale(2.*fEtaCut);

                        hDPhiReal2Dmixed1D[0][ic][iptt][ipta] = (TH1D*)hDPhiRaw2D[0][ic][iptt][ipta]->Clone();
                        hDPhiReal2Dmixed1D[0][ic][iptt][ipta]->Divide(hDPhiMix2D[0][ic][iptt][ipta] );
                        */

                    }
                }
            }
            fProjectedDEtaDPhi = true;
        }

        inline
        void DrawWingCorr()
        {
            TString fname;
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        MPlot * mwcorr = new MPlot(iplot++, "#Delta#eta", "wing correction", false);
                        hList = { hDEtaReal2DFar[ic][iptt][ipta] };
                        legList={""};
                        mwcorr->addHList(hList, legList, "pe");
                        mwcorr->SetLimitsX(-1.6, 1.6);
                        //mwcorr->SetRatioLimitsXY(-1.6, 1.6, 0, 2.5);
                        mwcorr->AddInfo( BuildInfo() );
                        if(fType != kPP) mwcorr->AddInfo( BuildCentTitle(ic) );
                        mwcorr->AddInfo( BuildPTtTitle(iptt) );
                        mwcorr->AddInfo( BuildPTaTitle(ipta) );
                        mwcorr->SetLimitsX(-1.3,1.3);
                        mwcorr->Draw();
						fname = Form("../figs/Corr/wingcorr_%s_%s", fTypeName.Data(), GetCTA(ic,iptt,ipta).Data());
                        mwcorr->Save( fname );
                    }
                }
            }
        }

        inline
        // -------------------------------------
        // Plot 2D histograms
        // -------------------------------------
        void Draw2DHistos()
        {
            TCanvas * c[100];
            int iplot2d = 1000;
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        c[iplot2d] = new TCanvas(Form("canvasfor2dnear%d",iplot2d),"c",800,700);
						hDEtaDPhiReal[ic][iptt][ipta]->SetAxisRange(-1.5, 1.5,"x");
						hDEtaDPhiReal[ic][iptt][ipta]->SetAxisRange(-0.45, 0.45,"y");
						hDEtaDPhiReal[ic][iptt][ipta]->Draw("surf5");

						c[iplot2d]->SaveAs(Form("../figs/Corr/2D_Near_%s_C%dT%dA%d.pdf",fTypeName.Data(),ic,iptt,ipta));
                        ++iplot2d;

                        c[iplot2d] = new TCanvas(Form("canvasfor2dfar%d",iplot2d),"c",800,700);
						hDEtaDPhiReal[ic][iptt][ipta]->SetAxisRange(0.45, 1.4,"y");
						hDEtaDPhiReal[ic][iptt][ipta]->Draw("surf5");
						c[iplot2d]->SaveAs(Form("../figs/Corr/2D_Far_%s_C%dT%dA%d.pdf",fTypeName.Data(),ic,iptt,ipta));
                        ++iplot2d;

                    }
                }
            }
        }
        inline
        void DrawCorr2D1D(bool save, int which)
        {
            if(!fLoadedDEta && !fLoadedDPhi)
                return;

            TString xtit, ytit, fname;
            double xlimit = 0;

            MTools * mt = new MTools();

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

                                hList   = {  hDEtaSig[0][ic][iptt][ipta], hDEtaSig2D[0][ic][iptt][ipta] };
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
            delete mt;
        }

        void LoadDEtaHistos() {
            std::cout << "MCorr::LoadDEtaHistos()...\n";

            //const int iptt_mixed = fPTt->GetBin( fsumTriggBinsForMixAbove );
            //const int ipta_mixed = fPTa->GetBin( fsumAssocBinsForMixAbove );

            const double dEta = fEta->At(fNumEta);
            const int PhiSkipPi = fPhi->GetBin(fPhiCut * TMath::Pi() );
			//cout << "PHI SUM BIN = " << PhiSkipPi << "\t" << fPhi->At(fNumPhi) <<  endl;
            double scaleMix = 1.;

            TString cta, cvta;
            MTools mt;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt) < fPTa->At(ipta))
                        continue; // PTa upper border should be smaller than PTt lower
                    for(int ic=0; ic<fNumCent; ic++)
                    {
                        cta = Form("C%.0fT%.0fA%.0f",fCent->At(ic),fPTt->At(iptt),fPTa->At(ipta));
                        hDEtaRaw[ic][iptt][ipta] = (TH1D*)fhDEtaNearRaw[ic][0][0][iptt][ipta]->Clone(Form("hDEtaRaw_%s",cta.Data()) );
                        hDEtaRaw[ic][iptt][ipta]->Reset();

                        if( fWhichMixed == 0 ) {
                            hDEtaMix[ic][iptt][ipta] = (TH1D*)fhDEtaNearMix[ic][0][0][iptt][ipta]->Clone(Form("hDEtaMix_%s",cta.Data()) );
                            hDEtaMix[ic][iptt][ipta]->Reset();
                        }
                        if( fWhichMixed == 1 ) {
							hDEtaMix[ic][iptt][ipta] = (TH1D*)fhDEtaNearMix[ic][0][iptt][ipta]->Clone(Form("hDEtaMix_%s",cta.Data()) );
							hDEtaMix[ic][iptt][ipta]->Reset();
                        }


                        // loading vertex bins
                        for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                        {
                            cvta = Form("C%.0fV%dT%.0fA%.0f",fCent->At(ic),(iv),fPTt->At(iptt),fPTa->At(ipta));

                            hDEtaRawVtx[ic][iv][iptt][ipta] = (TH1D*)fhDEtaNearRaw[ic][iv][0][iptt][ipta]->Clone(Form("hDEtaRawVtx_%s",cvta.Data()) );
                            hDEtaRawVtx[ic][iv][iptt][ipta]->Reset();

							hDEtaMixVtx[ic][iv][iptt][ipta] = (TH1D*)fhDEtaNearMix[ic][iv][0][iptt][ipta]->Clone(Form("hDEtaMixVtx_%s",cvta.Data()) );
							hDEtaMixVtx[ic][iv][iptt][ipta]->Reset();
							hDEtaMixVtxTmp[ic][iv][iptt][ipta] = (TH1D*)fhDEtaNearMix[ic][iv][0][iptt][ipta]->Clone(Form("hDEtaRawVtx_%s",cvta.Data()) );
							hDEtaMixVtxTmp[ic][iv][iptt][ipta]->Reset();

							if( fWhichMixed == 0 ) {
								for(int ip=0; ip<PhiSkipPi; ip++) {
									hDEtaMixVtx[ic][iv][iptt][ipta]->Add( (TH1D*)fhDEtaNearMix[ic][iv][ip][iptt][ipta] );
									hDEtaMixVtxTmp[ic][iv][iptt][ipta]->Add( (TH1D*)fhDEtaNearMix[ic][iv][ip][iptt][ipta] );
								}
							}
							if( fWhichMixed == 1 ) {
								hDEtaMixVtx[ic][iv][iptt][ipta]->Add( (TH1D*)fhDEtaNearMix[ic][iv][iptt][ipta] );
								hDEtaMixVtxTmp[ic][iv][iptt][ipta]->Add( (TH1D*)fhDEtaNearMix[ic][iv][iptt][ipta] );
							}

							// simply adding up PhiGap bins in range
							for(int ip=0; ip<PhiSkipPi; ip++)
                            {
                                hDEtaRawVtx[ic][iv][iptt][ipta]->Add( (TH1D*)fhDEtaNearRaw[ic][iv][ip][iptt][ipta] );
                                hDEtaRaw[ic][iptt][ipta]->Add( (TH1D*)fhDEtaNearRaw[ic][iv][ip][iptt][ipta] );
                            }
                            // prepare to correct with mixed event
                            hDEtaRealVtx[ic][iv][iptt][ipta] = (TH1D*)hDEtaRawVtx[ic][iv][iptt][ipta]->Clone( "hDEtaRealVtx_"+cvta );
                            hDEtaMix[ic][iptt][ipta]->Add( hDEtaMixVtx[ic][iv][iptt][ipta] );
                        }
                    }
                }
            }
			// Combine larger pT bins
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
                            // Add mixed (combine large pT bin)
                            if(fPTa->At(ipta)+0.1 < fsumAssocBinsForMixAbove || fPTt->At(iptt) + 0.1 < fsumTriggBinsForMixAbove )
                            {
                                // this part is done already
                            } else {
                                //hDEtaMixVtx[ic][iv][iptt][ipta]->Reset();
                                for(int it=fMinPTtBin; it<fNumPtt; it++ ) {
                                    for(int ia=0; ia<fNumPta; ia++)
                                    {
                                        if(fPTt->At(it) < fPTa->At(ia)) continue;
                                        if(fPTa->At(ia)+0.1 < fsumAssocBinsForMixAbove || fPTt->At(it) + 0.1 < fsumTriggBinsForMixAbove ) continue;
                                        //cout << fPTt->At(iptt) << "\t" << fPTa->At(ipta) << "\t" << fPTt->At(it) << "\t" << fPTa->At(ia) << endl;

                                        hDEtaMixVtx[ic][iv][iptt][ipta]->Add( hDEtaMixVtxTmp[ic][iv][it][ia] );
                                    }
                                }
                            }
                            scaleMix = 1./( GetMixedNorm1D(hDEtaMixVtx[ic][iv][iptt][ipta]) );
                            hDEtaMixVtx[ic][iv][iptt][ipta]->Scale(scaleMix);
                        }

                        scaleMix = 1./( GetMixedNorm1D(hDEtaMix[ic][iptt][ipta]) );
                        fFilipCorr_eta[ic][iptt][ipta] = scaleMix * 2 * dEta * hDEtaMix[ic][iptt][ipta]->Integral();
                        hDEtaMix[ic][iptt][ipta]->Scale( scaleMix );
                    }
                }

                // finish and restart loop so one can choose higher pta bin for mixed event
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt) < fPTa->At(ipta))
                        continue; // PTa upper border should be smaller than PTt lower
                    for(int ic=0; ic<fNumCent; ic++)
                    {
                        cta = Form("C%.0fT%.0fA%.0f",fCent->At(ic),fPTt->At(iptt),fPTa->At(ipta));
                        if( fVertexCorr )
                        // Correction in vertex bins
                        {
                            hDEtaReal[ic][iptt][ipta]=(TH1D*)hDEtaRawVtx[ic][fVertexSkip][iptt][ipta]->Clone( "hDEtaReal_"+cta );
                            hDEtaReal[ic][iptt][ipta]->Reset();
                            for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                            {
                                hDEtaRealVtx[ic][iv][iptt][ipta]->Divide( hDEtaMixVtx[ic][iv][iptt][ipta] );
                                hDEtaReal[ic][iptt][ipta]->Add( hDEtaRealVtx[ic][iv][iptt][ipta] );
                            }
                        } else {
                        // Correction after vertex summation
                            hDEtaReal[ic][iptt][ipta] = (TH1D*)hDEtaRaw[ic][iptt][ipta]->Clone( "hDEtaReal_"+cta );
                            hDEtaReal[ic][iptt][ipta]->Divide( hDEtaMix[ic][iptt][ipta] );
                        }
						hDEtaReal[ic][iptt][ipta]->Rebin(2);
                        hDEtaReal[ic][iptt][ipta]->Scale( 1./fNTrigg[ic][iptt], "width" );
                        //hDEtaReal[ic][iptt][ipta]->Multiply( hDEtaReal2DFar[ic][iptt][ipta] ); // wing correction
                        //hDEtaReal[ic][iptt][ipta]->Multiply( hDEtaReal2DFar[ic][iptt][ipta]);
                        hDEtaRealFlip[ic][iptt][ipta] = (TH1D*) mt.Flip( hDEtaReal[ic][iptt][ipta] );
                    }
                }
            }
            fLoadedDEta=true;
        }

        void FitDEtaHistos( TString opt="EMRNS", double fitmax = 1.6 )
        {
			FitResultsEta = new THashList( 50000 );
            CreateFitContainers( FitResultsEta );

            std::cout << "MCorr::FitDEtaHistos\n";

            TFitResultPtr r1[5][kC][kT][kA];
            TFitResultPtr r2[5][kC][kT][kA];

            double bg_1d=0;
            double bg_2d=0;
            double bgerr_1d=0;
            double bgerr_2d=0;

            MTools * mt = new MTools();

            for(int ic=0; ic<fNumCent; ic++)
            {
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                {
                    for(int ipta=0; ipta<fNumPta; ipta++)
                    {
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        int ifit = 0;

                        // Gauss + Constant
                        mfit_eta_1d[ifit][ic][iptt][ipta]  = new MFit(kOneGenGaussConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax, true);
                        r1[0][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, opt);
                        mfit_eta_2d[ifit][ic][iptt][ipta]  = new MFit(kOneGenGaussConst,kDEta,hDEtaReal2DFlip[ic][iptt][ipta], 0, fitmax, true);
                        r2[0][ic][iptt][ipta] = hDEtaReal2DFlip[ic][iptt][ipta]->Fit( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, opt);
                        ++ifit;

                        // Generalized Gauss + Constant
                        mfit_eta_1d[ifit][ic][iptt][ipta] = new MFit(kOneGenGaussConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax, false);
                        mfit_eta_1d[ifit][ic][iptt][ipta]->ffit->SetParameters( mfit_eta_1d[0][ic][iptt][ipta]->ffit->GetParameters() ); // init from previous fit instead
                        r1[1][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, opt );
                        mfit_eta_2d[ifit][ic][iptt][ipta] = new MFit(kOneGenGaussConst,kDEta,hDEtaReal2DFlip[ic][iptt][ipta], 0, fitmax, false);
                        mfit_eta_2d[ifit][ic][iptt][ipta]->ffit->SetParameters( mfit_eta_2d[0][ic][iptt][ipta]->ffit->GetParameters() ); // init from previous fit instead
                        r2[1][ic][iptt][ipta] = hDEtaReal2DFlip[ic][iptt][ipta]->Fit( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, opt );
                        ++ifit;

                        // Kaplan + Constant
                        mfit_eta_1d[ifit][ic][iptt][ipta] = new MFit(kKaplanConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax);
                        r1[2][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, opt ); //EMRNSQ
                        mfit_eta_2d[ifit][ic][iptt][ipta] = new MFit(kKaplanConst,kDEta,hDEtaReal2DFlip[ic][iptt][ipta], 0, fitmax);
                        r2[2][ic][iptt][ipta] = hDEtaReal2DFlip[ic][iptt][ipta]->Fit( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, opt ); //EMRNSQ
                        ++ifit;

                        // Cauchy + Constant
                        mfit_eta_1d[ifit][ic][iptt][ipta]  = new MFit(kCauchyConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax);
                        r1[3][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, "RNSQ" ); // EMRNSQ
                        mfit_eta_2d[ifit][ic][iptt][ipta]  = new MFit(kCauchyConst,kDEta,hDEtaReal2DFlip[ic][iptt][ipta], 0, fitmax);
                        r2[3][ic][iptt][ipta] = hDEtaReal2DFlip[ic][iptt][ipta]->Fit( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, "RNSQ" ); // EMRNSQ
                        ++ifit;

                        // Two Gauss + Constant
                        mfit_eta_1d[ifit][ic][iptt][ipta] = new MFit(kTwoGenGaussConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax, true);
                        r1[4][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, opt ); //EMRNSQ
                        mfit_eta_2d[ifit][ic][iptt][ipta] = new MFit(kTwoGenGaussConst,kDEta,hDEtaReal2DFlip[ic][iptt][ipta], 0, fitmax, true);
                        r2[4][ic][iptt][ipta] = hDEtaReal2DFlip[ic][iptt][ipta]->Fit( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, opt ); //EMRNSQ
                        ++ifit;

                        for(int i=0; i<ifit; i++)
                        {
                            if(ic==0 && iptt==fMinPTtBin && ipta==0) fFits.push_back(mfit_eta_1d[ifit][ic][iptt][ipta]->GetName());
                        //	FillIntHistos( hYield_eta_Int[i][ic][iptt], (TH1D*) heta[ic][iptt][ipta], ipta+1, mfit_eta[i][ic][iptt][ipta] );

                            bg_1d    = mfit_eta_1d[i][ic][iptt][ipta]->ffit->GetParameter(0);
                            bgerr_1d = mfit_eta_1d[i][ic][iptt][ipta]->ffit->GetParError(0);
                            bg_2d    = mfit_eta_2d[i][ic][iptt][ipta]->ffit->GetParameter(0);
                            bgerr_2d = mfit_eta_2d[i][ic][iptt][ipta]->ffit->GetParError(0);
							hDEtaSig[i][ic][iptt][ipta]   = mt->subtractConstTH1( hDEtaRealFlip[ic][iptt][ipta], bg_1d, bgerr_1d, false );
							hDEtaSig2D[i][ic][iptt][ipta] = mt->subtractConstTH1( hDEtaReal2DFlip[ic][iptt][ipta], bg_2d, bgerr_2d, false );
                        }
                    }
                }
            }
            FillAllFitHistos(1,FitResultsEta, mfit_eta_1d, r1);
            FillAllFitHistos(2,FitResultsEta, mfit_eta_2d, r2);

            delete mt;
            fFittedDEta=true;


        }

        void CreateFitContainers(THashList * results)
        {
            std::cout << "MCorr::CreateFitContainers()...\n";

            double ptaBO[10];
            for(int i=0;i<(fNumPta+1);i++) {
                ptaBO[i] = fPTa->At(i);
            }
            TString name;

            for(int id=1; id<3; id++) // 1D, 2D
            {
                for(int ic=0; ic<fNumCent; ic++)
                {
                    for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                    {
                        for(int ifit = 0;ifit<kF; ifit++)
                        {
                            name = Form("%dD_F%d_C%dT%d",id,ifit,ic,iptt);
                            results->Add( new TH1D(Form("hYield_eta_%s",name.Data()), "", fNumPta, ptaBO) );
                            results->Add( new TH1D(Form("hYield_eta_INT_%s",name.Data()),"", fNumPta, ptaBO) );
                        }
                    }
                }
            }
        }

        void FillAllFitHistos(int id, THashList * results, MFit * feta[kF][kC][kT][kA], TFitResultPtr r[5][kC][kT][kA])
        {
            TString name = "";
            double val = 0;
            double valerr = 0;

            TF1 * fUE;
            TH1D * htmp;

            for(int ifit = 0;ifit<kF; ifit++)
            {
                for(int ic=0; ic<fNumCent; ic++)
                {
                    for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                    {
                        name = Form("%dD_F%d_C%dT%d",id,ifit,ic,iptt);
                        TH1D * hyield = (TH1D*) results->FindObject( Form("hYield_eta_%s",name.Data()) );
                        TH1D * hyield_int = (TH1D*) results->FindObject( Form("hYield_eta_INT_%s",name.Data()) );

                        for(int ipta=0; ipta<fNumPta; ipta++)
                        {
                            if(fPTt->At(iptt) < fPTa->At(ipta))
                                continue; // PTa upper border should be smaller than PTt lower

                            val =  feta[ifit][ic][iptt][ipta]->GetYield();
                            valerr = feta[ifit][ic][iptt][ipta]->GetYieldError(r[ifit][ic][iptt][ipta]);
                            hyield->SetBinContent(ipta+1, val);
                            hyield->SetBinError(ipta+1, valerr);

                            htmp = (TH1D*) hDEtaReal[ic][iptt][ipta]->Clone();
                            // subtract constant before integrating
                            fUE = feta[ifit][ic][iptt][ipta]->GetUE();
                            htmp->Add(fUE, -1.0);
                            int int_binmin = htmp->FindBin(-0.6);
                            int int_binmax = htmp->FindBin( 0.6);
                            val = htmp->IntegralAndError(int_binmin, int_binmax, valerr, "width" );
                            hyield_int->SetBinContent(ipta+1, val);
                            hyield_int->SetBinError(ipta+1, valerr);
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
		void DrawDEta(int id)
        {
            if(!fLoadedDEta)
                return;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        if(id==1) hList   = { hDEtaRealFlip[ic][iptt][ipta] };
                        if(id==2) hList   = { hDEtaReal2DFlip[ic][iptt][ipta] };
                        legList = {  Form("%dD |#Delta#phi|<%.1f", id, fPhiCut) };
						MPlot * meta = new MPlot(iplot++, "#Delta#eta", "1/N_{trigg.}dN/d#Delta#eta", false);

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
                            for(int ifit=0; ifit<5; ifit++)
                            {
                                if(id==1)
                                {
                                    mfit_eta_1d[ifit][ic][iptt][ipta]->Draw();
                                    meta->AddThisEntry( mfit_eta_1d[ifit][ic][iptt][ipta]->ffit, mfit_eta_1d[ifit][ic][iptt][ipta]->GetName(), "l");
                                    meta->AddThisEntry((TObject*)0, Chi2NDF(mfit_eta_1d[ifit][ic][iptt][ipta]->ffit),"");
                                }
                                if(id==2)
                                {
                                    mfit_eta_2d[ifit][ic][iptt][ipta]->Draw();
                                    meta->AddThisEntry( mfit_eta_2d[ifit][ic][iptt][ipta]->ffit, mfit_eta_2d[ifit][ic][iptt][ipta]->GetName(), "l");
                                    meta->AddThisEntry((TObject*)0, Chi2NDF(mfit_eta_2d[ifit][ic][iptt][ipta]->ffit),"");
                                }
                            }
                        }
						meta->Save( Form("../figs/Corr/eta_%s_%dD_%s", fTypeName.Data(), id, GetCTA(ic,iptt,ipta).Data()) );
                    }
                }
            }
        }

        void DrawRaw2D1D()
        {
            TH1D * hDEtaTmp = nullptr;
			TH1D * hDEtaMixTmp = nullptr;
			TH1D * hDEtaMix2DTmp = nullptr;
            MTools mt;
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        //cout << ic << "\t" << iptt << "\t" << ipta << endl;
                        MPlot * mdeta = new MPlot(++iplot, "#Delta#eta", "counts",true);

                        hDEtaRaw2D[ic][iptt][ipta]->Scale(1., "width");
                        hDEtaRaw[ic][iptt][ipta]->Scale(1., "width");
						//hDEtaTmp = (TH1D*)mt.RebinHistoToOther( hDEtaRaw[ic][iptt][ipta], hDEtaRaw2D[ic][iptt][ipta] );

						hList = { hDEtaRaw2D[ic][iptt][ipta],  hDEtaRaw[ic][iptt][ipta] };
                        legList={ "2D", "1D" };
                        mdeta->addHList(hList, legList);
                        mdeta->AddInfo( BuildInfo() );
                        mdeta->AddInfo( BuildCentTitle(ic) );
                        mdeta->AddInfo( BuildPTtTitle(iptt) );
                        mdeta->AddInfo( BuildPTaTitle(ipta)  );
                        mdeta->Draw();
                        mdeta->SetLimitsX(-1.6,1.6);
                        mdeta->SetRatioLimits(0, 2.);
						mdeta->Save(Form("../figs/Corr/2D1D_Raw_%s_C0%dT0%dA0%d", fTypeName.Data(), ic,iptt,ipta));

                        MPlot * mdetam = new MPlot(++iplot, "#Delta#eta", "counts",true);
						//hDEtaMixTmp = (TH1D*)mt.RebinHistoToOther( hDEtaMix[ic][iptt][ipta], hDEtaMix2D[ic][iptt][ipta] );
						hDEtaMixTmp = (TH1D*) hDEtaMix[ic][iptt][ipta]->Clone();
						hDEtaMixTmp->Scale( 1./GetMixedNorm1D(hDEtaMixTmp) );
						hDEtaMix2DTmp = (TH1D*) hDEtaMix2D[ic][iptt][ipta]->Clone();
                        hDEtaMix2DTmp->Scale( 1./GetMixedNorm1D(hDEtaMix2DTmp) );
                        hList = { hDEtaMixTmp, hDEtaMix2DTmp };
                        mdetam->addHList(hList, legList);
                        mdetam->AddInfo( BuildInfo() );
                        mdetam->AddInfo( BuildCentTitle(ic) );
                        mdetam->AddInfo( BuildPTtTitle(iptt) );
                        mdetam->AddInfo( BuildPTaTitle(ipta)  );
                        mdetam->Draw();
                        mdetam->SetLimitsX(-1.6,1.6);
                        mdetam->SetRatioLimits(0, 2.);
						mdetam->Save(Form("../figs/Corr/2D1D_Mix_%s_C0%dT0%dA0%d", fTypeName.Data(), ic,iptt,ipta));

                    }
                }
            }
        }

        void DrawDPhi(bool save)
        {
            if(!fLoadedDPhi)
                return;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
                        for(int isub=0; isub<3; isub++){

                               MPlot * mphi = new MPlot(iplot++, "#Delta#phi", "1/N_{trigg.}dN/d#Delta#phi", false);

                               if( isub<2 ) {
                                   hList   = { hDPhiReal2D[0][ic][iptt][ipta], hDPhiReal2D[1][ic][iptt][ipta] };
                                   legList = { Form("|#Delta#eta|<%.1f",fEtaCut), Form("|#Delta#eta|#geq%.1f",fEtaCut) };
                               }
                               if( isub==2) {
                                   hList   = { hDPhiReal2D[2][ic][iptt][ipta] };
                                   legList = { "#eta-gap subtracted" };
                               }
                               mphi->addHList(hList, legList, "pe");
                               mphi->SetLimitsX(-0.45, 1.45);

                               mphi->AddInfo( BuildInfo() );
                               if( fType != kPP ) mphi->AddInfo( BuildCentTitle(ic) );
                               mphi->AddInfo( BuildPTtTitle(iptt) );
                               mphi->AddInfo( BuildPTaTitle(ipta) );
                               mphi->Draw();

                               if(fFittedDPhi)
                               {
                                   mphi->fPad->cd();

                                   mphi->AddThisEntry( mfit_phi_gc[ic][iptt][ipta][isub]->ffit, mfit_phi_gc[ic][iptt][ipta][isub]->GetName(), "l");
                                   mphi->AddThisEntry( mfit_phi_ggc[ic][iptt][ipta][isub]->ffit, mfit_phi_ggc[ic][iptt][ipta][isub]->GetName(), "l");
                                   mphi->AddThisEntry( mfit_phi_kc[ic][iptt][ipta][isub]->ffit, mfit_phi_kc[ic][iptt][ipta][isub]->GetName(), "l");
                                   mphi->AddThisEntry( mfit_phi_cc[ic][iptt][ipta][isub]->ffit, mfit_phi_cc[ic][iptt][ipta][isub]->GetName(), "l");
                                   mphi->AddThisEntry( mfit_phi_qgc[ic][iptt][ipta][isub]->ffit, mfit_phi_qgc[ic][iptt][ipta][isub]->GetName(), "l");
                                   mfit_phi_gc[ic][iptt][ipta][isub]->Draw();
                                   mfit_phi_ggc[ic][iptt][ipta][isub]->Draw();
                                   mfit_phi_kc[ic][iptt][ipta][isub]->Draw();
                                   mfit_phi_cc[ic][iptt][ipta][isub]->Draw();
                                   mfit_phi_qgc[ic][iptt][ipta][isub]->Draw();
                               }
							   if(save) mphi->Save(Form("../figs/Corr/phi_%s_SUB%d_%s", fTypeName.Data(), isub, GetCTA(ic,iptt,ipta).Data()));
                        }
                    }
                }
            }
        }

        void DrawFitQA(bool save)
        {
            //DrawFitWidth(save);
            DrawFitYield(save);
            //DrawFitQuality(save);
            DrawFitBackg(save);
        }
        void DrawDataFitRatios(bool save)
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
                        TH1D * htmp_ggc_ratio = (TH1D*) hDEtaReal[ic][iptt][ipta]->Clone();
                        TH1D * htmp_gc_ratio = (TH1D*) hDEtaReal[ic][iptt][ipta]->Clone();
                        // since cloning doesn't work for this constructor, copying things manually...
                        MFit * mfit_gc  = new MFit(0,0,htmp_ggc_ratio, fmin, fmax, true);
                        MFit * mfit_ggc = new MFit(0,0,htmp_ggc_ratio, fmin, fmax, true);
                        mfit_gc->ffit->SetParameters(mfit_eta_1d[0][ic][iptt][ipta]->ffit->GetParameters());
                        mfit_ggc->ffit->SetParameters(mfit_eta_1d[1][ic][iptt][ipta]->ffit->GetParameters());

                        htmp_gc_ratio->Divide((TF1*) mfit_gc->ffit);
                        htmp_ggc_ratio->Divide((TF1*) mfit_ggc->ffit);

                        MPlot * meta_r = new MPlot(iplot++, "#Delta#eta", "data/fit", false);

                        hList   = { htmp_gc_ratio, htmp_ggc_ratio };
                        legList = { "Gauss", "Generalized Gauss" };
                        meta_r->addHList(hList,legList);
                        meta_r->Draw();
						if(save) meta_r->Save( Form("../figs/Corr/eta_ratio_%s_%s", fTypeName.Data(), GetCTA(ic,iptt,ipta).Data()) );
                    }
                }
            }
            // TODO
        }


/*
        void DrawFitExponent(bool save)
        {
            if(!fFittedDEta || !fFittedDPhi)
                return;
            // -------------------------------------
            // Plot exponent(pta) of GGC of fit
            // -------------------------------------
            TLine *gauss_line = new TLine(0,2,10,2);
            MPlot * mexp = new MPlot(++iplot, "p_{T, assoc} [GeV]", "#alpha (fit)", false);

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    hList = { hExpo_eta[ic][iptt], hExpo_phi[0][ic][iptt], hExpo_phi[1][ic][iptt], hExpo_phi[2][ic][iptt] };
                    legList = { "#eta", "#phi (const.)", "#phi (flow)", "#phi (#eta-gap)" };
                    mexp->addHList(hList, legList);
                    mexp->SetLimitsY(0.5,3);
                    mexp->AddInfo( BuildCentTitle(ic) );
                    mexp->AddInfo( BuildPTtTitle(iptt) );
                    mexp->Draw(); mexp->fPad->cd(); gauss_line->Draw();

                    if(save) mexp->Save(Form("figs/Fit/Expo_%s_C0%dT0%d",fTypeName.Data(),ic,iptt));
                    mexp->resetHList();
                    mexp->resetInfo();
                }
            }
        }
*/

        // draws yield
        void DrawFitYield(bool save)
        {
            TString etacutString = Form("|#Delta#eta<%.1f", fEtaCut);
            TString phicutString = Form("|#Delta#phi<%.1f", fPhiCut);

            TString name = "";
            for(int ic=0; ic<fNumCent; ic++)
            {
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                {
                    for(int ifit=0; ifit<kF; ifit++)
                    {
                        MPlot * myp = new MPlot(++iplot, "p_{T, assoc} [GeV]", "yield", false);
                        name = Form("F%d_C%dT%d",ifit,ic,iptt);
                        hList   = { (TH1D*)FitResultsEta->FindObject(Form("hYield_eta_1D_%s",name.Data())), (TH1D*)FitResultsEta->FindObject(Form("hYield_eta_2D_%s",name.Data())) };

                        legList = { "#eta 1D"+phicutString, "#eta 3D "+phicutString };
                        myp->addHList(hList, legList, "PE");
                        //myp->SetLimitsXY(0, 10, 0, 5);
                        //myp->SetAppearance( hYield_eta[ifit][ic][iptt],1, 1, 24);
                        //myp->SetAppearance( hYield_phi[ifit][0][ic][iptt],1, 1, 24);
                        //myp->SetAppearance( hYield_phi[ifit][1][ic][iptt],1, 1, 25);
                        //myp->SetAppearance( hYield_phi[ifit][2][ic][iptt],1, 1, 26);
                        //hYield_phi[ifit][0][ic][iptt]->SetMarkerColor(kBlue+1);
                        //hYield_phi[ifit][1][ic][iptt]->SetMarkerColor(kRed+1);
                        //hYield_phi[ifit][2][ic][iptt]->SetMarkerColor(kGreen+1);
                        //hYield_phi[ifit][0][ic][iptt]->SetLineColor(kBlue+1);
                        //hYield_phi[ifit][1][ic][iptt]->SetLineColor(kRed+1);
                        //hYield_phi[ifit][2][ic][iptt]->SetLineColor(kGreen+1);
                        myp->Draw();
                        myp->AddInfo( BuildInfo() );
                        myp->AddInfo( fFits[ifit] );
                        if(fType == kPbPb) myp->AddInfo( Form("Cent: %.0f-%.0f %%",fCent->At(ic), fCent->At(ic+1)));
                        myp->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", fPTt->At(iptt), fPTt->At(iptt+1)) );
                        if(save) myp->Save(Form("figs/Fit/Yield_%s_FIT%d_T0%dC0%d", fTypeName.Data(), ifit, iptt, ic));
                    }

                    name = Form("C%dT%d",ic,iptt);
                    MPlot * myp_f = new MPlot(++iplot, "p_{T, assoc} [GeV]", "yield", false);
                    hList   = { (TH1D*)FitResultsEta->FindObject(Form("hYield_eta_F0%s",name.Data())), (TH1D*)FitResultsEta->FindObject(Form("hYield_eta_F1%s",name.Data())), (TH1D*)FitResultsEta->FindObject(Form("hYield_eta_F2%s",name.Data())) };
                    legList = { fFits[0], fFits[1], fFits[2] };
                    myp_f->addHList(hList, legList, "PE");
                    //myp_f->SetLimitsXY(0, 10, 0, 5);
                    myp_f->Draw();
					if(save) myp_f->Save(Form("../figs/Fit/Yield_%s_T0%dC0%d",fTypeName.Data(), iptt, ic));
                }
            }
        }
        // drawing constant background of the fit
        void DrawFitBackg(bool save)
        {
            double * bins = new double[10];
            for(int ipta=0; ipta<=fNumPta; ipta++) bins[ipta] = fPTa->At(ipta);

            TH1D * hback[kC][kT];

            for(int ic=0; ic<fNumCent; ic++)
            {
                MPlot * mfb = new MPlot(++iplot, "p_{T, assoc} [GeV]", "fitted background", false);
                hList.clear();
                legList.clear();
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                {
                    hback[ic][iptt] = new TH1D(Form("hback_%d%d",ic,iptt), "", fNumPta, bins);
                    hback[ic][iptt]->Print();
                    for(int ipta=0; ipta<fNumPta; ipta++) {
                        if(fPTt->At(iptt)<fPTa->At(ipta))
                            continue;
                        hback[ic][iptt]->SetBinContent(ipta+1, mfit_eta_1d[0][ic][iptt][ipta]->ffit->GetParameter(0));
                        hback[ic][iptt]->SetBinError(ipta+1, mfit_eta_1d[0][ic][iptt][ipta]->ffit->GetParError(0));
                    }
                    hList.push_back( hback[ic][iptt] );
                    legList.push_back( BuildPTtTitle(iptt) );
                }
                mfb->addHList(hList, legList, "PE", "P");
                mfb->Draw();
				if(save) mfb->Save(Form("../figs/Fit/Backg_%s_C0%d",fTypeName.Data(), ic));
            }
        }

        // chi2/ndf plots
        void DrawFitQuality() {
        /*
            for(int ic=0; ic<fNumCent; ic++) {
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++) {
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower
                        MPlot * mchi = new MPlot(++iplot, "p_{T, assoc} [GeV]", "yield", false);
                        hList   = { hYield_eta[ifit][ic][iptt], hYield_phi[ifit][0][ic][iptt], hYield_phi[ifit][1][ic][iptt], hYield_phi[ifit][2][ic][iptt] };
                        legList = {
        */
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
            return Form("p_{Ta}#in %.0f-%.0f GeV", fPTa->At(ipta), fPTa->At(ipta+1));
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
};

#endif /* MCORR_H */

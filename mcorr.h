#ifndef MCORR_H
#define MCORR_H

/* 
 * ***************************************************************************
 * This macro compares our results to published
 * http://hepdata.cedar.ac.uk/view/ins930312
 * correction in vertex bins can be chosen
 * ***************************************************************************
 */

#include <TFitResultPtr.h>
#include<iostream>
#include<fstream>

#include "AliJHistManagerROOT6.cxx"
#include "mfit.h"
#include "mtools.h"
#include "mplot.h"


// Overloaded helpers for fit histogram filling
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
    double val    = 0;
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

enum kType { kPP, kPbPb };


class MCorr
{
    private:
        // info strings related to I/O:
        int fType;            // 0=pp, 1=PbPb
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
        TH1D * hWidth_eta[kF][kC][kT], 
             * hExpo_eta[kC][kT],
             * hYield_eta[kF][kC][kT], 
             * hYield_eta_Int[kF][kC][kT], 
             * hWidth_phi[kF][kS][kC][kT],  
             * hExpo_phi[kS][kC][kT],
             * hYield_phi[kF][kS][kC][kT],
             * hYield_phi_Int[kF][kS][kC][kT]; 

        // histos for DEta
        TH1D * hDEtaRaw[kC][kT][kA], 
             * hDEtaMix[kC][kT][kA], 
             * hDEtaReal[kC][kT][kA],
             * hDEtaRealFlip[kC][kT][kA],
             * hDEtaSig[kF][kC][kT][kA],
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
        MFit * mfit_eta[kF][kC][kT][kA],
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
                case kPP:   fTypeName = "pp";     break;
                case kPbPb: fTypeName = "PbPb";   break;
                default:    fTypeName = "NOTYPE"; break;
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

        inline
        void Initialize()
        {
            LoadInputFile();

            LoadNtriggers();

            CreateOutHistos();

            std::cout <<"\nMCorr initialization done... \n";
        }

        inline
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

            fhChargedPt = fhst->GetTH1D("hChargedPt");   //! // 1D: C
            fhiCentr    = fhst->GetTH1D("hiCentr");      //! // 0D
            fhCentr     = fhst->GetTH1D("hCentr");       //! // 0D
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
                    delete hExpo_eta[ic][iptt];
                    for(int isub=0; isub<kS; isub++){
                        delete hExpo_phi[isub][ic][iptt];
                    }

                    for(int ifit=0; ifit<kF; ifit++)
                    {
                        delete hYield_eta[ifit][ic][iptt];
                        delete hYield_eta_Int[ifit][ic][iptt];
                        delete hWidth_eta[ifit][ic][iptt];

                        for(int isub=0; isub<kS; isub++){
                            delete hYield_phi[ifit][isub][ic][iptt];
                            delete hWidth_phi[ifit][isub][ic][iptt];
                        }
                    }
                }
            }
            fInFile->Close();
        }
        inline 
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
                cout << endl << ic << endl;
                fNEve[ic]=fhiCentr[0]->GetBinContent(ic+1);
                cout << fNEve[ic] << endl;
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
        inline
        void DoQA()
        {
            if(fType==kPbPb) DrawCentr();

            DrawIetaIphi();

            DrawInclPt();
        }
        inline
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
            mcentr->Save( "figs/QA/centrality" );

            //print event stat. table
            for(int ic=0; ic<fNumCent; ic++)
            {
                std::cout << fCent->BuildTitle(ic) << "\t" << fNEve[ic] << std::endl;
            }
        }

        inline
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
                mincl->Save( Form("figs/QA/ieta_%s_C%d",fTypeName.Data(), ic) );
            }
        }
        inline
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
                mptt->Save( Form("figs/QA/incl_pt_%s_C%d",fTypeName.Data(), ic));

            }
        }

        inline
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

        inline
        void CreateOutHistos() {
            double ptaBO[10], centBO[10];
            TString sub[] = {"Const", "Flow", "Etagap"};
            TString fit[] = {"Gauss", "GenGauss", "Kaplan", "Cauchy", "QGauss"};
            TString name;
            for(int ipta=0; ipta<=fNumPta; ipta++){ ptaBO[ipta] = fPTa->At(ipta); }
            for(int ic=0; ic<=fNumCent; ic++){ centBO[ic] = fCent->At(ic); }
            centBO[fNumCent]=100; // create an addition bin to store pp data

            for(int ic=0; ic<fNumCent; ic++){
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                    name = Form("C%.0fT%.0f",fCent->At(ic),fPTt->At(iptt));
                    hExpo_eta[ic][iptt] = new TH1D(Form("hExpo_eta_%s",name.Data()), "", fNumPta, ptaBO);
                    for(int isub=0; isub<kS; isub++){
                        hExpo_phi[isub][ic][iptt] = new TH1D(Form("hExpo_phi_%s_SUB%d",name.Data(),isub), "", fNumPta, ptaBO);
                    }

                    for(int ifit=0; ifit<kF; ifit++){
                        name = Form("%s_%s_C%.0fT%.0f",fit[ifit].Data(),sub[0].Data(),fCent->At(ic),fPTt->At(iptt));
                        hYield_eta[ifit][ic][iptt] = new TH1D(Form("hYield_eta_%s",name.Data()), "", fNumPta, ptaBO);
                        hYield_eta_Int[ifit][ic][iptt] = new TH1D(Form("hYield_eta_fromIntegral_%s",name.Data()), "", fNumPta, ptaBO);
                        hWidth_eta[ifit][ic][iptt] = new TH1D(Form("hWidth_eta_%s",name.Data()), "", fNumPta, ptaBO);

                        for(int isub=0; isub<kS; isub++){
                            name = Form("%s_%s_C%.0fT%.0f",fit[ifit].Data(),sub[isub].Data(),fCent->At(ic),fPTt->At(iptt));
                            hYield_phi[ifit][isub][ic][iptt] = new TH1D(Form("hYield_phi_%s",name.Data()), "", fNumPta, ptaBO);
                            hYield_phi_Int[ifit][isub][ic][iptt] = new TH1D(Form("hYield_phi_fromIntegral_%s",name.Data()), "", fNumPta, ptaBO);
                            hWidth_phi[ifit][isub][ic][iptt] = new TH1D(Form("hWidth_phi_%s",name.Data()), "", fNumPta, ptaBO);
                        }
                    }
                }
            }
        }
        inline
        void LoadDPhiHistos() 
        {
            std::cout << "MCorr::LoadDPhiHistos()...\n";

            const int ipta_mixed = fPTa->GetBin(3.);
            const double dPhi = 1.;
            double scaleMix = 1;
            double dEtaGapWidth = fEta->At(1)-fEta->At(0);
            TString cta;
            TH1D * htmp = nullptr;
            TH1 * htmp_mix = nullptr;
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
                            else                  { hDPhiReal[ie][ic][iptt][ipta]->Divide( hDPhiMix[ie][ic][iptt][ipta] ); }

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


        inline
        TString GetCTA(int ic, int iptt, int ipta)
        {
            return Form("C%.0fT%.0fA%.0f",fCent->At(ic),fPTt->At(iptt),fPTa->At(ipta));
        }
        inline
        TString GetCVTA(int ic, int iv, int iptt, int ipta)
        {
            return Form("C%.0fV%dT%.0fA%.0f",fCent->At(ic),(iv),fPTt->At(iptt),fPTa->At(ipta));
        }

        inline
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
        inline
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

        inline
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
                        hDEtaDPhiReal[ic][iptt][ipta]->Scale( 1./fNTrigg[ic][iptt], "width" );
                        //hDEtaDPhiRaw[ic][iptt][ipta]->Scale( 1., "width" );

                    }
                }
            }
            fLoadedDEtaDPhi=true;
        }
        inline
        void ProjectDEtaDPhi(Double_t Rcut)
        {

            std::cout << "MCorr::ProjectDEtaDPhi()...\n";
            TString cta;
            Int_t phi_firstbin, phi_lastbin;
            Int_t phi_firstbin_far, phi_lastbin_far;

            Int_t eta_firstbin, eta_midbin, eta_lastbin;
            const double dEta = fEta->At(fNumEta) - fEta->At(0);
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

                        cout << "\t" << fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->GetNbins() << "\t" << fhDphiDetaPta[0][ic][0][iptt][ipta]->GetYaxis()->GetBinWidth(2) << endl;
                        cout << "\t" << phi_firstbin << "\t" << phi_lastbin << "\t" << fPhiCut << endl;

                        hDEtaRaw2D[ic][iptt][ipta]  = (TH1D*) hDEtaDPhiRaw[ic][iptt][ipta]->ProjectionX( Form("hDEtaRaw2D_%s",cta.Data()), phi_firstbin, phi_lastbin,"e" );
                        hDEtaMix2D[ic][iptt][ipta]  = (TH1D*) hDEtaDPhiMix[ic][iptt][ipta]->ProjectionX( Form("hDEtaMix2D_%s",cta.Data()), phi_firstbin, phi_lastbin,"e" );
                        hDEtaReal2D[ic][iptt][ipta] = (TH1D*) hDEtaDPhiReal[ic][iptt][ipta]->ProjectionX( Form("hDEtaReal2D_%s",cta.Data()), phi_firstbin, phi_lastbin,"e" );
                        //hDEtaReal2D[ic][iptt][ipta]->Scale(1./double(phi_lastbin-phi_firstbin));
                        //hDEtaReal2D[ic][iptt][ipta]->Scale(2.*fPhiCut);
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
                        fname = Form("figs/Corr/wingcorr_%s_%s", fTypeName.Data(), GetCTA(ic,iptt,ipta).Data());
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

                        c[iplot2d]->SaveAs(Form("figs/Corr/2D_Near_%s_C%dT%dA%d.pdf",fTypeName.Data(),ic,iptt,ipta));
                        ++iplot2d;

                        c[iplot2d] = new TCanvas(Form("canvasfor2dfar%d",iplot2d),"c",800,700);
                        hDEtaDPhiReal[ic][iptt][ipta]->SetAxisRange(0.45, 1.4,"y");
                        hDEtaDPhiReal[ic][iptt][ipta]->Draw("surf5");
                        c[iplot2d]->SaveAs(Form("figs/Corr/2D_Far_%s_C%dT%dA%d.pdf",fTypeName.Data(),ic,iptt,ipta));
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

            TH1D * htmp, * htmp_2D1D, * htmp_2D;
            MTools * mt = new MTools();
            MFit * mfit_sub0;
            TString opt = "RNS";

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        switch( which )
                        {
                            case kDEta:
                                htmp = mt->RebinHistoToOther( (TH1D*)hDEtaReal[ic][iptt][ipta], (TH1D*)hDEtaReal2D[ic][iptt][ipta] );
                                xtit =  "#Delta#eta"; ytit = "1/N_{trigg.}dN/d#Delta#eta";
                                fname = Form("figs/Corr/eta_2D1D_%s_%s", fTypeName.Data(), GetCTA(ic,iptt,ipta).Data());

                                mfit_sub0 = new MFit(0,0,htmp, 0, 1.4, false);
                                htmp->Fit(mfit_sub0->ffit, opt);
                                mt->subtractConstTH1(htmp, mfit_sub0->ffit->GetParameter(0));

                                htmp_2D1D = (TH1D*) hDEtaReal2Dmixed1D[ic][iptt][ipta]->Clone();
                                htmp_2D1D->Fit(mfit_sub0->ffit, opt);
                                mt->subtractConstTH1(htmp_2D1D, mfit_sub0->ffit->GetParameter(0));

                                htmp_2D = (TH1D*) hDEtaReal2D[ic][iptt][ipta]->Clone();
                                htmp_2D->Fit(mfit_sub0->ffit, opt);
                                mt->subtractConstTH1(htmp_2D, mfit_sub0->ffit->GetParameter(0));

                                hList   = {  htmp, htmp_2D1D, htmp_2D };
                                legList = {  "1D histo", "2D histo (1D corr.)", "2D histo (2D corr.)" };
                                xlimit = 1.6;
                                //delete mfit_sub0;
                                break;
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
                        }

                        MPlot * mcorr = new MPlot(iplot++, xtit, ytit, true);
                        mcorr->addHList(hList, legList, "pe");
                        mcorr->SetLimitsX(-xlimit, xlimit);
                        mcorr->SetRatioLimitsXY(-xlimit, xlimit, 0, 2.5);
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


        inline
        void LoadDEtaHistos() {
            std::cout << "MCorr::LoadDEtaHistos()...\n";

            //const int iptt_mixed = fPTt->GetBin( fsumTriggBinsForMixAbove );
            //const int ipta_mixed = fPTa->GetBin( fsumAssocBinsForMixAbove );

            const double dEta = fEta->At(fNumEta);
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
                        }


                        // loading vertex bins
                        for(int iv=fVertexSkip; iv<(fNumVtx-fVertexSkip); iv++)
                        {
                            cvta = Form("C%.0fV%dT%.0fA%.0f",fCent->At(ic),(iv),fPTt->At(iptt),fPTa->At(ipta));

                            hDEtaRawVtx[ic][iv][iptt][ipta] = (TH1D*)fhDEtaNearRaw[ic][iv][0][iptt][ipta]->Clone(Form("hDEtaRawVtx_%s",cvta.Data()) );
                            hDEtaRawVtx[ic][iv][iptt][ipta]->Reset();

                            if( fWhichMixed == 0) {
                                hDEtaMixVtx[ic][iv][iptt][ipta] = (TH1D*)fhDEtaNearMix[ic][iv][0][iptt][ipta]->Clone(Form("hDEtaMixVtx_%s",cvta.Data()) );
                                hDEtaMixVtx[ic][iv][iptt][ipta]->Reset();
                                hDEtaMixVtxTmp[ic][iv][iptt][ipta] = (TH1D*)fhDEtaNearMix[ic][iv][0][iptt][ipta]->Clone(Form("hDEtaRawVtx_%s",cvta.Data()) );
                                hDEtaMixVtxTmp[ic][iv][iptt][ipta]->Reset();
                            }
                            // simply adding up PhiGap bins in range
                            for(int ip=0; ip<fPhiSkip; ip++)
                            {
                                hDEtaRawVtx[ic][iv][iptt][ipta]->Add( (TH1D*)fhDEtaNearRaw[ic][iv][ip][iptt][ipta] );
                                hDEtaRaw[ic][iptt][ipta]->Add( (TH1D*)fhDEtaNearRaw[ic][iv][ip][iptt][ipta] );
                                if( fWhichMixed == 0 ) {
                                    hDEtaMixVtx[ic][iv][iptt][ipta]->Add( (TH1D*)fhDEtaNearMix[ic][iv][ip][iptt][ipta] );
                                    hDEtaMixVtxTmp[ic][iv][iptt][ipta]->Add( (TH1D*)fhDEtaNearMix[ic][iv][ip][iptt][ipta] );
                                }
                            }
                            // prepare to correct with mixed event
                            hDEtaRealVtx[ic][iv][iptt][ipta] = (TH1D*)hDEtaRawVtx[ic][iv][iptt][ipta]->Clone( "hDEtaRealVtx_"+cvta );
                            hDEtaMix[ic][iptt][ipta]->Add( hDEtaMixVtx[ic][iv][iptt][ipta] );
                        }
                    }
                }
            }
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

                        hDEtaReal[ic][iptt][ipta]->Scale( 1./fNTrigg[ic][iptt], "width" );
                        //hDEtaReal[ic][iptt][ipta]->Multiply( hDEtaReal2DFar[ic][iptt][ipta] ); // wing correction
                        //hDEtaReal[ic][iptt][ipta]->Multiply( hDEtaReal2DFar[ic][iptt][ipta]);
                        hDEtaRealFlip[ic][iptt][ipta] = (TH1D*) mt.Flip( hDEtaReal[ic][iptt][ipta] );
                    }
                }
            }
            fLoadedDEta=true;
        }

        inline
        void FitDEtaHistos(TString opt="EMRNS", double fitmax=1.6)
        {
            std::cout << "MCorr::FitDEtaHistos\n";

            fFits = {"Gauss", "Gen.Gauss", "Kaplan", "Cauchy", "Q-Gauss"};
            TFitResultPtr r[5][kC][kT][kA];
            double bg = 0;
            double bgerr = 0;
            MTools * mt = new MTools();

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt) < fPTa->At(ipta)) 
                        continue; // PTa upper border should be smaller than PTt lower

                    for(int ic=0; ic<fNumCent; ic++)
                    {
                        int ifit = 0;
                        // Gauss + Constant 
                        mfit_eta[ifit][ic][iptt][ipta]  = new MFit(kOneGenGaussConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax, true);
                        r[0][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta[ifit][ic][iptt][ipta]->ffit, opt);
                        FillFitHistos( (TH1D*)hYield_eta[0][ic][iptt], (TH1D*)hWidth_eta[0][ic][iptt], ipta+1, mfit_eta[ifit][ic][iptt][ipta], r[0][ic][iptt][ipta] );
                        ++ifit;
                        // Generalized Gauss + Constant 
                        mfit_eta[ifit][ic][iptt][ipta] = new MFit(kOneGenGaussConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax, false);
                        mfit_eta[ifit][ic][iptt][ipta]->ffit->SetParameters( mfit_eta[0][ic][iptt][ipta]->ffit->GetParameters() ); // init from previous fit instead
                        r[1][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta[ifit][ic][iptt][ipta]->ffit, opt );
                        FillFitHistos( (TH1D*)hYield_eta[1][ic][iptt], (TH1D*)hWidth_eta[0][ic][iptt], hExpo_eta[ic][iptt], ipta+1, mfit_eta[ifit][ic][iptt][ipta], r[1][ic][iptt][ipta] );
                        ++ifit;
                        // Kaplan + Constant
                        mfit_eta[ifit][ic][iptt][ipta] = new MFit(kKaplanConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax);
                        r[2][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta[ifit][ic][iptt][ipta]->ffit, opt ); //EMRNSQ
                        FillFitHistos( (TH1D*)hYield_eta[2][ic][iptt], ipta+1, mfit_eta[ifit][ic][iptt][ipta], r[2][ic][iptt][ipta]);
                        ++ifit;
                        // Cauchy + Constant
                        mfit_eta[ifit][ic][iptt][ipta]  = new MFit(kCauchyConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax);
                        r[3][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta[ifit][ic][iptt][ipta]->ffit, "RNSQ" ); // EMRNSQ
                        FillFitHistos( (TH1D*)hYield_eta[3][ic][iptt], ipta+1, mfit_eta[ifit][ic][iptt][ipta], r[3][ic][iptt][ipta]);
                        ++ifit;
                        // Two Gauss + Constant
                        mfit_eta[ifit][ic][iptt][ipta] = new MFit(kTwoGenGaussConst,kDEta,hDEtaRealFlip[ic][iptt][ipta], 0, fitmax, true);
                        r[4][ic][iptt][ipta] = hDEtaRealFlip[ic][iptt][ipta]->Fit( mfit_eta[ifit][ic][iptt][ipta]->ffit, opt ); //EMRNSQ
                        FillFitHistos( (TH1D*)hYield_eta[4][ic][iptt], ipta+1, mfit_eta[ifit][ic][iptt][ipta], r[4][ic][iptt][ipta] );
                        ++ifit;

                        for(int i=0; i<ifit; i++)
                        {
                            FillIntHistos( hYield_eta_Int[i][ic][iptt], (TH1D*) hDEtaRealFlip[ic][iptt][ipta], ipta+1, mfit_eta[i][ic][iptt][ipta] );

                            bg = mfit_eta[i][ic][iptt][ipta]->ffit->GetParameter(0);
                            bgerr = mfit_eta[i][ic][iptt][ipta]->ffit->GetParError(0);
                            hDEtaSig[i][ic][iptt][ipta] = mt->subtractConstTH1( hDEtaRealFlip[ic][iptt][ipta], bg, bgerr );
                        }
                    }
                }
            }
            fFittedDEta=true;
            delete mt;
        }
        inline
        void FitDPhiHistos(TString opt="EMRNS", double fitrange=0.5)
        {
            std::cout << "MCorr::FitDPhiHistos()...\n";

            double ptt=-1;
            double pta=-1;

            int fit_index = 0;
            double fitmin=-fitrange, fitmax=fitrange;

            TFitResultPtr r[5][3][kC][kT][kA];
            int subtr_i=0;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
            {
                for(int ipta=0; ipta<fNumPta; ipta++)
                {
                    if(fPTt->At(iptt) < fPTa->At(ipta)) 
                        continue; // PTa upper border should be smaller than PTt lower

                    if(fType == kPbPb)
                    {
                        ptt = (fPTt->At(iptt)+fPTt->At(iptt+1))/2.;
                        pta = (fPTa->At(ipta)+fPTt->At(ipta+1))/2.;
                    }
                    for(int ic=0; ic<fNumCent; ic++)
                    {
                        for(int isub=0; isub<3; isub++)
                        {
                            // subtraction method is kFitType + [ 0 (flat), 1 (flow), 0 (etagap) ]
                            //   => kFitType + isub%2
                            subtr_i=0;
                            if(isub==2) subtr_i=2;

                            // Gauss
                            fit_index = kOneGenGaussConst+isub%2;
                            mfit_phi_gc[ic][iptt][ipta][isub] = new MFit(fit_index,kDPhi,hDPhiReal2D[subtr_i][ic][iptt][ipta], -fitmax, fitmax, true);
                            if(isub==1) mfit_phi_gc[ic][iptt][ipta][isub]->SetPtFlow(ptt, pta);

                            r[0][isub][ic][iptt][ipta] = hDPhiReal2D[subtr_i][ic][iptt][ipta]->Fit( mfit_phi_gc[ic][iptt][ipta][isub]->ffit, opt);
                            FillFitHistos( (TH1D*)hYield_phi[0][isub][ic][iptt], (TH1D*)hWidth_phi[0][isub][ic][iptt], ipta+1, mfit_phi_gc[ic][iptt][ipta][isub], r[0][isub][ic][iptt][ipta] );

                            // Generalized Gauss
                            fit_index = kOneGenGaussConst+isub%2;
                            mfit_phi_ggc[ic][iptt][ipta][isub] = new MFit(fit_index,kDPhi,hDPhiReal2D[subtr_i][ic][iptt][ipta], -fitmax, fitmax, false);
                            mfit_phi_ggc[ic][iptt][ipta][isub]->ffit->SetParameters( mfit_phi_gc[ic][iptt][ipta][isub]->ffit->GetParameters() ); // init from previous fit instead
                            r[1][isub][ic][iptt][ipta] = hDPhiReal2D[subtr_i][ic][iptt][ipta]->Fit( mfit_phi_ggc[ic][iptt][ipta][isub]->ffit, opt );
                            FillFitHistos( (TH1D*)hYield_phi[1][isub][ic][iptt], (TH1D*)hWidth_phi[0][isub][ic][iptt],
                                    hExpo_phi[isub][ic][iptt], ipta+1, mfit_phi_ggc[ic][iptt][ipta][isub], r[1][isub][ic][iptt][ipta] );

                            // Kaplan
                            fit_index = kKaplanConst+isub%2;
                            mfit_phi_kc[ic][iptt][ipta][isub] = new MFit(fit_index,kDPhi,hDPhiReal2D[subtr_i][ic][iptt][ipta], -fitmax, fitmax);
                            if(isub==1) mfit_phi_kc[ic][iptt][ipta][isub]->SetPtFlow(ptt, pta);
                            r[2][isub][ic][iptt][ipta] = hDPhiReal2D[subtr_i][ic][iptt][ipta]->Fit( mfit_phi_kc[ic][iptt][ipta][0]->ffit, "RNSQ" );
                            FillFitHistos( (TH1D*)hYield_phi[2][isub][ic][iptt], ipta+1, mfit_phi_kc[ic][iptt][ipta][isub], r[2][isub][ic][iptt][ipta]);

                            // Cauchy
                            fit_index = kCauchyConst+isub%2;
                            mfit_phi_cc[ic][iptt][ipta][isub]  = new MFit(fit_index,kDPhi,hDPhiReal2D[subtr_i][ic][iptt][ipta], -fitmax, fitmax);
                            if(isub==1) mfit_phi_cc[ic][iptt][ipta][isub]->SetPtFlow(ptt, pta);
                            r[3][isub][ic][iptt][ipta] = hDPhiReal2D[subtr_i][ic][iptt][ipta]->Fit( mfit_phi_cc[ic][iptt][ipta][isub]->ffit, "RNSQ" ); // EMRNSQ
                            FillFitHistos( (TH1D*)hYield_phi[3][isub][ic][iptt], ipta+1, mfit_phi_cc[ic][iptt][ipta][isub], r[3][isub][ic][iptt][ipta]);

                            // Q-Gaussian
                            fit_index = kQGaussConst+isub%2;
                            mfit_phi_qgc[ic][iptt][ipta][isub] = new MFit(fit_index,kDPhi,hDPhiReal2D[subtr_i][ic][iptt][ipta], -fitmax, fitmax);
                            if(isub==1) mfit_phi_qgc[ic][iptt][ipta][isub]->SetPtFlow(ptt, pta);
                            r[4][isub][ic][iptt][ipta] = hDPhiReal2D[subtr_i][ic][iptt][ipta]->Fit( mfit_phi_qgc[ic][iptt][ipta][isub]->ffit, "RNSQ" ); //EMRNSQ
                            FillFitHistos( (TH1D*)hYield_phi[4][isub][ic][iptt], ipta+1, mfit_phi_qgc[ic][iptt][ipta][isub], r[4][isub][ic][iptt][ipta] );

                            // Filling yield with integration
                            FillIntHistos( hYield_phi_Int[0][isub][ic][iptt], (TH1D*) hDPhiReal2D[isub][ic][iptt][ipta], ipta+1, mfit_phi_gc[ic][iptt][ipta][isub] );
                            FillIntHistos( hYield_phi_Int[1][isub][ic][iptt], (TH1D*) hDPhiReal2D[isub][ic][iptt][ipta], ipta+1, mfit_phi_ggc[ic][iptt][ipta][isub] );
                            FillIntHistos( hYield_phi_Int[2][isub][ic][iptt], (TH1D*) hDPhiReal2D[isub][ic][iptt][ipta], ipta+1, mfit_phi_kc[ic][iptt][ipta][isub] );
                            FillIntHistos( hYield_phi_Int[3][isub][ic][iptt], (TH1D*) hDPhiReal2D[isub][ic][iptt][ipta], ipta+1, mfit_phi_cc[ic][iptt][ipta][isub] );
                        }
                    }
                }
            }
            fFittedDPhi=true;
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

            // write eta/phi histograms and fits
            TString cta, cvta;
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ipta=0; ipta<fNumPta; ipta++){
                    if(fPTt->At(iptt) < fPTa->At(ipta)) 
                        continue; 
                    for(int ic=0; ic<fNumCent; ic++){
                        cta = Form("C%.0fT%.0fA%.0f",fCent->At(ic),fPTt->At(iptt),fPTa->At(ipta));
                        if( fLoadedDEta ) {
                            hDEtaReal[ic][iptt][ipta]->Write();
                            hDEtaRealFlip[ic][iptt][ipta]->Write();
                            hDEtaRaw[ic][iptt][ipta]->Write();
                            hDEtaMix[ic][iptt][ipta]->Write();
                            
                            if( fFittedDEta ) {
                                //mfit_eta_gc[ic][iptt][ipta]->ffit->Write("ffit_eta_GC_"+cta);
                                //mfit_eta_ggc[ic][iptt][ipta]->ffit->Write("ffit_eta_GGC_"+cta);
                                //mfit_eta_kc[ic][iptt][ipta]->ffit->Write("ffit_eta_KC_"+cta);
                            } // FitEta
                        } // Eta

                        if( fLoadedDPhi ) {
                            for(int ietagap=0; ietagap<2; ietagap++){
                                hDPhiRaw[ietagap][ic][iptt][ipta]->Write();
                                hDPhiMix[ietagap][ic][iptt][ipta]->Write();
                                hDPhiReal[ietagap][ic][iptt][ipta]->Write();
                            }
                            if( fFittedDPhi ) {
                                for(int isub=0; isub<kS; isub++){
                                    mfit_phi_gc[ic][iptt][ipta][isub]->ffit->Write(Form("ffit_phi_GC_S%d_%s",isub,cta.Data()));
                                    mfit_phi_ggc[ic][iptt][ipta][isub]->ffit->Write(Form("ffit_phi_GGC_S%d_%s",isub,cta.Data()));
                                    mfit_phi_kc[ic][iptt][ipta][isub]->ffit->Write(Form("ffit_phi_KC_S%d_%s",isub,cta.Data()));
                                }
                            } // FitPhi
                        } // Phi
                    } // Cent
                } // PTa
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ifit=0; ifit<kF; ifit++){
                        if(fLoadedDEta && fFittedDEta)
                            hYield_eta[ifit][ic][iptt]->Write();
                        if(fLoadedDPhi && fFittedDPhi) {
                            for(int isub=0; isub<kS; isub++){
                                hYield_phi[ifit][isub][ic][iptt]->Write();
                            }
                        }
                    }
                }
            } // PTt

            fhst->WriteConfigToNew( fOutFile );
            fOutFile->Write();
            fOutFile->Close();
            delete fOutFile;

        }

        // Basic plotting will take place here
        // more elaborate plots from other macro
        void DrawDEta(bool save)
        {
            if(!fLoadedDEta)
                return;

            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    for(int ipta=0; ipta<fNumPta; ipta++){
                        if(fPTt->At(iptt) < fPTa->At(ipta))
                            continue; // PTa upper border should be smaller than PTt lower

                        hList   = {  hDEtaRealFlip[ic][iptt][ipta] };
                        legList = {  Form("|#Delta#phi|<%.1f", fPhiCut) };
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
                                mfit_eta[ifit][ic][iptt][ipta]->Draw();
                                meta->AddThisEntry( mfit_eta[ifit][ic][iptt][ipta]->ffit, mfit_eta[ifit][ic][iptt][ipta]->GetName(), "l");
                                meta->AddThisEntry((TObject*)0, Chi2NDF(mfit_eta[ifit][ic][iptt][ipta]->ffit),"");
                            }
                        }
                        if(save) meta->Save( Form("figs/Corr/eta_%s_%s", fTypeName.Data(), GetCTA(ic,iptt,ipta).Data()) );
                    }
                }
            }
        }
        inline
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
                        hDEtaTmp = (TH1D*)mt.RebinHistoToOther( hDEtaRaw[ic][iptt][ipta], hDEtaRaw2D[ic][iptt][ipta] );

                        hList = { hDEtaRaw2D[ic][iptt][ipta],  hDEtaTmp };
                        legList={ "2D", "1D" };
                        mdeta->addHList(hList, legList);
                        mdeta->AddInfo( BuildInfo() );
                        mdeta->AddInfo( BuildCentTitle(ic) );
                        mdeta->AddInfo( BuildPTtTitle(iptt) );
                        mdeta->AddInfo( BuildPTaTitle(ipta)  );
                        mdeta->Draw();
                        mdeta->SetLimitsX(-1.6,1.6);
                        mdeta->SetRatioLimits(0, 2.);
                        mdeta->Save(Form("figs/Corr/2D1D_Raw_%s_C0%dT0%dA0%d", fTypeName.Data(), ic,iptt,ipta));

                        MPlot * mdetam = new MPlot(++iplot, "#Delta#eta", "counts",true);
                        hDEtaMixTmp = (TH1D*)mt.RebinHistoToOther( hDEtaMix[ic][iptt][ipta], hDEtaMix2D[ic][iptt][ipta] );
                        //hDEtaMixTmp = (TH1D*) hDEtaMix[ic][iptt][ipta]->Clone();
                        hDEtaMix2DTmp = (TH1D*) hDEtaMix2D[ic][iptt][ipta]->Clone();
                        hDEtaMix2DTmp->Scale( 1./GetMixedNorm1D(hDEtaMix2DTmp) );
                        hDEtaMixTmp->Scale( 1./GetMixedNorm1D(hDEtaMixTmp) );
                        hList = { hDEtaMixTmp, hDEtaMix2DTmp };
                        mdetam->addHList(hList, legList);
                        mdetam->AddInfo( BuildInfo() );
                        mdetam->AddInfo( BuildCentTitle(ic) );
                        mdetam->AddInfo( BuildPTtTitle(iptt) );
                        mdetam->AddInfo( BuildPTaTitle(ipta)  );
                        mdetam->Draw();
                        mdetam->SetLimitsX(-1.6,1.6);
                        mdetam->SetRatioLimits(0, 2.);
                        mdetam->Save(Form("figs/Corr/2D1D_Mix_%s_C0%dT0%dA0%d", fTypeName.Data(), ic,iptt,ipta));

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
                               if(save) mphi->Save(Form("figs/Corr/phi_%s_SUB%d_%s", fTypeName.Data(), isub, GetCTA(ic,iptt,ipta).Data()));
                        }
                    }
                }
            }
        }

        void DrawFitQA(bool save)
        {
            //DrawDataFitRatios(save);
            DrawFitExponent(save);
            DrawFitWidth(save);
            DrawFitExponent(save);
            DrawFitYield(save);
            //DrawFitQuality(save);
            DrawFitBackg(save);
        }
        void DrawDataFitRatios(bool save)
        {
            if(!fFittedDEta)
                return;

            double fmin = mfit_eta[0][0][fMinPTtBin][0]->fitmin;
            double fmax = mfit_eta[0][0][fMinPTtBin][0]->fitmax;

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
                        mfit_gc->ffit->SetParameters(mfit_eta[0][ic][iptt][ipta]->ffit->GetParameters());
                        mfit_ggc->ffit->SetParameters(mfit_eta[1][ic][iptt][ipta]->ffit->GetParameters());

                        htmp_gc_ratio->Divide((TF1*) mfit_gc->ffit);
                        htmp_ggc_ratio->Divide((TF1*) mfit_ggc->ffit);

                        MPlot * meta_r = new MPlot(iplot++, "#Delta#eta", "data/fit", false);

                        hList   = { htmp_gc_ratio, htmp_ggc_ratio };
                        legList = { "Gauss", "Generalized Gauss" };
                        meta_r->addHList(hList,legList);
                        meta_r->Draw();
                        if(save) meta_r->Save( Form("figs/Corr/eta_ratio_%s_%s", fTypeName.Data(), GetCTA(ic,iptt,ipta).Data()) );
                    }
                }
            }
            // TODO
        }

        void DrawFitWidth(bool save)
        {
            // -------------------------------------
            // Plot width(pta) from fit
            // -------------------------------------
            for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++){
                for(int ic=0; ic<fNumCent; ic++){
                    MPlot * mwidth_e = new MPlot(++iplot, "p_{T, assoc} [GeV]", "#sigma (fit)", false);
                    hList   = { hWidth_eta[0][ic][iptt], hWidth_eta[1][ic][iptt] };
                    legList = { "Gauss", "Gen.Gauss" };
                    mwidth_e->addHList(hList, legList);
                    mwidth_e->SetLimitsY(0,1);
                    mwidth_e->AddInfo( BuildCentTitle(ic) );
                    mwidth_e->AddInfo( BuildPTtTitle(iptt) );
                    mwidth_e->Draw();
                    if(save) mwidth_e->Save(Form("figs/Fit/Width_eta_%s_C0%dT0%d",fTypeName.Data(),ic,iptt));

                    MPlot * mwidth_p = new MPlot(++iplot, "p_{T, assoc} [GeV]", "#sigma (fit)", false);
                    hList   = { hWidth_phi[0][0][ic][iptt], hWidth_phi[0][1][ic][iptt] };
                    legList = { "Gauss", "Gen.Gauss" };
                    mwidth_p->addHList(hList, legList);
                    mwidth_p->SetLimitsY(0,1);
                    mwidth_p->AddInfo( BuildCentTitle(ic) );
                    mwidth_p->AddInfo( BuildPTtTitle(iptt) );
                    mwidth_p->Draw();
                    if(save) mwidth_p->Save(Form("figs/Fit/Width_phi_%s_C0%dT0%d",fTypeName.Data(),ic,iptt));

                    MPlot * mwidth_all = new MPlot(++iplot, "p_{T, assoc} [GeV]", "#sigma (fit)", false);
                    hList   = { hWidth_eta[1][ic][iptt], hWidth_phi[1][0][ic][iptt] };
                    legList = { "#eta Gen.Gauss", "#phi Gen.Gauss" };
                    mwidth_all->addHList(hList, legList);
                    mwidth_all->SetLimitsY(0,1);
                    mwidth_all->AddInfo( BuildCentTitle(ic) );
                    mwidth_all->AddInfo( BuildPTtTitle(iptt) );
                    mwidth_all->Draw();
                    if(save) mwidth_all->Save(Form("figs/Fit/Width_all_%s_C0%dT0%d",fTypeName.Data(),ic,iptt));
                }
            }
        }

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


        // draws yield
        void DrawFitYield(bool save)
        {
            TString etacutString = Form("|#Delta#eta<%.1f", fEtaCut);
            TString phicutString = Form("|#Delta#phi<%.1f", fPhiCut);

            for(int ic=0; ic<fNumCent; ic++)
            {
                for(int iptt=fMinPTtBin; iptt<fMaxPTtBin; iptt++)
                {
                    for(int ifit=0; ifit<kF; ifit++)
                    {
                        MPlot * myp = new MPlot(++iplot, "p_{T, assoc} [GeV]", "yield", false);
                        hList   = { hYield_eta[ifit][ic][iptt], hYield_phi[ifit][0][ic][iptt], hYield_phi[ifit][1][ic][iptt], hYield_phi[ifit][2][ic][iptt] };
                        hYield_eta[ifit][ic][iptt]->Print();
                        hYield_phi[ifit][0][ic][iptt]->Print();
                        hYield_phi[ifit][1][ic][iptt]->Print();
                        hYield_phi[ifit][2][ic][iptt]->Print();
                        legList = { "#eta (flat) "+phicutString, "#phi (const.) "+etacutString, "#phi (flow) "+etacutString, "#phi (#eta-gap) "+etacutString };
                        myp->addHList(hList, legList, "PE");
                        //myp->SetLimitsXY(0, 10, 0, 5);
                        myp->SetAppearance( hYield_eta[ifit][ic][iptt],1, 1, 24);
                        myp->SetAppearance( hYield_phi[ifit][0][ic][iptt],1, 1, 24);
                        myp->SetAppearance( hYield_phi[ifit][1][ic][iptt],1, 1, 25);
                        myp->SetAppearance( hYield_phi[ifit][2][ic][iptt],1, 1, 26);
                        hYield_phi[ifit][0][ic][iptt]->SetMarkerColor(kBlue+1);
                        hYield_phi[ifit][1][ic][iptt]->SetMarkerColor(kRed+1);
                        hYield_phi[ifit][2][ic][iptt]->SetMarkerColor(kGreen+1);
                        hYield_phi[ifit][0][ic][iptt]->SetLineColor(kBlue+1);
                        hYield_phi[ifit][1][ic][iptt]->SetLineColor(kRed+1);
                        hYield_phi[ifit][2][ic][iptt]->SetLineColor(kGreen+1);
                        myp->Draw();
                        myp->AddInfo( BuildInfo() );
                        myp->AddInfo( fFits[ifit] );
                        if(fType == kPbPb) myp->AddInfo( Form("Cent: %.0f-%.0f %%",fCent->At(ic), fCent->At(ic+1)));
                        myp->AddInfo( Form("p_{Tt}#in %.0f-%.0f GeV", fPTt->At(iptt), fPTt->At(iptt+1)) );
                        if(save) myp->Save(Form("figs/Fit/Yield_%s_FIT%d_T0%dC0%d", fTypeName.Data(), ifit, iptt, ic));
                    }

                    MPlot * myp_f = new MPlot(++iplot, "p_{T, assoc} [GeV]", "yield", false);
                    hList   = { hYield_eta[0][ic][iptt], hYield_eta[1][ic][iptt], hYield_eta[2][ic][iptt], hYield_eta[3][ic][iptt], hYield_eta[4][ic][iptt] };
                    legList = { fFits[0], fFits[1], fFits[2], fFits[3], fFits[4] };
                    myp_f->addHList(hList, legList, "PE");
                    //myp_f->SetLimitsXY(0, 10, 0, 5);
                    myp_f->Draw();
                    if(save) myp_f->Save(Form("figs/Fit/Yield_%s_T0%dC0%d",fTypeName.Data(), iptt, ic));
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
                        hback[ic][iptt]->SetBinContent(ipta+1, mfit_eta[0][ic][iptt][ipta]->ffit->GetParameter(0));
                        hback[ic][iptt]->SetBinError(ipta+1, mfit_eta[0][ic][iptt][ipta]->ffit->GetParError(0));
                    }
                    hList.push_back( hback[ic][iptt] );
                    legList.push_back( BuildPTtTitle(iptt) );
                }
                mfb->addHList(hList, legList, "PE", "P");
                mfb->Draw();
                if(save) mfb->Save(Form("figs/Fit/Backg_%s_C0%d",fTypeName.Data(), ic));
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

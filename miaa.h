#include "mtools.h"
#include "mcorr.h"

#ifndef MIAA_H
#define MIAA_H

class MIaa
{

private:
    MCorr * fPP;
    MCorr * fAA;
	TString ffilenamePP;
	TString ffilenameAA;
	TFile * ffilePP;
	TFile * ffileAA;
	const TString fdatadir;

public:
    static const int kF=3;
    static const int kC=5;
    static const int kS=3;
	static const int kT=4;
    static const int kA=6;

	std::vector<double> kTBo = {3, 4, 6, 8, 15};
    std::vector<double> kABo = {0.6, 1, 2, 3, 4, 6, 8};

    TH1D * hIAA_deta_1d[kF][kC][kT][kA];
    TH1D * hIAA_deta_2d[kF][kC][kT][kA];
    TH1D * hIAA_eta_1d[kF][kC][kT];
    TH1D * hIAA_eta_2d[kF][kC][kT];
    TH1D * hIAA_eta_INT_1d[kF][kC][kT];
    TH1D * hIAA_eta_INT_2d[kF][kC][kT];

	TH1D * hDEta_1D[2][kF][kC][kT][kA];
	TH1D * hDEta_2D[2][kF][kC][kT][kA];

	MIaa(TString filenamePP, TString filenameAA):
		ffilenamePP(filenamePP),
		ffilenameAA(filenameAA),
		fdatadir("../data/syst/")
	{
		ffilePP = TFile::Open(fdatadir+ffilenamePP, "READ");
		ffileAA = TFile::Open(fdatadir+ffilenameAA, "READ");
		//CreateIAApta_file(filePP, fileAA); // not used, for this one needs to save the fit results as well
		CreateIAAdeta_file(ffilePP, ffileAA);
		//SaveIAA_file();

	}
	~MIaa()
	{
		ffilePP->Close();
		ffileAA->Close();
	}

	void CreateIAApta_file(TFile * file_pp, TFile * file_aa)
	{
		std::cout << "MIaa::CreateIAApta_file()...\n";
		for(int ic=0; ic<kC; ic++)
		{
			for(int ic=0; ic<kA; ic++)
			{
				for(int iptt=0; iptt<kT; iptt++)
				{
					for(int ifit=0; ifit<kF; ifit++)
					{
						//hIAA_deta_1d[ifit][ic][iptt] = (TH1D*) file_aa->Get("hSig_1D_F%d_PbPb_%d",GetCTA(ic,iptt,ipta).Data());
					}
				}
			}
		}
	}
	void CreateIAAdeta_file(TFile * file_pp, TFile * file_aa)
	{
		std::cout << "MIaa::CreateIAAdeta_file()...\n";

		TString name;
		for(int ic=0; ic<kC; ic++)
		{
			for(int ipta=0; ipta<kA; ipta++)
			{
				for(int iptt=0; iptt<kT; iptt++)
				{
					if(kTBo[iptt] < kABo[ipta])
						continue;
					for(int ifit=0; ifit<kF; ifit++)
					{
						name = Form("F%d_PbPb_C%dT%dA%d",ifit,ic,iptt,ipta);
						hDEta_1D[1][ifit][ic][iptt][ipta] = (TH1D*) file_aa->Get(Form("hSig_1D_%s",name.Data()))->Clone();
						name = Form("F%d_pp_C%dT%dA%d",ifit,0,iptt,ipta);
						hDEta_1D[0][ifit][ic][iptt][ipta] = (TH1D*) file_pp->Get(Form("hSig_1D_%s",name.Data()))->Clone();

						hIAA_deta_1d[ifit][ic][iptt][ipta] = (TH1D*)hDEta_1D[1][ifit][ic][iptt][ipta]->Clone();
						hIAA_deta_1d[ifit][ic][iptt][ipta]->Divide( hDEta_1D[0][ifit][ic][iptt][ipta] );

						name = Form("F%d_PbPb_C%dT%dA%d",ifit,ic,iptt,ipta);
						hDEta_2D[1][ifit][ic][iptt][ipta] = (TH1D*) file_aa->Get(Form("hSig_2D_%s",name.Data()))->Clone();
						name = Form("F%d_pp_C%dT%dA%d",ifit,0,iptt,ipta);
						hDEta_2D[0][ifit][ic][iptt][ipta] = (TH1D*) file_pp->Get(Form("hSig_2D_%s",name.Data()))->Clone();

						hIAA_deta_2d[ifit][ic][iptt][ipta] = (TH1D*)hDEta_2D[1][ifit][ic][iptt][ipta]->Clone();
						hIAA_deta_2d[ifit][ic][iptt][ipta]->Divide( hDEta_2D[0][ifit][ic][iptt][ipta] );
					}
				}
			}
		}
	}
	void SaveIAA_file()
	{
		std::cout << "MIaa::SaveIAA_file()...\n";

		TString outname = Form("%s%s__%s.root",fdatadir.Data(), ffilenamePP.Data(),ffilenameAA.Data());
		std::cout << outname << std::endl;
		TFile * foutfile = new TFile(outname,"RECREATE", "title");
		TString cta = "";

		for(int ic=0; ic<kC; ic++)
		{
			for(int ipta=0; ipta<kA; ipta++)
			{
				for(int iptt=0; iptt<kT; iptt++)
				{
					if(kTBo[iptt] < kABo[ipta])
						continue;
					for(int ifit=0; ifit<kF; ifit++)
					{
						cta = Form("C%dT%dA%d",ic,iptt,ipta);
						hIAA_deta_1d[ifit][ic][iptt][ipta]->Write(Form("IAA_1D_%s_fit%d",cta.Data(),ifit));
						hIAA_deta_2d[ifit][ic][iptt][ipta]->Write(Form("IAA_2D_%s_fit%d",cta.Data(),ifit));
					}
				}
			}
		}
		foutfile->Write();
		foutfile->Close();
	}


	MIaa(MCorr *mcorrPP, MCorr *mcorrAA):
		fPP(mcorrPP),
		fAA(mcorrAA)
	{
		cout << "MIaa::MIaa() created from MCorr objects...\n";

		CreateIAApta();
		CreateIAAdeta();
		SaveIAA();
	}
    void CreateIAApta()
    {
        // -------------------------------------
        // create IAA(pta) histos from yields 
        // extracted by MCorr for various fits 
        // and back. subtraction
        // -------------------------------------
        TString name = "";
        for(int ifit=0; ifit<kF; ifit++)
        {
            for(int ic=0; ic<fAA->fNumCent; ic++)
            {
                for(int iptt=fPP->fMinPTtBin; iptt<fPP->fNumPtt; iptt++)
                {
					name = Form("hYield_%s",fAA->GetNotUniqueFitResultKeyCent(1,ifit,ic,iptt).Data());
					hIAA_eta_1d[ifit][ic][iptt] = (TH1D*)((TH1D*)fAA->FitResultsEta->FindObject(name))->Clone();
					name = Form("hYield_%s",fPP->GetNotUniqueFitResultKeyCent(1,ifit,0,iptt).Data());
                    hIAA_eta_1d[ifit][ic][iptt]->Divide( (TH1D*)fPP->FitResultsEta->FindObject(name) );
                    //hIAA_eta_1d[ifit][ic][iptt]->Scale(fPP->GetCorrToFilip(ic,iptt,ipta)/fAA->GetCorrToFilip(ic,iptt,ipta));

					name = Form("hYield_%s",fAA->GetNotUniqueFitResultKeyCent(2,ifit,ic,iptt).Data());
                    hIAA_eta_2d[ifit][ic][iptt] = (TH1D*)((TH1D*)fAA->FitResultsEta->FindObject(name))->Clone();
					name = Form("hYield_%s",fPP->GetNotUniqueFitResultKeyCent(2,ifit,0,iptt).Data());
                    hIAA_eta_2d[ifit][ic][iptt]->Divide( (TH1D*)fPP->FitResultsEta->FindObject(name) );
                    //hIAA_eta_2d[ifit][ic][iptt]->Scale(fPP->GetCorrToFilip(ic,iptt,ipta)/fAA->GetCorrToFilip(ic,iptt,ipta));

					name = Form("hYield_INT_%s",fAA->GetNotUniqueFitResultKeyCent(1,ifit,ic,iptt).Data());
                    hIAA_eta_INT_1d[ifit][ic][iptt] = (TH1D*)((TH1D*)fAA->FitResultsEta->FindObject(name))->Clone();
					name = Form("hYield_INT_%s",fPP->GetNotUniqueFitResultKeyCent(1,ifit,0,iptt).Data());
                    hIAA_eta_INT_1d[ifit][ic][iptt]->Divide( (TH1D*)fPP->FitResultsEta->FindObject(name) );
                    //hIAA_eta_INT_1d[ifit][ic][iptt]->Scale(fPP->GetCorrToFilip(ic,iptt,ipta)/fAA->GetCorrToFilip(ic,iptt,ipta));

					name = Form("hYield_INT_%s",fAA->GetNotUniqueFitResultKeyCent(2,ifit,ic,iptt).Data());
                    hIAA_eta_INT_2d[ifit][ic][iptt] = (TH1D*)((TH1D*)fAA->FitResultsEta->FindObject(name))->Clone();
					name = Form("hYield_INT_%s",fPP->GetNotUniqueFitResultKeyCent(2,ifit,0,iptt).Data());
                    hIAA_eta_INT_2d[ifit][ic][iptt]->Divide( (TH1D*)fPP->FitResultsEta->FindObject(name) );
                    //hIAA_eta_INT_2d[ifit][ic][iptt]->Scale(fPP->GetCorrToFilip(ic,iptt,ipta)/fAA->GetCorrToFilip(ic,iptt,ipta));

                }
            }
        }
        std::cout << "MIaa::CreateIAApta() done...\n";

    }
    void CreateIAAdeta()
    {
        // -------------------------------------
        // create IAA(DEta) histos from the
        // correlation histos themselves (bin-wise)
        // -------------------------------------

        for(int ifit=0; ifit<kF; ifit++)
        {
            for(int ic=0; ic<fAA->fNumCent; ic++)
            {
                for(int iptt=fPP->fMinPTtBin; iptt<fPP->fNumPtt; iptt++)
                {
                    for(int ipta=0; ipta<fPP->fNumPta; ipta++)
                    {
                        if(fPP->fPTt->At(iptt) < fPP->fPTa->At(ipta))
                            continue;

						hIAA_deta_1d[ifit][ic][iptt][ipta] = (TH1D*) fAA->hDEtaSig1D[ifit][ic][iptt][ipta]->Clone();
						hIAA_deta_1d[ifit][ic][iptt][ipta]->Divide( fPP->hDEtaSig1D[ifit][0][iptt][ipta] );

                        hIAA_deta_2d[ifit][ic][iptt][ipta] = (TH1D*) fAA->hDEtaSig2D[ifit][ic][iptt][ipta]->Clone();
                        hIAA_deta_2d[ifit][ic][iptt][ipta]->Divide( fPP->hDEtaSig2D[ifit][0][iptt][ipta] );
                    }
                }
            }
        }
        std::cout << "MIaa::CreateIAAdeta() done...\n";
    }
	void SaveIAA()
	{
		TString outname = Form("../data/syst/%s__%s.root",fPP->GenerateOutputName().Data(),fAA->GenerateOutputName().Data());
		TFile * foutfile = new TFile(outname,"RECREATE","title");
		TString cta = "";

		for(int ifit=0; ifit<kF; ifit++)
		{
			for(int ic=0; ic<fAA->fNumCent; ic++)
			{
				for(int iptt=fPP->fMinPTtBin; iptt<fPP->fNumPtt; iptt++)
				{
					for(int ipta=0; ipta<fPP->fNumPta; ipta++)
					{
						if(fPP->fPTt->At(iptt) < fPP->fPTa->At(ipta))
							continue;

						cta = fAA->GetNotUniqueCTA(ic,iptt,ipta);
						hIAA_deta_1d[ifit][ic][iptt][ipta]->Write(Form("IAA_1D_%s_fit%d",cta.Data(),ifit));
						hIAA_deta_2d[ifit][ic][iptt][ipta]->Write(Form("IAA_2D_%s_fit%d",cta.Data(),ifit));
					}
				}
			}
		}
		foutfile->Write();
		foutfile->Close();
	}


};

#endif /* MIAA_H */

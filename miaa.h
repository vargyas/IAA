#include "mtools.h"
#include "mcorr.h"

#ifndef MIAA_H
#define MIAA_H

class MIaa
{

private:
    MCorr * fPP;
    MCorr * fAA;

public:
	static const int kF=5;
	static const int kC=5;
	static const int kS=3;
	static const int kT=5;
	static const int kA=7;

	TH1D * hIAA_deta[kF][kC][kT][kA];
	TH1D * hIAA_eta[kF][kC][kT];
	TH1D * hIAA_eta_INT[kF][kC][kT];
	//        TH1D * hIAA_phi[kF][kS][kC][kT];
	//        TH1D * hIAA_phi_INT[kF][kS][kC][kT];

    MIaa(MCorr *mcorrPP, MCorr *mcorrAA):
        fPP(mcorrPP),
        fAA(mcorrAA)
    {
		cout << "MIaa::MIaa() created...\n";
    }
        
    void CreateIAApta()
    {
        // -------------------------------------
        // create IAA(pta) histos from yields 
        // extracted by MCorr for various fits 
        // and back. subtraction
        // -------------------------------------

        for(int ifit=0; ifit<kF; ifit++)
        {
			for(int ic=0; ic<fAA->fNumCent; ic++)
            {
				for(int iptt=fPP->fMinPTtBin; iptt<fPP->fNumPtt; iptt++)
                {
					hIAA_eta[ifit][ic][iptt] = (TH1D*)fAA->hYield_eta[ifit][ic][iptt]->Clone();
					hIAA_eta[ifit][ic][iptt]->Divide( (TH1D*)fPP->hYield_eta[ifit][0][iptt] );

					hIAA_eta_INT[ifit][ic][iptt] = (TH1D*)fAA->hYield_eta_Int[ifit][ic][iptt]->Clone();
					hIAA_eta_INT[ifit][ic][iptt]->Divide( (TH1D*)fPP->hYield_eta_Int[ifit][0][iptt] );

					/*
					for(int isub=0; isub<kS; isub++)
					{
						hIAA_phi[ifit][isub][ic][iptt] = (TH1D*)fAA->hYield_phi[ifit][isub][ic][iptt]->Clone();
						hIAA_phi[ifit][isub][ic][iptt]->Divide( (TH1D*)fPP->hYield_phi[ifit][isub][0][iptt] );

						hIAA_phi_INT[ifit][isub][ic][iptt] = (TH1D*)fAA->hYield_phi_Int[ifit][isub][ic][iptt]->Clone();
						hIAA_phi_INT[ifit][isub][ic][iptt]->Divide( (TH1D*)fPP->hYield_phi_Int[ifit][isub][0][iptt] );

					}
					*/
                }
            }
        }
		cout << "MIaa::CreateIAApta() done...\n";

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

						hIAA_deta[ifit][ic][iptt][ipta] = (TH1D*) fAA->hDEtaSig[ifit][ic][iptt][ipta]->Clone();
						hIAA_deta[ifit][ic][iptt][ipta]->Divide( fPP->hDEtaSig[ifit][0][iptt][ipta] );
					}
				}
			}
		}
		cout << "MIaa::CreateIAAdeta() done...\n";
	}


};

#endif /* MIAA_H */

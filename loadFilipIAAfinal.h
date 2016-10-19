#include "mtools.h"
#include <TPad.h>

// based on files:
// ======= I_AA ======
// IaaDetaStatAndSystErr_Date17Jul2012.root
// ----------------------------------------------------------


class FilipIAA
{
    private:
        TFile * infile;
        TString inName = "../data/others/IaaDetaStatAndSystErr_Date17Jul2012.root";
        TH1D * htmp, * htmp_pp;

    public:
        TH1D * hIAAstat[5][10][10]; 
        TH1D * hIAAsyst[5][10][10]; 
        TH1D * hIAAscale[5][10][10]; 
        TGraphAsymmErrors * gIAAstat[5][10][10];
        TGraphAsymmErrors * gIAAsyst[5][10][10]; 
        TGraphAsymmErrors * gIAAscale[5][10][10]; 

        FilipIAA(){
			// Common bins:
			const double pttBorders[] = {4, 6, 8, 15};
			const double ptaBorders[] = {2, 3, 4, 6, 8};
			const double centBorders[]= {0, 10, 20, 40, 60, 90};
			const int numCent = 5;
			const int numPtt  = 3;
			const int numPta  = 4;

			// PP bins:
			const double vertBorters[] = {-10, -5, 0, 5, 10};
			const double phiBorders[]  = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
			const int numVert = 4;
			const int numPhi  = 8;


            cout << "loading Filip IAA " << endl;
            infile = TFile::Open(inName);
			//MTools * mt = new MTools();

            for(int ic=0; ic<numCent; ic++)
            {
                for(int iptt=0; iptt<numPtt; iptt++)
                {
                    for(int ipta=0; ipta<numPta; ipta++)
                    {
                        if( ptaBorders[ipta] >= pttBorders[iptt])
                            continue;
                        gIAAstat[ic][iptt][ipta] = (TGraphAsymmErrors*)infile->Get(Form("grFinalIAAStat0%d0%d0%d",ic,iptt,ipta));
                        
                        //hIAAstat[ic][iptt][ipta] = (TH1D*)(mt->TGraph2TH1D(gIAAstat[ic][iptt][ipta]))->Clone();

                        gIAAsyst[ic][iptt][ipta] = (TGraphAsymmErrors*)infile->Get(Form("grFinalIAASyst0%d0%d0%d",ic,iptt,ipta));
                        //hIAAsyst[ic][iptt][ipta] = (TH1D*)(mt->TGraph2TH1D(gIAAsyst[ic][iptt][ipta]))->Clone();

                        gIAAscale[ic][iptt][ipta] = (TGraphAsymmErrors*)infile->Get(Form("grFinalIAAScaling0%d0%d0%d",ic,iptt,ipta));
                        //hIAAscale[ic][iptt][ipta] = (TH1D*)(mt->TGraph2TH1D(gIAAscale[ic][iptt][ipta]))->Clone();
                        
                    }
                }
            }

        }
};



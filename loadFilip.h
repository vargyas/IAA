// based on files:
// ======= A+A ======
// etaOut_Cut1_WingCorr1_EffRef_Skip0000_MixLast64_Fit1Range1.6_ppzbin_PPnoSDD_CentTPC.root
// assuming others have similar arrangement

// ======= P+P ======
// pp_ALL_CUT1_pp276_withSDD_list100.root
// pp_ALL_noSDDCut1_pp276_ESD_pass2_v2.13_withoutSDD_very_long_lists.root



// Flips histogram around 0, creates a new with half the bins
TH1D *Flip(TH1D* hin){
    int nb  = hin->GetNbinsX();
    double max = hin->GetBinLowEdge(nb+1);
    TString hname = hin->GetName();
    TString newName = Form("%s_flip",hname.Data());

    TH1D *hout = new TH1D(newName.Data(), newName.Data(), (int) nb/2, 0, max);
    int zero = hin->FindBin(0.00001);
    for(int ib=zero; ib<=nb; ib++){
        double valPos = hin->GetBinContent(ib);
        double errPos = hin->GetBinError(ib);
        double valNeg = hin->GetBinContent(nb - ib+1);
        double errNeg = hin->GetBinError(nb - ib+1);

        hout->SetBinContent(ib-zero+1, (valPos+valNeg));
        hout->SetBinError(ib-zero+1, sqrt( errPos*errPos + errNeg*errNeg ) );
    }
    return hout;
}
// ----------------------------------------------------------


class FilipHistos
{
    private:
        TString fCollType;
        TFile * infile;
        TString inName;
        TH1D * htmp, * htmp_pp;

    public:
        TH1D * hIAA[5][10][10]; 
        TH1D * hDEta[5][10][10]; // flipped histos

        FilipHistos(TString type){



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



            cout << "Filip histos being initialized for [" << type << "]" << endl;
            fCollType = type;

            // processing Pb+Pb if chosen
            if(fCollType=="AA" or fCollType=="PbPb"){
                cout << "loading Filip PB+PB " << endl;
                inName = "../data/others/etaOut_Cut1_WingCorr1_EffRef_Skip0000_MixLast64_Fit1Range1.6_ppzbin_PPnoSDD_CentTPC.root";
                infile = TFile::Open(inName);

                for(int ic=0; ic<numCent; ic++){
                    for(int iptt=0; iptt<numPtt; iptt++){
                        for(int ipta=0; ipta<numPta; ipta++){
                            if( ptaBorders[ipta] > pttBorders[iptt])
                                continue;
                            hDEta[ic][iptt][ipta] = (TH1D*)infile->Get(Form("hDEtaNearSide000%d0%d0%d_flip0_sig",ic,iptt,ipta));
                            hIAA[ic][iptt][ipta] = (TH1D*)infile->Get(Form("hRatioIAA0%d0%d0%d",ic,iptt,ipta));
                        }
                    }
                }
            }

            // processing p+p if chosen
            // this needs to be summed and normalized
            if(fCollType=="PP" or fCollType=="pp") {
                cout << "loading Filip P+P " << endl;
                inName = "../data/others/pp_ALL_noSDDCut1_pp276_ESD_pass2_v2.13_withoutSDD_very_long_lists.root";
                infile = TFile::Open(inName);
                TH1D * htrigg[10];


                for(int iptt=0; iptt<numPtt; iptt++){
                    htrigg[iptt] = (TH1D*) infile->Get(Form("hTriggPtBin000%d0%d",0,iptt));
                    htrigg[iptt]->Reset();
                    for(int iv=0; iv<numVert; iv++){
                        htrigg[iptt]->Add( (TH1D*) infile->Get(Form("hTriggPtBin000%d0%d",iv,iptt)) );
                    }
                    cout << "ntrigg loaded: " << htrigg[iptt]->GetEntries() << "\t" << htrigg[iptt]->Integral() << endl;

                    for(int ipta=0; ipta<numPta; ipta++){
                        if( ptaBorders[ipta] > pttBorders[iptt])
                            continue;

                        htmp_pp = (TH1D*)infile->Get(Form("hDEtaNear00000%d0%d0%d0%d", 0, 0, iptt,ipta));
                        htmp_pp->Reset();
                        for(int iv=0; iv<numVert; iv++){
                            for(int iphi=0; iphi<numPhi; iphi++){
                                htmp = (TH1D*)infile->Get(Form("hDEtaNear00000%d0%d0%d0%d", iv, iphi, iptt,ipta));
                                htmp_pp->Add( htmp );
                            }
                        }
                        hDEta[0][iptt][ipta] = Flip( htmp_pp );
                        hDEta[0][iptt][ipta]->Scale( 1./htrigg[iptt]->Integral(), "width" );
                        hDEta[0][iptt][ipta]->Rebin(4); hDEta[0][iptt][ipta]->Scale(1./4.);
                    }
                }
            }
            if(fCollType!="PbPb" or fCollType!="AA" or fCollType!="PP" or fCollType!="pp") {
                cerr << "Please choose one from options: [PP/pp] or [AA]/PbPb" << endl;
//                return 1;
            }
        }
};



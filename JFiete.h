//#include <TFile.h>
//#include <TString.h>
//#include <TH2D.h>
//#include <TH1D.h>

const double PttBinPiete[] = { 2, 3, 4, 8, 15 };
const int NPttBinPiete = 4;
const double PtaBinPiete[] = { 0, .5, 1, 2, 3, 4, 6, 8 };
static const int NPtaBinPiete = 7;
const double CentBinPiete[8][2]={ 
    {0,10}, {60,70},{0,100}, {20,40}, {10,20},{40,60},{70,80},{80,90}};
const int NCentBinPiete = 8; // loading only central here


class JFieteCF {

    public:

        TFile * fInRoot;
        TString fInRootName;
        TH2D * fDphi2D[NPttBinPiete][NPtaBinPiete][3];
        TH1D * fDphi1D[NPttBinPiete][NPtaBinPiete][3]; // All EtaGap
        TH1D * fDeta1D[NPttBinPiete][NPtaBinPiete][3]; // All EtaGap

        enum { kPP, kPA };
        enum { kCentral, kPeripheral };

        JFieteCF(TString rootname, double phiCut, double etaCut ):
            fInRootName(rootname)
        {
            LoadHists(phiCut, etaCut);
        }

        void LoadHists( double phiCut, double etaCut )
        {
            int iphiMin, iphiMax, ietaMin, ietaMax;
            TH1D * hphi, * heta;

            fInRoot = TFile::Open(fInRootName);

            for( int ic=0;ic<NCentBinPiete;ic++ ){
                for( int iptt=0;iptt<NPttBinPiete;iptt++ ){
                    for( int ipta=1;ipta<NPtaBinPiete;ipta++ ){
                        if( PttBinPiete[iptt]<PtaBinPiete[ipta] ) continue;
                        fDphi2D[iptt][ipta][ic] = (TH2D*) fInRoot->Get(Form("dphi_%d_%d_%d",iptt, ipta, ic));
                        if( fDphi2D[iptt][ipta][ic] ){
                            // project to dphi in range -etaCut:etaCut
                            ietaMin = fDphi2D[iptt][ipta][ic]->GetYaxis()->FindBin(-etaCut); 
                            ietaMax = fDphi2D[iptt][ipta][ic]->GetYaxis()->FindBin(etaCut); 
                            hphi = fDphi2D[iptt][ipta][ic]->ProjectionX(Form("phi%d%d%d",ic,iptt,ipta),ietaMin,ietaMax,"e");

                            // project to deta in range -phiCut:phiCut
                            iphiMin = fDphi2D[iptt][ipta][ic]->GetXaxis()->FindBin(-phiCut*TMath::Pi()); 
                            iphiMax = fDphi2D[iptt][ipta][ic]->GetXaxis()->FindBin(phiCut*TMath::Pi()); 
                            heta = fDphi2D[iptt][ipta][ic]->ProjectionY(Form("eta%d%d%d",ic,iptt,ipta), iphiMin, iphiMax, "e");

                            // rebin to have our units (phi/pi)
                            TArrayD newbin(*(hphi->GetXaxis()->GetXbins()));
                            for( int i=0;i<newbin.GetSize();i++ ){
                                newbin[i] = newbin[i]/TMath::Pi();
                            }
                            hphi->SetBins( newbin.GetSize()-1, newbin.GetArray() );
                            //h->Scale(pi);
                            // store to class
                            fDphi1D[iptt][ipta][ic]=(TH1D*) hphi->Clone();
                            fDeta1D[iptt][ipta][ic]=(TH1D*) heta->Clone();
                        }
                    }
                }
            }
            cout << "Jan-Fiete correlations loaded\n";
        }

        TH1D * Find1D(TString btyp, float c0, float ptt0, float pta0){
            int iCent= -1;
            int iPtt = -1;
            int iPta = -1;
            if( btyp.Contains("pp") ){ iCent=2; }

            else if ( btyp.Contains("pA") ){
                if( c0  < 10 ) { iCent=0; }
                if( c0 > 50 ){ iCent = 2; }
            }
            cout << iCent << "\t" << iPtt << "\t" << iPta << endl;

            if( iCent < 0 ) return 0x0;
            for( int i=0;i<NPttBinPiete;i++ ){
                if( ptt0 < PttBinPiete[i+1]-1e-4 ) {
                    iPtt = i;continue;
                }
            }
            cout << iCent << "\t" << iPtt << "\t" << iPta << endl;
            if( iPtt < 0 ) return 0x0;
            for( int i=0;i<NPtaBinPiete;i++ ){
                if( pta0 < PtaBinPiete[i+1]-1e-4 ) {
                    iPta = i;break;
                }
            }
            if( iPta < 0 ) return 0x0;
            cout << iCent << "\t" << iPtt << "\t" << iPta << endl;
             fDeta1D[iPtt][iPta][iCent]->Print();
            return fDeta1D[iPtt][iPta][iCent];
        }

};


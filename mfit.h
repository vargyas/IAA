#ifndef MFIT_H
#define MFIT_H

#include "TH1.h"
#include "TF1.h"
#include <iostream>
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"

/*
 * ***************************************************************************
 * 
 * Fitting class containing function definitions and initalization methods
 *
 * ***************************************************************************
 */

// numbering convention is that [fit flow] = [fit const +1 ]
enum FitFuncType {kOneGenGaussConst, kOneGenGaussFlow, kKaplanConst, kKaplanFlow, kCauchyConst, kCauchyFlow, kQGaussConst, kQGaussFlow, kTwoGenGaussConst};
enum HistType {kDEta, kDPhi};


// -------------------------------------
// published v2{SP,Deltaeta=1.0} up to 5GeV
// http://inspirehep.net/record/900651
// -------------------------------------
/*double getPublishedV2(double ptt) {

    int p8437_d14x1y1_numpoints = 22;
    double p8437_d14x1y1_xval[] = { 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.1, 
    1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8 };
    double p8437_d14x1y1_yval[] = { 0.039107, 0.053955, 0.068788, 0.083037, 0.095749, 0.10708, 0.119547, 0.127823, 0.143002, 0.160203, 0.180128, 0.191935, 0.201014, 0.210955, 0.220726, 0.22758, 0.24464, 0.23065, 0.234024, 0.220566, 0.21736, 0.208241 };
    double p8437_d14x1y1_yerr[] = { 7.46E-4, 7.44E-4, 7.9E-4, 8.67E-4, 9.6E-4, 0.001069, 0.001187, 0.001321, 0.001091, 
    0.001339, 0.001642, 0.002016, 0.002483, 0.00305, 0.003736, 0.004558, 0.004267, 0.006108, 0.008517, 
    0.011491, 0.017248, 0.017496 };

    double pubv2 = -1.0;

    for(int ib=0; ib<p8437_d14x1y1_numpoints; ib++){
        if(p8437_d14x1y1_xval[ib+1] >= ptt) {
            if(p8437_d14x1y1_xval[ib] < ptt) {
                pubv2 = p8437_d14x1y1_yval[ib];
                break;
            }
        }
    }
    // for pt>5 use the last bin available
    if( pubv2 < 0 ) {
        pubv2 = p8437_d14x1y1_yval[p8437_d14x1y1_numpoints-1];
    }
    // option to disable flow (useful for pp)
    if (ptt==-1) { pubv2=0; } 

    return pubv2;
}
*/

// -------------------------------------
// published V2{EP,ABS(DELTA(ETA))>2.0}
// http://aliceinfo.cern.ch/ArtSubmission/node/1860
// http://hepdata.cedar.ac.uk/view/ins1116150
// -------------------------------------
double getPublishedV2(double ptt) {

    static int p8508_d1x1y1_numpoints = 28;
    double p8508_d1x1y1_xval[] = { 0.45, 0.55, 0.6499999999999999, 0.75, 0.8500000000000001, 0.95, 1.1, 1.2999999999999998, 1.5, 
        1.7000000000000002, 1.9, 2.2, 2.5999999999999996, 3.0, 3.4000000000000004, 3.8, 4.25, 4.75, 5.25, 
        5.75, 6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 18.0 };

    // Cent: 0-5%
    double p8508_d1x1y1_yval[] = { 0.0187, 0.02243, 0.02569, 0.02897, 0.03162, 0.03474, 0.03808, 0.04235, 0.04671, 
        0.05058, 0.0537, 0.05903, 0.06404, 0.06668, 0.06857, 0.06342, 0.0588, 0.05085, 0.03913, 
        0.03577, 0.02937, 0.02045, 0.02344, 0.03872, 0.01858, 0.02262, 0.01621, 0.04911 };
    double p8508_d1x1y1_yerrminus[] = { 6.166036003787198E-4, 7.368174807915458E-4, 8.868483523128404E-4, 0.0010568348972285122, 0.0010977249200050075, 0.001229186723000212, 0.0013254433220624712, 0.0014875819305167697, 0.001709502851708648, 
        0.0019226284092356484, 0.0021172859986312662, 0.0023533168082517067, 0.002617728022541685, 0.002603171143048417, 0.0034052312696790507, 0.003590055709874152, 0.0038268786236304907, 0.004266626301892397, 0.004879856555268813, 
        0.0057246921314599965, 0.005876095642516381, 0.0075600595235751945, 0.009629792313440617, 0.011963987629548938, 0.012010849262229546, 0.016623387139809983, 0.02207280000362437, 0.023145412072374087 };
    double p8508_d1x1y1_yerrplus[] = { 4.98196748283246E-4, 5.984145720150872E-4, 6.88839603971781E-4, 7.792945527847606E-4, 8.500000000000001E-4, 9.321480569094161E-4, 0.0010171037311896955, 0.00113, 0.0012529964086141667, 
        0.0013579396157414364, 0.0014454411091428113, 0.0015799050604387594, 0.0017321951391226105, 0.00185, 0.0019954197553397127, 0.0020703864373589776, 0.0022132555207205517, 0.002658307732374113, 0.0033362403990120377, 
        0.004306100323959022, 0.00421283752357007, 0.00605413082118317, 0.008253320543878081, 0.010678824841713624, 0.010381570208788266, 0.015121907948403866, 0.020664474346084878, 0.021299708918198858 };

/*
    // Cent: 5-10%
    double p8508_d1x1y2_yval[] = { 0.03128, 0.03765, 0.04316, 0.04853, 0.05358, 0.05799, 0.0643, 0.07251, 0.07989, 
        0.08708, 0.09446, 0.1035, 0.1136, 0.1208, 0.1215, 0.119, 0.1117, 0.1011, 0.09417, 
        0.08098, 0.07163, 0.05977, 0.06111, 0.05712, 0.05211, 0.02767, 0.02995, 0.03031 };
    double p8508_d1x1y2_yerrminus[] = { 9.126883367283708E-4, 0.001102905254316979, 0.0012931357237351382, 0.0015033296378372909, 0.0016137533888422977, 0.0017640861656959958, 0.0019425756098540928, 0.002193285207172109, 0.0024739644298170495, 
        0.002745268657162719, 0.003017316688715323, 0.003306055050963308, 0.003712142238654117, 0.00382099463490856, 0.004455333881989093, 0.0045891175622335065, 0.004657252408878007, 0.004810405388322277, 0.005183830629949247, 
        0.005599946428315185, 0.005747660393586246, 0.006887045810795802, 0.00835259241194014, 0.010086000198294665, 0.010418493173199279, 0.013793262848216878, 0.017975897752268172, 0.01889832003115621 };
    double p8508_d1x1y2_yerrplus[] = { 8.32946576918352E-4, 0.0010031948963187562, 0.001143547113152755, 0.0012839003076563227, 0.001424254190796011, 0.0015346986674914392, 0.0017029386365926401, 0.0019237463450257676, 0.00211463944917331, 
        0.002306274051365102, 0.0025088044961694405, 0.0027073972741361767, 0.003014962686336267, 0.0032249030993194198, 0.00327566787083184, 0.0032280024783137946, 0.0032310988842807024, 0.0031906112267087636, 0.003451405510802809, 
        0.003808162286457866, 0.0036359317925395685, 0.004816523642628571, 0.006379130034730442, 0.008210054811022886, 0.008088590730158128, 0.011573046271401492, 0.01586967548502489, 0.016129851208241196 };

    // Cent: 10-20%
    double p8508_d1x1y3_yval[] = { 0.04525, 0.05433, 0.06282, 0.07061, 0.07768, 0.08448, 0.09339, 0.1052, 0.1163, 
        0.1273, 0.1376, 0.1504, 0.164, 0.1727, 0.1751, 0.1714, 0.1595, 0.1462, 0.1302, 
        0.1138, 0.1009, 0.08808, 0.07287, 0.0784, 0.06818, 0.07142, 0.05879, 0.04567 };
    double p8508_d1x1y3_yerrminus[] = { 0.0012809761902549165, 0.0015411683879446786, 0.001800999722376436, 0.0020711832367031166, 0.0022414281161795037, 0.002451305774480205, 0.002700907254979334, 0.003001666203960727, 0.0034014702703389897, 
        0.0038013155617496425, 0.004101219330881976, 0.004501110973970759, 0.005003998401278721, 0.005108815909777921, 0.005813776741499453, 0.0058309518948453, 0.005854912467321779, 0.00582494635168428, 0.005824087911424415, 
        0.005980802621722271, 0.0061351446600711864, 0.006838318506767581, 0.0076966551176468855, 0.00887405769645431, 0.009413442515891834, 0.011713214759407427, 0.014571153008598873, 0.01587052614124686 };
    double p8508_d1x1y3_yerrplus[] = { 0.0012010412149464312, 0.0014412494579357177, 0.0016610839834276892, 0.0018713097017864252, 0.0020615528128088306, 0.0022414281161795037, 0.0024709917037497313, 0.00280178514522438, 0.0031016124838541643, 
        0.0034014702703389897, 0.0036013886210738214, 0.004001249804748511, 0.004304648650006177, 0.004609772228646444, 0.004617358552246078, 0.004539823785126467, 0.004275511665286389, 0.004080441152620633, 0.003757658845611187, 
        0.003661966684720111, 0.0033600595232822885, 0.003790646382874562, 0.00452249930901045, 0.005736497189051869, 0.005568922696536557, 0.007880012690345111, 0.010823017139411728, 0.011066349895064767 };
*/

    // it is measured up to 50%, omitting the rest now...
    //double centralities[] = {0, 5, 10, 20};

    // implemented only for most central at the moment

    double pubv2 = -1.0;

    // determining the centrality: MISSING



    // determining the trigg bin:
    for(int ib=0; ib<p8508_d1x1y1_numpoints; ib++){
        if( p8508_d1x1y1_xval[ib+1] >= ptt ) {
            if( p8508_d1x1y1_xval[ib] < ptt ) {
                pubv2 = p8508_d1x1y1_yval[ib];
                break;
            }
        }
    }
    // option to disable flow (useful for pp)
    if (ptt==-1) { pubv2=0; } 

    return pubv2;
}


/*
 * FIT FUNCTIONS as listed above
 * (flow is only v2 at the moment)
 */
double OneGenGaussConst(double *x, double *par){
    // Fit generalized Gaussian + constant [0]
    // https://en.wikipedia.org/wiki/Generalized_normal_distribution
    // (yield should be divided by 2. if flipped histograms are used)
    double deta = x[0];
    double bg    = par[0];
    double yield = par[1];
    double alpha = par[2];
    double beta  = par[3];

//    if(alpha<0.0001) alpha = 0.0001;
//    double arg = pow(fabs(deta/alpha), beta);
    double arg = TMath::Power(TMath::Abs(deta/alpha), beta);

    return bg + yield*beta/2./alpha/TMath::Gamma(1.0/beta) * exp( -1.*arg );
}

double Flow(double *x, double *par){
    double dphi    = x[0];
    double bg      = par[0];
    double v2trigg = par[1];//getPublishedV2(ptt);
    double v2assoc = par[2];//getPublishedV2(pta);

    return bg*(1+v2trigg*v2assoc*cos(2.0*dphi*TMath::Pi())); 
}

double OneGenGaussFlow(double *x, double *par){
    // Fit generalized Gaussian + flow
    // https://en.wikipedia.org/wiki/Generalized_normal_distribution
    // needs additional parameters ptt and pta to get the published v2
    double dphi = x[0];
    double bg    = par[0];
    double yield = par[1];
    double alpha = par[2];
    double beta  = par[3];
    double v2trigg = par[4];
    double v2assoc = par[5];

    if(alpha<0.0001) alpha = 0.0001;
    double arg = pow(fabs(dphi/alpha), beta);
    double flow = bg*(1+v2trigg*v2assoc*cos(2.0*dphi*TMath::Pi())); // dphi is in rad/pi units, scaling up

    return flow + yield*beta/2./alpha/TMath::Gamma(1./beta) * exp( -1.*arg );
}

double TwoGenGaussConst(double *x, double *par) {
    //Fit 2 generalized Gaussian + constant
    double deta   = x[0];
    double bg     = par[0];
    double yield1 = par[1];
    double alpha1 = par[2];
    double beta1  = par[3];
    double yield2 = par[4];
    double alpha2 = par[5];
    double beta2  = par[6];

    double arg1 = pow(fabs(deta/alpha1), beta1);
    double arg2 = pow(fabs(deta/alpha2), beta2);

    double peak1 = yield1*beta1/2./alpha1/TMath::Gamma(1./beta1) * exp( -1.*arg1 );
    double peak2 = yield2*beta2/2./alpha2/TMath::Gamma(1./beta2) * exp( -1.*arg2 );

    return bg + peak1 + peak2;
}

double TwoGenGaussFlow(double *x, double *par) {
    //Fit 2 generalized Gaussian + constant
    //used for DPhi at the moment
    double dphi   = x[0];
    double bg     = par[0];
    double v2     = par[1];
    double yield1 = par[2];
    double alpha1 = par[3];
    double beta1  = par[4];
    double yield2 = par[5];
    double alpha2 = par[6];
    double beta2  = par[7];

    double flow = bg*(1+v2*v2*cos(2*dphi*TMath::Pi()));

    double arg_near = pow(fabs((dphi-0.)/alpha1), beta1);
    double arg_far  = pow(fabs((dphi-1.)/alpha2), beta2);

    double peak_near = yield1*beta1/2./alpha1/TMath::Gamma(1./beta1) * exp( -1.*arg_near );
    double peak_far  = yield2*beta2/2./alpha2/TMath::Gamma(1./beta2) * exp( -1.*arg_far );

    return flow + peak_near + peak_far;
}

double KaplanConst(double *x, double *par){
    //Fit Kaplan as used in Filip's analysis
    double deta = x[0];
    double bg   = par[0];
    double ampl = par[1];
    double k    = par[2];
    double n    = par[3];

    return bg + ampl*pow( 1 + k*deta*deta, -n );
}

double KaplanFlow(double *x, double *par){
    //Fit Kaplan as used in Filip's analysis
    double dphi = x[0];
    double bg   = par[0];
    double ampl = par[1];
    double k    = par[2];
    double n    = par[3];
    double v2trigg = par[4];
    double v2assoc = par[5];

    double flow = bg*(1+v2trigg*v2assoc*cos(2.0*dphi*TMath::Pi()));

    return flow + ampl*pow( 1 + k*dphi*dphi,-n);
}

double CauchyConst(double *x, double *par){
    // Cauchy distribution as in:
    // https://en.wikipedia.org/wiki/Cauchy_distribution
    //double dx    = x[0];
    //double bg    = par[0];
    //double yield = par[1];
    //double g     = par[2];
    //double peak  = yield * 1/pi/g * pow(g,2.0)/(pow(dx,2.0)+pow(g,2.0));

    // Lorentzian peak instead:
    // ftp://root.cern.ch/root/doc/5FittingHistograms.pdf
    double dx = x[0];
    double bg = par[0];
    double A  = par[1];
    double g  = par[2];
    double peak = 0.5*A*g/TMath::Pi() / TMath::Max(1.e-10, dx*dx + 0.25*g*g);

    return bg + peak;
}

double CauchyFlow(double *x, double *par){
    double dx = x[0];
    double bg = par[0];
    double A  = par[1];
    double g  = par[2];
    double v2trigg = par[3];
    double v2assoc = par[4];

    double flow = bg*(1+v2trigg*v2assoc*cos(2.0*dx*TMath::Pi()));

    double peak = 0.5*A*g/TMath::Pi() / TMath::Max(1.e-10, dx*dx + 0.25*g*g);

    return flow + peak;
}

// supplementary function for qgaussian
double e(double x, double q)
{
    double ret = 1.;
    if(q!=1.) ret = pow(1.+(1.-q)*x, 1./(1.-q));
    if(q==1.) ret = exp(x);
    return ret;
}
double QGaussConst(double *x, double *par)
{
    // Fit Q-Gaussian, as described in
    // https://en.wikipedia.org/wiki/Q-Gaussian_distribution
    double dx    = x[0];
    double bg    = par[0];
    double yield = par[1];
    double q     = par[2];
    double beta  = par[3];

//    double c = 1;
//    if(q<1.)
//        c = 2*sqrt(TMath::Pi())*TMath::Gamma(1./(1.-q)) / ( (3.-q)*sqrt(1.-q)*TMath::Gamma((3.-q)/2./(1.-q)) );
//    if(q==1.)
//        c = sqrt(TMath::Pi());
//    if(1.<q && q<3.)
    double c = sqrt(TMath::Pi())*TMath::Gamma((3.-q)/2./(q-1.)) / ( sqrt(q-1.)*TMath::Gamma(1./(q-1.)) );

    double peak = yield * sqrt(beta)/c * e( -beta*dx*dx, q );

    return bg + peak;
}
double QGaussFlow(double *x, double *par)
{
    // Fit Q-Gaussian, as described in
    // https://en.wikipedia.org/wiki/Q-Gaussian_distribution
    double dx    = x[0];
    double bg    = par[0];
    double yield = par[1];
    double q     = par[2];
    double beta  = par[3];
    double v2trigg = par[4];
    double v2assoc = par[5];

//    double c = 1;
//    if(q<1.)
//        c = 2*sqrt(TMath::Pi())*TMath::Gamma(1./(1.-q)) / ( (3.-q)*sqrt(1.-q)*TMath::Gamma((3.-q)/2./(1.-q)) );
//    if(q==1.)
//        c = sqrt(TMath::Pi());
//    if(1.<q && q<3.)
    double c = sqrt(TMath::Pi())*TMath::Gamma((3.-q)/2./(q-1.)) / ( sqrt(q-1.)*TMath::Gamma(1./(q-1.)) );

    double flow = bg*(1+v2trigg*v2assoc*cos(2.0*dx*TMath::Pi()));
    double peak = yield * sqrt(beta)/c * e( -beta*dx*dx, q );

    return flow + peak;
}

class MFit 
{  
private:
    int    binUnderLeft, binUnderRight;
    int    fFuncType, fHistType;
    bool   fIsGauss;
    double underLeft, underRight, nearRange, undereve;
    double deta, dphi;
    double v2;
    double fInt_min, fInt_max;
    TString myName;

public:

    TF1 * ffit;
    double fitmin, fitmax;

    /*
     * DEFAULT CONSTRUCTOR
     */
    MFit() :
        binUnderLeft(-1),
        binUnderRight(-1),
        fFuncType(-1),
        fHistType(-1),
        fIsGauss(true),
        underLeft(0),
        underRight(0),
        nearRange(0),
        undereve(0),
        deta(0),
        dphi(0),
        v2(0),
        fInt_min(0),
        fInt_max(0),
        myName(0),
        fitmin(0),
        fitmax(0)
    {
        std::cout << "MFit::Default constructor\n";
    }

    /*
     * CONSTRUCTOR
     */
    MFit( int fitType, int histType, TH1 * hidFor,  double fmin, double fmax, bool IsGauss = false  ) :
        binUnderLeft(-1),
        binUnderRight(-1),
        fFuncType(fitType),
        fHistType(histType),
        fIsGauss(IsGauss),
        underLeft(0),
        underRight(0),
        nearRange(0),
        undereve(0),
        deta(0),
        dphi(0),
        v2(0),
        fInt_min(0),
        fInt_max(0),
        myName(0),
        fitmin(fmin),
        fitmax(fmax)
    {
        switch( fFuncType )
        {
            case kOneGenGaussConst : InitOneGenGaussConst(hidFor); break;
            case kKaplanConst      : InitKaplanConst(hidFor);      break;
            case kOneGenGaussFlow  : InitOneGenGaussFlow(hidFor);  break;
            case kKaplanFlow       : InitKaplanFlow(hidFor);       break;
            case kCauchyConst      : InitCauchyConst(hidFor);      break;
            case kCauchyFlow       : InitCauchyFlow(hidFor);       break;
            case kQGaussConst      : InitQGaussConst(hidFor);      break;
            case kQGaussFlow       : InitQGaussFlow(hidFor);       break;
            case kTwoGenGaussConst : InitTwoGenGaussConst(hidFor); break;
            default: std::cerr << "MFIT case not recognized...\n";
        }
        switch( fHistType )
        {
            case kDEta:
                fInt_min = 0.0001; fInt_max = 1.4;
                break;
            case kDPhi:
                fInt_min = -0.45; fInt_max = 0.45;
                break;
        }
    }
    /*
     * DESTRUCTOR
     */
    ~MFit()
    {
        if(ffit) delete ffit;
    }

    TString GetName()
    {
        return myName;
    }

    /*
     * Gets yield for the different cases
     */
    double GetYield()
    {
        double yield=-1;
        switch( fFuncType )
        {
            case kKaplanConst : yield=GetYieldIntegral();  break;
            case kKaplanFlow  : yield=GetYieldIntegral();  break;
            case kCauchyConst : yield=GetYieldIntegral();  break;
            case kCauchyFlow  : yield=GetYieldIntegral();  break;
            default: yield=ffit->GetParameter(1); break;
        }
        return yield;
    }
    double GetYieldError( TFitResultPtr ptr )
    {
        double yielderr=-1;
        switch( fFuncType )
        {
            case kKaplanConst : yielderr=GetYiedIntegralError(ptr);  break;
            case kKaplanFlow  : yielderr=GetYiedIntegralError(ptr);  break;
            case kCauchyConst : yielderr=GetYiedIntegralError(ptr);  break;
            case kCauchyFlow  : yielderr=GetYiedIntegralError(ptr);  break;
            default: yielderr = ffit->GetParError(1); break;
        }
        return yielderr;
    }
    double GetYieldIntegral()
    {
        double yield = -1;
        double fitconst = ffit->GetParameter(0);
        ffit->SetParameter(0, 0); // set const. to 0 then integrate
        yield = ffit->Integral(fInt_min, fInt_max);
        ffit->SetParameter(0, fitconst);
        return yield;
    }
    double GetYiedIntegralError( TFitResultPtr ptr )
    {
        double yielderr = -1;
        double fitconst = ffit->GetParameter(0);
        ffit->SetParameter(0, 0); // set const. to 0 then integrate
        yielderr = ffit->IntegralError( fInt_min, fInt_max, ptr->GetParams(), ptr->GetCovarianceMatrix().GetMatrixArray() );
        ffit->SetParameter(0, fitconst);
        return yielderr;
    }
    /*
     * Gets exponent for generalized gaussian
     */
    double GetExpo()
    {
        double expo=-1;
        switch( fFuncType )
        {
            case kOneGenGaussConst : expo=ffit->GetParameter(3); break; // 0
            case kOneGenGaussFlow  : expo=ffit->GetParameter(3); break; // 1
            default: expo = -1; break;
        }
        return expo;
    }
    double GetExpoError()
    {
        double experr=-1;
        switch( fFuncType )
        {
            case kOneGenGaussConst : experr=ffit->GetParError(3); break; // 0
            case kOneGenGaussFlow  : experr=ffit->GetParError(3); break; // 1
            default: experr = -1; break;
        }
        return experr;
    }
    /*
     * Gets width for the different cases
     */
    double GetWidth(){
        double width=-1.;
        switch( fFuncType )
        {
            case kOneGenGaussConst : width = GetWidthOneGenGauss(); break; // 0
            case kOneGenGaussFlow  : width = GetWidthOneGenGauss(); break; // 1
            default: width = -1;
        }
        return width;
    }
    double GetWidthError(TFitResultPtr ptr)
    {
        double widtherr = -1;
        switch( fFuncType )
        {
            case kOneGenGaussConst : widtherr = GetWidthErrorOneGenGauss(ptr); break; // 0
            case kOneGenGaussFlow  : widtherr = GetWidthErrorOneGenGauss(ptr); break; // 1
            default: widtherr = -1;

        }
        return widtherr;
    }
    double GetWidthOneGenGauss() {
        double alpha = ffit->GetParameter(2);
        double beta  = ffit->GetParameter(3);
        return alpha*sqrt( TMath::Gamma(3./beta)/ TMath::Gamma(1./beta) );
    }
    double GetWidthErrorOneGenGauss( TFitResultPtr ptr )
    {
        TMatrixDSym cov = ptr->GetCovarianceMatrix();

        double alpha = ffit->GetParameter(2);
        double beta  = ffit->GetParameter(3);

        double alphaDer = TMath::Sqrt(TMath::Gamma(3./beta)/TMath::Gamma(1./beta));
        TF1* tmp = new TF1("tmp","TMath::Sqrt(TMath::Gamma(3./x)/TMath::Gamma(1./x))",1,2);
        double betaDer = alpha*tmp->Derivative(beta);
        double rmsError =
            TMath::Power(alphaDer * ffit->GetParError(2), 2) +
            TMath::Power(betaDer * ffit->GetParError(3), 2) +
            2. * alphaDer * betaDer * cov(2, 3);
        return TMath::Sqrt(rmsError);
    }


    /*
     * Drawing for the different cases (all formatting is done here)
     * sets also the name for TLegend()
     */
    void Draw(){
        switch( fFuncType ){
            case kOneGenGaussConst : DrawOneGenGaussConst(); break; // 0
            case kKaplanConst      : DrawKaplanConst();      break; // 2
            case kOneGenGaussFlow  : DrawOneGenGaussFlow();  break; // 3
            case kKaplanFlow       : DrawKaplanFlow();       break; // 5
            case kCauchyConst      : DrawCauchyConst();      break; // 6
            case kQGaussConst      : DrawQGaussConst();      break; //
            case kTwoGenGaussConst : DrawTwoGenGaussConst(); break;
        }
    }
    void DrawOneGenGaussConst(){
        ffit->SetLineWidth(2);
        //Sum
        if(fIsGauss)  { ffit->SetLineStyle(3); ffit->SetLineColor(kBlack); ffit->DrawClone("lsame"); }
        if(!fIsGauss) { ffit->SetLineStyle(1); ffit->SetLineColor(kRed+1); ffit->DrawClone("lsame"); }
        //Pedestal
        ffit->SetParameter(1, 0);
        if(fIsGauss)  { ffit->SetLineStyle(3); ffit->SetLineColor(kBlack); ffit->DrawClone("lsame"); }
        if(!fIsGauss) { ffit->SetLineStyle(1); ffit->SetLineColor(kRed+1); ffit->DrawClone("lsame"); }
    }
    void DrawOneGenGaussFlow(){
        ffit->SetLineWidth(2);
        //Sum
        if(fIsGauss)  { ffit->SetLineStyle(3); ffit->SetLineColor(kBlack); ffit->DrawClone("lsame"); }
        if(!fIsGauss) { ffit->SetLineStyle(1); ffit->SetLineColor(kRed+1); ffit->DrawClone("lsame"); }
        //Flow
        TF1 * flow = GetOneGenGaussFlow();
        if(fIsGauss)  { flow->SetLineStyle(3); flow->SetLineColor(kBlack); flow->DrawClone("lsame"); }
        if(!fIsGauss) { flow->SetLineStyle(1); flow->SetLineColor(kRed+1); flow->DrawClone("lsame"); }
        delete flow;
    }
    void DrawKaplanConst(){
        ffit->SetLineWidth(2);
        //Sum
        ffit->SetLineStyle(7); ffit->SetLineColor(kGreen+1); ffit->DrawClone("lsame");
        //Pedestal
        ffit->SetParameter(1, 0);
        ffit->SetLineStyle(7); ffit->SetLineColor(kGreen+1); ffit->DrawClone("lsame");
    }
    void DrawKaplanFlow(){
        ffit->SetLineWidth(2);
        //Sum
        ffit->SetLineStyle(7); ffit->SetLineColor(kGreen+1); ffit->DrawClone("lsame");
        //Flow
        TF1 * flow = GetKaplanFlow();
        flow->SetLineStyle(7); flow->SetLineColor(kGreen+1); flow->DrawClone("lsame");
        delete flow;
    }
    void DrawCauchyConst(){
        //Sum
        ffit->SetLineWidth(2); ffit->SetLineStyle(1); ffit->SetLineColor(kMagenta+1); ffit->DrawClone("lsame");
        //Pedestal
        ffit->SetParameter(1,0); ffit->DrawClone("lsame");
    }
    void DrawQGaussConst(){
        //Sum
        ffit->SetLineWidth(2); ffit->SetLineStyle(7); ffit->SetLineColor(kYellow+2); ffit->DrawClone("lsame");
        //Pedestal
        ffit->SetParameter(1,0); ffit->DrawClone("lsame");
    }
    void DrawTwoGenGaussConst(){
        //Sum
        ffit->SetLineWidth(2); ffit->SetLineStyle(9); ffit->SetLineColor(kYellow+2); ffit->DrawClone("lsame");
        //Pedestal
        ffit->SetParameter(1,0); ffit->SetParameter(4,0); ffit->DrawClone("lsame");
    }



    /*
     * Initialize one generalized gaussian + const
     */
    void InitOneGenGaussConst(TH1 * hidFor ) {
        myName = "Gen.Gauss+const.";
        double lr=fitmin, ur=fitmax;
        double yield, width;

        ffit = new TF1("ffit_GC", OneGenGaussConst, lr , ur, 4 );
        ffit->SetParNames("background", "yield", "width","power");

        if( fHistType==0 ) { // deta background region
            binUnderLeft  = hidFor->FindBin(1.);
            binUnderRight = hidFor->FindBin(1.5);
        }
        if( fHistType==1 ) { // dphi background region
            binUnderLeft  = hidFor->FindBin(0.45);
            binUnderRight = hidFor->FindBin(0.65);
        }
        undereve = hidFor->GetBinContent(binUnderLeft); // apart from stat. fluctuations this is correct
        width = hidFor->GetRMS()/sqrt(2.);
//            yield  = (hidFor->GetBinContent( hidFor->FindBin(0.0) ) -undereve ) * width/2.;
        yield = hidFor->Integral(hidFor->FindBin(0.0), hidFor->FindBin(width))-undereve;
//            if(fHistType==1) yield=2.0*yield; // phi is not flipped

        ffit->SetParameters( undereve, yield, width, 2. );

        //ffit->SetParLimits( 0, 0, 1e9 );
        //ffit->SetParLimits( 1, 0, 1e9 );
        ffit->SetParLimits( 2, 0, 2.0 );
        ffit->SetParLimits( 3, 0.5, 500. );

        if( fIsGauss ) {
            myName = "Gauss+const.";
            ffit->FixParameter(3, 2.);
        }
    }
    /*
     * Initialize one generalized gaussian + flow
     */
    void InitOneGenGaussFlow(TH1 * hidFor ) {
        myName = "Gen.Gauss+flow";
        double lr=fitmin, ur=fitmax;
        double yield, width;

        ffit = new TF1("ffit_GF", OneGenGaussFlow, lr , ur, 6 );
        ffit->SetParNames("background", "yield", "width", "power", "v2ptt", "v2pta");

        if( fHistType==0 ) { // deta background region
            binUnderLeft  = hidFor->FindBin(1.);
            binUnderRight = hidFor->FindBin(1.5);
        }
        if( fHistType==1 ) { // dphi background region
            binUnderLeft  = hidFor->FindBin(0.45);
            binUnderRight = hidFor->FindBin(0.65);
        }
        undereve = hidFor->GetBinContent(binUnderLeft);
        width = hidFor->GetRMS()/sqrt(2.);
//            yield  = (hidFor->GetBinContent( hidFor->FindBin(0.0) ) - undereve ) * width/2.;
        yield = hidFor->Integral(hidFor->FindBin(0.0), hidFor->FindBin(width))-undereve;
        if(fHistType==1) yield=2.0*yield; // phi is not flipped

        // setting arbritary pT, should be customized in setPt(double, double)
        ffit->SetParameters( undereve, yield, width, 2., 8, 4 );

        //ffit->SetParLimits( 0, 0, 1e9);
        //ffit->SetParLimits( 1, 0, 1e9);
        ffit->SetParLimits( 2, 0, 1. );
        ffit->SetParLimits( 3, 0.5, 5. );

        if( fIsGauss ) {
            myName = "Gauss+flow";
            ffit->FixParameter(3, 2.);
        }
    }
    /*
     * Fixing PTt and PTa values for v2 extraction
     */
    void SetPtFlow(double ptt, double pta)
    {
        switch( fFuncType )
        {
            case kOneGenGaussFlow : SetPtOneGenGaussFlow(ptt, pta); break;
            case kKaplanFlow      : SetPtKaplanFlow(ptt, pta); break;
            case kCauchyFlow      : SetPtCauchyFlow(ptt,pta); break;
            case kQGaussFlow      : SetPtQGaussianFlow(ptt,pta); break;
            default: std::cout << "SetPtFlow() case not defined\n"; break;
        }
    }
    void SetPtOneGenGaussFlow(double ptt, double pta){
        double v2trigg = getPublishedV2(ptt);
        double v2assoc = getPublishedV2(pta);
        //std::cout << "setting:ptt:v2ptt: " << ptt << " " << v2trigg << std::endl;
        //std::cout << "setting:pta:v2pta: " << pta << " " << v2assoc << std::endl;

        ffit->FixParameter(4, v2trigg);
        ffit->FixParameter(5, v2assoc);
    }
    void SetPtKaplanFlow(double ptt, double pta){
        double v2trigg = getPublishedV2(ptt);
        double v2assoc = getPublishedV2(pta);

        ffit->FixParameter(4, v2trigg);
        ffit->FixParameter(5, v2assoc);
    }
    void SetPtCauchyFlow(double ptt, double pta)
    {
        double v2trigg = getPublishedV2(ptt);
        double v2assoc = getPublishedV2(pta);
        ffit->FixParameter(3, v2trigg);
        ffit->FixParameter(4, v2assoc);
    }
    void SetPtQGaussianFlow(double ptt, double pta)
    {
        double v2trigg = getPublishedV2(ptt);
        double v2assoc = getPublishedV2(pta);
        ffit->FixParameter(4, v2trigg);
        ffit->FixParameter(5, v2assoc);
    }

    /*
     * Getting either flow component or constant of fit
     */
    TF1 * GetUE()
    {
        switch( fFuncType )
        {
            case kOneGenGaussFlow: return GetOneGenGaussFlow(); break;
            case kKaplanFlow:  return GetKaplanFlow(); break;
            default: return GetConst(); break;
        }
    }
    TF1 * GetKaplanFlow()
    {
        TF1 * kaplanflow = new TF1("kaplanflow", Flow, fitmin, fitmax, 3);
        kaplanflow->SetParameters(ffit->GetParameter(0), ffit->GetParameter(4), ffit->GetParameter(5));
        return kaplanflow;
    }
    TF1 * GetOneGenGaussFlow()
    {
        TF1 * gaussflow = new TF1("gaussflow", Flow, fitmin, fitmax, 3);
        gaussflow->SetParameters(ffit->GetParameter(0), ffit->GetParameter(4), ffit->GetParameter(5));
        return gaussflow;
    }
    TF1 * GetConst()
    {
        TF1 * fconst = new TF1("const", "pol0", fitmin, fitmax);
        fconst->SetParameter(0, ffit->GetParameter(0) );
        return fconst;
    }


    /*
     * Initialize two generalized gaussian + const
     */
    void InitTwoGenGaussConst(TH1 * hidFor ) {
        myName = "2Gen.Gauss+const.";

        double lr=fitmin, ur=fitmax;

        ffit = new TF1("ffit_GGC", TwoGenGaussConst, lr , ur, 7 );
        ffit->SetParNames("background", "yield1", "width1","power1", "yield2", "width2","power2");

        if( fHistType==0 ) { // deta background region
            binUnderLeft  = hidFor->FindBin(1.);
            binUnderRight = hidFor->FindBin(1.5);
        }
        if( fHistType==1 ) { // dphi background region
            binUnderLeft  = hidFor->FindBin(0.45);
            binUnderRight = hidFor->FindBin(0.65);
        }
        undereve = hidFor->Integral(binUnderLeft, binUnderRight)/(binUnderRight-binUnderLeft);

        double width = hidFor->GetRMS()/sqrt(2.);
        double yield  = hidFor->Integral( hidFor->FindBin(0.0), hidFor->FindBin(2.*width))-undereve;

        ffit->SetParameters( undereve, yield/2., width*2., 2., yield/2., width/2., 2. );

        //ffit->SetParLimits( 0, 0, 1e7 );
        //ffit->SetParLimits( 1, 0, 1e7 );
        ffit->SetParLimits( 2, 0, 3. );
        ffit->SetParLimits( 3, 0.5, 5 );
        //ffit->SetParLimits( 4, 0, 1e7 );
        ffit->SetParLimits( 5, 0, 3. );
        ffit->SetParLimits( 6, 0.5, 5 );

        if( fIsGauss ) {
            myName = "2Gauss+const.";
            ffit->FixParameter(3, 2.);
            ffit->FixParameter(6, 2.);
        }
    }
    /*
     * Initialize two generalized gaussian + flow
     */
    void InitTwoGenGaussFlow(TH1 * hidFor){
        myName = "2Gen.Gauss+flow";
        double lr=-0.5, ur=1.5;

        ffit = new TF1("ffit_GGF", TwoGenGaussFlow, lr , ur, 8 );
        ffit->SetParNames("background", "v2", "yield_near", "width_near","power_near", "yield_far", "width_far","power_far");

        binUnderLeft  = hidFor->FindBin(.45);
        binUnderRight = hidFor->FindBin(.65);
        undereve = hidFor->Integral(binUnderLeft, binUnderRight)/(binUnderRight-binUnderLeft);

        double maxVal = hidFor->GetMaximum();
        //double minVal = hidFor->GetMinimum();
        double yield  = maxVal-undereve;

        ffit->SetParameters( undereve, yield, 0.3, 2., yield, 0.3, 2. );

        ffit->SetParLimits(1, 0, 1);
        ffit->SetParLimits( 2, 0, 1e7 );
        ffit->SetParLimits( 3, 0, 1. );
        ffit->SetParLimits( 4, 0.5, 5 );
        ffit->SetParLimits( 5, 0, 1e7 );
        ffit->SetParLimits( 6, 0, 1. );
        ffit->SetParLimits( 7, 0.5, 5 );

        if( fIsGauss ) {
            myName = "2Gauss+flow.";
            ffit->FixParameter(4, 2.);
            ffit->FixParameter(7, 2.);
        }
    }
    /*
     * Initialize Kaplan + const
     */
    void InitKaplanConst(TH1 * hidFor) {
        myName = "Kaplan+const.";
        double lr=fitmin, ur=fitmax;

        ffit = new TF1("ffit_kaplanC", KaplanConst, lr, ur, 4);
        ffit->SetParNames("A", "b", "n", "k");

        if( fHistType==0 ) { // deta background region
            binUnderLeft  = hidFor->FindBin(1.);
            binUnderRight = hidFor->FindBin(1.5);
        }
        if( fHistType==1 ) { // dphi background region
            binUnderLeft  = hidFor->FindBin(0.45);
            binUnderRight = hidFor->FindBin(0.65);
        }
//            undereve = hidFor->Integral(binUnderLeft, binUnderRight)/(binUnderRight-binUnderLeft);
        undereve = hidFor->GetMinimum();

        double maxVal = hidFor->GetMaximum();
        //double minVal = hidFor->GetMinimum();
        double yield  = maxVal-undereve;

        ffit->SetParameters(undereve, yield, 10.0, 1.0);

        ffit->SetParLimits(0, 0, 1e9);
        ffit->SetParLimits(1, 2e-3, 100);
        ffit->SetParLimits(3, 2e-3, 10);
    }
    /*
     * Initialize Kaplan + const
     */
    void InitKaplanFlow(TH1 * hidFor) {
        myName = "Kaplan+flow";
        double lr=fitmin, ur=fitmax;

        ffit = new TF1("ffit_kaplanC", KaplanFlow, lr, ur, 6);
        ffit->SetParNames("A", "b", "n", "k", "v2ptt", "v2pta");

        if( fHistType==0 ) { // deta background region
            binUnderLeft  = hidFor->FindBin(1.);
            binUnderRight = hidFor->FindBin(1.5);
        }
        if( fHistType==1 ) { // dphi background region
            binUnderLeft  = hidFor->FindBin(0.45);
            binUnderRight = hidFor->FindBin(0.65);
        }
        undereve = hidFor->GetMinimum();

        double maxVal = hidFor->GetMaximum();
        //double minVal = hidFor->GetMinimum();
        double yield  = maxVal-undereve;

        ffit->SetParameters(undereve, yield, 20.0, 1.0, 1.0, 1.0);

        ffit->SetParLimits(0, 0, 1e9);
//            ffit->SetParLimits(1, 2e-3, 100);
//            ffit->SetParLimits(3, 2e-3, 10);
    }
    /*
     * Initialize Cauchy + const
     */
    void InitCauchyConst(TH1 * hidFor) {
        myName = "Lorentz+const.";
        double lr=fitmin, ur=fitmax;

        ffit = new TF1("ffit_cauchyC", CauchyConst, lr, ur, 3);
        ffit->SetParNames("background", "yield", "gamma");

        if( fHistType==0 ) { // deta background region
            binUnderLeft  = hidFor->FindBin(1.);
            binUnderRight = hidFor->FindBin(1.5);
        }
        if( fHistType==1 ) { // dphi background region
            binUnderLeft  = hidFor->FindBin(0.45);
            binUnderRight = hidFor->FindBin(0.65);
        }

        double maxVal = hidFor->GetMaximum();
        //double minVal = hidFor->GetMinimum();
        double undereve = hidFor->GetBinContent( binUnderLeft );
        double width = hidFor->GetRMS()/sqrt(2.);
        double yield = hidFor->Integral(hidFor->FindBin(0.0), hidFor->FindBin(width))-undereve;

        ffit->SetParameters(undereve, yield, width);

        //ffit->SetParLimits(0, 0, 1e9);
        //ffit->SetParLimits(1, 0, 10);
        ffit->SetParLimits(2, 0.01, 10);
    }
    /*
     * Initalize Cauchy + flow
     */
    void InitCauchyFlow(TH1 * hidFor){
        myName = "Lorentz+flow.";
        double lr=fitmin, ur=fitmax;

        ffit = new TF1("ffit_cauchyF", CauchyConst, lr, ur, 5);
        ffit->SetParNames("background", "yield", "gamma", "v2ptt", "v2pta");

        if( fHistType==0 ) { // deta background region
            binUnderLeft  = hidFor->FindBin(1.);
            binUnderRight = hidFor->FindBin(1.5);
        }
        if( fHistType==1 ) { // dphi background region
            binUnderLeft  = hidFor->FindBin(0.45);
            binUnderRight = hidFor->FindBin(0.65);
        }

        double maxVal = hidFor->GetMaximum();
        //double minVal = hidFor->GetMinimum();
        double undereve = hidFor->GetBinContent( binUnderLeft );
        double width = hidFor->GetRMS()/sqrt(2.);
        double yield = hidFor->Integral(hidFor->FindBin(0.0), hidFor->FindBin(width))-undereve;

        ffit->SetParameters(undereve, yield, width, 8, 4);

        //ffit->SetParLimits(0, 0, 1e9);
        //ffit->SetParLimits(1, 0, 10);
        ffit->SetParLimits(2, 0.01, 10);
    };
    /*
     * Initialize QGaussian + const
     */
    void InitQGaussConst( TH1 * hidFor )
    {
        myName = "QGaussian+const.";
        double lr=fitmin, ur=fitmax;

        ffit = new TF1("ffit_qgaussC", QGaussConst, lr, ur, 4);
        ffit->SetParNames("background", "yield", "q", "beta");

        if( fHistType==0 ) { // deta background region
            binUnderLeft  = hidFor->FindBin(1.);
            binUnderRight = hidFor->FindBin(1.5);
        }
        if( fHistType==1 ) { // dphi background region
            binUnderLeft  = hidFor->FindBin(0.45);
            binUnderRight = hidFor->FindBin(0.65);
        }

        undereve = hidFor->GetBinContent(binUnderLeft);
        double width = hidFor->GetRMS()/2.;
        double yield = hidFor->Integral(hidFor->FindBin(0.0), hidFor->FindBin(width))-undereve;

        ffit->SetParameters(undereve, yield, 1.1, width);

        ffit->SetParLimits(2, 1, 3);
        ffit->SetParLimits(3, 0, 5);
    }
    /*
     * Initialize QGaussian + flow
     */
    void InitQGaussFlow(TH1 * hidFor)
    {
        myName = "QGaussian+flow.";
        double lr=fitmin, ur=fitmax;

        ffit = new TF1("ffit_qgaussF", QGaussFlow, lr, ur, 6);
        ffit->SetParNames("background", "yield", "q", "beta", "v2ptt", "v2pta");

        if( fHistType==0 ) { // deta background region
            binUnderLeft  = hidFor->FindBin(1.);
            binUnderRight = hidFor->FindBin(1.5);
        }
        if( fHistType==1 ) { // dphi background region
            binUnderLeft  = hidFor->FindBin(0.45);
            binUnderRight = hidFor->FindBin(0.65);
        }

        undereve = hidFor->GetBinContent(binUnderLeft);
        double width = hidFor->GetRMS()/2.;
        double yield = hidFor->Integral(hidFor->FindBin(0.0), hidFor->FindBin(width))-undereve;

        ffit->SetParameters(undereve, yield, 1.1, width, 8, 4);

        ffit->SetParLimits(2, 1, 3);
        ffit->SetParLimits(3, 0, 5);
    }

    int GetFitIndex() { return fFuncType; }
    int GetHistIndex(){ return fHistType; }
};

#endif /* MFIT_H */

bool IsFilipPub(int iptt, int ipta) {
    if(iptt==2) {
        if(ipta==3) return true;
        else if(ipta==4) return true;
        else return false;
    }
    if(iptt==3) {
        if(ipta==4) return true;
        else if(ipta==5) return true;
        else return false;
    }
    return false;
}

// -------------------------------------
// Filip's prelimiary IAA(pta)
// obtained with plotdigitizer...
// from fig. IAA.png in talks/introfigs
// -------------------------------------
TGraphAsymmErrors * filipIAApta(bool stat)
{
	TGraphAsymmErrors * fiaa = nullptr;

	double xval[] = { 0.85, 1.5, 2.5, 3.5, 5., 7 };
	double yval[] = {1.53677, 1.41644, 1.32706, 1.12765, 1.2136, 1.1861};
	double xerrminus[]     = { 0.15, 0.5, 0.5, 0.5, 1.0, 1.0 }; // symm. bins, same as plus
	double ystaterrminus[] = { 0.28535, 0.10314, 0.06188, 0.03782, 0.03782, 0.04469 }; // symm. bins, same as plus
	double ysysterrplus[]  = { 0.56383, 0.13065, 0.11001, 0.08939, 0.09627, 0.09282 };
	double ysysterrminus[] = { 0.15127, 0.07563, 0.06533, 0.05844, 0.05157, 0.04126 };

	if(stat)

	{
		fiaa = new TGraphAsymmErrors(6, xval, yval, xerrminus, xerrminus, ystaterrminus, ystaterrminus);
		fiaa->SetName("gFilipPrelimIAAstat");
	}
	if(!stat)
	{
		fiaa = new TGraphAsymmErrors(6, xval, yval, xerrminus, xerrminus, ysysterrminus, ysysterrplus );
		fiaa->SetName("gFilipPrelimIAAsyst");
		fiaa->SetFillColorAlpha(41,0.2) ;
	}
	fiaa->SetTitle("Filip's preliminary I_{AA}(pta)");

	return fiaa;
}
// -------------------------------------
// published near-side I_AA, pTt 8-15GeV, C:0-5%
// http://arxiv.org/pdf/1110.0121.pdf
// -------------------------------------
TGraphAsymmErrors * publishedIAA(int icase)
{
	TGraphAsymmErrors * g = nullptr; // return graph

	double xval[] = {3.5, 5.0, 7.0, 9.0};
	double xerrminus[] = {0.5, 1.0, 1.0, 1.0};
	double xerrplus[] = {0.5, 1.0, 1.0, 1.0};

	double yval1[] = {1.4, 1.21, 1.12, 1.25};
	double yerrminus1[] = {0.1118033988749895, 0.09848857801796104, 0.08944271909999159, 0.12041594578792295};
	double yerrplus1[] = {0.1118033988749895, 0.09848857801796104, 0.08944271909999159, 0.12041594578792295};
	int numpoints = 4;
	TGraphAsymmErrors * p8088_d1x1y1 = new TGraphAsymmErrors(numpoints, xval, yval1, xerrminus, xerrplus, yerrminus1, yerrplus1);
	p8088_d1x1y1->SetName("/HepData/8088/d1x1y1");
	p8088_d1x1y1->SetTitle("PRL flat backg.");

	double yval2[] = {1.29, 1.19, 1.12, 1.25};
	double yerrminus2[] ={0.10295630140987, 0.08544003745317531, 0.08944271909999159, 0.12041594578792295};
	double yerrplus2[] = {0.10295630140987, 0.08544003745317531, 0.08944271909999159, 0.12041594578792295};
	TGraphAsymmErrors * p8088_d1x1y2 = new TGraphAsymmErrors(numpoints, xval, yval2, xerrminus, xerrplus, yerrminus2, yerrplus2);
	p8088_d1x1y2->SetName("/HepData/8088/d1x1y2");
	p8088_d1x1y2->SetTitle("PRL v2 backg.");

	double yval3[] = {1.25, 1.22, 1.09, 1.3};
	double yerrminus3[] = {0.1140175425099138, 0.09848857801796104, 0.09433981132056604, 0.12727922061357855};
	double yerrplus3[] = {0.1140175425099138, 0.09848857801796104, 0.09433981132056604, 0.12727922061357855};
	TGraphAsymmErrors * p8088_d1x1y3 = new TGraphAsymmErrors(numpoints, xval, yval3, xerrminus, xerrplus, yerrminus3, yerrplus3);
	p8088_d1x1y3->SetName("/HepData/8088/d1x1y3");
	p8088_d1x1y3->SetTitle("PRL #eta-gap backg.");

	if(icase == 0)
		g = new TGraphAsymmErrors( *p8088_d1x1y1 ); // 0: flat back.
	if(icase == 1)
		g = new TGraphAsymmErrors( *p8088_d1x1y2 ); // 1: v2 backg.
	if(icase == 2)
		g = new TGraphAsymmErrors( *p8088_d1x1y3 ); // 2: etaGap backg.

	return g;
}
// -------------------------------------
// Published $I_{\rm AA}$ from pi0-hadron
// correlation: arxiv:1608.07201
// -------------------------------------
TH1F * publishedIAA_pi0()
{
	// TODO handle syst errors
	TFile * pubFile = TFile::Open("../data/published/HEPData-ins1483164-1-root.root");
	TH1F * h = (TH1F*) pubFile->Get("Table 3/Hist1D_y1");
	return h;
}

std::vector<double> ptt = {3, 4, 6, 8, 15};
//std::vector<double> pta = {0.6, 1, 2, 3, 4, 6, 8, 10};
std::vector<double> pta = {0.6, 1, 2, 3, 4, 6, 8};
std::vector<double> cent = {0, 10, 20, 40, 60, 90};

TString CentTitle(int ic)
{
	return Form("C: %.0f-%.0f %%",cent.at(ic), cent.at(ic+1));
}
TString PTaTitle(int ipta)
{
	if(ipta==0) return Form("p_{Ta}#in %.1f-%.0f GeV", pta.at(ipta), pta.at(ipta+1));
	else return Form("p_{Ta}#in %.0f-%.0f GeV",  pta.at(ipta), pta.at(ipta+1));
}
TString PTtTitle(int iptt)
{
	return Form("p_{Tt}#in %.0f-%.0f GeV",  ptt.at(iptt), ptt.at(iptt+1));
}

// converts vector of histograms into one graph
// with systematic error. Base histogram is 0th element
TGraphAsymmErrors * systerr(std::vector<TH1D*> hsys)
{
    std::cout << "processing systerr() with " << hsys.size() << " parameters\n";

	double x=0;
	double y=0;
	double exl=0;
	double exh=0;
	double eyl=0;
	double eyh=0;

	const int nbins = hsys[0]->GetNbinsX();
	TString name = Form("%s_systerr",hsys.at(0)->GetName());
	TGraphAsymmErrors * g = new TGraphAsymmErrors(nbins);
	g->SetName(name.Data());
	std::vector<double> systerr;
	for(int ib=1; ib<=nbins; ib++)
	{
		systerr.clear();
		for(auto ih : hsys)
			systerr.push_back( ih->GetBinContent(ib) );

		x = hsys[0]->GetBinCenter(ib);
		y = hsys[0]->GetBinContent(ib);
		exl = x-hsys[0]->GetXaxis()->GetBinLowEdge(ib);
		exh = hsys[0]->GetXaxis()->GetBinUpEdge(ib)-x;
		auto err = std::minmax_element(systerr.begin(), systerr.end());
		eyl = y-TMath::Abs(*(err.first));
		eyh = *(err.second)-y;
		g->SetPoint(ib-1, x, y);
		g->SetPointError(ib-1, exl, exh, eyl, eyh);
	}
	return g;
}




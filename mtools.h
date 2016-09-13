#ifndef MTOOLS_H
#define MTOOLS_H

#include "mfit.h"

/* 
 * ***************************************************************************
 * This class is my toolbox mostly for ROOT-related operations
 * (e.g. shrink, flip histograms, etc) 
 * ***************************************************************************
 */

class MTools
{

private:
    // ensures unique naming
    // (problems e.g. flipping two histograms 
    // with same name in a loop... )
    int mToolsIndex;

public:
    MTools() { mToolsIndex=0; }

    // --------------------------------------------
    // Shrinks x axis of a histogram (h) with (sc)
    // --------------------------------------------
    TH1 * shrinkHist( TH1 * h, double sc){
        int nbins   = h->GetNbinsX();
        double xmin = h->GetXaxis()->GetXmin()/sc;
        double xmax = h->GetXaxis()->GetXmax()/sc;
        TH1D * hn = new TH1D(Form("%s_shrink",h->GetName()),"", nbins, xmin, xmax);
        for(int ib=0; ib<nbins; ib++) {
            hn->SetBinContent(ib+1, h->GetBinContent(ib+1));
            hn->SetBinError(ib+1, h->GetBinError(ib+1));
        }
        return hn;
    }
    // --------------------------------------------
    // Modifies given histogram (h) to have
    // Zero Yield at Minimum
    // --------------------------------------------
    void zyam(TH1 * h){
        double ymin = h->GetMinimum();
        double y, yerr;
        for(int ib=0; ib<h->GetNbinsX(); ib++){
            y = h->GetBinContent(ib+1)-ymin;
            yerr = h->GetBinError(ib+1)/h->GetBinContent(ib+1);

            h->SetBinContent(ib+1, y);
            h->SetBinError(ib+1, yerr*y);
        }
    }
    // --------------------------------------------
	// Substracts background (bg)
	// from a given histogram (h)
    // --------------------------------------------
	void subtractConstTH1(TH1 * h, double bg)
	{
        double y, yerr;
		for(int ib=1; ib<=h->GetNbinsX(); ib++){
			y = h->GetBinContent(ib)-bg;
			yerr = h->GetBinError(ib)/h->GetBinContent(ib);

			h->SetBinContent(ib, y);
			h->SetBinError(ib, yerr*y);
        }
    }
	TH1D* subtractConstTH1(TH1 *h, double bg, double ebg)
	{
		int nb =h->GetNbinsX();
		TString hname = h->GetName();

		TH1D *hsig  = (TH1D*) h->Clone(Form("%s_sig",hname.Data()));
		hsig->Reset();

		for(int ib=1; ib<=nb; ib++){
			double val = h->GetBinContent(ib);
			double err = h->GetBinError(ib);
			hsig->SetBinContent(ib,val-bg);
			hsig->SetBinError(ib,sqrt(err*err +  ebg * ebg));
		}
		return hsig;
	}
	void subtractConstTH1(TH1 * h) // fit and subtract constant from deta histogram
	{
		MFit * fit = new MFit(0, 0, h, 0, 1.6, false); // larger range for better backg. estimation
		h->Fit( fit->ffit, "IEMR");
		double bg = fit->ffit->GetParameter(0);
		double bgerr = fit->ffit->GetParError(0);
		subtractConstTH1(h, bg);
		//return h;
		delete fit;
	}



    // --------------------------------------------
    // Normalize given histogram (h) to its tail,
    // the tail is defined as range: xmin-xmax
    // --------------------------------------------
    void normalizeToTail(TH1 * h, double xmin, double xmax){
        int xBinMin = h->FindBin(xmin);
        int xBinMax = h->FindBin(xmax);
        double scale = h->Integral(xBinMin, xBinMax)/float(xBinMax-xBinMin);
        h->Scale(1./scale);
    }
    // --------------------------------------------
	// Shift histogram to another histogram's tail,
    // tail is defined similarly
    // --------------------------------------------
    void shiftToThisTail(TH1 * h2, TH1 * h1, double xmin, double xmax){
        int xBinMin1 = h1->FindBin(xmin);
        int xBinMax1 = h1->FindBin(xmax);
        int xBinMin2 = h2->FindBin(xmin);
        int xBinMax2 = h2->FindBin(xmax);

        double shift1 = h1->Integral(xBinMin1, xBinMax1)/float(xBinMax1-xBinMin1);
        double shift2 = h2->Integral(xBinMin2, xBinMax2)/float(xBinMax2-xBinMin2);
        for(int ib=0; ib<h2->GetNbinsX(); ib++){
            h2->SetBinContent(ib+1, h2->GetBinContent(ib+1)+shift1-shift2);
        }
    }
	// --------------------------------------------
	// Scale histogram to another histogram's tail,
	// tail is defined similarly
	// --------------------------------------------
	void scaleToThisTail(TH1 * h2, TH1 * h1, double xmin, double xmax) {
		int xBinMin1 = h1->FindBin(xmin);
		int xBinMax1 = h1->FindBin(xmax);
		int xBinMin2 = h2->FindBin(xmin);
		int xBinMax2 = h2->FindBin(xmax);

		double shift1 = h1->Integral(xBinMin1, xBinMax1)/double(xBinMax1-xBinMin1);
		double shift2 = h2->Integral(xBinMin2, xBinMax2)/double(xBinMax2-xBinMin2);
		h2->Scale(1./shift2*shift1);
	}

    // --------------------------------------------
    // Returns a histogram, which flips the given 
    // histogram (h) around 0 (and has half the bins)
    // --------------------------------------------
    TH1D* Flip( TH1D* hin ) {
        int nb  = hin->GetNbinsX();
        int zero = hin->FindBin(0.0000001);
        double max = hin->GetBinLowEdge(nb+1);
        double valPos, valNeg, errPos, errNeg, err;
        TString hname = hin->GetName();
        TString newName = Form("%s_flip_%d",hname.Data(), ++mToolsIndex);

        TH1D * hout = new TH1D(newName.Data(), newName.Data(), (int) nb/2, 0, max);

        for(int ib=zero; ib<=nb; ib++){
            valPos = hin->GetBinContent(ib);
            valNeg = hin->GetBinContent(nb - ib+1);
            errPos = hin->GetBinError(ib);
            errNeg = hin->GetBinError(nb - ib+1);

            if(valPos!=0 && valNeg!=0) err = sqrt( errPos/valPos*errPos/valPos + errNeg/valNeg*errNeg/valNeg );
            else err = 0;
            hout->SetBinContent(ib-zero+1, (valPos+valNeg));
            hout->SetBinError(ib-zero+1, (valPos+valNeg) * err);
        }
        return hout;
    }
    // --------------------------------------------
    // Gets formatted text Chi2/NDF from a given function
    // --------------------------------------------
    TString GetChi(TF1 * ffit){
        double chi2 = ffit->GetChisquare();
        int ndf     = ffit->GetNDF();

        TString text = Form("#chi^{2}/NDF = %.1f/%d", chi2, ndf);
        return text;
    }
    // --------------------------------------------
    // Divides two histogram with different binning
    // (why it is not defined in ROOT?!)
    // --------------------------------------------
    void DivideDiffBin(TH1* h1, TH1* h2){
        int nb1 = h1->GetNbinsX();
        int nb2 = h2->GetNbinsX();

        if(nb1 >= nb2) {
            for(int ib=0; ib<nb1; ib++){
                double x  = h1->GetBinCenter(ib);
                double y1 = h1->GetBinContent(ib);
                double e1 = h1->GetBinError(ib);
                double y2 = h2->GetBinContent(h2->FindBin(x));
                double e2 = h2->GetBinError(h2->FindBin(x));
                h1->SetBinContent(ib,y1/y2);
                if(y1>0 && y2>0){
                    h1->SetBinError(ib, y1/y2*sqrt(e1*e1/y1/y1 + e2*e2/y2/y2));
                }else{
                    h1->SetBinError(ib,0);
                }
            }
        }
        if(nb1 < nb2) {
            for(int ib=0; ib<nb2; ib++){
                double x  = h2->GetBinCenter(ib);
                double y2 = h2->GetBinContent(ib);
                double e2 = h2->GetBinError(ib);
                double y1 = h1->GetBinContent(h1->FindBin(x));
                double e1 = h1->GetBinError(h1->FindBin(x));
                h1->SetBinContent(ib,y1/y2);
                if(y1>0 && y2>0){
                    h1->SetBinError(ib, y1/y2*sqrt(e1*e1/y1/y1 + e2*e2/y2/y2));
                }else{
                    h1->SetBinError(ib,0);
                }
            }
        }
    }
    // --------------------------------------------
    // Creates a histogram with unique binning and 
    // copies the content of a given histogram to that
    // (New histogram should have fewer bins.)
    // Errors are averaged
    // --------------------------------------------
    TH1D * NewHistoWithUniqueBins(TH1D * h, int nbins, double xbins[] ){
//        std::cout << "DEBUG: setting histogram: y \t yerr \t yerr/y\n";

        int iFirstBin=0, iLastBin=0;
        double xl_old=0, xh_old=0, y_old=0, xl_new=0, xh_new=0, area=0;
        double y_new=0, y_err=0;
        TString newname = Form("%s_newbin_%d", h->GetName(), ++mToolsIndex);

        TH1D * hnew = new TH1D( newname, "", nbins, xbins );
        hnew->Sumw2();

        std::vector<int> ibinInRange;
        std::vector<double> xInRange;
        std::vector<double> yInRange;

        for(int inewbin=1; inewbin<=nbins; inewbin++)
        {
            y_err = 0; // reset errors
            xl_new = hnew->GetBinLowEdge(inewbin);
            xh_new = hnew->GetBinLowEdge(inewbin+1);

            ibinInRange.clear();
            xInRange.clear();
            yInRange.clear();
            
            // store bins which are inside or overlap with the new bin
            for(int ioldBin=1; ioldBin<=h->GetNbinsX(); ioldBin++)
            {
                xl_old = h->GetBinLowEdge(ioldBin);
                xh_old = h->GetBinLowEdge(ioldBin+1);
                y_old = h->GetBinContent(ioldBin);

                if( xl_new>xh_old || xh_new<xl_old ) 
                    continue;

                ibinInRange.push_back(ioldBin);
            }
            // calculate area for "rebin"...
            iFirstBin = ibinInRange[0];
            iLastBin  = ibinInRange[ibinInRange.size()-1];
            area = (h->GetBinLowEdge(iFirstBin+1)-xl_new)*h->GetBinContent(iFirstBin) + (xh_new-h->GetBinLowEdge(iLastBin))*h->GetBinContent(iLastBin);
            // ... and add bins which are contained in the new bin
            if( ibinInRange.size()>2 ){
                for(unsigned long ibin=1; ibin<ibinInRange.size()-1; ibin++) {
                    area += h->GetBinWidth( ibinInRange.at(ibin) ) * h->GetBinContent( ibinInRange.at(ibin) );
                }
            }
            y_new = area/hnew->GetBinWidth(inewbin);
            hnew->SetBinContent(inewbin, y_new);

            // error is average of bins of the bins contained 
            for(unsigned long ibin=1; ibin<ibinInRange.size(); ibin++){ 
                //yerr += h->GetBinError(ibin)/h->GetBinContent(ibin) / double( ibinInRange.size() ); 
                y_err += h->GetBinError( ibinInRange.at(ibin) ) / double( ibinInRange.size() ); 
            }
            hnew->SetBinError(inewbin, y_err);
//            cout << inewbin << "\t" << y_new << "\t" << y_err << "\t" << y_err/y_new << endl;
        }
        // divide with bin width to get back the y from the area
//        hnew->Scale(1., "width"); 

        return hnew;
    }

    // Rebins a histogram according to other histograms bins
    TH1D * RebinHistoToOther(TH1D * h_value, TH1D * h_bins) 
   {
        const int nbins = h_bins->GetNbinsX();

        double xbins[nbins+1];
        for(int ib=0; ib<nbins; ib++){ xbins[ib] = h_bins->GetXaxis()->GetBinLowEdge(ib); }
        xbins[nbins] = h_bins->GetXaxis()->GetBinLowEdge(nbins);
        
        int iFirstBin=0, iLastBin=0;
        double xl_old=0, xh_old=0, y_old=0, xl_new=0, xh_new=0, area=0;
        double y_new=0, y_err=0;
        TString newname = Form("%s_newbin_%d", h_value->GetName(), ++mToolsIndex);

        TH1D * hnew = new TH1D( newname, "", nbins, xbins );
        hnew->Sumw2();

        std::vector<int> ibinInRange;
        std::vector<double> xInRange;
        std::vector<double> yInRange;

        for(int inewbin=1; inewbin<=nbins; inewbin++)
        {
            y_err = 0; // reset errors
            xl_new = hnew->GetBinLowEdge(inewbin);
            xh_new = hnew->GetBinLowEdge(inewbin+1);

            ibinInRange.clear();
            xInRange.clear();
            yInRange.clear();
            
            // store bins which are inside or overlap with the new bin
            for(int ioldBin=1; ioldBin<=h_value->GetNbinsX(); ioldBin++)
            {
                xl_old = h_value->GetBinLowEdge(ioldBin);
                xh_old = h_value->GetBinLowEdge(ioldBin+1);
                y_old = h_value->GetBinContent(ioldBin);

                if( xl_new>xh_old || xh_new<xl_old ) 
                    continue;

                ibinInRange.push_back(ioldBin);
            }
            // calculate area for "rebin"...
            iFirstBin = ibinInRange[0];
            iLastBin  = ibinInRange[ibinInRange.size()-1];
            area = (h_value->GetBinLowEdge(iFirstBin+1)-xl_new)*h_value->GetBinContent(iFirstBin) 
                + (xh_new-h_value->GetBinLowEdge(iLastBin))*h_value->GetBinContent(iLastBin);
            // ... and add bins which are contained in the new bin
            if( ibinInRange.size()>2 )
            {
                for(unsigned long ibin=1; ibin<ibinInRange.size()-1; ibin++) 
                {
                    area += h_value->GetBinWidth( ibinInRange.at(ibin) ) * h_value->GetBinContent( ibinInRange.at(ibin) );
                }
            }
            y_new = area/hnew->GetBinWidth(inewbin);
            hnew->SetBinContent(inewbin, y_new);

            // error is average of bins of the bins contained 
            for(unsigned long ibin=1; ibin<ibinInRange.size(); ibin++)
            { 
                y_err += h_value->GetBinError( ibinInRange.at(ibin) ) / double( ibinInRange.size() ); 
            }
            hnew->SetBinError(inewbin, y_err);
        }
        return hnew;
    }
    // Project TH2 to either X or Y axis along with an ellipse defined by major axises 
    TH1D * DoProjectionCircle(TH2D * h2, bool onX, const char *name, Int_t firstbin, Int_t lastbin, Double_t R, Option_t *option) const
    {
       const char *expectedName = 0;
       Int_t inNbin;
       Int_t firstOutBin, lastOutBin;
       const TAxis* outAxis;
       const TAxis* inAxis;

       TString opt = option;
       TString cut;
       Int_t i1 = opt.Index("[");
       if (i1>=0) {
          Int_t i2 = opt.Index("]");
          cut = opt(i1,i2-i1+1);
       }
       opt.ToLower();  //must be called after having parsed the cut name
       bool originalRange = opt.Contains("o");

       if ( onX )
       {
          expectedName = "_px";
          inNbin = h2->GetYaxis()->GetNbins();
          outAxis = h2->GetXaxis();
          inAxis = h2->GetYaxis();
       }
       else
       {
          expectedName = "_py";
          inNbin = h2->GetXaxis()->GetNbins();
          outAxis = h2->GetYaxis();
          inAxis = h2->GetXaxis();
       }

       firstOutBin = outAxis->GetFirst();
       lastOutBin = outAxis->GetLast();
       if (firstOutBin == 0 && lastOutBin == 0) {
          firstOutBin = 1; lastOutBin = outAxis->GetNbins();
       }

       if ( lastbin < firstbin && inAxis->TestBit(TAxis::kAxisRange) ) {
          firstbin = inAxis->GetFirst();
          lastbin = inAxis->GetLast();
          // For special case of TAxis::SetRange, when first == 1 and last
          // = N and the range bit has been set, the TAxis will return 0
          // for both.
          if (firstbin == 0 && lastbin == 0)
          {
             firstbin = 1;
             lastbin = inAxis->GetNbins();
          }
       }
       if (firstbin < 0) firstbin = 0;
       if (lastbin  < 0) lastbin  = inNbin + 1;
       if (lastbin  > inNbin+1) lastbin  = inNbin + 1;

       // Create the projection histogram
       char *pname = (char*)name;
       if (name && strcmp(name,expectedName) == 0) {
          Int_t nch = strlen(h2->GetName()) + 4;
          pname = new char[nch];
          snprintf(pname,nch,"%s%s",h2->GetName(),name);
       }
       TH1D *h1=0;
       //check if histogram with identical name exist
       // if compatible reset and re-use previous histogram
       // (see https://savannah.cern.ch/bugs/?54340)
       TObject *h1obj = gROOT->FindObject(pname);
       if (h1obj && h1obj->InheritsFrom(TH1::Class())) {
          if (h1obj->IsA() != TH1D::Class() ) {
             Error("DoProjection","Histogram with name %s must be a TH1D and is a %s",name,h1obj->ClassName());
             return 0;
          }
          h1 = (TH1D*)h1obj;
          // reset the existing histogram and set always the new binning for the axis
          // This avoid problems when the histogram already exists and the histograms is rebinned or its range has changed
          // (see https://savannah.cern.ch/bugs/?94101 or https://savannah.cern.ch/bugs/?95808 )
          h1->Reset();
          const TArrayD *xbins = outAxis->GetXbins();
          if (xbins->fN == 0) {
             if ( originalRange )
                h1->SetBins(outAxis->GetNbins(),outAxis->GetXmin(),outAxis->GetXmax());
             else
                h1->SetBins(lastOutBin-firstOutBin+1,outAxis->GetBinLowEdge(firstOutBin),outAxis->GetBinUpEdge(lastOutBin));
          } else {
             // case variable bins
             if (originalRange )
                h1->SetBins(outAxis->GetNbins(),xbins->fArray);
             else
                h1->SetBins(lastOutBin-firstOutBin+1,&xbins->fArray[firstOutBin-1]);
          }
       }

       Int_t ncuts = 0;
       //if (opt.Contains("[")) {
       //   TVirtualHistPainter * painter  = ((TH2 *)h2)->GetPainter();
       //   if (painter) ncuts = painter->MakeCuts((char*)cut.Data());
       //}

       if (!h1) {
          const TArrayD *bins = outAxis->GetXbins();
          if (bins->fN == 0) {
             if ( originalRange )
                h1 = new TH1D(pname,h2->GetTitle(),outAxis->GetNbins(),outAxis->GetXmin(),outAxis->GetXmax());
             else
                h1 = new TH1D(pname,h2->GetTitle(),lastOutBin-firstOutBin+1,
                              outAxis->GetBinLowEdge(firstOutBin),outAxis->GetBinUpEdge(lastOutBin));
          } else {
             // case variable bins
             if (originalRange )
                h1 = new TH1D(pname,h2->GetTitle(),outAxis->GetNbins(),bins->fArray);
             else
                h1 = new TH1D(pname,h2->GetTitle(),lastOutBin-firstOutBin+1,&bins->fArray[firstOutBin-1]);
          }
          if (opt.Contains("e") || h1->GetSumw2N() ) h1->Sumw2();
       }
       if (pname != name)  delete [] pname;

       // Copy the axis attributes and the axis labels if needed.
       h1->GetXaxis()->ImportAttributes(outAxis);
//       THashList* labels=outAxis->GetLabels();
//       if (labels) {
//          TIter iL(labels);
//          TObjString* lb;
//          Int_t i = 1;
//          while ((lb=(TObjString*)iL())) {
//             h1->GetXaxis()->SetBinLabel(i,lb->String().Data());
//             i++;
//          }
//       }

       h1->SetLineColor(h2->GetLineColor());
       h1->SetFillColor(h2->GetFillColor());
       h1->SetMarkerColor(h2->GetMarkerColor());
       h1->SetMarkerStyle(h2->GetMarkerStyle());

       // Fill the projected histogram
       Double_t cont,err2;
       Double_t totcont = 0;
       Bool_t  computeErrors = h1->GetSumw2N();

       Double_t rx=0, ry=0;
       // implement filling of projected histogram
       // outbin is bin number of outAxis (the projected axis). Loop is done on all bin of TH2 histograms
       // inbin is the axis being integrated. Loop is done only on the selected bins
       for ( Int_t outbin = 0; outbin <= outAxis->GetNbins() + 1;  ++outbin) {
          err2 = 0;
          cont = 0;
          if (outAxis->TestBit(TAxis::kAxisRange) && ( outbin < firstOutBin || outbin > lastOutBin )) continue;

          for (Int_t inbin = firstbin ; inbin <= lastbin ; ++inbin) {
             Int_t binx, biny;
             if (onX) { binx = outbin; biny=inbin; }
             else     { binx = inbin;  biny=outbin; }

             // do the circular cut in R
             // phi needs to be multiplied by PI
             if(onX) {
                 rx = outAxis->GetBinCenter(outbin);
                 ry = inAxis->GetBinCenter(inbin) * TMath::Pi();
             }
             else {
                 rx = outAxis->GetBinCenter(outbin) * TMath::Pi();
                 ry = inAxis->GetBinCenter(inbin);
             }
             if( sqrt(rx*rx + ry*ry) > R)
                  continue;


//             if (ncuts) {
//                if (!fPainter->IsInside(binx,biny)) continue;
//             }
             // sum bin content and error if needed
             cont  += h2->GetBinContent(binx,biny);
             if (computeErrors) {
                Double_t exy = h2->GetBinError(binx,biny);
                err2  += exy*exy;
             }
          }
          // find corresponding bin number in h1 for outbin
          Int_t binOut = h1->GetXaxis()->FindBin( outAxis->GetBinCenter(outbin) );
          h1->SetBinContent(binOut ,cont);
          if (computeErrors) h1->SetBinError(binOut,TMath::Sqrt(err2));
          // sum  all content
          totcont += cont;
       }

       // check if we can re-use the original statistics from  the previous histogram
       bool reuseStats = false;
       if ( ( firstbin == 1 && lastbin == inNbin     ) ||
            ( firstbin == 0 && lastbin == inNbin + 1 ) )
          reuseStats = true;
       else {
          // also if total content match we can re-use
          double eps = 1.E-12;
          if (h2->IsA() == TH2F::Class() ) eps = 1.E-6;
          if (TMath::Abs( h2->GetEntries() - totcont) <  TMath::Abs(h2->GetEntries()) * eps)
             reuseStats = true;
       }
       if (ncuts) reuseStats = false;
       // retrieve  the statistics and set in projected histogram if we can re-use it
       bool reuseEntries = reuseStats;
       // can re-use entries if underflow/overflow are included
       reuseEntries &= (firstbin==0 && lastbin == inNbin+1);
       if (reuseStats) {
          Double_t stats[13];
          h2->GetStats(stats);
          if (!onX) {  // case of projection on Y
             stats[2] = stats[4];
             stats[3] = stats[5];
          }
          h1->PutStats(stats);
       }
       else {
          // the statistics is automatically recalulated since it is reset by the call to SetBinContent
          // we just need to set the entries since they have not been correctly calculated during the projection
          // we can only set them to the effective entries
          h1->SetEntries( h1->GetEffectiveEntries() );
       }
       if (reuseEntries) {
          h1->SetEntries(h2->GetEntries());
       }
       else {
          // re-compute the entries
          // in case of error calculation (i.e. when Sumw2() is set)
          // use the effective entries for the entries
          // since this  is the only way to estimate them
          Double_t entries =  TMath::Floor( totcont + 0.5); // to avoid numerical rounding
          if (h1->GetSumw2N()) entries = h1->GetEffectiveEntries();
          h1->SetEntries( entries );
       }

       if (opt.Contains("d")) {
          TVirtualPad *padsav = gPad;
          TVirtualPad *pad = gROOT->GetSelectedPad();
          if (pad) pad->cd();
          opt.Remove(opt.First("d"),1);
          // remove also other options
          if (opt.Contains("e")) opt.Remove(opt.First("e"),1);
          if (!gPad || !gPad->FindObject(h1)) {
             h1->Draw(opt);
          } else {
             h1->Paint(opt);
          }
          if (padsav) padsav->cd();
       }

       return h1;
    }

	void DivideWithX(TH1D * h)
	{
		const int n = h->GetNbinsX();
		double xval = 0;
		double yval = 0;
		for(int ib=1; ib<=n; ib++)
		{
			xval = h->GetBinCenter(ib);
			yval = h->GetBinContent(ib);
			h->SetBinContent(ib, yval/xval);
		}
	}
	const TH1D * TGraph2TH1D(const TGraphAsymmErrors * gr)
	{
		const int nbins = gr->GetN();
		double x[nbins+1], y[nbins], yerr[nbins];
		double xerrl=0., xerrh=0.;
		for(int ib=0; ib<nbins; ib++)
		{
			(*gr).GetPoint(ib, x[ib], y[ib]);
			xerrl = (*gr).GetErrorXlow(ib);
			xerrh = (*gr).GetErrorXhigh(ib);
			x[ib] = x[ib]-xerrl;
			yerr[ib] = (*gr).GetErrorY(ib);
		}
		x[nbins] = x[nbins-1]+xerrh;

		TH1D * h = new TH1D(Form("h%s",gr->GetName()), "", nbins-1, x);
		for(int ib=0; ib<nbins; ib++)
		{
			//cout << x[ib] << " - " << x[ib+1] << endl;
			h->SetBinContent(ib+1, y[ib]);
			h->SetBinError(ib+1, yerr[ib]);
		}
		h->SetTitle( gr->GetTitle() );
		return h;
	}

//    TH1D * ProjectionCircX(TH2D * h2, const char *name, Int_t firstbin, Int_t lastbin, Double_t R, Option_t * option) const
//    {
//        DoProjectionCircle(h2, true, name, firstbin, lastbin, R, option);
//    }
//
//   TH1D * ProjectionCircY(TH2D * h2, const char *name, Int_t firstbin, Int_t lastbin, Double_t R, Option_t * option) const
//    {
//        DoProjectionCircle(h2, false, name, firstbin, lastbin, R, option);
//    }

};

#endif /* MTOOLS_H */

#ifndef MPLOT_H
#define MPLOT_H

#include <string>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <array>

#include "TH1.h"
#include "TString.h"
#include "TLegend.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TStyle.h"
#include "TF1.h"

/*
 * ***************************************************************************
 * This class contains all style formatting and abbreviations for an easier plot. 
 * The goal is to set the most common practices as default for a shorter plotting macro. 
 * More specific plots requiring additional setup can be also set here.
 * ***************************************************************************
 */

//int Colors[] = { 634, 617, 601, 433, 417, 635, 619, 603, 435, 419 };
std::array<int,10> Colors { { 601, 634, 417, 433, 635, 619, 603, 636, 719, 435} };
std::array<int,13> Markers_full { { 20, 21, 25, 22, 26, 23, 27, 33, 27, 20, 24, 25, 26 } }; 
std::array<int,13> Markers_open { { 24, 25, 26, 27, 33, 27, 20, 24, 25, 26, 27, 28, 29 } }; 

void DivideDiffBin(TH1* h1, TH1* h2){
    int nb1 = h1->GetNbinsX();
    int nb2 = h2->GetNbinsX();

    if(nb1 >= nb2) {
        for(int ib=1; ib<=nb1; ib++){
            double x  = h1->GetBinCenter(ib);
            double y1 = h1->GetBinContent(ib);
            double e1 = h1->GetBinError(ib);
            double y2 = h2->GetBinContent(h2->FindBin(x));
            double e2 = h2->GetBinError(h2->FindBin(x));
            y2!=0 ? h1->SetBinContent(ib,y1/y2) :h1->SetBinContent(ib,0);
            if(y1!=0 && y2!=0){
                h1->SetBinError(ib, y1/y2*sqrt(e1*e1/y1/y1 + e2*e2/y2/y2));
            }else{
                h1->SetBinError(ib,0);
            }
        }
    }

    if(nb1 < nb2) {
        for(int ib=1; ib<=nb2; ib++){
            double x  = h2->GetBinCenter(ib);
            double y2 = h2->GetBinContent(ib);
            double e2 = h2->GetBinError(ib);
            double y1 = h1->GetBinContent(h1->FindBin(x));
            double e1 = h1->GetBinError(h1->FindBin(x));
            y2!=0 ? h1->SetBinContent(ib,y1/y2) :h1->SetBinContent(ib,0);
            if(y1!=0 && y2!=0){
                h1->SetBinError(ib, y1/y2*sqrt(e1*e1/y1/y1 + e2*e2/y2/y2));
            }else{
                h1->SetBinError(ib,0);
            }
        }
    }
}

class MPlot

{
    private:
        // fCanvas 
        double fMarginLeft;     // left margin for the canvas
        double fMarginRight;    // right margin for the canvas
        double fMarginBottom;   // bottom margin for the canvas
        double fMarginTop;      // top margin for the canvas

        int fTopLeftX;          // Number of pixels the canvas is out of the top left corner of the screen in x-direction
        int fTopLeftY;          // Number of pixels the the canvas is out of the top left corner of the screen in y-direction
        int fCanvasWidth;       // Width of the canvas
        int fCanvasHeight;      // Height of the canvas

        // Appearance settings for histogram
        TString fXTitle;
        TString fYTitle;

        // TPad appearance settings
        int fLogX;              // 1: Logarithmic x-axis, 0: linear x-axis
        int fLogY;              // 1: Logarithmic y-axis, 0: linear y-axis
        int fLogZ;              // 1: Logarithmic z-axis, 0: linear z-axis
        int fGridX;             // 1: Grid in x-axis, 0: no grid
        int fGridY;             // 1: Grid in y-axis, 0: no grid
        int fTickX;             // 1: Ticks in x-axis, 0: no ticks
        int fTickY;             // 1: Ticks in y-axis, 0: no ticks

        // Canvas displacement settings
        int fCanvasDisplacementX;  // Displacement of the canvas for one index in x-direction
        int fCanvasDisplacementY;  // Displacement of the canvas for one index in y-direction
        int fCanvasesInOneRow;     // Number of canvases in a row before a new row is started

        // Handles needed for custom canvas creation
        bool isSplit; // automatic split canvas if this is true

        // Set palette (open or full symbols at the moment
        int fpalette; // 0: open, 1: full
        std::array<int,13> fMarkers;

        // Plotting histos and graphs stored in TObjArray
        TObjArray * histList;
        TObjArray * histRatioList;
        TObjArray * functList;
        TObjArray * graphList;

        // Automatic legend and base histo creation is also provided
        TLegend *l;
        TLegend *lNFO;
        TH2F * fHfr; TH2F * fHfr_ratio;

        // Settings for canvas splitting
        float fSplitRatio;     // The percentage of the lower pad from the total canvas when splitting the canvas

        // Index for generating names
        int fNameIndex;      // Names for canvases and pads are generated using this index. Given to the class at creation

        TString fOptHisto;
        TString fOptRatioHisto;
        TString fOptLegend;

    public:
        TCanvas * fCanvas;    
        TPad * fPad;              
        TPad * fPadRatio;

		MPlot()
		{
			Reset();
			fSplitRatio=0.4;
			Initialize(0,"","",false);
		}

        MPlot(int index, TString xlab, TString ylab, bool sp, float splitratio=0.4){
            Reset();
            fSplitRatio = splitratio;
            Initialize(index, xlab, ylab, sp);
            fMarkers = Markers_open; //open default, set later with SetPalette
        }
        virtual ~MPlot(){
            if (fCanvas) delete fCanvas;
            if (fPad) delete fPad;
            if (fPadRatio) delete fPadRatio;
            if (histList) delete histList;
            if (histRatioList) delete histRatioList;
        }

        /*
         * Initializing: plot index, x-y labels and base histogram, assuming autoscale. 
         * Special scale then be applied later if needed.
         */
        void Initialize(int index, TString xtit, TString ytit, bool split){
            fNameIndex = index; 
            isSplit = split;
            fXTitle = xtit;
            fYTitle = ytit;
            
            histList  = new TObjArray();
            functList = new TObjArray();
            graphList = new TObjArray();

            l = new TLegend(0.55, 0.65, 0.92, 0.92, "", "brNDC");
            l->SetFillStyle(0); l->SetBorderSize(0); l->SetTextSize(0.03);

            lNFO = new TLegend(0.0, 0.8, 0.92, 0.92, "", "brNDC");
            lNFO->SetFillStyle(0); lNFO->SetBorderSize(0); lNFO->SetTextSize(0.03);
            
            if(split) {
                InitializeSplitCanvas();
                histRatioList = new TObjArray();
            }
            else {
                fSplitRatio = 1;
                InitializeCanvas();
            }
        }

        /* 
         * Initialize a single JYU-style canvas (based on mc), and base histogram. 
         * Range is set to [0,1] now, proper range will be set while Draw
         */
        void InitializeCanvas() {
            TString cname = Form("c%d",fNameIndex);
            TString pname = Form("p%d",fNameIndex);
            // mc-style formatting, brute force now.
            float sdx = 300;
            float sdy = 100;
            int ncolumns = 5;
            //int ic0 = fNameIndex-1;
            //int cx = sdx*(ic0%ncolumns)+10;
            //int cy = sdy*(ic0-ic0%ncolumns)/ncolumns+10;
            fCanvas = new TCanvas(cname, cname, 440, 440); // should be formatted here
            fPad = new TPad(pname, pname, 0.01, 0.01, 0.99, 0.99, 0,0,0);
            SetPadRange(fPad, fMarginLeft, fMarginRight, fMarginTop, fMarginBottom);
            fHfr = new TH2F(Form("hfr%d", fNameIndex), "tit", 800, -1, 1, 800, 0, 1000);
            SetHfrStyle( fHfr, fXTitle, fYTitle, 1);
        }

        /* 
         * Initialize a single JYU-style canvas (based on mc), and base histogram. 
         * Range is set to [0,1] now, proper range will be set while Draw
         */
        void InitializeSplitCanvas() {
            TString cname = Form("c%d",fNameIndex);
            TString pname = Form("p%d",fNameIndex);
            fCanvas = new TCanvas(cname, cname, 462, 660); // should be formatted here

            fPad = new TPad(pname+"U", pname+"U", 0.0, fSplitRatio, 1.0, 1.0, 0,0,0);
            SetPadRange(fPad, fMarginLeft, fMarginRight, fMarginTop, 0.0);

            fPadRatio = new TPad(pname+"L", pname+"L", 0.0, 0.0, 1.0, fSplitRatio, 0,0,0);
            SetPadRange(fPadRatio, fMarginLeft, fMarginRight, 0.0, fMarginBottom/fSplitRatio*(1.0-fSplitRatio));

            fHfr = new TH2F(Form("hfr%d", fNameIndex), "tit", 800, 0, 1, 800, 0, 1);
            SetHfrStyle( fHfr, fXTitle, fYTitle,1);

            fHfr_ratio = new TH2F(Form("hfr%d_ratio", fNameIndex), "tit", 800, 0, 1, 800, 0, 1);
            SetHfrStyle( fHfr_ratio, fXTitle, "ratio", fSplitRatio/(1.0-fSplitRatio));
        }
        /* 
         * Sets palette (for open or full) symbols
         */
        void SetPalette(int ipalette){
            fpalette = ipalette;
            if(ipalette==0) fMarkers = Markers_open;
            if(ipalette==1) fMarkers = Markers_full;
        }

        /*
         * get combined range of histograms of histList
         */
        std::vector<double> GetHRanges( TObjArray* hList ) {
            std::vector<double> xmin;
            std::vector<double> xmax;
            std::vector<double> ymin;
            std::vector<double> ymax;
            int N = hList->GetEntries(); 
            for(int index=0; index<N; index++){
                xmin.push_back( ((TH1*) (hList->At(index)))->GetXaxis()->GetXmin() );
                xmax.push_back( ((TH1*) (hList->At(index)))->GetXaxis()->GetXmax() );
                ymin.push_back( ((TH1*) (hList->At(index)))->GetMinimum() );
                ymax.push_back( ((TH1*) (hList->At(index)))->GetMaximum() );
            }
            std::vector<double> hranges;
            hranges.push_back( get_min(xmin) );
            hranges.push_back( get_max(xmax) );
            hranges.push_back( get_min(ymin) );
            hranges.push_back( get_max(ymax) );

            return hranges;
        }

        /*
         * get combined y-range of histograms in a given x-range
         */
        std::vector<double> GetYRanges( TObjArray* hList, double xmin, double xmax )
        {
            std::vector<double> y, yranges;
            int N = hList->GetEntries();
            for(int ihist=0; ihist<N; ihist++)
            {
                int Nbins = ((TH1*)(hList->At(ihist)))->GetNbinsX();
                int ixmin = ((TH1*)(hList->At(ihist)))->FindBin(xmin);
                int ixmax = ((TH1*)(hList->At(ihist)))->FindBin(xmax);
                for(int ibin=ixmin; ibin<ixmax; ibin++){
                    y.push_back( ((TH1*)(hList->At(ihist)))->GetBinContent(ibin) );
                }
            }
            yranges.push_back( get_min(y));
            yranges.push_back( get_max(y));
            return yranges;
        }
                

        double get_min( std::vector<double> v ){
            double min = v[0];
            for( auto iv : v ) { if(iv < min) { min = iv; } }
            return min;
        }
        double get_max( std::vector<double> v ) {
            double max = v[0];
            for( auto iv : v ) { if(iv >= max) { max = iv; } }
            return max;
        }

        /*
         * Initialize list of histograms with the corresponding labels. 
         * This is needed for automatic range determination
         * Also initializes all ratio histograms
         */
        void addHList(std::vector<TH1*> hList, std::vector<TString> lList, TString opt_h="PE", TString opt_l="P", TString opt_hr="l") {
            fOptHisto = opt_h; fOptLegend = opt_l; fOptRatioHisto = opt_hr;
            for( unsigned i=0; i< hList.size(); i++) {
                histList->Add( hList.at(i) );
                SetAppearance( hList.at(i), Colors[i], Colors[i], fMarkers[i]);
                l->AddEntry( hList.at(i), lList.at(i), opt_l );
            }
            SetHfrRange( fHfr, histList, "auto" );
            if(isSplit){
                InitHRatio();
                SetHfrRange( fHfr_ratio, histRatioList, "auto" );
            }
        }

        /*
         * Clears list of histograms (and the ratio histograms to if they 
         * exist, if isSplit==True) with their corresponding legends *l. 
         * Useful to recycle the plotting class MPlot.
         */
        void resetHList() {
            histList->RemoveAll();
            if(isSplit){ histRatioList->RemoveAll(); }
            l->Clear();
        }

        /* 
         * Clear the info legend
         */
        void resetInfo() { lNFO->Clear(); }

        /*
         * Sets colors, marker styles, and marker size of given histogram
         */
        void SetAppearance( TH1 * h, int mcol, int lcol, int msty, double msiz = 1.0 ){
            h->SetMarkerColor( mcol );   
            h->SetMarkerStyle( msty );
            h->SetMarkerSize( msiz );
            h->SetLineColor( lcol );
        }

        /*
         * Draws all objects of the class, with the predefined options
         */ 
        void Draw() {
            gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
            fCanvas->Draw(); fPad->Draw(); fPad->cd(); fHfr->Draw();
            l->Draw(); lNFO->Draw();
            for( Int_t i=0; i<histList->GetEntriesFast();i++) 
                histList->At(i)->Draw(fOptHisto+"same");

            if(isSplit){
                fCanvas->cd(); fPadRatio->Draw(); fPadRatio->cd(); fHfr_ratio->Draw();
                for( Int_t i=0; i<histRatioList->GetEntriesFast();i++) { 
                    histRatioList->At(i)->Draw(fOptRatioHisto+"same");
                }
            }
        }

        /*
         * Draws given object to the already created top pad
         * (should work, but usually does not)
         */ 
        void DrawThis( TObject * obj, TString opt){
            fPad->cd();
            obj->Draw(opt);
        }
        /*
         * Draws given TF1 to the pad
         */
        void DrawThisTF1( TF1 * f, TString opt){
            fPad->cd();
            f->Draw(opt.Data());
        }
        /*
         * Draws and formats given graph to the already created top pad
         */
        void DrawThisTGraph( TGraph * g, TString opt, int marker, int color){
            g->SetLineColor(color);
            g->SetMarkerStyle(marker);
            g->Draw(opt);
        }

        void DrawThisTGraphAsymmErrors( TGraphAsymmErrors * g, TString opt, int marker, int color){
            fPad->cd();
            g->SetMarkerColor(color);
            g->SetLineColor(color);
            g->SetMarkerStyle(marker);
            g->Draw(opt);
        }


        void DrawThisTH1( TH1 * h, TString opt, int marker, int color, double alpha=1){
            fPad->cd();
            h->SetMarkerStyle(marker);
            if(alpha==1){ h->SetMarkerColor(color); }
            else { h->SetFillColorAlpha(color, alpha); }
            h->SetLineColor(color);
            h->Draw(opt);
        }

        void Save(TString filename ){
            fCanvas->SaveAs( filename+".pdf" );
        }

        void AddInfo( TString textNFO ) {
            lNFO->AddEntry((TObject*)0, textNFO, "");
        }
        void AddThisEntry( TObject * obj, TString text, TString opt="PE" ) {
            l->AddEntry( obj, text, opt );
        }

        void InitHRatio() {
            TH1D * hnom; 
            TH1D * hdenom;
            hdenom = (TH1D*)(histList->At(0))->Clone();
            for(Int_t i=1; i<histList->GetEntriesFast(); i++){
                hnom = (TH1D*)(histList->At(i))->Clone();
//                hnom->Divide(hdenom);
                DivideDiffBin( hnom, hdenom);
                histRatioList->Add( hnom );
            }
        }
        void SetRatioLimits(double ymin, double ymax){
            SetHfrRange( fHfr_ratio, histRatioList, "", ymin, ymax);
        }
        void SetRatioLimitsXY(double xmin, double xmax, double ymin, double ymax){
            SetHfrRange( fHfr_ratio, histRatioList, "", ymin, ymax, xmin, xmax, "");
        }
        void SetLimitsX( double xmin, double xmax ){
            SetHfrRange( fHfr, histList, "auto", 0, 1, xmin, xmax, "");
            if(isSplit) SetHfrRange( fHfr_ratio, histRatioList, "auto", 0, 1, xmin, xmax, "");
            fPad->Update();
        }
        void SetLimitsY( double ymin, double ymax ){
            SetHfrRange( fHfr, histList, "", ymin, ymax);
        }
        void SetLimitsXY( double xmin, double xmax, double ymin, double ymax ){
            SetHfrRange( fHfr, histList, "", ymin, ymax, xmin, xmax, "");
            if(isSplit) SetHfrRange( fHfr_ratio, histRatioList, "auto", 0, 1, xmin, xmax, "");
        }
        void SetHfrRange( TH1*h, TObjArray * hList, TString auto_y, double ymin=0, double ymax=1,double xmin=0, double xmax=1, TString auto_x="auto"){
           if( auto_x == "auto"){
               std::vector<double> ranges;
               ranges = GetHRanges(hList);
               double xmin = ranges[0]; double xmax = ranges[1];
               h->GetXaxis()->SetLimits( xmin, xmax );
               h->GetYaxis()->SetLimits( ymin, ymax );
            }
           if( auto_y == "auto"){
               std::vector<double> ranges;
               ranges = GetYRanges(hList, xmin, xmax);
               double yw_low=0.1; double yw_high=0.2;

               double height = ranges[1]-ranges[0];
               double ymin = ranges[0]-height*yw_low;
               double ymax = ranges[1]+height*yw_high; 
               h->GetYaxis()->SetLimits( ymin, ymax );
            }
            if( auto_x!="auto"){
               h->GetXaxis()->SetLimits( xmin, xmax );
            }
            if( auto_y!="auto"){
               h->GetYaxis()->SetLimits( ymin, ymax );
            }
        }
        void SetPadRange(TPad * pad, double r1, double r2, double r3, double r4) {
            pad->SetLeftMargin(r1); pad->SetRightMargin(r2);
            pad->SetTopMargin(r3);  pad->SetBottomMargin(r4);
        }
        void SetRatioLabel(TString ytit){
            fHfr_ratio->GetYaxis()->SetTitle(ytit);
        }
        void SetLog(bool logx, bool logy){
            fPad->SetLogx(logx);
            fPad->SetLogy(logy);
            if(isSplit && logx){ fPadRatio->SetLogx( logx ); }
        }
        void SetGridY() { fPad->SetGridy(1); }
        void SetGridX() { 
            fPad->SetGridx(1);
            if(isSplit) fPadRatio->SetGridx(1);
        }
        void SetGridXY() { 
            fPad->SetGrid(1, 1);
            if(isSplit) fPadRatio->SetGrid(1, 1);
        }

        void SetHfrStyle( TH1 * histo, TString xtit, TString ytit, float ratio ){
            double fTitleOffsetX = 1.1;    // Offset of the x-axis title
            double fTitleOffsetY = 1.1;    // Offset of the y-axis title
            double fTitleSizeX = 0.06;     // Size of the x-axis title
            double fTitleSizeY = 0.06;     // Size of the y-axis title
            double fLabelOffsetX = 0.01;   // Offset of the x-axis label
            double fLabelOffsetY = 0.001;  // Offset of the y-axis label
            double fLabelSizeX = 0.04;     // Size of the x-axis label
            double fLabelSizeY = 0.04;     // Size of the y-axis label
            int fDivisionsX = 505;      // The number of divisions in x-axis
            int fDivisionsY = 505;      // The number of divisions in y-axis
            int fFont = 42;             // Font index for titles and labels
            
            histo->GetXaxis()->CenterTitle(1);    // Axis titles are centered
            histo->GetYaxis()->CenterTitle(1);    // Axis titles are centered

            histo->GetXaxis()->SetTitleOffset(fTitleOffsetX); // Give a small offset to the title so that it does overlap with axis
            histo->GetYaxis()->SetTitleOffset(fTitleOffsetY*ratio); // Give a small offset to the title so that it does overlap with axis

            histo->GetXaxis()->SetTitleSize(fTitleSizeX/ratio); // Define the sixe of the title
            histo->GetYaxis()->SetTitleSize(fTitleSizeY/ratio); // Define the sixe of the title

            histo->GetXaxis()->SetLabelOffset(fLabelOffsetX/ratio); // Give a small offset to the label so that it does overlap with axis
            histo->GetYaxis()->SetLabelOffset(fLabelOffsetY/ratio); // Give a small offset to the label so that it does overlap with axis

            histo->GetXaxis()->SetLabelSize(fLabelSizeX/ratio); // Define the sixe of the label
            histo->GetYaxis()->SetLabelSize(fLabelSizeY/ratio); // Define the sixe of the label

            histo->GetXaxis()->SetNdivisions(fDivisionsX); // Set the number of division markers
            histo->GetYaxis()->SetNdivisions(fDivisionsY); // Set the number of division markers

            histo->GetXaxis()->SetTitle(xtit); // Set the axis title
            histo->GetYaxis()->SetTitle(ytit); // Set the axis title

            histo->GetXaxis()->SetLabelFont(fFont); // Set the label font
            histo->GetYaxis()->SetLabelFont(fFont); // Set the label font
            histo->GetXaxis()->SetTitleFont(fFont); // Set the title font
            histo->GetYaxis()->SetTitleFont(fFont); // Set the title font
        }

        void PrintHlist(){
            std::cout << "List of histograms:" << std::endl;
            for(Int_t i=0; i<histList->GetEntriesFast(); i++) histList->At(i)->Print() ;
            if(isSplit){
                std::cout << "List of ratio histograms:" << std::endl;
                for(Int_t i=0; i<histRatioList->GetEntriesFast(); i++) histRatioList->At(i)->Print() ;
            }
        }

        void Reset(){

            // Set default values for margins
            fMarginLeft = 0.15;     // left margin
            fMarginRight = 0.06;    // right margin
            fMarginBottom = 0.15;   // bottom margin
            fMarginTop = 0.06;      // top margin

            // Set default values for the size and place of the canvas
            fTopLeftX = 10;         // Number of pixels the the canvas is out of the top left corner of the screen in x-direction
            fTopLeftY = 10;         // Number of pixels the the canvas is out of the top left corner of the screen in y-direction
            fCanvasWidth = 700;     // Width of the canvas
            fCanvasHeight = 500;    // Height of the canvas

            // Set default values for histogram appearance settings

            // Axes are linear by default
            fLogX = 0;              // 1: Logarithmic x-axis, 0: linear x-axis
            fLogY = 0;              // 1: Logarithmic y-axis, 0: linear y-axis
            fLogZ = 0;              // 1: Logarithmic z-axis, 0: linear z-axis

            // Set default values for canves displacement information
            fCanvasDisplacementX = 100;  // Displacement of the canvas for one index in x-direction
            fCanvasDisplacementY = 100;  // Displacement of the canvas for one index in y-direction
            fCanvasesInOneRow = 10;      // Number of canvases in a row before a new row is started

            // Set defauls values for pad properties
            fGridX = 0;             // 1: Grid in x-axis, 0: no grid
            fGridY = 0;             // 1: Grid in y-axis, 0: no grid
            fTickX = 1;             // 1: Ticks in x-axis, 0: no ticks
            fTickY = 1;             // 1: Ticks in y-axis, 0: no ticks

            // Default setting values for splitted canvas
            fSplitRatio = 0.4;      // The percentage of the lower pad from the total canvas when splitting the canvas


            // Set the canvas handles to null
//            fCanvas = NULL;
//            fUpperSplitPad = NULL;
//            fLowerSplitPad = NULL;
        }

		//ClassDef(MPlot,1)
};

#endif /* MPLOT_H */

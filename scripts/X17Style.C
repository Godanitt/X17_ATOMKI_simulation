#ifndef X17STYLE_C
#define X17STYLE_C

#include <TROOT.h>
#include <TStyle.h>
#include <TGaxis.h>

void SetX17Style()
{
    TStyle *x17Style = new TStyle("X17Style", "X17 publication style");

    // ----------------------------
    // Canvas
    // ----------------------------
    x17Style->SetCanvasBorderMode(0);
    x17Style->SetCanvasColor(kWhite);
    x17Style->SetCanvasDefH(700);
    x17Style->SetCanvasDefW(800);
    x17Style->SetCanvasDefX(0);
    x17Style->SetCanvasDefY(0);

    // ----------------------------
    // Pads
    // ----------------------------
    x17Style->SetPadBorderMode(0);
    x17Style->SetPadColor(kWhite);
    x17Style->SetPadGridX(false);
    x17Style->SetPadGridY(false);
    x17Style->SetGridColor(0);
    x17Style->SetGridStyle(3);
    x17Style->SetGridWidth(1);

    // ----------------------------
    // Frame
    // ----------------------------
    x17Style->SetFrameBorderMode(0);
    x17Style->SetFrameBorderSize(1);
    x17Style->SetFrameFillColor(kWhite);
    x17Style->SetFrameFillStyle(0);
    x17Style->SetFrameLineColor(kBlack);
    x17Style->SetFrameLineStyle(1);
    x17Style->SetFrameLineWidth(2);

    // ----------------------------
    // Histograms
    // ----------------------------
    x17Style->SetHistLineColor(kBlack);
    x17Style->SetHistLineStyle(1);
    x17Style->SetHistLineWidth(2);
    x17Style->SetEndErrorSize(0);
    x17Style->SetMarkerStyle(20);
    x17Style->SetMarkerSize(1.0);

    // ----------------------------
    // Fits
    // ----------------------------
    x17Style->SetFuncColor(kRed + 1);
    x17Style->SetFuncStyle(1);
    x17Style->SetFuncWidth(2);

    // ----------------------------
    // Statistics box
    // ----------------------------
    x17Style->SetOptStat(0);
    x17Style->SetOptFit(0);

    // ----------------------------
    // Margins
    // ----------------------------
    x17Style->SetPadTopMargin(0.07);
    x17Style->SetPadRightMargin(0.04);
    x17Style->SetPadBottomMargin(0.15);
    x17Style->SetPadLeftMargin(0.16);

    // ----------------------------
    // Titles
    // ----------------------------
    x17Style->SetOptTitle(0);
    x17Style->SetTitleFont(42, "XYZ");
    x17Style->SetTitleSize(0.055, "XYZ");
    x17Style->SetTitleXOffset(1.10);
    x17Style->SetTitleYOffset(1.25);

    // ----------------------------
    // Axis labels
    // ----------------------------
    x17Style->SetLabelFont(42, "XYZ");
    x17Style->SetLabelSize(0.045, "XYZ");
    x17Style->SetLabelOffset(0.010, "XYZ");

    x17Style->SetAxisColor(kBlack, "XYZ");
    x17Style->SetStripDecimals(kTRUE);
    x17Style->SetTickLength(0.03, "XYZ");
    x17Style->SetNdivisions(510, "XYZ");

    // Ticks on all sides
    x17Style->SetPadTickX(1);
    x17Style->SetPadTickY(1);

    // ----------------------------
    // Legend
    // ----------------------------
    x17Style->SetLegendBorderSize(0);
    x17Style->SetLegendFillColor(kWhite);
    x17Style->SetLegendFont(42);
    x17Style->SetLegendTextSize(0.040);

    // ----------------------------
    // Text
    // ----------------------------
    x17Style->SetTextFont(42);
    x17Style->SetTextSize(0.045);

    // ----------------------------
    // Palette, useful later for 2D
    // ----------------------------
    x17Style->SetPalette(kViridis);

    // ----------------------------
    // Apply
    // ----------------------------
    gROOT->SetStyle("X17Style");
    gROOT->ForceStyle();

    TGaxis::SetMaxDigits(3);
}

#endif
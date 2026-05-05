// ROOT macro to plot the angular distributions used by the X17-like generator.
//
// Usage from the project root:
//   root -l -q 'scripts/plot_angles_ROOT.C("data/data_pair_creation.txt", "plots")'
//
// It creates PDF files only:
//   plots/thetaee_distribution_ROOT.pdf
//   plots/angular_distributions_ROOT.pdf

#include "X17Style.C"

#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TString.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <initializer_list>

void plot_angles_ROOT(const char* infile = "data/data_pair_creation.txt",
                      const char* outdir = "plots")
{
    // Apply common paper-style ROOT configuration
    SetX17Style();

    // Create output directory if needed
    gSystem->mkdir(outdir, kTRUE);

    std::ifstream fin(infile);
    if (!fin.is_open()) {
        std::cerr << "[plot_angles_ROOT] ERROR: no se pudo abrir "
                  << infile << std::endl;
        return;
    }

    double theta_ee, theta_em, T_em, theta_ep, T_ep;

    const int nbins = 90;
    const double xmin = 0.0;
    const double xmax = 180.0;

    TH1D *hThetaEE = new TH1D("hThetaEE", "", nbins, xmin, xmax);
    TH1D *hThetaEm = new TH1D("hThetaEm", "", nbins, xmin, xmax);
    TH1D *hThetaEp = new TH1D("hThetaEp", "", nbins, xmin, xmax);

    hThetaEE->Sumw2();
    hThetaEm->Sumw2();
    hThetaEp->Sumw2();

    long nrows = 0;

    while (fin >> theta_ee >> theta_em >> T_em >> theta_ep >> T_ep) {
        hThetaEE->Fill(theta_ee);
        hThetaEm->Fill(theta_em);
        hThetaEp->Fill(theta_ep);
        nrows++;
    }

    std::cout << "[plot_angles_ROOT] Read " << nrows
              << " rows from " << infile << std::endl;

    // ----------------------------
    // Histogram style
    // ----------------------------
    hThetaEE->SetLineColor(kBlack);
    hThetaEE->SetLineWidth(3);

    hThetaEm->SetLineColor(kBlue + 1);
    hThetaEm->SetLineWidth(2);

    hThetaEp->SetLineColor(kRed + 1);
    hThetaEp->SetLineWidth(2);

    hThetaEE->GetXaxis()->SetTitle("#theta [deg]");
    hThetaEE->GetYaxis()->SetTitle("Entries / 2 deg");
    hThetaEE->GetXaxis()->SetRangeUser(xmin, xmax);

    hThetaEm->GetXaxis()->SetTitle("#theta [deg]");
    hThetaEm->GetYaxis()->SetTitle("Entries / 2 deg");
    hThetaEm->GetXaxis()->SetRangeUser(xmin, xmax);

    hThetaEp->GetXaxis()->SetTitle("#theta [deg]");
    hThetaEp->GetYaxis()->SetTitle("Entries / 2 deg");
    hThetaEp->GetXaxis()->SetRangeUser(xmin, xmax);

    // ----------------------------
    // Figure 1: theta_ee only
    // ----------------------------
    TCanvas *c1 = new TCanvas("c1", "thetaee", 800, 700);
    c1->cd();

    hThetaEE->SetMaximum(1.20 * hThetaEE->GetMaximum());
    hThetaEE->Draw("HIST");

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.045);
    label.DrawLatex(0.20, 0.86, "X17-like generator");
    label.SetTextSize(0.038);
    label.DrawLatex(0.20, 0.80, "Input angular distribution");

    c1->SaveAs(Form("%s/thetaee_distribution_ROOT.pdf", outdir));

    // ----------------------------
    // Figure 2: all angular distributions
    // ----------------------------
    TCanvas *c2 = new TCanvas("c2", "all_angles", 800, 700);
    c2->cd();

    double maxy = std::max({
        hThetaEE->GetMaximum(),
        hThetaEm->GetMaximum(),
        hThetaEp->GetMaximum()
    });

    hThetaEE->SetMaximum(1.25 * maxy);
    hThetaEE->Draw("HIST");
    hThetaEm->Draw("HIST SAME");
    hThetaEp->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.62, 0.70, 0.88, 0.88);
    leg->AddEntry(hThetaEE, "#theta_{ee}", "l");
    leg->AddEntry(hThetaEm, "#theta_{e^{-}}", "l");
    leg->AddEntry(hThetaEp, "#theta_{e^{+}}", "l");
    leg->Draw();

    TLatex label2;
    label2.SetNDC();
    label2.SetTextFont(42);
    label2.SetTextSize(0.045);
    label2.DrawLatex(0.20, 0.86, "X17-like generator");

    c2->SaveAs(Form("%s/angular_distributions_ROOT.pdf", outdir));

    std::cout << "[plot_angles_ROOT] Saved:\n"
              << "  " << outdir << "/thetaee_distribution_ROOT.pdf\n"
              << "  " << outdir << "/angular_distributions_ROOT.pdf"
              << std::endl;
}
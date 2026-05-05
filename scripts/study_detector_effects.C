// Detector-effects study for the X17-like simulation.
//
// Default usage:
//   root -l -b -q -I scripts 'scripts/study_detector_effects.C'
//
// Input:
//   x17_output.root
//
// Output ROOT ntuple:
//   x17_detector_effects.root
//
// Output PDFs:
//   plots/detector_acceptance_vs_theta.pdf
//   plots/detector_thetaee_distortion.pdf
//   plots/detector_thetaee_resolution.pdf
//   plots/detector_threshold_scan.pdf

#include "X17Style.C"

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TMath.h>

#include <iostream>
#include <cmath>
#include <algorithm>
#include <initializer_list>
#include <vector>

namespace
{
    double OpeningAngleDeg(double x1, double y1, double z1,
                           double x2, double y2, double z2)
    {
        const double r1 = std::sqrt(x1*x1 + y1*y1 + z1*z1);
        const double r2 = std::sqrt(x2*x2 + y2*y2 + z2*z2);

        if (r1 <= 0.0 || r2 <= 0.0) return -999.0;

        double c = (x1*x2 + y1*y2 + z1*z2) / (r1 * r2);
        c = std::max(-1.0, std::min(1.0, c));

        return std::acos(c) * 180.0 / TMath::Pi();
    }

    bool finite3(double x, double y, double z)
    {
        return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
    }

    double safe_eff_err(double passed, double total)
    {
        if (total <= 0.0) return 0.0;
        const double e = passed / total;
        return std::sqrt(e * (1.0 - e) / total);
    }
}

void study_detector_effects(const char* inputFile = "x17_output.root",
                            const char* outdir = "plots")
{
    SetX17Style();
    gSystem->mkdir(outdir, kTRUE);

    TFile* fin = TFile::Open(inputFile, "READ");

    if (!fin || fin->IsZombie())
    {
        std::cerr << "[study_detector_effects] ERROR: cannot open "
                  << inputFile << std::endl;
        return;
    }

    TTree* events = nullptr;
    fin->GetObject("events", events);

    if (!events)
    {
        std::cerr << "[study_detector_effects] ERROR: tree 'events' not found in "
                  << inputFile << std::endl;
        fin->ls();
        return;
    }

    // ------------------------------------------------------------
    // Input branches
    // ------------------------------------------------------------
    Int_t eventID = -1;

    Double_t thetaEE_gen_deg = -999.0;
    Double_t thetaEm_gen_deg = -999.0;
    Double_t thetaEp_gen_deg = -999.0;

    Double_t kinEm_gen_MeV = -999.0;
    Double_t kinEp_gen_MeV = -999.0;

    Int_t hasEmHit = 0;
    Int_t hasEpHit = 0;
    Int_t detEm = -1;
    Int_t detEp = -1;

    Double_t xEm_mm = 0.0;
    Double_t yEm_mm = 0.0;
    Double_t zEm_mm = 0.0;

    Double_t xEp_mm = 0.0;
    Double_t yEp_mm = 0.0;
    Double_t zEp_mm = 0.0;

    Double_t thetaEE_reco_stored_deg = -999.0;

    Double_t edepEm_MeV = 0.0;
    Double_t edepEp_MeV = 0.0;
    Double_t edepAll_MeV = 0.0;

    Int_t pairAccepted = 0;
    Int_t nHitRows = 0;

    events->SetBranchAddress("eventID", &eventID);

    events->SetBranchAddress("thetaEE_gen_deg", &thetaEE_gen_deg);
    events->SetBranchAddress("thetaEm_gen_deg", &thetaEm_gen_deg);
    events->SetBranchAddress("thetaEp_gen_deg", &thetaEp_gen_deg);

    events->SetBranchAddress("kinEm_gen_MeV", &kinEm_gen_MeV);
    events->SetBranchAddress("kinEp_gen_MeV", &kinEp_gen_MeV);

    events->SetBranchAddress("hasEmHit", &hasEmHit);
    events->SetBranchAddress("hasEpHit", &hasEpHit);

    events->SetBranchAddress("detEm", &detEm);
    events->SetBranchAddress("detEp", &detEp);

    events->SetBranchAddress("xEm_mm", &xEm_mm);
    events->SetBranchAddress("yEm_mm", &yEm_mm);
    events->SetBranchAddress("zEm_mm", &zEm_mm);

    events->SetBranchAddress("xEp_mm", &xEp_mm);
    events->SetBranchAddress("yEp_mm", &yEp_mm);
    events->SetBranchAddress("zEp_mm", &zEp_mm);

    events->SetBranchAddress("thetaEE_reco_deg", &thetaEE_reco_stored_deg);

    events->SetBranchAddress("edepEm_MeV", &edepEm_MeV);
    events->SetBranchAddress("edepEp_MeV", &edepEp_MeV);
    events->SetBranchAddress("edepAll_MeV", &edepAll_MeV);

    events->SetBranchAddress("pairAccepted", &pairAccepted);
    events->SetBranchAddress("nHitRows", &nHitRows);

    // ------------------------------------------------------------
    // Output ntuple only
    // ------------------------------------------------------------
    TFile* fout = TFile::Open("x17_detector_effects.root", "RECREATE");
    TTree* effects = new TTree("detector_effects",
                               "event-level detector-effects study");

    Double_t thetaEE_reco_deg = -999.0;
    Double_t thetaEE_smeared_deg = -999.0;
    Double_t thetaEE_res_deg = -999.0;
    Double_t thetaEE_smear_res_deg = -999.0;

    Int_t pass_hits = 0;
    Int_t pass_position = 0;
    Int_t pass_edep_005 = 0;
    Int_t pass_edep_010 = 0;
    Int_t pass_edep_020 = 0;
    Int_t pass_edep_050 = 0;
    Int_t pass_basic = 0;

    effects->Branch("eventID", &eventID, "eventID/I");

    effects->Branch("thetaEE_gen_deg", &thetaEE_gen_deg, "thetaEE_gen_deg/D");
    effects->Branch("thetaEE_reco_deg", &thetaEE_reco_deg, "thetaEE_reco_deg/D");
    effects->Branch("thetaEE_smeared_deg", &thetaEE_smeared_deg, "thetaEE_smeared_deg/D");

    effects->Branch("thetaEE_res_deg", &thetaEE_res_deg, "thetaEE_res_deg/D");
    effects->Branch("thetaEE_smear_res_deg", &thetaEE_smear_res_deg, "thetaEE_smear_res_deg/D");

    effects->Branch("thetaEm_gen_deg", &thetaEm_gen_deg, "thetaEm_gen_deg/D");
    effects->Branch("thetaEp_gen_deg", &thetaEp_gen_deg, "thetaEp_gen_deg/D");

    effects->Branch("kinEm_gen_MeV", &kinEm_gen_MeV, "kinEm_gen_MeV/D");
    effects->Branch("kinEp_gen_MeV", &kinEp_gen_MeV, "kinEp_gen_MeV/D");

    effects->Branch("detEm", &detEm, "detEm/I");
    effects->Branch("detEp", &detEp, "detEp/I");

    effects->Branch("edepEm_MeV", &edepEm_MeV, "edepEm_MeV/D");
    effects->Branch("edepEp_MeV", &edepEp_MeV, "edepEp_MeV/D");
    effects->Branch("edepAll_MeV", &edepAll_MeV, "edepAll_MeV/D");

    effects->Branch("pass_hits", &pass_hits, "pass_hits/I");
    effects->Branch("pass_position", &pass_position, "pass_position/I");
    effects->Branch("pass_edep_005", &pass_edep_005, "pass_edep_005/I");
    effects->Branch("pass_edep_010", &pass_edep_010, "pass_edep_010/I");
    effects->Branch("pass_edep_020", &pass_edep_020, "pass_edep_020/I");
    effects->Branch("pass_edep_050", &pass_edep_050, "pass_edep_050/I");
    effects->Branch("pass_basic", &pass_basic, "pass_basic/I");

    // ------------------------------------------------------------
    // Histograms for PDFs only
    // ------------------------------------------------------------
    TH1D* hGenAll = new TH1D("hGenAll", "", 90, 0.0, 180.0);
    TH1D* hGenAccepted = new TH1D("hGenAccepted", "", 90, 0.0, 180.0);
    TH1D* hRecoAccepted = new TH1D("hRecoAccepted", "", 90, 0.0, 180.0);
    TH1D* hRecoSmeared = new TH1D("hRecoSmeared", "", 90, 0.0, 180.0);

    TH1D* hResolution = new TH1D("hResolution", "", 100, -50.0, 50.0);
    TH1D* hResolutionSmeared = new TH1D("hResolutionSmeared", "", 100, -50.0, 50.0);

    TH1D* hAccNum = new TH1D("hAccNum", "", 18, 0.0, 180.0);
    TH1D* hAccDen = new TH1D("hAccDen", "", 18, 0.0, 180.0);

    hGenAll->SetLineColor(kGray + 2);
    hGenAll->SetLineWidth(2);

    hGenAccepted->SetLineColor(kBlack);
    hGenAccepted->SetLineWidth(3);

    hRecoAccepted->SetLineColor(kRed + 1);
    hRecoAccepted->SetLineWidth(3);

    hRecoSmeared->SetLineColor(kBlue + 1);
    hRecoSmeared->SetLineWidth(3);

    hResolution->SetLineColor(kBlack);
    hResolution->SetLineWidth(3);

    hResolutionSmeared->SetLineColor(kBlue + 1);
    hResolutionSmeared->SetLineWidth(3);

    // ------------------------------------------------------------
    // Smearing model
    //
    // It represents finite hit-position resolution.
    // For silicon-strip-like information, 1 mm is a simple conservative default.
    // ------------------------------------------------------------
    TRandom3 rng(12345);
    const double sigma_pos_mm = 1.0;

    Long64_t nAll = 0;
    Long64_t nAccepted = 0;

    Long64_t nThr005 = 0;
    Long64_t nThr010 = 0;
    Long64_t nThr020 = 0;
    Long64_t nThr050 = 0;

    const Long64_t nentries = events->GetEntries();

    for (Long64_t i = 0; i < nentries; ++i)
    {
        events->GetEntry(i);

        nAll++;
        hGenAll->Fill(thetaEE_gen_deg);
        hAccDen->Fill(thetaEE_gen_deg);

        thetaEE_reco_deg = OpeningAngleDeg(xEm_mm, yEm_mm, zEm_mm,
                                           xEp_mm, yEp_mm, zEp_mm);

        const double xEm_s = xEm_mm + rng.Gaus(0.0, sigma_pos_mm);
        const double yEm_s = yEm_mm + rng.Gaus(0.0, sigma_pos_mm);
        const double zEm_s = zEm_mm + rng.Gaus(0.0, sigma_pos_mm);

        const double xEp_s = xEp_mm + rng.Gaus(0.0, sigma_pos_mm);
        const double yEp_s = yEp_mm + rng.Gaus(0.0, sigma_pos_mm);
        const double zEp_s = zEp_mm + rng.Gaus(0.0, sigma_pos_mm);

        thetaEE_smeared_deg = OpeningAngleDeg(xEm_s, yEm_s, zEm_s,
                                              xEp_s, yEp_s, zEp_s);

        thetaEE_res_deg = thetaEE_reco_deg - thetaEE_gen_deg;
        thetaEE_smear_res_deg = thetaEE_smeared_deg - thetaEE_gen_deg;

        pass_hits = (hasEmHit == 1 && hasEpHit == 1 && pairAccepted == 1);

        const double rEm = std::sqrt(xEm_mm*xEm_mm + yEm_mm*yEm_mm + zEm_mm*zEm_mm);
        const double rEp = std::sqrt(xEp_mm*xEp_mm + yEp_mm*yEp_mm + zEp_mm*zEp_mm);

        pass_position =
            finite3(xEm_mm, yEm_mm, zEm_mm) &&
            finite3(xEp_mm, yEp_mm, zEp_mm) &&
            rEm > 1.0 &&
            rEp > 1.0 &&
            thetaEE_reco_deg >= 0.0 &&
            thetaEE_reco_deg <= 180.0;

        pass_edep_005 = (edepEm_MeV > 0.05 && edepEp_MeV > 0.05);
        pass_edep_010 = (edepEm_MeV > 0.10 && edepEp_MeV > 0.10);
        pass_edep_020 = (edepEm_MeV > 0.20 && edepEp_MeV > 0.20);
        pass_edep_050 = (edepEm_MeV > 0.50 && edepEp_MeV > 0.50);

        pass_basic = (pass_hits && pass_position && pass_edep_005);

        if (pass_hits && pass_position)
        {
            if (pass_edep_005) nThr005++;
            if (pass_edep_010) nThr010++;
            if (pass_edep_020) nThr020++;
            if (pass_edep_050) nThr050++;
        }

        if (pass_basic)
        {
            nAccepted++;

            hGenAccepted->Fill(thetaEE_gen_deg);
            hRecoAccepted->Fill(thetaEE_reco_deg);
            hRecoSmeared->Fill(thetaEE_smeared_deg);

            hResolution->Fill(thetaEE_res_deg);
            hResolutionSmeared->Fill(thetaEE_smear_res_deg);

            hAccNum->Fill(thetaEE_gen_deg);
        }

        effects->Fill();
    }

    fout->cd();
    effects->Write();
    fout->Close();

    std::cout << "\n[study_detector_effects] Summary\n";
    std::cout << "  All generated events      : " << nAll << "\n";
    std::cout << "  Basic accepted events     : " << nAccepted << "\n";

    if (nAll > 0)
    {
        std::cout << "  Global acceptance         : "
                  << 100.0 * double(nAccepted) / double(nAll)
                  << " %\n";
    }

    if (nAccepted > 0)
    {
        std::cout << "  Mean reco residual        : "
                  << hResolution->GetMean()
                  << " deg\n";
        std::cout << "  RMS reco residual         : "
                  << hResolution->GetRMS()
                  << " deg\n";
        std::cout << "  RMS smeared residual      : "
                  << hResolutionSmeared->GetRMS()
                  << " deg\n";
    }

    std::cout << "  Wrote ntuple              : x17_detector_effects.root\n";

    // ------------------------------------------------------------
    // Plot 1: acceptance vs theta_ee
    // ------------------------------------------------------------
    TGraphErrors* gAcc = new TGraphErrors();

    int ip = 0;

    for (int b = 1; b <= hAccDen->GetNbinsX(); ++b)
    {
        const double den = hAccDen->GetBinContent(b);
        const double num = hAccNum->GetBinContent(b);

        if (den <= 0.0) continue;

        const double x = hAccDen->GetBinCenter(b);
        const double ex = 0.5 * hAccDen->GetBinWidth(b);
        const double y = num / den;
        const double ey = safe_eff_err(num, den);

        gAcc->SetPoint(ip, x, y);
        gAcc->SetPointError(ip, ex, ey);
        ip++;
    }

    TCanvas* c1 = new TCanvas("c1", "acceptance", 800, 700);

    TH1D* frameAcc = new TH1D("frameAcc", "", 18, 0.0, 180.0);
    frameAcc->SetMinimum(0.0);
    frameAcc->SetMaximum(1.05);
    frameAcc->GetXaxis()->SetTitle("#theta_{ee}^{gen} [deg]");
    frameAcc->GetYaxis()->SetTitle("Acceptance");
    frameAcc->Draw("AXIS");

    gAcc->SetMarkerStyle(20);
    gAcc->SetMarkerSize(1.1);
    gAcc->SetLineWidth(2);
    gAcc->Draw("P SAME");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.20, 0.86, "Detector acceptance");
    latex.SetTextSize(0.036);
    latex.DrawLatex(0.20, 0.80, "Basic quality cuts");

    c1->SaveAs(Form("%s/detector_acceptance_vs_theta.pdf", outdir));

    // ------------------------------------------------------------
    // Plot 2: angular distortion
    // ------------------------------------------------------------
    TCanvas* c2 = new TCanvas("c2", "distortion", 800, 700);

    hGenAll->GetXaxis()->SetTitle("#theta_{ee} [deg]");
    hGenAll->GetYaxis()->SetTitle("Entries / 2 deg");
    hGenAll->SetMaximum(1.35 * std::max({
        hGenAll->GetMaximum(),
        hGenAccepted->GetMaximum(),
        hRecoAccepted->GetMaximum(),
        hRecoSmeared->GetMaximum()
    }));

    hGenAll->Draw("HIST");
    hGenAccepted->Draw("HIST SAME");
    hRecoAccepted->Draw("HIST SAME");
    hRecoSmeared->Draw("HIST SAME");

    TLegend* leg2 = new TLegend(0.50, 0.66, 0.88, 0.88);
    leg2->AddEntry(hGenAll, "Generated, all", "l");
    leg2->AddEntry(hGenAccepted, "Generated, accepted", "l");
    leg2->AddEntry(hRecoAccepted, "Reco, accepted", "l");
    leg2->AddEntry(hRecoSmeared, "Reco + 1 mm smearing", "l");
    leg2->Draw();

    latex.SetTextSize(0.045);
    latex.DrawLatex(0.20, 0.86, "Detector distortion");

    c2->SaveAs(Form("%s/detector_thetaee_distortion.pdf", outdir));

    // ------------------------------------------------------------
    // Plot 3: resolution
    // ------------------------------------------------------------
    TCanvas* c3 = new TCanvas("c3", "resolution", 800, 700);

    hResolution->GetXaxis()->SetTitle("#theta_{ee}^{reco} - #theta_{ee}^{gen} [deg]");
    hResolution->GetYaxis()->SetTitle("Entries");

    hResolution->SetMaximum(1.30 * std::max(hResolution->GetMaximum(),
                                             hResolutionSmeared->GetMaximum()));
    hResolution->Draw("HIST");
    hResolutionSmeared->Draw("HIST SAME");

    TLine* zero = new TLine(0.0, 0.0, 0.0, hResolution->GetMaximum());
    zero->SetLineStyle(2);
    zero->SetLineWidth(2);
    zero->Draw("SAME");

    TLegend* leg3 = new TLegend(0.56, 0.75, 0.88, 0.88);
    leg3->AddEntry(hResolution, "Reco", "l");
    leg3->AddEntry(hResolutionSmeared, "Reco + 1 mm smearing", "l");
    leg3->Draw();

    latex.SetTextSize(0.045);
    latex.DrawLatex(0.20, 0.86, "Angular resolution");

    c3->SaveAs(Form("%s/detector_thetaee_resolution.pdf", outdir));

    // ------------------------------------------------------------
    // Plot 4: threshold scan
    // ------------------------------------------------------------
    const std::vector<double> thresholds = {0.05, 0.10, 0.20, 0.50};
    const std::vector<double> passed = {
        double(nThr005),
        double(nThr010),
        double(nThr020),
        double(nThr050)
    };

    TGraphErrors* gThr = new TGraphErrors();

    for (int i = 0; i < int(thresholds.size()); ++i)
    {
        const double eff = (nAll > 0) ? passed[i] / double(nAll) : 0.0;
        const double err = safe_eff_err(passed[i], double(nAll));

        gThr->SetPoint(i, thresholds[i], eff);
        gThr->SetPointError(i, 0.0, err);
    }

    TCanvas* c4 = new TCanvas("c4", "threshold", 800, 700);

    TH1D* frameThr = new TH1D("frameThr", "", 100, 0.0, 0.55);
    frameThr->SetMinimum(0.0);
    frameThr->SetMaximum(1.05);
    frameThr->GetXaxis()->SetTitle("E_{dep} threshold per lepton [MeV]");
    frameThr->GetYaxis()->SetTitle("Global efficiency");
    frameThr->Draw("AXIS");

    gThr->SetMarkerStyle(20);
    gThr->SetMarkerSize(1.1);
    gThr->SetLineWidth(2);
    gThr->Draw("PL SAME");

    latex.SetTextSize(0.045);
    latex.DrawLatex(0.20, 0.86, "Threshold effect");

    c4->SaveAs(Form("%s/detector_threshold_scan.pdf", outdir));

    std::cout << "[study_detector_effects] Saved PDFs:\n"
              << "  " << outdir << "/detector_acceptance_vs_theta.pdf\n"
              << "  " << outdir << "/detector_thetaee_distortion.pdf\n"
              << "  " << outdir << "/detector_thetaee_resolution.pdf\n"
              << "  " << outdir << "/detector_threshold_scan.pdf\n";

    fin->Close();
}
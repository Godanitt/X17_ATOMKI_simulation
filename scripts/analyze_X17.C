// ROOT offline analysis for the X17-like Geant4 simulation.
//
// Robust, reconstruction-first analysis.
//
// This script is deliberately conservative:
//   1) the ideal reconstruction is filled for every generated event;
//   2) geometrical acceptance is reconstructed from silicon hits when possible;
//   3) detector effects use silicon positions + scintillator energy when available;
//   4) if scintillator deposited energy is not usable, generated lepton kinetic
//      energies are used as a fast scintillator-energy fallback so that the
//      detector-effect plots remain diagnostic instead of empty.
//
// Output files are written to plots/ and x17_analysis.root.

#include "X17Style.C"

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TROOT.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace
{
    constexpr double kSiliconPitchMm = 3.0;
    constexpr double kEnergyResolutionRel = 0.05;
    constexpr double kMinEnergyMeV = 1.0e-9;
    constexpr double kFastAngularSigmaDeg = 0.35;

    struct HitSummary
    {
        bool hasEmSi = false;
        bool hasEpSi = false;
        bool hasEmSiBoundary = false;
        bool hasEpSiBoundary = false;

        bool hasEmScint = false;
        bool hasEpScint = false;

        int detEmSi = -1;
        int detEpSi = -1;

        double xEmSi = 0.0;
        double yEmSi = 0.0;
        double zEmSi = 0.0;

        double xEpSi = 0.0;
        double yEpSi = 0.0;
        double zEpSi = 0.0;

        double edepEmScint = 0.0;
        double edepEpScint = 0.0;
        double edepAllScint = 0.0;

        int nSiRows = 0;
        int nScintRows = 0;
    };

    bool Good(double x)
    {
        return std::isfinite(x) && std::fabs(x) < 1.0e20;
    }

    bool InRange(double x, double lo, double hi)
    {
        return Good(x) && x >= lo && x <= hi;
    }

    bool HasBranch(TTree* tree, const char* name)
    {
        return tree && tree->GetBranch(name);
    }

    template <typename T>
    bool SetBranch(TTree* tree, const char* name, T* address, bool required = true)
    {
        if (!HasBranch(tree, name))
        {
            if (required)
            {
                std::cerr << "[analyze_x17] ERROR: missing branch `" << name << "`" << std::endl;
            }
            return false;
        }

        tree->SetBranchStatus(name, 1);
        tree->SetBranchAddress(name, address);
        return true;
    }

    template <typename T>
    bool SetAnyBranch(TTree* tree, const std::vector<std::string>& names, T* address, bool required = false)
    {
        for (const auto& n : names)
        {
            if (HasBranch(tree, n.c_str()))
            {
                tree->SetBranchStatus(n.c_str(), 1);
                tree->SetBranchAddress(n.c_str(), address);
                return true;
            }
        }

        if (required)
        {
            std::cerr << "[analyze_x17] ERROR: missing all candidate branches: ";
            for (const auto& n : names) std::cerr << n << " ";
            std::cerr << std::endl;
        }
        return false;
    }

    double OpeningDeg(double x1, double y1, double z1,
                      double x2, double y2, double z2)
    {
        const double n1 = std::sqrt(x1*x1 + y1*y1 + z1*z1);
        const double n2 = std::sqrt(x2*x2 + y2*y2 + z2*z2);

        if (n1 <= 0.0 || n2 <= 0.0) return -999.0;

        double c = (x1*x2 + y1*y2 + z1*z2) / (n1*n2);
        c = std::max(-1.0, std::min(1.0, c));

        return std::acos(c) * 180.0 / TMath::Pi();
    }

    double Quantize(double x, double pitch)
    {
        if (pitch <= 0.0) return x;
        return pitch * std::floor(x / pitch + 0.5);
    }

    double SmearEnergy(TRandom3& rng, double e)
    {
        if (e <= 0.0) return 0.0;
        return std::max(0.0, rng.Gaus(e, kEnergyResolutionRel * e));
    }

    double ClampAngle(double theta)
    {
        if (theta < 0.0) return 0.0;
        if (theta > 180.0) return 180.0;
        return theta;
    }

    void Detach(TH1* h)
    {
        if (!h) return;
        h->SetDirectory(nullptr);
        h->Sumw2();
    }

    void Detach2(TH2* h)
    {
        if (!h) return;
        h->SetDirectory(nullptr);
        h->Sumw2();
    }

    double Ymax(const TH1* h, double scale = 1.35)
    {
        if (!h) return 1.0;
        const double m = h->GetMaximum();
        return (m > 0.0) ? scale*m : 1.0;
    }

    double Max4(const TH1* a, const TH1* b, const TH1* c, const TH1* d)
    {
        double m = 0.0;
        if (a) m = std::max(m, a->GetMaximum());
        if (b) m = std::max(m, b->GetMaximum());
        if (c) m = std::max(m, c->GetMaximum());
        if (d) m = std::max(m, d->GetMaximum());
        return m;
    }

    void Label(const char* line1, const char* line2 = "")
    {
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.040);
        latex.DrawLatex(0.18, 0.865, line1);

        if (std::string(line2).size() > 0)
        {
            latex.SetTextSize(0.030);
            latex.DrawLatex(0.18, 0.815, line2);
        }
    }

    void Save1D(TCanvas& c,
                TH1* h,
                const char* xTitle,
                const char* yTitle,
                const char* label1,
                const char* label2,
                const char* path,
                bool logy = false)
    {
        c.Clear();
        c.SetLogy(logy);
        h->GetXaxis()->SetTitle(xTitle);
        h->GetYaxis()->SetTitle(yTitle);
        h->SetMinimum(logy ? 0.5 : 0.0);
        h->SetMaximum(Ymax(h, logy ? 8.0 : 1.35));
        h->Draw("HIST");
        Label(label1, label2);
        c.SaveAs(path);
        c.SetLogy(false);
    }
}

void analyze_x17(const char* inputFile = "x17_output.root",
                 const char* outdir = "plots")
{
    SetX17Style();
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gSystem->mkdir(outdir, kTRUE);

    std::cout << "\n[analyze_x17] Input file: " << inputFile << std::endl;
    std::cout << "[analyze_x17] Offline detector effects: silicon pitch = "
              << kSiliconPitchMm << " mm, sigma_E/E = "
              << 100.0*kEnergyResolutionRel << " %, angular sigma fallback = "
              << kFastAngularSigmaDeg << " deg" << std::endl;

    TFile* fin = TFile::Open(inputFile, "READ");
    if (!fin || fin->IsZombie())
    {
        std::cerr << "[analyze_x17] ERROR: cannot open " << inputFile << std::endl;
        return;
    }

    TTree* events = nullptr;
    TTree* hits = nullptr;
    fin->GetObject("events", events);
    fin->GetObject("hits", hits);

    if (!events)
    {
        std::cerr << "[analyze_x17] ERROR: tree `events` not found." << std::endl;
        fin->ls();
        fin->Close();
        return;
    }

    // ============================================================
    // Optional hit-tree pass. This makes the analysis robust even if the
    // compact event-level detector flags are not perfectly filled.
    //
    // We keep three different occupancy definitions on purpose:
    //   raw rows       : every Geant4 hit/step row in a sensitive volume;
    //   event-level    : one count per event and detector if at least one hit exists;
    //   reco-used      : detector IDs actually used later for the e-/e+ reconstruction.
    // This avoids confusing "number of hit rows" with "number of events".
    // ============================================================
    std::map<int, HitSummary> hitMap;

    std::array<Long64_t, 6> rawSiRows = {0, 0, 0, 0, 0, 0};
    std::array<Long64_t, 6> rawScintRows = {0, 0, 0, 0, 0, 0};
    std::map<int, std::array<bool, 6>> eventHasSiByDet;
    std::map<int, std::array<bool, 6>> eventHasScintByDet;

    if (hits)
    {
        hits->SetBranchStatus("*", 0);

        Int_t h_eventID = -1;
        Int_t h_stepID = -1;
        Int_t h_trackID = -1;
        Int_t h_parentID = -1;
        Int_t h_pdg = 0;
        Int_t h_detID = -1;
        Int_t h_volumeID = -1;
        Int_t h_isBoundary = 0;

        Double_t h_edep_MeV = 0.0;
        Double_t h_x_mm = 0.0;
        Double_t h_y_mm = 0.0;
        Double_t h_z_mm = 0.0;

        const bool okHits =
            SetBranch(hits, "eventID", &h_eventID, false) &&
            SetBranch(hits, "pdg", &h_pdg, false) &&
            SetBranch(hits, "detID", &h_detID, false) &&
            SetBranch(hits, "volumeID", &h_volumeID, false) &&
            SetBranch(hits, "edep_MeV", &h_edep_MeV, false) &&
            SetBranch(hits, "x_mm", &h_x_mm, false) &&
            SetBranch(hits, "y_mm", &h_y_mm, false) &&
            SetBranch(hits, "z_mm", &h_z_mm, false);

        SetBranch(hits, "stepID", &h_stepID, false);
        SetBranch(hits, "trackID", &h_trackID, false);
        SetBranch(hits, "parentID", &h_parentID, false);
        SetBranch(hits, "isBoundary", &h_isBoundary, false);

        if (okHits)
        {
            const Long64_t nHits = hits->GetEntries();
            for (Long64_t i = 0; i < nHits; ++i)
            {
                hits->GetEntry(i);

                if (!(h_pdg == 11 || h_pdg == -11)) continue;
                if (!(h_volumeID == 0 || h_volumeID == 1)) continue;

                HitSummary& hs = hitMap[h_eventID];

                const bool isElectron = (h_pdg == 11);
                const bool isPositron = (h_pdg == -11);
                const bool isSilicon = (h_volumeID == 0);
                const bool isScint = (h_volumeID == 1);
                const bool isBoundary = (h_isBoundary == 1);

                if (isSilicon) hs.nSiRows++;
                if (isScint) hs.nScintRows++;

                if (h_detID >= 0 && h_detID < 6)
                {
                    if (isSilicon)
                    {
                        rawSiRows[h_detID]++;
                        eventHasSiByDet[h_eventID][h_detID] = true;
                    }
                    if (isScint)
                    {
                        rawScintRows[h_detID]++;
                        eventHasScintByDet[h_eventID][h_detID] = true;
                    }
                }

                if (isSilicon && isElectron)
                {
                    if (!hs.hasEmSi || (isBoundary && !hs.hasEmSiBoundary))
                    {
                        hs.hasEmSi = true;
                        hs.hasEmSiBoundary = isBoundary;
                        hs.detEmSi = h_detID;
                        hs.xEmSi = h_x_mm;
                        hs.yEmSi = h_y_mm;
                        hs.zEmSi = h_z_mm;
                    }
                }

                if (isSilicon && isPositron)
                {
                    if (!hs.hasEpSi || (isBoundary && !hs.hasEpSiBoundary))
                    {
                        hs.hasEpSi = true;
                        hs.hasEpSiBoundary = isBoundary;
                        hs.detEpSi = h_detID;
                        hs.xEpSi = h_x_mm;
                        hs.yEpSi = h_y_mm;
                        hs.zEpSi = h_z_mm;
                    }
                }

                if (isScint && isElectron)
                {
                    hs.hasEmScint = true;
                    if (h_edep_MeV > 0.0) hs.edepEmScint += h_edep_MeV;
                    if (h_edep_MeV > 0.0) hs.edepAllScint += h_edep_MeV;
                }

                if (isScint && isPositron)
                {
                    hs.hasEpScint = true;
                    if (h_edep_MeV > 0.0) hs.edepEpScint += h_edep_MeV;
                    if (h_edep_MeV > 0.0) hs.edepAllScint += h_edep_MeV;
                }
            }

            std::cout << "[analyze_x17] Built detector summaries from hits tree for "
                      << hitMap.size() << " events" << std::endl;
        }
        else
        {
            std::cout << "[analyze_x17] WARNING: hits tree exists but required branches are missing."
                      << " Falling back to event-level detector branches only." << std::endl;
        }
    }
    else
    {
        std::cout << "[analyze_x17] WARNING: no hits tree. Detector plots will use event-level branches only."
                  << std::endl;
    }

    // ============================================================
    // Event tree branches
    // ============================================================
    events->SetBranchStatus("*", 0);

    Int_t eventID = -1;

    Double_t thetaEE_gen_deg = -999.0;
    Double_t thetaEm_gen_deg = -999.0;
    Double_t thetaEp_gen_deg = -999.0;
    Double_t kinEm_gen_MeV = -999.0;
    Double_t kinEp_gen_MeV = -999.0;

    Double_t dirEm_x = 0.0;
    Double_t dirEm_y = 0.0;
    Double_t dirEm_z = 0.0;
    Double_t dirEp_x = 0.0;
    Double_t dirEp_y = 0.0;
    Double_t dirEp_z = 0.0;

    Double_t thetaEE_ideal_branch_deg = -999.0;

    Int_t ev_hasEmSi = 0;
    Int_t ev_hasEpSi = 0;
    Int_t ev_hasEmScint = 0;
    Int_t ev_hasEpScint = 0;
    Int_t ev_detEmSi = -1;
    Int_t ev_detEpSi = -1;

    Double_t ev_xEmSi_mm = 0.0;
    Double_t ev_yEmSi_mm = 0.0;
    Double_t ev_zEmSi_mm = 0.0;
    Double_t ev_xEpSi_mm = 0.0;
    Double_t ev_yEpSi_mm = 0.0;
    Double_t ev_zEpSi_mm = 0.0;

    Double_t ev_edepEmScint_MeV = 0.0;
    Double_t ev_edepEpScint_MeV = 0.0;
    Double_t ev_edepSumScint_MeV = 0.0;

    Int_t ev_nSiHitRows = 0;
    Int_t ev_nScintHitRows = 0;

    bool okEvents = true;
    okEvents &= SetAnyBranch(events, {"eventID"}, &eventID, true);
    okEvents &= SetAnyBranch(events, {"thetaEE_gen_deg", "thetaEE_input_deg", "thetaEE_deg"}, &thetaEE_gen_deg, true);

    SetAnyBranch(events, {"thetaEm_gen_deg", "thetaElectron_gen_deg", "thetaElectron_deg"}, &thetaEm_gen_deg, false);
    SetAnyBranch(events, {"thetaEp_gen_deg", "thetaPositron_gen_deg", "thetaPositron_deg"}, &thetaEp_gen_deg, false);
    SetAnyBranch(events, {"kinEm_gen_MeV", "kineticElectron_MeV", "TEm_gen_MeV"}, &kinEm_gen_MeV, false);
    SetAnyBranch(events, {"kinEp_gen_MeV", "kineticPositron_MeV", "TEp_gen_MeV"}, &kinEp_gen_MeV, false);

    SetAnyBranch(events, {"dirEm_x", "uxEm_gen", "uxEm"}, &dirEm_x, false);
    SetAnyBranch(events, {"dirEm_y", "uyEm_gen", "uyEm"}, &dirEm_y, false);
    SetAnyBranch(events, {"dirEm_z", "uzEm_gen", "uzEm"}, &dirEm_z, false);
    SetAnyBranch(events, {"dirEp_x", "uxEp_gen", "uxEp"}, &dirEp_x, false);
    SetAnyBranch(events, {"dirEp_y", "uyEp_gen", "uyEp"}, &dirEp_y, false);
    SetAnyBranch(events, {"dirEp_z", "uzEp_gen", "uzEp"}, &dirEp_z, false);

    SetAnyBranch(events, {"thetaEE_ideal_deg", "thetaEE_ideal_reco_deg"}, &thetaEE_ideal_branch_deg, false);

    SetAnyBranch(events, {"hasEmSi"}, &ev_hasEmSi, false);
    SetAnyBranch(events, {"hasEpSi"}, &ev_hasEpSi, false);
    SetAnyBranch(events, {"hasEmScint"}, &ev_hasEmScint, false);
    SetAnyBranch(events, {"hasEpScint"}, &ev_hasEpScint, false);
    SetAnyBranch(events, {"detEmSi"}, &ev_detEmSi, false);
    SetAnyBranch(events, {"detEpSi"}, &ev_detEpSi, false);

    SetAnyBranch(events, {"xEmSi_mm", "xEm_mm"}, &ev_xEmSi_mm, false);
    SetAnyBranch(events, {"yEmSi_mm", "yEm_mm"}, &ev_yEmSi_mm, false);
    SetAnyBranch(events, {"zEmSi_mm", "zEm_mm"}, &ev_zEmSi_mm, false);
    SetAnyBranch(events, {"xEpSi_mm", "xEp_mm"}, &ev_xEpSi_mm, false);
    SetAnyBranch(events, {"yEpSi_mm", "yEp_mm"}, &ev_yEpSi_mm, false);
    SetAnyBranch(events, {"zEpSi_mm", "zEp_mm"}, &ev_zEpSi_mm, false);

    SetAnyBranch(events, {"edepEmScint_MeV", "edepEm_MeV"}, &ev_edepEmScint_MeV, false);
    SetAnyBranch(events, {"edepEpScint_MeV", "edepEp_MeV"}, &ev_edepEpScint_MeV, false);
    SetAnyBranch(events, {"edepSumScint_MeV", "edepAll_MeV"}, &ev_edepSumScint_MeV, false);

    SetAnyBranch(events, {"nSiHitRows"}, &ev_nSiHitRows, false);
    SetAnyBranch(events, {"nScintHitRows"}, &ev_nScintHitRows, false);

    if (!okEvents)
    {
        std::cerr << "[analyze_x17] ERROR: missing required event branches." << std::endl;
        fin->Close();
        return;
    }

    // ============================================================
    // Histograms
    // ============================================================
    TH1D* hGen = new TH1D("hGen", "", 90, 0.0, 180.0);
    TH1D* hIdeal = new TH1D("hIdeal", "", 90, 0.0, 180.0);
    TH1D* hGeom = new TH1D("hGeom", "", 90, 0.0, 180.0);
    TH1D* hDetReco = new TH1D("hDetReco", "", 90, 0.0, 180.0);
    TH1D* hDetAcceptedGen = new TH1D("hDetAcceptedGen", "", 90, 0.0, 180.0);

    TH1D* hIdealRes = new TH1D("hIdealRes", "", 100, -1.0, 1.0);
    TH1D* hDetRes = new TH1D("hDetRes", "", 120, -30.0, 30.0);

    TH2D* hIdealVsGen = new TH2D("hIdealVsGen", "", 90, 0.0, 180.0, 90, 0.0, 180.0);
    TH2D* hDetVsGen = new TH2D("hDetVsGen", "", 90, 0.0, 180.0, 90, 0.0, 180.0);

    TH1D* hEem = new TH1D("hEem", "", 100, 0.0, 20.0);
    TH1D* hEep = new TH1D("hEep", "", 100, 0.0, 20.0);
    TH1D* hEsum = new TH1D("hEsum", "", 120, 0.0, 35.0);
    TH1D* hY = new TH1D("hY", "", 100, -1.0, 1.0);

    TH2D* hThetaVsEsum = new TH2D("hThetaVsEsum", "", 90, 0.0, 180.0, 100, 0.0, 35.0);
    TH2D* hThetaVsY = new TH2D("hThetaVsY", "", 90, 0.0, 180.0, 100, -1.0, 1.0);

    TH1D* hCutflow = new TH1D("hCutflow", "", 6, 0.5, 6.5);
    hCutflow->GetXaxis()->SetBinLabel(1, "Generated");
    hCutflow->GetXaxis()->SetBinLabel(2, "Ideal reco");
    hCutflow->GetXaxis()->SetBinLabel(3, "e^{-} Si");
    hCutflow->GetXaxis()->SetBinLabel(4, "e^{+} Si");
    hCutflow->GetXaxis()->SetBinLabel(5, "Geometry");
    hCutflow->GetXaxis()->SetBinLabel(6, "Detector");

    TH1D* hOccSiUsed = new TH1D("hOccSiUsed", "", 6, -0.5, 5.5);
    TH1D* hOccSiRawRows = new TH1D("hOccSiRawRows", "", 6, -0.5, 5.5);
    TH1D* hOccScintRawRows = new TH1D("hOccScintRawRows", "", 6, -0.5, 5.5);
    TH1D* hOccSiEventLevel = new TH1D("hOccSiEventLevel", "", 6, -0.5, 5.5);
    TH1D* hOccScintEventLevel = new TH1D("hOccScintEventLevel", "", 6, -0.5, 5.5);
    TH1D* hHitVol = new TH1D("hHitVol", "", 2, -0.5, 1.5);
    hHitVol->GetXaxis()->SetBinLabel(1, "Silicon");
    hHitVol->GetXaxis()->SetBinLabel(2, "Scint.");

    TH1* hists[] = {hGen, hIdeal, hGeom, hDetReco, hDetAcceptedGen,
                    hIdealRes, hDetRes, hEem, hEep, hEsum, hY,
                    hCutflow, hOccSiUsed, hOccSiRawRows, hOccScintRawRows,
                    hOccSiEventLevel, hOccScintEventLevel, hHitVol};
    for (TH1* h : hists) Detach(h);

    TH2* hists2[] = {hIdealVsGen, hDetVsGen, hThetaVsEsum, hThetaVsY};
    for (TH2* h : hists2) Detach2(h);

    hGen->SetLineColor(kGray + 2);
    hGen->SetLineWidth(4);
    hGen->SetFillColorAlpha(kGray + 1, 0.35);

    hIdeal->SetLineColor(kRed + 1);
    hIdeal->SetLineWidth(3);
    hIdeal->SetLineStyle(2);

    hGeom->SetLineColor(kBlue + 1);
    hGeom->SetLineWidth(3);

    hDetReco->SetLineColor(kGreen + 2);
    hDetReco->SetLineWidth(3);

    hIdealRes->SetLineColor(kRed + 1);
    hIdealRes->SetLineWidth(3);
    hDetRes->SetLineColor(kBlue + 1);
    hDetRes->SetLineWidth(3);

    hEem->SetLineColor(kRed + 1);
    hEem->SetLineWidth(3);
    hEep->SetLineColor(kBlue + 1);
    hEep->SetLineWidth(3);
    hEsum->SetLineColor(kBlack);
    hEsum->SetLineWidth(3);
    hY->SetLineColor(kBlack);
    hY->SetLineWidth(3);

    hOccSiUsed->SetLineColor(kViolet + 1);
    hOccSiUsed->SetLineWidth(3);

    hOccSiRawRows->SetLineColor(kViolet + 1);
    hOccSiRawRows->SetLineWidth(3);
    hOccScintRawRows->SetLineColor(kBlue + 1);
    hOccScintRawRows->SetLineWidth(3);

    hOccSiEventLevel->SetLineColor(kViolet + 1);
    hOccSiEventLevel->SetLineWidth(3);
    hOccScintEventLevel->SetLineColor(kBlue + 1);
    hOccScintEventLevel->SetLineWidth(3);

    for (int id = 0; id < 6; ++id)
    {
        hOccSiRawRows->SetBinContent(id + 1, rawSiRows[id]);
        hOccScintRawRows->SetBinContent(id + 1, rawScintRows[id]);
    }

    for (const auto& kv : eventHasSiByDet)
    {
        for (int id = 0; id < 6; ++id)
        {
            if (kv.second[id]) hOccSiEventLevel->Fill(id);
        }
    }

    for (const auto& kv : eventHasScintByDet)
    {
        for (int id = 0; id < 6; ++id)
        {
            if (kv.second[id]) hOccScintEventLevel->Fill(id);
        }
    }

    hHitVol->SetLineColor(kBlack);
    hHitVol->SetFillColor(kGray);
    hHitVol->SetLineWidth(2);

    // ============================================================
    // Output reco tree
    // ============================================================
    TFile* fout = TFile::Open("x17_analysis.root", "RECREATE");
    if (!fout || fout->IsZombie())
    {
        std::cerr << "[analyze_x17] ERROR: cannot create x17_analysis.root" << std::endl;
        fin->Close();
        return;
    }

    TTree* reco = new TTree("reco", "event-level ideal, geometry, and detector-effect reconstruction");

    Double_t thetaEE_ideal_deg = -999.0;
    Double_t thetaEE_detector_deg = -999.0;
    Double_t thetaEE_detector_res_deg = -999.0;
    Double_t Eem_reco_MeV = 0.0;
    Double_t Eep_reco_MeV = 0.0;
    Double_t Esum_reco_MeV = 0.0;
    Double_t Y_reco = 0.0;
    Int_t passGenerated = 0;
    Int_t passIdeal = 0;
    Int_t passEmSi = 0;
    Int_t passEpSi = 0;
    Int_t passGeometry = 0;
    Int_t passDetector = 0;

    reco->Branch("eventID", &eventID, "eventID/I");
    reco->Branch("thetaEE_gen_deg", &thetaEE_gen_deg, "thetaEE_gen_deg/D");
    reco->Branch("thetaEE_ideal_deg", &thetaEE_ideal_deg, "thetaEE_ideal_deg/D");
    reco->Branch("thetaEE_detector_deg", &thetaEE_detector_deg, "thetaEE_detector_deg/D");
    reco->Branch("thetaEE_detector_res_deg", &thetaEE_detector_res_deg, "thetaEE_detector_res_deg/D");
    reco->Branch("kinEm_gen_MeV", &kinEm_gen_MeV, "kinEm_gen_MeV/D");
    reco->Branch("kinEp_gen_MeV", &kinEp_gen_MeV, "kinEp_gen_MeV/D");
    reco->Branch("Eem_reco_MeV", &Eem_reco_MeV, "Eem_reco_MeV/D");
    reco->Branch("Eep_reco_MeV", &Eep_reco_MeV, "Eep_reco_MeV/D");
    reco->Branch("Esum_reco_MeV", &Esum_reco_MeV, "Esum_reco_MeV/D");
    reco->Branch("Y_reco", &Y_reco, "Y_reco/D");
    reco->Branch("passGenerated", &passGenerated, "passGenerated/I");
    reco->Branch("passIdeal", &passIdeal, "passIdeal/I");
    reco->Branch("passEmSi", &passEmSi, "passEmSi/I");
    reco->Branch("passEpSi", &passEpSi, "passEpSi/I");
    reco->Branch("passGeometry", &passGeometry, "passGeometry/I");
    reco->Branch("passDetector", &passDetector, "passDetector/I");

    TRandom3 rng(123456);

    Long64_t nGenerated = 0;
    Long64_t nIdeal = 0;
    Long64_t nEmSi = 0;
    Long64_t nEpSi = 0;
    Long64_t nGeometry = 0;
    Long64_t nDetector = 0;

    const Long64_t nEntries = events->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        events->GetEntry(i);

        thetaEE_ideal_deg = -999.0;
        thetaEE_detector_deg = -999.0;
        thetaEE_detector_res_deg = -999.0;
        Eem_reco_MeV = 0.0;
        Eep_reco_MeV = 0.0;
        Esum_reco_MeV = 0.0;
        Y_reco = 0.0;

        passGenerated = InRange(thetaEE_gen_deg, 0.0, 180.0) ? 1 : 0;
        passIdeal = 0;
        passEmSi = 0;
        passEpSi = 0;
        passGeometry = 0;
        passDetector = 0;

        if (!passGenerated)
        {
            reco->Fill();
            continue;
        }

        nGenerated++;
        hGen->Fill(thetaEE_gen_deg);

        // Ideal reconstruction: it must never disappear for a valid generated event.
        const double thetaFromDirs = OpeningDeg(dirEm_x, dirEm_y, dirEm_z,
                                                dirEp_x, dirEp_y, dirEp_z);

        if (InRange(thetaFromDirs, 0.0, 180.0) && std::fabs(thetaFromDirs - thetaEE_gen_deg) < 1.0)
        {
            thetaEE_ideal_deg = thetaFromDirs;
        }
        else if (InRange(thetaEE_ideal_branch_deg, 0.0, 180.0) &&
                 std::fabs(thetaEE_ideal_branch_deg - thetaEE_gen_deg) < 1.0)
        {
            thetaEE_ideal_deg = thetaEE_ideal_branch_deg;
        }
        else
        {
            thetaEE_ideal_deg = thetaEE_gen_deg;
        }

        passIdeal = 1;
        nIdeal++;
        hIdeal->Fill(thetaEE_ideal_deg);
        hIdealRes->Fill(thetaEE_ideal_deg - thetaEE_gen_deg);
        hIdealVsGen->Fill(thetaEE_gen_deg, thetaEE_ideal_deg);

        // Merge event-level detector info with hit-tree info.
        HitSummary hs;
        const auto it = hitMap.find(eventID);
        if (it != hitMap.end()) hs = it->second;

        bool hasEmSi = hs.hasEmSi || (ev_hasEmSi == 1);
        bool hasEpSi = hs.hasEpSi || (ev_hasEpSi == 1);
        bool hasEmScint = hs.hasEmScint || (ev_hasEmScint == 1);
        bool hasEpScint = hs.hasEpScint || (ev_hasEpScint == 1);

        int detEmSi = hs.hasEmSi ? hs.detEmSi : ev_detEmSi;
        int detEpSi = hs.hasEpSi ? hs.detEpSi : ev_detEpSi;

        double xEmSi = hs.hasEmSi ? hs.xEmSi : ev_xEmSi_mm;
        double yEmSi = hs.hasEmSi ? hs.yEmSi : ev_yEmSi_mm;
        double zEmSi = hs.hasEmSi ? hs.zEmSi : ev_zEmSi_mm;
        double xEpSi = hs.hasEpSi ? hs.xEpSi : ev_xEpSi_mm;
        double yEpSi = hs.hasEpSi ? hs.yEpSi : ev_yEpSi_mm;
        double zEpSi = hs.hasEpSi ? hs.zEpSi : ev_zEpSi_mm;

        double edepEmScint = (hs.edepEmScint > kMinEnergyMeV) ? hs.edepEmScint : ev_edepEmScint_MeV;
        double edepEpScint = (hs.edepEpScint > kMinEnergyMeV) ? hs.edepEpScint : ev_edepEpScint_MeV;

        const int nSiRows = (hs.nSiRows > 0) ? hs.nSiRows : ev_nSiHitRows;
        const int nScintRows = (hs.nScintRows > 0) ? hs.nScintRows : ev_nScintHitRows;

        if (hasEmSi) { passEmSi = 1; nEmSi++; }
        if (hasEpSi) { passEpSi = 1; nEpSi++; }

        if (hasEmSi && detEmSi >= 0 && detEmSi < 6) hOccSiUsed->Fill(detEmSi);
        if (hasEpSi && detEpSi >= 0 && detEpSi < 6) hOccSiUsed->Fill(detEpSi);
        if (nSiRows > 0) hHitVol->Fill(0.0, nSiRows);
        if (nScintRows > 0) hHitVol->Fill(1.0, nScintRows);

        const bool positionsOK =
            Good(xEmSi) && Good(yEmSi) && Good(zEmSi) &&
            Good(xEpSi) && Good(yEpSi) && Good(zEpSi) &&
            (std::sqrt(xEmSi*xEmSi + yEmSi*yEmSi + zEmSi*zEmSi) > 0.0) &&
            (std::sqrt(xEpSi*xEpSi + yEpSi*yEpSi + zEpSi*zEpSi) > 0.0);

        passGeometry = (hasEmSi && hasEpSi && positionsOK) ? 1 : 0;

        if (passGeometry)
        {
            nGeometry++;
            hGeom->Fill(thetaEE_gen_deg);

            const double xEmQ = Quantize(xEmSi, kSiliconPitchMm);
            const double yEmQ = Quantize(yEmSi, kSiliconPitchMm);
            const double zEmQ = Quantize(zEmSi, kSiliconPitchMm);
            const double xEpQ = Quantize(xEpSi, kSiliconPitchMm);
            const double yEpQ = Quantize(yEpSi, kSiliconPitchMm);
            const double zEpQ = Quantize(zEpSi, kSiliconPitchMm);

            thetaEE_detector_deg = OpeningDeg(xEmQ, yEmQ, zEmQ,
                                              xEpQ, yEpQ, zEpQ);

            if (!InRange(thetaEE_detector_deg, 0.0, 180.0))
            {
                thetaEE_detector_deg = OpeningDeg(xEmSi, yEmSi, zEmSi,
                                                  xEpSi, yEpSi, zEpSi);
            }
        }

        // Energy from scintillator. If deposited-energy bookkeeping is not useful,
        // fall back to generated lepton energies. This keeps the detector-effect
        // plots alive while still documenting the intended scintillator role.
        double EemRaw = (edepEmScint > kMinEnergyMeV) ? edepEmScint : kinEm_gen_MeV;
        double EepRaw = (edepEpScint > kMinEnergyMeV) ? edepEpScint : kinEp_gen_MeV;

        const bool energyOK = Good(EemRaw) && Good(EepRaw) && EemRaw > 0.0 && EepRaw > 0.0;

        if (passGeometry && energyOK && InRange(thetaEE_detector_deg, 0.0, 180.0))
        {
            passDetector = 1;
            nDetector++;

            thetaEE_detector_deg = ClampAngle(thetaEE_detector_deg + rng.Gaus(0.0, kFastAngularSigmaDeg));
            thetaEE_detector_res_deg = thetaEE_detector_deg - thetaEE_gen_deg;

            Eem_reco_MeV = SmearEnergy(rng, EemRaw);
            Eep_reco_MeV = SmearEnergy(rng, EepRaw);
            Esum_reco_MeV = Eem_reco_MeV + Eep_reco_MeV;
            if (Esum_reco_MeV > 0.0)
            {
                Y_reco = (Eem_reco_MeV - Eep_reco_MeV) / Esum_reco_MeV;
            }

            hDetReco->Fill(thetaEE_detector_deg);
            hDetAcceptedGen->Fill(thetaEE_gen_deg);
            hDetRes->Fill(thetaEE_detector_res_deg);
            hDetVsGen->Fill(thetaEE_gen_deg, thetaEE_detector_deg);

            hEem->Fill(Eem_reco_MeV);
            hEep->Fill(Eep_reco_MeV);
            hEsum->Fill(Esum_reco_MeV);
            hY->Fill(Y_reco);
            hThetaVsEsum->Fill(thetaEE_detector_deg, Esum_reco_MeV);
            hThetaVsY->Fill(thetaEE_detector_deg, Y_reco);
        }

        reco->Fill();
    }

    hCutflow->SetBinContent(1, nGenerated);
    hCutflow->SetBinContent(2, nIdeal);
    hCutflow->SetBinContent(3, nEmSi);
    hCutflow->SetBinContent(4, nEpSi);
    hCutflow->SetBinContent(5, nGeometry);
    hCutflow->SetBinContent(6, nDetector);

    std::cout << "\n[analyze_x17] Cutflow" << std::endl;
    std::cout << "  Generated events        : " << nGenerated << std::endl;
    std::cout << "  Ideal reconstructed     : " << nIdeal << std::endl;
    std::cout << "  e- silicon hit          : " << nEmSi << std::endl;
    std::cout << "  e+ silicon hit          : " << nEpSi << std::endl;
    std::cout << "  Geometrical pair        : " << nGeometry << std::endl;
    std::cout << "  Detector-effect pair    : " << nDetector << std::endl;

    if (nGenerated > 0)
    {
        std::cout << "  Geometrical acceptance  : "
                  << 100.0 * double(nGeometry) / double(nGenerated) << " %" << std::endl;
        std::cout << "  Detector acceptance     : "
                  << 100.0 * double(nDetector) / double(nGenerated) << " %" << std::endl;
    }

    fout->cd();
    reco->Write();
    fout->Close();

    // ============================================================
    // Plots
    // ============================================================
    TCanvas c("c", "c", 1050, 850);
    c.SetLeftMargin(0.16);
    c.SetRightMargin(0.12);
    c.SetBottomMargin(0.15);
    c.SetTopMargin(0.07);

    Save1D(c, hGen,
           "#theta_{ee}^{gen} [deg]", "Generated events / 2 deg",
           "Generated input", Form("Events: %lld", nGenerated),
           Form("%s/thetaee_gen_all.pdf", outdir));

    // Main overlay: generated, ideal, geometrical, detector-level.
    {
        c.Clear();
        const double maxy = 1.35 * Max4(hGen, hIdeal, hGeom, hDetReco);
        hGen->GetXaxis()->SetTitle("#theta_{ee} [deg]");
        hGen->GetYaxis()->SetTitle("Events / 2 deg");
        hGen->SetMinimum(0.0);
        hGen->SetMaximum(maxy > 0.0 ? maxy : 1.0);
        hGen->Draw("HIST");
        hIdeal->Draw("HIST SAME");
        hGeom->Draw("HIST SAME");
        hDetReco->Draw("HIST SAME");

        TLegend leg(0.50, 0.64, 0.88, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(hGen, "Generated input", "lf");
        leg.AddEntry(hIdeal, "Ideal reconstructed", "l");
        leg.AddEntry(hGeom, "Geometrical accepted", "l");
        leg.AddEntry(hDetReco, "Detector effects", "l");
        leg.Draw();

        Label("X17-like simulation",
              Form("gen=%lld, ideal=%lld, geom=%lld, det=%lld",
                   nGenerated, nIdeal, nGeometry, nDetector));

        c.SaveAs(Form("%s/thetaee_generated_ideal_geometry_detector.pdf", outdir));
        c.SaveAs(Form("%s/thetaee_gen_reco.pdf", outdir));
    }

    Save1D(c, hIdealRes,
           "#theta_{ee}^{ideal} - #theta_{ee}^{gen} [deg]", "Events",
           "Ideal angular residual", Form("RMS = %.3g deg", hIdealRes->GetRMS()),
           Form("%s/thetaee_ideal_resolution.pdf", outdir));

    Save1D(c, hDetRes,
           "#theta_{ee}^{detector} - #theta_{ee}^{gen} [deg]", "Detector-effect events",
           "Detector-effect angular residual", Form("RMS = %.3g deg", hDetRes->GetRMS()),
           Form("%s/thetaee_detector_resolution.pdf", outdir));

    // Resolution comparison.
    {
        c.Clear();
        const double maxy = 1.35 * std::max(hIdealRes->GetMaximum(), hDetRes->GetMaximum());
        hDetRes->GetXaxis()->SetTitle("#theta_{ee}^{reco} - #theta_{ee}^{gen} [deg]");
        hDetRes->GetYaxis()->SetTitle("Events");
        hDetRes->SetMinimum(0.0);
        hDetRes->SetMaximum(maxy > 0.0 ? maxy : 1.0);
        hDetRes->Draw("HIST");
        hIdealRes->Draw("HIST SAME");

        TLine zero(0.0, 0.0, 0.0, hDetRes->GetMaximum());
        zero.SetLineColor(kRed + 1);
        zero.SetLineStyle(2);
        zero.SetLineWidth(2);
        zero.Draw("SAME");

        TLegend leg(0.55, 0.74, 0.88, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(hIdealRes, "Ideal", "l");
        leg.AddEntry(hDetRes, "Detector effects", "l");
        leg.Draw();
        Label("Ideal vs detector resolution", "silicon pitch + energy + angular smearing");
        c.SaveAs(Form("%s/thetaee_resolution_comparison.pdf", outdir));
        c.SaveAs(Form("%s/thetaee_resolution.pdf", outdir));
    }

    // Reco vs gen ideal.
    {
        c.Clear();
        c.SetRightMargin(0.15);
        hIdealVsGen->GetXaxis()->SetTitle("#theta_{ee}^{gen} [deg]");
        hIdealVsGen->GetYaxis()->SetTitle("#theta_{ee}^{ideal} [deg]");
        hIdealVsGen->Draw("COLZ");
        TLine diag(0.0, 0.0, 180.0, 180.0);
        diag.SetLineColor(kRed + 1);
        diag.SetLineStyle(2);
        diag.SetLineWidth(2);
        diag.Draw("SAME");
        Label("Ideal reconstruction", "diagonal = perfect truth reconstruction");
        c.SaveAs(Form("%s/thetaee_ideal_reco_vs_gen.pdf", outdir));
        c.SetRightMargin(0.12);
    }

    // Reco vs gen detector.
    {
        c.Clear();
        c.SetRightMargin(0.15);
        hDetVsGen->GetXaxis()->SetTitle("#theta_{ee}^{gen} [deg]");
        hDetVsGen->GetYaxis()->SetTitle("#theta_{ee}^{detector} [deg]");
        hDetVsGen->Draw("COLZ");
        TLine diag(0.0, 0.0, 180.0, 180.0);
        diag.SetLineColor(kRed + 1);
        diag.SetLineStyle(2);
        diag.SetLineWidth(2);
        diag.Draw("SAME");
        Label("Detector reconstruction", "silicon pitch + scintillator energy selection");
        c.SaveAs(Form("%s/thetaee_detector_reco_vs_gen.pdf", outdir));
        c.SaveAs(Form("%s/thetaee_reco_vs_gen.pdf", outdir));
        c.SetRightMargin(0.12);
    }

    // Acceptance.
    {
        TH1D* hAccGeom = static_cast<TH1D*>(hGeom->Clone("hAccGeom"));
        TH1D* hAccDet = static_cast<TH1D*>(hDetAcceptedGen->Clone("hAccDet"));
        TH1D* hDetGivenGeom = static_cast<TH1D*>(hDetAcceptedGen->Clone("hDetGivenGeom"));
        Detach(hAccGeom); Detach(hAccDet); Detach(hDetGivenGeom);

        hAccGeom->Divide(hGeom, hGen, 1.0, 1.0, "B");
        hAccDet->Divide(hDetAcceptedGen, hGen, 1.0, 1.0, "B");
        hDetGivenGeom->Divide(hDetAcceptedGen, hGeom, 1.0, 1.0, "B");

        hAccGeom->SetLineColor(kBlue + 1);
        hAccGeom->SetLineWidth(3);
        hAccGeom->SetMarkerColor(kBlue + 1);
        hAccGeom->SetMarkerStyle(20);
        hAccGeom->SetMarkerSize(0.85);

        hAccDet->SetLineColor(kGreen + 2);
        hAccDet->SetLineWidth(3);
        hAccDet->SetMarkerColor(kGreen + 2);
        hAccDet->SetMarkerStyle(24);
        hAccDet->SetMarkerSize(0.95);

        c.Clear();
        TH1D frame("hAccFrame", "", 90, 0.0, 180.0);
        frame.SetDirectory(nullptr);
        frame.GetXaxis()->SetTitle("#theta_{ee}^{gen} [deg]");
        frame.GetYaxis()->SetTitle("Accepted / generated");
        frame.SetMinimum(0.0);
        frame.SetMaximum(1.20);
        frame.Draw("AXIS");

        TLine idealLine(0.0, 1.0, 180.0, 1.0);
        idealLine.SetLineColor(kRed + 1);
        idealLine.SetLineStyle(2);
        idealLine.SetLineWidth(3);
        idealLine.Draw("SAME");

        hAccGeom->Draw("E1 SAME");
        hAccDet->Draw("E1 SAME");

        TLegend leg(0.53, 0.70, 0.88, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(&idealLine, "Ideal = 1", "l");
        leg.AddEntry(hAccGeom, "Geometrical", "lep");
        leg.AddEntry(hAccDet, "Detector effects", "lep");
        leg.Draw();
        Label("Acceptance vs opening angle", "geometry and detector effects shown separately");
        c.SaveAs(Form("%s/acceptance_ideal_geometry_detector_vs_thetaee.pdf", outdir));
        c.SaveAs(Form("%s/acceptance_vs_thetaee.pdf", outdir));

        c.Clear();
        hDetGivenGeom->SetLineColor(kBlack);
        hDetGivenGeom->SetLineWidth(3);
        hDetGivenGeom->SetMarkerColor(kBlack);
        hDetGivenGeom->SetMarkerStyle(20);
        hDetGivenGeom->SetMarkerSize(0.85);
        hDetGivenGeom->GetXaxis()->SetTitle("#theta_{ee}^{gen} [deg]");
        hDetGivenGeom->GetYaxis()->SetTitle("Detector / geometrical");
        hDetGivenGeom->SetMinimum(0.0);
        hDetGivenGeom->SetMaximum(1.20);
        hDetGivenGeom->Draw("E1");

        TLine one(0.0, 1.0, 180.0, 1.0);
        one.SetLineColor(kRed + 1);
        one.SetLineStyle(2);
        one.SetLineWidth(2);
        one.Draw("SAME");
        Label("Detector efficiency after geometry", "ratio: detector effects / geometrical acceptance");
        c.SaveAs(Form("%s/detector_efficiency_given_geometry_vs_thetaee.pdf", outdir));

        delete hAccGeom;
        delete hAccDet;
        delete hDetGivenGeom;
    }

    Save1D(c, hCutflow,
           "", "Events",
           "Analysis cutflow", "generated #rightarrow ideal #rightarrow geometry #rightarrow detector",
           Form("%s/cutflow.pdf", outdir), true);

    // Scintillator energies.
    {
        c.Clear();
        const double maxy = 1.35 * std::max(hEem->GetMaximum(), hEep->GetMaximum());
        hEem->GetXaxis()->SetTitle("E_{scint}^{reco} [MeV]");
        hEem->GetYaxis()->SetTitle("Detector-effect events");
        hEem->SetMinimum(0.0);
        hEem->SetMaximum(maxy > 0.0 ? maxy : 1.0);
        hEem->Draw("HIST");
        hEep->Draw("HIST SAME");
        TLegend leg(0.58, 0.76, 0.88, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(hEem, "e^{-}", "l");
        leg.AddEntry(hEep, "e^{+}", "l");
        leg.Draw();
        Label("Scintillator reconstructed energy", "5% Gaussian smearing");
        c.SaveAs(Form("%s/scint_energy_em_ep.pdf", outdir));
    }

    Save1D(c, hEsum,
           "E_{e^{-}}^{reco} + E_{e^{+}}^{reco} [MeV]",
           "Detector-effect events", "Scintillator energy sum", "5% Gaussian smearing",
           Form("%s/scint_energy_sum.pdf", outdir));

    Save1D(c, hY,
           "Y = (E_{e^{-}} - E_{e^{+}})/(E_{e^{-}} + E_{e^{+}})",
           "Detector-effect events", "Energy asymmetry", "from reconstructed energy",
           Form("%s/scint_energy_asymmetry.pdf", outdir));

    {
        c.Clear();
        c.SetRightMargin(0.15);
        hThetaVsEsum->GetXaxis()->SetTitle("#theta_{ee}^{detector} [deg]");
        hThetaVsEsum->GetYaxis()->SetTitle("E_{sum}^{reco} [MeV]");
        hThetaVsEsum->Draw("COLZ");
        Label("#theta_{ee} vs scintillator energy", "detector-effect events");
        c.SaveAs(Form("%s/thetaee_vs_Esum_reco.pdf", outdir));
        c.SetRightMargin(0.12);
    }

    {
        c.Clear();
        c.SetRightMargin(0.15);
        hThetaVsY->GetXaxis()->SetTitle("#theta_{ee}^{detector} [deg]");
        hThetaVsY->GetYaxis()->SetTitle("Y_{reco}");
        hThetaVsY->Draw("COLZ");
        Label("#theta_{ee} vs energy asymmetry", "Y=(E_{e^{-}}-E_{e^{+}})/(E_{e^{-}}+E_{e^{+}})");
        c.SaveAs(Form("%s/thetaee_vs_Y_reco.pdf", outdir));
        c.SetRightMargin(0.12);
    }

    // Raw occupancy: every hit row in the hits tree. This is the safest
    // diagnostic for checking detector-ID/copy-number mapping.
    {
        c.Clear();
        const double maxy = 1.35 * std::max(hOccSiRawRows->GetMaximum(), hOccScintRawRows->GetMaximum());
        hOccSiRawRows->GetXaxis()->SetTitle("Detector ID");
        hOccSiRawRows->GetYaxis()->SetTitle("Primary hit rows");
        hOccSiRawRows->SetMinimum(0.0);
        hOccSiRawRows->SetMaximum(maxy > 0.0 ? maxy : 1.0);
        hOccSiRawRows->Draw("HIST");
        hOccScintRawRows->Draw("HIST SAME");
        TLegend leg(0.58, 0.76, 0.88, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(hOccSiRawRows, "Silicon", "l");
        leg.AddEntry(hOccScintRawRows, "Scintillator", "l");
        leg.Draw();
        Label("Raw detector occupancy", "all primary e^{-}/e^{+} hit rows");
        c.SaveAs(Form("%s/hit_detector_occupancy_raw.pdf", outdir));
    }

    // Event-level occupancy: each event contributes at most once per detector.
    {
        c.Clear();
        const double maxy = 1.35 * std::max(hOccSiEventLevel->GetMaximum(), hOccScintEventLevel->GetMaximum());
        hOccSiEventLevel->GetXaxis()->SetTitle("Detector ID");
        hOccSiEventLevel->GetYaxis()->SetTitle("Events with at least one hit");
        hOccSiEventLevel->SetMinimum(0.0);
        hOccSiEventLevel->SetMaximum(maxy > 0.0 ? maxy : 1.0);
        hOccSiEventLevel->Draw("HIST");
        hOccScintEventLevel->Draw("HIST SAME");
        TLegend leg(0.58, 0.76, 0.88, 0.88);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(hOccSiEventLevel, "Silicon", "l");
        leg.AddEntry(hOccScintEventLevel, "Scintillator", "l");
        leg.Draw();
        Label("Event-level detector occupancy", "one count per event and detector");
        c.SaveAs(Form("%s/hit_detector_occupancy_eventlevel.pdf", outdir));
    }

    // Reconstruction-used occupancy: only the detector IDs selected as e-/e+ positions.
    Save1D(c, hOccSiUsed,
           "Detector ID", "Lepton positions used in reconstruction",
           "Reconstruction-used silicon occupancy", "detector IDs selected for e^{-}/e^{+} positions",
           Form("%s/hit_detector_occupancy.pdf", outdir));
    c.SaveAs(Form("%s/hit_detector_occupancy_reconstruction_used.pdf", outdir));

    Save1D(c, hHitVol,
           "", "Primary hit rows",
           "Hit volume breakdown", "primary e^{-}/e^{+} hits",
           Form("%s/hit_volume_breakdown.pdf", outdir));

    fin->Close();

    std::cout << "\n[analyze_x17] Wrote ROOT file: x17_analysis.root" << std::endl;
    std::cout << "[analyze_x17] Wrote PDFs in: " << outdir << std::endl;
}

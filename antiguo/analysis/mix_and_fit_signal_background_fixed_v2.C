// mix_and_fit_signal_background.C
//
// Combine reconstructed signal and background samples into a single pseudo-data sample
// and extract the signal/background composition with a template fit.
//
// Main idea:
//   1. Read signal detected_detector_effects.
//   2. Read background detected_detector_effects.
//   3. Build one mixed tree called "mixed_events".
//   4. Fit the mixed thetaEE_reco_deg distribution with signal and background templates.
//   5. Estimate N_signal, N_background, purity and contamination.
//   6. Save ROOT histograms and PDF plots.
//
// Usage:
//
//   root -l -q 'analysis/mix_and_fit_signal_background.C(
//      "analysis_hits_detector_effects.root",
//      "background_analysis_hits_detector_effects.root",
//      "signal_background_fit.root",
//      "plots_signal_background_fit",
//      -1,
//      -1,
//      120.0,
//      180.0,
//      12345
//   )'
//
// Arguments:
//   sigRecoFile      : signal ROOT file containing tree detected_detector_effects
//   bkgRecoFile      : background ROOT file containing tree detected_detector_effects
//   outRootFile      : output ROOT file
//   outPlotDir       : directory for PDF plots
//   nSignalToMix     : number of signal events injected into pseudo-data; -1 = all signal entries
//   nBackgroundToMix : number of background events injected into pseudo-data; -1 = all background entries
//   thetaRegionMin   : lower edge of signal-region window in thetaEE_reco_deg
//   thetaRegionMax   : upper edge of signal-region window in thetaEE_reco_deg
//   seed             : random seed for pseudo-data sampling
//
// Notes:
//   - The fit is a template fit in thetaEE_reco_deg using ROOT TFractionFitter.
//   - The mixed tree keeps a MC truth label only for validation:
//         isSignalTruth = 1 for injected signal
//         isSignalTruth = 0 for injected background
//     In real data this label would not exist.
//   - The contamination is reported as:
//         contamination = N_background / (N_signal + N_background)
//     both globally and inside the theta signal-region window.

#include <TFile.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TFractionFitter.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>

#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <random>

namespace sbfit {

constexpr double kElectronMassMeV = 0.51099895;

struct EventRow {
    int eventID = -1;
    int dataRow = -1;
    int isSignalTruth = -1;

    double thetaEE_reco_deg = -999.0;
    double thetaEe_reco_deg = -999.0;
    double thetaEp_reco_deg = -999.0;

    double energyEe_reco_MeV = -999.0;  // kinetic energy
    double energyEp_reco_MeV = -999.0;  // kinetic energy
    double energyEE_reco_MeV = -999.0;  // kinetic sum

    double massEE_reco_MeV = -999.0;
};

bool HasBranch(TTree* t, const char* name) {
    return t && t->GetBranch(name);
}

double SafeMassFromKinematics(double k1, double k2, double thetaDeg) {
    if (!std::isfinite(k1) || !std::isfinite(k2) || !std::isfinite(thetaDeg)) return -999.0;
    if (k1 < 0.0 || k2 < 0.0) return -999.0;

    const double e1 = k1 + kElectronMassMeV;
    const double e2 = k2 + kElectronMassMeV;

    const double p1sq = std::max(0.0, e1 * e1 - kElectronMassMeV * kElectronMassMeV);
    const double p2sq = std::max(0.0, e2 * e2 - kElectronMassMeV * kElectronMassMeV);

    const double p1 = std::sqrt(p1sq);
    const double p2 = std::sqrt(p2sq);
    const double c  = std::cos(thetaDeg * TMath::DegToRad());

    const double m2 = 2.0 * kElectronMassMeV * kElectronMassMeV
                    + 2.0 * (e1 * e2 - p1 * p2 * c);

    return (m2 > 0.0) ? std::sqrt(m2) : 0.0;
}

std::vector<EventRow> ReadRecoTree(TTree* tree, int truthLabel) {
    std::vector<EventRow> rows;

    if (!tree) {
        std::cerr << "[ERROR] Null tree passed to ReadRecoTree.\n";
        return rows;
    }

    if (!HasBranch(tree, "thetaEE_reco_deg")) {
        std::cerr << "[ERROR] Tree '" << tree->GetName()
                  << "' does not contain thetaEE_reco_deg.\n";
        return rows;
    }

    const bool hasEventID = HasBranch(tree, "eventID");
    const bool hasDataRow = HasBranch(tree, "dataRow");

    const bool hasThetaEe = HasBranch(tree, "thetaEe_reco_deg");
    const bool hasThetaEp = HasBranch(tree, "thetaEp_reco_deg");

    const bool hasEnergyEe = HasBranch(tree, "energyEe_reco_MeV");
    const bool hasEnergyEp = HasBranch(tree, "energyEp_reco_MeV");
    const bool hasEnergyEE = HasBranch(tree, "energyEE_reco_MeV");

    const bool hasMassEE = HasBranch(tree, "massEE_reco_MeV");

    TTreeFormula fThetaEE("fThetaEE", "thetaEE_reco_deg", tree);

    std::unique_ptr<TTreeFormula> fEventID;
    std::unique_ptr<TTreeFormula> fDataRow;
    std::unique_ptr<TTreeFormula> fThetaEe;
    std::unique_ptr<TTreeFormula> fThetaEp;
    std::unique_ptr<TTreeFormula> fEnergyEe;
    std::unique_ptr<TTreeFormula> fEnergyEp;
    std::unique_ptr<TTreeFormula> fEnergyEE;
    std::unique_ptr<TTreeFormula> fMassEE;

    if (hasEventID) fEventID = std::make_unique<TTreeFormula>("fEventID", "eventID", tree);
    if (hasDataRow) fDataRow = std::make_unique<TTreeFormula>("fDataRow", "dataRow", tree);
    if (hasThetaEe) fThetaEe = std::make_unique<TTreeFormula>("fThetaEe", "thetaEe_reco_deg", tree);
    if (hasThetaEp) fThetaEp = std::make_unique<TTreeFormula>("fThetaEp", "thetaEp_reco_deg", tree);

    if (hasEnergyEe) fEnergyEe = std::make_unique<TTreeFormula>("fEnergyEe", "energyEe_reco_MeV", tree);
    if (hasEnergyEp) fEnergyEp = std::make_unique<TTreeFormula>("fEnergyEp", "energyEp_reco_MeV", tree);
    if (hasEnergyEE) fEnergyEE = std::make_unique<TTreeFormula>("fEnergyEE", "energyEE_reco_MeV", tree);

    if (hasMassEE) fMassEE = std::make_unique<TTreeFormula>("fMassEE", "massEE_reco_MeV", tree);

    const Long64_t n = tree->GetEntries();
    rows.reserve(n);

    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);

        EventRow r;
        r.isSignalTruth = truthLabel;
        r.eventID = hasEventID ? static_cast<int>(fEventID->EvalInstance()) : static_cast<int>(i);
        r.dataRow = hasDataRow ? static_cast<int>(fDataRow->EvalInstance()) : -1;

        r.thetaEE_reco_deg = fThetaEE.EvalInstance();
        r.thetaEe_reco_deg = hasThetaEe ? fThetaEe->EvalInstance() : -999.0;
        r.thetaEp_reco_deg = hasThetaEp ? fThetaEp->EvalInstance() : -999.0;

        r.energyEe_reco_MeV = hasEnergyEe ? fEnergyEe->EvalInstance() : -999.0;
        r.energyEp_reco_MeV = hasEnergyEp ? fEnergyEp->EvalInstance() : -999.0;

        if (hasEnergyEE) {
            r.energyEE_reco_MeV = fEnergyEE->EvalInstance();
        } else if (hasEnergyEe && hasEnergyEp) {
            r.energyEE_reco_MeV = r.energyEe_reco_MeV + r.energyEp_reco_MeV;
        }

        if (hasMassEE) {
            r.massEE_reco_MeV = fMassEE->EvalInstance();
        } else if (hasEnergyEe && hasEnergyEp) {
            r.massEE_reco_MeV = SafeMassFromKinematics(
                r.energyEe_reco_MeV,
                r.energyEp_reco_MeV,
                r.thetaEE_reco_deg
            );
        }

        if (!std::isfinite(r.thetaEE_reco_deg)) continue;
        rows.push_back(r);
    }

    return rows;
}

void FillHistograms(const std::vector<EventRow>& rows,
                    TH1D* hTheta,
                    TH1D* hEnergy,
                    TH1D* hMass) {
    for (const auto& r : rows) {
        if (hTheta  && std::isfinite(r.thetaEE_reco_deg)  && r.thetaEE_reco_deg > -900) hTheta->Fill(r.thetaEE_reco_deg);
        if (hEnergy && std::isfinite(r.energyEE_reco_MeV) && r.energyEE_reco_MeV > -900) hEnergy->Fill(r.energyEE_reco_MeV);
        if (hMass   && std::isfinite(r.massEE_reco_MeV)   && r.massEE_reco_MeV > -900) hMass->Fill(r.massEE_reco_MeV);
    }
}

std::vector<EventRow> SampleRows(const std::vector<EventRow>& input,
                                 Long64_t requested,
                                 TRandom3& rng) {
    std::vector<EventRow> out;

    if (input.empty()) return out;

    if (requested < 0) {
        return input;
    }

    out.reserve(requested);
    for (Long64_t i = 0; i < requested; ++i) {
        const int idx = rng.Integer(static_cast<UInt_t>(input.size()));
        out.push_back(input[idx]);
    }

    return out;
}

double IntegralInRange(TH1D* h, double xmin, double xmax) {
    if (!h) return 0.0;

    const int b1 = h->GetXaxis()->FindFixBin(xmin + 1e-9);
    const int b2 = h->GetXaxis()->FindFixBin(xmax - 1e-9);

    const int lo = std::max(1, b1);
    const int hi = std::min(h->GetNbinsX(), b2);

    if (hi < lo) return 0.0;
    return h->Integral(lo, hi);
}

TH1D* CloneNormalized(TH1D* h, const char* name) {
    if (!h) return nullptr;
    TH1D* c = dynamic_cast<TH1D*>(h->Clone(name));
    c->SetDirectory(nullptr);
    const double integral = c->Integral();
    if (integral > 0) c->Scale(1.0 / integral);
    return c;
}

TH1D* ScaledTemplate(TH1D* hTemplate, const char* name, double yield) {
    TH1D* h = CloneNormalized(hTemplate, name);
    if (h) h->Scale(yield);
    return h;
}

void StyleHist(TH1D* h, int color, int width = 3) {
    if (!h) return;
    h->SetLineColor(color);
    h->SetLineWidth(width);
    h->SetMarkerColor(color);
}

void SaveTemplateFitPlot(TH1D* hData,
                         TH1D* hSigFit,
                         TH1D* hBkgFit,
                         TH1D* hTotalFit,
                         const std::string& outPdf,
                         const std::string& xTitle,
                         const std::string& canvasName,
                         double thetaRegionMin = -999,
                         double thetaRegionMax = -999) {
    if (!hData || !hSigFit || !hBkgFit || !hTotalFit) return;

    TCanvas c(canvasName.c_str(), canvasName.c_str(), 950, 720);
    c.SetTicks(1, 1);

    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.9);
    hData->SetLineColor(kBlack);
    hData->SetMarkerColor(kBlack);

    StyleHist(hSigFit, kBlue + 1, 3);
    StyleHist(hBkgFit, kRed + 1, 3);
    StyleHist(hTotalFit, kBlack, 4);

    hData->GetXaxis()->SetTitle(xTitle.c_str());
    hData->GetYaxis()->SetTitle("events");

    double ymax = hData->GetMaximum();
    ymax = std::max(ymax, hSigFit->GetMaximum());
    ymax = std::max(ymax, hBkgFit->GetMaximum());
    ymax = std::max(ymax, hTotalFit->GetMaximum());
    hData->SetMaximum(ymax > 0 ? 1.35 * ymax : 1.0);

    hData->Draw("E");
    hBkgFit->Draw("HIST SAME");
    hSigFit->Draw("HIST SAME");
    hTotalFit->Draw("HIST SAME");
    hData->Draw("E SAME");

    if (thetaRegionMin > -900 && thetaRegionMax > -900) {
        TLine l1(thetaRegionMin, 0, thetaRegionMin, hData->GetMaximum());
        TLine l2(thetaRegionMax, 0, thetaRegionMax, hData->GetMaximum());
        l1.SetLineStyle(2);
        l2.SetLineStyle(2);
        l1.SetLineWidth(2);
        l2.SetLineWidth(2);
        l1.Draw("SAME");
        l2.Draw("SAME");
    }

    TLegend leg(0.58, 0.68, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(hData, "mixed pseudo-data", "pe");
    leg.AddEntry(hSigFit, "signal template fit", "l");
    leg.AddEntry(hBkgFit, "background template fit", "l");
    leg.AddEntry(hTotalFit, "signal + background", "l");
    leg.Draw();

    c.SaveAs(outPdf.c_str());
}

void SaveOverlayNormalized(TH1D* hSig,
                           TH1D* hBkg,
                           TH1D* hData,
                           const std::string& outPdf,
                           const std::string& xTitle,
                           const std::string& canvasName) {
    if (!hSig || !hBkg || !hData) return;

    std::unique_ptr<TH1D> s(CloneNormalized(hSig, "tmp_s_norm"));
    std::unique_ptr<TH1D> b(CloneNormalized(hBkg, "tmp_b_norm"));
    std::unique_ptr<TH1D> d(CloneNormalized(hData, "tmp_d_norm"));

    TCanvas c(canvasName.c_str(), canvasName.c_str(), 950, 720);
    c.SetTicks(1, 1);

    StyleHist(s.get(), kBlue + 1, 3);
    StyleHist(b.get(), kRed + 1, 3);
    StyleHist(d.get(), kBlack, 4);

    double ymax = std::max({s->GetMaximum(), b->GetMaximum(), d->GetMaximum()});
    d->SetMaximum(ymax > 0 ? 1.35 * ymax : 1.0);
    d->GetXaxis()->SetTitle(xTitle.c_str());
    d->GetYaxis()->SetTitle("normalized entries");

    d->Draw("HIST");
    b->Draw("HIST SAME");
    s->Draw("HIST SAME");

    TLegend leg(0.58, 0.70, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(d.get(), "mixed", "l");
    leg.AddEntry(s.get(), "signal template", "l");
    leg.AddEntry(b.get(), "background template", "l");
    leg.Draw();

    c.SaveAs(outPdf.c_str());
}

} // namespace sbfit


void mix_and_fit_signal_background_fixed_v2(
    const char* sigRecoFile = "analysis_hits_detector_effects.root",
    const char* bkgRecoFile = "background_analysis_hits_detector_effects.root",
    const char* outRootFile = "signal_background_fit.root",
    const char* outPlotDir = "plots_signal_background_fit",
    Long64_t nSignalToMix = -1,
    Long64_t nBackgroundToMix = -1,
    double thetaRegionMin = 120.0,
    double thetaRegionMax = 180.0,
    unsigned int seed = 12345
) {
    using namespace sbfit;

    gStyle->SetOptStat(0);
    gSystem->mkdir(outPlotDir, kTRUE);

    std::unique_ptr<TFile> fSig(TFile::Open(sigRecoFile, "READ"));
    std::unique_ptr<TFile> fBkg(TFile::Open(bkgRecoFile, "READ"));

    if (!fSig || fSig->IsZombie()) {
        std::cerr << "[ERROR] Cannot open signal file: " << sigRecoFile << "\n";
        return;
    }
    if (!fBkg || fBkg->IsZombie()) {
        std::cerr << "[ERROR] Cannot open background file: " << bkgRecoFile << "\n";
        return;
    }

    TTree* tSig = dynamic_cast<TTree*>(fSig->Get("detected_detector_effects"));
    TTree* tBkg = dynamic_cast<TTree*>(fBkg->Get("detected_detector_effects"));

    if (!tSig) {
        std::cerr << "[ERROR] Missing tree detected_detector_effects in " << sigRecoFile << "\n";
        return;
    }
    if (!tBkg) {
        std::cerr << "[ERROR] Missing tree detected_detector_effects in " << bkgRecoFile << "\n";
        return;
    }

    std::vector<EventRow> sigRowsAll = ReadRecoTree(tSig, 1);
    std::vector<EventRow> bkgRowsAll = ReadRecoTree(tBkg, 0);

    std::cout << "\n=== Signal/background mixing and template fit ===\n";
    std::cout << "Signal detector-effects file     : " << sigRecoFile << "\n";
    std::cout << "Background detector-effects file : " << bkgRecoFile << "\n";
    std::cout << "Signal entries available         : " << sigRowsAll.size() << "\n";
    std::cout << "Background entries available     : " << bkgRowsAll.size() << "\n";
    std::cout << "Requested signal events          : " << nSignalToMix << "  (-1 = all)\n";
    std::cout << "Requested background events      : " << nBackgroundToMix << "  (-1 = all)\n";
    std::cout << "Theta signal region              : [" << thetaRegionMin << ", "
              << thetaRegionMax << "] deg\n";

    if (sigRowsAll.empty() || bkgRowsAll.empty()) {
        std::cerr << "[ERROR] Need non-empty signal and background detector-effects trees.\n";
        std::cerr << "        If one is empty, rerun with more generated events or looser selections.\n";
        return;
    }

    TRandom3 rng(seed);

    std::vector<EventRow> sigMixed = SampleRows(sigRowsAll, nSignalToMix, rng);
    std::vector<EventRow> bkgMixed = SampleRows(bkgRowsAll, nBackgroundToMix, rng);

    std::vector<EventRow> mixed;
    mixed.reserve(sigMixed.size() + bkgMixed.size());
    mixed.insert(mixed.end(), sigMixed.begin(), sigMixed.end());
    mixed.insert(mixed.end(), bkgMixed.begin(), bkgMixed.end());

    // Shuffle pseudo-data so the truth origin is not ordered in the output tree.
    std::mt19937 shuffleEngine(seed);
    std::shuffle(mixed.begin(), mixed.end(), shuffleEngine);

    std::cout << "Injected signal events           : " << sigMixed.size() << "\n";
    std::cout << "Injected background events       : " << bkgMixed.size() << "\n";
    std::cout << "Total mixed pseudo-data events   : " << mixed.size() << "\n";

    std::unique_ptr<TFile> fout(TFile::Open(outRootFile, "RECREATE"));
    if (!fout || fout->IsZombie()) {
        std::cerr << "[ERROR] Cannot create output file: " << outRootFile << "\n";
        return;
    }

    // Output mixed tree.
    TTree tMixed("mixed_events", "Mixed signal+background pseudo-data after detector effects");

    EventRow out;
    tMixed.Branch("eventID", &out.eventID, "eventID/I");
    tMixed.Branch("dataRow", &out.dataRow, "dataRow/I");
    tMixed.Branch("isSignalTruth", &out.isSignalTruth, "isSignalTruth/I");

    tMixed.Branch("thetaEE_reco_deg", &out.thetaEE_reco_deg, "thetaEE_reco_deg/D");
    tMixed.Branch("thetaEe_reco_deg", &out.thetaEe_reco_deg, "thetaEe_reco_deg/D");
    tMixed.Branch("thetaEp_reco_deg", &out.thetaEp_reco_deg, "thetaEp_reco_deg/D");

    tMixed.Branch("energyEe_reco_MeV", &out.energyEe_reco_MeV, "energyEe_reco_MeV/D");
    tMixed.Branch("energyEp_reco_MeV", &out.energyEp_reco_MeV, "energyEp_reco_MeV/D");
    tMixed.Branch("energyEE_reco_MeV", &out.energyEE_reco_MeV, "energyEE_reco_MeV/D");

    tMixed.Branch("massEE_reco_MeV", &out.massEE_reco_MeV, "massEE_reco_MeV/D");

    for (const auto& r : mixed) {
        out = r;
        tMixed.Fill();
    }

    // Histograms.
    TH1D hSigTheta("hSigTheta_template", "Signal template;#theta_{ee}^{reco} [deg];events", 90, 0, 180);
    TH1D hBkgTheta("hBkgTheta_template", "Background template;#theta_{ee}^{reco} [deg];events", 90, 0, 180);
    TH1D hDataTheta("hDataTheta_mixed", "Mixed pseudo-data;#theta_{ee}^{reco} [deg];events", 90, 0, 180);

    TH1D hSigEnergy("hSigEnergy_template", "Signal template;E_{ee}^{reco} [MeV];events", 80, 0, 25);
    TH1D hBkgEnergy("hBkgEnergy_template", "Background template;E_{ee}^{reco} [MeV];events", 80, 0, 25);
    TH1D hDataEnergy("hDataEnergy_mixed", "Mixed pseudo-data;E_{ee}^{reco} [MeV];events", 80, 0, 25);

    TH1D hSigMass("hSigMass_template", "Signal template;m_{ee}^{reco} [MeV/c^{2}];events", 80, 0, 25);
    TH1D hBkgMass("hBkgMass_template", "Background template;m_{ee}^{reco} [MeV/c^{2}];events", 80, 0, 25);
    TH1D hDataMass("hDataMass_mixed", "Mixed pseudo-data;m_{ee}^{reco} [MeV/c^{2}];events", 80, 0, 25);

    hSigTheta.Sumw2(); hBkgTheta.Sumw2(); hDataTheta.Sumw2();
    hSigEnergy.Sumw2(); hBkgEnergy.Sumw2(); hDataEnergy.Sumw2();
    hSigMass.Sumw2(); hBkgMass.Sumw2(); hDataMass.Sumw2();

    // Important in ROOT macros:
    // the output TFile is already open at this point, so newly created histograms
    // would otherwise be owned by that file/directory. If the file is closed before
    // the PDF plots are produced, those histograms can become dangling pointers and
    // ROOT may segfault inside TCanvas. Detach them explicitly.
    hSigTheta.SetDirectory(nullptr);
    hBkgTheta.SetDirectory(nullptr);
    hDataTheta.SetDirectory(nullptr);

    hSigEnergy.SetDirectory(nullptr);
    hBkgEnergy.SetDirectory(nullptr);
    hDataEnergy.SetDirectory(nullptr);

    hSigMass.SetDirectory(nullptr);
    hBkgMass.SetDirectory(nullptr);
    hDataMass.SetDirectory(nullptr);

    FillHistograms(sigRowsAll, &hSigTheta, &hSigEnergy, &hSigMass);
    FillHistograms(bkgRowsAll, &hBkgTheta, &hBkgEnergy, &hBkgMass);
    FillHistograms(mixed,     &hDataTheta, &hDataEnergy, &hDataMass);

    // Template fraction fit in theta.
    TObjArray templates(2);
    templates.Add(&hSigTheta);
    templates.Add(&hBkgTheta);

    TFractionFitter fitter(&hDataTheta, &templates);
    fitter.Constrain(0, 0.0, 1.0);
    fitter.Constrain(1, 0.0, 1.0);

    const int fitStatus = fitter.Fit();

    double fracSig = 0.0, fracSigErr = 0.0;
    double fracBkg = 0.0, fracBkgErr = 0.0;

    if (fitStatus == 0) {
        fitter.GetResult(0, fracSig, fracSigErr);
        fitter.GetResult(1, fracBkg, fracBkgErr);
    } else {
        std::cerr << "[WARNING] TFractionFitter returned status " << fitStatus
                  << ". Using truth fractions as fallback.\n";
        const double nTot = static_cast<double>(mixed.size());
        fracSig = nTot > 0 ? static_cast<double>(sigMixed.size()) / nTot : 0.0;
        fracBkg = nTot > 0 ? static_cast<double>(bkgMixed.size()) / nTot : 0.0;
    }

    const double nData = hDataTheta.Integral();
    const double nSigFit = fracSig * nData;
    const double nBkgFit = fracBkg * nData;

    std::unique_ptr<TH1D> hSigThetaFit(ScaledTemplate(&hSigTheta, "hSigTheta_fit", nSigFit));
    std::unique_ptr<TH1D> hBkgThetaFit(ScaledTemplate(&hBkgTheta, "hBkgTheta_fit", nBkgFit));
    std::unique_ptr<TH1D> hTotalThetaFit(dynamic_cast<TH1D*>(hSigThetaFit->Clone("hTotalTheta_fit")));
    hTotalThetaFit->Add(hBkgThetaFit.get());

    std::unique_ptr<TH1D> hSigEnergyFit(ScaledTemplate(&hSigEnergy, "hSigEnergy_fit", nSigFit));
    std::unique_ptr<TH1D> hBkgEnergyFit(ScaledTemplate(&hBkgEnergy, "hBkgEnergy_fit", nBkgFit));
    std::unique_ptr<TH1D> hTotalEnergyFit(dynamic_cast<TH1D*>(hSigEnergyFit->Clone("hTotalEnergy_fit")));
    hTotalEnergyFit->Add(hBkgEnergyFit.get());

    std::unique_ptr<TH1D> hSigMassFit(ScaledTemplate(&hSigMass, "hSigMass_fit", nSigFit));
    std::unique_ptr<TH1D> hBkgMassFit(ScaledTemplate(&hBkgMass, "hBkgMass_fit", nBkgFit));
    std::unique_ptr<TH1D> hTotalMassFit(dynamic_cast<TH1D*>(hSigMassFit->Clone("hTotalMass_fit")));
    hTotalMassFit->Add(hBkgMassFit.get());

    // Region estimates from fitted theta templates.
    const double sigRegionFit = IntegralInRange(hSigThetaFit.get(), thetaRegionMin, thetaRegionMax);
    const double bkgRegionFit = IntegralInRange(hBkgThetaFit.get(), thetaRegionMin, thetaRegionMax);
    const double totalRegionFit = sigRegionFit + bkgRegionFit;

    const double contaminationGlobal = (nSigFit + nBkgFit) > 0 ? nBkgFit / (nSigFit + nBkgFit) : 0.0;
    const double purityGlobal        = (nSigFit + nBkgFit) > 0 ? nSigFit / (nSigFit + nBkgFit) : 0.0;

    const double contaminationRegion = totalRegionFit > 0 ? bkgRegionFit / totalRegionFit : 0.0;
    const double purityRegion        = totalRegionFit > 0 ? sigRegionFit / totalRegionFit : 0.0;

    // Truth validation for pseudo-data.
    Long64_t truthSigRegion = 0;
    Long64_t truthBkgRegion = 0;
    for (const auto& r : mixed) {
        if (r.thetaEE_reco_deg >= thetaRegionMin && r.thetaEE_reco_deg < thetaRegionMax) {
            if (r.isSignalTruth == 1) ++truthSigRegion;
            if (r.isSignalTruth == 0) ++truthBkgRegion;
        }
    }

    const double truthContaminationGlobal =
        mixed.empty() ? 0.0 : static_cast<double>(bkgMixed.size()) / static_cast<double>(mixed.size());
    const double truthContaminationRegion =
        (truthSigRegion + truthBkgRegion) > 0
        ? static_cast<double>(truthBkgRegion) / static_cast<double>(truthSigRegion + truthBkgRegion)
        : 0.0;

    std::cout << "\n=== Template fit result using thetaEE_reco_deg ===\n";
    std::cout << "Fit status                       : " << fitStatus << "\n";
    std::cout << "Total pseudo-data events          : " << nData << "\n";
    std::cout << "Fitted signal fraction            : " << fracSig << " +/- " << fracSigErr << "\n";
    std::cout << "Fitted background fraction        : " << fracBkg << " +/- " << fracBkgErr << "\n";
    std::cout << "Estimated signal events           : " << nSigFit << "\n";
    std::cout << "Estimated background events       : " << nBkgFit << "\n";
    std::cout << "Global purity S/(S+B)             : " << purityGlobal << "\n";
    std::cout << "Global contamination B/(S+B)      : " << contaminationGlobal << "\n";

    std::cout << "\n=== Signal-region estimate ===\n";
    std::cout << "Theta region                      : [" << thetaRegionMin << ", " << thetaRegionMax << "] deg\n";
    std::cout << "Estimated signal in region        : " << sigRegionFit << "\n";
    std::cout << "Estimated background in region    : " << bkgRegionFit << "\n";
    std::cout << "Region purity S/(S+B)             : " << purityRegion << "\n";
    std::cout << "Region contamination B/(S+B)      : " << contaminationRegion << "\n";

    std::cout << "\n=== MC truth validation of mixed sample ===\n";
    std::cout << "Injected signal events            : " << sigMixed.size() << "\n";
    std::cout << "Injected background events        : " << bkgMixed.size() << "\n";
    std::cout << "Truth global contamination        : " << truthContaminationGlobal << "\n";
    std::cout << "Truth signal in theta region      : " << truthSigRegion << "\n";
    std::cout << "Truth background in theta region  : " << truthBkgRegion << "\n";
    std::cout << "Truth region contamination        : " << truthContaminationRegion << "\n";

    // Summary tree with scalar results.
    TTree tSummary("fit_summary", "Signal/background template-fit summary");

    int fitStatusOut = fitStatus;
    double nDataOut = nData;
    double nSigFitOut = nSigFit;
    double nBkgFitOut = nBkgFit;
    double fracSigOut = fracSig;
    double fracSigErrOut = fracSigErr;
    double fracBkgOut = fracBkg;
    double fracBkgErrOut = fracBkgErr;

    double purityGlobalOut = purityGlobal;
    double contaminationGlobalOut = contaminationGlobal;
    double purityRegionOut = purityRegion;
    double contaminationRegionOut = contaminationRegion;
    double sigRegionFitOut = sigRegionFit;
    double bkgRegionFitOut = bkgRegionFit;
    double thetaRegionMinOut = thetaRegionMin;
    double thetaRegionMaxOut = thetaRegionMax;
    double truthContaminationGlobalOut = truthContaminationGlobal;
    double truthContaminationRegionOut = truthContaminationRegion;

    tSummary.Branch("fitStatus", &fitStatusOut, "fitStatus/I");
    tSummary.Branch("nData", &nDataOut, "nData/D");
    tSummary.Branch("nSignalFit", &nSigFitOut, "nSignalFit/D");
    tSummary.Branch("nBackgroundFit", &nBkgFitOut, "nBackgroundFit/D");
    tSummary.Branch("fracSignalFit", &fracSigOut, "fracSignalFit/D");
    tSummary.Branch("fracSignalFitErr", &fracSigErrOut, "fracSignalFitErr/D");
    tSummary.Branch("fracBackgroundFit", &fracBkgOut, "fracBackgroundFit/D");
    tSummary.Branch("fracBackgroundFitErr", &fracBkgErrOut, "fracBackgroundFitErr/D");
    tSummary.Branch("purityGlobal", &purityGlobalOut, "purityGlobal/D");
    tSummary.Branch("contaminationGlobal", &contaminationGlobalOut, "contaminationGlobal/D");
    tSummary.Branch("thetaRegionMin", &thetaRegionMinOut, "thetaRegionMin/D");
    tSummary.Branch("thetaRegionMax", &thetaRegionMaxOut, "thetaRegionMax/D");
    tSummary.Branch("nSignalRegionFit", &sigRegionFitOut, "nSignalRegionFit/D");
    tSummary.Branch("nBackgroundRegionFit", &bkgRegionFitOut, "nBackgroundRegionFit/D");
    tSummary.Branch("purityRegion", &purityRegionOut, "purityRegion/D");
    tSummary.Branch("contaminationRegion", &contaminationRegionOut, "contaminationRegion/D");
    tSummary.Branch("truthContaminationGlobal", &truthContaminationGlobalOut, "truthContaminationGlobal/D");
    tSummary.Branch("truthContaminationRegion", &truthContaminationRegionOut, "truthContaminationRegion/D");
    tSummary.Fill();

    // Write output.
    fout->cd();

    tMixed.Write();
    tSummary.Write();

    hSigTheta.Write();
    hBkgTheta.Write();
    hDataTheta.Write();
    hSigThetaFit->Write();
    hBkgThetaFit->Write();
    hTotalThetaFit->Write();

    hSigEnergy.Write();
    hBkgEnergy.Write();
    hDataEnergy.Write();
    hSigEnergyFit->Write();
    hBkgEnergyFit->Write();
    hTotalEnergyFit->Write();

    hSigMass.Write();
    hBkgMass.Write();
    hDataMass.Write();
    hSigMassFit->Write();
    hBkgMassFit->Write();
    hTotalMassFit->Write();

    // Do not close fout yet. Produce the plots first, then close the file at the end.
    // Closing here can invalidate ROOT-owned objects in some interactive/cling sessions.
    fout->Write();

    const std::string dir(outPlotDir);

    SaveTemplateFitPlot(&hDataTheta, hSigThetaFit.get(), hBkgThetaFit.get(), hTotalThetaFit.get(),
                        dir + "/fit_thetaEE.pdf",
                        "#theta_{ee}^{reco} [deg]",
                        "cFitThetaEE",
                        thetaRegionMin,
                        thetaRegionMax);

    SaveTemplateFitPlot(&hDataEnergy, hSigEnergyFit.get(), hBkgEnergyFit.get(), hTotalEnergyFit.get(),
                        dir + "/fit_energyEE.pdf",
                        "E_{ee}^{reco} [MeV]",
                        "cFitEnergyEE");

    SaveTemplateFitPlot(&hDataMass, hSigMassFit.get(), hBkgMassFit.get(), hTotalMassFit.get(),
                        dir + "/fit_massEE.pdf",
                        "m_{ee}^{reco} [MeV/c^{2}]",
                        "cFitMassEE");

    SaveOverlayNormalized(&hSigTheta, &hBkgTheta, &hDataTheta,
                          dir + "/normalized_templates_thetaEE.pdf",
                          "#theta_{ee}^{reco} [deg]",
                          "cNormTheta");

    SaveOverlayNormalized(&hSigEnergy, &hBkgEnergy, &hDataEnergy,
                          dir + "/normalized_templates_energyEE.pdf",
                          "E_{ee}^{reco} [MeV]",
                          "cNormEnergy");

    SaveOverlayNormalized(&hSigMass, &hBkgMass, &hDataMass,
                          dir + "/normalized_templates_massEE.pdf",
                          "m_{ee}^{reco} [MeV/c^{2}]",
                          "cNormMass");

    fout->Close();

    std::cout << "\n=== Outputs ===\n";
    std::cout << "Output ROOT file : " << outRootFile << "\n";
    std::cout << "Output plot dir  : " << outPlotDir << "\n";
    std::cout << "Main tree        : mixed_events\n";
    std::cout << "Summary tree     : fit_summary\n";
    std::cout << "Main plots       : fit_thetaEE.pdf, fit_energyEE.pdf, fit_massEE.pdf\n";
    std::cout << "\nDone.\n";
}

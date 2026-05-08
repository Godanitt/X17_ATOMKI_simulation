// final_signal_background_analysis_fixed.C
//
// Robust final signal/background analysis.
// It reads:
//   analysis_hits_detector_effects.root
//   background_analysis_hits_detector_effects.root
//
// It expects the tree:
//   detected_detector_effects
//
// It compares signal/background in:
//   thetaEE_reco_deg
//   energyEE_reco_MeV
//   massEE_reco_MeV
//
// If energyEE_reco_MeV is missing, it computes:
//   energyEE_reco_MeV = energyEe_reco_MeV + energyEp_reco_MeV
//
// If massEE_reco_MeV is missing, it computes it from:
//   energyEe_reco_MeV, energyEp_reco_MeV, thetaEE_reco_deg
//
// Usage:
//   root -l -q 'analysis/final_signal_background_analysis_fixed.C(
//      "analysis_hits_detector_effects.root",
//      "background_analysis_hits_detector_effects.root",
//      "final_signal_background.root",
//      1.0,
//      1.0
//   )'

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>

#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

namespace fsb {

constexpr double me = 0.51099895; // MeV

bool HasBranch(TTree* t, const char* name) {
    return t && t->GetBranch(name);
}

double MassEEFromKinetic(double kEe, double kEp, double thetaDeg) {
    if (!std::isfinite(kEe) || !std::isfinite(kEp) || !std::isfinite(thetaDeg)) return -999.0;
    if (kEe < 0.0 || kEp < 0.0) return -999.0;

    const double EeTot = kEe + me;
    const double EpTot = kEp + me;

    const double pEe = std::sqrt(std::max(0.0, EeTot * EeTot - me * me));
    const double pEp = std::sqrt(std::max(0.0, EpTot * EpTot - me * me));

    const double c = std::cos(thetaDeg * TMath::DegToRad());
    const double m2 = 2.0 * me * me + 2.0 * (EeTot * EpTot - pEe * pEp * c);

    return m2 > 0.0 ? std::sqrt(m2) : 0.0;
}

struct BranchReader {
    bool ok = false;

    bool hasTheta = false;
    bool hasEnergyEe = false;
    bool hasEnergyEp = false;
    bool hasEnergyEE = false;
    bool hasMassEE = false;

    double thetaEE = -999.0;
    double energyEe = -999.0;
    double energyEp = -999.0;
    double energyEE = -999.0;
    double massEE = -999.0;
};

bool Attach(TTree* t, BranchReader& r, const std::string& label) {
    if (!t) {
        std::cerr << "[ERROR] Null tree for " << label << "\n";
        return false;
    }

    r.hasTheta    = HasBranch(t, "thetaEE_reco_deg");
    r.hasEnergyEe = HasBranch(t, "energyEe_reco_MeV");
    r.hasEnergyEp = HasBranch(t, "energyEp_reco_MeV");
    r.hasEnergyEE = HasBranch(t, "energyEE_reco_MeV");
    r.hasMassEE   = HasBranch(t, "massEE_reco_MeV");

    if (!r.hasTheta) {
        std::cerr << "[ERROR] " << label << " does not contain thetaEE_reco_deg\n";
        return false;
    }

    if (!r.hasEnergyEE && !(r.hasEnergyEe && r.hasEnergyEp)) {
        std::cerr << "[WARNING] " << label
                  << " has no energyEE_reco_MeV and cannot compute it because "
                  << "energyEe_reco_MeV/energyEp_reco_MeV are missing.\n";
    }

    if (!r.hasMassEE && !(r.hasEnergyEe && r.hasEnergyEp)) {
        std::cerr << "[WARNING] " << label
                  << " has no massEE_reco_MeV and cannot compute it because "
                  << "energyEe_reco_MeV/energyEp_reco_MeV are missing.\n";
    }

    t->SetBranchStatus("*", 0);

    t->SetBranchStatus("thetaEE_reco_deg", 1);
    t->SetBranchAddress("thetaEE_reco_deg", &r.thetaEE);

    if (r.hasEnergyEe) {
        t->SetBranchStatus("energyEe_reco_MeV", 1);
        t->SetBranchAddress("energyEe_reco_MeV", &r.energyEe);
    }

    if (r.hasEnergyEp) {
        t->SetBranchStatus("energyEp_reco_MeV", 1);
        t->SetBranchAddress("energyEp_reco_MeV", &r.energyEp);
    }

    if (r.hasEnergyEE) {
        t->SetBranchStatus("energyEE_reco_MeV", 1);
        t->SetBranchAddress("energyEE_reco_MeV", &r.energyEE);
    }

    if (r.hasMassEE) {
        t->SetBranchStatus("massEE_reco_MeV", 1);
        t->SetBranchAddress("massEE_reco_MeV", &r.massEE);
    }

    r.ok = true;
    return true;
}

Long64_t FillFromTree(TTree* t,
                      const std::string& label,
                      TH1D* hTheta,
                      TH1D* hEnergy,
                      TH1D* hMass,
                      double scale) {
    BranchReader r;
    if (!Attach(t, r, label)) return 0;

    Long64_t nFilled = 0;
    const Long64_t n = t->GetEntries();

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);

        const double theta = r.thetaEE;

        double energyEE = -999.0;
        if (r.hasEnergyEE) {
            energyEE = r.energyEE;
        } else if (r.hasEnergyEe && r.hasEnergyEp) {
            energyEE = r.energyEe + r.energyEp;
        }

        double massEE = -999.0;
        if (r.hasMassEE) {
            massEE = r.massEE;
        } else if (r.hasEnergyEe && r.hasEnergyEp) {
            massEE = MassEEFromKinetic(r.energyEe, r.energyEp, theta);
        }

        if (std::isfinite(theta) && theta > -900.0) {
            hTheta->Fill(theta, scale);
        }

        if (std::isfinite(energyEE) && energyEE > -900.0) {
            hEnergy->Fill(energyEE, scale);
        }

        if (std::isfinite(massEE) && massEE > -900.0) {
            hMass->Fill(massEE, scale);
        }

        ++nFilled;
    }

    t->SetBranchStatus("*", 1);
    return nFilled;
}

TH1D* SumHist(TH1D* s, TH1D* b, const char* name, const char* title) {
    TH1D* h = (TH1D*) b->Clone(name);
    h->SetTitle(title);
    h->Add(s);
    h->SetLineColor(kBlack);
    h->SetLineWidth(3);
    return h;
}

TH1D* SignificanceHist(TH1D* s, TH1D* b, const char* name, const char* title) {
    TH1D* h = (TH1D*) s->Clone(name);
    h->Reset("ICES");
    h->SetTitle(title);

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double S = s->GetBinContent(i);
        const double B = b->GetBinContent(i);
        h->SetBinContent(i, B > 0.0 ? S / std::sqrt(B) : 0.0);
        h->SetBinError(i, 0.0);
    }

    h->SetLineColor(kBlack);
    h->SetLineWidth(3);
    return h;
}

void Style(TH1D* sig, TH1D* bkg, TH1D* sum) {
    sig->SetLineColor(kBlue + 1);
    sig->SetLineWidth(3);

    bkg->SetLineColor(kRed + 1);
    bkg->SetLineWidth(3);

    sum->SetLineColor(kBlack);
    sum->SetLineWidth(3);
}

void SaveOverlay(TH1D* sig,
                 TH1D* bkg,
                 TH1D* sum,
                 const char* pdf,
                 const char* cname,
                 const char* xtitle) {
    TCanvas c(cname, cname, 900, 700);
    c.SetTicks(1, 1);

    Style(sig, bkg, sum);

    const double ymax = std::max({sig->GetMaximum(), bkg->GetMaximum(), sum->GetMaximum()});
    sum->SetMaximum(ymax > 0.0 ? 1.25 * ymax : 1.0);
    sum->GetXaxis()->SetTitle(xtitle);
    sum->GetYaxis()->SetTitle("entries");

    sum->Draw("HIST");
    bkg->Draw("HIST SAME");
    sig->Draw("HIST SAME");

    TLegend leg(0.60, 0.70, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(sig, "signal", "l");
    leg.AddEntry(bkg, "background", "l");
    leg.AddEntry(sum, "signal + background", "l");
    leg.Draw();

    c.SaveAs(pdf);
}

void SaveSingle(TH1D* h,
                const char* pdf,
                const char* cname,
                const char* xtitle) {
    TCanvas c(cname, cname, 900, 700);
    c.SetTicks(1, 1);

    h->SetLineColor(kBlack);
    h->SetLineWidth(3);
    h->GetXaxis()->SetTitle(xtitle);
    h->Draw("HIST");

    c.SaveAs(pdf);
}

} // namespace fsb


void final_signal_background_analysis_fixed(
    const char* signalRecoFile = "analysis_hits_detector_effects.root",
    const char* backgroundRecoFile = "background_analysis_hits_detector_effects.root",
    const char* outputFile = "final_signal_background.root",
    double signalScale = 1.0,
    double backgroundScale = 1.0
) {
    using namespace fsb;

    gStyle->SetOptStat(0);
    gSystem->mkdir("plots_final_signal_background", kTRUE);

    TFile* fSig = TFile::Open(signalRecoFile, "READ");
    TFile* fBkg = TFile::Open(backgroundRecoFile, "READ");

    if (!fSig || fSig->IsZombie()) {
        std::cerr << "[ERROR] Cannot open signal file: " << signalRecoFile << "\n";
        return;
    }

    if (!fBkg || fBkg->IsZombie()) {
        std::cerr << "[ERROR] Cannot open background file: " << backgroundRecoFile << "\n";
        return;
    }

    TTree* tSig = (TTree*) fSig->Get("detected_detector_effects");
    TTree* tBkg = (TTree*) fBkg->Get("detected_detector_effects");

    if (!tSig) {
        std::cerr << "[ERROR] Missing detected_detector_effects in " << signalRecoFile << "\n";
        return;
    }

    if (!tBkg) {
        std::cerr << "[ERROR] Missing detected_detector_effects in " << backgroundRecoFile << "\n";
        return;
    }

    TH1D* hThetaSig = new TH1D("hThetaEE_reco_signal",
                               "Reco opening angle;#theta_{ee}^{reco} [deg];entries",
                               90, 0, 180);
    TH1D* hThetaBkg = new TH1D("hThetaEE_reco_background",
                               "Reco opening angle;#theta_{ee}^{reco} [deg];entries",
                               90, 0, 180);

    TH1D* hEnergySig = new TH1D("hEnergyEE_reco_signal",
                                "Reco pair energy;E_{ee}^{reco} [MeV];entries",
                                80, 0, 25);
    TH1D* hEnergyBkg = new TH1D("hEnergyEE_reco_background",
                                "Reco pair energy;E_{ee}^{reco} [MeV];entries",
                                80, 0, 25);

    TH1D* hMassSig = new TH1D("hMassEE_reco_signal",
                              "Reco invariant mass;m_{ee}^{reco} [MeV/c^{2}];entries",
                              80, 0, 25);
    TH1D* hMassBkg = new TH1D("hMassEE_reco_background",
                              "Reco invariant mass;m_{ee}^{reco} [MeV/c^{2}];entries",
                              80, 0, 25);

    hThetaSig->Sumw2(); hThetaBkg->Sumw2();
    hEnergySig->Sumw2(); hEnergyBkg->Sumw2();
    hMassSig->Sumw2(); hMassBkg->Sumw2();

    const Long64_t nSigFilled = FillFromTree(tSig, "signal", hThetaSig, hEnergySig, hMassSig, signalScale);
    const Long64_t nBkgFilled = FillFromTree(tBkg, "background", hThetaBkg, hEnergyBkg, hMassBkg, backgroundScale);

    TH1D* hThetaSum = SumHist(hThetaSig, hThetaBkg,
                              "hThetaEE_reco_signalPlusBackground",
                              "Reco opening angle;#theta_{ee}^{reco} [deg];entries");

    TH1D* hEnergySum = SumHist(hEnergySig, hEnergyBkg,
                               "hEnergyEE_reco_signalPlusBackground",
                               "Reco pair energy;E_{ee}^{reco} [MeV];entries");

    TH1D* hMassSum = SumHist(hMassSig, hMassBkg,
                             "hMassEE_reco_signalPlusBackground",
                             "Reco invariant mass;m_{ee}^{reco} [MeV/c^{2}];entries");

    TH1D* hSigTheta = SignificanceHist(hThetaSig, hThetaBkg,
                                       "hSignificance_thetaEE",
                                       "Per-bin S/#sqrt{B};#theta_{ee}^{reco} [deg];S/#sqrt{B}");

    TH1D* hSigMass = SignificanceHist(hMassSig, hMassBkg,
                                      "hSignificance_massEE",
                                      "Per-bin S/#sqrt{B};m_{ee}^{reco} [MeV/c^{2}];S/#sqrt{B}");

    TFile* fout = TFile::Open(outputFile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "[ERROR] Cannot create output file: " << outputFile << "\n";
        return;
    }

    hThetaSig->Write();
    hThetaBkg->Write();
    hThetaSum->Write();

    hEnergySig->Write();
    hEnergyBkg->Write();
    hEnergySum->Write();

    hMassSig->Write();
    hMassBkg->Write();
    hMassSum->Write();

    hSigTheta->Write();
    hSigMass->Write();

    fout->Close();

    SaveOverlay(hThetaSig, hThetaBkg, hThetaSum,
                "plots_final_signal_background/thetaEE_signal_background.pdf",
                "cThetaEE",
                "#theta_{ee}^{reco} [deg]");

    SaveOverlay(hEnergySig, hEnergyBkg, hEnergySum,
                "plots_final_signal_background/energyEE_signal_background.pdf",
                "cEnergyEE",
                "E_{ee}^{reco} [MeV]");

    SaveOverlay(hMassSig, hMassBkg, hMassSum,
                "plots_final_signal_background/massEE_signal_background.pdf",
                "cMassEE",
                "m_{ee}^{reco} [MeV/c^{2}]");

    SaveSingle(hSigTheta,
               "plots_final_signal_background/significance_thetaEE.pdf",
               "cSignificanceThetaEE",
               "#theta_{ee}^{reco} [deg]");

    SaveSingle(hSigMass,
               "plots_final_signal_background/significance_massEE.pdf",
               "cSignificanceMassEE",
               "m_{ee}^{reco} [MeV/c^{2}]");

    std::cout << "\n=== Final signal/background analysis summary ===\n";
    std::cout << "Signal detector-effects file     : " << signalRecoFile << "\n";
    std::cout << "Background detector-effects file : " << backgroundRecoFile << "\n";
    std::cout << "Output file                      : " << outputFile << "\n";
    std::cout << "Signal reconstructed entries     : " << tSig->GetEntries()
              << "  scale = " << signalScale << "\n";
    std::cout << "Background reconstructed entries : " << tBkg->GetEntries()
              << "  scale = " << backgroundScale << "\n";
    std::cout << "Signal rows processed            : " << nSigFilled << "\n";
    std::cout << "Background rows processed        : " << nBkgFilled << "\n";
    std::cout << "Scaled signal yield              : " << hThetaSig->Integral() << "\n";
    std::cout << "Scaled background yield          : " << hThetaBkg->Integral() << "\n";
    std::cout << "Global contamination B/(S+B)     : "
              << ((hThetaSig->Integral() + hThetaBkg->Integral()) > 0.0
                  ? hThetaBkg->Integral() / (hThetaSig->Integral() + hThetaBkg->Integral())
                  : 0.0)
              << "\n";
    std::cout << "Plots directory                  : plots_final_signal_background\n";
    std::cout << "\nDone.\n";
}

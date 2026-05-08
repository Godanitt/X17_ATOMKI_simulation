// fit_ipc_excess.C
//
// Ajuste final fenomenológico:
//
//   D(x) = fondo IPC-like suave + exceso tipo X17
//
// con:
//
//   f(x) = A exp(-lambda x) + H exp[-0.5 ((x-mu)/sigma)^2]
//
// Lee:
//   analysis_hits_detector_effects.root
//   background_analysis_hits_detector_effects.root
//
// Tree esperado:
//   detected_detector_effects
//
// Genera:
//   ipc_excess_fit.root
//   plots_ipc_excess_fit/

#include <TFile.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TH1.h>

#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <algorithm>
#include <cmath>
#include <string>

namespace ipcfit {

constexpr double me = 0.51099895; // MeV

bool valid(double x) {
    return (x == x) && std::abs(x) < 1e20;
}

bool hasBranch(TTree* t, const char* name) {
    return t && t->GetBranch(name);
}

double massEEFromKinetic(double k1, double k2, double thetaDeg) {
    if (!valid(k1) || !valid(k2) || !valid(thetaDeg)) return -999.0;
    if (k1 < 0.0 || k2 < 0.0) return -999.0;

    const double E1 = k1 + me;
    const double E2 = k2 + me;

    const double p1 = std::sqrt(std::max(0.0, E1 * E1 - me * me));
    const double p2 = std::sqrt(std::max(0.0, E2 * E2 - me * me));

    const double c = std::cos(thetaDeg * TMath::DegToRad());
    const double m2 = 2.0 * me * me + 2.0 * (E1 * E2 - p1 * p2 * c);

    return m2 > 0.0 ? std::sqrt(m2) : 0.0;
}

struct Event {
    int isSignalTruth = -1;
    double thetaEE = -999.0;
    double energyEE = -999.0;
    double massEE = -999.0;
};

std::vector<Event> readReco(TTree* t, int truthLabel) {
    std::vector<Event> rows;

    if (!t) {
        std::cerr << "[ERROR] Null tree\n";
        return rows;
    }

    if (!hasBranch(t, "thetaEE_reco_deg")) {
        std::cerr << "[ERROR] Missing thetaEE_reco_deg in tree "
                  << t->GetName() << "\n";
        return rows;
    }

    const bool hasE1 = hasBranch(t, "energyEe_reco_MeV");
    const bool hasE2 = hasBranch(t, "energyEp_reco_MeV");
    const bool hasEE = hasBranch(t, "energyEE_reco_MeV");
    const bool hasM  = hasBranch(t, "massEE_reco_MeV");

    TTreeFormula fTheta("fTheta", "thetaEE_reco_deg", t);

    std::unique_ptr<TTreeFormula> fE1;
    std::unique_ptr<TTreeFormula> fE2;
    std::unique_ptr<TTreeFormula> fEE;
    std::unique_ptr<TTreeFormula> fM;

    if (hasE1) fE1 = std::make_unique<TTreeFormula>("fE1", "energyEe_reco_MeV", t);
    if (hasE2) fE2 = std::make_unique<TTreeFormula>("fE2", "energyEp_reco_MeV", t);
    if (hasEE) fEE = std::make_unique<TTreeFormula>("fEE", "energyEE_reco_MeV", t);
    if (hasM)  fM  = std::make_unique<TTreeFormula>("fM",  "massEE_reco_MeV", t);

    const Long64_t n = t->GetEntries();
    rows.reserve(n);

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);

        Event ev;
        ev.isSignalTruth = truthLabel;
        ev.thetaEE = fTheta.EvalInstance();

        double e1 = hasE1 ? fE1->EvalInstance() : -999.0;
        double e2 = hasE2 ? fE2->EvalInstance() : -999.0;

        if (hasEE) {
            ev.energyEE = fEE->EvalInstance();
        } else if (hasE1 && hasE2) {
            ev.energyEE = e1 + e2;
        }

        if (hasM) {
            ev.massEE = fM->EvalInstance();
        } else if (hasE1 && hasE2) {
            ev.massEE = massEEFromKinetic(e1, e2, ev.thetaEE);
        }

        if (valid(ev.thetaEE)) rows.push_back(ev);
    }

    return rows;
}

std::vector<Event> sampleRows(const std::vector<Event>& input,
                              Long64_t requested,
                              TRandom3& rng) {
    if (input.empty()) return {};

    if (requested < 0) return input;

    std::vector<Event> out;
    out.reserve(requested);

    for (Long64_t i = 0; i < requested; ++i) {
        int idx = rng.Integer((UInt_t) input.size());
        out.push_back(input[idx]);
    }

    return out;
}

void fillHists(const std::vector<Event>& data,
               TH1D* hTheta,
               TH1D* hMass,
               TH1D* hEnergy) {
    for (const auto& ev : data) {
        if (valid(ev.thetaEE) && ev.thetaEE > -900.0)
            hTheta->Fill(ev.thetaEE);

        if (valid(ev.massEE) && ev.massEE > -900.0)
            hMass->Fill(ev.massEE);

        if (valid(ev.energyEE) && ev.energyEE > -900.0)
            hEnergy->Fill(ev.energyEE);
    }
}

struct FitResult {
    int status = -1;

    double A = 0.0;
    double lambda = 0.0;
    double H = 0.0;
    double mu = 0.0;
    double sigma = 0.0;

    double nSignal = 0.0;
    double nBackground = 0.0;
    double purity = 0.0;
    double contamination = 0.0;
};

double gaussianYield(double H, double sigma, double binWidth) {
    if (sigma <= 0.0 || binWidth <= 0.0) return 0.0;
    return H * sigma * std::sqrt(2.0 * TMath::Pi()) / binWidth;
}

double expYield(double A, double lambda, double xmin, double xmax, double binWidth) {
    if (binWidth <= 0.0) return 0.0;

    if (std::abs(lambda) < 1e-12) {
        return A * (xmax - xmin) / binWidth;
    }

    return A * (std::exp(-lambda * xmin) - std::exp(-lambda * xmax))
           / lambda / binWidth;
}

FitResult fitOne(TH1D* hData,
                 const char* tag,
                 double xmin,
                 double xmax,
                 double mu0,
                 double sigma0,
                 TH1D*& hBkg,
                 TH1D*& hSig,
                 TH1D*& hModel) {
    FitResult r;

    if (!hData || hData->Integral() <= 0.0) {
        std::cerr << "[WARNING] Empty histogram for " << tag << "\n";
        return r;
    }

    const double maxY = std::max(1.0, hData->GetMaximum());

    TF1* model = new TF1(
        Form("f_%s_model", tag),
        "[0]*TMath::Exp(-[1]*x) + [2]*TMath::Exp(-0.5*((x-[3])/[4])^2)",
        xmin,
        xmax
    );

    model->SetParNames("A_IPC", "lambda_IPC", "H_X17", "mu_X17", "sigma_X17");

    model->SetParameters(
        maxY,
        0.01,
        0.25 * maxY,
        mu0,
        sigma0
    );

    model->SetParLimits(0, 0.0, 100.0 * maxY);
    model->SetParLimits(1, 0.0, 10.0);
    model->SetParLimits(2, 0.0, 100.0 * maxY);
    model->SetParLimits(3, xmin, xmax);
    model->SetParLimits(4, 0.05, xmax - xmin);

    r.status = hData->Fit(model, "RQ0");

    r.A      = model->GetParameter(0);
    r.lambda = model->GetParameter(1);
    r.H      = model->GetParameter(2);
    r.mu     = model->GetParameter(3);
    r.sigma  = std::abs(model->GetParameter(4));

    const double binWidth = hData->GetXaxis()->GetBinWidth(1);

    r.nSignal     = gaussianYield(r.H, r.sigma, binWidth);
    r.nBackground = expYield(r.A, r.lambda, xmin, xmax, binWidth);

    const double total = r.nSignal + r.nBackground;

    r.purity        = total > 0.0 ? r.nSignal / total : 0.0;
    r.contamination = total > 0.0 ? r.nBackground / total : 0.0;

    hBkg   = (TH1D*) hData->Clone(Form("h%s_background_smooth", tag));
    hSig   = (TH1D*) hData->Clone(Form("h%s_signal_excess", tag));
    hModel = (TH1D*) hData->Clone(Form("h%s_total_fit", tag));

    hBkg->Reset("ICES");
    hSig->Reset("ICES");
    hModel->Reset("ICES");

    hBkg->SetDirectory(nullptr);
    hSig->SetDirectory(nullptr);
    hModel->SetDirectory(nullptr);

    hBkg->SetLineColor(kRed + 1);
    hBkg->SetLineWidth(3);

    hSig->SetLineColor(kBlue + 1);
    hSig->SetLineWidth(3);

    hModel->SetLineColor(kBlack);
    hModel->SetLineWidth(4);

    for (int i = 1; i <= hData->GetNbinsX(); ++i) {
        const double x = hData->GetBinCenter(i);

        const double b = r.A * std::exp(-r.lambda * x);
        const double s = r.H * std::exp(-0.5 * std::pow((x - r.mu) / r.sigma, 2));
        const double m = b + s;

        const double d = hData->GetBinContent(i);

        hBkg->SetBinContent(i, b);
        hSig->SetBinContent(i, d - b);
        hSig->SetBinError(i, hData->GetBinError(i));
        hModel->SetBinContent(i, m);
    }

    model->SetLineColor(kBlack);
    model->SetLineWidth(4);

    return r;
}

void saveFitPlot(TH1D* hData,
                 TH1D* hBkg,
                 TH1D* hSig,
                 TH1D* hModel,
                 const std::string& pdf,
                 const std::string& xtitle) {
    if (!hData || !hBkg || !hSig || !hModel) return;

    TCanvas c("c", "c", 950, 720);
    c.SetTicks(1, 1);

    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(0.9);
    hData->SetLineColor(kBlack);
    hData->SetMarkerColor(kBlack);

    double ymax = std::max({
        hData->GetMaximum(),
        hBkg->GetMaximum(),
        hSig->GetMaximum(),
        hModel->GetMaximum()
    });

    hData->SetMaximum(ymax > 0.0 ? 1.35 * ymax : 1.0);
    hData->GetXaxis()->SetTitle(xtitle.c_str());
    hData->GetYaxis()->SetTitle("events");

    hData->Draw("E");
    hBkg->Draw("HIST SAME");
    hSig->Draw("HIST SAME");
    hModel->Draw("HIST SAME");
    hData->Draw("E SAME");

    TLegend leg(0.55, 0.66, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(hData, "D = signal + background", "pe");
    leg.AddEntry(hBkg, "smooth IPC-like background", "l");
    leg.AddEntry(hSig, "extracted excess D - background", "l");
    leg.AddEntry(hModel, "background + Gaussian excess", "l");
    leg.Draw();

    c.SaveAs(pdf.c_str());
}

void saveSignalPlot(TH1D* hSig,
                    const std::string& pdf,
                    const std::string& xtitle) {
    if (!hSig) return;

    TCanvas c("cSig", "cSig", 900, 700);
    c.SetTicks(1, 1);

    hSig->SetLineColor(kBlue + 1);
    hSig->SetLineWidth(3);
    hSig->GetXaxis()->SetTitle(xtitle.c_str());
    hSig->GetYaxis()->SetTitle("excess events");
    hSig->Draw("HIST E");

    c.SaveAs(pdf.c_str());
}

} // namespace ipcfit


void fit_ipc_excess(
    const char* signalRecoFile = "analysis_hits_detector_effects.root",
    const char* backgroundRecoFile = "background_analysis_hits_detector_effects.root",
    const char* outputRootFile = "ipc_excess_fit.root",
    const char* outputPlotDir = "plots_ipc_excess_fit",
    Long64_t nSignalToMix = -1,
    Long64_t nBackgroundToMix = -1,
    unsigned int seed = 12345,
    double thetaMu0 = 140.0,
    double thetaSigma0 = 12.0,
    double massMu0 = 17.0,
    double massSigma0 = 1.0,
    double energyMu0 = 17.5,
    double energySigma0 = 0.8
) {
    using namespace ipcfit;

    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    gSystem->mkdir(outputPlotDir, kTRUE);

    std::unique_ptr<TFile> fSig(TFile::Open(signalRecoFile, "READ"));
    std::unique_ptr<TFile> fBkg(TFile::Open(backgroundRecoFile, "READ"));

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
        std::cerr << "[ERROR] Missing tree detected_detector_effects in "
                  << signalRecoFile << "\n";
        return;
    }

    if (!tBkg) {
        std::cerr << "[ERROR] Missing tree detected_detector_effects in "
                  << backgroundRecoFile << "\n";
        return;
    }

    auto sigAll = readReco(tSig, 1);
    auto bkgAll = readReco(tBkg, 0);

    if (sigAll.empty() || bkgAll.empty()) {
        std::cerr << "[ERROR] Need non-empty signal and background samples.\n";
        std::cerr << "        Generate more events or loosen selections.\n";
        return;
    }

    TRandom3 rng(seed);

    auto sig = sampleRows(sigAll, nSignalToMix, rng);
    auto bkg = sampleRows(bkgAll, nBackgroundToMix, rng);

    std::vector<Event> data;
    data.reserve(sig.size() + bkg.size());
    data.insert(data.end(), sig.begin(), sig.end());
    data.insert(data.end(), bkg.begin(), bkg.end());

    std::mt19937 engine(seed);
    std::shuffle(data.begin(), data.end(), engine);

    std::cout << "\n=== IPC smooth background + Gaussian excess fit ===\n";
    std::cout << "Signal entries available      : " << sigAll.size() << "\n";
    std::cout << "Background entries available  : " << bkgAll.size() << "\n";
    std::cout << "Injected signal events        : " << sig.size() << "\n";
    std::cout << "Injected background events    : " << bkg.size() << "\n";
    std::cout << "Pseudo-data D events          : " << data.size() << "\n";

    TH1D hTheta("hThetaEE_data",
                "Pseudo-data opening angle;#theta_{ee}^{reco} [deg];events",
                90, 0, 180);

    TH1D hMass("hMassEE_data",
               "Pseudo-data invariant mass;m_{ee}^{reco} [MeV/c^{2}];events",
               80, 0, 25);

    TH1D hEnergy("hEnergyEE_data",
                 "Pseudo-data pair energy;E_{ee}^{reco} [MeV];events",
                 80, 0, 25);

    hTheta.Sumw2();
    hMass.Sumw2();
    hEnergy.Sumw2();

    fillHists(data, &hTheta, &hMass, &hEnergy);

    TH1D* hThetaBkg = nullptr;
    TH1D* hThetaSig = nullptr;
    TH1D* hThetaModel = nullptr;

    TH1D* hMassBkg = nullptr;
    TH1D* hMassSig = nullptr;
    TH1D* hMassModel = nullptr;

    TH1D* hEnergyBkg = nullptr;
    TH1D* hEnergySig = nullptr;
    TH1D* hEnergyModel = nullptr;

    FitResult thetaFit = fitOne(
        &hTheta,
        "ThetaEE",
        0.0,
        180.0,
        thetaMu0,
        thetaSigma0,
        hThetaBkg,
        hThetaSig,
        hThetaModel
    );

    FitResult massFit = fitOne(
        &hMass,
        "MassEE",
        0.0,
        25.0,
        massMu0,
        massSigma0,
        hMassBkg,
        hMassSig,
        hMassModel
    );

    FitResult energyFit = fitOne(
        &hEnergy,
        "EnergyEE",
        0.0,
        25.0,
        energyMu0,
        energySigma0,
        hEnergyBkg,
        hEnergySig,
        hEnergyModel
    );

    std::cout << "\n=== Fit results ===\n";

    std::cout << "ThetaEE: Nsig = " << thetaFit.nSignal
              << ", Nbkg = " << thetaFit.nBackground
              << ", purity = " << thetaFit.purity
              << ", contamination = " << thetaFit.contamination
              << ", mu = " << thetaFit.mu
              << ", sigma = " << thetaFit.sigma
              << ", status = " << thetaFit.status << "\n";

    std::cout << "MassEE : Nsig = " << massFit.nSignal
              << ", Nbkg = " << massFit.nBackground
              << ", purity = " << massFit.purity
              << ", contamination = " << massFit.contamination
              << ", mu = " << massFit.mu
              << ", sigma = " << massFit.sigma
              << ", status = " << massFit.status << "\n";

    std::cout << "Energy : Nsig = " << energyFit.nSignal
              << ", Nbkg = " << energyFit.nBackground
              << ", purity = " << energyFit.purity
              << ", contamination = " << energyFit.contamination
              << ", mu = " << energyFit.mu
              << ", sigma = " << energyFit.sigma
              << ", status = " << energyFit.status << "\n";

    std::string outDir(outputPlotDir);

    saveFitPlot(&hTheta, hThetaBkg, hThetaSig, hThetaModel,
                outDir + "/fit_thetaEE.pdf",
                "#theta_{ee}^{reco} [deg]");

    saveFitPlot(&hMass, hMassBkg, hMassSig, hMassModel,
                outDir + "/fit_massEE.pdf",
                "m_{ee}^{reco} [MeV/c^{2}]");

    saveFitPlot(&hEnergy, hEnergyBkg, hEnergySig, hEnergyModel,
                outDir + "/fit_energyEE.pdf",
                "E_{ee}^{reco} [MeV]");

    saveSignalPlot(hThetaSig,
                   outDir + "/signal_excess_thetaEE.pdf",
                   "#theta_{ee}^{reco} [deg]");

    saveSignalPlot(hMassSig,
                   outDir + "/signal_excess_massEE.pdf",
                   "m_{ee}^{reco} [MeV/c^{2}]");

    saveSignalPlot(hEnergySig,
                   outDir + "/signal_excess_energyEE.pdf",
                   "E_{ee}^{reco} [MeV]");

    std::unique_ptr<TFile> fout(TFile::Open(outputRootFile, "RECREATE"));

    if (!fout || fout->IsZombie()) {
        std::cerr << "[ERROR] Cannot create output ROOT file: "
                  << outputRootFile << "\n";
        return;
    }

    fout->cd();

    TTree tSummary("fit_summary", "IPC-like smooth background + Gaussian excess summary");

    double nTruthSignal = (double) sig.size();
    double nTruthBackground = (double) bkg.size();
    double truthContamination =
        (nTruthSignal + nTruthBackground) > 0.0
        ? nTruthBackground / (nTruthSignal + nTruthBackground)
        : 0.0;

    int thetaStatus = thetaFit.status;
    int massStatus = massFit.status;
    int energyStatus = energyFit.status;

    tSummary.Branch("nTruthSignal", &nTruthSignal, "nTruthSignal/D");
    tSummary.Branch("nTruthBackground", &nTruthBackground, "nTruthBackground/D");
    tSummary.Branch("truthContamination", &truthContamination, "truthContamination/D");

    tSummary.Branch("thetaFitStatus", &thetaStatus, "thetaFitStatus/I");
    tSummary.Branch("thetaSignalYield", &thetaFit.nSignal, "thetaSignalYield/D");
    tSummary.Branch("thetaBackgroundYield", &thetaFit.nBackground, "thetaBackgroundYield/D");
    tSummary.Branch("thetaPurity", &thetaFit.purity, "thetaPurity/D");
    tSummary.Branch("thetaContamination", &thetaFit.contamination, "thetaContamination/D");
    tSummary.Branch("thetaSignalMean", &thetaFit.mu, "thetaSignalMean/D");
    tSummary.Branch("thetaSignalSigma", &thetaFit.sigma, "thetaSignalSigma/D");

    tSummary.Branch("massFitStatus", &massStatus, "massFitStatus/I");
    tSummary.Branch("massSignalYield", &massFit.nSignal, "massSignalYield/D");
    tSummary.Branch("massBackgroundYield", &massFit.nBackground, "massBackgroundYield/D");
    tSummary.Branch("massPurity", &massFit.purity, "massPurity/D");
    tSummary.Branch("massContamination", &massFit.contamination, "massContamination/D");
    tSummary.Branch("massSignalMean", &massFit.mu, "massSignalMean/D");
    tSummary.Branch("massSignalSigma", &massFit.sigma, "massSignalSigma/D");

    tSummary.Branch("energyFitStatus", &energyStatus, "energyFitStatus/I");
    tSummary.Branch("energySignalYield", &energyFit.nSignal, "energySignalYield/D");
    tSummary.Branch("energyBackgroundYield", &energyFit.nBackground, "energyBackgroundYield/D");
    tSummary.Branch("energyPurity", &energyFit.purity, "energyPurity/D");
    tSummary.Branch("energyContamination", &energyFit.contamination, "energyContamination/D");
    tSummary.Branch("energySignalMean", &energyFit.mu, "energySignalMean/D");
    tSummary.Branch("energySignalSigma", &energyFit.sigma, "energySignalSigma/D");

    tSummary.Fill();

    hTheta.Write();
    hMass.Write();
    hEnergy.Write();

    if (hThetaBkg) hThetaBkg->Write();
    if (hThetaSig) hThetaSig->Write();
    if (hThetaModel) hThetaModel->Write();

    if (hMassBkg) hMassBkg->Write();
    if (hMassSig) hMassSig->Write();
    if (hMassModel) hMassModel->Write();

    if (hEnergyBkg) hEnergyBkg->Write();
    if (hEnergySig) hEnergySig->Write();
    if (hEnergyModel) hEnergyModel->Write();

    tSummary.Write();

    fout->Close();

    std::cout << "\n=== Outputs ===\n";
    std::cout << "Output ROOT file : " << outputRootFile << "\n";
    std::cout << "Plot directory   : " << outputPlotDir << "\n";
    std::cout << "Summary tree     : fit_summary\n";
    std::cout << "\nDone.\n";
}

// extract_root_summaries.C
//
// ROOT-enabled helper for the Typst report.  It extracts the official summary
// TTrees produced by analysis/analyze_hits.C and
// analysis/mix_and_fit_signal_background_fixed_v2.C and writes compact CSV
// tables under memoria/tables/.
//
// Usage from the project root:
//   root -l -q 'memoria/scripts/extract_root_summaries.C()'
//
// Optional:
//   root -l -q 'memoria/scripts/extract_root_summaries.C("results", "memoria/tables")'

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace memoria_extract {

struct EfficiencyRow {
    std::string geometry;
    std::string sample;
    Long64_t nGenerated = 0;
    Long64_t nDetectedElectron = 0;
    Long64_t nDetectedPositron = 0;
    Long64_t nDetectedCoincidence = 0;
    double effCoincidence = 0.0;
    double effCoincidenceErr = 0.0;
};

bool read_efficiency(const std::string& filename,
                     const std::string& geometryFallback,
                     const std::string& sample,
                     EfficiencyRow& row)
{
    TFile f(filename.c_str(), "READ");
    if (f.IsZombie()) {
        std::cerr << "[warning] cannot open " << filename << "\n";
        return false;
    }

    auto* t = dynamic_cast<TTree*>(f.Get("summary"));
    if (!t) {
        std::cerr << "[warning] no summary tree in " << filename << "\n";
        return false;
    }

    std::string* geometry = nullptr;
    int selectedVolumeID = 0;
    Long64_t nHitRows = 0;
    double effElectron = 0.0;
    double effPositron = 0.0;

    t->SetBranchAddress("geometry", &geometry);
    t->SetBranchAddress("selectedVolumeID", &selectedVolumeID);
    t->SetBranchAddress("nGenerated", &row.nGenerated);
    t->SetBranchAddress("nHitRows", &nHitRows);
    t->SetBranchAddress("nDetectedElectron", &row.nDetectedElectron);
    t->SetBranchAddress("nDetectedPositron", &row.nDetectedPositron);
    t->SetBranchAddress("nDetectedCoincidence", &row.nDetectedCoincidence);
    t->SetBranchAddress("effElectron", &effElectron);
    t->SetBranchAddress("effPositron", &effPositron);
    t->SetBranchAddress("effCoincidence", &row.effCoincidence);

    t->GetEntry(0);

    row.geometry = (geometry && !geometry->empty()) ? *geometry : geometryFallback;
    row.sample = sample;
    row.effCoincidenceErr = row.nGenerated > 0
        ? std::sqrt(row.effCoincidence * (1.0 - row.effCoincidence) / double(row.nGenerated))
        : 0.0;
    return true;
}

void write_efficiency_typst(const std::vector<EfficiencyRow>& rows,
                            const std::string& outTyp)
{
    std::ofstream out(outTyp);
    out << "#table(\n";
    out << "  columns: (1.15fr, 1.20fr, 1.20fr),\n";
    out << "  inset: 5pt,\n";
    out << "  align: (left, right, right),\n";
    out << "  [*Geometry*], [*Signal*], [*Background*],\n";

    const std::vector<std::string> geoms = {"current", "2pi", "4pi", "padplane"};
    for (const auto& g : geoms) {
        const EfficiencyRow* s = nullptr;
        const EfficiencyRow* b = nullptr;
        for (const auto& r : rows) {
            if (r.geometry == g && r.sample == "signal") s = &r;
            if (r.geometry == g && r.sample == "background") b = &r;
        }
        if (!s || !b) continue;
        out << "  [" << g << "], ["
            << std::fixed << std::setprecision(2) << 100.0 * s->effCoincidence
            << " ± " << 100.0 * s->effCoincidenceErr << "%], ["
            << 100.0 * b->effCoincidence << " ± " << 100.0 * b->effCoincidenceErr
            << "%],\n";
    }
    out << ")\n";
}

void write_fit_csv_and_typst(const std::string& resultsDir, const std::string& outDir)
{
    const std::vector<std::string> geoms = {"current", "2pi", "4pi", "padplane"};
    std::ofstream csv(outDir + "/root_fit_summary.csv");
    csv << "geometry,fitStatus,nSignalFit,nBackgroundFit,purityGlobal,contaminationGlobal,"
        << "signalRegionFit,backgroundRegionFit,purityRegion,contaminationRegion,"
        << "truthSignalGlobal,truthBackgroundGlobal,truthSignalRegion,truthBackgroundRegion\n";

    std::ofstream typ(outDir + "/root_fit_summary.typ");
    typ << "#table(\n";
    typ << "  columns: (1.05fr, 1.10fr, 1.10fr, 1.05fr),\n";
    typ << "  inset: 4.5pt,\n";
    typ << "  align: (left, right, right, right),\n";
    typ << "  [*Geometry*], [*$N_X^fit$*], [*$N_b^fit$*], [*purity*],\n";

    for (const auto& g : geoms) {
        const std::string filename = resultsDir + "/" + g + "/signal_background_fit.root";
        TFile f(filename.c_str(), "READ");
        if (f.IsZombie()) {
            std::cerr << "[warning] cannot open " << filename << "\n";
            continue;
        }
        auto* t = dynamic_cast<TTree*>(f.Get("fit_summary"));
        if (!t) {
            std::cerr << "[warning] no fit_summary tree in " << filename << "\n";
            continue;
        }

        int fitStatus = -1;
        double nSignalFit = 0.0, nBackgroundFit = 0.0;
        double purityGlobal = 0.0, contaminationGlobal = 0.0;
        double signalRegionFit = 0.0, backgroundRegionFit = 0.0;
        double purityRegion = 0.0, contaminationRegion = 0.0;
        double truthSignalGlobal = 0.0, truthBackgroundGlobal = 0.0;
        double truthSignalRegion = 0.0, truthBackgroundRegion = 0.0;

        t->SetBranchAddress("fitStatus", &fitStatus);
        t->SetBranchAddress("nSignalFit", &nSignalFit);
        t->SetBranchAddress("nBackgroundFit", &nBackgroundFit);
        t->SetBranchAddress("purityGlobal", &purityGlobal);
        t->SetBranchAddress("contaminationGlobal", &contaminationGlobal);
        t->SetBranchAddress("signalRegionFit", &signalRegionFit);
        t->SetBranchAddress("backgroundRegionFit", &backgroundRegionFit);
        t->SetBranchAddress("purityRegion", &purityRegion);
        t->SetBranchAddress("contaminationRegion", &contaminationRegion);
        t->SetBranchAddress("truthSignalGlobal", &truthSignalGlobal);
        t->SetBranchAddress("truthBackgroundGlobal", &truthBackgroundGlobal);
        t->SetBranchAddress("truthSignalRegion", &truthSignalRegion);
        t->SetBranchAddress("truthBackgroundRegion", &truthBackgroundRegion);
        t->GetEntry(0);

        csv << g << "," << fitStatus << "," << nSignalFit << "," << nBackgroundFit << ","
            << purityGlobal << "," << contaminationGlobal << ","
            << signalRegionFit << "," << backgroundRegionFit << ","
            << purityRegion << "," << contaminationRegion << ","
            << truthSignalGlobal << "," << truthBackgroundGlobal << ","
            << truthSignalRegion << "," << truthBackgroundRegion << "\n";

        typ << "  [" << g << "], [" << std::fixed << std::setprecision(0) << nSignalFit
            << "], [" << nBackgroundFit << "], [" << std::setprecision(1) << 100.0 * purityGlobal
            << "%],\n";
    }
    typ << ")\n";
}

} // namespace memoria_extract

void extract_root_summaries(const char* resultsDir = "results",
                            const char* outDir = "memoria/tables")
{
    using namespace memoria_extract;
    gSystem->mkdir(outDir, true);

    const std::vector<std::string> geoms = {"current", "2pi", "4pi", "padplane"};
    std::vector<EfficiencyRow> rows;

    for (const auto& g : geoms) {
        EfficiencyRow sig, bkg;
        if (read_efficiency(std::string(resultsDir) + "/" + g + "/analysis_hits.root", g, "signal", sig)) {
            rows.push_back(sig);
        }
        if (read_efficiency(std::string(resultsDir) + "/" + g + "/background_analysis_hits.root", g, "background", bkg)) {
            rows.push_back(bkg);
        }
    }

    std::ofstream csv(std::string(outDir) + "/root_geometrical_efficiencies.csv");
    csv << "geometry,sample,nGenerated,nDetectedElectron,nDetectedPositron,"
        << "nDetectedCoincidence,effCoincidence,effCoincidenceErr\n";
    for (const auto& r : rows) {
        csv << r.geometry << "," << r.sample << "," << r.nGenerated << ","
            << r.nDetectedElectron << "," << r.nDetectedPositron << ","
            << r.nDetectedCoincidence << "," << std::setprecision(10)
            << r.effCoincidence << "," << r.effCoincidenceErr << "\n";
    }
    write_efficiency_typst(rows, std::string(outDir) + "/root_geometrical_efficiencies.typ");
    write_fit_csv_and_typst(resultsDir, outDir);

    std::cout << "Wrote ROOT-derived CSV/Typst tables under " << outDir << "\n";
}

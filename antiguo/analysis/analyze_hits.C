#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TError.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>

namespace
{
struct GenRow
{
    int eventID = -1;
    int dataRow = -1;
    double thetaEE_txt_deg = -999.0;
    double thetaEe_txt_deg = -999.0;
    double energyEe_txt_MeV = -999.0;
    double thetaEp_txt_deg = -999.0;
    double energyEp_txt_MeV = -999.0;
};

struct HitRow
{
    int eventID = -1;
    int dataRow = -1;
    int trackID = -1;
    int pdg = 0;
    int charge = 0;
    int detID = -1;
    int volumeID = -1;
    double x_hit_mm = -999.0;
    double y_hit_mm = -999.0;
    double z_hit_mm = -999.0;
    double theta_hit_deg = -999.0;
    double phi_hit_deg = -999.0;
    double kinetic_hit_MeV = -999.0;
    double dir_hit_x = 0.0;
    double dir_hit_y = 0.0;
    double dir_hit_z = 0.0;
    double theta_gen_deg = -999.0;
    double phi_gen_deg = -999.0;
    double kinetic_gen_MeV = -999.0;
};

struct PairHits
{
    bool hasEe = false;
    bool hasEp = false;
    HitRow Ee;
    HitRow Ep;
};

double opening_angle_deg(const HitRow& a, const HitRow& b)
{
    const double ma = std::sqrt(a.dir_hit_x*a.dir_hit_x + a.dir_hit_y*a.dir_hit_y + a.dir_hit_z*a.dir_hit_z);
    const double mb = std::sqrt(b.dir_hit_x*b.dir_hit_x + b.dir_hit_y*b.dir_hit_y + b.dir_hit_z*b.dir_hit_z);
    if (ma <= 0.0 || mb <= 0.0) return -999.0;

    double c = (a.dir_hit_x*b.dir_hit_x + a.dir_hit_y*b.dir_hit_y + a.dir_hit_z*b.dir_hit_z) / (ma * mb);
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c) * 180.0 / TMath::Pi();
}

const char* volume_name(int volumeID)
{
    if (volumeID == 0) return "SiliconStripLV";
    if (volumeID == 1) return "ScintillatorLV";
    return "all detector volumes";
}
}

void analyze_hits(const char* inputFile = "sampled.root",
                  const char* outputFile = "analysis_hits.root",
                  int selectedVolumeID = 0)
{
    // selectedVolumeID = 0 -> use silicon hits for the reconstructed pair.
    // selectedVolumeID = 1 -> use scintillator hits.
    // selectedVolumeID < 0 -> allow both, taking the first hit per particle/event.

    TFile fin(inputFile, "READ");
    if (fin.IsZombie()) {
        std::cerr << "ERROR: cannot open input ROOT file: " << inputFile << std::endl;
        return;
    }

    auto* generated = dynamic_cast<TTree*>(fin.Get("generated"));
    auto* hits = dynamic_cast<TTree*>(fin.Get("hits"));
    if (!generated || !hits) {
        std::cerr << "ERROR: input file must contain trees named generated and hits." << std::endl;
        return;
    }

    std::map<int, GenRow> genByEvent;

    GenRow g;
    generated->SetBranchAddress("eventID", &g.eventID);
    generated->SetBranchAddress("dataRow", &g.dataRow);
    generated->SetBranchAddress("thetaEE_txt_deg", &g.thetaEE_txt_deg);
    generated->SetBranchAddress("thetaEe_txt_deg", &g.thetaEe_txt_deg);
    generated->SetBranchAddress("energyEe_txt_MeV", &g.energyEe_txt_MeV);
    generated->SetBranchAddress("thetaEp_txt_deg", &g.thetaEp_txt_deg);
    generated->SetBranchAddress("energyEp_txt_MeV", &g.energyEp_txt_MeV);

    for (Long64_t i = 0; i < generated->GetEntries(); ++i) {
        generated->GetEntry(i);
        genByEvent[g.eventID] = g;
    }

    std::map<int, PairHits> pairByEvent;
    std::map<int, int> hitMultiplicityByEvent;
    std::set<int> eventsWithEe;
    std::set<int> eventsWithEp;

    HitRow h;
    hits->SetBranchAddress("eventID", &h.eventID);
    hits->SetBranchAddress("dataRow", &h.dataRow);
    hits->SetBranchAddress("trackID", &h.trackID);
    hits->SetBranchAddress("pdg", &h.pdg);
    hits->SetBranchAddress("charge", &h.charge);
    hits->SetBranchAddress("detID", &h.detID);
    hits->SetBranchAddress("volumeID", &h.volumeID);
    hits->SetBranchAddress("x_hit_mm", &h.x_hit_mm);
    hits->SetBranchAddress("y_hit_mm", &h.y_hit_mm);
    hits->SetBranchAddress("z_hit_mm", &h.z_hit_mm);
    hits->SetBranchAddress("theta_hit_deg", &h.theta_hit_deg);
    hits->SetBranchAddress("phi_hit_deg", &h.phi_hit_deg);
    hits->SetBranchAddress("kinetic_hit_MeV", &h.kinetic_hit_MeV);
    hits->SetBranchAddress("dir_hit_x", &h.dir_hit_x);
    hits->SetBranchAddress("dir_hit_y", &h.dir_hit_y);
    hits->SetBranchAddress("dir_hit_z", &h.dir_hit_z);
    hits->SetBranchAddress("theta_gen_deg", &h.theta_gen_deg);
    hits->SetBranchAddress("phi_gen_deg", &h.phi_gen_deg);
    hits->SetBranchAddress("kinetic_gen_MeV", &h.kinetic_gen_MeV);

    for (Long64_t i = 0; i < hits->GetEntries(); ++i) {
        hits->GetEntry(i);
        ++hitMultiplicityByEvent[h.eventID];

        if (selectedVolumeID >= 0 && h.volumeID != selectedVolumeID) continue;

        auto& pair = pairByEvent[h.eventID];
        if (h.pdg == 11) {
            eventsWithEe.insert(h.eventID);
            if (!pair.hasEe) { pair.Ee = h; pair.hasEe = true; }
        } else if (h.pdg == -11) {
            eventsWithEp.insert(h.eventID);
            if (!pair.hasEp) { pair.Ep = h; pair.hasEp = true; }
        }
    }

    TFile fout(outputFile, "RECREATE");

    TH1D hThetaEEtxt("hThetaEEtxt", "Generated opening angle;#theta_{ee}^{txt} [deg];Events", 180, 0, 180);
    TH1D hThetaEEhit("hThetaEEhit", "Reconstructed opening angle from hit directions;#theta_{ee}^{hit} [deg];Coincidences", 180, 0, 180);
    TH1D hThetaEeHit("hThetaEeHit", "Electron hit polar angle;#theta_{e-}^{hit} [deg];Hits", 180, 0, 180);
    TH1D hThetaEpHit("hThetaEpHit", "Positron hit polar angle;#theta_{e+}^{hit} [deg];Hits", 180, 0, 180);
    TH1D hEnergyEeHit("hEnergyEeHit", "Electron hit kinetic energy;E_{e-}^{hit} [MeV];Hits", 200, 0, 20);
    TH1D hEnergyEpHit("hEnergyEpHit", "Positron hit kinetic energy;E_{e+}^{hit} [MeV];Hits", 200, 0, 20);
    TH1D hDetEe("hDetEe", "Electron detector ID;detID;Hits", 6, -0.5, 5.5);
    TH1D hDetEp("hDetEp", "Positron detector ID;detID;Hits", 6, -0.5, 5.5);
    TH1D hHitMultiplicity("hHitMultiplicity", "Hit rows per generated event;hit rows/event;Events", 12, -0.5, 11.5);
    TH2D hThetaEEtxtVsHit("hThetaEEtxtVsHit", "Opening angle: generated vs hit;#theta_{ee}^{txt} [deg];#theta_{ee}^{hit} [deg]", 180, 0, 180, 180, 0, 180);

    int eventID = -1;
    int dataRow = -1;
    int detIDEe = -1;
    int detIDEp = -1;
    double thetaEE_txt_deg = -999.0;
    double thetaEE_hit_deg = -999.0;
    double thetaEe_txt_deg = -999.0;
    double thetaEp_txt_deg = -999.0;
    double energyEe_txt_MeV = -999.0;
    double energyEp_txt_MeV = -999.0;
    double thetaEe_hit_deg = -999.0;
    double thetaEp_hit_deg = -999.0;
    double phiEe_hit_deg = -999.0;
    double phiEp_hit_deg = -999.0;
    double energyEe_hit_MeV = -999.0;
    double energyEp_hit_MeV = -999.0;
    double xEe_hit_mm = -999.0, yEe_hit_mm = -999.0, zEe_hit_mm = -999.0;
    double xEp_hit_mm = -999.0, yEp_hit_mm = -999.0, zEp_hit_mm = -999.0;

    TTree detected("detected", "one row per event with both e- and e+ detected in selected volume");
    detected.Branch("eventID", &eventID);
    detected.Branch("dataRow", &dataRow);
    detected.Branch("selectedVolumeID", &selectedVolumeID);
    detected.Branch("detIDEe", &detIDEe);
    detected.Branch("detIDEp", &detIDEp);
    detected.Branch("thetaEE_txt_deg", &thetaEE_txt_deg);
    detected.Branch("thetaEE_hit_deg", &thetaEE_hit_deg);
    detected.Branch("thetaEe_txt_deg", &thetaEe_txt_deg);
    detected.Branch("thetaEp_txt_deg", &thetaEp_txt_deg);
    detected.Branch("energyEe_txt_MeV", &energyEe_txt_MeV);
    detected.Branch("energyEp_txt_MeV", &energyEp_txt_MeV);
    detected.Branch("thetaEe_hit_deg", &thetaEe_hit_deg);
    detected.Branch("thetaEp_hit_deg", &thetaEp_hit_deg);
    detected.Branch("phiEe_hit_deg", &phiEe_hit_deg);
    detected.Branch("phiEp_hit_deg", &phiEp_hit_deg);
    detected.Branch("energyEe_hit_MeV", &energyEe_hit_MeV);
    detected.Branch("energyEp_hit_MeV", &energyEp_hit_MeV);
    detected.Branch("xEe_hit_mm", &xEe_hit_mm);
    detected.Branch("yEe_hit_mm", &yEe_hit_mm);
    detected.Branch("zEe_hit_mm", &zEe_hit_mm);
    detected.Branch("xEp_hit_mm", &xEp_hit_mm);
    detected.Branch("yEp_hit_mm", &yEp_hit_mm);
    detected.Branch("zEp_hit_mm", &zEp_hit_mm);

    for (const auto& [eid, gen] : genByEvent) {
        hThetaEEtxt.Fill(gen.thetaEE_txt_deg);
        const int mult = hitMultiplicityByEvent.count(eid) ? hitMultiplicityByEvent[eid] : 0;
        hHitMultiplicity.Fill(mult);
    }

    Long64_t nCoincidences = 0;
    for (const auto& [eid, pair] : pairByEvent) {
        if (!pair.hasEe || !pair.hasEp) continue;

        const auto genIt = genByEvent.find(eid);
        if (genIt == genByEvent.end()) continue;
        const auto& gen = genIt->second;

        eventID = eid;
        dataRow = gen.dataRow;
        detIDEe = pair.Ee.detID;
        detIDEp = pair.Ep.detID;

        thetaEE_txt_deg = gen.thetaEE_txt_deg;
        thetaEe_txt_deg = gen.thetaEe_txt_deg;
        thetaEp_txt_deg = gen.thetaEp_txt_deg;
        energyEe_txt_MeV = gen.energyEe_txt_MeV;
        energyEp_txt_MeV = gen.energyEp_txt_MeV;

        thetaEe_hit_deg = pair.Ee.theta_hit_deg;
        thetaEp_hit_deg = pair.Ep.theta_hit_deg;
        phiEe_hit_deg = pair.Ee.phi_hit_deg;
        phiEp_hit_deg = pair.Ep.phi_hit_deg;
        energyEe_hit_MeV = pair.Ee.kinetic_hit_MeV;
        energyEp_hit_MeV = pair.Ep.kinetic_hit_MeV;
        xEe_hit_mm = pair.Ee.x_hit_mm;
        yEe_hit_mm = pair.Ee.y_hit_mm;
        zEe_hit_mm = pair.Ee.z_hit_mm;
        xEp_hit_mm = pair.Ep.x_hit_mm;
        yEp_hit_mm = pair.Ep.y_hit_mm;
        zEp_hit_mm = pair.Ep.z_hit_mm;
        thetaEE_hit_deg = opening_angle_deg(pair.Ee, pair.Ep);

        hThetaEEhit.Fill(thetaEE_hit_deg);
        hThetaEeHit.Fill(thetaEe_hit_deg);
        hThetaEpHit.Fill(thetaEp_hit_deg);
        hEnergyEeHit.Fill(energyEe_hit_MeV);
        hEnergyEpHit.Fill(energyEp_hit_MeV);
        hDetEe.Fill(detIDEe);
        hDetEp.Fill(detIDEp);
        hThetaEEtxtVsHit.Fill(thetaEE_txt_deg, thetaEE_hit_deg);

        detected.Fill();
        ++nCoincidences;
    }

    hThetaEEtxt.Write();
    hThetaEEhit.Write();
    hThetaEeHit.Write();
    hThetaEpHit.Write();
    hEnergyEeHit.Write();
    hEnergyEpHit.Write();
    hDetEe.Write();
    hDetEp.Write();
    hHitMultiplicity.Write();
    hThetaEEtxtVsHit.Write();
    detected.Write();

    fout.Close();

    const auto nGen = generated->GetEntries();
    const auto nHitRows = hits->GetEntries();
    const double effEe = nGen > 0 ? double(eventsWithEe.size()) / double(nGen) : 0.0;
    const double effEp = nGen > 0 ? double(eventsWithEp.size()) / double(nGen) : 0.0;
    const double effCoin = nGen > 0 ? double(nCoincidences) / double(nGen) : 0.0;

    std::cout << "\n=== X17 hit analysis summary ===\n";
    std::cout << "Input ROOT file              : " << inputFile << "\n";
    std::cout << "Output analysis file         : " << outputFile << "\n";
    std::cout << "Selected hit volumeID        : " << selectedVolumeID << " (" << volume_name(selectedVolumeID) << ")\n";
    std::cout << "Generated events             : " << nGen << "\n";
    std::cout << "Raw hit rows in hits tree    : " << nHitRows << "\n";
    std::cout << "Events with detected e-      : " << eventsWithEe.size() << "  efficiency = " << effEe << "\n";
    std::cout << "Events with detected e+      : " << eventsWithEp.size() << "  efficiency = " << effEp << "\n";
    std::cout << "Events with e- AND e+        : " << nCoincidences << "  efficiency = " << effCoin << "\n";
    std::cout << "\nRemember: hits entries are particle/volume crossings, not detected events.\n";
    std::cout << "A complete event requires one pdg=11 hit and one pdg=-11 hit with the same eventID.\n";
    std::cout << "================================\n";
}

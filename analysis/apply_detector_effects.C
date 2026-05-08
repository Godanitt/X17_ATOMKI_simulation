#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <iostream>

namespace
{
double clamp(double x, double lo, double hi)
{
    return std::max(lo, std::min(hi, x));
}

double wrap_phi_deg(double phi)
{
    while (phi <= -180.0) phi += 360.0;
    while (phi >   180.0) phi -= 360.0;
    return phi;
}

void direction_from_angles(double thetaDeg, double phiDeg, double& x, double& y, double& z)
{
    const double th = thetaDeg * TMath::Pi() / 180.0;
    const double ph = phiDeg   * TMath::Pi() / 180.0;
    x = std::sin(th) * std::cos(ph);
    y = std::sin(th) * std::sin(ph);
    z = std::cos(th);
}

double opening_angle_deg_from_angles(double theta1Deg, double phi1Deg,
                                     double theta2Deg, double phi2Deg)
{
    double x1, y1, z1;
    double x2, y2, z2;
    direction_from_angles(theta1Deg, phi1Deg, x1, y1, z1);
    direction_from_angles(theta2Deg, phi2Deg, x2, y2, z2);

    double c = x1*x2 + y1*y2 + z1*z2;
    c = clamp(c, -1.0, 1.0);
    return std::acos(c) * 180.0 / TMath::Pi();
}

// Energies stored in the simulation are kinetic energies. For invariant mass
// use total energies: E_total = T + m_e and |p| = sqrt(T^2 + 2 T m_e).
double mass_ee_mev_from_kinetic(double kineticEeMeV,
                                double kineticEpMeV,
                                double thetaEEDeg)
{
    constexpr double me = 0.51099895; // MeV/c^2
    const double e1 = kineticEeMeV + me;
    const double e2 = kineticEpMeV + me;
    const double p1 = std::sqrt(std::max(0.0, kineticEeMeV * kineticEeMeV + 2.0 * kineticEeMeV * me));
    const double p2 = std::sqrt(std::max(0.0, kineticEpMeV * kineticEpMeV + 2.0 * kineticEpMeV * me));
    const double c12 = std::cos(thetaEEDeg * TMath::Pi() / 180.0);
    const double m2 = 2.0 * me * me + 2.0 * (e1 * e2 - p1 * p2 * c12);
    return std::sqrt(std::max(0.0, m2));
}
}

void apply_detector_effects(const char* inputFile  = "analysis_hits.root",
                            const char* outputFile = "analysis_hits_detector_effects.root",
                            double particleEfficiency = 0.90,
                            double sigmaThetaDeg = 2.0,
                            double sigmaPhiDeg = 2.0,
                            double relativeEnergyResolution = 0.05,
                            double energyThresholdMeV = 1.0,
                            unsigned int seed = 12345)
{
    // Input:  analysis_hits.root with tree `detected`.
    // Output: analysis_hits_detector_effects.root with tree `detected_detector_effects`.
    //
    // The detector model is deliberately simple and analysis-level:
    //   1) per-particle detection efficiency;
    //   2) Gaussian angular smearing on theta and phi;
    //   3) Gaussian relative energy smearing;
    //   4) threshold on reconstructed energies.

    if (particleEfficiency < 0.0 || particleEfficiency > 1.0) {
        std::cerr << "ERROR: particleEfficiency must be between 0 and 1.\n";
        return;
    }
    if (sigmaThetaDeg < 0.0 || sigmaPhiDeg < 0.0 ||
        relativeEnergyResolution < 0.0 || energyThresholdMeV < 0.0) {
        std::cerr << "ERROR: resolutions and threshold must be non-negative.\n";
        return;
    }

    TFile fin(inputFile, "READ");
    if (fin.IsZombie()) {
        std::cerr << "ERROR: cannot open input ROOT file: " << inputFile << std::endl;
        return;
    }

    auto* detectedIn = dynamic_cast<TTree*>(fin.Get("detected"));
    if (!detectedIn) {
        std::cerr << "ERROR: input file must contain a tree named detected.\n";
        std::cerr << "Run scripts/analyze_hits.sh before applying detector effects.\n";
        return;
    }

    int eventID = -1;
    int dataRow = -1;
    int selectedVolumeID = -1;
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

    detectedIn->SetBranchAddress("eventID", &eventID);
    detectedIn->SetBranchAddress("dataRow", &dataRow);
    detectedIn->SetBranchAddress("selectedVolumeID", &selectedVolumeID);
    detectedIn->SetBranchAddress("detIDEe", &detIDEe);
    detectedIn->SetBranchAddress("detIDEp", &detIDEp);
    detectedIn->SetBranchAddress("thetaEE_txt_deg", &thetaEE_txt_deg);
    detectedIn->SetBranchAddress("thetaEE_hit_deg", &thetaEE_hit_deg);
    detectedIn->SetBranchAddress("thetaEe_txt_deg", &thetaEe_txt_deg);
    detectedIn->SetBranchAddress("thetaEp_txt_deg", &thetaEp_txt_deg);
    detectedIn->SetBranchAddress("energyEe_txt_MeV", &energyEe_txt_MeV);
    detectedIn->SetBranchAddress("energyEp_txt_MeV", &energyEp_txt_MeV);
    detectedIn->SetBranchAddress("thetaEe_hit_deg", &thetaEe_hit_deg);
    detectedIn->SetBranchAddress("thetaEp_hit_deg", &thetaEp_hit_deg);
    detectedIn->SetBranchAddress("phiEe_hit_deg", &phiEe_hit_deg);
    detectedIn->SetBranchAddress("phiEp_hit_deg", &phiEp_hit_deg);
    detectedIn->SetBranchAddress("energyEe_hit_MeV", &energyEe_hit_MeV);
    detectedIn->SetBranchAddress("energyEp_hit_MeV", &energyEp_hit_MeV);
    detectedIn->SetBranchAddress("xEe_hit_mm", &xEe_hit_mm);
    detectedIn->SetBranchAddress("yEe_hit_mm", &yEe_hit_mm);
    detectedIn->SetBranchAddress("zEe_hit_mm", &zEe_hit_mm);
    detectedIn->SetBranchAddress("xEp_hit_mm", &xEp_hit_mm);
    detectedIn->SetBranchAddress("yEp_hit_mm", &yEp_hit_mm);
    detectedIn->SetBranchAddress("zEp_hit_mm", &zEp_hit_mm);

    TRandom3 rng(seed);

    TFile fout(outputFile, "RECREATE");

    TH1D hCutFlow("hCutFlow", "Detector effects cut flow;step;events", 4, 0.5, 4.5);
    hCutFlow.GetXaxis()->SetBinLabel(1, "ideal detected");
    hCutFlow.GetXaxis()->SetBinLabel(2, "pass e-/e+ efficiency");
    hCutFlow.GetXaxis()->SetBinLabel(3, "pass E thresholds");
    hCutFlow.GetXaxis()->SetBinLabel(4, "stored");

    TH1D hThetaEEtxt("hThetaEEtxt", "Generated opening angle for input detected events;#theta_{ee}^{txt} [deg];Events", 180, 0, 180);
    TH1D hThetaEEhit("hThetaEEhit", "Ideal hit opening angle;#theta_{ee}^{hit} [deg];Events", 180, 0, 180);
    TH1D hThetaEEreco("hThetaEEreco", "Opening angle after detector effects;#theta_{ee}^{reco} [deg];Events", 180, 0, 180);
    TH1D hDeltaThetaEE("hDeltaThetaEE", "Opening angle residual;#theta_{ee}^{reco}-#theta_{ee}^{hit} [deg];Events", 160, -20, 20);
    TH2D hThetaEEhitVsReco("hThetaEEhitVsReco", "Opening angle response;#theta_{ee}^{hit} [deg];#theta_{ee}^{reco} [deg]", 180, 0, 180, 180, 0, 180);

    TH1D hEnergyEeReco("hEnergyEeReco", "Electron reconstructed kinetic energy;T_{e-}^{reco} [MeV];Events", 200, 0, 20);
    TH1D hEnergyEpReco("hEnergyEpReco", "Positron reconstructed kinetic energy;T_{e+}^{reco} [MeV];Events", 200, 0, 20);
    TH1D hEnergyEEReco("hEnergyEEReco", "Reconstructed kinetic-energy sum;T_{e-}^{reco}+T_{e+}^{reco} [MeV];Events", 220, 0, 22);
    TH1D hMassEEReco("hMassEEReco", "Reconstructed invariant mass;m_{ee}^{reco} [MeV/c^{2}];Events", 220, 0, 22);
    TH1D hThetaEeReco("hThetaEeReco", "Electron reconstructed polar angle;#theta_{e-}^{reco} [deg];Events", 180, 0, 180);
    TH1D hThetaEpReco("hThetaEpReco", "Positron reconstructed polar angle;#theta_{e+}^{reco} [deg];Events", 180, 0, 180);

    double thetaEE_reco_deg = -999.0;
    double thetaEe_reco_deg = -999.0;
    double thetaEp_reco_deg = -999.0;
    double phiEe_reco_deg = -999.0;
    double phiEp_reco_deg = -999.0;
    double energyEe_reco_MeV = -999.0;
    double energyEp_reco_MeV = -999.0;
    double energyEE_hit_MeV = -999.0;
    double energyEE_reco_MeV = -999.0;
    double massEE_hit_MeV = -999.0;
    double massEE_reco_MeV = -999.0;

    bool passedEfficiency = false;
    bool passedThreshold = false;

    TTree detectedOut("detected_detector_effects", "one row per e-/e+ event after analysis-level detector effects");
    detectedOut.Branch("eventID", &eventID);
    detectedOut.Branch("dataRow", &dataRow);
    detectedOut.Branch("selectedVolumeID", &selectedVolumeID);
    detectedOut.Branch("detIDEe", &detIDEe);
    detectedOut.Branch("detIDEp", &detIDEp);

    detectedOut.Branch("particleEfficiency", &particleEfficiency);
    detectedOut.Branch("sigmaThetaDeg", &sigmaThetaDeg);
    detectedOut.Branch("sigmaPhiDeg", &sigmaPhiDeg);
    detectedOut.Branch("relativeEnergyResolution", &relativeEnergyResolution);
    detectedOut.Branch("energyThresholdMeV", &energyThresholdMeV);

    detectedOut.Branch("passedEfficiency", &passedEfficiency);
    detectedOut.Branch("passedThreshold", &passedThreshold);

    detectedOut.Branch("thetaEE_txt_deg", &thetaEE_txt_deg);
    detectedOut.Branch("thetaEE_hit_deg", &thetaEE_hit_deg);
    detectedOut.Branch("thetaEE_reco_deg", &thetaEE_reco_deg);

    detectedOut.Branch("thetaEe_txt_deg", &thetaEe_txt_deg);
    detectedOut.Branch("thetaEp_txt_deg", &thetaEp_txt_deg);
    detectedOut.Branch("energyEe_txt_MeV", &energyEe_txt_MeV);
    detectedOut.Branch("energyEp_txt_MeV", &energyEp_txt_MeV);

    detectedOut.Branch("thetaEe_hit_deg", &thetaEe_hit_deg);
    detectedOut.Branch("thetaEp_hit_deg", &thetaEp_hit_deg);
    detectedOut.Branch("phiEe_hit_deg", &phiEe_hit_deg);
    detectedOut.Branch("phiEp_hit_deg", &phiEp_hit_deg);
    detectedOut.Branch("energyEe_hit_MeV", &energyEe_hit_MeV);
    detectedOut.Branch("energyEp_hit_MeV", &energyEp_hit_MeV);

    detectedOut.Branch("thetaEe_reco_deg", &thetaEe_reco_deg);
    detectedOut.Branch("thetaEp_reco_deg", &thetaEp_reco_deg);
    detectedOut.Branch("phiEe_reco_deg", &phiEe_reco_deg);
    detectedOut.Branch("phiEp_reco_deg", &phiEp_reco_deg);
    detectedOut.Branch("energyEe_reco_MeV", &energyEe_reco_MeV);
    detectedOut.Branch("energyEp_reco_MeV", &energyEp_reco_MeV);
    detectedOut.Branch("energyEE_hit_MeV", &energyEE_hit_MeV);
    detectedOut.Branch("energyEE_reco_MeV", &energyEE_reco_MeV);
    detectedOut.Branch("massEE_hit_MeV", &massEE_hit_MeV);
    detectedOut.Branch("massEE_reco_MeV", &massEE_reco_MeV);

    detectedOut.Branch("xEe_hit_mm", &xEe_hit_mm);
    detectedOut.Branch("yEe_hit_mm", &yEe_hit_mm);
    detectedOut.Branch("zEe_hit_mm", &zEe_hit_mm);
    detectedOut.Branch("xEp_hit_mm", &xEp_hit_mm);
    detectedOut.Branch("yEp_hit_mm", &yEp_hit_mm);
    detectedOut.Branch("zEp_hit_mm", &zEp_hit_mm);

    const Long64_t nIdealDetected = detectedIn->GetEntries();
    Long64_t nPassEfficiency = 0;
    Long64_t nPassThreshold = 0;
    Long64_t nStored = 0;

    for (Long64_t i = 0; i < nIdealDetected; ++i) {
        detectedIn->GetEntry(i);

        hCutFlow.Fill(1);
        hThetaEEtxt.Fill(thetaEE_txt_deg);
        hThetaEEhit.Fill(thetaEE_hit_deg);
        energyEE_hit_MeV = energyEe_hit_MeV + energyEp_hit_MeV;
        massEE_hit_MeV = mass_ee_mev_from_kinetic(energyEe_hit_MeV, energyEp_hit_MeV, thetaEE_hit_deg);

        passedEfficiency = false;
        passedThreshold = false;

        const bool keepEe = rng.Uniform() <= particleEfficiency;
        const bool keepEp = rng.Uniform() <= particleEfficiency;
        if (!keepEe || !keepEp) continue;

        passedEfficiency = true;
        ++nPassEfficiency;
        hCutFlow.Fill(2);

        thetaEe_reco_deg = clamp(thetaEe_hit_deg + rng.Gaus(0.0, sigmaThetaDeg), 0.0, 180.0);
        thetaEp_reco_deg = clamp(thetaEp_hit_deg + rng.Gaus(0.0, sigmaThetaDeg), 0.0, 180.0);
        phiEe_reco_deg = wrap_phi_deg(phiEe_hit_deg + rng.Gaus(0.0, sigmaPhiDeg));
        phiEp_reco_deg = wrap_phi_deg(phiEp_hit_deg + rng.Gaus(0.0, sigmaPhiDeg));

        energyEe_reco_MeV = energyEe_hit_MeV + rng.Gaus(0.0, relativeEnergyResolution * energyEe_hit_MeV);
        energyEp_reco_MeV = energyEp_hit_MeV + rng.Gaus(0.0, relativeEnergyResolution * energyEp_hit_MeV);
        energyEe_reco_MeV = std::max(0.0, energyEe_reco_MeV);
        energyEp_reco_MeV = std::max(0.0, energyEp_reco_MeV);

        if (energyEe_reco_MeV < energyThresholdMeV || energyEp_reco_MeV < energyThresholdMeV) continue;

        passedThreshold = true;
        ++nPassThreshold;
        hCutFlow.Fill(3);

        thetaEE_reco_deg = opening_angle_deg_from_angles(thetaEe_reco_deg, phiEe_reco_deg,
                                                        thetaEp_reco_deg, phiEp_reco_deg);
        energyEE_reco_MeV = energyEe_reco_MeV + energyEp_reco_MeV;
        massEE_reco_MeV = mass_ee_mev_from_kinetic(energyEe_reco_MeV, energyEp_reco_MeV, thetaEE_reco_deg);

        hThetaEEreco.Fill(thetaEE_reco_deg);
        hDeltaThetaEE.Fill(thetaEE_reco_deg - thetaEE_hit_deg);
        hThetaEEhitVsReco.Fill(thetaEE_hit_deg, thetaEE_reco_deg);
        hEnergyEeReco.Fill(energyEe_reco_MeV);
        hEnergyEpReco.Fill(energyEp_reco_MeV);
        hEnergyEEReco.Fill(energyEE_reco_MeV);
        hMassEEReco.Fill(massEE_reco_MeV);
        hThetaEeReco.Fill(thetaEe_reco_deg);
        hThetaEpReco.Fill(thetaEp_reco_deg);

        detectedOut.Fill();
        ++nStored;
        hCutFlow.Fill(4);
    }

    hCutFlow.Write();
    hThetaEEtxt.Write();
    hThetaEEhit.Write();
    hThetaEEreco.Write();
    hDeltaThetaEE.Write();
    hThetaEEhitVsReco.Write();
    hEnergyEeReco.Write();
    hEnergyEpReco.Write();
    hEnergyEEReco.Write();
    hMassEEReco.Write();
    hThetaEeReco.Write();
    hThetaEpReco.Write();
    detectedOut.Write();

    fout.Close();

    const double effFromGeneratedPlaceholder = 0.0; // Printed summary is relative to ideal detected input.
    (void)effFromGeneratedPlaceholder;

    std::cout << "\n=== X17 detector-effects analysis summary ===\n";
    std::cout << "Input analysis file          : " << inputFile << "\n";
    std::cout << "Output detector-effects file : " << outputFile << "\n";
    std::cout << "particleEfficiency           : " << particleEfficiency << " per particle\n";
    std::cout << "sigmaThetaDeg                : " << sigmaThetaDeg << " deg\n";
    std::cout << "sigmaPhiDeg                  : " << sigmaPhiDeg << " deg\n";
    std::cout << "relativeEnergyResolution     : " << relativeEnergyResolution << "\n";
    std::cout << "energyThresholdMeV           : " << energyThresholdMeV << " MeV\n";
    std::cout << "seed                         : " << seed << "\n";
    std::cout << "Ideal detected coincidences  : " << nIdealDetected << "\n";
    std::cout << "Pass e-/e+ efficiency        : " << nPassEfficiency
              << "  fraction = " << (nIdealDetected > 0 ? double(nPassEfficiency)/double(nIdealDetected) : 0.0) << "\n";
    std::cout << "Pass energy thresholds       : " << nPassThreshold
              << "  fraction = " << (nIdealDetected > 0 ? double(nPassThreshold)/double(nIdealDetected) : 0.0) << "\n";
    std::cout << "Stored reconstructed events  : " << nStored
              << "  fraction = " << (nIdealDetected > 0 ? double(nStored)/double(nIdealDetected) : 0.0) << "\n";
    std::cout << "Tree written                 : detected_detector_effects\n";
    std::cout << "=============================================\n";
}

{
    gSystem->AddIncludePath("-Iscripts");

    Int_t error = 0;

    gROOT->ProcessLine(".L scripts/plot_angles_ROOT.C", &error);

    if (error != 0) {
        std::cerr << "ERROR: could not load scripts/plot_angles_ROOT.C" << std::endl;
        gSystem->Exit(30);
    }

    gROOT->ProcessLine("plot_angles_ROOT(\"data/data_pair_creation.txt\", \"plots\");", &error);

    if (error != 0) {
        std::cerr << "ERROR: plot_angles_ROOT failed" << std::endl;
        gSystem->Exit(31);
    }
}

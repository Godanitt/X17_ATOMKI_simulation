{
    gSystem->AddIncludePath("-Iscripts");

    Int_t error = 0;

    gROOT->ProcessLine(".L scripts/study_detector_effects.C", &error);

    if (error != 0) {
        std::cerr << "ERROR: could not load scripts/study_detector_effects.C" << std::endl;
        gSystem->Exit(20);
    }

    gROOT->ProcessLine("study_detector_effects(\"x17_output.root\", \"plots\");", &error);

    if (error != 0) {
        std::cerr << "ERROR: detector-effects macro failed" << std::endl;
        gSystem->Exit(21);
    }
}

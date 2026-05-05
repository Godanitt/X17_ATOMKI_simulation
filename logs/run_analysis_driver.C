{
    gSystem->AddIncludePath("-Iscripts");

    Int_t error = 0;
    gROOT->ProcessLine(".L scripts/analyze_X17.C", &error);

    if (error != 0) {
        std::cerr << "ERROR: could not load scripts/analyze_X17.C" << std::endl;
        gSystem->Exit(10);
    }

    gROOT->ProcessLine("analyze_x17(\"x17_output.root\", \"plots\");", &error);

    if (error != 0) {
        std::cerr << "ERROR: analyze_x17 failed" << std::endl;
        gSystem->Exit(11);
    }
}

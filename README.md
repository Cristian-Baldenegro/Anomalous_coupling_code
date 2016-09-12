# Anomalous_coupling_code
C++ code to analyze FPMC output root files for AA->XY anomalous coupling (XY=AZ, ZZ,WW,AA)

To run:
-make clean
-make
-./aaAnom_signal (Outputs ROOT file with histograms that you can use for your analysis)

Workflow:
-You need an FPMC tmpntuple.ntp transformed to a root file with the command h2root tmpntuple.ntp yourfilename.root . For instructions to run FPMC see our Twiki page https://twiki.cern.ch/twiki/bin/view/Sandbox/ExclusiveDiphotonAnalysis#Forward_Physics_Monte_Carlo_FPMC

-aaAnom_signal.cpp reads the FPMC output rootfile, calls the relevant routines in aaAnom_analysis.h and feeds in the TTree with aAnom_tree.h.
Please note that you should select the integrated luminosity you're using for your analysis and the number of events you used for your simulation in FPMC. This information is needed for proper the normalization weight for your simulation.
-aaAnom_analysis.h is where the analysis is done. Its output is a rootfile with histograms containing Pt, invariant mass,
eta, phi distributions, what have you. You can add more observables as you see fit.
-The output root file can be used for further analysis (e.g., find correlations between observables or compare different FPMC simulations with different parameters) or you can directly browse it.


outPileup.root is Matthias Saimpert's previous Pileup simulation rootfile with ATLAS Forward Physics (APS). This needs to be corrected for the TOTEM case, but will be used provisionally. Also, the reconstruction efficiency is assumed for the ATLAS central detectors. We need to change this for CMS central detector (We could treat a reconstructino efficiency of almost 100% for muons and electrons for a Pt>10).

In the case of jets, the fastjet library is used and the anti-kt algorithm is used. You can learn more about fastjet here http://fastjet.fr/ . You don't need to use this for AA->AA, AA->WW analysis, since they don't decay hadronically.

In case you want to run the analysis code on your local computer, you'll need:
-ROOT
-Fastjet (If you're analyzing jet decays)

If you want to run FPMC on your local computer, you'll need gfortran and CERNLIB. My personal recommendation is to run FPMC on Lxplus, transform the NTUPLE file to a ROOT file and do the analysis in your local computer (It's faster to edit your code locally for the analysis), feel free to do as you wish.

Any questions regarding the code: Cristian Baldenegro cbaldenegro@ku.edu

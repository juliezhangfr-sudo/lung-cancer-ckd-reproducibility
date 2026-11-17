# Run from the repository root with Rscript analysis/run_all.R
source('analysis/additional_checks.R')
source('analysis/round2_sensitivities.R')
source('analysis/round2_figures.R')
source('analysis/continue_native_checks.R')
cat('R workflows completed. Fresh NCI execution is documented separately in native exports.\n')

# Lung cancer mortality and CKD co-recording in the United States

Reproducibility materials for the 1999–2024 CDC WONDER analysis. The materials contain aggregate data, not individual death records.

- [Publication figures](publication/figure-captions.md), [result tables](publication/table-index.md), and [supplementary results](publication/supplementary-results.md) were extracted from the supplied November 16 manuscript and updated supplementary Word document, using the final visible text of tracked revisions. Drafts, editorial comments, and reviewer correspondence are not included.
- [Analysis inputs](analysis/input) contain counts, rates, standard errors, and comparisons required by the scripts.
- [CDC WONDER query requests and responses](analysis/source/new_cdc_queries) preserve the saved XML specifications and responses. Additional [original exports](analysis/source/wonder_desktop) and [urbanicity exports](analysis/source/urbanicity_wonder) retain source reporting and query notes. See the [query guide](query-guide.md) for scope. Treat suppressed or unavailable observations as missing, not zero.
- [Joinpoint sessions and exports](analysis/source/joinpoint_native), [age-specific sessions](analysis/source/age_audit), and [regional verification sessions](analysis/source/native_continued) preserve the native inputs and output. Use NCI Joinpoint Regression Program 5.4.0.0; if prompted for an original local path, point Joinpoint to the matching bundled rate-and-SE input.
- [Saved analysis results](analysis/results) contain machine-readable outputs and regenerated figures, which are separate from the manuscript figures under `publication/`.

From the repository root, run `Rscript analysis/run_all.R` (tested with R 4.5.2). The R checks do not rerun native NCI model selection or confidence intervals. `python analysis/prepare_inputs.py` optionally regenerates selected intermediate inputs from bundled source exports; the ready-to-use inputs are included. Run that step in a separate copy if you want to preserve the provided input files.

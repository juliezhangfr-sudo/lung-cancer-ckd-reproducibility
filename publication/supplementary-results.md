# Supplementary results

Supplementary Table S1. National age-standardized CKD co-recording analysis, 2013–2024. Panel A. Change in CKD co-recording, 2013–2024

[Table S1, block 1](../publication/tables/table-S1-1.csv)

Panel B. Adjusted temporal trend

[Table S1, block 2](../publication/tables/table-S1-2.csv)

Year × age interaction: global Wald P = .094. This test was not statistically significant, so differences between age groups in the temporal trend are not supported by formal inference.

Grouped quasi-binomial logit model: logit[P(CKD co-recorded)] = year_c + age group, where year_c=year−2013. To retain all known ages despite suppression in finer younger-age CKD cells, ages <55 years were collapsed into one stratum; the remaining strata were 55–64, 65–74, 75–84, and ≥85 years. Decedents with age not stated were not assignable to an age stratum and were excluded (17 of 1,704,558 lung-cancer deaths and 1 of 18,366 Lung+CKD deaths during 2013–2024). The age-standardized model therefore included 1,704,541 lung-cancer deaths, of which 18,365 had CKD co-recorded. Standardization used the pooled 2013–2024 distribution of lung-cancer deaths across these five strata (weights 0.05308, 0.19469, 0.32871, 0.29447, and 0.12906, respectively). Pearson dispersion =

2.426. Annual aggregate Pearson residual lag-1 correlation r=.590 (P = .056); the lag-1 HAC sensitivity interval retained statistical significance. A categorical-year sensitivity model used indicator terms for calendar year (2013 reference) with the same categorical age groups and fixed pooled age weights; it yielded standardized predictions of 0.727% in 2013 and 1.374% in 2024, with modest year-to-year fluctuations, and was used to assess model-form sensitivity without imposing a linear year effect. Model structure: the model contained 60 grouped observations (12 calendar years × 5 age strata), 6 estimated coefficients (intercept, year_c and four age contrasts) and 54 residual degrees of freedom; the quasi-binomial dispersion estimate (2.426) scales the primary covariance matrix. The estimand is the age-standardized (marginally standardized) CKD co-recording proportion among lung-cancer decedents, that is, a conditional probability standardized to a descriptive reference distribution; it is not a conventional population-based age-adjusted mortality rate. The pooled 2013–2024 decedent distribution was used as a neutral descriptive reference, and the 2013 and 2024 reference distributions are shown as sensitivity analyses in the table above. The lag-1 and lag-2 HAC analyses are sensitivity analyses based on the same fitted model and are not independent replications of the trend; a one-year lag was used as the shortest dependence structure available with 12 annual points.

Revision sensitivity analyses

[Table S1, block 3](../publication/tables/table-S1-3.csv)

HAC covariances use year-aggregated binomial score contributions, Bartlett lags 1 or 2, no prewhitening and no finite-sample correction. They are sensitivity estimates from the same model. National source counts were rechecked with official queries; the complete 120-count reconciliation is included in the analysis inputs. The absolute changes obtained with the 2013 reference distribution (+0.687 pp) and the 2024 reference distribution (+0.726 pp) are close to the pooled-reference change (+0.703 pp), indicating limited dependence of the descriptive result on the chosen standard.

Supplementary Table S2. Age-standardized CKD co-recording sensitivity analysis among lung-cancer decedents aged ≥55 years by six-level 2023 NCHS urbanicity, 2018–2024. “2023 NCHS” denotes the urban-rural classification-scheme vintage; mortality follow-up extends through 2024.

[Table S2, block 1](../publication/tables/table-S2-1.csv)

[Table S2, block 2](../publication/tables/table-S2-2.csv)

Age-stratified grouped quasi-binomial model with continuous year centered at 2018 (year_c=year−2018), year × urbanicity interaction, and categorical age group. Common pooled 2018–2024 age weights: 55–64=0.196996, 65–74=0.349229, 75–84=0.318113, ≥85=0.135663. Pearson dispersion=1.131 and was used as the variance scale.

Global year × urbanicity Wald χZ(5)=8.31, P = .140. Records coded as urbanicity Not Available were excluded. This is a sensitivity analysis; no pairwise subgroup superiority is inferred. Urbanicity categories are listed in NCHS classification order and are not ranked. All urbanicity analyses are exploratory, and no multiplicity adjustment was applied.

Supplementary Table S3. NCI Joinpoint audit and retained model metadata.

Panel A. Retained overall Joinpoint estimates for the age-adjusted mortality rate of lung-cancer deaths with CKD co-recorded

[Table S3, block 1](../publication/tables/table-S3-1.csv)

Panel B. Documented NCI Joinpoint model settings and evidence status

[Table S3, block 2](../publication/tables/table-S3-2.csv)

Note: APC, annual percent change; AAPC, average annual percent change; CI, confidence interval; SE, standard error. The retained analysis record documents the software version, model form, weighting, permutation-test selection, and empirical-quantile confidence-interval procedure. Native .JPO/.JPS files for the retained overall Lung+CKD model summarized in this table have been located. Item-by-item verification against the saved native outputs (Run Summary, Model Selection, APC and AAPC exports) is complete for the retained overall Lung+CKD model: Joinpoint V5.4.0.0; log-linear model weighted by CDC WONDER-provided standard errors with uncorrelated errors; grid search; maximum 2 joinpoints; minimum 2 observations from a joinpoint to either end of the data and between two joinpoints; permutation test with overall α = .05 and 4,499 permutations; empirical-quantile confidence intervals with 5,001 resamples. The same configuration was verified for the primary age, sex, race/ethnicity, region and urbanicity Joinpoint sessions (native_session_audit.csv); the maximum number of joinpoints was 4 for the 1999–2024 series, 2 for the 2013–2024 series and 1 for the 2013–2020 urbanicity Lung+CKD series. For the overall all-lung 1999–2024 series the original saved full-range AAPC display was absent, and that AAPC was verified by a fresh native rerun. This table concerns the retained overall model, not the updated age-specific models. For the updated age analyses, the retained input CSVs, native .JPO/.JPS files, Run Summary, Model Selection, APC/AAPC exports, and selected-model outputs confirmed the reported AAPCs and confidence intervals, maximum allowed joinpoint counts of 4 for all-lung mortality and 2 for Lung+CKD mortality, minimum-observation constraints, and handling of the unreliable 2015 observation in the 45–54-year Lung+CKD series; no discrepancy was identified. The NCI Help System documents software methods and options but does not establish which settings were used in a specific analysis session.

The +2.43% estimate reported in sensitivity analysis is a one-segment 2013–2024 APC and is not the same estimand as the +2.1350% AAPC from the selected segmented model with a 2020 joinpoint. The reported primary APC/AAPC values are retained-record NCI Joinpoint results; no unsupported session settings were reconstructed.

Model selection and robustness of the 2020 joinpoint (overall Lung+CKD series, native Model Selection export): the permutation test of 0 versus 2 joinpoints gave P < .001 (significance level .025), and the test of 1 versus 2 joinpoints gave P = .068 (level .050), so one joinpoint (2020) was selected. The 2020–2024 segment contains five annual observations including the joinpoint year, and its APC (−1.17%; 95% CI −4.95 to 0.74; P = .202) is compatible with continued increase, a plateau or a decline; it is treated as inconclusive. The one-segment 2013–2024 OLS APC (+2.43%; 95% CI 1.43 to 3.43) is shown beside the segmented model so that the dependence of the 2020 break on model choice is visible; a fixed-2020 hinge model reproduced the segment slopes (+4.07% before and −1.17% after 2020). The data-driven joinpoint is a descriptive change in slope. It is not evidence that a biological or epidemiologic transition occurred in 2020, and no causal attribution to the COVID-19 pandemic is made; any such explanation would be post hoc. A comparison of 2013–2019 versus 2021–2024 has not been performed because too few post-2020 observations are available.

Supplementary Table S4. Race-era sensitivity analyses and 2018–2020 overlap quality control. (unweighted OLS sensitivity analysis)

Unweighted ordinary least-squares (OLS) sensitivity analysis. These APCs are not directly comparable with the standard-error-weighted NCI Joinpoint AAPCs reported in Table 3 and Supplementary Tables S3 and S6, and the two eras use different race definitions. The estimates indicate direction and continuity only and are not a single harmonized longitudinal effect.

[Table S4, block 1](../publication/tables/table-S4-1.csv)

Race-era APCs were estimated by unweighted ordinary least-squares regression of log-transformed AAMRs on calendar year, with t-based confidence intervals and P values using residual degrees of freedom. The Bridged Race analysis covered 1999–2020 and therefore spans the 2011–2012 coding transition; the Single Race analysis covered 2018–2024. These sensitivity models are separate from the SE-weighted NCI Joinpoint analyses and should not be interpreted as proving absence of a coding or database-transition effect. Overall Lung+CKD death counts and AAMRs were identical between database vintages in 2018–2020. Bridged Race Other comprised non-Hispanic American Indian or Alaska Native plus non-Hispanic Asian or Pacific Islander; Single Race Other comprised non-Hispanic American Indian or Alaska Native, Asian, Native Hawaiian or Other Pacific Islander, and More than one race. The continuous cross-era Other AAPC is not interpreted as a primary longitudinal estimate. The positive estimates for non-Hispanic White persons in both eras do not demonstrate that the underlying racial trend is unaffected by the database change; the within-era estimates are the only inferential results, and the Bridged Race estimate additionally includes the 2011–2012 coding-transition period.

Supplementary Table S5. Updated 2024 age-specific descriptive endpoints.

[Table S5, block 1](../publication/tables/table-S5-1.csv)

These are 2024 endpoints from final Single Race CDC WONDER data. Formal age-specific NCI Joinpoint analyses were updated through 2024 (all-lung, 1999–2024; Lung+CKD, 2013–2024). The archived 2015 Lung+CKD rate for ages 45–54 was based on 17 deaths and was flagged unreliable; that observation was excluded from formal trend inference without interpolation. The CKD co-recording percentages are crude within each age group, are not age-standardized, and are descriptive only. Sensitivity of the 45–54-year result to the unreliable 2015 point: the Joinpoint input for this age group contained 11 valid annual observations (2013, 2014 and 2016–2024; the 2015 row was omitted from the input file, so the software received a shortened series and not an imputed or missing-value placeholder), whereas each other age group contained 12 observations. The primary Joinpoint AAPC for ages 45–54 years was −0.39% (95% CI, −6.54% to 6.27%; P = .842). In a separate exploratory single-segment WLS analysis, the APC changed from −0.39% with 2015 excluded to +0.44% with the count-derived 2015 rate included. The reversal shows that the direction of the point estimate was sensitive to inclusion of the 2015 observation; it does not establish an increase or a decline.

Supplementary Table S6. Race, region, and historical urbanicity trend estimates for the age-adjusted mortality rate of lung-cancer deaths with CKD co-recorded.

[Table S6, block 1](../publication/tables/table-S6-1.csv)

Note: The aggregated Other race/ethnicity category is not shown because a continuous cross-era AAPC is not compositionally comparable; era-specific sensitivity estimates are provided in Supplementary Table S4.

Subgroup AAPCs are reported individually. Formal parallelism or pairwise tests were not performed, so these estimates should not be interpreted as statistically different from one another.   All subgroup analyses are exploratory and no multiplicity adjustment was applied. A statistically significant AAPC in one subgroup alongside a non-significant AAPC in another does not establish that the subgroup trends differ, and the subgroup with the highest point estimate should not be read as having the largest true trend.

Supplementary Figures

Supplementary Figure S1. Age-standardized predicted CKD co-recording proportions by six-level 2023 NCHS urbanicity among lung-cancer decedents aged ≥55 years, 2018–2024. Predictions were standardized to a common pooled 2018–2024 age distribution (55–64, 65–74, 75–84, and ≥85 years). The global year-by-urbanicity interaction was not statistically significant (Wald P = .140). This analysis is presented as a sensitivity analysis; the lines should not be interpreted as formal pairwise superiority tests. “2023 NCHS” denotes the classification-scheme vintage, not the mortality-data endpoint. No formal evidence of differential trends was found (global year × urbanicity Wald P = .140). The lines are plotted in NCHS classification order, and their vertical ordering in 2024 must not be read as a ranking of urbanicity categories.

Supplementary Table S7  Revision comparator and renal code sensitivities

All outcomes use C34 as the underlying cause, ages ≥55 years and four age strata. Rows represent separate grouped quasi-binomial models over 2013–2024, with 48 observations, five coefficients and 43 residual degrees of freedom when complete. Proportions use the same pooled C34-decedent age weights. These analyses were added during revision and are exploratory. N18.5 is a mortality-code stage-5/end-stage category; N18.6 is not an available clinical-code field in these mortality data [18,23]. A joint renal-code union was queried once per record, avoiding addition of potentially overlapping code-specific counts. Co-recording of diabetes (E10–E14) also rose over the same period (annual OR 1.0390), so a general increase in how comorbidities are recorded on certificates cannot be excluded as an explanation for the rise in N18 co-recording.

[Table S7, block 1](../publication/tables/table-S7-1.csv)

HAC intervals, Pearson dispersion, event totals, query requests, complete XML responses, release-overlap checks and count reconciliation are included in the saved analysis outputs. A count below 20 is a descriptive sparse-count flag here, not an assertion of the current database rate-reliability rule. Models use disclosed counts rather than population rates. Different event definitions share decedents, so separate ORs do not constitute a test of differing slopes.

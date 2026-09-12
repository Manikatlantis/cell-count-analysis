# Contract

Frozen before implementation. Do not rename anything here without telling me.

## Source columns, verbatim from cell-count.csv

project, subject, condition, age, sex, treatment, response, sample,
sample_type, time_from_treatment_start, b_cell, cd8_t_cell, cd4_t_cell,
nk_cell, monocyte

## Confirmed value spellings, case sensitive

project      prj1, prj2, prj3
condition    melanoma, carcinoma, healthy
sex          M, F
treatment    miraclib, phauximab, none
response     yes, no, or NULL (1422 rows: all healthy + treatment none)
sample_type  PBMC, WB
time         0, 7, 14

## Database

File: cell_count.db in repo root.
Schema: sql/schema.sql.
Tables: projects, conditions, treatments, subjects, treatment_courses,
        samples, cell_populations, sample_counts.
Views:  v_sample_frequencies, v_sample_annotated, v_analysis.

Loader invariants, assert all of these in tests:
- sample_counts row count == 10500 * 5 == 52500
- samples == 10500, subjects == 3500, treatment_courses == 3500
- response is SQL NULL for untreated, never the string 'nan' or 'None'
- every subject has exactly 3 samples
- running load_data.py twice leaves all counts unchanged

## Function signatures

src/analysis/frequency.py
    get_frequencies(conn) -> DataFrame
        columns exactly: sample, total_count, population, count, percentage

src/analysis/stats.py
    responder_comparison(conn, baseline_only: bool = False) -> DataFrame
        cohort: condition=melanoma, treatment=miraclib, sample_type=PBMC
        baseline_only=True adds time_from_treatment_start = 0
        columns: population, n_yes, n_no, n_subj_yes, n_subj_no,
                 median_yes, median_no, median_diff, cliffs_delta,
                 u_statistic, p_value, p_adj, significant

src/analysis/subsets.py
    baseline_cohort(conn) -> DataFrame
        melanoma + miraclib + PBMC + time 0
    samples_per_project(conn) -> DataFrame
        Counts samples, one per row of v_sample_annotated. Not COUNT(*)
        over v_analysis, which holds five rows per sample.
        Projects with no sample in the cohort are kept at 0.
    subjects_by_response(conn) -> DataFrame     # COUNT(DISTINCT subject)
    subjects_by_sex(conn) -> DataFrame          # COUNT(DISTINCT subject)
    mean_bcell_melanoma_male_responders_baseline(conn) -> float
        Returns the mean RAW b_cell COUNT, rounded to 2 decimals.
        FILTERS: condition=melanoma, sex=M, response=yes, time=0
        NO sample_type filter. NO treatment filter.
        Expected slice: 485 rows, 485 subjects. Expected value: 10206.15

    mean_bcell_percentage_melanoma_male_responders_baseline(conn) -> float
        Same slice, mean relative frequency. Expected value: 9.99
        Reported as a secondary figure because "average number of B cells"
        does not state whether count or frequency is meant.

src/analysis/plots.py
    boxplot_responders(df, outpath) -> Path

## Outputs, written by run_analysis.py

outputs/summary_frequencies.csv
outputs/responder_stats_all_timepoints.csv
outputs/responder_stats_baseline.csv
outputs/boxplot_responders.png
outputs/part4_baseline_cohort.csv
outputs/part4_summary.json
    Keys, both rounded to 2 decimals:
      mean_bcell_count_melanoma_male_responders_baseline       10206.15
      mean_bcell_percentage_melanoma_male_responders_baseline  9.99
    Plus baseline_cohort_samples, baseline_cohort_subjects,
    samples_per_project, subjects_by_response, subjects_by_sex.

## Known cohort sizes, verified in profiling. Treat as regression tests.

melanoma                                  5175 rows
  + miraclib                              2655
  + PBMC                                  1968   (656 subjects)
  + time 0                                 656   (656 subjects)

Part 3 arms, all timepoints:  993 yes / 975 no samples, 331 / 325 subjects
Part 3 arms, baseline only:   331 yes / 325 no samples, 331 / 325 subjects
Part 4 final slice:           485 rows, 485 subjects

If any implementation produces a number other than these, it is a bug.

## Additional verified counts

courses with NULL response  : 474   (subjects: healthy + treatment 'none')
samples with NULL response  : 1422  (474 x 3)

## Long-format gotcha

v_analysis and v_sample_frequencies are 5 rows per sample.
Counting samples there requires COUNT(DISTINCT sample), never COUNT(*).
Counting subjects requires COUNT(DISTINCT subject).
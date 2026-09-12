"""Independent recomputation of every published number, from cell-count.csv only.

Run as: python verify.py

This exists to catch a bug that the pipeline would otherwise reproduce
consistently and therefore hide. It recomputes the assignment from the CSV by a
completely separate route and diffs the result against the committed files in
outputs/.

Constraints, held on purpose so the two routes share no code:
  - Imports only pandas, numpy and scipy. Nothing from src/.
  - Reads cell-count.csv directly.
  - Never opens cell_count.db.
  - Never reads any .sql file.
  - Never reads load_data.py or run_analysis.py.

Where a check could be satisfied by copying the shipped approach, it is done a
different way:
  - Cliff's delta is counted pair by pair, not derived from the U statistic.
  - Benjamini-Hochberg is implemented here rather than called from statsmodels.
scipy's mannwhitneyu supplies the p value because the assignment names that exact
test, so reimplementing its tie handling would test scipy rather than this repo.

The tests here are fed full precision percentages computed from the CSV. The
pipeline feeds its tests the percentages as the summary table reports them,
which the assignment specifies and which are rounded to 4 decimals. That is a
real difference in input, so small differences in continuous values are expected
and are not failures. Anything categorical or countable must still agree exactly.

A difference fails this script only if it is one of:
  - any count, size or integer answer differing at all
  - any magnitude label differing
  - any significant flag differing
  - any p_value, p_adj, cliffs_delta or median differing by more than TOLERANCE

Everything else is reported under "within tolerance". The same tests are also
run on 4 decimal input as a secondary check, which prints but does not affect
the exit code.

Reports discrepancies. Fixes nothing.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent  # verify.py lives in the repo root
CSV = ROOT / "cell-count.csv"
OUT = ROOT / "outputs"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

# Half of the last digit the summary table keeps is 5e-5. A thousandth leaves
# room for that to propagate through a median or a rank test without hiding a
# difference that would matter.
TOLERANCE = 1e-3

checked = 0
failures: list[str] = []
tolerated: list[str] = []
rank_notes: list[str] = []


def compare(label: str, mine, theirs, kind: str = "exact") -> None:
    """Record one comparison.

    kind "exact"      : any difference is a failure.
    kind "continuous" : a difference up to TOLERANCE is reported, not failed.
    kind "rank"       : reported, never failed. See the u_statistic note below.
    """
    global checked
    checked += 1
    if kind == "exact":
        if isinstance(mine, (int, np.integer)) and isinstance(theirs, (int, np.integer)):
            ok = int(mine) == int(theirs)
        else:
            ok = mine == theirs
        if not ok:
            failures.append(f"{label}\n    verify.py: {mine!r}\n    outputs/ : {theirs!r}")
        return

    delta = abs(float(mine) - float(theirs))
    if kind == "rank":
        # U is a rank sum on a scale of n_yes * n_no, and a tie moves it by half
        # a rank, so an absolute tolerance set for a probability does not apply.
        # Cliff's delta is the same statistic rescaled to [-1, 1] and it is
        # checked against TOLERANCE above, which is the meaningful test.
        if delta > 0:
            rank_notes.append(f"{label:<44} {mine!r} vs {theirs!r}  (delta {delta:g})")
        return
    if delta > TOLERANCE:
        failures.append(
            f"{label}  (differs by {delta:.3e}, over tolerance)\n"
            f"    verify.py: {mine!r}\n    outputs/ : {theirs!r}"
        )
    elif delta > 0:
        tolerated.append(f"{label:<44} {mine!r} vs {theirs!r}  (delta {delta:.2e})")


# --- Independent computation -------------------------------------------------

raw = pd.read_csv(CSV)

long = raw.melt(
    id_vars=[c for c in raw.columns if c not in POPULATIONS],
    value_vars=POPULATIONS,
    var_name="population",
    value_name="count",
)
long["total_count"] = long.groupby("sample")["count"].transform("sum")
long["percentage"] = 100.0 * long["count"] / long["total_count"]
# Only for the secondary check further down. The primary path never uses this.
long["percentage_4dp"] = long["percentage"].round(4)

cohort = long[
    (long["condition"] == "melanoma")
    & (long["treatment"] == "miraclib")
    & (long["sample_type"] == "PBMC")
    & (long["response"].isin(["yes", "no"]))
]
baseline_cohort_long = cohort[cohort["time_from_treatment_start"] == 0]


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """P(x > y) - P(x < y), counted directly rather than via the U statistic."""
    gt = lt = 0
    for value in x:
        gt += int((value > y).sum())
        lt += int((value < y).sum())
    return (gt - lt) / (x.size * y.size)


def magnitude(delta: float) -> str:
    d = abs(delta)
    if d < 0.147:
        return "negligible"
    if d < 0.33:
        return "small"
    if d < 0.474:
        return "medium"
    return "large"


def benjamini_hochberg(pvals: list[float]) -> list[float]:
    """Step up BH, written out rather than imported."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adjusted = [0.0] * m
    running = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        position = m - rank + 1
        running = min(running, pvals[idx] * m / position)
        adjusted[idx] = running
    return adjusted


def comparison(frame: pd.DataFrame, column: str = "percentage") -> pd.DataFrame:
    rows = []
    for population in POPULATIONS:
        sub = frame[frame["population"] == population]
        yes = sub[sub["response"] == "yes"]
        no = sub[sub["response"] == "no"]
        x = yes[column].to_numpy(float)
        y = no[column].to_numpy(float)
        u = mannwhitneyu(x=x, y=y, alternative="two-sided")
        delta = cliffs_delta(x, y)
        rows.append(
            {
                "population": population,
                "n_yes": int(yes["sample"].nunique()),
                "n_no": int(no["sample"].nunique()),
                "n_subj_yes": int(yes["subject"].nunique()),
                "n_subj_no": int(no["subject"].nunique()),
                "median_yes": float(np.median(x)),
                "median_no": float(np.median(y)),
                "median_diff": float(np.median(x)) - float(np.median(y)),
                "cliffs_delta": delta,
                "magnitude": magnitude(delta),
                "u_statistic": float(u.statistic),
                "p_value": float(u.pvalue),
            }
        )
    out = pd.DataFrame(rows)
    out["p_adj"] = benjamini_hochberg(out["p_value"].tolist())
    out["significant"] = out["p_adj"] < 0.05
    return out


mine_all = comparison(cohort)
mine_base = comparison(baseline_cohort_long)

base_samples = raw[
    (raw["condition"] == "melanoma")
    & (raw["treatment"] == "miraclib")
    & (raw["sample_type"] == "PBMC")
    & (raw["time_from_treatment_start"] == 0)
]
mine_per_project = {p: 0 for p in sorted(raw["project"].unique())}
mine_per_project.update(base_samples.groupby("project")["sample"].count().to_dict())
mine_by_response = base_samples.groupby("response")["subject"].nunique().to_dict()
mine_by_sex = base_samples.groupby("sex")["subject"].nunique().to_dict()

q_slice = raw[
    (raw["condition"] == "melanoma")
    & (raw["sex"] == "M")
    & (raw["response"] == "yes")
    & (raw["time_from_treatment_start"] == 0)
]
mine_mean_count = round(float(q_slice["b_cell"].mean()), 2)
q_long = long[
    (long["condition"] == "melanoma")
    & (long["sex"] == "M")
    & (long["response"] == "yes")
    & (long["time_from_treatment_start"] == 0)
    & (long["population"] == "b_cell")
]
mine_mean_pct = round(float(q_long["percentage"].mean()), 2)


# --- Diff against outputs/ ---------------------------------------------------

print("Recomputed from cell-count.csv on full precision percentages.")
print(f"Comparing against outputs/. Tolerance for continuous values: {TOLERANCE:g}.\n")

theirs_freq = pd.read_csv(OUT / "summary_frequencies.csv")
mine_freq = long[["sample", "population", "count", "total_count", "percentage"]]
merged = theirs_freq.merge(
    mine_freq, on=["sample", "population"], how="outer",
    suffixes=("_theirs", "_mine"), indicator=True,
)
compare("part2.row_count", len(mine_freq), len(theirs_freq))
compare("part2.rows_matched_on_key", 0, int((merged["_merge"] != "both").sum()))

for col in ("count", "total_count"):
    diff = (merged[f"{col}_mine"] - merged[f"{col}_theirs"]).abs()
    checked += len(merged)
    if int((diff > 0).sum()):
        row = merged[diff > 0].iloc[0]
        failures.append(
            f"part2.{col}: {int((diff > 0).sum())} rows differ, first "
            f"{row['sample']}/{row['population']}\n"
            f"    verify.py: {row[f'{col}_mine']!r}\n    outputs/ : {row[f'{col}_theirs']!r}"
        )

# Bulk column, so reported as an aggregate rather than 52,500 separate lines.
pct_diff = (merged["percentage_mine"] - merged["percentage_theirs"]).abs()
checked += len(merged)
over = int((pct_diff > TOLERANCE).sum())
if over:
    failures.append(f"part2.percentage: {over} rows differ by more than {TOLERANCE:g}")
elif float(pct_diff.max()) > 0:
    tolerated.append(
        f"{'part2.percentage (all 52,500 rows)':<44} max delta {pct_diff.max():.2e}, "
        f"{int((pct_diff > 0).sum())} rows nonzero"
    )

# Gated at TOLERANCE. u_statistic is handled separately, see compare().
CONTINUOUS = {"median_yes", "median_no", "median_diff", "cliffs_delta", "p_value", "p_adj"}
for name, mine in (("all_timepoints", mine_all), ("baseline", mine_base)):
    theirs = pd.read_csv(OUT / f"responder_stats_{name}.csv").set_index("population")
    m = mine.set_index("population")
    compare(f"{name}.population_set", sorted(m.index), sorted(theirs.index))
    for population in POPULATIONS:
        a, b = m.loc[population], theirs.loc[population]
        for col in ("n_yes", "n_no", "n_subj_yes", "n_subj_no"):
            compare(f"{name}.{population}.{col}", int(a[col]), int(b[col]))
        for col in sorted(CONTINUOUS):
            if col in theirs.columns:
                compare(f"{name}.{population}.{col}", float(a[col]), float(b[col]), "continuous")
        if "u_statistic" in theirs.columns:
            compare(f"{name}.{population}.u_statistic",
                    float(a["u_statistic"]), float(b["u_statistic"]), "rank")
        if "magnitude" in theirs.columns:
            compare(f"{name}.{population}.magnitude", a["magnitude"], b["magnitude"])
        compare(f"{name}.{population}.significant", bool(a["significant"]), bool(b["significant"]))

summary = json.loads((OUT / "part4_summary.json").read_text())
compare("part4.baseline_cohort_samples", int(base_samples["sample"].nunique()),
        int(summary["baseline_cohort_samples"]))
compare("part4.baseline_cohort_subjects", int(base_samples["subject"].nunique()),
        int(summary["baseline_cohort_subjects"]))

theirs_projects = summary["samples_per_project"]
compare("part4.samples_per_project.keys", sorted(mine_per_project), sorted(theirs_projects))
for project in sorted(set(mine_per_project) | set(theirs_projects)):
    compare(f"part4.samples_per_project[{project}]",
            int(mine_per_project.get(project, 0)), int(theirs_projects.get(project, 0)))

for label, mine_map, theirs_map in (
    ("subjects_by_response", mine_by_response, summary["subjects_by_response"]),
    ("subjects_by_sex", mine_by_sex, summary["subjects_by_sex"]),
):
    compare(f"part4.{label}.keys", sorted(mine_map), sorted(theirs_map))
    for key in sorted(set(mine_map) | set(theirs_map)):
        compare(f"part4.{label}[{key}]", int(mine_map.get(key, 0)), int(theirs_map.get(key, 0)))

compare("part4.slice_rows", int(len(q_slice)), 485)
compare("part4.mean_bcell_count", mine_mean_count,
        float(summary["mean_bcell_count_melanoma_male_responders_baseline"]), "continuous")
compare("part4.mean_bcell_percentage", mine_mean_pct,
        float(summary["mean_bcell_percentage_melanoma_male_responders_baseline"]), "continuous")

theirs_base_csv = pd.read_csv(OUT / "part4_baseline_cohort.csv")
compare("part4_baseline_cohort.row_count", int(len(base_samples)), int(len(theirs_base_csv)))
compare("part4_baseline_cohort.sample_set",
        sorted(base_samples["sample"]), sorted(theirs_base_csv["sample"]))


# --- Secondary check: the same tests on 4 decimal input ----------------------

print("Secondary check. The same tests fed 4 decimal percentages, which is what")
print("the pipeline uses. This prints only and does not affect the exit code.")
secondary_exact = 0
secondary_total = 0
for name, frame in (("all_timepoints", cohort), ("baseline", baseline_cohort_long)):
    theirs = pd.read_csv(OUT / f"responder_stats_{name}.csv").set_index("population")
    rounded = comparison(frame, column="percentage_4dp").set_index("population")
    for population in POPULATIONS:
        for col in ("median_yes", "median_no", "cliffs_delta", "u_statistic", "p_value", "p_adj"):
            if col in theirs.columns:
                secondary_total += 1
                if abs(float(rounded.loc[population, col]) - float(theirs.loc[population, col])) <= 1e-9:
                    secondary_exact += 1
print(f"  {secondary_exact}/{secondary_total} values reproduce outputs/ exactly at 1e-9")
print("  A shortfall here would mean the difference is not explained by rounding.\n")


# --- Report ------------------------------------------------------------------

print(f"values checked      : {checked:,}")
print(f"failures            : {len(failures)}")
print(f"within tolerance    : {len(tolerated)}")
print(f"rank sum notes      : {len(rank_notes)}")

if tolerated:
    print(f"\nWithin tolerance, expected from 4 decimal rounding in the summary table.")
    print(f"Reported for completeness. None exceeds {TOLERANCE:g} absolute.\n")
    for line in tolerated:
        print(f"  {line}")

if rank_notes:
    print("\nU statistic differences. Not gated on an absolute tolerance, because U")
    print("is a rank sum on a scale of n_yes * n_no and each tie moves it half a")
    print("rank. Cliff's delta is the same statistic rescaled to [-1, 1] and it is")
    print("checked against the tolerance above.\n")
    for line in rank_notes:
        print(f"  {line}")

if failures:
    print(f"\nFAILURES. These are substantive and not explained by rounding.\n")
    for line in failures:
        print(line)
    raise SystemExit(1)

print("\nEvery count, label and flag agrees exactly.")
print("Continuous values agree within tolerance.")

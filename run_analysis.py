"""Run the whole analysis and write everything under outputs/.

No arguments. Exits non-zero if any cohort size drifts from the sizes frozen in
CONTRACT.md, so a silent filter change fails the pipeline instead of the review.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

from src.analysis.db import connect
from src.analysis.frequency import get_frequencies
from src.analysis.plots import boxplot_responders
from src.analysis.stats import cohort_frame, responder_comparison
from src.analysis.subsets import (
    baseline_cohort,
    mean_bcell_melanoma_male_responders_baseline,
    mean_bcell_percentage_melanoma_male_responders_baseline,
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)

REPO_ROOT = Path(__file__).resolve().parent
OUTPUTS = REPO_ROOT / "outputs"

_MELANOMA = "condition = 'melanoma'"
_MIRACLIB = f"{_MELANOMA} AND treatment = 'miraclib'"
_PBMC = f"{_MIRACLIB} AND sample_type = 'PBMC'"
_BASE = f"{_PBMC} AND time_from_treatment_start = 0"


class CountMismatch(Exception):
    pass


def _scalar(conn, sql: str) -> int:
    return int(conn.execute(sql).fetchone()[0])


def check_cohort_sizes(conn) -> list[str]:
    """Re-derive every frozen count straight from the views. Returns failures."""
    checks: list[tuple[str, int, int]] = []

    def add(label: str, sql: str, expected: int) -> None:
        checks.append((label, _scalar(conn, sql), expected))

    add("melanoma samples", f"SELECT COUNT(*) FROM v_sample_annotated WHERE {_MELANOMA}", 5175)
    add("  + miraclib", f"SELECT COUNT(*) FROM v_sample_annotated WHERE {_MIRACLIB}", 2655)
    add("  + PBMC", f"SELECT COUNT(*) FROM v_sample_annotated WHERE {_PBMC}", 1968)
    add("  + PBMC subjects", f"SELECT COUNT(DISTINCT subject) FROM v_sample_annotated WHERE {_PBMC}", 656)
    add("  + time 0", f"SELECT COUNT(*) FROM v_sample_annotated WHERE {_BASE}", 656)
    add("  + time 0 subjects", f"SELECT COUNT(DISTINCT subject) FROM v_sample_annotated WHERE {_BASE}", 656)

    arms = {("all timepoints", _PBMC): {"yes": (993, 331), "no": (975, 325)},
            ("baseline", _BASE): {"yes": (331, 331), "no": (325, 325)}}
    for (scope, where), expected in arms.items():
        for response, (n_samples, n_subjects) in expected.items():
            add(f"part 3 {scope} {response} samples",
                f"SELECT COUNT(*) FROM v_sample_annotated WHERE {where} AND response = '{response}'",
                n_samples)
            add(f"part 3 {scope} {response} subjects",
                f"SELECT COUNT(DISTINCT subject) FROM v_sample_annotated WHERE {where} AND response = '{response}'",
                n_subjects)

    part4 = ("condition = 'melanoma' AND sex = 'M' AND response = 'yes' "
             "AND time_from_treatment_start = 0")
    add("part 4 slice rows", f"SELECT COUNT(*) FROM v_sample_annotated WHERE {part4}", 485)
    add("part 4 slice subjects", f"SELECT COUNT(DISTINCT subject) FROM v_sample_annotated WHERE {part4}", 485)

    failures = []
    for label, actual, expected in checks:
        status = "ok" if actual == expected else "MISMATCH"
        print(f"  {label:38s} {actual:6d}  expected {expected:6d}  {status}")
        if actual != expected:
            failures.append(f"{label}: got {actual}, contract says {expected}")
    return failures


def main() -> int:
    OUTPUTS.mkdir(exist_ok=True)
    conn = connect()

    print("Cohort size checks against CONTRACT.md:")
    failures = check_cohort_sizes(conn)

    freq = get_frequencies(conn)
    freq.to_csv(OUTPUTS / "summary_frequencies.csv", index=False)
    print(f"\nPart 2: {len(freq)} frequency rows over {freq['sample'].nunique()} samples")

    stats_all = responder_comparison(conn, baseline_only=False)
    stats_base = responder_comparison(conn, baseline_only=True)
    stats_all.to_csv(OUTPUTS / "responder_stats_all_timepoints.csv", index=False)
    stats_base.to_csv(OUTPUTS / "responder_stats_baseline.csv", index=False)

    for label, frame, expected in (
        ("all timepoints", stats_all, (993, 975, 331, 325)),
        ("baseline", stats_base, (331, 325, 331, 325)),
    ):
        row = frame.iloc[0]
        actual = (int(row.n_yes), int(row.n_no), int(row.n_subj_yes), int(row.n_subj_no))
        status = "ok" if actual == expected else "MISMATCH"
        print(f"  part 3 frame {label:15s} {actual}  expected {expected}  {status}")
        if actual != expected:
            failures.append(f"part 3 frame {label}: got {actual}, contract says {expected}")

    print("\nPart 3, all timepoints:")
    print(stats_all.to_string(index=False))
    print("\nPart 3, baseline only:")
    print(stats_base.to_string(index=False))

    for label, frame in (("all timepoints", stats_all), ("baseline", stats_base)):
        for message in frame.attrs.get("sign_disagreements", []):
            print(f"  note ({label}): {message}")

    boxplot_responders(cohort_frame(conn), OUTPUTS / "boxplot_responders.png")

    base = baseline_cohort(conn)
    base.to_csv(OUTPUTS / "part4_baseline_cohort.csv", index=False)
    per_project = samples_per_project(conn)
    by_response = subjects_by_response(conn)
    by_sex = subjects_by_sex(conn)
    mean_bcell_count = mean_bcell_melanoma_male_responders_baseline(conn)
    mean_bcell_pct = mean_bcell_percentage_melanoma_male_responders_baseline(conn)

    summary = {
        "baseline_cohort_samples": int(len(base)),
        "baseline_cohort_subjects": int(base["subject"].nunique()),
        "samples_per_project": per_project.set_index("project")["n_samples"].astype(int).to_dict(),
        "subjects_by_response": by_response.set_index("response")["n_subjects"].astype(int).to_dict(),
        "subjects_by_sex": by_sex.set_index("sex")["n_subjects"].astype(int).to_dict(),
        "mean_bcell_count_melanoma_male_responders_baseline": mean_bcell_count,
        "mean_bcell_percentage_melanoma_male_responders_baseline": mean_bcell_pct,
    }
    with open(OUTPUTS / "part4_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print("\nPart 4:")
    for key, value in sorted(summary.items()):
        print(f"  {key}: {value}")

    conn.close()

    if failures:
        print("\nFAILED cohort checks:", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        return 1
    print("\nAll cohort checks passed. Outputs written to outputs/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Part 3. Responders vs non-responders, one test per cell population."""

from __future__ import annotations

import sqlite3
import warnings

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from .db import population_order

COLUMNS = [
    "population",
    "n_yes",
    "n_no",
    "n_subj_yes",
    "n_subj_no",
    "median_yes",
    "median_no",
    "median_diff",
    "cliffs_delta",
    "magnitude",
    "u_statistic",
    "p_value",
    "p_adj",
    "significant",
]

_COHORT_QUERY = """
SELECT sample, subject, response, population, percentage, time_from_treatment_start
FROM v_analysis
WHERE condition = 'melanoma'
  AND treatment = 'miraclib'
  AND sample_type = 'PBMC'
  AND response IN ('yes', 'no')
  {time_clause}
ORDER BY sample, population
"""


def cohort_frame(conn: sqlite3.Connection, baseline_only: bool = False) -> pd.DataFrame:
    """Long frame for the Part 3 cohort. Five rows per sample."""
    time_clause = "AND time_from_treatment_start = 0" if baseline_only else ""
    df = pd.read_sql_query(_COHORT_QUERY.format(time_clause=time_clause), conn)
    if df.empty:
        raise ValueError("Part 3 cohort query returned no rows")
    return df


def magnitude_label(delta: float) -> str:
    """Romano thresholds on |Cliff's delta|."""
    d = abs(delta)
    if d < 0.147:
        return "negligible"
    if d < 0.33:
        return "small"
    if d < 0.474:
        return "medium"
    return "large"


def compare_frame(df: pd.DataFrame, populations: list[str]) -> pd.DataFrame:
    """Run the comparison over an already loaded cohort frame.

    Split out from responder_comparison so the figure can annotate its panels
    without a second database round trip, and without a second copy of the math.
    """
    rows: list[dict] = []
    disagreements: list[str] = []
    for population in populations:
        sub = df[df["population"] == population]
        if sub.empty:
            raise ValueError(f"no cohort rows for population {population}")
        yes = sub[sub["response"] == "yes"]
        no = sub[sub["response"] == "no"]
        x = yes["percentage"].to_numpy(dtype=float)
        y = no["percentage"].to_numpy(dtype=float)
        if x.size == 0 or y.size == 0:
            raise ValueError(f"empty arm for population {population}")

        u = mannwhitneyu(x=x, y=y, alternative="two-sided")
        # U is computed for x, so 2U/(nx*ny) - 1 is positive when responders sit higher.
        delta = 2.0 * u.statistic / (x.size * y.size) - 1.0

        median_yes = float(np.median(x))
        median_no = float(np.median(y))
        median_diff = median_yes - median_no
        # The median difference compares one point of each distribution. Cliff's delta
        # compares every pair. When the two distributions cross, they can point opposite
        # ways and both still be right, so a mismatch is flagged rather than raised.
        # Verified by brute force pair counting: the delta here matches
        # (#x>y - #x<y) / (nx*ny) exactly for every population.
        if median_diff != 0.0 and delta != 0.0 and np.sign(delta) != np.sign(median_diff):
            disagreements.append(
                f"{population}: cliffs_delta {delta:+.6f} and median_diff "
                f"{median_diff:+.6f} disagree in sign (p_value {u.pvalue:.4f})"
            )

        rows.append(
            {
                "population": population,
                "n_yes": int(yes["sample"].nunique()),
                "n_no": int(no["sample"].nunique()),
                "n_subj_yes": int(yes["subject"].nunique()),
                "n_subj_no": int(no["subject"].nunique()),
                "median_yes": median_yes,
                "median_no": median_no,
                "median_diff": median_diff,
                "cliffs_delta": float(delta),
                "magnitude": magnitude_label(delta),
                "u_statistic": float(u.statistic),
                "p_value": float(u.pvalue),
            }
        )

    out = pd.DataFrame(rows)
    _, p_adj, _, _ = multipletests(out["p_value"].to_numpy(), method="fdr_bh")
    out["p_adj"] = p_adj
    out["significant"] = out["p_adj"] < 0.05
    # Carried on attrs so the CSV keeps exactly the contracted columns.
    out.attrs["sign_disagreements"] = disagreements
    for message in disagreements:
        warnings.warn(f"sign check: {message}", stacklevel=2)
    return out[COLUMNS]


def responder_comparison(
    conn: sqlite3.Connection, baseline_only: bool = False
) -> pd.DataFrame:
    return compare_frame(
        cohort_frame(conn, baseline_only=baseline_only), population_order(conn)
    )

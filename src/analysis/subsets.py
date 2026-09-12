"""Part 4. Baseline cohort and the breakdowns asked for on top of it."""

from __future__ import annotations

import sqlite3

import pandas as pd

_BASELINE_WHERE = """
WHERE condition = 'melanoma'
  AND treatment = 'miraclib'
  AND sample_type = 'PBMC'
  AND time_from_treatment_start = 0
"""


def baseline_cohort(conn: sqlite3.Connection) -> pd.DataFrame:
    """Melanoma, miraclib, PBMC, time 0. One row per sample."""
    df = pd.read_sql_query(
        "SELECT sample, subject, project, condition, treatment, sample_type, "
        "time_from_treatment_start, response, sex, age "
        f"FROM v_sample_annotated {_BASELINE_WHERE} ORDER BY sample",
        conn,
    )
    if df.empty:
        raise ValueError("baseline cohort query returned no rows")
    return df


def samples_per_project(conn: sqlite3.Connection) -> pd.DataFrame:
    """Counts samples, one per row of v_sample_annotated.

    Projects with no baseline sample are kept at 0. An omitted row looks like a
    filter nobody noticed. prj2 lands here because its melanoma miraclib subjects
    are all WB, never PBMC.
    """
    return pd.read_sql_query(
        "SELECT p.project_code AS project, COUNT(v.sample) AS n_samples "
        "FROM projects p "
        f"LEFT JOIN (SELECT sample, project FROM v_sample_annotated {_BASELINE_WHERE}) v "
        "  ON v.project = p.project_code "
        "GROUP BY p.project_code ORDER BY p.project_code",
        conn,
    )


def subjects_by_response(conn: sqlite3.Connection) -> pd.DataFrame:
    return pd.read_sql_query(
        "SELECT response, COUNT(DISTINCT subject) AS n_subjects "
        f"FROM v_sample_annotated {_BASELINE_WHERE} GROUP BY response ORDER BY response",
        conn,
    )


def subjects_by_sex(conn: sqlite3.Connection) -> pd.DataFrame:
    return pd.read_sql_query(
        "SELECT sex, COUNT(DISTINCT subject) AS n_subjects "
        f"FROM v_sample_annotated {_BASELINE_WHERE} GROUP BY sex ORDER BY sex",
        conn,
    )


_BCELL_SLICE = """
SELECT count, percentage
FROM v_analysis
WHERE condition = 'melanoma'
  AND sex = 'M'
  AND response = 'yes'
  AND time_from_treatment_start = 0
  AND population = 'b_cell'
"""


def _bcell_slice(conn: sqlite3.Connection) -> pd.DataFrame:
    """The 485 row slice behind both Part 4 b_cell means.

    The question says all sample and treatment types, so neither is filtered here.
    """
    df = pd.read_sql_query(_BCELL_SLICE, conn)
    if df.empty:
        raise ValueError("b_cell slice returned no rows")
    return df


def mean_bcell_melanoma_male_responders_baseline(conn: sqlite3.Connection) -> float:
    """Mean b_cell count for melanoma males who responded, at time 0.

    "Average number of B cells" is a count, not a frequency. Part 3 is the part
    that works off relative frequencies. Part 4 filters the database directly, so
    it reads the raw counts.
    """
    return round(float(_bcell_slice(conn)["count"].mean()), 2)


def mean_bcell_percentage_melanoma_male_responders_baseline(
    conn: sqlite3.Connection,
) -> float:
    """Same slice as above, reported as a relative frequency instead of a count."""
    return round(float(_bcell_slice(conn)["percentage"].mean()), 2)

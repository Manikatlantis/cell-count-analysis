#!/usr/bin/env python3
"""Build cell_count.db from cell-count.csv.

Run as: python load_data.py

No arguments. The grader invokes this file directly from the repo root, so
every path below is resolved relative to this file and not to the caller's
working directory.

The load is a full rebuild. sql/schema.sql drops every object it owns before
recreating it, and the database file itself is removed first, so a second run
produces the same database as the first. That is what makes the pipeline
idempotent without needing upserts or delete-then-insert bookkeeping.
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
SCHEMA_PATH = ROOT / "sql" / "schema.sql"
DB_PATH = ROOT / "cell_count.db"

# Verbatim from the CSV header. Order matters: it fixes the population ids.
POPULATION_COLUMNS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

METADATA_COLUMNS = [
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
]

EXPECTED_COLUMNS = METADATA_COLUMNS + POPULATION_COLUMNS

# Columns the schema stores once per subject rather than once per sample.
SUBJECT_LEVEL_COLUMNS = ["project", "condition", "age", "sex"]


class LoadError(RuntimeError):
    """Raised when the CSV does not match what the schema can represent."""


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise LoadError(f"missing input file: {path}")
    df = pd.read_csv(path)
    if list(df.columns) != EXPECTED_COLUMNS:
        raise LoadError(
            "unexpected CSV header.\n"
            f"  expected: {EXPECTED_COLUMNS}\n"
            f"  found:    {list(df.columns)}"
        )
    return df


def validate(df: pd.DataFrame) -> None:
    """Fail loudly on anything the schema would silently mangle.

    A CHECK violation halfway through an insert leaves a confusing error, so
    the same conditions are tested here where the message can name the column.
    """
    if df.empty:
        raise LoadError("CSV has a header but no rows")

    dup_samples = df["sample"][df["sample"].duplicated()].unique()
    if len(dup_samples):
        raise LoadError(f"duplicate sample codes: {list(dup_samples[:10])}")

    for col in METADATA_COLUMNS:
        if col == "response":
            continue  # NULL is meaningful here: untreated or not yet assessed
        if df[col].isna().any():
            raise LoadError(f"nulls in required column {col!r}")

    for col in POPULATION_COLUMNS:
        if df[col].isna().any():
            raise LoadError(f"nulls in count column {col!r}")
        if (df[col] < 0).any():
            raise LoadError(f"negative counts in {col!r}")

    # subjects holds one row per subject, so these must not vary within one.
    for col in SUBJECT_LEVEL_COLUMNS:
        varying = df.groupby("subject")[col].nunique(dropna=False)
        bad = varying[varying > 1]
        if len(bad):
            raise LoadError(
                f"{col!r} varies within subject, cannot store it per subject: "
                f"{list(bad.index[:10])}"
            )

    # treatment_courses is keyed on (subject, treatment) and carries one
    # response, so a subject on one drug cannot have two different outcomes.
    resp = df.groupby(["subject", "treatment"])["response"].nunique(dropna=False)
    bad = resp[resp > 1]
    if len(bad):
        raise LoadError(f"conflicting response values for: {list(bad.index[:10])}")

    bad_resp = set(df["response"].dropna().unique()) - {"yes", "no"}
    if bad_resp:
        raise LoadError(f"unexpected response values: {sorted(bad_resp)}")

    bad_sex = set(df["sex"].dropna().unique()) - {"M", "F"}
    if bad_sex:
        raise LoadError(f"unexpected sex values: {sorted(bad_sex)}")

    # The schema declares UNIQUE (course_id, sample_type, time).
    key = ["subject", "treatment", "sample_type", "time_from_treatment_start"]
    dup_key = df[df.duplicated(key, keep=False)]
    if len(dup_key):
        raise LoadError(
            f"{len(dup_key)} rows collide on {key}, first: "
            f"{dup_key.iloc[0][key].to_dict()}"
        )


def _lookup_rows(values) -> list[tuple[int, str]]:
    """Number distinct values from 1 in sorted order.

    Sorting keeps the surrogate ids stable across runs, which is what lets a
    rebuilt database compare equal to the previous one.
    """
    return [(i, v) for i, v in enumerate(sorted(values), start=1)]


def _as_null(value):
    """pandas gives NaN for a blank cell. SQLite needs None to store NULL."""
    return None if pd.isna(value) else value


def load(conn: sqlite3.Connection, df: pd.DataFrame) -> None:
    projects = _lookup_rows(df["project"].unique())
    conditions = _lookup_rows(df["condition"].unique())
    treatments = _lookup_rows(df["treatment"].unique())
    populations = list(enumerate(POPULATION_COLUMNS, start=1))

    conn.executemany("INSERT INTO projects VALUES (?, ?)", projects)
    conn.executemany("INSERT INTO conditions VALUES (?, ?)", conditions)
    conn.executemany("INSERT INTO treatments VALUES (?, ?)", treatments)
    conn.executemany("INSERT INTO cell_populations VALUES (?, ?)", populations)

    project_id = {code: i for i, code in projects}
    condition_id = {name: i for i, name in conditions}
    treatment_id = {name: i for i, name in treatments}

    subject_meta = (
        df.sort_values("subject")
        .drop_duplicates("subject")
        .set_index("subject")[SUBJECT_LEVEL_COLUMNS]
    )
    subject_id = {code: i for i, code in enumerate(subject_meta.index, start=1)}
    conn.executemany(
        "INSERT INTO subjects VALUES (?, ?, ?, ?, ?, ?)",
        [
            (
                subject_id[code],
                code,
                project_id[row.project],
                condition_id[row.condition],
                int(row.age),
                row.sex,
            )
            for code, row in subject_meta.iterrows()
        ],
    )

    courses = (
        df[["subject", "treatment", "response"]]
        .drop_duplicates(["subject", "treatment"])
        .sort_values(["subject", "treatment"])
    )
    course_id = {}
    course_rows = []
    for i, row in enumerate(courses.itertuples(index=False), start=1):
        course_id[(row.subject, row.treatment)] = i
        course_rows.append(
            (i, subject_id[row.subject], treatment_id[row.treatment], _as_null(row.response))
        )
    conn.executemany("INSERT INTO treatment_courses VALUES (?, ?, ?, ?)", course_rows)

    samples = df.sort_values("sample")
    sample_rows = []
    count_rows = []
    for i, row in enumerate(samples.itertuples(index=False), start=1):
        sample_rows.append(
            (
                i,
                row.sample,
                course_id[(row.subject, row.treatment)],
                row.sample_type,
                int(row.time_from_treatment_start),
            )
        )
        for population_id, column in populations:
            count_rows.append((i, population_id, int(getattr(row, column))))

    conn.executemany("INSERT INTO samples VALUES (?, ?, ?, ?, ?)", sample_rows)
    conn.executemany("INSERT INTO sample_counts VALUES (?, ?, ?)", count_rows)


def report(conn: sqlite3.Connection) -> None:
    tables = [
        "projects",
        "conditions",
        "treatments",
        "cell_populations",
        "subjects",
        "treatment_courses",
        "samples",
        "sample_counts",
    ]
    print("\nrow counts")
    for table in tables:
        n = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"  {table:<20} {n:>7}")

    print("\nviews")
    for view in ["v_sample_frequencies", "v_sample_annotated", "v_analysis"]:
        n = conn.execute(f"SELECT COUNT(*) FROM {view}").fetchone()[0]
        print(f"  {view:<20} {n:>7}")

    null_courses = conn.execute(
        "SELECT COUNT(*) FROM treatment_courses WHERE response IS NULL"
    ).fetchone()[0]
    print(f"\ncourses with NULL response: {null_courses}")


def main() -> int:
    print(f"reading  {CSV_PATH}")
    df = read_csv(CSV_PATH)
    print(f"  {len(df)} rows, {len(df.columns)} columns")

    validate(df)
    print("  validation passed")

    # Removing the file is what guarantees a rerun cannot accumulate rows.
    if DB_PATH.exists():
        print(f"removing {DB_PATH} (rebuilding from scratch)")
        DB_PATH.unlink()

    print(f"applying {SCHEMA_PATH}")
    schema_sql = SCHEMA_PATH.read_text()

    conn = sqlite3.connect(DB_PATH)
    try:
        # executescript commits first, which would drop a pragma set inside it.
        conn.executescript(schema_sql)
        conn.execute("PRAGMA foreign_keys = ON")
        with conn:
            load(conn, df)
        violations = conn.execute("PRAGMA foreign_key_check").fetchall()
        if violations:
            raise LoadError(f"foreign key violations: {violations[:5]}")
        report(conn)
    finally:
        conn.close()

    print(f"\nwrote {DB_PATH}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except LoadError as exc:
        print(f"load failed: {exc}", file=sys.stderr)
        sys.exit(1)

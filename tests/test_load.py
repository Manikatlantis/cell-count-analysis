"""Loader invariants from CONTRACT.md.

Every number asserted here is copied from the "Loader invariants" and "Known
cohort sizes" sections of CONTRACT.md. If one of these fails, either the
loader is wrong or the contract needs renegotiating. Do not edit a number
here to make a test pass.
"""

from __future__ import annotations

import hashlib
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "cell_count.db"
LOADER = ROOT / "load_data.py"

TABLES = [
    "projects",
    "conditions",
    "treatments",
    "cell_populations",
    "subjects",
    "treatment_courses",
    "samples",
    "sample_counts",
]


def run_loader() -> subprocess.CompletedProcess:
    """Invoke the loader the way the grader does, from the repo root."""
    proc = subprocess.run(
        [sys.executable, "load_data.py"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, f"load_data.py failed:\n{proc.stdout}\n{proc.stderr}"
    return proc


@pytest.fixture(scope="session", autouse=True)
def database():
    if not DB_PATH.exists():
        run_loader()
    return DB_PATH


@pytest.fixture
def conn():
    c = sqlite3.connect(DB_PATH)
    try:
        yield c
    finally:
        c.close()


def scalar(conn: sqlite3.Connection, sql: str, params=()) -> int:
    return conn.execute(sql, params).fetchone()[0]


def samples_in(conn: sqlite3.Connection, where: str, params=()) -> int:
    """Distinct sample count from v_analysis.

    v_analysis is long format, five population rows per sample, so a bare
    COUNT(*) there is five times the CSV row count.
    """
    return scalar(conn, f"SELECT COUNT(DISTINCT sample) FROM v_analysis WHERE {where}", params)


def subjects_in(conn: sqlite3.Connection, where: str, params=()) -> int:
    return scalar(conn, f"SELECT COUNT(DISTINCT subject) FROM v_analysis WHERE {where}", params)


MELANOMA = "condition = 'melanoma'"
MIRACLIB = MELANOMA + " AND treatment = 'miraclib'"
PBMC = MIRACLIB + " AND sample_type = 'PBMC'"
BASELINE = PBMC + " AND time_from_treatment_start = 0"


# --- Loader invariants, CONTRACT.md -----------------------------------------


def test_database_file_exists():
    assert DB_PATH.exists(), "load_data.py must create cell_count.db in the repo root"


def test_sample_counts_is_ten_thousand_five_hundred_times_five(conn):
    assert scalar(conn, "SELECT COUNT(*) FROM sample_counts") == 10500 * 5 == 52500


def test_core_table_counts(conn):
    assert scalar(conn, "SELECT COUNT(*) FROM samples") == 10500
    assert scalar(conn, "SELECT COUNT(*) FROM subjects") == 3500
    assert scalar(conn, "SELECT COUNT(*) FROM treatment_courses") == 3500


def test_five_populations_named_as_in_the_csv(conn):
    names = [r[0] for r in conn.execute("SELECT population_name FROM cell_populations ORDER BY population_id")]
    assert names == ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def test_response_is_sql_null_not_a_string(conn):
    values = {r[0] for r in conn.execute("SELECT DISTINCT response FROM treatment_courses")}
    assert values == {"yes", "no", None}

    # The failure mode this guards is a NaN reaching SQLite as text.
    strays = scalar(
        conn,
        "SELECT COUNT(*) FROM treatment_courses "
        "WHERE response IN ('nan', 'NaN', 'None', 'NULL', '')",
    )
    assert strays == 0

    assert scalar(conn, "SELECT COUNT(*) FROM treatment_courses WHERE response IS NULL") == 474
    assert samples_in(conn, "response IS NULL") == 1422


def test_untreated_subjects_are_the_ones_with_null_response(conn):
    # CONTRACT.md: 1422 rows, all healthy and treatment none.
    assert samples_in(conn, "response IS NULL AND condition = 'healthy' AND treatment = 'none'") == 1422
    assert samples_in(conn, "response IS NOT NULL AND treatment = 'none'") == 0


def test_every_subject_has_exactly_three_samples(conn):
    lo, hi, n = conn.execute(
        "SELECT MIN(n), MAX(n), COUNT(*) FROM ("
        "  SELECT subject, COUNT(DISTINCT sample) AS n FROM v_analysis GROUP BY subject"
        ")"
    ).fetchone()
    assert (lo, hi, n) == (3, 3, 3500)


def test_no_duplicate_sample_codes(conn):
    assert scalar(conn, "SELECT COUNT(DISTINCT sample_code) FROM samples") == 10500


def test_loading_twice_changes_nothing():
    def snapshot() -> dict[str, tuple[int, str]]:
        c = sqlite3.connect(DB_PATH)
        try:
            out = {}
            for table in TABLES:
                rows = c.execute(f"SELECT * FROM {table} ORDER BY 1, 2").fetchall()
                digest = hashlib.sha256(repr(rows).encode()).hexdigest()
                out[table] = (len(rows), digest)
            return out
        finally:
            c.close()

    run_loader()
    first = snapshot()
    run_loader()
    second = snapshot()

    assert first == second, "a second load_data.py run changed the database"
    assert first["sample_counts"][0] == 52500
    assert first["samples"][0] == 10500


# --- Known cohort sizes, CONTRACT.md, queried through v_analysis ------------


def test_v_analysis_is_five_rows_per_sample(conn):
    assert scalar(conn, "SELECT COUNT(*) FROM v_analysis") == 52500
    assert scalar(conn, "SELECT COUNT(DISTINCT sample) FROM v_analysis") == 10500


def test_cohort_filter_chain(conn):
    assert samples_in(conn, MELANOMA) == 5175
    assert samples_in(conn, MIRACLIB) == 2655
    assert samples_in(conn, PBMC) == 1968
    assert samples_in(conn, BASELINE) == 656

    assert subjects_in(conn, PBMC) == 656
    assert subjects_in(conn, BASELINE) == 656


def test_part3_arms_all_timepoints(conn):
    assert samples_in(conn, PBMC + " AND response = 'yes'") == 993
    assert samples_in(conn, PBMC + " AND response = 'no'") == 975
    assert subjects_in(conn, PBMC + " AND response = 'yes'") == 331
    assert subjects_in(conn, PBMC + " AND response = 'no'") == 325


def test_part3_arms_baseline_only(conn):
    assert samples_in(conn, BASELINE + " AND response = 'yes'") == 331
    assert samples_in(conn, BASELINE + " AND response = 'no'") == 325
    assert subjects_in(conn, BASELINE + " AND response = 'yes'") == 331
    assert subjects_in(conn, BASELINE + " AND response = 'no'") == 325


def test_part4_final_slice_has_no_sample_type_or_treatment_filter(conn):
    where = (
        "condition = 'melanoma' AND sex = 'M' AND response = 'yes' "
        "AND time_from_treatment_start = 0"
    )
    assert samples_in(conn, where) == 485
    assert subjects_in(conn, where) == 485


# --- Frequency view sanity ---------------------------------------------------


def test_percentages_sum_to_one_hundred_per_sample(conn):
    off = scalar(
        conn,
        "SELECT COUNT(*) FROM ("
        "  SELECT sample, SUM(percentage) AS pct FROM v_sample_frequencies GROUP BY sample"
        ") WHERE ABS(pct - 100.0) > 0.01",
    )
    assert off == 0


def test_total_count_matches_the_sum_of_the_five_populations(conn):
    mismatched = scalar(
        conn,
        "SELECT COUNT(*) FROM ("
        "  SELECT sample, total_count, SUM(count) AS s FROM v_sample_frequencies"
        "  GROUP BY sample, total_count"
        ") WHERE s != total_count",
    )
    assert mismatched == 0

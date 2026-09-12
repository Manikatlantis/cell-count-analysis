"""Streamlit dashboard over cell_count.db.

Parts 2, 3 and 4 of the assignment, one tab each. The database is the only
source of data here. The CSV is never read.
"""

import sys
from pathlib import Path

# Streamlit runs this file as a script from inside app/, so the repo root is not
# on sys.path and "import src.analysis..." would fail. Put the root in front.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import sqlite3
from contextlib import closing

import matplotlib

# No display in a Codespace, so a GUI backend would fail at import time.
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

DB_PATH = REPO_ROOT / "cell_count.db"
OUTPUTS_DIR = REPO_ROOT / "outputs"

# Populations are fixed by the schema. Keeping the order stable keeps the
# boxplot panels in the same place between reruns.
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
RESPONSE_ORDER = ["yes", "no"]
RESPONSE_COLORS = {"yes": "#4c78a8", "no": "#f58518"}

# Part 3 cohort, locked: melanoma, miraclib, PBMC.
COHORT_SQL = """
SELECT sample, subject, response, population, percentage, count, total_count,
       time_from_treatment_start
FROM v_analysis
WHERE condition = 'melanoma'
  AND treatment = 'miraclib'
  AND sample_type = 'PBMC'
  AND response IS NOT NULL
"""

# The analysis package may not be importable in a half-built checkout. Each
# module is tried on its own so one missing file only disables its own tab, and
# the page shows a message instead of a stack trace.
import importlib  # noqa: E402

ANALYSIS = {}
IMPORT_ERRORS = {}
for _name in ("frequency", "stats", "subsets", "plots"):
    try:
        ANALYSIS[_name] = importlib.import_module(f"src.analysis.{_name}")
    except Exception as exc:  # noqa: BLE001
        IMPORT_ERRORS[_name] = exc


st.set_page_config(page_title="Cell count dashboard", layout="wide")


def connect():
    """Fresh connection per read. Streamlit reruns across threads."""
    return closing(sqlite3.connect(DB_PATH, check_same_thread=False))


def analysis_ready(part: str, module: str) -> bool:
    exc = IMPORT_ERRORS.get(module)
    if exc is None:
        return True
    st.error(
        f"{part} needs src/analysis/{module}.py, which did not import: "
        f"{type(exc).__name__}: {exc}. Run: make pipeline"
    )
    return False


def read_output(name: str, required: set) -> pd.DataFrame | None:
    """Load a precomputed CSV from outputs/, or None if it is unusable.

    outputs/ is written by run_analysis.py and may be absent, partial, or stale.
    Any problem means the caller recomputes from the database instead.
    """
    path = OUTPUTS_DIR / name
    if not path.exists():
        return None
    try:
        df = pd.read_csv(path)
    except Exception:  # noqa: BLE001
        return None
    if not required.issubset(df.columns):
        return None
    return df


@st.cache_data(show_spinner=False)
def load_frequencies() -> tuple[pd.DataFrame, str]:
    cols = {"sample", "total_count", "population", "count", "percentage"}
    cached = read_output("summary_frequencies.csv", cols)
    if cached is not None:
        return cached, "outputs/summary_frequencies.csv"
    with connect() as conn:
        return ANALYSIS["frequency"].get_frequencies(conn), "database"


@st.cache_data(show_spinner=False)
def load_responder_stats(baseline_only: bool) -> tuple[pd.DataFrame, str]:
    name = (
        "responder_stats_baseline.csv"
        if baseline_only
        else "responder_stats_all_timepoints.csv"
    )
    cached = read_output(name, {"population", "p_value", "p_adj"})
    if cached is not None:
        return cached, f"outputs/{name}"
    with connect() as conn:
        fn = ANALYSIS["stats"].responder_comparison
        return fn(conn, baseline_only=baseline_only), "database"


@st.cache_data(show_spinner=False)
def load_cohort(baseline_only: bool) -> pd.DataFrame:
    sql = COHORT_SQL
    if baseline_only:
        sql += "  AND time_from_treatment_start = 0\n"
    with connect() as conn:
        return pd.read_sql_query(sql, conn)


@st.cache_data(show_spinner=False)
def load_sample_index() -> pd.DataFrame:
    sql = """
    SELECT sample, subject, project, condition, treatment, response,
           sample_type, time_from_treatment_start, sex, age
    FROM v_sample_annotated
    ORDER BY sample
    """
    with connect() as conn:
        return pd.read_sql_query(sql, conn)


@st.cache_data(show_spinner=False)
def load_part4() -> dict:
    mod = ANALYSIS["subsets"]
    with connect() as conn:
        return {
            "baseline_cohort": mod.baseline_cohort(conn),
            "samples_per_project": mod.samples_per_project(conn),
            "subjects_by_response": mod.subjects_by_response(conn),
            "subjects_by_sex": mod.subjects_by_sex(conn),
            "mean_bcell": mod.mean_bcell_melanoma_male_responders_baseline(conn),
            "mean_bcell_pct": (
                mod.mean_bcell_percentage_melanoma_male_responders_baseline(conn)
            ),
        }


def format_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Order the contract columns first and make the p columns readable."""
    order = [
        "population",
        "n_yes",
        "n_no",
        "n_subj_yes",
        "n_subj_no",
        "median_yes",
        "median_no",
        "median_diff",
        "cliffs_delta",
        "u_statistic",
        "p_value",
        "p_adj",
        "significant",
    ]
    known = [c for c in order if c in df.columns]
    rest = [c for c in df.columns if c not in known]
    out = df[known + rest].copy()
    for col in ("median_yes", "median_no", "median_diff", "cliffs_delta"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce").round(3)
    for col in ("p_value", "p_adj"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce").map(
                lambda v: "" if pd.isna(v) else f"{v:.3g}"
            )
    return out


def two_column_chart(df: pd.DataFrame):
    """Bar chart for a small two column summary. Skipped if the shape differs."""
    if df.shape[1] != 2:
        return
    label, value = df.columns[0], df.columns[1]
    if not pd.api.types.is_numeric_dtype(df[value]):
        return
    st.bar_chart(df.set_index(label)[value])


if not DB_PATH.exists():
    st.title("Cell count dashboard")
    st.error(f"{DB_PATH.name} not found in the repo root. Run: make pipeline")
    st.stop()

st.title("Cell count dashboard")
st.caption(
    "All figures are read from cell_count.db. Precomputed CSVs under outputs/ "
    "are used when present, otherwise the database is queried live."
)

if IMPORT_ERRORS:
    st.warning(
        "These analysis modules did not import, so their tabs are unavailable: "
        + ", ".join(sorted(IMPORT_ERRORS))
    )

tab2, tab3, tab4 = st.tabs(
    ["Part 2 relative frequencies", "Part 3 responders", "Part 4 subsets"]
)


# Each tab body is a function so a failure inside the analysis layer shows up as
# a message on the page instead of a stack trace.
def render_part2():
    freqs, source = load_frequencies()
    index = load_sample_index()

    left, mid, right = st.columns(3)
    n_samples = freqs["sample"].nunique()
    left.metric("Samples", f"{n_samples:,}")
    mid.metric("Populations", f"{freqs['population'].nunique():,}")
    right.metric("Rows", f"{len(freqs):,}")
    st.caption(f"Source: {source}")

    f1, f2, f3 = st.columns([2, 2, 2])
    text = f1.text_input("Sample code contains", "")
    pops = f2.multiselect("Populations", POPULATIONS, default=POPULATIONS)
    joinable = index.set_index("sample")
    conditions = ["(any)"] + sorted(joinable["condition"].dropna().unique())
    cond = f3.selectbox("Condition", conditions)

    view = freqs
    if text.strip():
        view = view[view["sample"].str.contains(text.strip(), case=False)]
    if pops:
        view = view[view["population"].isin(pops)]
    if cond != "(any)":
        keep = set(joinable.index[joinable["condition"] == cond])
        view = view[view["sample"].isin(keep)]

    st.write(
        f"{view['sample'].nunique():,} samples, {len(view):,} rows after filtering"
    )
    if view.empty:
        st.warning("No rows match those filters.")
    else:
        st.dataframe(view.head(2000), width="stretch", hide_index=True)
        if len(view) > 2000:
            st.caption("Showing the first 2000 rows. Download for the full table.")
        st.download_button(
            "Download filtered table as CSV",
            view.to_csv(index=False).encode(),
            file_name="relative_frequencies.csv",
            mime="text/csv",
        )

    picked = st.selectbox(
        "Inspect one sample",
        ["(none)"] + sorted(view["sample"].unique().tolist())[:5000],
    )
    if picked != "(none)":
        st.dataframe(
            view[view["sample"] == picked], width="stretch", hide_index=True
        )
        if picked in joinable.index:
            st.write(joinable.loc[[picked]])


def render_part3():
    scope = st.radio(
        "Timepoints",
        ["All timepoints", "Baseline only (time 0)"],
        horizontal=True,
    )
    baseline_only = scope.startswith("Baseline")

    table, source = load_responder_stats(baseline_only)
    cohort = load_cohort(baseline_only)

    samples = cohort.groupby("response")["sample"].nunique()
    subjects = cohort.groupby("response")["subject"].nunique()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Responder samples", int(samples.get("yes", 0)))
    c2.metric("Non responder samples", int(samples.get("no", 0)))
    c3.metric("Responder subjects", int(subjects.get("yes", 0)))
    c4.metric("Non responder subjects", int(subjects.get("no", 0)))

    if table is None or table.empty:
        st.error("The comparison returned no rows.")
    else:
        st.dataframe(format_stats(table), width="stretch", hide_index=True)
        st.caption(f"Source: {source}")
        st.download_button(
            "Download stats table as CSV",
            table.to_csv(index=False).encode(),
            file_name=(
                "responder_stats_baseline.csv"
                if baseline_only
                else "responder_stats_all_timepoints.csv"
            ),
            mime="text/csv",
        )

    if cohort.empty:
        st.error("The cohort query returned no rows.")
    elif analysis_ready("The figure", "plots"):
        # The same builder run_analysis.py saves to outputs/, so the committed
        # PNG and this page cannot drift. It draws both runs, so it is fed the
        # unfiltered cohort regardless of the timepoint toggle above.
        fig = ANALYSIS["plots"]._build_boxplot(load_cohort(False))
        st.pyplot(fig)
        plt.close(fig)
        st.caption(
            "Each point is one sample. Both runs are shown, so this figure does "
            "not follow the timepoint selector. Outlier markers are off because "
            "every point is already drawn."
        )


def render_part4():
    part4 = load_part4()

    baseline = part4["baseline_cohort"]
    st.markdown("Baseline cohort: melanoma, miraclib, PBMC, time 0")
    b1, b2, b3 = st.columns(3)
    b1.metric("Rows", f"{len(baseline):,}")
    if "sample" in baseline.columns:
        b2.metric("Samples", f"{baseline['sample'].nunique():,}")
    if "subject" in baseline.columns:
        b3.metric("Subjects", f"{baseline['subject'].nunique():,}")
    if baseline.empty:
        st.error("The baseline cohort is empty.")
    else:
        st.dataframe(baseline.head(1000), width="stretch", hide_index=True)
        if len(baseline) > 1000:
            st.caption("Showing the first 1000 rows.")
        st.download_button(
            "Download baseline cohort as CSV",
            baseline.to_csv(index=False).encode(),
            file_name="part4_baseline_cohort.csv",
            mime="text/csv",
        )

    st.divider()
    col_a, col_b, col_c = st.columns(3)
    with col_a:
        st.markdown("Samples per project")
        # A project with no baseline sample must read 0. A blank cell would look
        # like a filter nobody noticed.
        per_project = part4["samples_per_project"].copy()
        if per_project.shape[1] == 2:
            counts = per_project.columns[1]
            per_project[counts] = (
                pd.to_numeric(per_project[counts], errors="coerce")
                .fillna(0)
                .astype(int)
            )
        st.dataframe(per_project, width="stretch", hide_index=True)
        two_column_chart(per_project)
    with col_b:
        st.markdown("Subjects by response")
        st.dataframe(
            part4["subjects_by_response"], width="stretch", hide_index=True
        )
        two_column_chart(part4["subjects_by_response"])
    with col_c:
        st.markdown("Subjects by sex")
        st.dataframe(part4["subjects_by_sex"], width="stretch", hide_index=True)
        two_column_chart(part4["subjects_by_sex"])
    st.caption(
        "Samples per project counts sample rows. The subject tables count "
        "distinct subjects."
    )

    st.divider()
    st.markdown("Melanoma males, responders, time 0, all sample and treatment types")
    m1, m2 = st.columns(2)
    m1.metric("Mean b_cell count", f"{part4['mean_bcell']:,.2f}")
    m2.metric("Mean b_cell relative frequency (%)", f"{part4['mean_bcell_pct']:.2f}")
    st.caption(
        "The wording \"average number of B cells\" is ambiguous. The count is "
        "the answer here. Over the same 485 row slice the mean relative "
        f"frequency is {part4['mean_bcell_pct']:.2f} percent."
    )

with tab2:
    st.subheader("Relative frequency per sample and population")
    if analysis_ready("Part 2", "frequency"):
        try:
            render_part2()
        except Exception as exc:  # noqa: BLE001
            st.error(
                "Part 2 could not be rendered: "
                f"{type(exc).__name__}: {exc}"
            )

with tab3:
    st.subheader("Responders versus non responders")
    st.caption(
        "Cohort: melanoma, miraclib, PBMC. Test: two sided Mann Whitney U per "
        "population. Correction: Benjamini Hochberg across the five populations."
    )
    if analysis_ready("Part 3", "stats"):
        try:
            render_part3()
        except Exception as exc:  # noqa: BLE001
            st.error(
                "Part 3 could not be rendered: "
                f"{type(exc).__name__}: {exc}"
            )

with tab4:
    st.subheader("Data subsets")
    if analysis_ready("Part 4", "subsets"):
        try:
            render_part4()
        except Exception as exc:  # noqa: BLE001
            st.error(
                "Part 4 could not be rendered: "
                f"{type(exc).__name__}: {exc}"
            )

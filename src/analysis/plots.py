"""Part 3 figure. Two rows of five panels, responders against non-responders.

One implementation of this figure exists. run_analysis.py saves it to a PNG and
the dashboard renders the same Figure object, so the committed image and the
screen never drift apart.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # No display in Codespaces or CI.

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .stats import compare_frame  # noqa: E402

# CSV column order, which is the order the loader numbers the populations.
# Hardcoded because the figure builder takes a frame, not a connection.
POPULATION_ORDER = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

ARMS = ["yes", "no"]
ARM_LABELS = {"yes": "responder", "no": "non-responder"}
ARM_COLORS = {"yes": "#2b7bba", "no": "#d1704a"}
JITTER_SEED = 0  # Fixed so reruns produce a byte identical PNG.

REQUIRED = {"population", "response", "percentage", "time_from_treatment_start"}


def _panel(ax, sub: pd.DataFrame, stats_row, rng, show_ylabel: bool) -> None:
    arms = [sub.loc[sub["response"] == a, "percentage"].to_numpy(float) for a in ARMS]

    # Points first, boxes after, so the summary reads on top of the data.
    for pos, arm, values in zip([1, 2], ARMS, arms):
        x = pos + rng.uniform(-0.18, 0.18, size=values.size)
        ax.scatter(
            x, values, s=6, alpha=0.25, linewidths=0, color=ARM_COLORS[arm], zorder=2
        )
    ax.boxplot(
        arms,
        positions=[1, 2],
        widths=0.55,
        showfliers=False,
        patch_artist=False,  # Unfilled, so the points stay visible through the box.
        medianprops={"color": "black", "linewidth": 1.7},
        boxprops={"color": "black", "linewidth": 1.2},
        whiskerprops={"color": "black", "linewidth": 1.0},
        capprops={"color": "black", "linewidth": 1.0},
        zorder=3,
    )

    ax.set_xticks([1, 2])
    ax.set_xticklabels(
        [f"{ARM_LABELS[a]}\nn={v.size}" for a, v in zip(ARMS, arms)], fontsize=8
    )
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)
    if show_ylabel:
        ax.set_ylabel("relative frequency (%)", fontsize=9)

    # The reader should not have to hold the stats table alongside the figure.
    note = (
        f"p_adj = {stats_row.p_adj:.3f}\n"
        f"d = {stats_row.cliffs_delta:+.3f} ({stats_row.magnitude})"
    )
    ax.text(
        0.5,
        0.97,
        note,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=8,
        color="#333333",
        bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.85, "pad": 2.2},
    )


def _build_boxplot(df: pd.DataFrame) -> plt.Figure:
    """Build the two row responder figure from the all timepoints cohort frame.

    df is the long Part 3 cohort: one row per sample per population. The baseline
    row is derived from it rather than passed in, so both rows are guaranteed to
    come from the same query.
    """
    missing = REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"_build_boxplot needs columns {sorted(missing)}")

    baseline = df[df["time_from_treatment_start"] == 0]
    if baseline.empty:
        raise ValueError("no time 0 rows in the cohort frame")

    runs = [
        ("All timepoints", df),
        ("Baseline only (time 0)", baseline),
    ]
    populations = [p for p in POPULATION_ORDER if p in set(df["population"])]
    if not populations:
        raise ValueError("no populations to plot")

    rng = np.random.default_rng(JITTER_SEED)
    fig, axes = plt.subplots(
        len(runs),
        len(populations),
        figsize=(3.0 * len(populations), 4.4 * len(runs)),
        sharey=True,
        squeeze=False,
    )

    for row, (label, frame) in enumerate(runs):
        # run_analysis.py already reports sign disagreements. Annotating the
        # figure should not print them a second time.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            stats = compare_frame(frame, populations).set_index("population")
        for col, population in enumerate(populations):
            ax = axes[row][col]
            _panel(
                ax,
                frame[frame["population"] == population],
                stats.loc[population],
                rng,
                show_ylabel=(col == 0),
            )
            if row == 0:
                ax.set_title(population, fontsize=11)
        axes[row][0].annotate(
            label,
            xy=(0, 0.5),
            xytext=(-52, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            rotation=90,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
        )

    fig.suptitle(
        "Melanoma, miraclib, PBMC: cell population frequency by response", fontsize=13
    )
    fig.tight_layout(rect=(0.02, 0, 1, 0.98))
    return fig


def boxplot_responders(df: pd.DataFrame, outpath) -> Path:
    """Save the responder figure. Contracted entry point for run_analysis.py."""
    fig = _build_boxplot(df)
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    # Suppress the matplotlib version stamp so repeat runs hash the same.
    fig.savefig(outpath, dpi=150, metadata={"Software": None})
    plt.close(fig)
    return outpath

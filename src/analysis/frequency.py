"""Part 2. Per sample relative frequency of each cell population."""

from __future__ import annotations

import sqlite3

import pandas as pd

COLUMNS = ["sample", "total_count", "population", "count", "percentage"]

_QUERY = """
SELECT sample, total_count, population, count, percentage
FROM v_sample_frequencies
ORDER BY sample, population
"""


def get_frequencies(conn: sqlite3.Connection) -> pd.DataFrame:
    """Return the Part 2 table. The view already does the total and the percentage."""
    df = pd.read_sql_query(_QUERY, conn)
    if df.empty:
        raise ValueError("v_sample_frequencies returned no rows")
    return df[COLUMNS]

"""Connection helper so every module points at the same database file."""

from __future__ import annotations

import sqlite3
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DB_PATH = REPO_ROOT / "cell_count.db"


def connect(db_path: str | Path = DB_PATH) -> sqlite3.Connection:
    db_path = Path(db_path)
    if not db_path.exists():
        raise FileNotFoundError(f"{db_path} not found. Run: python load_data.py")
    conn = sqlite3.connect(db_path)
    # The views join across every table, so referential checks stay on.
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def population_order(conn: sqlite3.Connection) -> list[str]:
    """Population names in load order, so report rows do not shuffle run to run."""
    rows = conn.execute(
        "SELECT population_name FROM cell_populations ORDER BY population_id"
    ).fetchall()
    return [r[0] for r in rows]

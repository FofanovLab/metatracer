"""Count per-read combinations of selected annotation columns."""

from __future__ import annotations

import csv
import re
import sqlite3
import tempfile
from pathlib import Path
from typing import Iterable


def _normalize_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.strip().lower())


def _delimiter(path: Path) -> str:
    with path.open(encoding="utf-8", errors="replace") as handle:
        first = handle.readline()
    return "\t" if "\t" in first else ","


def _find_column(fields: list[str], name: str) -> str | None:
    wanted = _normalize_header(name)
    return next((field for field in fields if _normalize_header(field) == wanted), None)


def _values(value: str) -> set[str]:
    """Return distinct nonempty alternatives stored in one table cell."""
    text = (value or "").strip()
    if not text or text.upper() in {"NA", "N/A", "NONE", "-"}:
        return set()
    return {part.strip() for part in re.split(r"[;,]", text) if part.strip()}


def _ingest(
    connection: sqlite3.Connection,
    inputs: Iterable[Path],
    columns: list[str],
) -> tuple[int, int]:
    """Store each unique value observed for a read and selected column."""
    inserted = 0
    input_rows = 0
    batch: list[tuple[str, str, int, str]] = []
    read_batch: list[tuple[str, str]] = []
    statement = "INSERT OR IGNORE INTO candidates VALUES (?, ?, ?, ?)"

    for path in inputs:
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter=_delimiter(path))
            fields = reader.fieldnames or []
            read_column = _find_column(fields, "read_id")
            sample_column = _find_column(fields, "sample_id") or _find_column(fields, "sample")
            selected = [_find_column(fields, column) for column in columns]
            missing = [column for column, actual in zip(columns, selected) if actual is None]
            if read_column is None:
                missing.insert(0, "read_id")
            if missing:
                raise ValueError(
                    f"{path}: missing required count columns: {', '.join(missing)}"
                )

            fallback_sample = path.stem
            for row in reader:
                input_rows += 1
                read_id = (row.get(read_column) or "").strip()
                sample_id = (
                    (row.get(sample_column) or "").strip() if sample_column else fallback_sample
                )
                if not read_id:
                    continue
                read_batch.append((sample_id, read_id))
                for index, actual in enumerate(selected):
                    for value in _values(row.get(actual, "")):
                        batch.append((sample_id, read_id, index, value))
                if len(batch) >= 50_000:
                    before = connection.total_changes
                    connection.executemany(statement, batch)
                    inserted += connection.total_changes - before
                    batch.clear()
                    connection.executemany("INSERT OR IGNORE INTO reads VALUES (?, ?)", read_batch)
                    read_batch.clear()

    if batch:
        before = connection.total_changes
        connection.executemany(statement, batch)
        inserted += connection.total_changes - before
    if read_batch:
        connection.executemany("INSERT OR IGNORE INTO reads VALUES (?, ?)", read_batch)
    connection.commit()
    return inserted, input_rows


def _aggregate(connection: sqlite3.Connection, column_count: int) -> tuple[int, int]:
    """Collapse every read to one multi-value group and count those groups."""
    counts: dict[tuple[str, ...], int] = {}
    skipped_reads = 0
    current_read: tuple[str, str] | None = None
    values = [set() for _ in range(column_count)]

    def flush() -> None:
        nonlocal skipped_reads
        if current_read is None:
            return
        if any(not column_values for column_values in values):
            skipped_reads += 1
            return
        sample_id, _read_id = current_read
        grouped = tuple(";".join(sorted(column_values)) for column_values in values)
        key = (sample_id, *grouped)
        counts[key] = counts.get(key, 0) + 1

    rows = connection.execute(
        "SELECT reads.sample_id, reads.read_id, candidates.column_index, candidates.value "
        "FROM reads LEFT JOIN candidates USING (sample_id, read_id) "
        "ORDER BY reads.sample_id, reads.read_id, candidates.column_index, candidates.value"
    )
    for sample_id, read_id, column_index, value in rows:
        key = (sample_id, read_id)
        if current_read is not None and key != current_read:
            flush()
            values = [set() for _ in range(column_count)]
        current_read = key
        if column_index is not None:
            values[column_index].add(value)
    flush()

    connection.execute(
        "CREATE TABLE counts (group_key TEXT PRIMARY KEY, sample_id TEXT, "
        + ", ".join(f"column_{index} TEXT" for index in range(column_count))
        + ", read_count INTEGER)"
    )
    connection.executemany(
        "INSERT INTO counts VALUES (" + ",".join("?" for _ in range(column_count + 3)) + ")",
        [
            ("\x1f".join(key), *key, count)
            for key, count in counts.items()
        ],
    )
    connection.commit()
    return skipped_reads, len(counts)


def run(
    inputs: list[Path],
    output: Path,
    columns: list[str],
    tmpdir: Path | None = None,
) -> dict[str, int]:
    if not inputs:
        raise ValueError("At least one annotation input is required")
    if not columns:
        raise ValueError("At least one count column is required")
    normalized = [_normalize_header(column) for column in columns]
    if any(not column for column in normalized):
        raise ValueError("Count column names may not be empty")
    if len(set(normalized)) != len(normalized):
        raise ValueError("Count columns must be unique")

    paths = [Path(path) for path in inputs]
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temp_parent = str(tmpdir) if tmpdir else None
    with tempfile.TemporaryDirectory(prefix="metatracer-count-", dir=temp_parent) as work:
        connection = sqlite3.connect(Path(work) / "counts.sqlite")
        try:
            connection.execute("PRAGMA journal_mode=OFF")
            connection.execute("PRAGMA synchronous=OFF")
            connection.execute("PRAGMA temp_store=FILE")
            connection.execute(
                "CREATE TABLE reads (sample_id TEXT NOT NULL, read_id TEXT NOT NULL, "
                "PRIMARY KEY (sample_id, read_id)) WITHOUT ROWID"
            )
            connection.execute(
                "CREATE TABLE candidates (sample_id TEXT NOT NULL, read_id TEXT NOT NULL, "
                "column_index INTEGER NOT NULL, value TEXT NOT NULL, "
                "PRIMARY KEY (sample_id, read_id, column_index, value)) WITHOUT ROWID"
            )
            unique_values, input_rows = _ingest(connection, paths, columns)
            skipped_reads, count_rows = _aggregate(connection, len(columns))

            selected_sql = ", ".join(f"column_{index}" for index in range(len(columns)))
            order_sql = ", ".join(
                ["sample_id", *(f"column_{index}" for index in range(len(columns)))]
            )
            temporary = output.with_name(output.name + ".tmp")
            with temporary.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample_id", *columns, "count"])
                for row in connection.execute(
                    f"SELECT sample_id, {selected_sql}, read_count FROM counts ORDER BY {order_sql}"
                ):
                    writer.writerow(row)
            temporary.replace(output)
        finally:
            connection.close()

    return {
        "input_rows": input_rows,
        "unique_values": unique_values,
        "skipped_reads": skipped_reads,
        "count_rows": count_rows,
    }

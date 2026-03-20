#!/usr/bin/env python3
"""Monitor thermo output in a LAMMPS log while the run is still active.

Examples
--------
./watch_lammps_log.py run_mpi4_omp4.log -c Temp KinEng
./watch_lammps_log.py log.lammps -c temp,v_h2_count --history 20
./watch_lammps_log.py run_mpi4_omp4.log --list-columns
./watch_lammps_log.py run_mpi4_omp4.log -c Temp Press --no-follow
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


DEFAULT_PREFERRED_COLUMNS = [
    "Temp",
    "Press",
    "PotEng",
    "KinEng",
    "TotEng",
    "v_h2_count",
    "Atoms",
]


def normalize_name(name: str) -> str:
    return name.strip().lower()


def is_number(token: str) -> bool:
    try:
        float(token.replace("D", "E").replace("d", "e"))
        return True
    except ValueError:
        return False


def looks_like_header(tokens: list[str]) -> bool:
    return len(tokens) >= 2 and tokens[0] == "Step" and not any(is_number(tok) for tok in tokens[1:])


@dataclass
class ThermoSection:
    columns: list[str] = field(default_factory=list)
    rows: list[dict[str, str]] = field(default_factory=list)

    @property
    def has_data(self) -> bool:
        return bool(self.columns) and bool(self.rows)

    @property
    def lookup(self) -> dict[str, str]:
        return {normalize_name(name): name for name in self.columns}


class ThermoParser:
    def __init__(self) -> None:
        self.latest = ThermoSection()
        self.section_count = 0

    def reset(self) -> None:
        self.latest = ThermoSection()
        self.section_count = 0

    def feed_line(self, line: str) -> list[tuple[str, object]]:
        stripped = line.strip()
        if not stripped:
            return []

        tokens = stripped.split()
        if looks_like_header(tokens):
            self.latest = ThermoSection(columns=tokens, rows=[])
            self.section_count += 1
            return [("header", tokens)]

        if not self.latest.columns:
            return []

        if len(tokens) != len(self.latest.columns):
            return []

        if not all(is_number(token) for token in tokens):
            return []

        row = dict(zip(self.latest.columns, tokens))
        self.latest.rows.append(row)
        return [("row", row)]


def parse_column_args(raw_values: list[str] | None) -> list[str]:
    if not raw_values:
        return []

    requested: list[str] = []
    for value in raw_values:
        for item in value.split(","):
            item = item.strip()
            if item:
                requested.append(item)
    return requested


def choose_default_columns(section: ThermoSection) -> list[str]:
    selected = ["Step"]
    lookup = section.lookup
    for candidate in DEFAULT_PREFERRED_COLUMNS:
        name = lookup.get(normalize_name(candidate))
        if name and name not in selected:
            selected.append(name)

    if len(selected) == 1:
        for name in section.columns:
            if name not in selected:
                selected.append(name)
            if len(selected) >= min(6, len(section.columns)):
                break
    return selected


def resolve_columns(section: ThermoSection, requested: list[str]) -> list[str]:
    if not section.columns:
        return []

    if not requested:
        return choose_default_columns(section)

    lookup = section.lookup
    available_normalized = list(lookup.keys())
    resolved = ["Step"] if "Step" in section.columns else []

    for item in requested:
        key = normalize_name(item)
        exact = lookup.get(key)
        if exact:
            if exact not in resolved:
                resolved.append(exact)
            continue

        matches = [lookup[name] for name in available_normalized if key in name]
        matches = list(dict.fromkeys(matches))
        if len(matches) == 1:
            if matches[0] not in resolved:
                resolved.append(matches[0])
            continue

        available = ", ".join(section.columns)
        if matches:
            joined = ", ".join(matches)
            raise SystemExit(
                f"Ambiguous column '{item}'. Matches: {joined}\nAvailable columns: {available}"
            )
        raise SystemExit(f"Unknown column '{item}'. Available columns: {available}")

    return resolved


def format_table(rows: list[dict[str, str]], columns: list[str]) -> str:
    if not rows:
        return "No thermo rows parsed yet."

    widths: dict[str, int] = {}
    for column in columns:
        widths[column] = len(column)
        for row in rows:
            widths[column] = max(widths[column], len(row.get(column, "")))

    header = "  ".join(column.rjust(widths[column]) for column in columns)
    separator = "  ".join("-" * widths[column] for column in columns)
    body = [
        "  ".join(row.get(column, "").rjust(widths[column]) for column in columns)
        for row in rows
    ]
    return "\n".join([header, separator, *body])


def format_stats(section: ThermoSection, columns: list[str]) -> str:
    stats_lines: list[str] = []
    for column in columns:
        if column == "Step":
            continue
        values: list[float] = []
        for row in section.rows:
            token = row.get(column)
            if token is None:
                continue
            try:
                value = float(token.replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
            if math.isfinite(value):
                values.append(value)
        if not values:
            continue
        avg = sum(values) / len(values)
        stats_lines.append(
            f"{column}: min={min(values):.6g} avg={avg:.6g} max={max(values):.6g}"
        )
    return "\n".join(stats_lines)


def render_screen(
    path: Path,
    parser: ThermoParser,
    columns: list[str],
    history: int,
    interval: float,
    show_stats: bool,
) -> str:
    section = parser.latest
    lines = [
        f"File: {path}",
        f"Thermo sections seen: {parser.section_count} | Rows in latest section: {len(section.rows)} | Poll interval: {interval:.2f}s",
    ]
    if section.rows:
        latest = section.rows[-1]
        latest_summary = "  ".join(f"{col}={latest.get(col, '')}" for col in columns)
        lines.append(f"Latest row: {latest_summary}")
    else:
        lines.append("Waiting for thermo rows...")

    lines.append("")
    lines.append(format_table(section.rows[-history:], columns))
    if show_stats and section.rows:
        stats_block = format_stats(section, columns)
        if stats_block:
            lines.extend(["", "Stats over latest thermo section:", stats_block])
    return "\n".join(lines)


def scan_entire_file(path: Path, parser: ThermoParser) -> tuple[int, int]:
    parser.reset()
    inode = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parser.feed_line(line)
        position = handle.tell()
        inode = os.fstat(handle.fileno()).st_ino
    return position, inode


def render_once(
    path: Path,
    parser: ThermoParser,
    columns: list[str],
    history: int,
    interval: float,
    show_stats: bool,
    line_mode: bool,
) -> None:
    if line_mode:
        print(render_screen(path, parser, columns, history, interval, show_stats), flush=True)
        return

    print("\033[2J\033[H", end="")
    print(render_screen(path, parser, columns, history, interval, show_stats))
    sys.stdout.flush()


def resolve_default_logfile() -> Path:
    for candidate in ("log.lammps", "run_mpi4_omp4.log"):
        path = Path(candidate)
        if path.exists():
            return path
    raise SystemExit("No logfile argument provided and neither 'log.lammps' nor 'run_mpi4_omp4.log' exists in the current directory.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Monitor the latest LAMMPS thermo section in a log file."
    )
    parser.add_argument(
        "logfile",
        nargs="?",
        help="LAMMPS log to monitor. Defaults to ./log.lammps or ./run_mpi4_omp4.log if present.",
    )
    parser.add_argument(
        "-c",
        "--columns",
        nargs="+",
        help="Columns to display. Case-insensitive. Accepts space-separated or comma-separated names.",
    )
    parser.add_argument(
        "-n",
        "--history",
        type=int,
        default=12,
        help="Number of recent thermo rows to keep on screen.",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=1.0,
        help="Polling interval in seconds while following the log.",
    )
    parser.add_argument(
        "--list-columns",
        action="store_true",
        help="Print the column names from the latest thermo section and exit.",
    )
    parser.add_argument(
        "--no-follow",
        action="store_true",
        help="Parse the current file once and exit.",
    )
    parser.add_argument(
        "--line-mode",
        action="store_true",
        help="Do not clear the terminal between refreshes.",
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Show min/avg/max for selected numeric columns in the latest thermo section.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    path = Path(args.logfile) if args.logfile else resolve_default_logfile()
    path = path.expanduser().resolve()

    if not path.exists():
        raise SystemExit(f"Log file does not exist: {path}")

    if args.history <= 0:
        raise SystemExit("--history must be >= 1")

    if args.interval <= 0:
        raise SystemExit("--interval must be > 0")

    parser = ThermoParser()
    position, inode = scan_entire_file(path, parser)

    if not parser.latest.columns:
        raise SystemExit(f"No LAMMPS thermo header found in {path}")

    if args.list_columns:
        print("\n".join(parser.latest.columns))
        return 0

    requested = parse_column_args(args.columns)
    selected_columns = resolve_columns(parser.latest, requested)

    render_once(
        path=path,
        parser=parser,
        columns=selected_columns,
        history=args.history,
        interval=args.interval,
        show_stats=args.stats,
        line_mode=args.line_mode or args.no_follow or not sys.stdout.isatty(),
    )

    if args.no_follow:
        return 0

    partial = ""
    while True:
        try:
            stat = path.stat()
        except FileNotFoundError:
            time.sleep(args.interval)
            continue

        rotated = stat.st_ino != inode or stat.st_size < position
        if rotated:
            position, inode = scan_entire_file(path, parser)
            selected_columns = resolve_columns(parser.latest, requested)
            render_once(
                path=path,
                parser=parser,
                columns=selected_columns,
                history=args.history,
                interval=args.interval,
                show_stats=args.stats,
                line_mode=args.line_mode or not sys.stdout.isatty(),
            )
            time.sleep(args.interval)
            continue

        with path.open("r", encoding="utf-8", errors="replace") as handle:
            handle.seek(position)
            chunk = handle.read()
            position = handle.tell()

        if chunk:
            partial += chunk
            lines = partial.splitlines(keepends=True)
            partial = ""
            if lines and not lines[-1].endswith(("\n", "\r")):
                partial = lines.pop()

            updated = False
            for line in lines:
                events = parser.feed_line(line)
                if events:
                    updated = True

            if updated:
                selected_columns = resolve_columns(parser.latest, requested)
                render_once(
                    path=path,
                    parser=parser,
                    columns=selected_columns,
                    history=args.history,
                    interval=args.interval,
                    show_stats=args.stats,
                    line_mode=args.line_mode or not sys.stdout.isatty(),
                )

        time.sleep(args.interval)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nStopped.", file=sys.stderr)
        raise SystemExit(130)

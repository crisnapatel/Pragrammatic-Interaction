#!/usr/bin/env python3
"""Live terminal plots for LAMMPS thermo output.

Examples
--------
python3 plot_lammps_log_tty.py run_mpi4_omp4.log -x Step -y Temp
python3 plot_lammps_log_tty.py log.lammps -x Time -y temp,v_h2_count -n 120
python3 plot_lammps_log_tty.py run_mpi4_omp4.log -x Step -y Temp --width 70 --height 12
python3 plot_lammps_log_tty.py run_mpi4_omp4.log --list-columns

Sizing
------
Use ``--width`` to limit the plot width in terminal columns.
Use ``--height`` to limit the total plot height in terminal rows.
Example: ``--width 70 --height 12``
"""

from __future__ import annotations

import argparse
import math
import shutil
import sys
import time
from pathlib import Path

from watch_lammps_log import (
    ThermoParser,
    normalize_name,
    parse_column_args,
    resolve_default_logfile,
    scan_entire_file,
)

try:
    import plotext as plt
except ImportError:
    plt = None


def coerce_float(token: str | None) -> float | None:
    if token is None:
        return None
    try:
        value = float(token.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def resolve_one_column(section, requested: str | None, default: str) -> str:
    lookup = section.lookup
    key = normalize_name(requested or default)

    exact = lookup.get(key)
    if exact:
        return exact

    matches = [actual for name, actual in lookup.items() if key in name]
    matches = list(dict.fromkeys(matches))
    if len(matches) == 1:
        return matches[0]
    if matches:
        raise SystemExit(f"Ambiguous x-column '{requested}'. Matches: {', '.join(matches)}")
    raise SystemExit(f"Unknown x-column '{requested or default}'. Available: {', '.join(section.columns)}")


def resolve_y_columns(section, requested: list[str]) -> list[str]:
    lookup = section.lookup
    if not requested:
        raise SystemExit("At least one y-column is required.")

    resolved: list[str] = []
    for item in requested:
        key = normalize_name(item)
        exact = lookup.get(key)
        if exact:
            if exact not in resolved:
                resolved.append(exact)
            continue

        matches = [actual for name, actual in lookup.items() if key in name]
        matches = list(dict.fromkeys(matches))
        if len(matches) == 1:
            if matches[0] not in resolved:
                resolved.append(matches[0])
            continue
        if matches:
            raise SystemExit(f"Ambiguous y-column '{item}'. Matches: {', '.join(matches)}")
        raise SystemExit(f"Unknown y-column '{item}'. Available: {', '.join(section.columns)}")

    return resolved


def format_number(value: float | None) -> str:
    if value is None:
        return "n/a"
    abs_value = abs(value)
    if abs_value >= 1e4 or (abs_value > 0 and abs_value < 1e-3):
        return f"{value:.3e}"
    return f"{value:.6f}".rstrip("0").rstrip(".")


def build_canvas(width: int, height: int) -> list[list[str]]:
    return [[" " for _ in range(width)] for _ in range(height)]


def canvas_to_lines(canvas: list[list[str]]) -> list[str]:
    return ["".join(row).rstrip() for row in canvas]


def draw_line_plot(
    x_values: list[float | None],
    y_values: list[float | None],
    width: int,
    height: int,
) -> tuple[list[str], float | None, float | None]:
    if width < 8 or height < 4:
        return ["terminal too small"], None, None

    finite_pairs = [
        (x, y)
        for x, y in zip(x_values, y_values)
        if x is not None and y is not None and math.isfinite(x) and math.isfinite(y)
    ]
    if not finite_pairs:
        return ["no numeric data yet"], None, None

    x_min = min(x for x, _ in finite_pairs)
    x_max = max(x for x, _ in finite_pairs)
    y_min = min(y for _, y in finite_pairs)
    y_max = max(y for _, y in finite_pairs)

    if x_max == x_min:
        x_max = x_min + 1.0
    if y_max == y_min:
        bump = abs(y_min) * 0.05 or 1.0
        y_min -= bump
        y_max += bump

    plot_w = width - 8
    plot_h = height - 2
    canvas = build_canvas(width, height)

    for row in range(plot_h):
        canvas[row][7] = "|"
    for col in range(7, width):
        canvas[plot_h][col] = "-"
    canvas[plot_h][7] = "+"

    for idx, (x, y) in enumerate(finite_pairs):
        px = 8 + int(round((plot_w - 1) * (x - x_min) / (x_max - x_min)))
        py = int(round((plot_h - 1) * (y - y_min) / (y_max - y_min)))
        py = (plot_h - 1) - py
        px = max(8, min(width - 1, px))
        py = max(0, min(plot_h - 1, py))
        canvas[py][px] = "*" if idx == len(finite_pairs) - 1 else "."

    top_label = f"{format_number(y_max):>7}"
    bot_label = f"{format_number(y_min):>7}"
    for i, ch in enumerate(top_label[:7]):
        canvas[0][i] = ch
    for i, ch in enumerate(bot_label[:7]):
        canvas[plot_h - 1][i] = ch

    x0 = format_number(x_min)
    x1 = format_number(x_max)
    x_line = [" " for _ in range(width)]
    for i, ch in enumerate(x0[: max(0, width // 2 - 1)]):
        x_line[8 + i] = ch
    start = max(8, width - len(x1))
    for i, ch in enumerate(x1):
        if start + i < width:
            x_line[start + i] = ch
    canvas[height - 1] = x_line

    return canvas_to_lines(canvas), finite_pairs[-1][1], finite_pairs[0][1]


def render(
    path: Path,
    parser: ThermoParser,
    x_column: str,
    y_columns: list[str],
    history: int,
    interval: float,
    plot_width: int | None,
) -> str:
    terminal = shutil.get_terminal_size((120, 40))
    width = min(plot_width or terminal.columns, terminal.columns)
    width = max(40, width)
    section = parser.latest
    rows = section.rows[-history:]
    x_values = [coerce_float(row.get(x_column)) for row in rows]

    lines = [
        f"File: {path}",
        f"Latest thermo section rows: {len(section.rows)} | x={x_column} | poll={interval:.2f}s | history={len(rows)}",
    ]

    if rows:
        latest = rows[-1]
        summary_fields = list(dict.fromkeys(["Step", x_column, *y_columns]))
        summary = "  ".join(f"{name}={latest.get(name, '')}" for name in summary_fields if name in latest)
        lines.append(f"Latest row: {summary}")
    else:
        lines.append("Waiting for thermo rows...")

    for column in y_columns:
        y_values = [coerce_float(row.get(column)) for row in rows]
        finite = [value for value in y_values if value is not None]
        lines.append("")
        if finite:
            lines.append(
                f"{column}  latest={format_number(finite[-1])}  min={format_number(min(finite))}  max={format_number(max(finite))}"
            )
        else:
            lines.append(f"{column}  no numeric data yet")
        plot_lines, _, _ = draw_line_plot(x_values, y_values, width, 10)
        lines.extend(plot_lines)

    return "\n".join(lines)


def render_with_plotext(
    path: Path,
    parser: ThermoParser,
    x_column: str,
    y_columns: list[str],
    history: int,
    interval: float,
    line_mode: bool,
    plot_width: int | None,
    plot_height: int | None,
) -> None:
    if plt is None:
        text = render(path, parser, x_column, y_columns, history, interval, plot_width)
        render_once(text, line_mode)
        if not line_mode:
            print(
                "\nplotext is not available in this Python environment. "
                "Run with `conda run -n Analysis python3 plot_lammps_log_tty.py ...` for better terminal plots.",
                file=sys.stderr,
            )
        return

    terminal = shutil.get_terminal_size((120, 40))
    width = min(plot_width or terminal.columns, terminal.columns)
    width = max(60, width)
    header_lines = 4
    total_plot_height = min(plot_height or (terminal.lines - header_lines), terminal.lines - header_lines)
    total_plot_height = max(8, total_plot_height)

    section = parser.latest
    rows = section.rows[-history:]
    x_values = [coerce_float(row.get(x_column)) for row in rows]

    if not line_mode:
        plt.clear_terminal()

    print(f"File: {path}")
    print(
        f"Latest thermo section rows: {len(section.rows)} | x={x_column} | poll={interval:.2f}s | history={len(rows)}"
    )
    if rows:
        latest = rows[-1]
        summary_fields = list(dict.fromkeys(["Step", x_column, *y_columns]))
        summary = "  ".join(f"{name}={latest.get(name, '')}" for name in summary_fields if name in latest)
        print(f"Latest row: {summary}")
    else:
        print("Waiting for thermo rows...")
    print("")

    plt.clear_data()
    plt.clear_figure()
    plt.clf()
    plt.theme("pro")
    plt.plotsize(width, total_plot_height)
    use_subplots = len(y_columns) > 1
    if use_subplots:
        plt.subplots(len(y_columns), 1)

    colors = ["cyan", "orange", "green", "magenta", "red", "blue", "yellow"]
    for index, column in enumerate(y_columns, start=1):
        y_values = [coerce_float(row.get(column)) for row in rows]
        finite = [value for value in y_values if value is not None]
        if use_subplots:
            plt.subplot(index, 1)
        plt.title(
            f"{column}  latest={format_number(finite[-1]) if finite else 'n/a'}  "
            f"min={format_number(min(finite)) if finite else 'n/a'}  "
            f"max={format_number(max(finite)) if finite else 'n/a'}"
        )
        if finite:
            xs = []
            ys = []
            for x, y in zip(x_values, y_values):
                if x is None or y is None:
                    continue
                xs.append(x)
                ys.append(y)
            if xs and ys:
                plt.plot(xs, ys, color=colors[(index - 1) % len(colors)], marker="braille")
        else:
            plt.plot([0], [0], color=colors[(index - 1) % len(colors)])
        plt.ylabel(column)
        if index == len(y_columns):
            plt.xlabel(x_column)

    plt.show()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot LAMMPS thermo columns directly in the terminal.")
    parser.add_argument(
        "logfile",
        nargs="?",
        help="LAMMPS log to monitor. Defaults to ./log.lammps or ./run_mpi4_omp4.log if present.",
    )
    parser.add_argument(
        "-x",
        "--x-column",
        help="Column for the x-axis.",
    )
    parser.add_argument(
        "-y",
        "--y-columns",
        nargs="+",
        help="Thermo columns to plot on the y-axis. Case-insensitive. Accepts space-separated or comma-separated names.",
    )
    parser.add_argument(
        "-n",
        "--history",
        type=int,
        default=100,
        help="Number of recent points to plot.",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=1.0,
        help="Polling interval in seconds.",
    )
    parser.add_argument(
        "--list-columns",
        action="store_true",
        help="Print thermo column names from the latest section and exit.",
    )
    parser.add_argument(
        "--no-follow",
        action="store_true",
        help="Render the current plot once and exit.",
    )
    parser.add_argument(
        "--line-mode",
        action="store_true",
        help="Do not clear the terminal between refreshes.",
    )
    parser.add_argument(
        "--width",
        type=int,
        help="Plot width in terminal columns. Defaults to the current terminal width.",
    )
    parser.add_argument(
        "--height",
        type=int,
        help="Total plot height in terminal rows. Defaults to available terminal height.",
    )
    return parser


def render_once(text: str, line_mode: bool) -> None:
    if line_mode:
        print(text, flush=True)
        return
    print("\033[2J\033[H", end="")
    print(text)
    sys.stdout.flush()


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
    if args.width is not None and args.width < 40:
        raise SystemExit("--width must be >= 40")
    if args.height is not None and args.height < 8:
        raise SystemExit("--height must be >= 8")

    parser = ThermoParser()
    position, inode = scan_entire_file(path, parser)

    if not parser.latest.columns:
        raise SystemExit(f"No LAMMPS thermo header found in {path}")

    if args.list_columns:
        print("\n".join(parser.latest.columns))
        return 0

    if not args.x_column or not args.y_columns:
        raise SystemExit("Explicit axes are required. Use --x-column <name> and --y-columns <name> [more names].")

    requested_columns = parse_column_args(args.y_columns)
    resolved_y = resolve_y_columns(parser.latest, requested_columns)
    x_column = resolve_one_column(parser.latest, args.x_column, "Step")
    if x_column in resolved_y:
        raise SystemExit(f"x-column '{x_column}' cannot also appear in y-columns.")

    if not resolved_y:
        raise SystemExit("No y-columns left to plot after resolving x-column.")

    line_mode = args.line_mode or args.no_follow or not sys.stdout.isatty()
    render_with_plotext(
        path,
        parser,
        x_column,
        resolved_y,
        args.history,
        args.interval,
        line_mode,
        args.width,
        args.height,
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
            x_column = resolve_one_column(parser.latest, args.x_column, "Step")
            resolved_y = resolve_y_columns(parser.latest, requested_columns)
            render_with_plotext(
                path,
                parser,
                x_column,
                resolved_y,
                args.history,
                args.interval,
                line_mode,
                args.width,
                args.height,
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
                if parser.feed_line(line):
                    updated = True

            if updated:
                x_column = resolve_one_column(parser.latest, args.x_column, "Step")
                resolved_y = resolve_y_columns(parser.latest, requested_columns)
                render_with_plotext(
                    path,
                    parser,
                    x_column,
                    resolved_y,
                    args.history,
                    args.interval,
                    line_mode,
                    args.width,
                    args.height,
                )

        time.sleep(args.interval)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nStopped.", file=sys.stderr)
        raise SystemExit(130)

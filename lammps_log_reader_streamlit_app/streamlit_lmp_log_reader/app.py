
# LAMMPS Log Lite v2.5.2 — Streamlit (hotfix)

import io, os, re, gzip, json, hashlib, time
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import plotly.io as pio

st.set_page_config(page_title="LAMMPS Log Lite v2.5.2", layout="wide")

# ----- Parsers -----
HAS_LL = False
HAS_OFFICIAL = False
LL_MSG, OFF_MSG = "", ""
try:
    import lammps_logfile as ll
    HAS_LL = True
    LL_MSG = "lammps-logfile OK"
except Exception as e:
    LL_MSG = f"lammps-logfile import failed: {e}"

try:
    from lammps.formats import LogFile as OfficialLogFile
    HAS_OFFICIAL = True
    OFF_MSG = "lammps.formats.LogFile OK"
except Exception as e:
    OFF_MSG = f"official parser not available: {e}"

# ----- Regexes & helpers -----
HEADER_RE = re.compile(r'^\s*Step(\s+\S+)+\s*$')
LOOP_RE = re.compile(r'^\s*Loop time', re.IGNORECASE)
WARN_RE = re.compile(r'^\s*WARNING:', re.IGNORECASE)
ERR_RE = re.compile(r'^\s*ERROR:', re.IGNORECASE)
BREAKDOWN_RE_OLD = re.compile(r'^\s*(Pair|Neigh|Comm|Output|Modify|Other)\s+time\s*\(%\)\s*[:=]\s*([+\-.\deE]+)\s*\(([\d.\-eE]+)\)', re.IGNORECASE)
BREAKDOWN_TABLE_HEADER = re.compile(r'^\s*MPI task timing breakdown', re.IGNORECASE)
BREAKDOWN_TABLE_ROW = re.compile(r'^\s*(\w+)\s*\|\s*([+\-.\deE]*)\s*\|\s*([+\-.\deE]*)\s*\|\s*([+\-.\deE]*)\s*\|\s*([+\-.\deE]*)\s*\|\s*([+\-.\deE]*)')

def is_gzip(bytez: bytes) -> bool:
    return len(bytez) >= 2 and bytez[0] == 0x1F and bytez[1] == 0x8B

def to_text(bytez: bytes) -> str:
    if isinstance(bytez, str):
        return bytez
    if is_gzip(bytez):
        data = gzip.decompress(bytez)
    else:
        data = bytez
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return data.decode("latin-1", errors="ignore")

def sha1(data: bytes) -> str:
    return hashlib.sha1(data).hexdigest()

@dataclass
class Scope:
    file_key: str
    display: str
    run_num: int  # -1 means "all runs"

@dataclass
class ParsedFile:
    name: str
    text: str
    hash: str
    ll_file: Optional[object]
    official: Optional[object]
    runs: List[Tuple[List[str], int, int]]
    lines: List[str]

def split_runs(text: str) -> List[Tuple[List[str], int, int]]:
    lines = text.splitlines()
    runs = []
    collecting = False
    header = None
    start = None
    for i, line in enumerate(lines):
        if not collecting and HEADER_RE.match(line):
            header = line.split()
            start = i + 1
            collecting = True
            continue
        if collecting:
            if LOOP_RE.match(line) or (not line.strip()):
                runs.append((header, start, i))
                collecting = False
    if collecting:
        runs.append((header, start, len(lines)-1))
    return runs

def parse_file(name: str, data: bytes) -> ParsedFile:
    text = to_text(data)
    h = sha1(data if isinstance(data, (bytes, bytearray)) else text.encode("utf-8", errors="ignore"))
    ll_file = None
    if HAS_LL:
        try:
            ll_file = ll.File(io.StringIO(text))
        except Exception:
            ll_file = None
    official_file = None
    if HAS_OFFICIAL:
        try:
            official_file = OfficialLogFile(io.StringIO(text))
        except Exception:
            official_file = None
    runs = split_runs(text)
    lines = text.splitlines()
    return ParsedFile(name=name, text=text, hash=h, ll_file=ll_file, official=official_file, runs=runs, lines=lines)

def get_keywords_union_ll(pf: ParsedFile) -> List[str]:
    """Union of keywords across runs without an unbounded loop."""
    cols = set()
    n_runs = max(0, len(pf.runs))
    # hard failsafe cap to avoid pathological files
    cap = min(n_runs if n_runs else 32, 256)
    if pf.ll_file is not None:
        for i in range(cap):
            try:
                ks = pf.ll_file.get_keywords(run_num=i)
            except Exception:
                break
            if not ks:
                # If the parser returns empty for an out-of-range run, stop.
                if i >= n_runs:
                    break
                # else keep going (this run may be empty header-only)
            cols.update(ks or [])
    if not cols:
        for h, _, _ in pf.runs:
            cols.update(h or [])
    return sorted(list(cols))

def get_keywords(pf: ParsedFile, run_num: int) -> List[str]:
    if run_num == -1:
        return get_keywords_union_ll(pf)
    if pf.ll_file is not None:
        try:
            ks = pf.ll_file.get_keywords(run_num=run_num)
            if ks: return ks
        except Exception:
            pass
    if pf.official is not None:
        try:
            return list(pf.official.runs[run_num].keys())
        except Exception:
            pass
    return pf.runs[run_num][0] if (0 <= run_num < len(pf.runs)) else []

def get_series(pf: ParsedFile, key: str, run_num: int) -> Optional[np.ndarray]:
    if pf.ll_file is not None:
        try:
            arr = pf.ll_file.get(key, run_num=run_num)
            return None if arr is None else np.asarray(arr, dtype=float)
        except Exception:
            pass
    if pf.official is not None:
        try:
            if run_num == -1:
                vals = []
                for r in pf.official.runs:
                    if key in r:
                        vals.extend(r[key])
                return np.asarray(vals, dtype=float) if len(vals) else None
            else:
                run = pf.official.runs[run_num]
                if key in run:
                    return np.asarray(run[key], dtype=float)
                return None
        except Exception:
            pass
    return None

# ---------- Formula helpers ----------
IDENT_RE = re.compile(r'\b[A-Za-z_][A-Za-z0-9_]*\b')
PY_FUNCS = {"abs","sqrt","log","exp","sin","cos","tan","arcsin","arccos","arctan","where","clip","mean","std","min","max"}

def ensure_columns_for_formula(pf: ParsedFile, run_num: int, df: pd.DataFrame, expr: str, all_cols: List[str]) -> pd.DataFrame:
    names = set(IDENT_RE.findall(expr))
    want = [n for n in names if n in all_cols and n not in df.columns and n not in PY_FUNCS]
    for n in want:
        arr = get_series(pf, n, run_num)
        if arr is not None:
            df[n] = arr
    return df

# ----- Session state -----
if "parsed" not in st.session_state:
    st.session_state.parsed = {}
if "formulas" not in st.session_state:
    st.session_state.formulas = []
if "live_active" not in st.session_state:
    st.session_state.live_active = False
if "live_last" not in st.session_state:
    st.session_state.live_last = 0.0

st.title("LAMMPS Log Lite v2.5.2")
st.caption("Multi-file overlays, formulas, LTTB, dual axis, live tailing, performance & errors dashboards.")

with st.expander("Parsers diagnostics", expanded=False):
    st.write(LL_MSG); st.write(OFF_MSG)

# ---------- Upload & Live tail ----------
st.subheader("Upload logs")
uploads = st.file_uploader("Drag & drop one or more log files (.lammps/.log/.txt/.gz)",
                           type=["lammps","log","txt","gz"], accept_multiple_files=True)

with st.expander("Live tailing (path on this machine)", expanded=False):
    st.caption("Note: Uploading a file **does not** reveal its path. Enter an absolute path here.")
    c1, c2, c3, c4 = st.columns([3,1,1,1])
    live_path = c1.text_input("Path", value="", placeholder="/absolute/path/to/log.lammps")
    interval = c2.number_input("Refresh (s)", min_value=2, max_value=60, value=5, step=1)
    start = c3.button("Start"); stop = c4.button("Stop")
    if start and live_path:
        st.session_state.live_active = True; st.session_state.live_path = live_path
    if stop:
        st.session_state.live_active = False
    if st.session_state.live_active and st.session_state.get("live_path"):
        try:
            with open(st.session_state.live_path, "rb") as f: live_data = f.read()
            class SimpleFile: 
                def __init__(self, name, data): self.name = name; self._data = data
                def read(self): return self._data
            uploads = list(uploads) if uploads else []
            uploads.append(SimpleFile(os.path.basename(st.session_state.live_path), live_data))
            st.success(f"Live tail OK: {st.session_state.live_path}")
            now = time.time()
            if now - st.session_state.live_last >= interval:
                st.session_state.live_last = now
                time.sleep(interval); st.experimental_rerun()
        except Exception as e:
            st.error(f"Live read failed: {e}")

# Parse all uploaded
files: Dict[str, ParsedFile] = {}
if uploads:
    for up in uploads:
        data = up.read()
        pf = parse_file(up.name, data)
        files[pf.hash] = pf
        st.session_state.parsed[pf.hash] = pf
elif st.session_state.parsed:
    files = st.session_state.parsed

if not files:
    st.info("Upload at least one log to continue.")
    st.stop()

def scope_list_from_files(files: Dict[str, ParsedFile]) -> List[Scope]:
    scopes = []
    for h, pf in files.items():
        scopes.append(Scope(file_key=h, display=f"{pf.name} — All runs", run_num=-1))
        for i in range(len(pf.runs)):
            scopes.append(Scope(file_key=h, display=f"{pf.name} — Run #{i+1}", run_num=i))
    return scopes

scopes = scope_list_from_files(files)

tab_plot, tab_errs, tab_perf, tab_search, tab_presets = st.tabs(["Plot", "Errors & Warnings", "Performance", "Prefix search", "View presets"])

# ================= Plot Tab ==================
with tab_plot:
    r1c1, r1c2, r1c3, r1c4 = st.columns([1.4, 1.4, 1.2, 1.0])
    mode = r1c1.radio("Mode", ["Single scope (multi Y)", "Overlay (one Y)"])
    down_method = r1c2.selectbox("Downsampling", ["None", "Every Nth", "LTTB"])
    avg_N = r1c3.number_input("Running avg N", min_value=1, value=1, step=1)
    if down_method == "Every Nth":
        ds_every = r1c4.number_input("Every Nth", min_value=1, value=1, step=1); lttb_pts = 0
    else:
        lttb_pts = r1c4.number_input("Target pts (LTTB)", min_value=500, max_value=20000, value=5000, step=500) if down_method == "LTTB" else 0
        ds_every = 1

    if mode == "Single scope (multi Y)":
        scol1, scol2 = st.columns([2,1])
        sel_scope = scol1.selectbox("Scope", options=list(range(len(scopes))), format_func=lambda i: scopes[i].display)
        dpolicy = scol2.selectbox("Dup. timestep policy", ["keep_all","keep_first","keep_last"])

        pf = files[scopes[sel_scope].file_key]
        run_num = scopes[sel_scope].run_num
        cols = get_keywords(pf, run_num)
        if not cols:
            st.warning("No columns detected for this scope.")
            st.stop()

        c2a, c2b = st.columns([1,1])
        x_default = "Time" if "Time" in cols else ("Step" if "Step" in cols else cols[0])
        x_col = c2a.selectbox("X", cols, index=cols.index(x_default))
        y_choices = [c for c in cols if c != x_col]
        y_defaults = [c for c in ["Temp","Density","Press","PotEng","TotEng","KinEng"] if c in y_choices] or y_choices[:1]
        y_left = c2b.multiselect("Y (left)", y_choices, default=y_defaults)

        c3a, c3b = st.columns([1,1])
        y_right = c3a.selectbox("Y (right)", ["(none)"] + y_choices, index=0)

        # Formulas
        st.markdown("**Derived series** (use column names; e.g., `Density-200`, `(PotEng+KinEng)`):")
        fcol1, fcol2, fcol3 = st.columns([3,1,1])
        f_expr = fcol1.text_input("Formula", value="", placeholder="e.g., (PotEng+KinEng)/Step")
        if fcol2.button("Add"):
            if f_expr.strip():
                st.session_state.formulas.append(f_expr.strip())
        if fcol3.button("Clear all"):
            st.session_state.formulas = []
        if st.session_state.formulas:
            st.caption("Active: " + ", ".join([f"`{f}`" for f in st.session_state.formulas]))

        # Build dataframe
        x = get_series(pf, x_col, run_num)
        if x is None:
            st.error(f"Cannot read X column '{x_col}'.")
            st.stop()
        df = pd.DataFrame({x_col: x})
        for y in set(y_left + ([] if y_right == "(none)" else [y_right])):
            arr = get_series(pf, y, run_num)
            if arr is not None:
                df[y] = arr

        # Evaluate formulas: autoload needed columns referenced
        all_cols = get_keywords(pf, run_num)
        for expr in st.session_state.formulas:
            try:
                df = ensure_columns_for_formula(pf, run_num, df, expr, all_cols)
                series = df.eval(expr, engine="python")
                df[f"f:{expr}"] = pd.Series(series)
            except Exception as e:
                st.warning(f"Formula `{expr}` failed: {e}")

        if run_num == -1:
            df = df.dropna()
            if dpolicy != "keep_all":
                df = df.drop_duplicates(subset=[x_col], keep=("first" if dpolicy=="keep_first" else "last"))

        if avg_N > 1:
            for c in df.columns:
                if c != x_col:
                    df[c] = pd.Series(df[c]).rolling(window=int(avg_N), min_periods=1, center=True).mean()

        if down_method == "Every Nth" and ds_every > 1:
            df = df.iloc[::int(ds_every), :].reset_index(drop=True)
        elif down_method == "LTTB" and lttb_pts > 0:
            new = {x_col: None}
            xx = df[x_col].to_numpy()
            for c in df.columns:
                if c == x_col: continue
                yy = pd.to_numeric(df[c], errors="coerce").to_numpy()
                mask = ~np.isnan(yy); x = xx[mask]; y = yy[mask]
                if len(x) > 2:
                    n, T = len(x), int(lttb_pts)
                    if 3 <= T < n:
                        every = (n - 2) / (T - 2)
                        a = 0
                        xds = [x[0]]; yds = [y[0]]
                        for i in range(0, T-2):
                            ars = int(np.floor((i + 1) * every) + 1)
                            are = int(np.floor((i + 2) * every) + 1)
                            are = min(are, n)
                            avg_x = np.mean(x[ars:are]) if are > ars else x[ars-1]
                            avg_y = np.mean(y[ars:are]) if are > ars else y[ars-1]
                            ro = int(np.floor(i * every) + 1)
                            rt = int(np.floor((i + 1) * every) + 1)
                            rt = min(rt, n-1)
                            bx = x[ro:rt]; by = y[ro:rt]; ax = x[a]; ay = y[a]
                            area = np.abs((ax - avg_x) * (by - ay) - (ax - bx) * (avg_y - ay))
                            if len(area) == 0:
                                a = ro; continue
                            idx = np.argmax(area); a = ro + idx
                            xds.append(x[a]); yds.append(y[a])
                        xds.append(x[-1]); yds.append(y[-1])
                        new[x_col] = np.asarray(xds); new[c] = np.asarray(yds)
            if new.get(x_col) is not None:
                df = pd.DataFrame(new)

        fig = go.Figure()
        for y in y_left:
            if y in df:
                fig.add_trace(go.Scattergl(x=df[x_col], y=df[y], mode="lines", name=y, yaxis="y1"))
        for col in df.columns:
            if col.startswith("f:"):
                fig.add_trace(go.Scattergl(x=df[x_col], y=df[col], mode="lines", name=col, yaxis="y1"))
        if y_right != "(none)" and y_right in df:
            fig.add_trace(go.Scattergl(x=df[x_col], y=df[y_right], mode="lines", name=y_right, yaxis="y2"))
            fig.update_layout(yaxis2=dict(title=y_right, overlaying="y", side="right"))
        fig.update_layout(xaxis_title=x_col, yaxis_title="Value", margin=dict(l=10,r=10,t=30,b=40), height=520, legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0))
        st.plotly_chart(fig, use_container_width=True)

        num_cols = [c for c in df.columns if c != x_col]
        if num_cols:
            stats = df[num_cols].describe().transpose()[["mean","std","min","max"]]
            stats.columns = ["mean","std","min","max"]
            st.caption("Quick stats (after smoothing/downsampling):")
            st.dataframe(stats)

        c1, c2, c3 = st.columns(3)
        c1.download_button("Download CSV (all plotted columns)", df.to_csv(index=False).encode("utf-8"), file_name=f"{pf.name}_run{run_num}_selection.csv", mime="text/csv")
        y_quick = c2.selectbox("Two-column export: Y", [c for c in df.columns if c != x_col])
        c2.download_button("Download X,Y CSV", df[[x_col, y_quick]].to_csv(index=False).encode("utf-8"),
                           file_name=f"{pf.name}_run{run_num}_{x_col}_{y_quick}.csv", mime="text/csv")
        try:
            png_bytes = pio.to_image(fig, format="png", width=1400, height=600, scale=2)
            c3.download_button("Download PNG", data=png_bytes, file_name=f"{pf.name}_plot.png", mime="image/png")
        except Exception:
            st.info("PNG export requires `kaleido` (pip install kaleido).")
        html_bytes = pio.to_html(fig, include_plotlyjs="cdn").encode("utf-8")
        st.download_button("Download HTML (interactive)", data=html_bytes, file_name=f"{pf.name}_plot.html", mime="text/html")

        st.subheader("Data preview")
        st.dataframe(df.head(1000))

    else:
        choices = st.multiselect("Scopes (file/run)", options=list(range(len(scopes))), format_func=lambda i: scopes[i].display)
        if not choices:
            st.info("Select scopes above."); st.stop()
        common_cols = None
        for idx in choices:
            pf = files[scopes[idx].file_key]
            cols = set(get_keywords(pf, scopes[idx].run_num))
            common_cols = cols if common_cols is None else (common_cols & cols)
        common_cols = sorted(list(common_cols)) if common_cols else []
        if not common_cols:
            st.error("No common columns across selected scopes."); st.stop()
        cox, coy = st.columns([1,1])
        x_default = "Time" if "Time" in common_cols else ("Step" if "Step" in common_cols else common_cols[0])
        x_col = cox.selectbox("X", common_cols, index=common_cols.index(x_default))
        y_col = coy.selectbox("Y", [c for c in common_cols if c != x_col])
        dpolicy = st.selectbox("Dup. timestep policy (All runs)", ["keep_all","keep_first","keep_last"])

        fig = go.Figure()
        long_frames = []
        for idx in choices:
            sc = scopes[idx]; pf = files[sc.file_key]
            x = get_series(pf, x_col, sc.run_num); y = get_series(pf, y_col, sc.run_num)
            if x is None or y is None: st.warning(f"Skipping {sc.display} — missing data."); continue
            df = pd.DataFrame({x_col: x, y_col: y}).dropna()
            if sc.run_num == -1 and dpolicy != "keep_all":
                df = df.drop_duplicates(subset=[x_col], keep=("first" if dpolicy=="keep_first" else "last"))
            if avg_N > 1:
                df[y_col] = pd.Series(df[y_col]).rolling(window=int(avg_N), min_periods=1, center=True).mean()
            if down_method == "Every Nth" and ds_every > 1:
                df = df.iloc[::int(ds_every), :].reset_index(drop=True)
            elif down_method == "LTTB" and lttb_pts > 0 and len(df) > 3:
                xx = df[x_col].to_numpy(); yy = df[y_col].to_numpy(); n = len(xx); T = int(lttb_pts)
                if 3 <= T < n:
                    every = (n - 2) / (T - 2); a = 0; xds = [xx[0]]; yds = [yy[0]]
                    for i in range(0, T-2):
                        ars = int(np.floor((i + 1) * every) + 1); are = int(np.floor((i + 2) * every) + 1); are = min(are, n)
                        avg_x = np.mean(xx[ars:are]) if are > ars else xx[ars-1]
                        avg_y = np.mean(yy[ars:are]) if are > ars else yy[ars-1]
                        ro = int(np.floor(i * every) + 1); rt = int(np.floor((i + 1) * every) + 1); rt = min(rt, n-1)
                        bx = xx[ro:rt]; by = yy[ro:rt]; ax = xx[a]; ay = yy[a]
                        area = np.abs((ax - avg_x) * (by - ay) - (ax - bx) * (avg_y - ay))
                        if len(area) == 0: a = ro; continue
                        idxmax = np.argmax(area); a = ro + idxmax; xds.append(xx[a]); yds.append(yy[a])
                    xds.append(xx[-1]); yds.append(yy[-1]); df = pd.DataFrame({x_col: np.asarray(xds), y_col: np.asarray(yds)})
            fig.add_trace(go.Scattergl(x=df[x_col], y=df[y_col], mode="lines", name=sc.display))
            df["scope"] = sc.display
            long_frames.append(df)
        fig.update_layout(xaxis_title=x_col, yaxis_title=y_col, margin=dict(l=10,r=10,t=30,b=40), height=520, legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0))
        st.plotly_chart(fig, use_container_width=True)

        if long_frames:
            out = pd.concat(long_frames, ignore_index=True)
            st.download_button("Download CSV (overlay long format)", out.to_csv(index=False).encode("utf-8"),
                               file_name=f"overlay_{x_col}_{y_col}.csv", mime="text/csv")

# ================= Errors & Warnings ==================
with tab_errs:
    st.subheader("Errors & Warnings")
    rows = []
    for h, pf in files.items():
        for i, line in enumerate(pf.lines):
            if WARN_RE.match(line) or ERR_RE.match(line):
                rows.append({"file": pf.name, "line_no": i+1, "type": "ERROR" if ERR_RE.match(line) else "WARNING", "text": line.strip()})
    if rows:
        df = pd.DataFrame(rows)
        st.dataframe(df, use_container_width=True, height=350)
        st.write("Counts:", df.groupby(["file","type"]).size().rename("count").reset_index())
    else:
        st.info("No WARNING/ERROR lines found.")

# ================= Performance ==================
with tab_perf:
    st.subheader("Performance summary (per run)")
    for h, pf in files.items():
        st.markdown(f"### {pf.name}")
        for ridx, (hdr, s, e) in enumerate(pf.runs):
            endwin = min(e+80, len(pf.lines))
            block = pf.lines[s:endwin]
            loop_line = next((ln for ln in block if LOOP_RE.match(ln)), None)
            if not loop_line: continue
            m = re.search(r'Loop time of\s+([+\-.\deE]+)\s+on\s+(\d+)\s+procs(?:\s+for\s+(\d+)\s+steps)?', loop_line, re.IGNORECASE)
            if not m: continue
            loop_time = float(m.group(1)); procs = int(m.group(2))
            steps = int(m.group(3)) if m.group(3) else None
            if steps is None:
                x = get_series(pf, "Step", ridx)
                if x is not None and len(x) > 0:
                    steps = int(x[-1] - x[0]) if len(x) > 1 else int(x[0])
            steps_s = (steps / loop_time) if (steps and loop_time > 0) else None
            st.write(f"**Run #{ridx+1}** — loop time: {loop_time:.4g} s on {procs} procs; steps: {steps if steps is not None else '—'}; steps/s: {steps_s:.3g}" )
            cats_old = {}
            for ln in block:
                m2 = BREAKDOWN_RE_OLD.match(ln)
                if m2:
                    cats_old[m2.group(1).capitalize()] = float(m2.group(2))
            cats_tbl = {}
            for j, ln in enumerate(block):
                if BREAKDOWN_TABLE_HEADER.search(ln):
                    k = j + 3  # skip header lines
                    while k < len(block):
                        row = block[k]
                        if not row.strip(): break
                        mrow = BREAKDOWN_TABLE_ROW.match(row)
                        if mrow:
                            name = mrow.group(1).capitalize(); avg_time = mrow.group(3)
                            try: cats_tbl[name] = float(avg_time)
                            except: pass
                        k += 1
                    break
            cats = cats_tbl if cats_tbl else cats_old
            if cats:
                order = ["Pair","Neigh","Comm","Output","Modify","Other","Bond"]
                xs = [k for k in order if k in cats] + [k for k in cats.keys() if k not in order]
                ys = [cats[k] for k in xs]
                fig = go.Figure([go.Bar(x=xs, y=ys)])
                fig.update_layout(height=260, margin=dict(l=10,r=10,t=10,b=30), yaxis_title="avg time (s)")
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.caption("No timing breakdown found for this run.")

# ================= Prefix Search ==================
with tab_search:
    st.subheader("Prefix search")
    defaults = ["fix", "variable", "print", "pair_style", "bond_style", "kspace_style", "compute", "thermo", "run", "minimize"]
    cols = st.columns(2)
    custom = cols[0].text_input("Additional prefixes (comma-separated)", value="")
    use_def = cols[1].toggle("Include defaults", value=True)
    prefixes = [p.strip() for p in custom.split(",") if p.strip()]
    if use_def: prefixes = defaults + prefixes
    if not prefixes:
        st.info("Enter at least one prefix.")
    else:
        pattern = re.compile(rf'^\s*(?:{"|".join([re.escape(p) for p in prefixes])})\b', re.IGNORECASE)
        hits = []
        for h, pf in files.items():
            for i, line in enumerate(pf.lines):
                if pattern.match(line):
                    hits.append({"file": pf.name, "line_no": i+1, "text": line.rstrip()})
        if hits:
            outdf = pd.DataFrame(hits)
            st.dataframe(outdf, use_container_width=True, height=400)
            st.download_button("Download matches (CSV)", outdf.to_csv(index=False).encode("utf-8"),
                               file_name="prefix_matches.csv", mime="text/csv")
        else:
            st.info("No matches found for the chosen prefixes.")

# ================= View presets ==================
with tab_presets:
    st.subheader("Save current view")
    snapshot = {"formulas": st.session_state.formulas}
    st.download_button("Download preset (JSON)", json.dumps(snapshot, indent=2).encode("utf-8"),
                       file_name="lmp_loglite_preset.json", mime="application/json")
    st.subheader("Load preset")
    jp = st.file_uploader("Upload preset JSON", type=["json"], key="preset_uploader")
    if jp is not None:
        try:
            cfg = json.loads(jp.read().decode("utf-8"))
            st.session_state.formulas = cfg.get("formulas", [])
            st.success("Preset loaded (reselect scopes if needed).")
        except Exception as e:
            st.error(f"Invalid preset: {e}")

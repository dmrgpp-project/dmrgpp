#!/usr/bin/env python3
"""Create a dependency-free SVG comparison of impurity Green functions."""
from __future__ import annotations

import argparse
import html
import math
import re
from pathlib import Path

COLORS = {"1.5": "#1769aa", "2.5": "#d1495b", "3.5": "#2a9d55"}


def table(path: Path) -> list[tuple[float, float, float]]:
    rows = []
    if not path.is_file():
        raise ValueError(f"missing data file: {path}")
    for n, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"{path}:{n}: expected exactly three columns")
        try:
            row = tuple(float(v) for v in fields)
        except ValueError as exc:
            raise ValueError(f"{path}:{n}: invalid numeric value") from exc
        if not all(math.isfinite(v) for v in row):
            raise ValueError(f"{path}:{n}: values must be finite")
        rows.append(row)
    if not rows:
        raise ValueError(f"no data rows in {path}")
    return rows


def input_value(text: str, key: str) -> str:
    m = re.search(rf"^\s*(?:\w+\s+)?{re.escape(key)}\s*=\s*([^;]+);", text, re.M)
    if not m:
        raise ValueError(f"input.ain does not define {key}")
    return m.group(1).strip().strip('"')


def nice_ticks(lo: float, hi: float, count: int = 6) -> list[float]:
    if lo == hi:
        lo, hi = lo - 1, hi + 1
    raw = (hi - lo) / (count - 1)
    power = 10 ** math.floor(math.log10(raw))
    step = min((1, 2, 2.5, 5, 10), key=lambda m: abs(m * power - raw)) * power
    start = math.ceil(lo / step) * step
    ticks = []
    v = start
    while v <= hi + step * 1e-9 and len(ticks) < 30:
        ticks.append(v)
        v += step
    return ticks


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("runs", nargs="+", type=Path, help="completed run directories")
    ap.add_argument("--output", required=True, type=Path, help="destination SVG path")
    args = ap.parse_args()
    try:
        runs = []
        for directory in args.runs:
            inp = (directory / "input.ain").read_text()
            u = input_value(inp, "HubbardU")
            eta = input_value(inp, "OmegaDelta")
            solver_match = re.search(r"^\s*ImpuritySolver\s*=\s*([^;]+);", inp, re.M)
            solver = solver_match.group(1).strip().strip('"').lower() if solver_match else "dmrg"
            mats = [r for r in table(directory / f"gimp_{solver}.txt") if r[0] > 0]
            real = table(directory / f"gimp_{solver}_real.txt")
            if not mats:
                raise ValueError(f"{directory}: no positive Matsubara frequencies")
            logpath = directory / "cincuenta.log"
            log = logpath.read_text(errors="replace") if logpath.exists() else ""
            runs.append((directory, u, eta, mats, real, "Converged after" in log, solver))
        if len({r[2] for r in runs}) != 1:
            raise ValueError("runs use different OmegaDelta broadenings")
        eta = runs[0][2]
        # Fixed, readable two-panel canvas; shared vertical scale aids comparison.
        W, H = 1280, 600
        panels = [(90, 82, 510, 385), (690, 82, 510, 385)]
        all_y = [-row[2] for r in runs for rows in (r[3], r[4]) for row in rows]
        ymin, ymax = min(0.0, min(all_y)), max(0.0, max(all_y))
        pad = max((ymax - ymin) * .07, 1e-12)
        ymin -= pad
        ymax += pad
        out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
               '<rect width="100%" height="100%" fill="white"/>',
               '<style>text{font-family:Arial,sans-serif;fill:#222}.grid{stroke:#ddd;stroke-width:1}.axis{stroke:#333;stroke-width:1.4}.curve{fill:none;stroke-width:2.4}</style>',
               '<text x="640" y="30" text-anchor="middle" font-size="21" font-weight="bold">Impurity Green functions</text>',
               f'<text x="640" y="53" text-anchor="middle" font-size="14">CPU impurity solvers · real-axis broadening η={html.escape(eta)}</text>']
        datasets = [([r[3] for r in runs], "Matsubara frequency ωₙ", "-Im G(iωₙ)"),
                    ([r[4] for r in runs], "Real frequency ω", f"-Im G(ω+iη), η={eta}")]
        for pi, (x0, y0, pw, ph) in enumerate(panels):
            group = datasets[pi][0]
            xs = [row[0] for seq in group for row in seq]
            xmin, xmax = min(xs), max(xs)
            if pi == 0:
                xmin = 0.0
            if xmin == xmax: xmax = xmin + 1
            X = lambda v: x0 + (v-xmin)/(xmax-xmin)*pw
            Y = lambda v: y0 + ph - (v-ymin)/(ymax-ymin)*ph
            out.append(f'<text x="{x0+pw/2}" y="76" text-anchor="middle" font-size="17" font-weight="bold">{"Positive Matsubara axis" if pi==0 else "Broadened real axis"}</text>')
            for t in nice_ticks(xmin, xmax):
                xx=X(t); out.append(f'<line class="grid" x1="{xx:.2f}" y1="{y0}" x2="{xx:.2f}" y2="{y0+ph}"/><text x="{xx:.2f}" y="{y0+ph+21}" text-anchor="middle" font-size="12">{t:.4g}</text>')
            for t in nice_ticks(ymin, ymax):
                yy=Y(t); out.append(f'<line class="grid" x1="{x0}" y1="{yy:.2f}" x2="{x0+pw}" y2="{yy:.2f}"/><text x="{x0-9}" y="{yy+4:.2f}" text-anchor="end" font-size="12">{t:.4g}</text>')
            if ymin <= 0 <= ymax:
                yy=Y(0); out.append(f'<line class="axis" x1="{x0}" y1="{yy:.2f}" x2="{x0+pw}" y2="{yy:.2f}"/>')
            out.append(f'<line class="axis" x1="{x0}" y1="{y0+ph}" x2="{x0+pw}" y2="{y0+ph}"/><line class="axis" x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0+ph}"/>')
            out.append(f'<text x="{x0+pw/2}" y="{y0+ph+48}" text-anchor="middle" font-size="14">{datasets[pi][1]}</text>')
            out.append(f'<text transform="translate({x0-58} {y0+ph/2}) rotate(-90)" text-anchor="middle" font-size="14">{datasets[pi][2]}</text>')
            for ri, seq in enumerate(group):
                _, u, _, _, _, conv, solver = runs[ri]
                color=COLORS.get(u, "#7950a1")
                points=" ".join(f"{X(a):.2f},{Y(-c):.2f}" for a,_,c in seq)
                is_ed = solver in ("ed", "exactdiag")
                dash = ' stroke-dasharray="8 5"' if is_ed else ""
                if not conv:
                    dash = (' stroke-dasharray="2 4"' if is_ed
                            else ' stroke-dasharray="8 3 2 3"')
                out.append(f'<polyline class="curve" stroke="{color}"{dash} points="{points}"/>')
        # Shared two-row, three-column legend below panels.
        lx, ly = 250, 550
        for i, r in enumerate(runs):
            _, u, _, _, _, conv, solver = r
            xx=lx+(i%3)*300; yy=ly+(i//3)*28; color=COLORS.get(u,"#7950a1")
            is_ed = solver in ("ed", "exactdiag")
            dash = ' stroke-dasharray="8 5"' if is_ed else ""
            if not conv:
                dash = (' stroke-dasharray="2 4"' if is_ed
                        else ' stroke-dasharray="8 3 2 3"')
            label = f"{'ED' if is_ed else solver.upper()} U={u}" + (" (not converged)" if not conv else "")
            out.append(f'<line x1="{xx}" y1="{yy-5}" x2="{xx+34}" y2="{yy-5}" stroke="{color}" stroke-width="3"{dash}/><text x="{xx+42}" y="{yy}" font-size="14">{html.escape(label)}</text>')
        out.append('</svg>')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text("\n".join(out)+"\n")
    except (OSError, ValueError) as exc:
        ap.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

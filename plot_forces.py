import argparse
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from cp_integration import (alpha_from_name, integrate, naca4, read_foam_wall_cp,
                            read_xflr5_cp, xflr5_contour_cp)


def rotate(pts, alpha_deg):
    al = np.radians(alpha_deg)
    d = np.atleast_2d(np.asarray(pts, float)) - np.array([1.0, 0.0])
    x = d[:, 0] * np.cos(al) + d[:, 1] * np.sin(al)
    y = -d[:, 0] * np.sin(al) + d[:, 1] * np.cos(al)
    return np.column_stack((1.0 + x, y))


def load_results(args):
    res = {"OpenFOAM": {}, "XFLR5": {}}
    for case in glob.glob(args.cases):
        m = re.search(r"af_a(-?\d+(?:\.\d+)?)$", case)
        if not m:
            continue
        a = float(m.group(1))
        try:
            r, cp_e, _ = read_foam_wall_cp(case, args.patch, args.umag)
        except Exception as e:
            print(f"[OpenFOAM alpha={a}] bỏ qua: {e}")
            continue
        res["OpenFOAM"][a] = integrate(r, cp_e, a, cp_on="edges")
    for f in glob.glob(os.path.join(args.xflr5_dir, "*")):
        if os.path.isdir(f):
            continue
        a = alpha_from_name(f)
        if a is None:
            continue
        x, cp = read_xflr5_cp(f, args.xflr5_col)
        r, cpn = xflr5_contour_cp(x, cp, args.code, a)
        res["XFLR5"][a] = integrate(r, cpn, a, cp_on="nodes")
    return res


def cp_point(o, a):
    """Vị trí tâm áp suất trên dây cung; nếu không xác định thì dùng c/4."""
    xcp = o["Xcp"]
    ok = np.isfinite(xcp)
    return (xcp if ok else 0.25), ok


def draw_panel(ax, code, a, o, name, scale, gain, ylim):
    upper, lower = naca4(code)
    ru, rl = rotate(upper, a), rotate(lower, a)
    le, te = rotate([[0.0, 0.0]], a)[0], np.array([1.0, 0.0])
    xcp, ok = cp_point(o, a)
    p = rotate([[xcp, 0.0]], a)[0]

    ax.plot(ru[:, 0], ru[:, 1], "b-", lw=2)
    ax.plot(rl[:, 0], rl[:, 1], "r-", lw=2)
    ax.plot([le[0], te[0]], [le[1], te[1]], "k--", lw=1.5)
    ax.plot(*le, "ko", ms=6)
    ax.plot(*te, "ko", ms=6)
    ax.annotate("LE", le, xytext=(6, 8), textcoords="offset points", weight="bold")
    ax.annotate("TE", te, xytext=(6, -3), textcoords="offset points", weight="bold")
    ax.axhline(0, color="gray", lw=0.6, ls="--")

    ax.axvline(p[0], color="magenta", ls=":", lw=1.5)
    ax.plot(*p, "o", color="magenta", ms=10, zorder=5)
    hw = 0.025
    ax.arrow(p[0], p[1], 0, scale * o["CL"], width=0.006, head_width=hw, head_length=hw,
             length_includes_head=True, color="green", zorder=4)
    ax.arrow(p[0], p[1], gain * scale * o["CDp"], 0, width=0.006, head_width=hw, head_length=hw,
             length_includes_head=True, color="orange", zorder=4)

    ax.set_title(f"{name}", fontsize=13)
    ax.set_xlabel("Toạ độ x")
    ax.set_ylabel("Toạ độ y")
    ax.set_xlim(-0.05, 1.15)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(ls=":", alpha=0.6)

    xcp_label = f"$x_{{cp}}$ = {xcp:.3f}" if ok else "$x_{cp}$ không xác định (vẽ tại c/4)"
    handles = [Line2D([0], [0], color="b", lw=2, label="Upper Surface"),
               Line2D([0], [0], color="r", lw=2, label="Lower Surface"),
               Line2D([0], [0], color="k", ls="--", label="Chord Line"),
               Line2D([0], [0], marker="o", color="w", markerfacecolor="magenta", ms=10,
                      label=f"Center of Pressure ({xcp_label})"),
               Patch(color="green", label=f"Lift (L), CL = {o['CL']:.4f}"),
               Patch(color="orange", label=f"Drag áp suất (×{gain:g}), CDp = {o['CDp']:.4f}")]
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.95)
    ax.text(0.02, 0.03, f"$C_m$(c/4) = {o['Cm_c4']:.4f}", transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--code", default="2412", help="mã NACA 4 chữ số")
    ap.add_argument("--cases", default=os.path.join(here, "cases", "af_a*"))
    ap.add_argument("--xflr5-dir", default=os.path.join(here, "xflr5_cp"))
    ap.add_argument("--xflr5-col", type=int, default=None)
    ap.add_argument("--umag", type=float, default=26.0)
    ap.add_argument("--patch", default="walls")
    ap.add_argument("--outdir", default=os.path.join(here, "hinh_luc"))
    ap.add_argument("--scale", type=float, default=0.3, help="độ dài mũi tên lift = scale * CL (đơn vị chord)")
    ap.add_argument("--drag-gain", type=float, default=10.0, help="hệ số phóng đại mũi tên drag so với lift")
    args = ap.parse_args()

    res = load_results(args)
    sources = [s for s in res if res[s]]
    angles = sorted(set().union(*[set(res[s]) for s in sources])) if sources else []
    if not angles:
        raise SystemExit("Không đọc được dữ liệu OpenFOAM hoặc XFLR5.")
    os.makedirs(args.outdir, exist_ok=True)

    # cùng thang trục y cho mọi hình để so sánh được
    ys = [0.0, 0.25]
    for s in sources:
        for a, o in res[s].items():
            py = rotate([[cp_point(o, a)[0], 0.0]], a)[0][1]
            ys += [py + args.scale * o["CL"], py, rotate([[0.0, 0.0]], a)[0][1] + 0.12]
    ylim = (min(-0.22, min(ys) - 0.05), max(ys) + 0.08)

    for a in angles:
        srcs = [s for s in sources if a in res[s]]
        fig, axs = plt.subplots(1, len(srcs), figsize=(7.5 * len(srcs), 5.0), squeeze=False)
        for ax, s in zip(axs[0], srcs):
            draw_panel(ax, args.code, a, res[s][a], s, args.scale, args.drag_gain, ylim)
        fig.suptitle(f"NACA {args.code}: hình học, tâm áp suất và lực khí động ($\\alpha$ = {a:g}°)", fontsize=14)
        fig.tight_layout()
        path = os.path.join(args.outdir, f"forces_a{a:g}.png")
        fig.savefig(path, dpi=200)
        plt.close(fig)
        print("Đã lưu", path)


if __name__ == "__main__":
    main()

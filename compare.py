import argparse
import csv
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cp_integration import (aerodynamic_center, alpha_from_name, integrate, read_foam_wall_cp,
                            read_xflr5_cp, xflr5_contour_cp)

here = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--code", default="2412", help="mã NACA 4 chữ số")
ap.add_argument("--cases", default=os.path.join(here, "cases", "af_a*"), help="các case OpenFOAM theo góc tấn")
ap.add_argument("--xflr5-dir", default=os.path.join(here, "xflr5_cp"), help="thư mục file Cp XFLR5 (a0.txt, a5.txt, ...)")
ap.add_argument("--xflr5-col", type=int, default=None, help="chỉ số cột Cp (mặc định: cột Cpv)")
ap.add_argument("--umag", type=float, default=26.0, help="|U| của OpenFOAM (m/s)")
ap.add_argument("--patch", default="walls")
args = ap.parse_args()

results = {"OpenFOAM": {}, "XFLR5": {}}
cp_data = {"OpenFOAM": {}, "XFLR5": {}}
check = {}

for case in glob.glob(args.cases):
    m = re.search(r"af_a(-?\d+(?:\.\d+)?)$", case)
    if not m:
        continue
    a = float(m.group(1))
    try:
        r, cp_e, cen = read_foam_wall_cp(case, args.patch, args.umag)
    except Exception as e:
        print(f"[OpenFOAM alpha={a}] bỏ qua: {e}")
        continue
    results["OpenFOAM"][a] = integrate(r, cp_e, a, cp_on="edges")
    cp_data["OpenFOAM"][a] = (cen[:, 0], cp_e)
    fc = sorted(glob.glob(os.path.join(case, "postProcessing", "forceCoeffsDict", "*", "forceCoeffs.dat")))
    if fc:
        last = [l.split() for l in open(fc[-1]) if l.strip() and not l.startswith("#")][-1]
        check[a] = (-float(last[1]), float(last[3]), float(last[2]))  # Cm (đổi dấu), CL, CD từ forceCoeffs

for f in glob.glob(os.path.join(args.xflr5_dir, "*")):
    a = alpha_from_name(f)
    if a is None or os.path.isdir(f):
        continue
    x, cp = read_xflr5_cp(f, args.xflr5_col)
    r, cpn = xflr5_contour_cp(x, cp, args.code, a)
    results["XFLR5"][a] = integrate(r, cpn, a, cp_on="nodes")
    cp_data["XFLR5"][a] = (x, cp)

# ---- bảng kết quả
rows = []
print(f"{'nguồn':9s} {'alpha':>6s} {'CL':>8s} {'CDp':>8s} {'Cm_c4':>8s} {'Xcp/c':>7s}")
for src in results:
    for a in sorted(results[src]):
        o = results[src][a]
        rows.append([src, a, o["CL"], o["CDp"], o["Cm_c4"], o["Xcp"], o["Cy"]])
        print(f"{src:9s} {a:6.1f} {o['CL']:8.4f} {o['CDp']:8.4f} {o['Cm_c4']:8.4f} {o['Xcp']:7.3f}")
with open("ket_qua_tich_phan.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["source", "alpha", "CL", "CDp", "Cm_c4", "Xcp", "Cn"])
    w.writerows(rows)

xac = {}
for src in results:
    al = sorted(results[src])
    if len(al) >= 2:
        xac[src] = aerodynamic_center(al, [results[src][a]["Cy"] for a in al], [results[src][a]["Cm_c4"] for a in al])
        print(f"Xac/c {src} = {xac[src]:.3f}  (hồi quy trên alpha <= 10 deg)")

if check:
    print("\nĐối chiếu OpenFOAM: tích phân Cp so với forceCoeffs (forceCoeffs gồm cả ma sát, CD sẽ lớn hơn CDp)")
    for a in sorted(check):
        o = results["OpenFOAM"][a]
        cm, cl, cd = check[a]
        print(f"  alpha={a:5.1f}  CL {o['CL']:.4f} vs {cl:.4f} | Cm_c4 {o['Cm_c4']:.4f} vs {cm:.4f} | CDp {o['CDp']:.4f} vs CD {cd:.4f}")

# ---- đồ thị hệ số
style = {"OpenFOAM": "o--", "XFLR5": "s-"}
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
for src in results:
    al = sorted(results[src])
    if not al:
        continue
    g = lambda k: [results[src][a][k] for a in al]
    ax[0, 0].plot(al, g("CL"), style[src], label=src)
    ax[0, 1].plot(al, g("Cm_c4"), style[src], label=src)
    ax[1, 0].plot(al, g("Xcp"), style[src], label=src)
    ax[1, 1].plot(al, g("CDp"), style[src], label=src)
for a_, (xl, yl, t) in zip(ax.ravel(), [("alpha (deg)", "CL", "CL - alpha"), ("alpha (deg)", "Cm (c/4, ngóc mũi dương)", "Cm - alpha"),
                                         ("alpha (deg)", "Xcp/c", "Tâm áp suất"), ("alpha (deg)", "CD do áp suất", "Lực cản áp suất")]):
    a_.set(xlabel=xl, ylabel=yl, title=t); a_.grid(alpha=0.3); a_.legend()
plt.tight_layout(); plt.savefig("so_sanh.png", dpi=200)

# ---- Cp(x/c) chồng hai nguồn
angles = sorted(set(cp_data["OpenFOAM"]) | set(cp_data["XFLR5"]))
if angles:
    nc = min(2, len(angles)); nr = int(np.ceil(len(angles) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(6 * nc, 4.2 * nr), squeeze=False)
    for a_, a in zip(axs.ravel(), angles):
        if a in cp_data["XFLR5"]:
            a_.plot(*cp_data["XFLR5"][a], "-", label="XFLR5")
        if a in cp_data["OpenFOAM"]:
            a_.plot(*cp_data["OpenFOAM"][a], ".", ms=3, label="OpenFOAM")
        a_.invert_yaxis(); a_.grid(alpha=0.3); a_.legend()
        a_.set(title=f"alpha = {a:g} deg", xlabel="x/c", ylabel="Cp")
    plt.tight_layout(); plt.savefig("cp_so_sanh.png", dpi=200)
print("\nĐã lưu so_sanh.png, cp_so_sanh.png, ket_qua_tich_phan.csv")

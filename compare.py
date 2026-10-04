"""So sánh OpenFOAM (định tính) và XFLR5 (định lượng): CL, Cm, Xcp, Xac theo alpha.

Cách dùng:
    python3 compare.py --xflr5 polar_xflr5.txt
(XFLR5: xuất polar bằng File > Export (Polars) dạng .txt, Cm lấy quanh c/4)
"""
import argparse, glob, os, re
import numpy as np
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--xflr5", required=True, help="file polar XFLR5 (.txt)")
ap.add_argument("--cases", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "cases", "af_a*"), help="mẫu đường dẫn các case")
args = ap.parse_args()


def read_foam(case):
    """Đọc dòng cuối của forceCoeffs.dat -> (Cm, Cd, Cl)."""
    files = sorted(glob.glob(f"{case}/postProcessing/forceCoeffsDict/*/forceCoeffs.dat"))
    if not files:
        return None
    header, last = None, None
    for line in open(files[-1]):
        if line.startswith("#"):
            if "Cm" in line and "Cl" in line:
                header = line.lstrip("#").split()
        elif line.strip():
            last = line.split()
    col = {n: i for i, n in enumerate(header)}
    return float(last[col["Cm"]]), float(last[col["Cd"]]), float(last[col["Cl"]])


def read_xflr5(path):
    """Cột polar XFLR5: alpha CL CD CDp Cm ..."""
    rows = []
    for line in open(path, errors="ignore"):
        p = line.split()
        if len(p) >= 5 and re.fullmatch(r"-?\d+\.?\d*", p[0]):
            try:
                rows.append([float(x) for x in p[:5]])
            except ValueError:
                pass
    return np.array(rows)  # alpha, CL, CD, CDp, Cm


# --- OpenFOAM ---
alphas, foam = [], []
for c in glob.glob(args.cases):
    m = re.search(r"af_a(-?\d+)$", c)
    r = read_foam(c) if m else None
    if r:
        alphas.append(int(m.group(1))); foam.append(r)
order = np.argsort(alphas)
a_f = np.array(alphas)[order]
foam = np.array(foam)[order]
Cm_f, Cd_f, Cl_f = foam[:, 0], foam[:, 1], foam[:, 2]

# --- XFLR5 ---
x = read_xflr5(args.xflr5)
a_x, Cl_x, Cd_x, Cm_x = x[:, 0], x[:, 1], x[:, 2], x[:, 4]

# --- Xcp và Xac (moment quanh c/4): Xcp/c = 0.25 - Cm/CL ; Xac/c = 0.25 - dCm/dCL ---
def xcp(Cm, Cl):
    Cl = np.where(np.abs(Cl) < 0.05, np.nan, Cl)  # tránh chia cho CL ~ 0
    return 0.25 - Cm / Cl

def xac(Cm, Cl, mask):
    return 0.25 - np.polyfit(Cl[mask], Cm[mask], 1)[0]

lin_f = a_f <= 5            # vùng tuyến tính để tính tâm khí động
lin_x = a_x <= 5
print("alpha_OpenFOAM:", a_f)
print(f"Xac/c OpenFOAM = {xac(Cm_f, Cl_f, lin_f):.3f}" if lin_f.sum() >= 2 else "Cần >= 2 góc <= 5 deg để tính Xac")
print(f"Xac/c XFLR5    = {xac(Cm_x, Cl_x, lin_x):.3f}")

# --- Vẽ ---
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
ax[0, 0].plot(a_x, Cl_x, "-", label="XFLR5"); ax[0, 0].plot(a_f, Cl_f, "o--", label="OpenFOAM")
ax[0, 0].set(xlabel="alpha (deg)", ylabel="CL", title="CL - alpha")
ax[0, 1].plot(a_x, Cm_x, "-", label="XFLR5"); ax[0, 1].plot(a_f, Cm_f, "o--", label="OpenFOAM")
ax[0, 1].set(xlabel="alpha (deg)", ylabel="Cm (c/4)", title="Cm - alpha")
ax[1, 0].plot(a_x, xcp(Cm_x, Cl_x), "-", label="XFLR5"); ax[1, 0].plot(a_f, xcp(Cm_f, Cl_f), "o--", label="OpenFOAM")
ax[1, 0].set(xlabel="alpha (deg)", ylabel="Xcp/c", title="Tam ap suat")
ax[1, 1].plot(Cd_x, Cl_x, "-", label="XFLR5"); ax[1, 1].plot(Cd_f, Cl_f, "o--", label="OpenFOAM")
ax[1, 1].set(xlabel="CD", ylabel="CL", title="Drag polar")
for a in ax.ravel():
    a.grid(True, alpha=0.3); a.legend()
plt.tight_layout(); plt.savefig("so_sanh.png", dpi=200)
print("Đã lưu so_sanh.png")

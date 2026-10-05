import glob
import os
import re

import numpy as np

NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


# ---------------------------------------------------------------- hình học
def naca4(code, n=200):

    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n)))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
    if m == 0:
        yc, dyc = np.zeros_like(x), np.zeros_like(x)
    else:
        yc = np.where(x < p, m / p**2 * (2 * p * x - x**2), m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x**2))
        dyc = np.where(x < p, 2 * m / p**2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    th = np.arctan(dyc)
    upper = np.column_stack((x - yt * np.sin(th), yc + yt * np.cos(th)))
    lower = np.column_stack((x + yt * np.sin(th), yc - yt * np.cos(th)))
    return upper, lower


def closed_contour(upper, lower):
    """Contour kín không lặp điểm: mặt trên mũi->đuôi, rồi mặt dưới đuôi->mũi."""
    return np.vstack((upper, lower[-2:0:-1]))


# ---------------------------------------------------------------- tích phân
def integrate(r, cp, alpha_deg, cp_on="nodes"):
  
    r = np.asarray(r, float)
    cp = np.asarray(cp, float)
    a, b = r, np.roll(r, -1, axis=0)
    ca = cp
    cb = np.roll(cp, -1) if cp_on == "nodes" else cp
    edge = b - a
    ds = np.hypot(edge[:, 0], edge[:, 1])
    area2 = np.sum(a[:, 0] * b[:, 1] - b[:, 0] * a[:, 1])
    normal = np.sign(area2) * np.column_stack((edge[:, 1], -edge[:, 0])) / ds[:, None]
    force = np.zeros(2)
    moment = 0.0  # mô-men z quanh mũi (0,0), chiều ngược kim đồng hồ dương
    for t in (0.5 - 0.5 / np.sqrt(3), 0.5 + 0.5 / np.sqrt(3)):
        rq = (1 - t) * a + t * b
        d = -0.5 * ((1 - t) * ca + t * cb)[:, None] * normal * ds[:, None]
        force += d.sum(axis=0)
        moment += np.sum(rq[:, 0] * d[:, 1] - rq[:, 1] * d[:, 0])
    cx, cy = force
    al = np.radians(alpha_deg)
    cd = np.cos(al) * cx + np.sin(al) * cy
    cl = -np.sin(al) * cx + np.cos(al) * cy
    out = dict(Cx=cx, Cy=cy, CL=cl, CDp=cd, MzLE=moment)
    out["Xcp"] = moment / cy if abs(cy) > 0.05 else np.nan
    # Cm quy ước hàng không (ngóc mũi dương) quanh c/4: ngược dấu mô-men z
    out["Cm_c4"] = -(moment - 0.25 * cy)
    return out


def aerodynamic_center(alpha, cn, cm_c4, amax=10.0):
    """Xac/c = 0.25 - (dCm/dalpha)/(dCn/dalpha), hồi quy tuyến tính trên alpha <= amax."""
    alpha, cn, cm_c4 = map(np.asarray, (alpha, cn, cm_c4))
    k = alpha <= amax
    if k.sum() < 2:
        return np.nan
    return 0.25 - np.polyfit(alpha[k], cm_c4[k], 1)[0] / np.polyfit(alpha[k], cn[k], 1)[0]


 #XFLR5
def read_xflr5_cp(path, col=None):
    """File Cp xuất từ XFLR5: cột x, (Cpi), Cpv. Trả về (x, cp) theo thứ tự trong file."""
    rows, header = [], None
    for line in open(path, errors="ignore"):
        tok = [t for t in re.split(r"[,\s;]+", line.strip()) if t]
        if not tok:
            continue
        try:
            rows.append([float(t) for t in tok])
        except ValueError:
            if "Cp" in line:
                header = tok
    n = max(set(len(r) for r in rows), key=[len(r) for r in rows].count)
    arr = np.array([r for r in rows if len(r) == n])
    if col is None:
        col = header.index("Cpv") if header and "Cpv" in header and len(header) == n else n - 1
    return arr[:, 0], arr[:, col]


def xflr5_contour_cp(x, cp, code, alpha):
    """Tách 2 nhánh tại x nhỏ nhất, nội suy Cp lên contour NACA, trả về (r, cp_nodes)."""
    i0 = int(np.argmin(x))
    first, second = (x[: i0 + 1], cp[: i0 + 1]), (x[i0:], cp[i0:])
    if alpha >= 2 and first[1].mean() > second[1].mean():
        print(f"  [xflr5 alpha={alpha}] nhánh đầu có Cp trung bình cao hơn -> coi là mặt dưới, đã đổi chỗ")
        first, second = second, first
    upper, lower = naca4(code)
    cps = []
    for (bx, bc), geom in ((first, upper), (second, lower)):
        o = np.argsort(bx)
        cps.append(np.interp(geom[:, 0], bx[o], bc[o]))
    return closed_contour(upper, lower), np.r_[cps[0], cps[1][-2:0:-1]]


# OpenFOAM
def _body(path):
    s = open(path, errors="ignore").read()
    return s[s.index("}") + 1:]  # bỏ phần FoamFile


def _patch_range(case, name):
    s = _body(os.path.join(case, "constant", "polyMesh", "boundary"))
    blk = re.search(name + r"\s*\{(.*?)\}", s, re.S).group(1)
    return int(re.search(r"startFace\s+(\d+)", blk).group(1)), int(re.search(r"nFaces\s+(\d+)", blk).group(1))


def _latest_time(case):
    ts = [d for d in os.listdir(case) if re.fullmatch(r"\d+(\.\d+)?", d) and os.path.isfile(os.path.join(case, d, "p"))]
    return max(ts, key=float)


def read_foam_wall_cp(case, patch="walls", umag=26.0):
    """Đọc lưới và p của bước thời gian cuối; trả về contour kín và Cp hằng số trên từng mặt."""
    pm = os.path.join(case, "constant", "polyMesh")
    pts = np.array(re.findall(r"\(\s*(%s)\s+(%s)\s+(%s)\s*\)" % (NUM, NUM, NUM), _body(os.path.join(pm, "points"))), float)
    faces = [[int(v) for v in f.split()] for _, f in re.findall(r"(\d+)\(([\d\s]+)\)", _body(os.path.join(pm, "faces")))]
    owner = np.array(re.findall(r"\b\d+\b", _body(os.path.join(pm, "owner")))[1:], int)
    start, n = _patch_range(case, patch)

    s = open(os.path.join(case, _latest_time(case), "p"), errors="ignore").read()
    pin = np.array(re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\d+\s*\(([^)]*)\)", s, re.S).group(1).split(), float)
    bf = s[s.index("boundaryField"):]
    blk = re.search(r"\b%s\s*\{([^}]*)\}" % patch, bf, re.S).group(1)
    mv = re.search(r"value\s+nonuniform\s+List<scalar>\s*\d+\s*\(([^)]*)\)", blk, re.S)
    p_face = np.array(mv.group(1).split(), float) if mv else pin[owner[start:start + n]]  # zeroGradient: p mặt = p ô kề
    cp_face = p_face / (0.5 * umag**2)  # p động học, p_inf = 0

    zmin = pts[:, 2].min()
    edges = []
    for k in range(n):
        v = [i for i in faces[start + k] if abs(pts[i, 2] - zmin) < 1e-9]
        if len(v) != 2:
            raise RuntimeError("mặt biên dạng không phải tứ giác đùn theo z")
        edges.append((v[0], v[1]))
    adj = {}
    for k, (a, b) in enumerate(edges):
        adj.setdefault(a, []).append(k)
        adj.setdefault(b, []).append(k)
    cur, k, seq, nodes, used = edges[0][0], 0, [], [edges[0][0]], set()
    while True:
        a, b = edges[k]
        nxt = b if a == cur else a
        used.add(k); seq.append(k); nodes.append(nxt); cur = nxt
        nz = [j for j in adj[cur] if j not in used]
        if not nz:
            break
        k = nz[0]
    if len(seq) != n:
        raise RuntimeError("contour biên dạng không khép kín một vòng")
    r = pts[nodes[:-1], :2]
    cp_edges = cp_face[seq]
    centers = 0.5 * (r + np.roll(r, -1, axis=0))
    return r, cp_edges, centers


def alpha_from_name(name):
    m = re.search(r"(?<![A-Za-z])(?:alpha|a)[_=]?(-?\d+(?:\.\d+)?)", os.path.basename(name), re.I)
    return float(m.group(1)) if m else None

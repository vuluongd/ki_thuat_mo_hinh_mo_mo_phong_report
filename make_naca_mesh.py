"""Sinh lưới cho NACA 4 chữ số (chord = 1) bằng Gmsh, xuất naca.msh cho gmshToFoam.

Cách dùng:  python3 make_naca_mesh.py 2412
Patch: inlet (nửa cung trước), outlet (nửa cung sau), walls (biên dạng), frontAndBack (empty).
"""
import sys, math
import gmsh

code = sys.argv[1] if len(sys.argv) > 1 else "2412"
m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
R = 20.0          # bán kính miền (số dây cung)
Y1 = 3e-4         # bề dày lớp đầu tiên (y+ ~ 30 ở Re ~ 2.6e6); đổi nếu tutorial dùng y+ ~ 1
GROWTH = 1.2
NLAYERS = 25
THICK_Z = 0.1     # bề dày theo z  -> AREF = 1 * 0.1
NPTS = 60         # số điểm mỗi mặt


def yt(x):  # bề dày, đuôi đóng
    return 5 * t * (0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)


def yc(x):
    if m == 0:
        return 0.0, 0.0
    if x < p:
        return m / p**2 * (2 * p * x - x * x), 2 * m / p**2 * (p - x)
    return m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x * x), 2 * m / (1 - p) ** 2 * (p - x)


up, lo = [], []
for i in range(NPTS + 1):
    b = math.pi * i / NPTS
    x = 0.5 * (1 - math.cos(b))           # tập trung điểm ở mũi và đuôi
    c, dyc = yc(x)
    th = math.atan(dyc)
    y_t = yt(x)
    up.append((x - y_t * math.sin(th), c + y_t * math.cos(th)))
    lo.append((x + y_t * math.sin(th), c - y_t * math.cos(th)))

gmsh.initialize()
gmsh.model.add("naca")
g = gmsh.model.geo
lc_far = 2.0
# điểm biên dạng: từ đuôi (x=1) ngược lên mũi theo mặt trên, rồi mặt dưới về đuôi
pu = [g.addPoint(x, y, 0, 0.01) for (x, y) in reversed(up)]          # TE -> LE (upper)
pl = [g.addPoint(x, y, 0, 0.01) for (x, y) in lo[1:-1]]              # LE+ -> TE- (lower)
te = pu[0]
le = pu[-1]
pl_full = [le] + pl + [te]
su = g.addSpline(pu)
sl = g.addSpline(pl_full)

# miền ngoài: hai cung tròn quanh (0.5, 0)
cx = 0.5
c0 = g.addPoint(cx, 0, 0, lc_far)
pw = g.addPoint(cx - R, 0, 0, lc_far)    # xa phía trước
pe = g.addPoint(cx + R, 0, 0, lc_far)    # xa phía sau
pn = g.addPoint(cx, R, 0, lc_far)
ps = g.addPoint(cx, -R, 0, lc_far)
a1 = g.addCircleArc(pn, c0, pw)
a2 = g.addCircleArc(pw, c0, ps)
a3 = g.addCircleArc(ps, c0, pe)
a4 = g.addCircleArc(pe, c0, pn)
loop_out = g.addCurveLoop([a1, a2, a3, a4])
loop_af = g.addCurveLoop([su, sl])
surf = g.addPlaneSurface([loop_out, loop_af])
g.synchronize()

# lớp biên
bl = gmsh.model.mesh.field.add("BoundaryLayer")
gmsh.model.mesh.field.setNumbers(bl, "CurvesList", [su, sl])
gmsh.model.mesh.field.setNumber(bl, "Size", Y1)
gmsh.model.mesh.field.setNumber(bl, "Ratio", GROWTH)
gmsh.model.mesh.field.setNumber(bl, "NbLayers", NLAYERS)
gmsh.model.mesh.field.setNumber(bl, "Quads", 1)
gmsh.model.mesh.field.setNumber(bl, "IntersectMetrics", 0)
gmsh.model.mesh.field.setAsBoundaryLayer(bl)

# kích thước ô theo khoảng cách tới biên dạng
d = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(d, "CurvesList", [su, sl])
gmsh.model.mesh.field.setNumber(d, "Sampling", 200)
th_f = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(th_f, "InField", d)
gmsh.model.mesh.field.setNumber(th_f, "SizeMin", 0.01)
gmsh.model.mesh.field.setNumber(th_f, "SizeMax", lc_far)
gmsh.model.mesh.field.setNumber(th_f, "DistMin", 0.2)
gmsh.model.mesh.field.setNumber(th_f, "DistMax", R * 0.8)
gmsh.model.mesh.field.setAsBackgroundMesh(th_f)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

# đùn 1 lớp theo z (lăng trụ / hex), recombine
out = g.extrude([(2, surf)], 0, 0, THICK_Z, numElements=[1], recombine=True)
g.synchronize()
top = out[0][1]
vol = out[1][1]
inlet, outlet, walls = [], [], []
for dim, tag in out[2:]:
    bb = gmsh.model.getBoundingBox(dim, tag)
    # phân loại theo bounding box
    if bb[3] - bb[0] < 3 and bb[4] - bb[1] < 3:
        walls.append(tag)
    elif (bb[0] + bb[3]) / 2 < cx:
        inlet.append(tag)
    else:
        outlet.append(tag)
gmsh.model.addPhysicalGroup(2, inlet, name="inlet")
gmsh.model.addPhysicalGroup(2, outlet, name="outlet")
gmsh.model.addPhysicalGroup(2, walls, name="walls")
gmsh.model.addPhysicalGroup(2, [surf, top], name="frontAndBack")
gmsh.model.addPhysicalGroup(3, [vol], name="fluid")
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.model.mesh.generate(3)
gmsh.write("naca.msh")
types, tags, _ = gmsh.model.mesh.getElements(3)
print("Phan tu 3D:", {int(t): len(g_) for t, g_ in zip(types, tags)})
for nm, tg in [("inlet", inlet), ("outlet", outlet), ("walls", walls)]:
    print(nm, "so mat extrude:", len(tg))
gmsh.finalize()

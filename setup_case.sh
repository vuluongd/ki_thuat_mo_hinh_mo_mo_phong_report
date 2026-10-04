#!/bin/bash
# Tạo case gốc cho NACA bất kỳ từ tutorial airFoil2D, đặt trong ./cases/ cạnh script.
# Dùng:  bash setup_case.sh 2412
DIR=$(cd "$(dirname "$0")" && pwd)
CODE=${1:-2412}
NEW=$DIR/cases/naca${CODE}_base
mkdir -p "$DIR/cases"
rm -rf "$NEW"
cp -r $FOAM_TUTORIALS/incompressibleFluid/airFoil2D "$NEW"
cd "$NEW" || exit 1
rm -rf constant/polyMesh [1-9]* postProcessing log.*

python3 "$DIR/make_naca_mesh.py" "$CODE" || exit 1   # tạo naca.msh
gmshToFoam naca.msh > log.gmshToFoam 2>&1 || { echo "gmshToFoam lỗi, xem $NEW/log.gmshToFoam"; exit 1; }

# walls -> wall, frontAndBack -> empty
python3 - <<'PY'
import re
f = "constant/polyMesh/boundary"
s = open(f).read()
def fix(name, typ):
    global s
    s = re.sub(r"(%s\s*\{[^}]*?type\s+)\w+;" % name, r"\g<1>%s;" % typ, s, count=1)
fix("walls", "wall"); fix("frontAndBack", "empty")
open(f, "w").write(s)
PY

checkMesh > log.checkMesh 2>&1
tail -8 log.checkMesh
echo "Case gốc: $NEW"

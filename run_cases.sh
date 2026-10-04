#!/bin/bash
# Chạy airFoil2D (foamRun) ở nhiều góc tấn.
# Cách dùng:
#   1) cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/airFoil2D ~/airfoil_base
#   2) cd ~/airfoil_base && ./Allrun   (chạy 1 lần để chắc chắn tutorial hoạt động)
#   3) ./Allclean (nếu có) rồi chạy script này:  bash run_cases.sh

DIR=$(cd "$(dirname "$0")" && pwd)
CODE=${CODE:-2412}
BASE=${BASE:-$DIR/cases/naca${CODE}_base}
ANGLES=${ANGLES:-"0 5 10 15"}
UMAG=26.0          # |U| của tutorial (25.75, 3.62) ~ 26 m/s
PATCH=walls        # patch biên dạng (đã kiểm tra từ constant/polyMesh/boundary)
AREF=${AREF:-0.1}  # = chord * bề dày z. Lưới make_naca_mesh.py dùng 0.1; lưới airFoil2D gốc: xem checkMesh

for A in $ANGLES; do
  CASE=$DIR/cases/af_a$A
  rm -rf "$CASE"; cp -r "$BASE" "$CASE"; cd "$CASE" || exit 1
  # xóa kết quả cũ, giữ thư mục 0
  find . -maxdepth 1 -type d -regex './[0-9.]+' ! -name 0 -exec rm -rf {} +
  rm -rf postProcessing log.*

  UX=$(python3 -c "import math;print(round($UMAG*math.cos(math.radians($A)),5))")
  UY=$(python3 -c "import math;print(round($UMAG*math.sin(math.radians($A)),5))")
  LX=$(python3 -c "import math;print(round(-math.sin(math.radians($A)),6))")
  LY=$(python3 -c "import math;print(round(math.cos(math.radians($A)),6))")
  DX=$(python3 -c "import math;print(round(math.cos(math.radians($A)),6))")
  DY=$(python3 -c "import math;print(round(math.sin(math.radians($A)),6))")

  # đổi vector vận tốc dòng tới (thay mọi chỗ có (25.75 3.62 0))
  sed -i "s/(25.75 3.62 0)/($UX $UY 0)/g" 0/U

  cat > system/forceCoeffsDict <<EOF
type            forceCoeffs;
libs            ("libforces.so");
patches         ($PATCH);
rho             rhoInf;
rhoInf          1;
liftDir         ($LX $LY 0);
dragDir         ($DX $DY 0);
CofR            (0.25 0 0);
pitchAxis       (0 0 1);
magUInf         $UMAG;
lRef            1;
Aref            $AREF;
writeControl    timeStep;
writeInterval   1;
EOF

  # cho forceCoeffs chạy cùng solver (ghi vào postProcessing/forceCoeffsDict/0/)
  cat >> system/controlDict <<EOF2

functions
{
    forceCoeffsDict
    {
        #include "forceCoeffsDict"
    }
}
EOF2
  foamRun > log.foamRun 2>&1
  touch case.foam
  echo "Xong alpha=$A  ->  $CASE"
done

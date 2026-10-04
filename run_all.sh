#!/bin/bash
# Chạy một lượt: lưới -> case gốc -> các góc tấn -> so sánh với XFLR5.
# Đặt chung thư mục với: make_naca_mesh.py setup_case.sh run_cases.sh compare.py
# và file polar XFLR5 tên  polar_xflr5.txt  (nếu có).
# Dùng:  bash run_all.sh            (mặc định NACA 2412, góc 0 5 10 15)
#        CODE=4412 ANGLES="0 4 8 12" bash run_all.sh
DIR=$(cd "$(dirname "$0")" && pwd)
export CODE=${CODE:-2412}
cd "$DIR" || exit 1

echo "== 1/4 Tạo lưới và case gốc NACA $CODE"
bash setup_case.sh "$CODE" || exit 1

echo "== 2/4 Chạy thử case gốc"
BASE=$DIR/cases/naca${CODE}_base
( cd "$BASE" && foamRun > log.foamRun 2>&1; tail -4 log.foamRun )
grep -q "converged\|End" "$BASE/log.foamRun" || echo "CẢNH BÁO: case gốc có thể chưa hội tụ, xem $BASE/log.foamRun"

echo "== 3/4 Chạy các góc tấn"
bash run_cases.sh

echo "== 4/4 So sánh với XFLR5"
if [ -f polar_xflr5.txt ]; then
  python3 compare.py --xflr5 polar_xflr5.txt
else
  echo "Chưa có polar_xflr5.txt. Xuất polar từ XFLR5, lưu vào $DIR rồi chạy: python3 compare.py --xflr5 polar_xflr5.txt"
fi

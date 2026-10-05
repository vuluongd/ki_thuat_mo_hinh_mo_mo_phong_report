
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

echo "== 4/4 Tích phân Cp và so sánh với XFLR5"
if ls xflr5_cp/* >/dev/null 2>&1; then
  python3 compare.py --code "$CODE"
else
  echo "Chưa có thư mục xflr5_cp/. Xuất Cp từ XFLR5 (mỗi góc một file: a0.txt, a5.txt, a10.txt, a15.txt), rồi chạy: python3 compare.py --code $CODE"
  python3 compare.py --code "$CODE"
fi

# Khảo sát khí động biên dạng NACA 2412: đối sánh XFLR5 (định lượng) và OpenFOAM (định tính)

Báo cáo giữa kỳ. Thư mục này chứa toàn bộ code và hướng dẫn tái lập kết quả.
Các mục có dấu **[điền]** cần bổ sung số liệu sau khi chạy XFLR5 và chụp hình ParaView.

## 1. Mục tiêu

1. Xác định hệ số khí động (CL, CD, Cm), **tâm áp suất** (Xcp) và **tâm khí động** (Xac) của biên dạng NACA 2412 theo góc tấn α.
2. Đối sánh hai nguồn dữ liệu:
   - **Định lượng:** XFLR5 (XFoil, dòng thế kết hợp lớp biên nhớt).
   - **Định tính:** OpenFOAM (RANS, mô hình rối Spalart-Allmaras), thay cho Ansys Fluent do không có giấy phép.
3. Dùng ảnh trường dòng (contour áp suất, vận tốc, streamline) để nhận xét hiện tượng **tách dòng** ở góc tấn lớn.
4. Kiểm thử code tính Xcp và Xac trên NACA 0012, vì biên dạng này có đáp án biết trước (Xac ≈ 0,25c, Cm quanh c/4 ≈ 0).

Điều kiện khảo sát: chord c = 1 m, vận tốc dòng tới |U| = 26 m/s, độ nhớt động học ν = 1×10⁻⁵ m²/s (giá trị của tutorial, cần đối chiếu `constant/transportProperties`), suy ra **Re ≈ 2,6×10⁶**. Góc tấn: 0°, 5°, 10°, 15°.

## 2. Cấu trúc thư mục

| File | Vai trò |
|---|---|
| `run_all.sh` | Chạy một lượt: lưới, case gốc, các góc tấn, so sánh |
| `setup_case.sh` | Tạo lưới cho NACA bất kỳ và dựng case gốc từ tutorial `airFoil2D` |
| `make_naca_mesh.py` | Sinh biên dạng NACA 4 chữ số và lưới bằng Gmsh, xuất `naca.msh` |
| `run_cases.sh` | Với mỗi α: đổi hướng dòng tới, khai báo `forceCoeffs`, chạy solver |
| `compare.py` | Đọc kết quả OpenFOAM và polar XFLR5, tính Xcp, Xac, vẽ `so_sanh.png` |
| `polar_xflr5.txt` | Polar xuất từ XFLR5 (tự tạo, xem mục 5) |
| `cases/` | Sinh ra khi chạy: case gốc `naca2412_base` và `af_a0`, `af_a5`, `af_a10`, `af_a15` |

## 3. Phương pháp và giải thích code

### 3.1 Hình học biên dạng và lưới (`make_naca_mesh.py`)

- **Hình học:** mã NACA 4 chữ số `MPTT` (2412) cho độ cong lớn nhất m = 2%, vị trí cong lớn nhất p = 40% dây cung, độ dày t = 12%. Đường camber và phân bố độ dày tính theo công thức chuẩn NACA, mỗi điểm được dịch vuông góc với đường camber. Hệ số cuối của phân bố độ dày dùng −0,1036 để **đuôi đóng kín**. Điểm phân bố theo cosine để dày hơn ở mũi và đuôi (60 điểm mỗi mặt).
- **Miền tính:** hình tròn bán kính 20c tâm (0,5; 0), chia thành 4 cung. Hai cung phía trước là `inlet`, hai cung phía sau là `outlet`.
- **Lớp biên:** 25 lớp ô tứ giác, ô đầu tiên dày 3×10⁻⁴ m (y⁺ cỡ 30 ở Re ≈ 2,6×10⁶, phù hợp wall function), tỉ lệ tăng 1,2.
- **Kích thước ô:** nhỏ (0,01) gần biên dạng, tăng dần đến 2,0 ở biên xa.
- **Đùn 3D:** lưới 2D được kéo dày 0,1 theo z đúng một lớp ô (OpenFOAM luôn giải 3D). Hai mặt trước và sau được đặt là `empty` để bài toán thành 2D. Vì vậy diện tích tham chiếu **Aref = c × 0,1 = 0,1**.
- **Gán tên patch:** `inlet`, `outlet`, `walls` (biên dạng), `frontAndBack`. Tên khớp với các file điều kiện biên của tutorial nên không phải sửa `0/`.
- **Kết quả:** khoảng 25.000 ô (hex ở vùng lớp biên, lăng trụ ở vùng ngoài), xuất `naca.msh` định dạng 2.2 để `gmshToFoam` đọc.

### 3.2 Thiết lập OpenFOAM (`setup_case.sh`, `run_cases.sh`)

- **Phiên bản:** OpenFOAM Foundation, nhánh `dev`. Solver là module `incompressibleFluid`, chạy bằng lệnh `foamRun` (thuật toán SIMPLE, dừng khi đạt tiêu chí hội tụ hoặc sau 500 vòng).
- **Case gốc:** sao chép tutorial `incompressibleFluid/airFoil2D` (điều kiện biên, schemes, mô hình rối Spalart-Allmaras với biến `nuTilda`), thay lưới bằng lưới NACA 2412.
- **Điều kiện biên:** `inlet` và `outlet` kiểu `freestreamVelocity`/`freestream`, thành biên dạng `noSlip` với wall function cho độ nhớt rối.
- **Đổi góc tấn:** không xoay lưới mà xoay vector vận tốc dòng tới trong `0/U`: Ux = |U|cos α, Uy = |U|sin α.
- **Hệ số lực:** function `forceCoeffs` (thư viện `libforces.so`) chạy cùng solver trên patch `walls`:
  - `liftDir` = (−sin α, cos α, 0), `dragDir` = (cos α, sin α, 0) bám theo hướng dòng tới,
  - moment lấy quanh `CofR` = (0,25; 0; 0), tức **điểm 1/4 dây cung**, trục `pitchAxis` = (0, 0, 1),
  - `lRef` = 1, `Aref` = 0,1, `rhoInf` = 1 (áp suất trong solver là áp suất động học nên hệ số vẫn đúng).
- **Đầu ra:** `cases/af_aX/postProcessing/forceCoeffsDict/0/forceCoeffs.dat`, các cột `Time Cm Cd Cl Cl(f) Cl(r)`; script lấy dòng cuối.

### 3.3 Dữ liệu định lượng XFLR5

Biên dạng NACA 2412 tạo trong *Direct Foil Design*; phân tích trong *XFoil Direct Analysis* với polar Type 1, Re = 2,6×10⁶, Mach = 0, Ncrit = 9, quét α từ 0° đến 16° bước 1°. Xuất bằng *File > Export > Polars* thành `polar_xflr5.txt` (cột `alpha CL CD CDp Cm`, Cm quanh c/4). XFoil không mô phỏng được vùng tách dòng lớn, nên ở α cao dữ liệu chỉ có giá trị tham khảo.

### 3.4 Xử lý số liệu (`compare.py`)

- **Tâm áp suất:** Xcp/c = 0,25 − Cm/CL (Cm quanh c/4, gần đúng CN ≈ CL ở góc nhỏ). Bỏ các điểm |CL| < 0,05 vì phép chia mất ý nghĩa.
- **Tâm khí động:** Xac/c = 0,25 − dCm/dCL, với độ dốc dCm/dCL là hệ số góc của đường hồi quy tuyến tính Cm theo CL trên các góc ≤ 10° (vùng tuyến tính).
- **Quy ước dấu:** OpenFOAM trả Cm dương cho biên dạng camber (nose-up), ngược quy ước thường dùng (nose-down âm, như XFLR5). Code **đổi dấu Cm của OpenFOAM**; có thể tắt bằng `--keep-cm-sign`.
- **Đầu ra:** `so_sanh.png` gồm 4 đồ thị CL–α, Cm–α, Xcp–α, CL–CD; và hai giá trị Xac in ra màn hình.

## 4. Cách chạy

```bash
pip install gmsh
bash run_all.sh                              # mặc định NACA 2412, góc 0 5 10 15
CODE=4412 ANGLES="0 4 8 12" bash run_all.sh  # đổi biên dạng hoặc góc
python3 compare.py --xflr5 polar_xflr5.txt   # nếu polar XFLR5 có sau
```

Hình trường dòng: `paraview cases/af_a10/case.foam`, chọn thời điểm cuối, tô màu theo `p` và `U`, dùng *Stream Tracer* hoặc tô màu theo `U_X` (vùng âm ở mặt hút là dòng hồi lưu).

## 5. Nguồn dữ liệu

| Dữ liệu | Nguồn |
|---|---|
| Hình học biên dạng | Công thức NACA 4 chữ số, tính trong `make_naca_mesh.py` |
| Hệ số CL, CD, Cm (định tính, OpenFOAM) | Chạy `run_cases.sh`, đọc `forceCoeffs.dat` |
| Hệ số CL, CD, Cm (định lượng, XFLR5) | Người làm tự chạy XFLR5, file `polar_xflr5.txt` |
| Trường áp suất, vận tốc, streamline | ParaView mở `cases/af_aX/case.foam` |
| Xcp, Xac | Tính từ hai nguồn trên bằng `compare.py` |
| Template case OpenFOAM | Tutorial `airFoil2D` đi kèm OpenFOAM |

## 6. Kết quả

### 6.1 OpenFOAM (đã chạy, residual cuối cỡ 10⁻⁶ đến 10⁻⁵, 500 vòng, chưa đạt tiêu chí hội tụ)

| α (°) | CL | CD | Cm (c/4, sau khi đổi dấu) | Xcp/c |
|---|---|---|---|---|
| 0 | 0,2016 | 0,0238 | −0,0491 | 0,49 |
| 5 | 0,6724 | 0,0365 | −0,0464 | 0,32 |
| 10 | 1,0825 | 0,0699 | −0,0410 | 0,29 |
| 15 | 1,3648 | 0,1239 | −0,0388 | 0,28 |

Tâm khí động từ độ dốc Cm–CL trên 0° đến 10°: Xac/c ≈ 0,24.

### 6.2 XFLR5 **[điền]**

| α (°) | CL | CD | Cm (c/4) | Xcp/c |
|---|---|---|---|---|
| 0 | [điền] | [điền] | [điền] | [điền] |
| 5 | [điền] | [điền] | [điền] | [điền] |
| 10 | [điền] | [điền] | [điền] | [điền] |
| 15 | [điền] | [điền] | [điền] | [điền] |

Xac/c (XFLR5): **[điền]** (in ra từ `compare.py`).

### 6.3 Bảng đối sánh định lượng và định tính

| Nội dung | Định lượng (XFLR5) | Định tính (OpenFOAM) |
|---|---|---|
| Phân bố áp suất | Đồ thị Cp theo x/c **[chèn hình]** | Contour áp suất tại 0°, 5°, 10°, 15° **[chèn hình]** |
| Tâm áp suất Xcp | **[điền]** | 0,49 → 0,32 → 0,29 → 0,28 |
| Tâm khí động Xac | **[điền]** | ≈ 0,24c |
| Tách dòng | Chỉ gián tiếp: CL–α cong xuống, Cp mặt hút phẳng | Contour U, streamline, vùng U_X < 0 ở mặt hút **[chèn hình, nhất là 10° và 15°]** |
| Nhận xét | **[điền]** | **[điền]** |

## 7. Nhận xét và kết luận (dự kiến, kiểm lại theo số liệu cuối cùng)

1. **Tâm khí động:** cả hai nguồn cho Xac gần 0,25c như lý thuyết biên dạng mỏng, cho thấy quy trình tính Xac đúng. Phép kiểm trên NACA 0012 nên cho Xac ≈ 0,25c và Cm ≈ 0.
2. **Tâm áp suất:** Xcp dịch dần về phía mũi khi α tăng (từ khoảng 0,49c về khoảng 0,28c trong OpenFOAM), đúng quy luật của biên dạng có camber, và tiến dần đến Xac khi CL lớn.
3. **Đường CL–α:** OpenFOAM cho độ dốc khoảng 0,09 mỗi độ và CL thấp hơn lý thuyết cỡ 10 đến 20%, đường cong chậm lại ở 10° đến 15°, cho thấy tiệm cận trạng thái CLmax và tách dòng bắt đầu ở mặt hút. XFLR5 cho độ dốc gần 0,1 mỗi độ ở vùng tuyến tính **[kiểm lại bằng số liệu XFLR5]**.
4. **Tách dòng:** XFLR5 không mô tả được vùng hồi lưu; chỉ OpenFOAM cho thấy trực tiếp dòng tách ở mặt hút khi α lớn **[xác nhận bằng ảnh streamline]**. Hai nguồn lệch nhau nhiều nhất ở 10° đến 15°.
5. **Kết luận chung:** hai phương pháp bổ sung cho nhau: XFLR5 nhanh và cho số liệu định lượng tin cậy ở góc nhỏ, OpenFOAM cho bức tranh trường dòng và dấu hiệu tách dòng ở góc lớn. Với lưới hiện tại OpenFOAM chỉ nên dùng để nhận xét xu hướng và hiện tượng, không dùng làm số liệu định lượng tuyệt đối.

## 8. Hạn chế và hướng phát triển

- **CD của OpenFOAM cao bất thường** (≈ 0,024 ở 0°, trong khi giá trị dòng rối hoàn toàn cỡ 0,008 đến 0,01). Nguyên nhân có thể: lưới và y⁺ chưa tối ưu, đuôi đóng kín, chưa hội tụ hẳn. Cần kiểm tra y⁺, tăng số vòng lặp hoặc làm mịn lưới.
- Mô hình rối Spalart-Allmaras không dự đoán chính xác CLmax và vùng tách dòng; có thể thử k-ω SST.
- XFLR5 giả thiết dòng gần dính, không áp dụng tin cậy sau góc tách.
- Chưa có so sánh phân bố Cp chồng hai nguồn; cần xuất Cp từ XFLR5 và từ OpenFOAM (ParaView, *Plot On Intersection Curves*).
- Hướng làm tiếp cuối kỳ: kiểm định độc lập lưới, thêm góc tấn gần vùng tách (12° đến 18°), so sánh với số liệu thực nghiệm công bố cho NACA 2412.

## 9. Tài liệu tham khảo

1. I. H. Abbott, A. E. von Doenhoff, *Theory of Wing Sections*, Dover, 1959 (công thức NACA 4 chữ số).
2. M. Drela, "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils", 1989; tài liệu XFLR5.
3. OpenFOAM Foundation, *OpenFOAM User Guide* và tutorial `incompressibleFluid/airFoil2D`, openfoam.org.
4. C. Geuzaine, J.-F. Remacle, "Gmsh: a three-dimensional finite element mesh generator", 2009; tài liệu tại gmsh.info.
5. P. R. Spalart, S. R. Allmaras, "A one-equation turbulence model for aerodynamic flows", 1992.

*Ghi chú:* các script được xây dựng với sự hỗ trợ của công cụ AI (Claude); người thực hiện đã chạy, kiểm tra và chịu trách nhiệm về kết quả. Khai báo mục này nếu giảng viên yêu cầu.

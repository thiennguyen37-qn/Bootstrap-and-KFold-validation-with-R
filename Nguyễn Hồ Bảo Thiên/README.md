# Hướng Dẫn Chạy Mã Nguồn Tiểu Luận

Thư mục này chứa các tệp cần thiết để chạy lại mã nguồn phân tích trong R và biên dịch lại báo cáo khi cần.

## 1. Danh Sách Tệp

- `tieu_luan.pdf`: báo cáo hoàn chỉnh.
- `tieu_luan_code.R`: mã nguồn R dùng để chạy lại các kết quả trong báo cáo.
- `tieu_luan.Rnw`: file nguồn dùng để biên dịch lại báo cáo PDF.
- `Indian Liver Patient Dataset (ILPD).csv`: dữ liệu đầu vào.

## 2. Yêu Cầu Trước Khi Chạy

Máy cần cài đặt:

- R
- RStudio hoặc VS Code

Khi chạy mã nguồn, cần giữ `tieu_luan_code.R` và `Indian Liver Patient Dataset (ILPD).csv` trong cùng một thư mục. Không đổi tên file dữ liệu nếu chưa chỉnh lại đường dẫn trong code.

## 3. Chạy Bằng RStudio

1. Mở RStudio.
2. Mở file `tieu_luan_code.R`.
3. Đặt working directory về thư mục đang chứa file:

   Chọn `Session` -> `Set Working Directory` -> `To Source File Location`.

4. Chạy toàn bộ file bằng một trong các cách sau:

- Nhấn `Ctrl + Shift + Enter`.
- Chọn `Code` -> `Run Region` -> `Run All`.
- Hoặc chạy lệnh sau trong Console:

```r
source("tieu_luan_code.R")
```

Kết quả dạng chữ sẽ được in trong Console. Các biểu đồ sẽ xuất hiện trong tab `Plots`; dùng các nút mũi tên trong tab này để xem lần lượt những biểu đồ đã tạo. Khi chạy theo cách này, chương trình không tự lưu biểu đồ thành file.

## 4. Chạy Bằng VS Code

1. Mở VS Code.
2. Chọn `File` -> `Open Folder...`.
3. Chọn thư mục đang chứa `tieu_luan_code.R` và file dữ liệu CSV.
4. Mở Terminal trong VS Code bằng `Terminal` -> `New Terminal`.
5. Chạy lệnh:

```text
Rscript tieu_luan_code.R
```

Kết quả dạng chữ sẽ được in trực tiếp trong Terminal. Sau khi chương trình chạy xong, toàn bộ biểu đồ sẽ được lưu vào file `Rplots.pdf` trong chính thư mục đang mở. Đây là một file PDF nhiều trang, mỗi trang chứa một biểu đồ. Nếu chạy lại chương trình, file `Rplots.pdf` cũ sẽ bị ghi đè.

Nếu Terminal không nhận lệnh `Rscript`, cần kiểm tra lại việc cài đặt R hoặc thêm R vào biến môi trường `PATH`.

## 5. Các Biểu Đồ Được Tạo

Khi chạy `tieu_luan_code.R`, chương trình tạo các biểu đồ phục vụ cho phần phân tích, gồm phân bố biến mục tiêu, ma trận tương quan, boxplot trước và sau biến đổi log, phân phối bootstrap, forest plot, đường đánh giá ngưỡng phân loại và ma trận nhầm lẫn.

Vị trí nhận biểu đồ phụ thuộc vào cách chạy:

- Chạy bằng RStudio: biểu đồ xuất hiện trong tab `Plots`.
- Chạy lệnh `Rscript tieu_luan_code.R` trong Terminal của VS Code hoặc RStudio: biểu đồ được gom vào file `Rplots.pdf`.
- Biên dịch `tieu_luan.Rnw` bằng `knitr`: các biểu đồ được lưu riêng trong thư mục `figure` với tên bắt đầu bằng `tl-`, sau đó được chèn vào báo cáo.

## 6. Biên Dịch Lại Báo Cáo PDF

Thông thường không cần biên dịch lại PDF vì file `tieu_luan.pdf` đã có sẵn. Chỉ thực hiện bước này khi muốn tạo lại báo cáo từ file nguồn `tieu_luan.Rnw`.

Cần cài thêm:

- Gói R `knitr`
- Một bản phân phối LaTeX có `xelatex`

Trong R hoặc RStudio, chạy:

```r
knitr::knit("tieu_luan.Rnw")
```

Lệnh trên tạo file `tieu_luan.tex` và các file biểu đồ trong thư mục `figure`. Sau đó biên dịch `tieu_luan.tex` bằng `xelatex` để tạo lại báo cáo PDF.

## 7. Ghi Chu Ve Pipeline Model

Phan mo hinh trong `tieu_luan_code.R` goi `python_reference_pipeline.py` de dung cung pipeline voi notebook Python `03_tuning_and_bootstrap_ci.ipynb`: sklearn `train_test_split`, imblearn `SMOTE`, `GridSearchCV`, `LogisticRegression` va nguong phan loai chon tu out-of-fold tren tap train.

Vi vay khi chay R can co Python va cac goi `pandas`, `numpy`, `scikit-learn`, `imbalanced-learn`. Cach cai nhanh:

```text
pip install pandas numpy scikit-learn imbalanced-learn
```

Cach lam nay giup ket qua R va Python nam tren cung he quy chieu, khong con dung pipeline Logistic/SMOTE tu viet trong R nhu truoc.

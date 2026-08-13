# Hướng Dẫn Chạy Lại Báo Cáo

Thư mục này chứa bản nộp cuối cùng của tiểu luận. Toàn bộ pipeline mô hình trong
`tieu_luan.Rnw` và `tieu_luan_code.R` được thực hiện hoàn toàn bằng R, không
gọi script tham chiếu bên ngoài.

## Danh Sách Tệp

- `tieu_luan.pdf`: báo cáo hoàn chỉnh sau khi biên dịch.
- `tieu_luan.Rnw`: nguồn chính để tạo lại báo cáo PDF.
- `tieu_luan_code.R`: mã R để chạy lại các bước phân tích và mô hình.
- `Indian Liver Patient Dataset (ILPD).csv`: dữ liệu đầu vào.
- `logo_qnu.png`: logo dùng trên trang bìa.

## Yêu Cầu

Máy cần cài:

- R
- RStudio hoặc VS Code
- Một bản LaTeX có `xelatex`, ví dụ MiKTeX hoặc TeX Live
- Gói R `knitr` và `jsonlite`

Cài gói R nếu máy chưa có:

```r
install.packages(c("knitr", "jsonlite"))
```

## Cách 1: Chạy Bằng VS Code Hoặc Terminal

Mở terminal tại thư mục `FINAL REPORT`.

Chạy lại toàn bộ mã phân tích:

```text
Rscript tieu_luan_code.R
```

Tạo file LaTeX từ nguồn `.Rnw`:

```text
Rscript -e "knitr::knit('tieu_luan.Rnw')"
```

Biên dịch PDF bằng `xelatex` ít nhất hai lần:

```text
xelatex tieu_luan.tex
xelatex tieu_luan.tex
```

Sau bước này, file `tieu_luan.pdf` sẽ được tạo lại trong thư mục `FINAL REPORT`.

## Cách 2: Chạy Bằng RStudio

Mở RStudio, sau đó mở thư mục `FINAL REPORT` hoặc đặt working directory về thư
mục này.

Trong Console của RStudio, chạy:

```r
source("tieu_luan_code.R")
```

Tiếp theo, tạo file LaTeX từ nguồn `.Rnw`:

```r
knitr::knit("tieu_luan.Rnw")
```

Sau đó mở terminal hoặc tab Terminal trong RStudio và chạy:

```text
xelatex tieu_luan.tex
xelatex tieu_luan.tex
```

Nếu RStudio có nút Knit/Compile PDF và đã nhận đúng engine LaTeX, cũng có thể
dùng nút đó để biên dịch file `tieu_luan.Rnw`.

## File Có Thể Phát Sinh Khi Build

Khi build lại báo cáo, có thể phát sinh thêm các file và thư mục phụ như:

- `figure/`: thư mục chứa hình do `knitr` xuất ra.
- `tieu_luan.tex`: file LaTeX được tạo từ `tieu_luan.Rnw`.
- `tieu_luan.aux`, `tieu_luan.log`, `tieu_luan.out`, `tieu_luan.toc`: file phụ
  của LaTeX.
- `Rplots.pdf`: file plot phụ nếu R sinh hình ngoài thiết bị của `knitr`.


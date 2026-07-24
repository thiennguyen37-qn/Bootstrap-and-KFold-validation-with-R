# ILPD (Indian Liver Patient Dataset)

## 1. Giới thiệu chung

Bộ dữ liệu này được thu thập từ khu vực Đông Bắc bang Andhra Pradesh, Ấn Độ, gồm hồ sơ của các bệnh nhân được xét nghiệm để xác định có mắc xơ gan hay không. Mục tiêu của bộ dữ liệu là phục vụ bài toán **phân loại nhị phân**: dự đoán một bệnh nhân có mắc bệnh xơ gan hay không, dựa trên các chỉ số sinh hóa liên quan đến chức năng gan (như bilirubin, albumin, các enzyme chuyển hóa).

Bộ dữ liệu được công bố trên UCI Machine Learning Repository, đóng góp bởi Bendi Ramana và N. Venkateswarlu, ngày 20/05/2012.

## 2. Thông tin tổng quan

- **Lĩnh vực:** Y tế / Sức khỏe (Health and Medicine)
- **Số quan sát:** 583
- **Số biến:** 10 đặc trưng + 1 biến mục tiêu
- **Kiểu dữ liệu biến:** Số nguyên và số thực
- **Giấy phép:** Creative Commons Attribution 4.0 International (CC BY 4.0)
- **DOI:** 10.24432/C5D02C

## 3. Bảng mô tả các biến

| Tên biến | Mô tả |
|---|---|
| Age | Tuổi của bệnh nhân. Những bệnh nhân trên 89 tuổi được ghi nhận là 90. |
| Gender | Giới tính của bệnh nhân (Male/Female). |
| TB | Nồng độ bilirubin toàn phần trong máu. |
| DB | Nồng độ bilirubin trực tiếp (bilirubin liên hợp). |
| Alkphos | Nồng độ enzyme phosphatase kiềm (ALP). |
| Sgpt | Nồng độ enzyme Alanine Aminotransferase (ALT). |
| Sgot | Nồng độ enzyme Aspartate Aminotransferase (AST). |
| TP | Tổng lượng protein trong huyết thanh. |
| ALB | Nồng độ albumin trong huyết thanh. |
| A/G Ratio | Tỷ lệ giữa albumin và globulin trong máu. |
| Selector (target) | Nhãn phân loại: có mắc bệnh gan hay không |

## 4. Cách truy cập dữ liệu

Tải trực tiếp file CSV (khoảng 23 KB) từ trang dữ liệu gốc trên UCI Machine Learning Repository:

[https://archive.ics.uci.edu/dataset/225/ilpd+indian+liver+patient+dataset](https://archive.ics.uci.edu/dataset/225/ilpd+indian+liver+patient+dataset)

Sau khi tải về, đọc file CSV trực tiếp trong R bằng:

```r
df <- read.csv("Indian Liver Patient Dataset (ILPD).csv", header = FALSE)
colnames(df) <- c("Age", "Gender", "TB", "DB", "Alkphos", "Sgpt", "Sgot", "TP", "ALB", "A/G_Ratio", "Selector")
```

## 5. Trích dẫn nguồn

> Ramana, B. & Venkateswarlu, N. (2022). ILPD (Indian Liver Patient Dataset) [Dataset]. UCI Machine Learning Repository. https://doi.org/10.24432/C5D02C



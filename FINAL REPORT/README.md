# Huong Dan Chay Lai Bao Cao

Thu muc nay chua ban nop cuoi cung cua tieu luan. Toan bo pipeline mo hinh trong
`tieu_luan.Rnw` va `tieu_luan_code.R` duoc thuc hien hoan toan bang R, khong
goi script tham chieu ben ngoai.

## Danh Sach Tep

- `tieu_luan.pdf`: bao cao hoan chinh sau khi bien dich.
- `tieu_luan.Rnw`: nguon chinh de tao lai bao cao PDF.
- `tieu_luan_code.R`: ma R de chay lai cac buoc phan tich va mo hinh.
- `Indian Liver Patient Dataset (ILPD).csv`: du lieu dau vao.
- `logo_qnu.png`: logo dung tren trang bia.

## Yeu Cau

May can cai:

- R
- RStudio hoac VS Code
- Mot ban phan phoi LaTeX co `xelatex`
- Goi R `knitr` va `jsonlite`

Cai goi R neu may chua co:

```r
install.packages(c("knitr", "jsonlite"))
```

## Chay Lai Ma Phan Tich

Mo terminal tai thu muc nay va chay:

```text
Rscript tieu_luan_code.R
```

Lenh nay chay lai EDA, lua chon dac trung, chia train-test, SMOTE, logistic
regression co dieu chuan, grid search, chon nguong, danh gia test va bootstrap.

## Tao Lai PDF

Trong R hoac RStudio, chay:

```r
knitr::knit("tieu_luan.Rnw")
```

Sau do bien dich file `tieu_luan.tex` bang `xelatex` it nhat hai lan:

```text
xelatex tieu_luan.tex
xelatex tieu_luan.tex
```

Neu dung VS Code hoac terminal, co the chay cac lenh tren ngay trong thu muc nay.

## Ghi Chu Ve Pipeline R

Pipeline R trong ban nay tu cai dat cac buoc sau:

- Dien khuyet bang median tinh tren tap train cua tung fold.
- Chuan hoa theo cong thuc kieu `StandardScaler`, dung do lech chuan tong the
  `sqrt(mean((x - mean)^2))`.
- SMOTE chi thuc hien tren tap train cua tung fold, theo cach lay mau cap
  lang gieng va noi suy giong `imblearn` nhat co the trong R.
- Logistic regression toi uu truc tiep ham mat mat logistic co dieu chuan theo
  he so `C`.
- Grid search gom `C = {0.1, 1, 10}`, dang phat `L1/L2`, va
  `SMOTE k = {3, 5, 7}`.

Khong can them bat ky script ngoai R nao de tao lai ket qua trong thu muc
FINAL REPORT.

# =============================================================================
# PHÂN TÍCH BỘ DỮ LIỆU INDIAN LIVER PATIENT (ILPD)
# Bootstrap và K-Fold Cross-Validation
#
# Bản R gộp từ ba notebook Python:
#   01_eda.ipynb                    -> PHẦN 1, 2
#   02_bootstrap_kfold_validation   -> PHẦN 4, 5
#   03_tuning_and_bootstrap_ci      -> PHẦN 6, 7, 8
#
# CHỈ DÙNG BASE R + stats. Không cần caret / dplyr / glmnet / ggplot2 -
# các gói biên dịch đang bị Smart App Control chặn trên máy này. Việc tự cài đặt
# SMOTE, VIF, K-Fold, Bootstrap và hồi quy ridge cũng làm rõ cơ chế bên trong
# hơn là gọi một hàm có sẵn.
# =============================================================================

set.seed(42)

# =============================================================================
# PHẦN 1. NẠP DỮ LIỆU VÀ TIỀN XỬ LÝ
# =============================================================================

col_names <- c(
  "Age", "Gender", "Total_Bilirubin", "Direct_Bilirubin",
  "Alkaline_Phosphotase", "Alamine_Aminotransferase", "Aspartate_Aminotransferase",
  "Total_Protiens", "Albumin", "Albumin_and_Globulin_Ratio", "target"
)

df <- read.csv("Indian Liver Patient Dataset (ILPD).csv",
               header = FALSE, col.names = col_names, stringsAsFactors = FALSE)

# Nhãn gốc: 1 = có bệnh, 2 = không bệnh. Đổi về 1 / 0 cho tiện tính toán.
df$target <- ifelse(df$target == 2, 0, 1)

cat("Kích thước dữ liệu:", nrow(df), "dòng x", ncol(df), "cột\n")
cat("Phân bố target:", sum(df$target == 1), "có bệnh /",
    sum(df$target == 0), "không bệnh\n\n")

str(df)
summary(df)

# Số ô khuyết theo từng cột
cat("\nSố giá trị khuyết theo cột:\n")
print(colSums(is.na(df)))


# =============================================================================
# PHẦN 2. PHÂN TÍCH KHÁM PHÁ (EDA)
# =============================================================================

# ---- 2.1 Phân bố target -----------------------------------------------------
counts <- table(df$target)
pct <- round(100 * counts / sum(counts), 1)

bp <- barplot(counts,
              names.arg = c("Không bệnh (0)", "Có bệnh (1)"),
              col = c("#66c2a5", "#fc8d62"),
              main = "Phân bố target", ylab = "Số lượng",
              ylim = c(0, max(counts) * 1.15))
text(bp, counts, labels = paste0(counts, "\n(", pct, "%)"), pos = 3)


# ---- 2.2 Quan hệ Gender và target: kiểm định Chi-square ---------------------
# H0: tỉ lệ mắc bệnh không khác nhau giữa hai giới tính
contingency <- table(df$Gender, df$target)
cat("\n--- Bảng tần số Gender x target ---\n")
print(contingency)

chi_res <- chisq.test(contingency)
cat("\nChi-square =", round(chi_res$statistic, 4),
    "| p-value =", round(chi_res$p.value, 4), "\n")
cat("Bảng tần số kỳ vọng (cần >= 5 để kiểm định hợp lệ):\n")
print(round(chi_res$expected, 2))

if (chi_res$p.value < 0.05) {
  cat("=> Bác bỏ H0: có khác biệt về tỉ lệ mắc bệnh giữa hai giới.\n")
} else {
  cat("=> Không đủ bằng chứng bác bỏ H0: chưa thấy khác biệt giữa hai giới.\n")
}


# ---- 2.3 Ma trận tương quan -------------------------------------------------
numeric_cols <- c(
  "Age", "Total_Bilirubin", "Direct_Bilirubin",
  "Alkaline_Phosphotase", "Alamine_Aminotransferase", "Aspartate_Aminotransferase",
  "Total_Protiens", "Albumin", "Albumin_and_Globulin_Ratio"
)

corr_mat <- cor(df[, c(numeric_cols, "target")], use = "complete.obs")
cat("\n--- Ma trận tương quan ---\n")
print(round(corr_mat, 2))

# Heatmap bằng base R: image() vẽ ma trận, trục y phải đảo cho đúng chiều đọc
par(mar = c(9, 9, 3, 2))
image(seq_len(ncol(corr_mat)), seq_len(nrow(corr_mat)),
      t(corr_mat[nrow(corr_mat):1, ]),
      col = colorRampPalette(c("#3288bd", "white", "#d53e4f"))(64),
      zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "",
      main = "Ma trận tương quan")
axis(1, seq_len(ncol(corr_mat)), colnames(corr_mat), las = 2, cex.axis = 0.7)
axis(2, seq_len(nrow(corr_mat)), rev(rownames(corr_mat)), las = 2, cex.axis = 0.7)
par(mar = c(5, 4, 4, 2) + 0.1)


# ---- 2.4 VIF: kiểm tra đa cộng tuyến ----------------------------------------
# VIF của biến j = 1 / (1 - R^2_j), với R^2_j lấy từ hồi quy biến j lên TẤT CẢ
# các biến còn lại. Tự tính bằng lm() nên không cần gói car.
compute_vif <- function(data) {
  vars <- colnames(data)
  sapply(vars, function(v) {
    others <- setdiff(vars, v)
    fml <- as.formula(paste(v, "~", paste(others, collapse = " + ")))
    r2 <- summary(lm(fml, data = data))$r.squared
    1 / (1 - r2)
  })
}

# Chuẩn bị ma trận thiết kế: Gender -> nhị phân, điền khuyết bằng median
X_vif <- df[, c(numeric_cols, "Gender")]
X_vif$Gender <- ifelse(X_vif$Gender == "Male", 1, 0)
X_vif$Albumin_and_Globulin_Ratio[is.na(X_vif$Albumin_and_Globulin_Ratio)] <-
  median(X_vif$Albumin_and_Globulin_Ratio, na.rm = TRUE)

vif_raw <- compute_vif(X_vif)
cat("\n--- VIF trên thang đo gốc ---\n")
print(round(sort(vif_raw, decreasing = TRUE), 3))
cat("Ngưỡng cảnh báo thường dùng: VIF > 5 đáng lo, VIF > 10 nghiêm trọng.\n")


# ---- 2.5 Boxplot 9 biến số theo target --------------------------------------
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
for (col in numeric_cols) {
  boxplot(df[[col]] ~ df$target, col = c("#66c2a5", "#fc8d62"),
          main = col, xlab = "target", ylab = "")
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)


# ---- 2.6 log1p cho các biến lệch phải nặng ----------------------------------
# Outlier ở đây là tín hiệu lâm sàng thật (AST ~ 5000 = suy gan cấp) và tập trung
# ở nhóm bệnh, nên winsorization sẽ kẹp mất khả năng phân biệt. log1p đơn điệu
# tăng nên BẢO TOÀN THỨ HẠNG - p-value Mann-Whitney không đổi - mà vẫn nén đuôi.
skewed_cols <- c("Total_Bilirubin", "Direct_Bilirubin", "Alkaline_Phosphotase",
                 "Alamine_Aminotransferase", "Aspartate_Aminotransferase")

df_raw <- df                                  # giữ bản gốc để tra thang đo lâm sàng
for (col in skewed_cols) df[[col]] <- log1p(df[[col]])

# Hệ số bất đối xứng, tự tính vì base R không có hàm skewness
skewness <- function(x) {
  x <- x[!is.na(x)]
  mean((x - mean(x))^3) / (mean((x - mean(x))^2))^1.5
}
cat("\n--- Skewness trước / sau log1p ---\n")
print(round(data.frame(
  goc = sapply(df_raw[skewed_cols], skewness),
  sau_log1p = sapply(df[skewed_cols], skewness)
), 2))

# Boxplot lại sau khi biến đổi để so sánh trực quan
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
for (col in numeric_cols) {
  boxplot(df[[col]] ~ df$target, col = c("#66c2a5", "#fc8d62"),
          main = if (col %in% skewed_cols) paste0("log1p(", col, ")") else col,
          xlab = "target", ylab = "")
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# VIF tính lại trên thang log: đây mới là ma trận thiết kế mà mô hình thực sự dùng.
# Pearson/VIF KHÔNG bất biến với log (khác Mann-Whitney) nên con số sẽ lệch.
X_vif_log <- X_vif
for (col in skewed_cols) X_vif_log[[col]] <- log1p(df_raw[[col]])
vif_log <- compute_vif(X_vif_log)

cat("\n--- VIF: thang gốc vs thang log ---\n")
print(round(data.frame(
  VIF_goc = vif_raw, VIF_log = vif_log[names(vif_raw)],
  chenh_lech = vif_log[names(vif_raw)] - vif_raw
)[order(-vif_log[names(vif_raw)]), ], 2))


# ---- 2.7 Thống kê theo nhóm (dùng df_raw để giữ đơn vị lâm sàng) ------------
cat("\n--- Thống kê theo nhóm target (thang đo gốc) ---\n")
for (col in numeric_cols) {
  s <- tapply(df_raw[[col]], df_raw$target, function(x) {
    x <- x[!is.na(x)]
    c(min = min(x), max = max(x), mean = mean(x), median = median(x), sd = sd(x))
  })
  cat("\n", col, "\n")
  print(round(do.call(rbind, s), 2))
}


# ---- 2.8 Mann-Whitney U + hiệu chỉnh FDR ------------------------------------
# Dùng Mann-Whitney thay t-test vì phần lớn biến lệch phải nặng, mà t-test giả
# định phân phối chuẩn. Kiểm định này dựa trên thứ hạng nên p-value không đổi
# dù đã log1p ở trên.
fs <- data.frame(feature = character(), p_value = numeric(),
                 effect_size = numeric(), stringsAsFactors = FALSE)

for (col in numeric_cols) {
  g0 <- df[[col]][df$target == 0]
  g1 <- df[[col]][df$target == 1]
  w <- wilcox.test(g0, g1, exact = FALSE)
  # Rank-biserial correlation: effect size cho Mann-Whitney, nằm trong [-1, 1]
  n0 <- sum(!is.na(g0)); n1 <- sum(!is.na(g1))
  rb <- 1 - 2 * w$statistic / (n0 * n1)
  fs <- rbind(fs, data.frame(feature = col, p_value = w$p.value,
                             effect_size = as.numeric(rb)))
}

# Gender: dùng Cramer's V làm effect size
n_total <- sum(contingency)
cramers_v <- sqrt(chi_res$statistic / (n_total * (min(dim(contingency)) - 1)))
fs <- rbind(fs, data.frame(feature = "Gender", p_value = chi_res$p.value,
                           effect_size = as.numeric(cramers_v)))

# Hiệu chỉnh Benjamini-Hochberg vì chạy 10 kiểm định cùng lúc ở alpha = 0.05
fs$p_fdr <- p.adjust(fs$p_value, method = "BH")
fs$significant <- fs$p_fdr < 0.05
fs <- fs[order(fs$p_value), ]

cat("\n--- Feature selection (Mann-Whitney + FDR) ---\n")
print(fs, row.names = FALSE)
cat("\nGiữ lại:", paste(fs$feature[fs$significant], collapse = ", "), "\n")
cat("Loại bỏ:", paste(fs$feature[!fs$significant], collapse = ", "), "\n")


# =============================================================================
# PHẦN 3. CÁC HÀM DÙNG CHUNG CHO PHẦN MÔ HÌNH
# =============================================================================

# Bỏ 3 biến:
#   Direct_Bilirubin  - trùng thông tin với Total_Bilirubin (r = 0.87, VIF ~ 20)
#   Gender, Total_Protiens - FDR p >= 0.05, không có quan hệ với target
model_cols <- c("Age", "Total_Bilirubin", "Alkaline_Phosphotase",
                "Alamine_Aminotransferase", "Aspartate_Aminotransferase",
                "Albumin", "Albumin_and_Globulin_Ratio")

X_all <- as.matrix(df[, model_cols])
y_all <- df$target


# ---- 3.1 Chia phân tầng -----------------------------------------------------
# Giữ nguyên tỉ lệ 71/29 của target ở cả hai tập.
stratified_split <- function(y, test_frac = 0.2, seed = 42) {
  set.seed(seed)
  test_idx <- integer(0)
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    test_idx <- c(test_idx, sample(idx, round(length(idx) * test_frac)))
  }
  list(train = setdiff(seq_along(y), test_idx), test = sort(test_idx))
}


# ---- 3.2 Gán fold phân tầng -------------------------------------------------
# Trả về vector số hiệu fold (1..k) cho từng quan sát. Xáo trong nội bộ từng lớp
# rồi chia đều nên mỗi fold giữ đúng tỉ lệ target.
stratified_folds <- function(y, k = 5, seed = 42) {
  set.seed(seed)
  fold <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- sample(which(y == cls))
    fold[idx] <- rep_len(seq_len(k), length(idx))
  }
  fold
}


# ---- 3.3 SMOTE tự cài đặt ---------------------------------------------------
# Sinh mẫu tổng hợp cho lớp thiểu số bằng nội suy tuyến tính giữa một mẫu và
# một trong k láng giềng gần nhất CÙNG LỚP của nó, cho tới khi hai lớp cân bằng.
smote_oversample <- function(X, y, k = 5, seed = 42) {
  set.seed(seed)
  tab <- table(y)
  min_lab <- as.numeric(names(tab)[which.min(tab)])
  n_new <- as.integer(max(tab) - min(tab))
  if (n_new <= 0) return(list(X = X, y = y))

  X_min <- X[y == min_lab, , drop = FALSE]
  if (nrow(X_min) <= k) k <- max(1, nrow(X_min) - 1)

  # Ma trận khoảng cách trong nội bộ lớp thiểu số; chặn đường chéo để không tự chọn mình
  D <- as.matrix(dist(X_min))
  diag(D) <- Inf

  synth <- matrix(0, nrow = n_new, ncol = ncol(X))
  for (i in seq_len(n_new)) {
    a <- sample.int(nrow(X_min), 1)
    neighbours <- order(D[a, ])[seq_len(k)]      # k láng giềng gần nhất
    b <- neighbours[sample.int(k, 1)]            # chọn ngẫu nhiên 1 trong k
    gap <- runif(1)                              # điểm mới nằm trên đoạn a-b
    synth[i, ] <- X_min[a, ] + gap * (X_min[b, ] - X_min[a, ])
  }
  colnames(synth) <- colnames(X)
  list(X = rbind(X, synth), y = c(y, rep(min_lab, n_new)))
}


# ---- 3.4 Hồi quy logistic có điều chuẩn ridge (IRLS) ------------------------
# glm() của base R không hỗ trợ phạt L2, mà glmnet là gói biên dịch nên không
# dùng được ở đây. Tự cài IRLS: mỗi vòng lặp giải một bài toán bình phương tối
# thiểu có trọng số, cộng thêm ma trận phạt lambda (không phạt hệ số chặn).
fit_logreg_ridge <- function(X, y, lambda = 0, max_iter = 50, tol = 1e-8) {
  Xd <- cbind(1, X)                       # thêm cột hệ số chặn
  p_dim <- ncol(Xd)
  beta <- rep(0, p_dim)

  P <- diag(lambda, p_dim)
  P[1, 1] <- 0                            # không điều chuẩn hệ số chặn

  for (it in seq_len(max_iter)) {
    eta <- as.vector(Xd %*% beta)
    mu <- 1 / (1 + exp(-eta))             # xác suất dự đoán
    W <- mu * (1 - mu)                    # trọng số IRLS
    W[W < 1e-10] <- 1e-10                 # chặn dưới để ma trận không suy biến
    z <- eta + (y - mu) / W               # biến phản hồi hiệu chỉnh

    beta_new <- tryCatch(
      solve(t(Xd) %*% (Xd * W) + P, t(Xd) %*% (W * z)),
      error = function(e) beta
    )
    if (max(abs(beta_new - beta)) < tol) { beta <- beta_new; break }
    beta <- as.vector(beta_new)
  }
  as.vector(beta)
}

predict_proba <- function(beta, X) {
  eta <- as.vector(cbind(1, X) %*% beta)
  1 / (1 + exp(-eta))
}


# ---- 3.5 Quy trình huấn luyện hoàn chỉnh ------------------------------------
# Gói impute -> scale -> SMOTE -> ridge logistic vào một hàm. Cả ba bước đầu đều
# HỌC THAM SỐ từ dữ liệu, nên phải nằm trong đây để chỉ nhìn thấy tập huấn luyện.
# Nếu chạy chúng trên toàn bộ dữ liệu trước khi chia fold thì tập kiểm tra bị rò rỉ.
fit_pipeline <- function(X_tr, y_tr, lambda = 1, smote_k = 5, seed = 42) {
  # 1. Điền khuyết bằng median CỦA TẬP TRAIN
  med <- apply(X_tr, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(X_tr))) X_tr[is.na(X_tr[, j]), j] <- med[j]

  # 2. Chuẩn hoá theo mean/sd CỦA TẬP TRAIN
  mu <- colMeans(X_tr)
  sdv <- apply(X_tr, 2, sd)
  sdv[sdv == 0] <- 1
  X_tr_s <- scale(X_tr, center = mu, scale = sdv)

  # 3. SMOTE chỉ trên tập train
  res <- smote_oversample(X_tr_s, y_tr, k = smote_k, seed = seed)

  # 4. Khớp mô hình
  beta <- fit_logreg_ridge(res$X, res$y, lambda = lambda)

  list(beta = beta, median = med, mean = mu, sd = sdv)
}

# Áp đúng chuỗi biến đổi đã học từ train lên dữ liệu mới
pipeline_proba <- function(fit, X_new) {
  for (j in seq_len(ncol(X_new))) X_new[is.na(X_new[, j]), j] <- fit$median[j]
  X_s <- scale(X_new, center = fit$mean, scale = fit$sd)
  predict_proba(fit$beta, X_s)
}


# ---- 3.6 Bộ chỉ số ----------------------------------------------------------
compute_metrics <- function(y_true, y_pred) {
  tp <- sum(y_true == 1 & y_pred == 1)
  tn <- sum(y_true == 0 & y_pred == 0)
  fp <- sum(y_true == 0 & y_pred == 1)
  fn <- sum(y_true == 1 & y_pred == 0)

  recall1 <- if ((tp + fn) > 0) tp / (tp + fn) else NA      # độ nhạy
  recall0 <- if ((tn + fp) > 0) tn / (tn + fp) else NA      # độ đặc hiệu
  prec1 <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  prec0 <- if ((tn + fn) > 0) tn / (tn + fn) else 0
  f1_1 <- if ((prec1 + recall1) > 0) 2 * prec1 * recall1 / (prec1 + recall1) else 0
  f1_0 <- if ((prec0 + recall0) > 0) 2 * prec0 * recall0 / (prec0 + recall0) else 0

  c(accuracy = (tp + tn) / length(y_true),
    recall1 = recall1, precision1 = prec1,
    recall0 = recall0, precision0 = prec0,
    macro_f1 = (f1_1 + f1_0) / 2)
}

print_confusion <- function(y_true, y_pred, title = "CONFUSION MATRIX") {
  tp <- sum(y_true == 1 & y_pred == 1); tn <- sum(y_true == 0 & y_pred == 0)
  fp <- sum(y_true == 0 & y_pred == 1); fn <- sum(y_true == 1 & y_pred == 0)

  cat("\n", title, "\n", sep = "")
  print(matrix(c(tn, fp, fn, tp), nrow = 2, byrow = TRUE,
               dimnames = list(`Thực tế` = c("Không bệnh", "Có bệnh"),
                               `Dự đoán` = c("Không bệnh", "Có bệnh"))))
  cat(sprintf("  TN = %3d   không bệnh, đoán đúng\n", tn))
  cat(sprintf("  FP = %3d   không bệnh nhưng báo có bệnh (báo động nhầm)\n", fp))
  cat(sprintf("  FN = %3d   CÓ BỆNH nhưng bỏ sót        <-- sai lầm tốn kém nhất\n", fn))
  cat(sprintf("  TP = %3d   có bệnh, bắt được\n", tp))
  invisible(compute_metrics(y_true, y_pred))
}


# =============================================================================
# PHẦN 4. BOOTSTRAP - KHOẢNG TIN CẬY CHO RECALL LỚP CÓ BỆNH
# =============================================================================
# Recall chỉ là MỘT con số tính từ hơn trăm bệnh nhân trong tập test. Lấy một
# nhóm khác thì con số đó khác. Bootstrap ước lượng mức dao động ấy.

sp <- stratified_split(y_all, test_frac = 0.2, seed = 42)
X_train <- X_all[sp$train, , drop = FALSE]; y_train <- y_all[sp$train]
X_test  <- X_all[sp$test,  , drop = FALSE]; y_test  <- y_all[sp$test]

cat("\n\n=== PHẦN 4: BOOTSTRAP ===\n")
cat("Train:", nrow(X_train), "mẫu | Test:", nrow(X_test), "mẫu\n")

# Huấn luyện MỘT lần, dự đoán MỘT lần. Mô hình cố định suốt phần bootstrap.
fit_80 <- fit_pipeline(X_train, y_train, lambda = 1, smote_k = 5)
proba_test <- pipeline_proba(fit_80, X_test)
pred_test <- as.integer(proba_test >= 0.5)

recall_point <- compute_metrics(y_test, pred_test)["recall1"]

CONF <- 0.90        # mức tin cậy - đổi mỗi số này là đủ
B <- 10000
n_test <- length(y_test)

recalls <- numeric(0)
for (b in seq_len(B)) {
  idx <- sample.int(n_test, n_test, replace = TRUE)   # rút CÓ HOÀN LẠI, cùng cỡ mẫu
  if (sum(y_test[idx] == 1) == 0) next                # không có ca bệnh -> recall vô nghĩa
  recalls <- c(recalls, compute_metrics(y_test[idx], pred_test[idx])["recall1"])
}

# Khoảng tin cậy percentile: cắt (1-CONF)/2 ở MỖI đuôi
tail_pct <- (1 - CONF) / 2
ci <- quantile(recalls, c(tail_pct, 1 - tail_pct))

cat(sprintf("\nRecall lớp có bệnh   = %.3f\n", recall_point))
cat(sprintf("Khoảng tin cậy %.0f%%   = [%.3f, %.3f]  (phân vị %.1f%% và %.1f%%)\n",
            CONF * 100, ci[1], ci[2], tail_pct * 100, (1 - tail_pct) * 100))
cat(sprintf("Bề rộng khoảng       = %.3f\n", ci[2] - ci[1]))
cat(sprintf("Sai số chuẩn         = %.3f\n", sd(recalls)))

hist(recalls, breaks = 40, col = "#66c2a5", border = "white",
     main = sprintf("Phân phối bootstrap của Recall (B = %d)", length(recalls)),
     xlab = "Recall lớp có bệnh", ylab = "Số lần xuất hiện")
abline(v = recall_point, col = "#d53e4f", lwd = 2)
abline(v = ci, col = "gray40", lty = 2, lwd = 1.5)
legend("topright", bty = "n", cex = 0.85,
       legend = c(sprintf("Điểm ước lượng = %.3f", recall_point),
                  sprintf("KTC %.0f%%", CONF * 100)),
       col = c("#d53e4f", "gray40"), lty = c(1, 2), lwd = 2)


# =============================================================================
# PHẦN 5. STRATIFIED K-FOLD (k = 5)
# =============================================================================
# Chia dữ liệu thành 5 phần bằng nhau, lần lượt lấy một phần làm tập kiểm tra.
# Mỗi quan sát được dự đoán đúng một lần.
#
# Chọn k = 5 vì mỗi fold kiểm tra đúng 20% dữ liệu, TRÙNG với tỉ lệ 80:20 ở
# Phần 4 -> hai phần cùng đánh giá mô hình huấn luyện trên ~466 mẫu.

cat("\n\n=== PHẦN 5: STRATIFIED 5-FOLD ===\n")

K <- 5
folds <- stratified_folds(y_all, k = K, seed = 42)

cv_results <- data.frame()
for (i in seq_len(K)) {
  te <- which(folds == i)
  tr <- which(folds != i)

  fit_i <- fit_pipeline(X_all[tr, , drop = FALSE], y_all[tr],
                        lambda = 1, smote_k = 5)
  pred_i <- as.integer(pipeline_proba(fit_i, X_all[te, , drop = FALSE]) >= 0.5)

  m <- compute_metrics(y_all[te], pred_i)
  cv_results <- rbind(cv_results, data.frame(
    Fold = i, n_train = length(tr), n_test = length(te),
    Accuracy = m["accuracy"], Recall = m["recall1"], Macro_F1 = m["macro_f1"]
  ))
}
rownames(cv_results) <- NULL

cat("\n")
print(round(cv_results, 4), row.names = FALSE)

# Báo cáo dạng mean +/- sd, không phải một con số đơn lẻ.
# sd cho biết kết quả dao động bao nhiêu khi đổi sang tập kiểm tra khác.
cat(sprintf("\nKết quả %d-Fold (mean +/- sd giữa các fold):\n", K))
for (m in c("Accuracy", "Recall", "Macro_F1")) {
  cat(sprintf("  %-9s = %.3f +/- %.3f\n",
              m, mean(cv_results[[m]]), sd(cv_results[[m]])))
}


# =============================================================================
# PHẦN 6. GRID SEARCH + 5-FOLD TRÊN TẬP TRAIN
# =============================================================================
# Quy trình chuẩn: K-Fold để TINH CHỈNH THAM SỐ trên tập train, Bootstrap để
# ước lượng KHOẢNG TIN CẬY trên tập test. Tập test bị niêm phong tới Phần 8.
#
# Chấm điểm bằng macro F1, KHÔNG dùng recall thuần: tối đa hoá recall lớp 1 có
# nghiệm tầm thường là dự đoán tất cả "có bệnh" -> recall = 1 nhưng vô dụng.

cat("\n\n=== PHẦN 6: GRID SEARCH ===\n")

param_grid <- expand.grid(
  lambda = c(0.01, 0.1, 1, 10, 100),   # cường độ điều chuẩn ridge
  smote_k = c(3, 5, 7)                 # số láng giềng SMOTE dùng để nội suy
)

folds_inner <- stratified_folds(y_train, k = 5, seed = 42)

grid_mean <- grid_sd <- numeric(nrow(param_grid))
for (g in seq_len(nrow(param_grid))) {
  scores <- numeric(5)
  for (i in seq_len(5)) {
    te <- which(folds_inner == i); tr <- which(folds_inner != i)
    fit_g <- fit_pipeline(X_train[tr, , drop = FALSE], y_train[tr],
                          lambda = param_grid$lambda[g],
                          smote_k = param_grid$smote_k[g])
    pred_g <- as.integer(pipeline_proba(fit_g, X_train[te, , drop = FALSE]) >= 0.5)
    scores[i] <- compute_metrics(y_train[te], pred_g)["macro_f1"]
  }
  grid_mean[g] <- mean(scores)
  grid_sd[g] <- sd(scores)          # dao động giữa 5 fold, dùng để đo mức nhiễu
}

param_grid$macro_f1 <- grid_mean
param_grid$sd_fold <- grid_sd

# ---- Quy tắc phá thế hoà -----------------------------------------------------
# Lấy thẳng which.max là sai lầm khi các bộ tham số cho điểm chênh nhau ít hơn
# chính sai số của phép đo: lúc đó "bộ thắng" chỉ là bộ gặp may với cách chia fold
# cụ thể này, chọn nó không có cơ sở thống kê nào.
#
# Cách xử lý: coi mọi bộ nằm trong PHẠM VI MỘT SAI SỐ CHUẨN quanh điểm cao nhất
# là ngang tài nhau, rồi phá hoà bằng một tiêu chí có lý do rõ ràng.
#
# Ở đây tiêu chí là CHỌN LAMBDA NHỎ NHẤT. Lưu ý điều này NGƯỢC với quy tắc 1-SE
# quen thuộc (vốn ưu tiên mô hình điều chuẩn mạnh hơn cho gọn). Lý do đảo chiều:
# bước tiếp theo của quy trình là chỉnh ngưỡng phân loại, mà việc đó cần xác suất
# dự đoán trải rộng. Lambda lớn nén hệ số lại, kéo xác suất co về quanh 0.5, khiến
# ngưỡng trở nên cực kỳ nhạy - dịch một chút là recall của hai lớp đảo lộn.
# Khi các bộ đã ngang điểm nhau, chọn bộ cho xác suất trải rộng hơn là hợp lý.
i_top <- which.max(param_grid$macro_f1)
se_top <- param_grid$sd_fold[i_top] / sqrt(5)              # sai số chuẩn của trung bình
tied <- which(param_grid$macro_f1 >= param_grid$macro_f1[i_top] - se_top)

cand <- tied[param_grid$lambda[tied] == min(param_grid$lambda[tied])]
best_idx <- cand[which.max(param_grid$macro_f1[cand])]     # cùng lambda thì lấy điểm cao nhất
best <- param_grid[best_idx, ]

cat("\nTop 5 bộ tham số:\n")
print(head(param_grid[order(-param_grid$macro_f1), ], 5), row.names = FALSE)

cat(sprintf("\nĐiểm cao nhất: %.4f (lambda = %g), sai số chuẩn = %.4f\n",
            param_grid$macro_f1[i_top], param_grid$lambda[i_top], se_top))
cat(sprintf("Có %d/%d bộ nằm trong phạm vi một sai số chuẩn -> coi là ngang nhau\n",
            length(tied), nrow(param_grid)))
cat(sprintf("Phá hoà bằng lambda nhỏ nhất -> chọn lambda = %g, smote_k = %g (macro F1 = %.4f)\n",
            best$lambda, best$smote_k, best$macro_f1))
cat(sprintf("Đã thử %d bộ x 5 fold = %d lần khớp mô hình\n",
            nrow(param_grid), nrow(param_grid) * 5))


# =============================================================================
# PHẦN 7. CHỌN NGƯỠNG PHÂN LOẠI
# =============================================================================
# predict mặc định cắt ở 0.5, tối ưu cho độ chính xác tổng thể chứ không ưu tiên
# lớp nào. Bài toán bệnh tật thì bỏ sót ca bệnh (FN) tốn kém hơn báo động nhầm
# (FP), nên hạ ngưỡng xuống để đẩy recall lớp có bệnh lên.
#
# NGƯỠNG PHẢI DÒ TRÊN TẬP TRAIN. Dò trên tập test thì tập test đã tham gia vào
# một quyết định, con số báo cáo sẽ lạc quan giả tạo.

cat("\n\n=== PHẦN 7: CHỌN NGƯỠNG ===\n")

# Dự đoán out-of-fold: mỗi mẫu train được chấm bởi mô hình chưa từng thấy nó
oof_proba <- numeric(length(y_train))
for (i in seq_len(5)) {
  te <- which(folds_inner == i); tr <- which(folds_inner != i)
  fit_o <- fit_pipeline(X_train[tr, , drop = FALSE], y_train[tr],
                        lambda = best$lambda, smote_k = best$smote_k)
  oof_proba[te] <- pipeline_proba(fit_o, X_train[te, , drop = FALSE])
}

# Quét ngưỡng, theo dõi recall CẢ HAI LỚP:
#   recall lớp 1 = độ nhạy (sensitivity)    - bắt được bao nhiêu phần ca bệnh
#   recall lớp 0 = độ đặc hiệu (specificity) - loại đúng bao nhiêu phần người khoẻ
grid_thr <- seq(0.05, 0.95, by = 0.005)
rec1 <- rec0 <- numeric(length(grid_thr))
for (i in seq_along(grid_thr)) {
  pred_t <- as.integer(oof_proba >= grid_thr[i])
  m <- compute_metrics(y_train, pred_t)
  rec1[i] <- m["recall1"]; rec0[i] <- m["recall0"]
}

plot(grid_thr, rec1, type = "l", col = "#d53e4f", lwd = 2, ylim = c(0, 1),
     xlab = "Ngưỡng phân loại", ylab = "Recall",
     main = "Recall của hai lớp theo ngưỡng (out-of-fold trên tập train)")
lines(grid_thr, rec0, col = "#3288bd", lwd = 2)
abline(v = 0.5, col = "gray50", lty = 3)
grid(col = "gray90")
legend("right", bty = "n", cex = 0.85, lwd = 2,
       legend = c("Recall lớp 1 - có bệnh (độ nhạy)",
                  "Recall lớp 0 - không bệnh (độ đặc hiệu)",
                  "Ngưỡng mặc định 0.5"),
       col = c("#d53e4f", "#3288bd", "gray50"), lty = c(1, 1, 3))

# Bảng tra: muốn đạt mức recall lớp 1 nào thì phải hạ ngưỡng tới đâu,
# và recall lớp 0 phải trả giá bao nhiêu
cat("\nBảng tra ngưỡng (out-of-fold trên tập train):\n")
lookup <- data.frame()
for (target in c(0.70, 0.75, 0.80, 0.85, 0.90, 0.95)) {
  ok <- which(rec1 >= target)
  if (length(ok) == 0) next
  i <- ok[which.max(grid_thr[ok])]     # ngưỡng CAO NHẤT vẫn đạt mục tiêu
  lookup <- rbind(lookup, data.frame(
    recall_muc_tieu = target, nguong = grid_thr[i],
    recall_lop1 = rec1[i], recall_lop0 = rec0[i]
  ))
}
print(round(lookup, 3), row.names = FALSE)

TARGET_RECALL <- 0.70        # mốc chọn dựa trên biểu đồ và bảng tra ở trên
ok <- which(rec1 >= TARGET_RECALL)
threshold <- grid_thr[ok[which.max(grid_thr[ok])]]

cat(sprintf("\nNgưỡng mặc định                      : 0.500\n"))
cat(sprintf("Ngưỡng chọn từ out-of-fold trên train: %.3f  (mục tiêu recall >= %.2f)\n",
            threshold, TARGET_RECALL))


# =============================================================================
# PHẦN 8. ĐÁNH GIÁ TRÊN TẬP TEST + BOOTSTRAP KHOẢNG TIN CẬY
# =============================================================================
# Đây là lần ĐẦU TIÊN tập test được đụng tới sau khi niêm phong ở Phần 6.

cat("\n\n=== PHẦN 8: ĐÁNH GIÁ CUỐI CÙNG ===\n")

fit_best <- fit_pipeline(X_train, y_train,
                         lambda = best$lambda, smote_k = best$smote_k)
proba_final <- pipeline_proba(fit_best, X_test)
pred_final <- as.integer(proba_final >= threshold)
pred_05 <- as.integer(proba_final >= 0.5)

m_final <- print_confusion(y_test, pred_final,
                           sprintf("CONFUSION MATRIX (ngưỡng %.3f)", threshold))

cat("\nCác chỉ số trên tập test:\n")
print(round(m_final, 3))

# Cái giá phải trả khi hạ ngưỡng
cat("\nĐánh đổi khi hạ ngưỡng:\n")
cmp <- rbind(
  `0.500 (mặc định)` = compute_metrics(y_test, pred_05),
  `đã chỉnh`         = compute_metrics(y_test, pred_final)
)
print(round(cmp, 3))

# ---- Bootstrap khoảng tin cậy cho các chỉ số trên tập test ------------------
CONF2 <- 0.95
B2 <- 2000
boot_mat <- matrix(NA, nrow = B2, ncol = length(m_final),
                   dimnames = list(NULL, names(m_final)))

for (b in seq_len(B2)) {
  idx <- sample.int(n_test, n_test, replace = TRUE)
  if (length(unique(y_test[idx])) < 2) next        # mẫu suy biến, bỏ qua
  boot_mat[b, ] <- compute_metrics(y_test[idx], pred_final[idx])
}
boot_mat <- boot_mat[complete.cases(boot_mat), , drop = FALSE]

tail2 <- (1 - CONF2) / 2
ci_tbl <- data.frame(
  diem_uoc_luong = m_final,
  ci_duoi = apply(boot_mat, 2, quantile, probs = tail2),
  ci_tren = apply(boot_mat, 2, quantile, probs = 1 - tail2),
  sai_so_chuan = apply(boot_mat, 2, sd)
)

cat(sprintf("\nBootstrap B = %d, mức tin cậy %.0f%%\n", nrow(boot_mat), CONF2 * 100))
print(round(ci_tbl, 4))

cat(sprintf("\nKết quả chính: Recall lớp có bệnh = %.3f, KTC %.0f%% [%.3f, %.3f]\n",
            ci_tbl["recall1", "diem_uoc_luong"], CONF2 * 100,
            ci_tbl["recall1", "ci_duoi"], ci_tbl["recall1", "ci_tren"]))

# Forest plot: điểm ước lượng kèm thanh khoảng tin cậy
par(mar = c(5, 9, 4, 2))
n_m <- nrow(ci_tbl)
plot(ci_tbl$diem_uoc_luong, seq_len(n_m), pch = 19, col = "#3288bd", cex = 1.3,
     xlim = c(0, 1), yaxt = "n", ylab = "", xlab = "Giá trị",
     main = sprintf("Khoảng tin cậy bootstrap %.0f%%", CONF2 * 100))
arrows(ci_tbl$ci_duoi, seq_len(n_m), ci_tbl$ci_tren, seq_len(n_m),
       length = 0.05, angle = 90, code = 3, col = "gray40")
axis(2, seq_len(n_m), rownames(ci_tbl), las = 2, cex.axis = 0.85)
abline(v = 0.5, col = "#d53e4f", lty = 3)
par(mar = c(5, 4, 4, 2) + 0.1)

cat("\n=== HOÀN TẤT ===\n")

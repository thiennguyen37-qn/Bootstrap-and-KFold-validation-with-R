## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE,
                      fig.align = "center", fig.width = 6, fig.height = 4,
                      fig.path = "figure/tl-", dev = "cairo_pdf",
                      fig.pos = "ht", size = "scriptsize")
options(width = 72)
set.seed(42)

if (!interactive()) {
  pdf_ok <- try(
    grDevices::cairo_pdf("Rplots.pdf", width = 8, height = 6, onefile = TRUE),
    silent = TRUE
  )
  if (inherits(pdf_ok, "try-error")) {
    grDevices::pdf("Rplots.pdf", width = 8, height = 6, onefile = TRUE)
  }
  on.exit({
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
}


## ----load-data----------------------------------------------------------------
col_names <- c(
  "Age", "Gender", "Total_Bilirubin", "Direct_Bilirubin",
  "Alkaline_Phosphotase", "Alamine_Aminotransferase",
  "Aspartate_Aminotransferase", "Total_Protiens", "Albumin",
  "Albumin_and_Globulin_Ratio", "target"
)

df <- read.csv("Indian Liver Patient Dataset (ILPD).csv",
               header = FALSE, col.names = col_names,
               stringsAsFactors = FALSE)

# Nhãn gốc: 1 = có bệnh, 2 = không bệnh. Đổi về 1 / 0.
df$target <- ifelse(df$target == 2, 0, 1)

cat("Kích thước dữ liệu:", nrow(df), "dòng x", ncol(df), "cột\n")
cat("Phân bố target:", sum(df$target == 1), "có bệnh /",
    sum(df$target == 0), "không bệnh\n")


## ----data-structure-----------------------------------------------------------
str(df)


## ----summary-missing----------------------------------------------------------
summary(df)
colSums(is.na(df))


## ----target-dist, echo=FALSE, fig.width=6, fig.height=4.2, fig.cap="Phân bố biến mục tiêu"----
counts <- table(df$target)
pct <- round(100 * counts / sum(counts), 1)

bp <- barplot(counts,
              names.arg = c("Không bệnh (0)", "Có bệnh (1)"),
              col = c("#66c2a5", "#fc8d62"),
              main = "Phân bố biến mục tiêu", ylab = "Số lượng",
              ylim = c(0, max(counts) * 1.15))
text(bp, counts, labels = paste0(counts, "\n(", pct, "%)"), pos = 3)


## ----corr-matrix, echo=FALSE--------------------------------------------------
numeric_cols <- c(
  "Age", "Total_Bilirubin", "Direct_Bilirubin", "Alkaline_Phosphotase",
  "Alamine_Aminotransferase", "Aspartate_Aminotransferase",
  "Total_Protiens", "Albumin", "Albumin_and_Globulin_Ratio"
)

corr_mat <- cor(df[, c(numeric_cols, "target")], use = "complete.obs")


## ----corr-heatmap, echo=FALSE, fig.width=6.5, fig.height=6, fig.cap="Ma trận tương quan giữa các biến sinh hoá và biến mục tiêu"----
n_c <- nrow(corr_mat)
z <- t(corr_mat[n_c:1, ])   # đảo hàng để trục y đọc từ trên xuống

par(mar = c(9, 9, 3, 2))
image(seq_len(n_c), seq_len(n_c), z,
      col = colorRampPalette(c("#3288bd", "white", "#d53e4f"))(64),
      zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "",
      main = "Ma trận tương quan")
axis(1, seq_len(n_c), colnames(corr_mat), las = 2, cex.axis = 0.7)
axis(2, seq_len(n_c), rev(rownames(corr_mat)), las = 2, cex.axis = 0.7)

# Ghi hệ số r lên từng ô; chữ trắng trên ô đậm, chữ đen trên ô nhạt
xs <- rep(seq_len(n_c), times = n_c)
ys <- rep(seq_len(n_c), each = n_c)
vals <- z[cbind(xs, ys)]
text(xs, ys, sprintf("%.2f", vals), cex = 0.55,
     col = ifelse(abs(vals) > 0.55, "white", "black"))
par(mar = c(5, 4, 4, 2) + 0.1)


## ----boxplot-raw, echo=FALSE, fig.pos="H", fig.width=8, fig.height=7.5, fig.cap="Phân bố các biến sinh hoá theo nhóm bệnh trên thang đo gốc"----
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
for (col in numeric_cols) {
  boxplot(df[[col]] ~ df$target, col = c("#66c2a5", "#fc8d62"),
          main = col, xlab = "target", ylab = "")
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)


## ----log1p-transform, echo=FALSE, results="hide"------------------------------
skewed_cols <- c("Total_Bilirubin", "Direct_Bilirubin",
                 "Alkaline_Phosphotase", "Alamine_Aminotransferase",
                 "Aspartate_Aminotransferase")

df_raw <- df                       # giữ bản gốc để tra thang đo lâm sàng
for (col in skewed_cols) df[[col]] <- log1p(df[[col]])

# Hệ số bất đối xứng (skewness)
skewness <- function(x) {
  x <- x[!is.na(x)]
  mean((x - mean(x))^3) / (mean((x - mean(x))^2))^1.5
}

print(round(data.frame(
  goc       = sapply(df_raw[skewed_cols], skewness),
  sau_log1p = sapply(df[skewed_cols], skewness)
), 2))


## ----boxplot-log, echo=FALSE, fig.pos="H", fig.width=8, fig.height=7.5, fig.cap="Phân bố các biến sinh hoá theo nhóm bệnh sau biến đổi log"----
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
for (col in numeric_cols) {
  boxplot(df[[col]] ~ df$target, col = c("#66c2a5", "#fc8d62"),
          main = if (col %in% skewed_cols) paste0("log1p(", col, ")") else col,
          xlab = "target", ylab = "")
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)


## ----vif-compare--------------------------------------------------------------
compute_vif <- function(data) {
  vars <- colnames(data)
  sapply(vars, function(v) {
    others <- setdiff(vars, v)
    fml <- as.formula(paste(v, "~", paste(others, collapse = " + ")))
    r2 <- summary(lm(fml, data = data))$r.squared
    1 / (1 - r2)
  })
}

# Thang gốc: dùng df_raw (bản trước biến đổi log)
X_vif <- df_raw[, c(numeric_cols, "Gender")]
X_vif$Gender <- ifelse(X_vif$Gender == "Male", 1, 0)
X_vif$Albumin_and_Globulin_Ratio[is.na(X_vif$Albumin_and_Globulin_Ratio)] <-
  median(X_vif$Albumin_and_Globulin_Ratio, na.rm = TRUE)
vif_raw <- compute_vif(X_vif)

# Thang log: log1p cho các biến lệch phải
X_vif_log <- X_vif
for (col in skewed_cols) X_vif_log[[col]] <- log1p(df_raw[[col]])
vif_log <- compute_vif(X_vif_log)

print(round(data.frame(
  VIF_goc    = vif_raw,
  VIF_log    = vif_log[names(vif_raw)],
  chenh_lech = vif_log[names(vif_raw)] - vif_raw
)[order(-vif_log[names(vif_raw)]), ], 2))


## ----group-stats, size="tiny"-------------------------------------------------
for (col in numeric_cols) {
  s <- tapply(df_raw[[col]], df_raw$target, function(x) {
    x <- x[!is.na(x)]
    c(min = min(x), max = max(x), mean = mean(x),
      median = median(x), sd = sd(x))
  })
  cat("\n", col, "\n")
  print(round(do.call(rbind, s), 2))
}


## ----feature-select-----------------------------------------------------------
fs <- data.frame(feature = character(), p_value = numeric(),
                 effect_size = numeric(), stringsAsFactors = FALSE)

for (col in numeric_cols) {
  g0 <- df[[col]][df$target == 0]
  g1 <- df[[col]][df$target == 1]
  w <- wilcox.test(g0, g1, exact = FALSE)
  n0 <- sum(!is.na(g0)); n1 <- sum(!is.na(g1))
  rb <- 1 - 2 * w$statistic / (n0 * n1)     # rank-biserial correlation
  fs <- rbind(fs, data.frame(feature = col, p_value = w$p.value,
                             effect_size = as.numeric(rb)))
}

# Gender là biến phân loại: dùng Chi-square, kích thước ảnh hưởng Cramer's V
contingency <- table(df$Gender, df$target)
chi_res <- chisq.test(contingency)
n_total <- sum(contingency)
cramers_v <- sqrt(chi_res$statistic /
                    (n_total * (min(dim(contingency)) - 1)))
fs <- rbind(fs, data.frame(feature = "Gender", p_value = chi_res$p.value,
                           effect_size = as.numeric(cramers_v)))

fs$p_fdr <- p.adjust(fs$p_value, method = "BH")   # Benjamini-Hochberg
fs$significant <- fs$p_fdr < 0.05
fs <- fs[order(fs$p_value), ]

fs_disp <- fs
fs_disp$p_value <- signif(fs_disp$p_value, 3)
fs_disp$effect_size <- round(fs_disp$effect_size, 3)
fs_disp$p_fdr <- signif(fs_disp$p_fdr, 3)
print(fs_disp, row.names = FALSE)


## ----model-cols---------------------------------------------------------------
model_cols <- c("Age", "Total_Bilirubin", "Alkaline_Phosphotase",
                "Alamine_Aminotransferase", "Aspartate_Aminotransferase",
                "Albumin", "Albumin_and_Globulin_Ratio")

X_all <- as.matrix(df[, model_cols])
y_all <- df$target


## ----split-fold-fns-----------------------------------------------------------
stratified_split <- function(y, test_frac = 0.2, seed = 42) {
  set.seed(seed)
  test_idx <- integer(0)
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    test_idx <- c(test_idx, sample(idx, round(length(idx) * test_frac)))
  }
  list(train = setdiff(seq_along(y), test_idx), test = sort(test_idx))
}

stratified_folds <- function(y, k = 5, seed = 42) {
  set.seed(seed)
  fold <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- sample(which(y == cls))
    fold[idx] <- rep_len(seq_len(k), length(idx))
  }
  fold
}


## ----smote-fn-----------------------------------------------------------------
smote_oversample <- function(X, y, k = 5, seed = 42) {
  set.seed(seed)
  tab <- table(y)
  min_lab <- as.numeric(names(tab)[which.min(tab)])
  n_new <- as.integer(max(tab) - min(tab))
  if (n_new <= 0) return(list(X = X, y = y))

  X_min <- X[y == min_lab, , drop = FALSE]
  if (nrow(X_min) <= k) k <- max(1, nrow(X_min) - 1)

  # Khoảng cách trong nội bộ lớp thiểu số; chặn đường chéo để không tự chọn mình
  D <- as.matrix(dist(X_min))
  diag(D) <- Inf

  synth <- matrix(0, nrow = n_new, ncol = ncol(X))
  for (i in seq_len(n_new)) {
    a <- sample.int(nrow(X_min), 1)
    neighbours <- order(D[a, ])[seq_len(k)]
    b <- neighbours[sample.int(k, 1)]
    gap <- runif(1)
    synth[i, ] <- X_min[a, ] + gap * (X_min[b, ] - X_min[a, ])
  }
  colnames(synth) <- colnames(X)
  list(X = rbind(X, synth), y = c(y, rep(min_lab, n_new)))
}


## ----ridge-fn-----------------------------------------------------------------
fit_logreg_ridge <- function(X, y, lambda = 0, max_iter = 50, tol = 1e-8) {
  Xd <- cbind(1, X)                     # thêm cột hệ số chặn
  p_dim <- ncol(Xd)
  beta <- rep(0, p_dim)

  P <- diag(lambda, p_dim)
  P[1, 1] <- 0                          # không điều chuẩn hệ số chặn

  for (it in seq_len(max_iter)) {
    eta <- as.vector(Xd %*% beta)
    mu <- 1 / (1 + exp(-eta))           # xác suất dự đoán
    W <- mu * (1 - mu)                  # trọng số IRLS
    W[W < 1e-10] <- 1e-10               # chặn dưới để ma trận không suy biến
    z <- eta + (y - mu) / W             # biến phản hồi hiệu chỉnh

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
  1 / (1 + exp(-as.vector(cbind(1, X) %*% beta)))
}


## ----pipeline-fn--------------------------------------------------------------
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


## ----metrics-fn---------------------------------------------------------------
compute_metrics <- function(y_true, y_pred) {
  tp <- sum(y_true == 1 & y_pred == 1)
  tn <- sum(y_true == 0 & y_pred == 0)
  fp <- sum(y_true == 0 & y_pred == 1)
  fn <- sum(y_true == 1 & y_pred == 0)

  recall1 <- if ((tp + fn) > 0) tp / (tp + fn) else NA   # độ nhạy
  recall0 <- if ((tn + fp) > 0) tn / (tn + fp) else NA   # độ đặc hiệu
  prec1 <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  prec0 <- if ((tn + fn) > 0) tn / (tn + fn) else 0
  f1_1 <- if ((prec1 + recall1) > 0) 2*prec1*recall1/(prec1+recall1) else 0
  f1_0 <- if ((prec0 + recall0) > 0) 2*prec0*recall0/(prec0+recall0) else 0

  c(accuracy = (tp + tn) / length(y_true),
    recall1 = recall1, precision1 = prec1,
    recall0 = recall0, precision0 = prec0,
    macro_f1 = (f1_1 + f1_0) / 2)
}

print_confusion <- function(y_true, y_pred, title = "MA TRAN NHAM LAN") {
  tp <- sum(y_true == 1 & y_pred == 1); tn <- sum(y_true == 0 & y_pred == 0)
  fp <- sum(y_true == 0 & y_pred == 1); fn <- sum(y_true == 1 & y_pred == 0)

  cat("\n", title, "\n", sep = "")
  print(matrix(c(tn, fp, fn, tp), nrow = 2, byrow = TRUE,
               dimnames = list(`Thực tế` = c("Không bệnh", "Có bệnh"),
                               `Dự đoán` = c("Không bệnh", "Có bệnh"))))
  cat(sprintf("  TN = %3d   không bệnh, dự đoán đúng\n", tn))
  cat(sprintf("  FP = %3d   không bệnh, dự đoán nhầm là có bệnh\n", fp))
  cat(sprintf("  FN = %3d   có bệnh, bị bỏ sót\n", fn))
  cat(sprintf("  TP = %3d   có bệnh, phát hiện được\n", tp))
  invisible(compute_metrics(y_true, y_pred))
}


## ----python-reference---------------------------------------------------------
find_reference_script <- function() {
  candidates <- c("python_reference_pipeline.py",
                  file.path("..", "python_reference_pipeline.py"))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("Khong tim thay python_reference_pipeline.py.")
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

run_python_reference <- function(data_file) {
  py <- Sys.which("python")
  if (!nzchar(py)) stop("Khong tim thay lenh python tren PATH.")

  out_dir <- file.path(tempdir(), "ilpd_python_reference")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)

  cmd_out <- system2(
    py,
    c(shQuote(find_reference_script()), "--data", shQuote(data_file),
      "--out", shQuote(out_dir)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(cmd_out, "status")
  if (!is.null(status) && status != 0) {
    cat(paste(cmd_out, collapse = "\n"), "\n")
    stop("Pipeline Python tham chieu chay that bai.")
  }
  out_dir
}

ref_dir <- run_python_reference("Indian Liver Patient Dataset (ILPD).csv")
summary_ref <- jsonlite::fromJSON(file.path(ref_dir, "summary.json"))
pred_df <- read.csv(file.path(ref_dir, "test_predictions.csv"))
grid_top <- read.csv(file.path(ref_dir, "grid_top.csv"))
threshold_curve <- read.csv(file.path(ref_dir, "threshold_curve.csv"))
threshold_lookup <- read.csv(file.path(ref_dir, "threshold_lookup.csv"))
final_metrics <- read.csv(file.path(ref_dir, "final_metrics.csv"))
tradeoff <- read.csv(file.path(ref_dir, "tradeoff.csv"))
boot_samples <- read.csv(file.path(ref_dir, "bootstrap_samples.csv"))
ci_tbl <- read.csv(file.path(ref_dir, "bootstrap_ci.csv"))
rownames(ci_tbl) <- ci_tbl$metric
ci_tbl$metric <- NULL

y_test <- pred_df$y_test
pred_05 <- pred_df$pred_05
pred_final <- pred_df$pred_final
proba_final <- pred_df$proba_test
n_test <- length(y_test)
threshold <- summary_ref$threshold

cat("Train:", summary_ref$train_n, "mau | Test:", summary_ref$test_n,
    "mau (theo sklearn train_test_split)\n")
cat("Ca co benh - train:", summary_ref$train_pos,
    ", test:", summary_ref$test_pos, "\n")


## ----boot-hist, echo=FALSE, fig.width=6.5, fig.height=4.2, fig.cap="Phan phoi bootstrap cua recall lop co benh"----
recalls <- boot_samples$recall1
recall_point <- final_metrics$recall1
ci <- as.numeric(ci_tbl["recall1", c("ci_duoi", "ci_tren")])
CONF <- 0.95

cat(sprintf("Recall lop co benh   = %.3f\n", recall_point))
cat(sprintf("Khoang tin cay %.0f%%   = [%.3f, %.3f]\n",
            CONF * 100, ci[1], ci[2]))
cat(sprintf("Be rong khoang       = %.3f\n", ci[2] - ci[1]))
cat(sprintf("Sai so chuan         = %.3f\n", sd(recalls)))

hist(recalls, breaks = 40, col = "#66c2a5", border = "white",
     main = sprintf("Phan phoi bootstrap cua recall (B = %d)",
                    length(recalls)),
     xlab = "Recall lop co benh", ylab = "Tan so")
abline(v = recall_point, col = "#d53e4f", lwd = 2)
abline(v = ci, col = "gray40", lty = 2, lwd = 1.5)
legend("topright", bty = "n", cex = 0.85,
       legend = c(sprintf("Uoc luong diem = %.3f", recall_point),
                  sprintf("KTC %.0f%%", CONF * 100)),
       col = c("#d53e4f", "gray40"), lty = c(1, 2), lwd = 2)


## ----grid-search--------------------------------------------------------------
cat("\nBo tham so tot nhat theo notebook Python:\n")
for (nm in names(summary_ref$best_params)) {
  cat(sprintf("  %-20s = %s\n", nm, summary_ref$best_params[[nm]]))
}
cat(sprintf("\nDiem f1_macro (5-Fold tren train) = %.4f\n",
            summary_ref$best_score))
cat(sprintf("Tong so mo hinh da thu: %d bo x 5 fold = %d lan khop\n",
            summary_ref$n_param_sets, summary_ref$n_param_sets * 5))
print(grid_top, row.names = FALSE)


## ----threshold-plot, fig.width=6.8, fig.height=4.5, fig.cap="Danh doi giua do nhay va do dac hieu khi thay doi nguong phan loai"----
plot(threshold_curve$threshold, threshold_curve$recall1,
     type = "l", col = "#d53e4f", lwd = 2, ylim = c(0, 1),
     xlab = "Nguong phan loai", ylab = "Recall",
     main = "Recall hai lop theo nguong (out-of-fold tren tap train)")
lines(threshold_curve$threshold, threshold_curve$recall0,
      col = "#3288bd", lwd = 2)
abline(v = 0.5, col = "gray50", lty = 3)
abline(v = threshold, col = "#d53e4f", lty = 2)
grid(col = "gray90")
legend("right", bty = "n", cex = 0.85, lwd = 2,
       legend = c("Recall lop co benh",
                  "Recall lop khong benh",
                  "Nguong mac dinh 0.5",
                  sprintf("Nguong da chon %.3f", threshold)),
       col = c("#d53e4f", "#3288bd", "gray50", "#d53e4f"),
       lty = c(1, 1, 3, 2))


## ----threshold-lookup---------------------------------------------------------
print(data.frame(
  recall_muc_tieu = threshold_lookup$target_recall,
  nguong = round(threshold_lookup$threshold, 3),
  recall_lop1 = round(threshold_lookup$recall1, 3),
  recall_lop0 = round(threshold_lookup$recall0, 3)
), row.names = FALSE)

cat(sprintf("\nNguong mac dinh : 0.500\n"))
cat(sprintf("Nguong da chon  : %.3f (muc tieu recall >= %.2f)\n",
            threshold, summary_ref$target_recall))


## ----final-eval---------------------------------------------------------------
m_final <- print_confusion(y_test, pred_final,
                           sprintf("MA TRAN NHAM LAN (nguong %.3f)",
                                   threshold))
cat("\nCac chi so tren tap kiem tra:\n")
print(round(m_final, 3))


## ----final-tradeoff-----------------------------------------------------------
cmp <- as.data.frame(tradeoff)
rownames(cmp) <- cmp$threshold_label
cmp$threshold_label <- NULL
print(round(cmp, 3))


## ----final-boot---------------------------------------------------------------
print(round(ci_tbl, 4))


## ----forest-plot, echo=FALSE, fig.width=6.8, fig.height=4.5, fig.cap="Khoang tin cay bootstrap 95\\% cho sau chi so tren tap kiem tra"----
par(mar = c(5, 9, 4, 2))
n_m <- nrow(ci_tbl)
plot(ci_tbl$diem_uoc_luong, seq_len(n_m), pch = 19, col = "#3288bd",
     cex = 1.3, xlim = c(0, 1), yaxt = "n", ylab = "", xlab = "Gia tri",
     main = "Khoang tin cay bootstrap 95%")
arrows(ci_tbl$ci_duoi, seq_len(n_m), ci_tbl$ci_tren, seq_len(n_m),
       length = 0.05, angle = 90, code = 3, col = "gray40")
axis(2, seq_len(n_m), rownames(ci_tbl), las = 2, cex.axis = 0.85)
abline(v = 0.5, col = "#d53e4f", lty = 3)
par(mar = c(5, 4, 4, 2) + 0.1)

## ----old-r-pipeline-disabled--------------------------------------------------
if (FALSE) {
## ----boot-setup---------------------------------------------------------------
sp <- stratified_split(y_all, test_frac = 0.2, seed = 42)
X_train <- X_all[sp$train, , drop = FALSE]; y_train <- y_all[sp$train]
X_test  <- X_all[sp$test,  , drop = FALSE]; y_test  <- y_all[sp$test]

cat("Train:", nrow(X_train), "mẫu | Test:", nrow(X_test), "mẫu\n")

# Huấn luyện MỘT lần, dự đoán MỘT lần
fit_80 <- fit_pipeline(X_train, y_train, lambda = 1, smote_k = 5)
proba_test <- pipeline_proba(fit_80, X_test)
pred_test <- as.integer(proba_test >= 0.5)

recall_point <- compute_metrics(y_test, pred_test)["recall1"]


## ----boot-loop----------------------------------------------------------------
CONF <- 0.95
B <- 10000
n_test <- length(y_test)

set.seed(42)
recalls <- numeric(B)
n_valid <- 0
for (b in seq_len(B)) {
  # Rút có hoàn lại, cùng cỡ mẫu
  idx <- sample.int(n_test, n_test, replace = TRUE)
  if (sum(y_test[idx] == 1) == 0) next   # không có ca bệnh -> recall không xác định
  n_valid <- n_valid + 1
  recalls[n_valid] <- compute_metrics(y_test[idx], pred_test[idx])["recall1"]
}
recalls <- recalls[seq_len(n_valid)]

# Khoảng tin cậy percentile: cắt (1 - CONF) / 2 ở mỗi đuôi
tail_pct <- (1 - CONF) / 2
ci <- quantile(recalls, c(tail_pct, 1 - tail_pct))

cat(sprintf("Recall lớp có bệnh   = %.3f\n", recall_point))
cat(sprintf("Khoảng tin cậy %.0f%%   = [%.3f, %.3f]\n",
            CONF * 100, ci[1], ci[2]))
cat(sprintf("Bề rộng khoảng       = %.3f\n", ci[2] - ci[1]))
cat(sprintf("Sai số chuẩn         = %.3f\n", sd(recalls)))


## ----boot-hist, echo=FALSE, fig.width=6.5, fig.height=4.2, fig.cap="Phân phối bootstrap của recall lớp có bệnh"----
hist(recalls, breaks = 40, col = "#66c2a5", border = "white",
     main = sprintf("Phân phối bootstrap của recall (B = %d)",
                    length(recalls)),
     xlab = "Recall lớp có bệnh", ylab = "Tần số")
abline(v = recall_point, col = "#d53e4f", lwd = 2)
abline(v = ci, col = "gray40", lty = 2, lwd = 1.5)
legend("topright", bty = "n", cex = 0.85,
       legend = c(sprintf("Ước lượng điểm = %.3f", recall_point),
                  sprintf("KTC %.0f%%", CONF * 100)),
       col = c("#d53e4f", "gray40"), lty = c(1, 2), lwd = 2)


## ----kfold--------------------------------------------------------------------
K <- 5
folds <- stratified_folds(y_all, k = K, seed = 42)

cv_results <- data.frame()
for (i in seq_len(K)) {
  te <- which(folds == i)
  tr <- which(folds != i)

  fit_i <- fit_pipeline(X_all[tr, , drop = FALSE], y_all[tr],
                        lambda = 1, smote_k = 5)
  pred_i <- as.integer(pipeline_proba(fit_i,
                                      X_all[te, , drop = FALSE]) >= 0.5)

  m <- compute_metrics(y_all[te], pred_i)
  cv_results <- rbind(cv_results, data.frame(
    Fold = i, n_train = length(tr), n_test = length(te),
    Accuracy = m["accuracy"], Recall = m["recall1"],
    Macro_F1 = m["macro_f1"]
  ))
}
rownames(cv_results) <- NULL
print(round(cv_results, 4), row.names = FALSE)

for (m in c("Accuracy", "Recall", "Macro_F1")) {
  cat(sprintf("  %-9s = %.3f +/- %.3f\n",
              m, mean(cv_results[[m]]), sd(cv_results[[m]])))
}


## ----grid-search--------------------------------------------------------------
param_grid <- expand.grid(
  lambda  = c(0.1, 1, 10),   # cường độ điều chuẩn Ridge
  smote_k = c(3, 5, 7)                  # số láng giềng SMOTE
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
    pred_g <- as.integer(pipeline_proba(
      fit_g, X_train[te, , drop = FALSE]) >= 0.5)
    scores[i] <- compute_metrics(y_train[te], pred_g)["macro_f1"]
  }
  grid_mean[g] <- mean(scores)
  grid_sd[g] <- sd(scores)     # dao động giữa 5 fold, đo mức ổn định
}

param_grid$macro_f1 <- grid_mean
param_grid$sd_fold <- grid_sd

print(head(param_grid[order(-param_grid$macro_f1), ], 5), row.names = FALSE)


## ----grid-best----------------------------------------------------------------
i_top <- which.max(param_grid$macro_f1)
se_top <- param_grid$sd_fold[i_top] / sqrt(5)   # sai số chuẩn của trung bình

# Nhóm cấu hình tương đương: cách điểm cao nhất không quá một sai số chuẩn
tied <- which(param_grid$macro_f1 >= param_grid$macro_f1[i_top] - se_top)

# Trong nhóm đó, ưu tiên lambda nhỏ nhất
cand <- tied[param_grid$lambda[tied] == min(param_grid$lambda[tied])]
best_idx <- cand[which.max(param_grid$macro_f1[cand])]
best <- param_grid[best_idx, ]

cat(sprintf("\nSai số chuẩn của trung bình 5 fold: %.4f\n", se_top))
cat(sprintf("Số cấu hình tương đương: %d/%d\n",
            length(tied), nrow(param_grid)))
cat(sprintf("Cấu hình được chọn: lambda = %g, smote_k = %g (macro F1 = %.4f)\n",
            best$lambda, best$smote_k, best$macro_f1))
cat(sprintf("Đã thử %d cấu hình x 5 fold = %d lần khớp mô hình\n",
            nrow(param_grid), nrow(param_grid) * 5))


## ----threshold-oof------------------------------------------------------------
oof_proba <- numeric(length(y_train))
for (i in seq_len(5)) {
  te <- which(folds_inner == i); tr <- which(folds_inner != i)
  fit_o <- fit_pipeline(X_train[tr, , drop = FALSE], y_train[tr],
                        lambda = best$lambda, smote_k = best$smote_k)
  oof_proba[te] <- pipeline_proba(fit_o, X_train[te, , drop = FALSE])
}

grid_thr <- seq(0.05, 0.95, by = 0.005)
rec1 <- rec0 <- n_fn <- n_fp <- numeric(length(grid_thr))
for (i in seq_along(grid_thr)) {
  pred_t <- as.integer(oof_proba >= grid_thr[i])
  m <- compute_metrics(y_train, pred_t)
  rec1[i] <- m["recall1"]; rec0[i] <- m["recall0"]
  n_fn[i] <- sum(y_train == 1 & pred_t == 0)   # ca bệnh bị bỏ sót
  n_fp[i] <- sum(y_train == 0 & pred_t == 1)   # người khoẻ bị báo nhầm
}


## ----threshold-plot, fig.width=6.8, fig.height=4.5, fig.cap="Đánh đổi giữa độ nhạy và độ đặc hiệu khi thay đổi ngưỡng phân loại"----
plot(grid_thr, rec1, type = "l", col = "#d53e4f", lwd = 2, ylim = c(0, 1),
     xlab = "Ngưỡng phân loại", ylab = "Recall",
     main = "Recall hai lớp theo ngưỡng (out-of-fold trên tập train)")
lines(grid_thr, rec0, col = "#3288bd", lwd = 2)
abline(v = 0.5, col = "gray50", lty = 3)
grid(col = "gray90")
legend("right", bty = "n", cex = 0.85, lwd = 2,
       legend = c("Recall lớp có bệnh (độ nhạy)",
                  "Recall lớp không bệnh (độ đặc hiệu)",
                  "Ngưỡng mặc định 0.5"),
       col = c("#d53e4f", "#3288bd", "gray50"), lty = c(1, 1, 3))


## ----threshold-lookup---------------------------------------------------------
sel <- which(round(grid_thr, 3) %in% c(0.50, 0.45, 0.40, 0.375, 0.35,
                                       0.33, 0.30, 0.25))
print(data.frame(
  nguong = grid_thr[sel], FN = n_fn[sel], FP = n_fp[sel],
  tong_sai = n_fn[sel] + n_fp[sel],
  recall_lop1 = round(rec1[sel], 3), recall_lop0 = round(rec0[sel], 3)
), row.names = FALSE)

threshold <- 0.350
i_thr <- which(grid_thr == threshold)

cat(sprintf("\nNgưỡng mặc định : 0.500\n"))
cat(sprintf("Ngưỡng đã chọn  : %.3f\n", threshold))
cat(sprintf("  Trên out-of-fold: FN = %d, FP = %d, tổng = %d\n",
            n_fn[i_thr], n_fp[i_thr], n_fn[i_thr] + n_fp[i_thr]))


## ----final-eval---------------------------------------------------------------
fit_best <- fit_pipeline(X_train, y_train,
                         lambda = best$lambda, smote_k = best$smote_k)
proba_final <- pipeline_proba(fit_best, X_test)
pred_final <- as.integer(proba_final >= threshold)
pred_05 <- as.integer(proba_final >= 0.5)

m_final <- print_confusion(y_test, pred_final,
                           sprintf("MA TRAN NHAM LAN (nguong %.3f)",
                                   threshold))
cat("\nCác chỉ số trên tập kiểm tra:\n")
print(round(m_final, 3))


## ----final-tradeoff-----------------------------------------------------------
cmp <- rbind(
  `0.500 (mặc định)` = compute_metrics(y_test, pred_05),
  `đã chỉnh`         = compute_metrics(y_test, pred_final)
)
print(round(cmp, 3))


## ----final-boot---------------------------------------------------------------
CONF2 <- 0.95
B2 <- 10000
boot_mat <- matrix(NA, nrow = B2, ncol = length(m_final),
                   dimnames = list(NULL, names(m_final)))

for (b in seq_len(B2)) {
  idx <- sample.int(n_test, n_test, replace = TRUE)
  if (length(unique(y_test[idx])) < 2) next     # mẫu suy biến, bỏ qua
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
print(round(ci_tbl, 4))


## ----forest-plot, echo=FALSE, fig.width=6.8, fig.height=4.5, fig.cap="Khoảng tin cậy bootstrap 95\\% cho sáu chỉ số trên tập kiểm tra"----
par(mar = c(5, 9, 4, 2))
n_m <- nrow(ci_tbl)
plot(ci_tbl$diem_uoc_luong, seq_len(n_m), pch = 19, col = "#3288bd",
     cex = 1.3, xlim = c(0, 1), yaxt = "n", ylab = "", xlab = "Giá trị",
     main = sprintf("Khoảng tin cậy bootstrap %.0f%%", CONF2 * 100))
arrows(ci_tbl$ci_duoi, seq_len(n_m), ci_tbl$ci_tren, seq_len(n_m),
       length = 0.05, angle = 90, code = 3, col = "gray40")
axis(2, seq_len(n_m), rownames(ci_tbl), las = 2, cex.axis = 0.85)
abline(v = 0.5, col = "#d53e4f", lty = 3)
par(mar = c(5, 4, 4, 2) + 0.1)

## ----end-old-r-pipeline-disabled----------------------------------------------
}

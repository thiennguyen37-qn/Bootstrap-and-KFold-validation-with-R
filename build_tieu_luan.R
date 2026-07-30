# =============================================================================
# BUILD tieu_luan.Rnw -> tieu_luan.pdf
#
#   Rscript build_tieu_luan.R
#
# Vì sao cần script này thay vì gọi thẳng "knitr -> xelatex x2":
#
#   1. xelatex trên máy này thỉnh thoảng chết giữa chừng với access violation
#      (0xC0000005), khoảng 1/3 số lần chạy, không phụ thuộc tài liệu hay shell.
#      Lần chạy hỏng để lại .aux dở dang, khiến lần kế tiếp báo "undefined
#      references" — nên chỉ đếm số lượt SẠCH LIÊN TIẾP, và xoá .aux mỗi khi
#      gặp crash để bắt đầu lại từ trạng thái biết trước.
#
#   2. MiKTeX trên máy đã hỏng hoàn toàn (mọi engine thoát ngay với 0xC0E90002),
#      mà "xelatex" trên PATH vẫn trỏ về MiKTeX. Script gọi thẳng binary của
#      TinyTeX để không phụ thuộc thứ tự PATH.
# =============================================================================

RNW <- "tieu_luan.Rnw"
TEX <- sub("\\.Rnw$", ".tex", RNW)
BASE <- sub("\\.Rnw$", "", RNW)

# Cần 3 lượt xelatex sạch liên tiếp: lượt 1 dựng .aux/.toc, lượt 2 chốt mục lục,
# lượt 3 chốt tham chiếu chéo tới các hình.
NEED_STREAK <- 3
MAX_ATTEMPTS <- 40

# Tự dò trong TinyTeX thay vì ghi cứng đường dẫn: thư mục con của bin/ khác nhau
# theo hệ điều hành (windows, x86_64-darwin, ...). Không dùng "xelatex" trần vì
# trên máy này PATH vẫn trỏ về bản MiKTeX đã hỏng.
xelatex <- list.files(file.path(tinytex::tinytex_root(), "bin"),
                      pattern = "^xelatex(\\.exe)?$",
                      recursive = TRUE, full.names = TRUE)[1]
if (is.na(xelatex)) stop("Không tìm thấy xelatex trong TinyTeX. Cài bằng: ",
                         "tinytex::install_tinytex()")

aux_files <- paste0(BASE, c(".aux", ".toc", ".out"))
clean_aux <- function() unlink(aux_files)

# Các cảnh báo cho biết tài liệu chưa hội tụ, cần chạy thêm lượt nữa
WARN_PAT <- paste("Rerun to get", "Label\\(s\\) may have changed",
                  "undefined references", "Overfull \\\\hbox", sep = "|")

message("== knitr ==")
knitr::knit(RNW)

message("\n== xelatex ==")
clean_aux()
streak <- 0
attempt <- 0

while (streak < NEED_STREAK && attempt < MAX_ATTEMPTS) {
  attempt <- attempt + 1
  log_txt <- suppressWarnings(system2(
    xelatex,
    c("-interaction=nonstopmode", "-file-line-error", TEX),
    stdout = TRUE, stderr = TRUE
  ))
  crashed <- !is.null(attr(log_txt, "status")) && attr(log_txt, "status") != 0

  if (crashed) {
    # .aux có thể đã bị ghi dở -> bỏ đi để lượt sau không đọc phải rác
    streak <- 0
    clean_aux()
    message(sprintf("  lần %d: crash, xoá .aux và chạy lại", attempt))
  } else {
    streak <- streak + 1
    n_warn <- sum(grepl(WARN_PAT, log_txt))
    message(sprintf("  lần %d: OK (chuỗi %d/%d), %d cảnh báo",
                    attempt, streak, NEED_STREAK, n_warn))
    if (streak >= NEED_STREAK && n_warn == 0) break
    if (streak >= NEED_STREAK) streak <- NEED_STREAK - 1   # chưa hội tụ, chạy thêm
  }
}

pdf <- paste0(BASE, ".pdf")
if (!file.exists(pdf)) stop("Build thất bại sau ", attempt, " lượt.")

message(sprintf("\nXong sau %d lượt xelatex -> %s (%.0f KB)",
                attempt, pdf, file.size(pdf) / 1024))

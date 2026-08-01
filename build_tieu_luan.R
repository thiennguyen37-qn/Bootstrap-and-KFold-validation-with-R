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

# Dò engine xelatex theo thứ tự ưu tiên. TinyTeX trên máy này không dùng được:
# tlmgr của nó không chạy (thư mục dự án có dấu tiếng Việt khiến perl báo "unable
# to translate to a wide string"), nên không cài nổi các gói còn thiếu như
# extarticle.cls. MiKTeX thì compile bình thường, nên ưu tiên MiKTeX; vòng lặp
# retry bên dưới đã đủ để né các lần MiKTeX crash giữa chừng.
find_xelatex <- function() {
  miktex <- file.path(Sys.getenv("LOCALAPPDATA"),
                      "Programs", "MiKTeX", "miktex", "bin", "x64", "xelatex.exe")
  if (file.exists(miktex)) return(miktex)

  tt_root <- tryCatch(tinytex::tinytex_root(), error = function(e) "")
  if (nzchar(tt_root)) {
    tt <- list.files(file.path(tt_root, "bin"), pattern = "^xelatex(\\.exe)?$",
                     recursive = TRUE, full.names = TRUE)
    if (length(tt) && !is.na(tt[1])) return(tt[1])
  }

  p <- Sys.which("xelatex")
  if (nzchar(p)) return(unname(p))

  stop("Không tìm thấy xelatex (đã thử MiKTeX, TinyTeX, PATH).")
}
xelatex <- find_xelatex()
message("Dùng engine: ", xelatex)

aux_files <- paste0(BASE, c(".aux", ".toc", ".out"))
clean_aux <- function() unlink(aux_files)

# Các cảnh báo cho biết tài liệu CHƯA HỘI TỤ, tức chạy thêm lượt nữa thì hết.
# Chỉ ba loại này mới được dùng làm điều kiện dừng.
RERUN_PAT <- paste("Rerun to get", "Label\\(s\\) may have changed",
                   "undefined references", sep = "|")

# Overfull hbox là lỗi TRÌNH BÀY (dòng tràn lề), chạy lại bao nhiêu lượt cũng
# không mất. Trước đây nó nằm chung trong điều kiện dừng, nên chỉ cần tài liệu
# có một dòng tràn lề là vòng lặp chạy đủ MAX_ATTEMPTS lượt. Giờ chỉ đếm để
# báo lại, không chặn việc dừng.
OVERFULL_PAT <- "Overfull \\\\hbox"

message("== knitr ==")
knitr::knit(RNW)

message("\n== xelatex ==")
clean_aux()
streak <- 0
attempt <- 0
n_over <- 0        # số dòng tràn lề ở lượt chạy sạch gần nhất

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
    n_rerun <- sum(grepl(RERUN_PAT, log_txt))
    n_over <- sum(grepl(OVERFULL_PAT, log_txt))
    message(sprintf("  lần %d: OK (chuỗi %d/%d), %d cần chạy lại, %d tràn lề",
                    attempt, streak, NEED_STREAK, n_rerun, n_over))
    if (streak >= NEED_STREAK && n_rerun == 0) break
    if (streak >= NEED_STREAK) streak <- NEED_STREAK - 1   # chưa hội tụ, chạy thêm
  }
}

pdf <- paste0(BASE, ".pdf")
if (!file.exists(pdf)) stop("Build thất bại sau ", attempt, " lượt.")

message(sprintf("\nXong sau %d lượt xelatex -> %s (%.0f KB)",
                attempt, pdf, file.size(pdf) / 1024))

if (n_over > 0) {
  message(sprintf("Lưu ý: %d dòng tràn lề (Overfull hbox). Không ảnh hưởng",
                  n_over))
  message("việc build, nhưng nên xem lại. Tìm 'Overfull' trong ",
          BASE, ".log để biết dòng nào.")
}

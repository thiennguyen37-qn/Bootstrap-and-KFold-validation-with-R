$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
$reportDir = Join-Path $repoRoot "FINAL REPORT"

Set-Location -LiteralPath $reportDir

Rscript -e "knitr::knit('tieu_luan.Rnw', quiet = TRUE)"
xelatex -interaction=batchmode -halt-on-error -file-line-error tieu_luan.tex
xelatex -interaction=batchmode -halt-on-error -file-line-error tieu_luan.tex

Write-Host "Built FINAL REPORT/tieu_luan.pdf"

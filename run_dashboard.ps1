# Launch the IBD miRNA-seq Streamlit dashboard.
#   pwsh -File run_dashboard.ps1
# Finds a Python interpreter that has Streamlit installed and serves the app.

$ErrorActionPreference = "Stop"
$app = Join-Path $PSScriptRoot "pipelines/08_report/dashboard/app.py"

$candidates = @(
    "python",
    "$env:USERPROFILE\anaconda3\python.exe",
    "$env:USERPROFILE\miniconda3\python.exe",
    "py -3.11"
)

$py = $null
foreach ($c in $candidates) {
    try {
        & cmd /c "$c -c ""import streamlit""" 2>$null
        if ($LASTEXITCODE -eq 0) { $py = $c; break }
    } catch {}
}

if (-not $py) {
    Write-Error "No Python with Streamlit found. Install with: pip install -r requirements.txt"
    exit 1
}

Write-Host "Launching dashboard with: $py" -ForegroundColor Green
& cmd /c "$py -m streamlit run ""$app"""

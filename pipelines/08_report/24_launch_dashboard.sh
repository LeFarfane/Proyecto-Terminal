#!/usr/bin/env bash
# Stage 08 / step 21 — Launch the IBD miRNA-seq Streamlit dashboard (WSL/Linux).
#
#   bash 24_launch_dashboard.sh
#   (or: chmod +x 24_launch_dashboard.sh && ./24_launch_dashboard.sh)
#
# Finds a Python interpreter that has Streamlit installed and serves the app at
# dashboard/Dashboard.py. Tries the project conda envs first, then any python on PATH.
set -euo pipefail

# Resolve this script's directory so it works from any cwd.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP="$SCRIPT_DIR/dashboard/Dashboard.py"

# Locate conda (WSL miniconda is the project default).
CONDA_BIN=""
for c in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" "$(command -v conda || true)"; do
    if [[ -n "$c" && -x "$c" ]]; then CONDA_BIN="$c"; break; fi
done

has_streamlit() { "$@" -c "import streamlit" >/dev/null 2>&1; }

PY=()
# 1) Try the project conda envs via `conda run`.
if [[ -n "$CONDA_BIN" ]]; then
    for env in dash_env py_env r_env base; do
        if has_streamlit "$CONDA_BIN" run -n "$env" python; then
            PY=("$CONDA_BIN" run --no-capture-output -n "$env" python)
            break
        fi
    done
fi
# 2) Fall back to any python on PATH.
if [[ ${#PY[@]} -eq 0 ]]; then
    for p in python3 python; do
        if command -v "$p" >/dev/null 2>&1 && has_streamlit "$p"; then
            PY=("$p")
            break
        fi
    done
fi

if [[ ${#PY[@]} -eq 0 ]]; then
    echo "ERROR: No Python with Streamlit found." >&2
    echo "Install it into an env, e.g.:  pip install -r requirements.txt" >&2
    exit 1
fi

echo "Launching dashboard with: ${PY[*]}"
echo "Open http://localhost:8501 in your browser (Ctrl+C to stop)."
# --server.headless avoids the "gio: Operation not supported" browser-launch
# error under WSL, where there is no default browser handler.
exec "${PY[@]}" -m streamlit run "$APP" --server.headless true "$@"

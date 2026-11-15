#!/usr/bin/env bash
#
# Install bacteria-area into a local virtual environment.
# Works on Linux and macOS (no conda, no root, no system packages touched).
#
#   ./install.sh          # runtime dependencies + the package
#   ./install.sh --dev    # ... plus pytest, for running the test suite
#
# The environment is created in ./.venv; activate it with
#   source .venv/bin/activate
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="${HERE}/.venv"
PY_MIN_MAJOR=3
PY_MIN_MINOR=9

usage() {
    cat <<EOF
usage: ./install.sh [--dev] [--help]

  --dev    also install the test dependencies (pytest)
EOF
}

DEV=0
while [ $# -gt 0 ]; do
    case "$1" in
        --dev) DEV=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# ---------------------------------------------------------------- interpreter
find_python() {
    for c in python3.13 python3.12 python3.11 python3.10 python3.9 python3 python; do
        if command -v "$c" >/dev/null 2>&1; then
            if "$c" -c "import sys; sys.exit(0 if sys.version_info[:2] >= (${PY_MIN_MAJOR}, ${PY_MIN_MINOR}) else 1)" 2>/dev/null; then
                command -v "$c"
                return 0
            fi
        fi
    done
    return 1
}

if ! PY="$(find_python)"; then
    echo "ERROR: no Python >= ${PY_MIN_MAJOR}.${PY_MIN_MINOR} found on PATH." >&2
    echo "       install one (macOS: 'brew install python@3.12' or the python.org" >&2
    echo "       installer; Debian/Ubuntu: 'sudo apt install python3-venv') and re-run." >&2
    exit 1
fi

echo "==> interpreter: $PY ($("$PY" --version 2>&1))"

# ------------------------------------------------------------------- venv
if [ ! -x "${VENV}/bin/python" ]; then
    echo "==> creating virtual environment in ${VENV}"
    if ! "$PY" -m venv "$VENV"; then
        echo "ERROR: could not create the virtual environment." >&2
        echo "       on Debian/Ubuntu: sudo apt install python3-venv" >&2
        echo "       on macOS: install Python from python.org or 'brew install python'" >&2
        exit 1
    fi
else
    echo "==> reusing the virtual environment in ${VENV}"
fi

VPY="${VENV}/bin/python"

echo "==> upgrading pip / setuptools / wheel"
"$VPY" -m pip install --upgrade pip setuptools wheel

echo "==> installing bacteria-area (editable) and its dependencies"
"$VPY" -m pip install -e "$HERE"

if [ "$DEV" = "1" ]; then
    echo "==> installing the test dependencies"
    "$VPY" -m pip install pytest
fi

echo
echo "==> done."
cat <<EOF

Run the analysis (no activation needed):
    ${VENV}/bin/bacteria-area --input path/to/images --outdir results

Or activate the environment first:
    source ${VENV}/bin/activate
    bacteria-area --input path/to/images --outdir results

Try it on synthetic data:
    ${VPY} examples/make_synthetic.py
    ${VPY} -m bacteria_area --input examples/data --outdir results/example

Run the tests:
    ${VPY} -m pytest tests -q

Note: importing Leica .lif files needs Fiji/ImageJ with Bio-Formats
      (see lif_export.py); TIFF/PNG/JPEG inputs need no extra software.
EOF

#!/bin/bash
set -eo pipefail

# Change to script directory to ensure .env is found
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Load environment variables from .env file
if [[ -f .env ]]; then
  set -a
  # shellcheck disable=SC1091
  source <(grep -v '^#' .env | grep -v '^[[:space:]]*$')
  set +a
else
  echo "Error: .env file not found in $SCRIPT_DIR" >&2
  exit 1
fi

# Check if PFS_APP_HOSTNAME is set and matches current hostname
CURRENT_HOSTNAME=$(hostname -f)
if [[ -z "$PFS_APP_HOSTNAME" ]]; then
  echo "Error: PFS_APP_HOSTNAME is not set in .env file" >&2
  exit 1
fi

if [[ "$PFS_APP_HOSTNAME" != "$CURRENT_HOSTNAME" ]]; then
  echo "Hostname mismatch: PFS_APP_HOSTNAME=$PFS_APP_HOSTNAME, current hostname=$CURRENT_HOSTNAME"
  echo "Skipping application launch."
  exit 0
fi

echo "Hostname matches: $CURRENT_HOSTNAME"

# Server settings
ADDRESS="0.0.0.0"
ORIGIN="$PFS_APP_HOSTNAME"
PREFIX="fluxstd-lookup"

# Build command-line options based on PREFIX
if [[ -n "$PREFIX" ]]; then
  PREFIX_OPT="--prefix $PREFIX"
else
  PREFIX_OPT=""
fi

# Determine mode from first argument
MODE="${1:-production}"

if [[ "$MODE" == "dev" ]]; then
  echo "Running in development mode with auto-reload"
  PORT="5207"
  # shellcheck disable=SC2086
  uv run panel serve app.py \
    --address "$ADDRESS" \
    --port "$PORT" \
    --allow-websocket-origin="$ORIGIN:$PORT" \
    --allow-websocket-origin="localhost:$PORT" \
    --allow-websocket-origin="127.0.0.1:$PORT" \
    $PREFIX_OPT \
    --static-dirs "assets"="./assets" \
    --dev
else
  echo "Running in production mode with multi-threading"
  PORT="5107"
  # shellcheck disable=SC2086
  uv run panel serve app.py \
    --address "$ADDRESS" \
    --port "$PORT" \
    --allow-websocket-origin="$ORIGIN:$PORT" \
    --allow-websocket-origin="localhost:$PORT" \
    --allow-websocket-origin="127.0.0.1:$PORT" \
    $PREFIX_OPT \
    --static-dirs "assets"="./assets" \
    --num-threads 8
fi

echo "Application started successfully on $ORIGIN:$PORT"


# Notes
# python3 -m pip install --target "$LSST_PYTHON_USERLIB" panel watchfiles loguru ipywidgets_bokeh ipympl python-dotenv joblib datashader "holoviews[recommended]"


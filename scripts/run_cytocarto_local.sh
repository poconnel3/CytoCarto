#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV="$ROOT/.artifacts/cytocarto-venv"
API_PORT="${CYTOCARTO_API_PORT:-8000}"
WEB_PORT="${CYTOCARTO_WEB_PORT:-3000}"

if [[ ! -x "$VENV/bin/python" ]]; then
  python3 -m venv "$VENV"
fi
"$VENV/bin/pip" install --quiet -r "$ROOT/web_api/requirements.txt"

cleanup() {
  kill "${API_PID:-}" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

cd "$ROOT"
"$VENV/bin/python" -m uvicorn web_api.app:app --host 127.0.0.1 --port "$API_PORT" &
API_PID=$!

cd "$ROOT/web"
NEXT_PUBLIC_CYTOCARTO_API_URL="${NEXT_PUBLIC_CYTOCARTO_API_URL:-http://127.0.0.1:$API_PORT}" \
  npm run dev -- --port "$WEB_PORT"

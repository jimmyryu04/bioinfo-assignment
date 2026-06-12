#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

if command -v typst >/dev/null 2>&1; then
  typst compile final_report.typ final_report.pdf
else
  python build_report_pdf.py final_report.md final_report.pdf
fi

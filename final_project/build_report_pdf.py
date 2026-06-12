#!/usr/bin/env python3
"""Small fallback PDF builder for the final-project manuscript.

This is used only when Typst/Pandoc/LaTeX are unavailable on the local PATH.
It renders a readable report PDF directly from final_report.md and embeds the
first page of each figure PDF referenced by Markdown image links.
"""

from __future__ import annotations

import html
import os
import re
import sys
from pathlib import Path

import fitz


PAGE_W = 595.28
PAGE_H = 841.89
MARGIN_X = 54
MARGIN_TOP = 54
MARGIN_BOTTOM = 54
CONTENT_W = PAGE_W - (2 * MARGIN_X)

FONT_REG = os.environ.get("REPORT_FONT_REG", "tiro")
FONT_BOLD = os.environ.get("REPORT_FONT_BOLD", "tibo")
FONT_ITAL = os.environ.get("REPORT_FONT_ITAL", "tiit")
FONT_MONO = os.environ.get("REPORT_FONT_MONO", "cour")
FONT_FILE = os.environ.get("REPORT_FONT_FILE")
CUSTOM_FONT_NAMES = {FONT_REG, FONT_BOLD, FONT_ITAL, FONT_MONO} if FONT_FILE else set()
CUSTOM_FONT = fitz.Font(fontfile=FONT_FILE) if FONT_FILE else None

FIGURE_WIDTH_FRACTIONS = {
    "fig1_baseline_total_clip_vs_delta_rd.pdf": 0.82,
    "fig4_dose_response.pdf": 0.78,
}


REPLACEMENTS = {
    "\u2010": "-",
    "\u2011": "-",
    "\u2012": "-",
    "\u2013": "-",
    "\u2014": "-",
    "\u2018": "'",
    "\u2019": "'",
    "\u201c": '"',
    "\u201d": '"',
    "\u03b1": "alpha",
    "\u03b2": "beta",
    "\u03b4": "delta",
    "\u0394": "Delta",
    "\u03c1": "rho",
    "\u00d7": "x",
    "\u00b2": "2",
    "\u2032": "'",
    "\u2033": '"',
    "\u2264": "<=",
    "\u2265": ">=",
}


def clean_inline(text: str) -> str:
    text = html.unescape(text)
    text = re.sub(r"<sup>(.*?)</sup>", r"\1", text)
    text = re.sub(r"<[^>]+>", "", text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1 (\2)", text)
    text = text.replace("`", "")
    text = text.replace("**", "")
    text = text.replace("*", "")
    for old, new in REPLACEMENTS.items():
        text = text.replace(old, new)
    return re.sub(r"\s+", " ", text).strip()


def text_width(text: str, font: str, size: float) -> float:
    try:
        if CUSTOM_FONT is not None and font in CUSTOM_FONT_NAMES:
            return CUSTOM_FONT.text_length(text, fontsize=size)
        return fitz.get_text_length(text, fontname=font, fontsize=size)
    except Exception:
        return len(text) * size * 0.48


def split_long_word(word: str, max_width: float, font: str, size: float) -> list[str]:
    parts: list[str] = []
    current = ""
    for char in word:
        candidate = current + char
        if current and text_width(candidate, font, size) > max_width:
            parts.append(current)
            current = char
        else:
            current = candidate
    if current:
        parts.append(current)
    return parts


def wrap_text(text: str, max_width: float, font: str, size: float) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        if text_width(word, font, size) > max_width:
            fragments = split_long_word(word, max_width, font, size)
        else:
            fragments = [word]

        for fragment in fragments:
            candidate = fragment if not current else f"{current} {fragment}"
            if current and text_width(candidate, font, size) > max_width:
                lines.append(current)
                current = fragment
            else:
                current = candidate
    if current:
        lines.append(current)
    return lines


class ReportPdf:
    def __init__(self) -> None:
        self.doc = fitz.open()
        self.page: fitz.Page | None = None
        self.y = MARGIN_TOP
        self.new_page()

    def new_page(self) -> None:
        self.page = self.doc.new_page(width=PAGE_W, height=PAGE_H)
        if FONT_FILE:
            assert self.page is not None
            for font in CUSTOM_FONT_NAMES:
                self.page.insert_font(fontname=font, fontfile=FONT_FILE)
        self.y = MARGIN_TOP

    def ensure(self, height: float) -> None:
        if self.y + height > PAGE_H - MARGIN_BOTTOM:
            self.new_page()

    def draw_lines(
        self,
        lines: list[str],
        *,
        font: str = FONT_REG,
        size: float = 10.2,
        leading: float = 13.0,
        color: tuple[float, float, float] = (0, 0, 0),
        align: str = "left",
        space_before: float = 0,
        space_after: float = 5,
    ) -> None:
        if not lines:
            return
        height = space_before + (len(lines) * leading) + space_after
        self.ensure(height)
        self.y += space_before
        assert self.page is not None
        for line in lines:
            if align == "center":
                x = MARGIN_X + max(0, (CONTENT_W - text_width(line, font, size)) / 2)
            else:
                x = MARGIN_X
            self.page.insert_text(
                (x, self.y),
                line,
                fontname=font,
                fontsize=size,
                color=color,
            )
            self.y += leading
        self.y += space_after

    def paragraph(
        self,
        text: str,
        *,
        font: str = FONT_REG,
        size: float = 10.2,
        leading: float = 13.0,
        align: str = "left",
        space_before: float = 0,
        space_after: float = 5,
    ) -> None:
        text = clean_inline(text)
        if not text:
            return
        lines = wrap_text(text, CONTENT_W, font, size)
        self.draw_lines(
            lines,
            font=font,
            size=size,
            leading=leading,
            align=align,
            space_before=space_before,
            space_after=space_after,
        )

    def heading(self, text: str, level: int) -> None:
        text = clean_inline(text)
        if level == 1:
            lines = wrap_text(text, CONTENT_W, FONT_BOLD, 17.5)
            self.draw_lines(
                lines,
                font=FONT_BOLD,
                size=17.5,
                leading=21.5,
                align="center",
                space_after=14,
            )
        elif level == 2:
            self.ensure(42)
            self.y += 11
            lines = wrap_text(text, CONTENT_W, FONT_BOLD, 13.3)
            self.draw_lines(
                lines,
                font=FONT_BOLD,
                size=13.3,
                leading=16,
                space_after=6,
            )
        else:
            self.ensure(34)
            self.y += 7
            lines = wrap_text(text, CONTENT_W, FONT_BOLD, 11.2)
            self.draw_lines(
                lines,
                font=FONT_BOLD,
                size=11.2,
                leading=14,
                space_after=4,
            )

    def image_from_pdf(self, figure_path: Path) -> None:
        if not figure_path.exists():
            self.paragraph(f"[Missing figure: {figure_path}]", font=FONT_BOLD)
            return

        width_fraction = FIGURE_WIDTH_FRACTIONS.get(figure_path.name, 1.0)
        fig_doc = fitz.open(figure_path)
        try:
            fig_page = fig_doc[0]
            rect = fig_page.rect
            target_width = CONTENT_W * width_fraction
            scale = target_width / rect.width
            max_h = PAGE_H - MARGIN_TOP - MARGIN_BOTTOM
            if rect.height * scale > max_h * 0.72:
                scale = (max_h * 0.72) / rect.height
            width = rect.width * scale
            height = rect.height * scale
            self.ensure(height + 10)
            assert self.page is not None
            pix = fig_page.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
            x0 = MARGIN_X + (CONTENT_W - width) / 2
            target = fitz.Rect(x0, self.y, x0 + width, self.y + height)
            self.page.insert_image(target, stream=pix.tobytes("png"))
            self.y += height + 8
        finally:
            fig_doc.close()

    def finish(self, output_path: Path) -> None:
        for index, page in enumerate(self.doc, start=1):
            label = str(index)
            width = text_width(label, FONT_REG, 9)
            page.insert_text(
                ((PAGE_W - width) / 2, PAGE_H - 28),
                label,
                fontname=FONT_REG,
                fontsize=9,
                color=(0.25, 0.25, 0.25),
            )
        self.doc.save(output_path, deflate=True, garbage=4)
        self.doc.close()


def flush_paragraph(builder: ReportPdf, lines: list[str]) -> None:
    if not lines:
        return
    text = " ".join(line.strip() for line in lines)
    cleaned = clean_inline(text)
    if cleaned.startswith("Dongmin Ryu"):
        builder.paragraph(cleaned, size=10.5, leading=13, align="center", space_after=1)
    elif cleaned.startswith("Bioinformatics Final Project"):
        builder.paragraph(cleaned, size=10.5, leading=13, align="center", space_after=7)
    elif text.startswith("**") and text.endswith("**"):
        builder.paragraph(
            cleaned,
            font=FONT_BOLD,
            size=10.6,
            leading=13.4,
            space_before=2,
            space_after=8,
        )
    elif re.match(r"^\d+\.\s+", cleaned):
        builder.paragraph(cleaned, size=9.8, leading=12.4, space_after=4)
    else:
        builder.paragraph(cleaned)
    lines.clear()


def render_markdown(markdown_path: Path, output_path: Path) -> None:
    root = markdown_path.parent
    builder = ReportPdf()
    paragraph_lines: list[str] = []

    for raw in markdown_path.read_text(encoding="utf-8").splitlines():
        line = raw.rstrip()
        image_match = re.match(r"!\[[^\]]*\]\(([^)]+)\)", line.strip())
        heading_match = re.match(r"^(#{1,3})\s+(.*)$", line)

        if image_match:
            flush_paragraph(builder, paragraph_lines)
            builder.image_from_pdf(root / image_match.group(1))
            continue

        if heading_match:
            flush_paragraph(builder, paragraph_lines)
            level = len(heading_match.group(1))
            builder.heading(heading_match.group(2), level)
            continue

        if not line.strip():
            flush_paragraph(builder, paragraph_lines)
            continue

        paragraph_lines.append(line)

    flush_paragraph(builder, paragraph_lines)
    builder.finish(output_path)


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("usage: build_report_pdf.py final_report.md final_report.pdf", file=sys.stderr)
        return 2

    markdown_path = Path(argv[1]).resolve()
    output_path = Path(argv[2]).resolve()
    render_markdown(markdown_path, output_path)
    print(f"wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

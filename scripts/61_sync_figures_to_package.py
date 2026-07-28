#!/usr/bin/env python3
"""Sync the regenerated (uppercase-label, Arial) figures into the Scientific
Reports submission package.

Source of truth: results/figures_editable_v2_20260724/<name>.pdf  (the eight
editable cairo/Arial masters verified on 2026-07-24). For each figure this
rebuilds the three package representations from that single master:
  - figure_jpg/<name>.jpg          300-dpi white-background JPEG (upload asset)
  - figures_pdf_review/<name>.pdf  self-contained raster review PDF
  - figures_pdf_editable/<name>.pdf editable vector master (copied)
Raw data and analysis outputs are not touched.
"""
import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results" / "figures_editable_v2_20260724"
PKG = ROOT / "submission" / "Scientific_Reports_20260715"
PDFTOCAIRO = "/opt/homebrew/bin/pdftocairo"

NAMES = ["Fig1", "Fig2", "Fig3", "Fig4", "FigS1", "FigS2", "FigS3", "FigS4"]

(PKG / "figure_jpg").mkdir(parents=True, exist_ok=True)
(PKG / "figures_pdf_review").mkdir(parents=True, exist_ok=True)
(PKG / "figures_pdf_editable").mkdir(parents=True, exist_ok=True)


def jpg_from_pdf(pdf: Path, jpg: Path):
    with tempfile.TemporaryDirectory() as td:
        stem = Path(td) / "p"
        subprocess.run([PDFTOCAIRO, "-png", "-singlefile", "-r", "300",
                        str(pdf), str(stem)], check=True)
        im = Image.open(stem.with_suffix(".png")).convert("RGB")
        im.save(jpg, "JPEG", quality=95, subsampling=0, dpi=(300, 300))


def pdf_from_jpg(jpg: Path, pdf: Path):
    Image.open(jpg).convert("RGB").save(pdf, "PDF", resolution=300.0, quality=100, subsampling=0)


for name in NAMES:
    master = SRC / f"{name}.pdf"
    assert master.exists(), f"missing master {master}"
    jpg = PKG / "figure_jpg" / f"{name}.jpg"
    jpg_from_pdf(master, jpg)
    pdf_from_jpg(jpg, PKG / "figures_pdf_review" / f"{name}.pdf")
    shutil.copy2(master, PKG / "figures_pdf_editable" / f"{name}.pdf")
    print(f"{name}: jpg + review pdf + editable pdf synced")

print("Done. All 8 figures synced from", SRC.name)

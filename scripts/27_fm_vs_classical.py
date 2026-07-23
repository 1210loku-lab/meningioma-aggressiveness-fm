#!/usr/bin/env python3
"""Compatibility launcher for the corrected patient-level FM benchmark.

The historical implementation paired a tokenisation-reordered embedding with
cell metadata by row position. The canonical analysis is now scripts/52, which
uses the audited patient mapping from scripts/51 and repeated nested LOPO.
"""

from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
mapped = ROOT / "results/scrna/gf_out/tumor_emb_patient_mapped.csv"
if not mapped.exists():
    subprocess.run(["Rscript", "scripts/51_repair_geneformer_patient_mapping.R"], cwd=ROOT, check=True)
subprocess.run(["Rscript", "scripts/52_geneformer_patient_level_audit.R"], cwd=ROOT, check=True)

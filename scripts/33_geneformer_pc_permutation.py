#!/usr/bin/env python3
"""Compatibility launcher for the corrected patient-level permutation audit.

The historical script shuffled 1,800 cell labels even though WHO grade is a
patient-level label. Script 52 performs 100,000 patient-label permutations and
repeats best-of-top-10 component selection in every permutation.
"""

from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
mapped = ROOT / "results/scrna/gf_out/tumor_emb_patient_mapped.csv"
if not mapped.exists():
    subprocess.run(["Rscript", "scripts/51_repair_geneformer_patient_mapping.R"], cwd=ROOT, check=True)
subprocess.run(["Rscript", "scripts/52_geneformer_patient_level_audit.R"], cwd=ROOT, check=True)

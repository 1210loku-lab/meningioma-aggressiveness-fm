#!/usr/bin/env python3
"""Curate clinically exposed L1000 hits and integrate meningioma evidence."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
src = ROOT / "results/drug_repurposing/l1000_opposite_signatures.csv"
out = ROOT / "docs/Table_S4_translational_drug_candidates.csv"
rows = list(csv.DictReader(src.open()))

ann = {
 "Etoposide": ("marketed oncology drug", "meningioma cell-line evidence includes surviving drug-resistant stem-like cells", "deprioritise", "strongest current reversal but single MCF7 signature; nonspecific cytotoxicity/resistance concern", "https://pmc.ncbi.nlm.nih.gov/articles/PMC5521079/"),
 "Amsacrine": ("clinically used/marketed in some jurisdictions for acute leukaemia", "no direct meningioma validation identified", "deprioritise", "strongest current reversal but single MCF7 signature and substantial systemic toxicity", "https://maayanlab.cloud/dmoa/report/BRD-K98490050"),
 "itraconazole": ("marketed antifungal; oncology repurposing studied clinically", "no direct meningioma validation identified", "exploratory", "weak single-VCAP reversal; mechanistic and tumour-exposure validation required", "https://clinicaltrials.gov/study/NCT02357836"),
 "Perhexiline maleate": ("marketed for angina in limited jurisdictions", "no direct meningioma validation identified", "exploratory", "weak single-MCF7 reversal; metabolic mechanism is conceptually relevant to OXPHOS but hepatotoxicity and exposure are concerns", "https://pubchem.ncbi.nlm.nih.gov/compound/Perhexiline"),
 "Retinoic acid": ("marketed as tretinoin/all-trans retinoic acid", "in vitro meningioma evidence: reduced invasion/migration and increased matrix adhesion", "cross-validated candidate", "weak single-MCF7 reversal, but orthogonal meningioma-specific experimental support makes it testable", "https://doi.org/10.1038/sj.bjc.6690705"),
 "Danazol": ("marketed hormonal drug", "no therapeutic meningioma validation identified", "exclude", "hormonal pharmacology and meningioma context make repurposing unattractive despite reversal", ""),
 "ESTRADIOL BENZOATE": ("marketed hormonal drug in some jurisdictions", "no therapeutic meningioma validation identified", "exclude", "estrogenic exposure is not a rational unselected meningioma therapy", ""),
 "Norethindrone": ("marketed progestin", "no therapeutic meningioma validation identified", "exclude", "progestin exposure is a safety concern in meningioma", ""),
 "oxymetholone": ("marketed androgen/anabolic steroid", "no therapeutic meningioma validation identified", "exclude", "hormonal toxicity and absent meningioma evidence", ""),
 "AMIODARONE HYDROCHLORIDE": ("marketed antiarrhythmic", "no direct meningioma validation identified", "deprioritise", "weak single-HT29 reversal and major chronic toxicity", ""),
 "Cyclosporin A": ("marketed immunosuppressant", "no direct meningioma validation identified", "exclude", "immunosuppression is translationally unattractive for an unvalidated anticancer signal", ""),
 "mirtazapine": ("marketed antidepressant", "no direct meningioma validation identified", "deprioritise", "weakest listed single-VCAP reversal and no mechanistic support", ""),
}

selected=[]
for r in rows:
    drug=r.get("pert_desc","")
    if drug not in ann: continue
    status, evidence, decision, caveat, source = ann[drug]
    selected.append({"candidate":drug,"l1000_opposite_rank":int(r["opposite_rank"]),
      "l1000_score":r["scores"],"cell_line":r["cell_id"],"time_h":r["pert_time"],
      "dose":r["pert_dose"],"market_or_clinical_status":status,
      "meningioma_cross_validation":evidence,"decision":decision,"key_caveat":caveat,"source":source})
selected.sort(key=lambda x:x["l1000_opposite_rank"])
with out.open("w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=selected[0].keys()); w.writeheader(); w.writerows(selected)
print(f"wrote {len(selected)} clinically exposed hits to {out}")

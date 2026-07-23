#!/usr/bin/env python3
"""Query L1000FWD for compounds that oppose the 200-gene programme.

The raw response is retained for audit. Aggregation is deliberately descriptive:
L1000 signatures come from non-meningioma cell lines and are hypothesis-generating.
"""
from __future__ import annotations

import json
import csv
import hashlib
import time
from pathlib import Path
from statistics import median
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
PROGRAM = ROOT / "docs" / "Table_S2_aggressiveness_program_genes.csv"
OUTDIR = ROOT / "results" / "drug_repurposing"
OUTDIR.mkdir(parents=True, exist_ok=True)


def clean_symbols(rows: list[dict[str, str]], direction: str) -> list[str]:
    values = [row.get("symbol", "") for row in rows if row.get("direction") == direction]
    return [x for x in values if x and not x.startswith("ENSG") and "-AS" not in x and x != "nan"]


with PROGRAM.open(encoding="utf-8", newline="") as handle:
    program = list(csv.DictReader(handle))
up = clean_symbols(program, "up_in_WHO_II_vs_I")
down = clean_symbols(program, "down_in_WHO_II_vs_I")

payload = {"up_genes": up, "down_genes": down}
# The endpoint is case-sensitive. The upper-case legacy path currently returns a
# 307 redirect, which urllib intentionally does not replay for POST requests.
base = "https://maayanlab.cloud/l1000fwd"
query_file = OUTDIR / "l1000_query.json"
raw_file = OUTDIR / "l1000_topn_raw.json"
if query_file.exists() and raw_file.exists():
    saved_query = json.loads(query_file.read_text(encoding="utf-8"))
    if saved_query.get("up_genes") == up and saved_query.get("down_genes") == down:
        result_id = saved_query["result_id"]
        raw = json.loads(raw_file.read_text(encoding="utf-8"))
    else:
        raise RuntimeError("Saved L1000 query uses a different signature; move it before re-querying")
else:
    post = Request(
        f"{base}/sig_search",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "User-Agent": "Meningioma-AI/1.0"},
        method="POST",
    )
    with urlopen(post, timeout=120) as response:
        result_id = json.load(response)["result_id"]
    get = Request(f"{base}/result/topn/{result_id}", headers={"User-Agent": "Meningioma-AI/1.0"})
    with urlopen(get, timeout=120) as response:
        raw = json.load(response)

query_file.write_text(
    json.dumps({"result_id": result_id, "up_genes": up, "down_genes": down}, indent=2),
    encoding="utf-8",
)
raw_file.write_text(json.dumps(raw, indent=2), encoding="utf-8")

# API naming has varied across deployments; use the opposing/reversing collection.
hits = raw.get("opposite") or raw.get("reverse") or raw.get("reversers") or []
if not hits:
    raise RuntimeError(f"No opposing signatures found; response keys={list(raw)}")

# topn returns only score + signature ID. Resolve each signature so the exported
# table contains compound, cell line, dose and time rather than opaque BRD IDs.
resolved_hits = []
cache_file = OUTDIR / "l1000_signature_metadata_cache.json"
metadata_cache = json.loads(cache_file.read_text(encoding="utf-8")) if cache_file.exists() else {}
for rank, hit in enumerate(hits, start=1):
    record = dict(hit)
    record["opposite_rank"] = rank
    sig_id = record.get("sig_id")
    if sig_id:
        metadata = metadata_cache.get(str(sig_id))
        if metadata is None:
            last_error = None
            for attempt in range(5):
                try:
                    sig_request = Request(
                        f"{base}/sig/{quote(str(sig_id), safe='')}",
                        headers={"User-Agent": "Meningioma-AI/1.0"},
                    )
                    with urlopen(sig_request, timeout=45) as response:
                        metadata = json.load(response)
                    metadata_cache[str(sig_id)] = metadata
                    cache_file.write_text(json.dumps(metadata_cache, indent=2), encoding="utf-8")
                    break
                except (HTTPError, URLError, TimeoutError) as error:
                    last_error = error
                    if attempt < 4:
                        time.sleep(min(2 ** attempt, 8))
            if metadata is None:
                record["resolution_error"] = repr(last_error)
        if metadata:
            for key, value in metadata.items():
                record.setdefault(key, value)
    resolved_hits.append(record)
hits = resolved_hits

def flatten(value: dict, prefix: str = "") -> dict:
    out = {}
    for key, item in value.items():
        name = f"{prefix}.{key}" if prefix else key
        if isinstance(item, dict):
            out.update(flatten(item, name))
        elif isinstance(item, list):
            out[name] = ";".join(map(str, item))
        else:
            out[name] = item
    return out


flat = [flatten(hit) for hit in hits]
# L1000FWD uses -666 as a missing-value sentinel. Never aggregate all such rows
# as one compound. Resolve the one clinically relevant high-ranking BRD ID found
# in this query; retain other unresolved IDs explicitly for audit/exclusion.
known_brd_names = {"BRD-K98490050": "Amsacrine"}
for record in flat:
    if str(record.get("pert_desc", "")) == "-666":
        pert_id = str(record.get("pert_id", ""))
        record["pert_desc"] = known_brd_names.get(pert_id, f"UNRESOLVED[{pert_id}]")
cols = sorted({key for row in flat for key in row})
with (OUTDIR / "l1000_opposite_signatures.csv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=cols, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(flat)


def first_existing(columns: list[str], candidates: list[str]) -> str | None:
    for candidate in candidates:
        if candidate in columns:
            return candidate
    return None


drug_col = first_existing(cols, ["pert_desc", "pert_name", "drug_name", "name", "Name", "meta.pert_name", "sig_id"])
score_col = first_existing(cols, ["scores", "score", "similarity", "combined_score", "zscore", "scores.z-down"])
cell_col = first_existing(cols, ["cell_id", "cell_line", "meta.cell_id", "meta.cell_line"])
time_col = first_existing(cols, ["pert_time", "time", "meta.pert_time"])
phase_col = first_existing(cols, ["clinical_phase", "phase", "meta.clinical_phase"])

if drug_col is None:
    raise RuntimeError(f"Cannot identify drug-name field; columns={cols}")

rows = []
groups: dict[str, list[dict]] = {}
for record in flat:
    groups.setdefault(str(record.get(drug_col, "")), []).append(record)
for drug, group in groups.items():
    def unique_values(column: str | None) -> list[str]:
        if not column:
            return []
        return sorted({str(record[column]) for record in group if record.get(column) not in (None, "")})

    row = {
        "drug": drug,
        "n_opposing_signatures": len(group),
        "n_cell_lines": len(unique_values(cell_col)) if cell_col else "",
        "cell_lines": ";".join(unique_values(cell_col)),
        "times": ";".join(unique_values(time_col)),
        "clinical_phase_api": ";".join(unique_values(phase_col)),
    }
    if score_col:
        scores = []
        for record in group:
            try:
                scores.append(float(record.get(score_col)))
            except (TypeError, ValueError):
                pass
        row["median_api_score"] = median(scores) if scores else ""
        row["best_api_score"] = min(scores) if scores else ""
    rows.append(row)

summary = sorted(
    rows,
    key=lambda row: (row["n_opposing_signatures"], row["n_cell_lines"] if row["n_cell_lines"] != "" else -1),
    reverse=True,
)
summary_cols = list(summary[0]) if summary else ["drug"]
with (OUTDIR / "l1000_drug_reversal_summary.csv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_cols, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(summary)

print(f"result_id={result_id}")
print(f"query genes: up={len(up)} down={len(down)}")
print(f"opposing signatures={len(flat)} drugs={len(groups)}")
for row in summary[:25]:
    print(row)

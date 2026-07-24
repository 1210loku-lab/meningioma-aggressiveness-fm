#!/usr/bin/env python3
"""Refresh manuscript text after a Results-only structural revision.

This deliberately does not call the full package builder: the author-supplied
JPEG and editable-PDF figure assets in the submission package are read-only.
"""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "submission" / "Scientific_Reports_20260715"
BUILDER = ROOT / "scripts" / "44_build_scirep_submission_package.py"


def load_builder():
    spec = spec_from_file_location("scirep_builder", BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {BUILDER}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    builder = load_builder()
    markdown = builder.final_markdown()
    cited, entries, order = builder.validate_citations(markdown)
    if cited != entries or order != list(range(1, 31)):
        raise RuntimeError("Citation order or reference coverage changed unexpectedly")

    md_path = PACKAGE / "manuscript" / "Meningioma_Scientific_Reports_submission.md"
    docx_path = PACKAGE / "manuscript" / "Meningioma_Scientific_Reports_submission.docx"
    md_path.write_text(markdown, encoding="utf-8")
    builder.render_doc(md_path, docx_path)
    print(md_path)
    print(docx_path)


if __name__ == "__main__":
    main()

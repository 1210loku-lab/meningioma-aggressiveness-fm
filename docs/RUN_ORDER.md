# Reproducible Run Order

This file fixes the execution order for the submission manuscript. Raw GEO files under `data/raw/` are read-only inputs. Run commands from the repository root.

## Environment Capture

Before submission, capture the R and Python environments:

```bash
Rscript -e 'sink("docs/sessionInfo_R.txt"); sessionInfo(); sink()'
python -m pip freeze > docs/pip_freeze_python.txt
```

The Geneformer-specific Python environment used for the completed analysis is already frozen in `docs/geneformer_requirements.txt`.

## Bulk Discovery and Validation

1. `scripts/04_GSE136661_counts_labels.R`
   - Input: `data/raw/GSE136661_ex/*.txt.gz`, `data/raw/GSE136661_series_matrix.txt.gz`
   - Output: `data/raw/GSE136661_assembled.rds`
   - Purpose: assemble GSE136661 HTSeq counts and phenotype labels.

2. `scripts/05_aggressiveness_program.R`
   - Input: `data/raw/GSE136661_assembled.rds`
   - Outputs: `results/deg/GSE136661_aggressiveness_program.rds`, `results/deg/GSE136661_aggressiveness_WHO_II_vs_I.csv`
   - Used in: Figure 1A, Table S2, program definition.

3. `scripts/06_crosscohort_validation.R`
   - Inputs: `results/deg/GSE136661_aggressiveness_program.rds`, `data/raw/GSE16581_series_matrix.txt.gz`, `data/raw/GPL570.soft.gz`
   - Output: `results/deg/GSE16581_program_validation.rds`
   - Used in: Figure 1A, 1E, 1F; recurrence, grade, classifier and survival validation in GSE16581.

4. `scripts/08_gsea.R`
   - Input: `results/deg/GSE136661_aggressiveness_WHO_II_vs_I.csv`
   - Outputs: `results/deg/GSEA_GOBP_aggressiveness.csv`, `results/deg/GO_BP_up_program.csv`, `results/deg/GO_BP_dn_program.csv`
   - Used in: Figure 1D and pathway interpretation.

5. `scripts/12_classical_baseline.R`
   - Inputs: `results/deg/GSE136661_aggressiveness_program.rds`, `results/deg/GSE16581_program_validation.rds`
   - Output: `results/deg/classical_baseline.rds`
   - Used in: Figure 1E; LASSO internal and external AUC.

5a. `scripts/34_nested_lasso_validation.R`
   - Input: `data/raw/GSE136661_assembled.rds`
   - Outputs: `results/deg/R5_nested_lasso_validation.txt`, `results/deg/R5_nested_lasso_predictions.csv`
   - Purpose: independent internal validation using repeated five-fold outer CV, inner-CV lambda selection and training-fold-only z-scaling.
   - Used in: the primary internal LASSO AUC reported in the Abstract, Results and Figure 1 legend; supersedes the original `cv.glmnet` curve estimate for internal discrimination.

6. `scripts/13_bulk_figures.R`
   - Inputs: outputs from scripts 05, 06, 08, 12, and 16 plus local GEO phenotype files
   - Outputs: `results/figures_pub/Figure_bulk_program.pdf`, `results/figures_pub/Figure_bulk_program.png`, `results/figures_pub/Figure1_complete.pdf`, `results/figures_pub/Figure1_complete.png`
   - Used in: integrated Figure 1.

7. `scripts/16_GSE74385_recurrence.R` and `scripts/16b_fig_GSE74385.R`
   - Inputs: `data/raw/GSE74385/GSE74385_normalized.txt.gz`, `data/raw/GSE74385/GSE74385_series_matrix.txt.gz`, `results/deg/GSE136661_aggressiveness_program.rds`
   - Outputs: `results/deg/GSE74385_program_recurrence.csv`, `results/deg/GSE74385_program_recurrence.rds`, `results/figures_pub/fig_GSE74385_recurrence.pdf`
   - Used in: Figure 1B-C source data; independent recurrence validation and batch-confounded grade-I check.

8. `scripts/24_P0_robustness.R`
   - Inputs: `results/deg/GSE16581_program_validation.rds`, `results/deg/GSE74385_program_recurrence.rds`
   - Output: `results/deg/P0_robustness_results.txt`
   - Used in: Firth recurrence models, batch adjustment, and limitations.

9. `scripts/25_signature_benchmark.R`
   - Inputs: GSE16581 and GSE74385 validation objects plus program genes
   - Output: `results/deg/P1_signature_benchmark.txt`
   - Used in: proliferation/OXPHOS/random-signature benchmark.

10. `scripts/30_withingrade_proliferation_LRT.R`
    - Inputs: GSE16581 and GSE74385 validation objects plus proliferation scores
    - Output: `results/deg/R3A_withingrade_proliferation_LRT.txt`
    - Used in: within-grade recurrence checks and proliferation incremental LRT.
    - Important: GSE16581 Affymetrix intensities are log2-transformed before probe collapse and gene-wise z-scaling; the submission-stage audit identified and corrected the older untransformed implementation.

11. `scripts/28_nomogram_calibration.R`
    - Input: GSE16581 validation object
    - Outputs: `results/deg/P2_nomogram_calibration.txt`, `results/figures_pub/fig_nomogram_calibration.pdf`
    - Historical/retired: not used in the revised manuscript or submission package because 10 recurrence events do not support a calibrated clinical model.

## Single-Cell Analyses

12. `scripts/10_scrna_process.R`
    - Input: `data/raw/10x_full/`
    - Output: `results/scrna/GSE183655_processed.rds`
    - Purpose: process GSE183655 atlas.

13. `scripts/11b_reannotate.R`
    - Input: `results/scrna/GSE183655_processed.rds`
    - Outputs: `results/scrna/GSE183655_annotated.rds`, `results/scrna/aggr_score_by_celltype_v2.csv`, `results/scrna/fig_umap_celltype_v2.pdf`
    - Used in: Figure 2.

14. `scripts/15_GSE206647_grade.R` and `scripts/15b_fig_from_rds.R`
    - Input: `data/raw/GSE206647_ex/`
    - Outputs: `results/scrna/GSE206647_processed.rds`, `results/scrna/GSE206647_tumor_aggr_by_grade_sample.csv`, `results/scrna/fig_GSE206647_grade_program.pdf`
    - Used in: Figure 3 and Table S3.

15. `scripts/26_purity_control.R`
    - Input: `results/scrna/GSE206647_processed.rds`
    - Output: `results/deg/P1_purity_control.txt`
    - Used in: depth/library-size sensitivity checks.

16. `scripts/31_pseudobulk_cellcycle.R`
    - Input: `results/scrna/GSE206647_processed.rds`
    - Outputs: `results/deg/R3B_pseudobulk_cellcycle.txt`, `results/deg/R3B_pseudobulk_score_per_patient.csv`
    - Used in: primary per-patient pseudobulk grade association and mixed model control.

## Foundation-Model and Attribution Analyses

17. `scripts/17c_download_geneformer_assets.py`
    - Output: `results/scrna/gf_model_V2-104M/`
    - Purpose: download Geneformer V2-104M model assets.

18. `scripts/17a_export_tumor_for_geneformer.R`
    - Input: `results/scrna/GSE206647_processed.rds`
    - Outputs: `results/scrna/gf_export/cells.csv`, `results/scrna/gf_export/genes.csv`, `results/scrna/gf_export/counts.mtx`
    - Purpose: export grade-balanced tumour/tumour-like cells for Geneformer.

19. `scripts/17b_geneformer_embed.py`
    - Inputs: `results/scrna/gf_export/`, `results/scrna/gf_model_V2-104M/`
    - Outputs: `results/scrna/gf_out/geneformer_summary.txt`, `results/scrna/gf_out/tumor_emb.csv`
    - Used in: Figure 5 zero-shot embedding source data.

20. `scripts/33_geneformer_pc_permutation.py`
    - Compatibility launcher for `scripts/51` and `scripts/52`.
    - The retired implementation shuffled cell labels and must not be used.

21. `scripts/27_fm_vs_classical.py`
    - Compatibility launcher for the corrected `scripts/51` + `scripts/52` patient-level workflow.
    - Historical outputs are retired because embedding rows were not aligned to cell metadata.

22. `scripts/29_anchored_perturbation.py`
    - Inputs: Geneformer tumour-state dataset and selected expressed program genes
    - Outputs: `results/deg/P2_anchored_perturbation.txt`, `results/deg/P2_anchored_shift_stats.csv`, `results/scrna/gf_anchor/`
    - Used in: transparent perturbation result.

23. `scripts/17c_geneformer_composite.R`
    - Inputs: `results/scrna/gf_out/tumor_emb.csv`, `results/deg/R4_geneformer_pc_permutation.txt`, `results/deg/P2_fm_vs_classical.txt`, `results/deg/P2_anchored_shift_stats.csv`
    - Outputs: `results/figures_pub/fig_geneformer_embedding.pdf`, `results/figures_pub/fig_geneformer_embedding.png`
    - Historical Figure 5 only; superseded by `scripts/56_rebuild_review_response_figures.R`.

24. `scripts/18_classical_drivers.R`
    - Inputs: program genes and GSE136661 discovery-cohort expression
    - Outputs: `results/deg/classical_drivers_ranked.csv`, `results/figures_pub/fig_classical_drivers.pdf`
    - Used in: Figure 4 candidate gene ranking.

## Translational Candidate Prioritisation

25. `scripts/36_target_evidence_matrix.R`
    - Inputs: `docs/Table_S2_aggressiveness_program_genes.csv`, `results/scrna/GSE206647_processed.rds`, and saved GSE16581/GSE74385 validation gene sets.
    - Outputs: `results/drug_repurposing/target_evidence_matrix.csv`, `results/drug_repurposing/target_evidence_top40.csv`.
    - Purpose: patient-level scRNA grade association, tumour specificity, cross-cohort recovery, and transparent internal target/biomarker evidence scoring.
    - Important: pseudobulk CPM denominators use the full transcriptome library size. `scripts/42_finalize_corrected_target_matrix.R` applies the independently recomputed statistics when the offloaded historical Seurat object is unavailable.

## Submission-stage independent audit and build

30. `scripts/39_fm_classical_leakage_audit.R`
    - Outputs: `results/audit_submission/fm_classical_leakage_audit_*.csv` and `.txt`.
    - Purpose: 10 repeats of five-fold patient-grouped outer CV with fold-contained preprocessing and patient aggregation.

31. `scripts/40_independent_headline_audit.R`
    - Reads assembled/raw bulk inputs, raw GSE206647 10x matrices and saved embeddings; it does not read historical headline result tables.
    - Outputs: `results/audit_submission/independent_headline_audit.*` and independently recomputed score/evidence tables.
    - **Partially superseded (2026-07-17):** its `gse74385_firth_*` (old composite endpoint) and `geneformer_*` PC/permutation fields (old cell-level) are pre-revision; current values come from `scripts/53` (recurrence-only, OR 5.18) and `scripts/51`+`52` (patient-level Geneformer, p=0.0516/0.00070). See `results/audit_submission/independent_headline_audit_SUPERSEDED_NOTE.md`. GSE16581/GSE206647/target-gene fields remain current.

32. `scripts/41_build_citation_audit.py`
    - Outputs: citation verification table, semantic-support notes and Zotero manual insertion map in `results/audit_submission/`.

33. `scripts/42_finalize_corrected_target_matrix.R`
    - Replaces the old program-only-library CPM statistics with the independently verified full-library-normalised statistics and refreshes Table S5.

34. `scripts/45_export_submission_source_data.R`
    - Exports compact source-data CSV files for the submitted figures.

35. `scripts/44_build_scirep_submission_package.py`
    - Builds the target-specific manuscript DOCX/Markdown, supplementary files, figure JPG/PDF copies, source data, citation maps and checksummed package manifest.

36. `scripts/51_repair_geneformer_patient_mapping.R`
    - Inputs: `results/scrna/gf_out/tumor_emb.csv`, `results/scrna/gf_export/cells.csv`.
    - Outputs: `results/scrna/gf_out/tumor_emb_patient_mapped.csv`, `results/audit_submission/geneformer_patient_mapping_audit.csv`, `.txt`.
    - Purpose: recover and audit the one-to-one patient/cell mapping for the historical embedding without assuming row order. Future embeddings retain `gsm` and `cell_id` directly via `scripts/17b_geneformer_embed.py`.

37. `scripts/52_geneformer_patient_level_audit.R`
    - Input: `results/scrna/gf_out/tumor_emb_patient_mapped.csv` and matched program-gene expression.
    - Outputs: `results/audit_submission/geneformer_patient_level_*`.
    - Purpose: patient-first and aggregation-order sensitivity PCA, 100,000 patient-label permutations, and repeated nested leave-one-patient-out Geneformer/classical benchmarks. This supersedes the inferential outputs of historical scripts 27, 33 and 39.

38. `scripts/53_GSE74385_endpoint_sensitivity.R`
    - Outputs: `results/audit_submission/GSE74385_endpoint_sensitivity.csv`, `.txt`.
    - Purpose: separate recurrence from malignant progression and fit score + grade + batch Firth models.

39. `scripts/54_signature_scoring_sensitivity.R`
    - Outputs: `results/audit_submission/signature_scoring_sensitivity.csv`, `signature_gene_set_sensitivity_summary.csv`, `.txt`.
    - Purpose: common-platform genes, within-sample rank scores, program-size and direction-balance sensitivity analyses.

40. `scripts/55_copykat_cnv_sensitivity.R`
    - Inputs: raw GSE206647 10x matrices and `results/scrna/GSE206647_cellmeta.csv`.
    - Outputs: `results/audit_submission/copykat_GSE206647/` compact patient predictions; `results/audit_submission/copykat_all_patient_summary.csv`, `copykat_cnv_confirmed_pseudobulk_scores.csv`, and `copykat_cnv_sensitivity_summary.txt`; `docs/Table_S6_CopyKAT_CNV_sensitivity.csv`.
    - Purpose: patient-wise CopyKAT sensitivity analysis of tumour/tumour-like annotations using immune cells as same-patient diploid references.

41. `scripts/56_rebuild_review_response_figures.R`
    - Inputs: revised Geneformer, endpoint and scoring audit tables.
    - Outputs: `results/figures_pub/fig_geneformer_patient_level.*`, `fig_review_sensitivity.*`.
    - Used in: revised Figure 5 and Figure S2.

26. `scripts/35_l1000_drug_reversal.py`
    - Input: `docs/Table_S2_aggressiveness_program_genes.csv`.
    - Planned outputs: `results/drug_repurposing/l1000_query.json`, `l1000_topn_raw.json`, `l1000_opposite_signatures.csv`, `l1000_drug_reversal_summary.csv`.
    - Purpose: query official L1000FWD for compounds opposing the fixed program. Requires network access; absence of the planned outputs means that present-signature reversal has not been completed.

27. `scripts/37_integrate_translational_candidates.py`
    - Input: `results/drug_repurposing/l1000_opposite_signatures.csv`.
    - Output: `docs/Table_S4_translational_drug_candidates.csv`.
    - Purpose: restrict present-query hits to marketed/clinically exposed compounds and adjudicate them against independent meningioma evidence and safety/biological plausibility.

## Supplementary Descriptive Analyses

28. `scripts/19_bulk_deconvolution.R`
    - Input: GSE136661 VST matrix
    - Output: `results/descriptive/19_GSE136661_mcpcounter_scores.csv`
    - Used in: cautious compositional interpretation.

29. `scripts/23_cellchat_analysis.R`
    - Input: `results/scrna/GSE183655_annotated.rds`
    - Outputs: `results/descriptive/23_cellchat_results.rds`, `results/descriptive/23_fig_cellchat_network.pdf`, `results/descriptive/23_fig_cellchat_pathways.pdf`
    - Used in: descriptive supplementary cell-cell communication analysis.

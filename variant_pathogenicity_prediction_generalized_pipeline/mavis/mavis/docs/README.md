# Documentation

## Files

- **`technical_methods.docx`** — Comprehensive pipeline technical document suitable for inclusion in a Materials & Methods section. Covers all 12 sections of the pipeline architecture including the three-axis ΔΔG framework, 16-category mechanism classification, four-way concordance, and output column definitions.

## Pipeline Overview

MAVIS evaluates missense variants across four stages:

### Stage 1: Structural Feature Extraction
- Per-residue pLDDT confidence from AlphaFold B-factors
- Intra-chain contact counts (5.0 Å cutoff, sequence separation ≥ 3)
- Inter-chain contacts at protein–protein interfaces
- Relative solvent accessibility (Shrake–Rupley)
- Secondary structure (DSSP)
- Grantham physicochemical distance
- Neighborhood ±3 weighted contact scores

### Stage 2: FoldX Thermodynamic Calculations
- **Axis 1 (ddg_monomer):** Intrinsic fold stability via BuildModel on monomer
- **Axis 2 (ddg_fold_{partner}):** Fold stability in complex context via BuildModel on complex
- **Axis 3 (ddg_binding_{partner}):** Interaction energy change via AnalyseComplex (ΔIE)
- All structures pre-processed with RepairPDB; results cached per structure
- pLDDT confidence gating: monomer DDG gated by monomer_plddt ≥ 50; per-partner DDG gated by partner pLDDT ≥ 50

### Stage 3: Scoring and Classification
- Composite disruption score from contact disruption + interface bonus + burial bonus
- pLDDT multiplier (×0.4 / ×0.7 / ×1.0) on final score
- Tier assignment: T1 ≥ 5.0, T2 ≥ 3.0, T3 ≥ 1.5, T4 < 1.5
- 16-category mechanism classification based on fold/PPI/interface/burial

### Stage 4: Concordance (Full pipeline only)
- Four-way concordance: structure × ΔΔG × AlphaMissense × Franklin/ClinVar
- Adjusted X/N denominators (only evaluable axes contribute)
- Strict, relaxed, and T3-inclusive flavors
- Dual-pipeline agreement tracking (single-residue vs. neighborhood ±3)

## ΔΔG Classification Thresholds

| Range (kcal/mol) | Category |
|---|---|
| > 2.0 | Highly destabilizing |
| 1.0 – 2.0 | Destabilizing |
| 0.5 – 1.0 | Mildly destabilizing |
| −0.5 – 0.5 | Neutral |
| −1.0 – −0.5 | Mildly stabilizing |
| −2.0 – −1.0 | Stabilizing |
| < −2.0 | Highly stabilizing |

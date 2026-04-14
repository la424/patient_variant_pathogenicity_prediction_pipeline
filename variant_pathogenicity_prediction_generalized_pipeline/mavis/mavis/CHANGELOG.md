# Changelog

## v1.0.0 (April 2026)

Initial public release based on CHD Pipeline v6.0 + Addendum A3-6.

### Features
- Three-axis ΔΔG framework (monomer fold, complex fold, binding energy)
- 16-category pathogenic mechanism classification
- Dual scoring pipelines (single-residue + neighborhood ±3)
- Four-way concordance (structure × ΔΔG × AlphaMissense × Franklin/ClinVar)
- Integrated FoldX runner with RepairPDB pre-processing
- pLDDT confidence gating at both data inclusion and scoring layers
- Auto-detection of monomer/multimer structures from directory
- Support for construct offset corrections (truncated multimer constructs)
- Bidirectional chain assignment for variants in genes appearing as chain B
- Idempotent re-run protection for FoldX cells
- DRY_RUN mode for testing pipeline logic without FoldX

### Pipelines
- `mavis_structural.ipynb` — Structure + FoldX only (no external annotations)
- `mavis_full.ipynb` — Full four-way concordance with AlphaMissense and Franklin/ClinVar

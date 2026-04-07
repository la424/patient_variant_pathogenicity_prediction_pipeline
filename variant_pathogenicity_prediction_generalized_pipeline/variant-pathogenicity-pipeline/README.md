# Variant Pathogenicity Pipeline (VarPat)

A config-driven computational pipeline for assessing missense variant pathogenicity using AlphaFold-predicted structures, FoldX thermodynamics, AlphaMissense predictions, and clinical annotations.

## Overview

VarPat integrates four independent lines of evidence into a unified scoring, tiering, and mechanism-classification framework for missense variants in protein-coding genes:

1. **Structural analysis** — pLDDT confidence, intra-chain contacts, inter-chain (interface) contacts, solvent accessibility, and burial classification from AlphaFold monomer and multimer predictions
2. **Thermodynamic stability** — Three-axis FoldX ΔΔG framework: monomer fold stability, partner-context fold stability, and binding energy change
3. **Machine learning** — AlphaMissense pathogenicity predictions
4. **Clinical annotation** — Franklin/ClinVar classifications

The pipeline produces a comprehensive variant-level table with disruption scores, tier assignments (1–4), 16-category mechanism classifications, and four-way concordance scores.

## Installation

```bash
git clone https://github.com/your-username/variant-pathogenicity-pipeline.git
cd variant-pathogenicity-pipeline
pip install -e .
```

### Dependencies

| Software | Version | Required | Purpose |
|----------|---------|----------|---------|
| Python | ≥3.10 | Yes | Runtime |
| Biopython | ≥1.82 | Yes | PDB/CIF parsing, SASA |
| pandas | ≥2.0 | Yes | Data manipulation |
| numpy | ≥1.24 | Yes | Numerical operations |
| PyYAML | ≥6.0 | Yes | Config parsing |
| openpyxl | ≥3.1 | Yes | XLSX output |
| FoldX | 5.x | Optional | ΔΔG calculations |
| mkdssp | 4.x | Optional | Secondary structure |
| AlphaFold | 2.3+ | External | Structure prediction |

## Quick Start

### 1. Prepare your structures

Generate AlphaFold predictions for your proteins:
- **Monomer** predictions for each gene (PDB and/or mmCIF format)
- **Multimer** predictions for each protein complex of interest

### 2. Prepare your variant list

Create a CSV with columns: `gene`, `position`, `ref_aa`, `alt_aa`, plus optional annotation columns:

```csv
gene,position,ref_aa,alt_aa,alphamissense,franklin
SHROOM3,35,G,V,0.8766,VUS (mid)
ZIC3,350,R,G,0.998,VUS (mid)
```

### 3. Create a configuration file

Copy `config/example_config.yaml` and edit for your project:

```yaml
paths:
  working_dir: "/path/to/your/project"
  variants_file: "my_variants.csv"
  foldx_binary: "/path/to/foldx"

monomers:
  mygene:
    pdb: "fold_mygene_model_0.pdb"
    cif: "fold_mygene_model_0.cif"
    plddt_source: cif

multimers:
  - gene: mygene
    partner: mypartner
    pdb: "fold_mygene_mypartner_model_0.pdb"
    chain_gene: A
    chain_partner: B
    primary: true
```

### 4. Run the pipeline

```bash
varpat config/my_config.yaml
# or
python -m varpat config/my_config.yaml
```

## Configuration Reference

### Structure Definitions

**Monomers** — One entry per gene. Set `plddt_source: cif` if PDB B-factors are zeroed (common with ColabFold outputs).

**Multimers** — One entry per complex. Specify chain assignments (`chain_gene`, `chain_partner`) matching the PDB file.

### Construct Offsets

For genes where a multimer uses a truncated construct (e.g., a mature protein without signal peptide), specify offsets so variant positions can be mapped:

```yaml
construct_offsets:
  cdh2:
    cdh2_truncated: 159  # construct starts at full-length pos 160
```

Variants at positions ≤ offset are outside the construct and receive null structural values.

### Partner Normalization

FoldX labels partners by protein name. When running reciprocal complexes (e.g., a ROCK2 variant in a SHROOM3–ROCK2 complex), FoldX calls the partner "shroom3". The normalization section maps these back:

```yaml
partner_normalization:
  shroom3:
    rock2: rock2     # "shroom3" partner for ROCK2 variants → rock2
    dvl2: dvl2
```

### Thresholds

All scoring thresholds are configurable:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `plddt_confident` | 50 | Minimum pLDDT for DDG confidence |
| `ddg_destabilizing` | 1.0 | \|ΔΔG\| threshold for significance |
| `ddg_highly_destabilizing` | 2.0 | \|ΔΔG\| threshold for high significance |
| `am_pathogenic` | 0.564 | AlphaMissense pathogenic threshold |
| `contact_distance` | 5.0 Å | Residue-residue contact cutoff |

## Output Files

| File | Contents |
|------|----------|
| `variant_comprehensive.csv` | Complete table with all columns |
| `variant_comprehensive.xlsx` | Multi-sheet workbook (Summary, Monomer, Multimer, DDG, Concordance) |
| `high_priority_variants.csv` | Tier 1 and Tier 2 variants only |
| `foldx_ddg_expanded_results.csv` | DDG summary with flags |

## Three-Axis ΔΔG Framework

The pipeline computes three independent thermodynamic measures per variant:

| Axis | Measure | FoldX Command | What It Measures |
|------|---------|---------------|------------------|
| 1 | `ddg_monomer` | BuildModel (monomer) | Intrinsic fold stability |
| 2 | `ddg_fold_{partner}` | BuildModel (complex) | Fold stability in complex context |
| 3 | `ddg_binding_{partner}` | AnalyseComplex | Interaction energy change |

These can diverge: a mutation may destabilize the fold while paradoxically strengthening a specific protein–protein interaction ("conflicting" mechanism).

## Mechanism Categories

The pipeline classifies each evaluable variant into one of 16 mechanism categories based on fold stability × PPI effect:

**Fold disrupted:** Fold + PPI destabilization, Fold destab. at interface, Fold destabilization, Fold destab. + PPI stabilization (conflicting)

**Fold stabilized:** Stabilizing + PPI (GoF), Stabilizing at interface (GoF), Stabilizing (GoF), Stabilizing + PPI destab. (conflicting)

**Fold neutral:** PPI destabilization, PPI stabilization (GoF), Interface variant (DDG neutral), contact-driven, PPI conflicting, burial-driven, Benign, Structure unevaluable

## Four-Way Concordance

Agreement across four evidence axes (structure tier, ΔΔG, AlphaMissense, Franklin/ClinVar), reported as votes/evaluable_axes. Three flavors: strict, relaxed, and T3-inclusive.

## Citation

If you use this pipeline, please cite:

> [Your paper reference here]

## License

MIT License

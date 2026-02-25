# CHD Variant Structural Analysis Pipeline (v5.2)

A computational pipeline for assessing the structural and functional impact of missense variants in congenital heart defect (CHD) candidate genes using AlphaFold-predicted protein structures, FoldX thermodynamic modeling, and multi-evidence concordance analysis.

## Overview

This pipeline integrates four independent evidence lines to prioritize potentially pathogenic variants:

1. **Structural disruption scoring** — contact-based disruption, interface involvement, and burial depth from AlphaFold monomer and multimer structures
2. **FoldX ΔΔG thermodynamic modeling** — predicted destabilization of protein fold (monomer) and protein-protein interactions (multimer), with pLDDT-based confidence gating
3. **AlphaMissense** — deep learning pathogenicity predictions
4. **Franklin (Genoox)** — clinical variant classification

Variants are scored, tiered, and classified through a concordance framework that counts how many independent evidence lines support pathogenicity.

## Requirements

### Software

- Python ≥ 3.9
- [Biopython](https://biopython.org/) ≥ 1.80 (`pip install biopython`)
- pandas, numpy
- [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) (`mkdssp` binary for secondary structure assignment)
  - Install via `conda install -c salilab dssp` or `brew install dssp`
- [FoldX 5](https://foldxsuite.crg.eu/) (for ΔΔG calculations; requires academic license)

### Environment Variables (optional)

| Variable | Description | Default |
|---|---|---|
| `CHD_WORKING_DIR` | Root project directory | Current directory (`.`) |
| `DSSP_PATH` | Path to `mkdssp` binary | Auto-detected via `which mkdssp` |
| `FOLDX_PATH` | Path to FoldX binary | Auto-detected via `which foldx` |

## Directory Structure

```
project_root/
├── CHD_pipeline_v5_2.ipynb          # Pipeline notebook
├── variants_with_alphamissense_and_franklin_expanded.csv  # Input variants + annotations
├── foldx_ddg_monomer_results_all.csv   # Pre-computed monomer ΔΔG (optional)
├── foldx_ddg_multimer_results.csv      # Pre-computed multimer ΔΔG (optional)
├── structures/
│   ├── monomers/                    # AlphaFold monomer .pdb/.cif files
│   └── multimers/                   # AlphaFold multimer .pdb/.cif files
├── results/                         # Pipeline outputs (auto-created)
│   ├── variant_comprehensive_v5_2.csv
│   ├── high_priority_variants_v5_2.csv
│   └── variant_pipeline_results_summary_v5_2.csv
└── foldx_expanded/                  # FoldX working directory (auto-created)
```

## Pipeline Cells

| Cell | Name | Description |
|---|---|---|
| 0 | Header | Markdown documentation and changelog |
| 1 | Configuration | Paths, structure definitions, amino acid data, helper functions |
| 2 | Structure Functions | PDB/CIF loading, pLDDT extraction, contact counting, SASA, DSSP |
| 3 | Load Variants | Read input CSV, compute Grantham distances and severity |
| 4 | Monomer Metrics | Per-variant pLDDT, contacts, burial, secondary structure from monomers |
| 5 | Multimer Extraction | Per-complex metrics for each variant across all multimer partners |
| 6 | FoldX ΔΔG Loading | Load pre-computed ΔΔG, add confidence gating, classify DDG categories |
| 6a | FoldX Monomer ΔΔG | Run FoldX BuildModel for missing monomer ΔΔG values |
| 6b | FoldX Multimer ΔΔG | Run FoldX BuildModel + AnalyseComplex for missing multimer ΔΔG |
| 7 | Scoring & Tiers | Compute structural disruption scores, assign tiers, classify mechanisms |
| 8 | Annotations & Concordance | Normalize annotations, compute sub-scores and concordance |
| 9 | Save Results | Column ordering, output files, validation checks |

## Scoring Formula

```
final_score = (disruption_pts + interface_pts + burial_pts) × pLDDT_multiplier
```

| Component | Thresholds | Points |
|---|---|---|
| **Contact disruption** (severity × total contacts) | ≥20 | +4 |
| | ≥10 | +3 |
| | ≥4 | +2 |
| | ≥1 | +1 |
| **Interface partners** | ≥2 partners | +2 |
| | 1 partner | +1.5 |
| **Burial** (best across monomer + confident multimer) | buried_core | +2 |
| | partially_buried | +1 |
| **pLDDT multiplier** | ≥70 | ×1.0 |
| | 50–69 | ×0.7 |
| | <50 | ×0.4 |

**Maximum possible score: 8.0**

## Tier Assignments

| Tier | Score Range | Interpretation |
|---|---|---|
| Tier 1 — High confidence pathogenic | ≥ 5.0 | Strong structural evidence for disruption |
| Tier 2 — Likely pathogenic | 3.0 – 4.99 | Moderate structural evidence |
| Tier 3 — VUS with evidence | 1.5 – 2.99 | Weak structural signal, uncertain |
| Tier 4 — Likely benign | < 1.5 | No significant structural disruption |

## Concordance Framework

Four independent evidence lines are combined into concordance scores (0–4):

### Evidence Lines

| Line | Standard (strict) | Relaxed |
|---|---|---|
| **Structural tier** | Tier 1 or 2 → vote = 1 | Same |
| **FoldX ΔΔG** | max(\|mono\|, \|multi_max\|, \|multi_min\|) ≥ 2.0 AND `ddg_confidence` = high | ≥ 1.0 AND `ddg_confidence` ≠ low |
| **AlphaMissense** | likely_pathogenic → vote = 1 | likely_pathogenic OR ambiguous → vote = 1 |
| **Franklin** | pathogenic, LP, or VUS(high) → vote = 1 | + VUS(mid) → vote = 1 |

### DDG Confidence Gating

FoldX predictions on low-confidence AlphaFold structures are unreliable. The pipeline gates ΔΔG votes by pLDDT:

| `ddg_confidence` | pLDDT Range | Standard concordance | Relaxed concordance |
|---|---|---|---|
| high | ≥ 70 | ΔΔG vote counted | ΔΔG vote counted |
| moderate | 50–69 | ΔΔG vote **excluded** | ΔΔG vote counted |
| low | < 50 | ΔΔG vote **excluded** | ΔΔG vote **excluded** |

### Sub-Scores

The concordance is decomposed into two sub-scores for transparency:

- **`structure_strict`** / **`structure_relaxed`** (0–2): tier vote + ΔΔG vote
- **`external_strict`** / **`external_relaxed`** (0–2): AlphaMissense vote + Franklin vote

The four-way concordance = structure sub-score + external sub-score.

T3-inclusive variants (`_t3` suffix) also count Tier 3 as a structural vote.

## Integrated Classification

| Class | Criteria |
|---|---|
| **A — Concordant pathogenic** | Tier 1/2 AND AlphaMissense pathogenic |
| **B — Likely pathogenic (structural)** | Tier 1/2 AND AlphaMissense ambiguous AND interface variant |
| **C — Interface disruptor** | Tier 1/2 AND interface AND high confidence (no AM support) |
| **D — Structural only (high conf)** | Tier 1/2 AND high confidence (no AM or interface) |
| **E — Structural only (low conf)** | Tier 1/2 AND low confidence |
| **F — AM pathogenic (weak structural)** | Tier 3 AND AlphaMissense pathogenic |
| **G — VUS (mixed evidence)** | Tier 1/2 AND AlphaMissense ambiguous (not interface) |
| **H — VUS (weak evidence)** | Tier 1/2 with no AM support |
| **I — AlphaMissense only** | AlphaMissense pathogenic, Tier 3/4 |
| **J — Likely benign** | None of the above |

## Mechanism Classification

For variants with FoldX data, the pipeline classifies the predicted mechanism of disruption:

| Mechanism | Criteria |
|---|---|
| Dual mechanism (fold + PPI) | Tier 1/2, high monomer ΔΔG, interface, high multimer ΔΔG |
| Dual mechanism (fold + interface) | Tier 1/2, high monomer ΔΔG, interface present |
| Complex destabilization (PPI-specific) | Tier 1/2, neutral monomer ΔΔG, high multimer ΔΔG |
| Interface disruption (DDG-neutral) | Tier 1/2, interface present, low monomer ΔΔG |
| Fold destabilization | Tier 1/2, high monomer ΔΔG, no interface |
| Structural concern (moderate DDG) | Tier 1/2, moderate ΔΔG (0.5–1.5 kcal/mol) |
| Structural tier (low DDG) | Tier 1/2, ΔΔG below moderate threshold |
| High DDG only (Tier 3/4) | Low structural tier but high ΔΔG |
| Likely benign | Low tier, low ΔΔG |

## Output Column Reference

### Variant Identification

| Column | Description |
|---|---|
| `gene` | Gene symbol (lowercase) |
| `position` | Amino acid position (1-indexed) |
| `ref_aa` | Reference amino acid (single letter) |
| `alt_aa` | Alternate amino acid (single letter) |

### Amino Acid Properties

| Column | Description |
|---|---|
| `grantham_distance` | Grantham distance between ref and alt amino acids (0–215) |
| `grantham_class` | conservative / moderately_conservative / moderately_radical / radical |
| `substitution_severity` | Normalized Grantham score (0–4.0), used in disruption calculation |
| `property_changes` | Changed biophysical properties (size, charge, hydrophobicity) |

### Monomer Structural Metrics

| Column | Description |
|---|---|
| `monomer_plddt` | AlphaFold pLDDT at variant position (0–100) |
| `monomer_plddt_category` | very_high (≥90) / confident (≥70) / low (≥50) / very_low (<50) |
| `monomer_n_contacts` | Number of intra-chain residue contacts within 5Å (sequence separation ≥3) |
| `monomer_contact_category` | low_contact / medium_contact / high_contact |
| `monomer_aa` | Amino acid observed in the structure at this position |
| `monomer_accessibility` | Relative solvent accessible surface area (0–1.0) |
| `monomer_burial` | buried_core (<0.05) / partially_buried (<0.25) / surface_exposed |
| `monomer_secondary_structure` | DSSP assignment: H (helix), E (sheet), T (turn), S (bend), - (coil), etc. |
| `monomer_contact_disruption` | = monomer_n_contacts (used in disruption formula) |

### Multimer Complex Metrics

For each interaction partner, 8 columns are produced with prefix `multi_{partner}_`:

| Suffix | Description |
|---|---|
| `_plddt` | pLDDT at variant position in the multimer complex |
| `_n_contacts` | Total residue contacts (intra + inter chain) |
| `_inter_contacts` | Inter-chain contacts with the partner |
| `_is_interface` | True if inter_contacts > 0 |
| `_accessibility` | Relative SASA in the complex |
| `_burial` | Burial classification in the complex |
| `_sec_struct` | Secondary structure in the complex |
| `_disruption` | = n_contacts (used in disruption formula) |

### Multimer Summary

| Column | Description |
|---|---|
| `n_multimer_complexes` | Number of multimer complexes analyzed |
| `multimer_partners` | Semicolon-separated list of partner proteins |
| `is_interface_any` | True if variant is at interface in any complex |
| `interface_partners` | Partners where variant is at the interface |
| `n_interface_partners` | Number of interface partners |
| `multimer_plddt_max` / `_avg` | Max/mean pLDDT across all complexes |
| `multimer_contacts_max` / `_avg` | Max/mean contacts across all complexes |
| `multimer_disruption_max` / `_avg` | Max/mean disruption across all complexes |

### Structural Scoring

| Column | Description |
|---|---|
| `contact_disruption` | substitution_severity × total_contacts |
| `total_contacts` | monomer_contacts + sum of all inter-chain contacts |
| `inter_contacts_sum` | Total inter-chain contacts across all complexes |
| `best_burial` | Most buried classification across monomer + confident multimers |
| `best_burial_source` | Which structure provided the best burial |
| `final_score` | Composite structural disruption score (0–8.0) |
| `confidence` | high (pLDDT ≥ 70) or low (< 70) |
| `best_plddt` | Highest pLDDT across monomer and all multimers |
| `score_evidence` | Human-readable breakdown of scoring components |
| `tier` | Tier 1–4 classification (see Tier Assignments) |

### External Annotations

| Column | Description |
|---|---|
| `AlphaMissense` | Normalized category: likely_pathogenic / ambiguous / likely_benign |
| `AlphaMissense_pathogenicity` | Raw AlphaMissense score (0–1) |
| `AlphaMissense_raw` | Original value before normalization |
| `franklin` | Normalized Franklin classification |
| `franklin_raw` | Original Franklin value before normalization |
| `integrated_class` | Classes A–J integrating tier + AlphaMissense (see Integrated Classification) |

### FoldX ΔΔG

| Column | Description |
|---|---|
| `ddg_monomer` | Monomer ΔΔG in kcal/mol (positive = destabilizing) |
| `ddg_category` | Confidence-aware category (may append `_unreliable` at low pLDDT) |
| `ddg_category_raw` | Category before confidence adjustment |
| `ddg_confidence` | high (pLDDT ≥ 70) / moderate (50–69) / low (< 50) |
| `ddg_multimer_max` | Largest (most positive) multimer ΔΔG across complexes |
| `ddg_multimer_min` | Smallest (most negative) multimer ΔΔG across complexes |
| `ddg_multimer_mean` | Mean multimer ΔΔG across tested complexes |
| `n_complexes_tested` | Number of multimer complexes with ΔΔG data |
| `partners_tested` | Which partners had ΔΔG computed |

### Mechanism

| Column | Description |
|---|---|
| `pathogenic_mechanism` | Predicted disruption mechanism (see Mechanism Classification) |
| `final_mechanism` | Same as pathogenic_mechanism |
| `evidence_summary` | Compact evidence string (e.g., `DDG_high(5.2); interface(rock2); buried; AM_path`) |

### Concordance Sub-Scores (new in v5.2)

| Column | Range | Description |
|---|---|---|
| `structure_strict` | 0–2 | Tier vote (T1/T2 = 1) + DDG strict vote |
| `structure_relaxed` | 0–2 | Tier vote (T1/T2 = 1) + DDG relaxed vote |
| `structure_strict_t3` | 0–2 | Tier vote (T1/T2/T3 = 1) + DDG strict vote |
| `structure_relaxed_t3` | 0–2 | Tier vote (T1/T2/T3 = 1) + DDG relaxed vote |
| `external_strict` | 0–2 | AM strict vote + Franklin strict vote |
| `external_relaxed` | 0–2 | AM relaxed vote + Franklin relaxed vote |

### Concordance Totals

| Column | Range | Description |
|---|---|---|
| `three_way` | 0–3 | Tier + AM strict + Franklin strict (no DDG) |
| `three_way_ambiguous` | 0–3 | Tier + AM relaxed + Franklin relaxed (no DDG) |
| `three_way_t3` | 0–3 | T3-inclusive three-way |
| `four_way` | 0–4 | structure_strict + external_strict |
| `four_way_ambiguous_as_pathogenic_ddg_1_threshold` | 0–4 | structure_relaxed + external_relaxed |
| `four_way_t3` | 0–4 | T3-inclusive four-way strict |
| `four_way_t3_ambiguous` | 0–4 | T3-inclusive four-way relaxed |

## Version History

| Version | Date | Changes |
|---|---|---|
| **v5.2** | Feb 2026 | DDG concordance bug fix (multi_max now included in threshold check); new structure/external sub-score columns |
| v5.1 | Feb 2026 | pLDDT confidence gating on DDG concordance; ddg_multimer_min column; T3-inclusive concordance; ddg_category reclassification at low pLDDT |
| v5 | Feb 2026 | Disruption = severity × contacts; graduated pLDDT multiplier; best burial across confident multimers; 4-tier system |

## Citation

If you use this pipeline, please cite the relevant tools:

- **AlphaFold**: Jumper, J. et al. (2021). Nature, 596, 583–589.
- **FoldX**: Schymkowitz, J. et al. (2005). Nucleic Acids Res., 33, W382–W388.
- **AlphaMissense**: Cheng, J. et al. (2023). Science, 381, eadg7492.
- **Biopython**: Cock, P.J.A. et al. (2009). Bioinformatics, 25, 1422–1423.

## License

This pipeline is provided for academic research purposes.

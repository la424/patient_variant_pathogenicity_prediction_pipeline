# MAVIS — Multimer-Aware Variant Impact Scoring

**Systematic structural evaluation of missense variant pathogenicity with protein–protein interaction interface analysis.**

MAVIS takes AlphaFold-predicted structures and a list of missense variants, runs FoldX thermodynamic calculations, and produces a multi-axis pathogenicity assessment that accounts for effects at protein–protein interaction interfaces — something monomer-only tools miss.

---

## Why MAVIS?

Most variant pathogenicity tools evaluate proteins in isolation. But many disease-relevant mutations exert their effects at protein–protein interaction interfaces — disrupting binding, altering complex stability, or producing gain-of-function through interface strengthening. MAVIS addresses this gap by:

- **Evaluating variants in both monomer and multimer structural contexts**
- **Separating three independent thermodynamic axes:** intrinsic fold stability, complex fold stability, and binding energy
- **Detecting conflicting signals** (e.g., fold destabilization + interface stabilization) that indicate biologically meaningful mechanisms
- **Gating all assessments by AlphaFold confidence (pLDDT)** to prevent artifacts in disordered regions

## Two Pipelines

| Pipeline | Notebook | What it does |
|----------|----------|-------------|
| **Structural** | `mavis_structural.ipynb` | Structure + FoldX only. Takes a variant CSV and AlphaFold structures, runs FoldX, and produces structural tier assignments with 16-category mechanism classification. No external databases needed. |
| **Full** | `mavis_full.ipynb` | Everything in Structural, plus four-way concordance with AlphaMissense and Franklin/ClinVar annotations, dual-pipeline scoring (single-residue + neighborhood ±3), and pipeline agreement tracking. |

**Start with Structural** if you just want to assess your variants structurally. Use **Full** if you also have AlphaMissense scores and clinical annotations and want concordance analysis.

---

## Quick Start

### 1. Install dependencies

```bash
pip install biopython pandas numpy openpyxl
```

You also need:
- **FoldX 5** ([academic license](https://foldxsuite.crg.eu/)) — required for ΔΔG calculations
- **mkdssp** (optional) — for secondary structure assignment

### 2. Prepare your input

**Variant CSV** (`variants.csv`):
```csv
gene,position,ref_aa,alt_aa
SHROOM3,35,G,V
SHROOM3,60,G,V
ZIC3,402,S,P
```

Or with a combined variant column:
```csv
gene,variant
SHROOM3,G35V
ZIC3,S402P
```

**Structure files** — Place AlphaFold PDB/CIF files in a `structures/` directory:
```
structures/
├── fold_shroom3_model_0.pdb          # monomer
├── fold_shroom3_model_0.cif          # monomer (CIF for pLDDT)
├── fold_shroom3_rock2_model_0.pdb    # multimer complex
└── fold_zic3_gli3_model_0.pdb        # multimer complex
```

### 3. Configure and run

Open the notebook, edit Cell 1:

**Option A — Auto-detection** (easiest):
```python
STRUCTURE_DIR = Path("structures/")
AUTO_DETECT = True
```

**Option B — Manual configuration** (most control):
```python
AUTO_DETECT = False
MONOMER_STRUCTURES = {
    'shroom3': ('fold_shroom3_model_0.cif', 'fold_shroom3_model_0.pdb'),
    'zic3':    ('fold_zic3_model_0.cif',    'fold_zic3_model_0.pdb'),
}
MULTIMER_STRUCTURES = [
    ('shroom3', 'rock2', None, 'fold_shroom3_rock2_model_0.pdb', 'A', 'B', True),
    ('zic3',    'gli3',  None, 'fold_zic3_gli3_model_0.pdb',     'A', 'B', True),
]
```

Run all cells. Results appear in `results/`.

### 4. Interpret results

Key output columns:

| Column | Meaning |
|--------|---------|
| `v6_tier` | Structural impact tier (1 = highest, 4 = minimal) |
| `v6_mechanism` | Pathogenic mechanism (16 categories) |
| `ddg_monomer` | Monomer fold ΔΔG (kcal/mol) |
| `ddg_binding_{partner}` | Binding energy change per interaction partner |
| `ddg_fold_{partner}` | Complex fold ΔΔG per partner |
| `concordance_strict` | Four-way concordance score (Full pipeline only) |

---

## Pipeline Architecture

### Three-Axis ΔΔG Framework

```
Axis 1: ddg_monomer          Intrinsic fold stability (BuildModel on monomer)
Axis 2: ddg_fold_{partner}   Fold stability in complex context (BuildModel on complex)  
Axis 3: ddg_binding_{partner} Interaction energy change (AnalyseComplex: ΔIE)
```

The separation of Axis 2 and 3 is biologically critical — a variant can destabilize the local fold while strengthening interface contacts, producing a conflicting signal that indicates a distinct mechanism.

### 16-Category Mechanism Classification

Variants are classified by the interaction of fold stability, PPI energy, interface status, and structural context:

**Fold destabilized** (ΔΔG > +1.0): Fold + PPI destabilization, Fold destab. at interface, Fold destabilization, Fold destab. + PPI stabilization (conflicting)

**Fold stabilized** (ΔΔG < −1.0): Stabilizing + PPI (potential GoF), Stabilizing at interface, Stabilizing (potential GoF), Stabilizing + PPI destabilization (conflicting)

**Fold neutral**: PPI destabilization, PPI stabilization (potential GoF), Interface variant (DDG neutral), PPI conflicting (mixed partners), Structural variant – contact/burial-driven, Benign (structurally evaluated)

**Unevaluable**: Structure unevaluable (pLDDT < 50)

### Dual Scoring Pipelines (Full pipeline)

- **Pipeline 1:** Single-residue structural disruption score
- **Pipeline 2:** Neighborhood ±3 weighted contact score (robust to single-position artifacts)
- **Pipeline agreement** tracks concordance between the two views

### Four-Way Concordance (Full pipeline)

Each variant is assessed across four independent axes: structural tier, FoldX ΔΔG, AlphaMissense, and Franklin/ClinVar. The concordance score (X/N) uses adjusted denominators reflecting only evaluable axes.

---

## Output Files

| File | Contents |
|------|----------|
| `variant_comprehensive.csv` | All variants × all columns |
| `high_priority_variants.csv` | Tier 1 and 2 variants only |
| `foldx_ddg_expanded_results.csv` | ΔΔG summary per variant |
| `variant_comprehensive.xlsx` | Multi-sheet workbook (Summary, Monomer Detail, Multimer Detail, DDG Detail, Concordance Detail) |

---

## Structure File Conventions

### Auto-detection naming

MAVIS auto-detects structures using AlphaFold naming conventions:

- **Monomer:** `fold_{gene}_model_N.pdb` or `{gene}.pdb`
- **Multimer:** `fold_{gene1}_{gene2}_model_N.pdb`
- CIF files (`.cif`) are preferred for pLDDT extraction when PDB B-factors are zeroed

### Construct offsets

If a multimer construct uses a truncated protein (e.g., starting at residue 160 of the full-length sequence), define the offset in `CONSTRUCT_OFFSETS`:

```python
CONSTRUCT_OFFSETS = {
    'partner_label': 159,  # full-length pos 160 = construct pos 1
}
```

Variants at positions ≤ offset are automatically nullified for that construct.

---

## Configuration Reference

All parameters are set in Cell 1 of either notebook:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `WORKING_DIR` | `.` | Root project directory |
| `FOLDX_BINARY` | `foldx` | Path to FoldX 5 executable |
| `DSSP_PATH` | `mkdssp` | Path to mkdssp (or `None`) |
| `CONTACT_DISTANCE` | 5.0 | Å cutoff for residue contacts |
| `PLDDT_GATE_STRICT` | 70 | pLDDT threshold for strict ΔΔG concordance |
| `PLDDT_GATE_RELAXED` | 50 | pLDDT threshold for evaluability |
| `DDG_DESTAB` | 1.0 | ΔΔG threshold for destabilizing (kcal/mol) |
| `DDG_HIGHLY` | 2.0 | ΔΔG threshold for highly destabilizing |
| `FOLDX_N_RUNS` | 3 | Number of FoldX BuildModel replicates |
| `DRY_RUN` | `False` | Skip FoldX execution (for testing) |

---

## Citation

If you use MAVIS in your research, please cite:

> [Citation details to be added upon publication]

---

## License

MIT License. See [LICENSE](LICENSE) for details.

FoldX requires a separate academic license from the [FoldX Suite](https://foldxsuite.crg.eu/).

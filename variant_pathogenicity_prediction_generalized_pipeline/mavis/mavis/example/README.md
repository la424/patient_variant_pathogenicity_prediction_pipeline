# Example Data

This directory contains an example variant CSV to demonstrate the expected input format.

## Obtaining Structure Files

MAVIS requires AlphaFold-predicted structures in PDB or CIF format. These are not included in the repository due to file size.

### Monomer structures

Download from the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/):

1. Search for your protein by UniProt ID or gene name
2. Download the PDB or CIF file
3. Rename to match the expected convention: `fold_{gene}_model_0.pdb`

### Multimer structures

Generate predictions using the [AlphaFold Server](https://alphafoldserver.com/):

1. Submit your protein complex (two or more chains)
2. Download the predicted structure
3. Place in the `structures/` directory with naming: `fold_{gene1}_{gene2}_model_0.pdb`

### Directory setup

```
your_project/
├── variants.csv
├── structures/
│   ├── fold_geneA_model_0.pdb
│   ├── fold_geneA_model_0.cif
│   ├── fold_geneA_geneB_model_0.pdb
│   └── fold_geneA_geneB_model_0.cif
├── mavis_structural.ipynb    (copy from repo)
└── mavis_full.ipynb          (copy from repo)
```

## Running the example

1. Download AlphaFold structures for the genes in `variants_example.csv`
2. Place them in a `structures/` subdirectory
3. Copy a MAVIS notebook into the same directory
4. Edit Cell 1 to set `VARIANTS_FILE = Path("variants_example.csv")`
5. Run all cells

## Input format options

MAVIS accepts two CSV formats:

**Separate columns (recommended):**
```csv
gene,position,ref_aa,alt_aa
SHROOM3,35,G,V
```

**Combined variant column:**
```csv
gene,variant
SHROOM3,G35V
```

Any additional columns (e.g., `alphamissense`, `franklin`) are preserved and used by the full pipeline for concordance analysis.

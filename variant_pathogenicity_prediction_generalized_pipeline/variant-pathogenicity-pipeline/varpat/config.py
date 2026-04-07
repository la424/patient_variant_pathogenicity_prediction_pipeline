"""Configuration loader and validator."""

from pathlib import Path
from typing import Dict, List, Optional, Any
import yaml


class PipelineConfig:
    """Structured access to pipeline configuration."""

    def __init__(self, config_path: str):
        with open(config_path) as f:
            self._raw = yaml.safe_load(f)
        self._validate()
        self._resolve_paths()

    def _validate(self):
        required = ['paths', 'monomers', 'multimers']
        for key in required:
            if key not in self._raw:
                raise ValueError(f"Missing required config section: '{key}'")

    def _resolve_paths(self):
        p = self._raw['paths']
        self.working_dir = Path(p['working_dir']).expanduser().resolve()
        self.results_dir = self._resolve(p.get('results_dir', 'results'))
        self.results_dir.mkdir(parents=True, exist_ok=True)

        self.structures_dir = self._resolve(p.get('structures_dir', 'structures'))
        self.variants_file = self._resolve(p['variants_file'])

        self.monomer_ddg_file = self._resolve_opt(p.get('monomer_ddg_file'))
        self.multimer_ddg_file = self._resolve_opt(p.get('multimer_ddg_file'))
        self.per_partner_ddg_file = self._resolve_opt(p.get('per_partner_ddg_file'))

        self.foldx_binary = Path(p['foldx_binary']) if p.get('foldx_binary') else None
        self.dssp_binary = p.get('dssp_binary', 'mkdssp')

        # Build search dirs for structure files
        self.search_dirs = [
            self.working_dir,
            self.structures_dir,
            self.structures_dir / 'monomers',
            self.structures_dir / 'multimers',
        ]

    def _resolve(self, rel_path) -> Path:
        p = Path(rel_path)
        if p.is_absolute():
            return p
        return self.working_dir / p

    def _resolve_opt(self, rel_path) -> Optional[Path]:
        if rel_path is None:
            return None
        p = self._resolve(rel_path)
        return p if p.exists() else None

    def find_file(self, filename: Optional[str]) -> Optional[Path]:
        """Search configured directories for a file."""
        if filename is None:
            return None
        for d in self.search_dirs:
            p = d / filename
            if p.exists():
                return p
        return None

    # ─── Monomer definitions ──────────────────────────────────────────────
    @property
    def monomer_genes(self) -> List[str]:
        return list(self._raw.get('monomers', {}).keys())

    def monomer_files(self, gene: str) -> dict:
        """Return {cif: path_or_None, pdb: path_or_None, plddt_source: str}."""
        m = self._raw.get('monomers', {}).get(gene, {})
        return {
            'cif': self.find_file(m.get('cif')),
            'pdb': self.find_file(m.get('pdb')),
            'plddt_source': m.get('plddt_source', 'pdb'),
        }

    # ─── Multimer definitions ─────────────────────────────────────────────
    @property
    def multimer_complexes(self) -> List[dict]:
        return self._raw.get('multimers', [])

    @property
    def partner_labels(self) -> List[str]:
        """All unique partner labels across multimer definitions."""
        labels = []
        seen = set()
        for m in self.multimer_complexes:
            pl = m['partner']
            if pl not in seen:
                labels.append(pl)
                seen.add(pl)
        return labels

    # ─── Construct offsets ────────────────────────────────────────────────
    def get_offset(self, gene: str, partner_label: str) -> int:
        """Return construct offset or 0 if none defined."""
        offsets = self._raw.get('construct_offsets', {})
        gene_offsets = offsets.get(gene.lower(), {})
        return gene_offsets.get(partner_label.lower(), 0)

    def apply_offset(self, gene: str, partner_label: str, position: int):
        """Apply construct offset. Returns (corrected_pos, is_valid)."""
        offset = self.get_offset(gene, partner_label)
        if offset == 0:
            return position, True
        corrected = position - offset
        return (corrected, True) if corrected > 0 else (None, False)

    # ─── Partner normalization ────────────────────────────────────────────
    def normalize_partner(self, foldx_partner_name: str, variant_gene: str) -> str:
        """Map FoldX partner name to structural extraction label."""
        norms = self._raw.get('partner_normalization', {})
        p = foldx_partner_name.strip().lower()
        g = variant_gene.strip().lower()

        if p in norms:
            mapping = norms[p]
            if g in mapping:
                return mapping[g]
            if '_default' in mapping:
                return mapping['_default']
        return p

    def normalize_partners_tested(self, partners_str: str, gene: str) -> str:
        """Normalize a semicolon-separated partners_tested string."""
        if not partners_str or str(partners_str).strip() in ('', 'nan'):
            return partners_str
        parts = [p.strip() for p in str(partners_str).split(';') if p.strip()]
        normalized = [self.normalize_partner(p, gene) for p in parts]
        seen = set()
        result = []
        for p in normalized:
            if p not in seen:
                seen.add(p)
                result.append(p)
        return ';'.join(result)

    # ─── Column mapping ──────────────────────────────────────────────────
    @property
    def column_map(self) -> Dict[str, Optional[str]]:
        defaults = {
            'gene': 'gene', 'position': 'position',
            'ref_aa': 'ref_aa', 'alt_aa': 'alt_aa',
        }
        cm = self._raw.get('columns', {})
        defaults.update({k: v for k, v in cm.items() if v is not None})
        return defaults

    # ─── Thresholds ───────────────────────────────────────────────────────
    @property
    def thresholds(self) -> dict:
        defaults = {
            'plddt_confident': 50,
            'plddt_high': 70,
            'plddt_very_high': 90,
            'ddg_destabilizing': 1.0,
            'ddg_highly_destabilizing': 2.0,
            'disruption_points': [[20, 4.0], [10, 3.0], [4, 2.0], [1, 1.0]],
            'tiers': [[5.0, 'Tier 1'], [3.0, 'Tier 2'], [1.5, 'Tier 3']],
            'default_tier': 'Tier 4',
            'contact_driven_threshold': 6,
            'am_pathogenic': 0.564,
            'am_benign': 0.340,
            'neighborhood_radius': 3,
            'neighborhood_weights': {0: 1.0, 1: 0.75, 2: 0.50, 3: 0.25},
            'contact_distance': 5.0,
            'sequence_separation': 3,
        }
        defaults.update(self._raw.get('thresholds', {}))
        return defaults

    # ─── FoldX settings ──────────────────────────────────────────────────
    @property
    def foldx_settings(self) -> dict:
        defaults = {
            'n_runs': 3,
            'timeout': 300,
            'run_missing_monomer': True,
            'run_missing_multimer': True,
        }
        defaults.update(self._raw.get('foldx', {}))
        return defaults

    # ─── Output settings ─────────────────────────────────────────────────
    @property
    def output_settings(self) -> dict:
        defaults = {
            'csv_filename': 'variant_comprehensive.csv',
            'xlsx_filename': 'variant_comprehensive.xlsx',
            'high_priority_filename': 'high_priority_variants.csv',
            'ddg_expanded_filename': 'foldx_ddg_expanded_results.csv',
            'xlsx_sheets': True,
            'include_legacy_columns': False,
        }
        defaults.update(self._raw.get('output', {}))
        return defaults

    # ─── Reverse partner map ──────────────────────────────────────────────
    @property
    def partner_gene_map(self) -> Dict[str, str]:
        """Map partner labels → gene names for reciprocal extraction."""
        result = {}
        for m in self.multimer_complexes:
            partner = m['partner'].lower()
            gene = m['gene'].lower()
            # The partner label typically IS the gene name, but for special
            # constructs (cdh2_truncated, cdh2_cyto, actb_no_bind) we need
            # to map to the base gene.
            base = partner.split('_')[0]  # cdh2_truncated → cdh2
            if base in [mg.lower() for mg in self.monomer_genes]:
                result[partner] = base
            elif partner in [mg.lower() for mg in self.monomer_genes]:
                result[partner] = partner
        return result

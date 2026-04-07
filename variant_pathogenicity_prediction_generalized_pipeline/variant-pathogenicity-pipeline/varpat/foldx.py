"""FoldX thermodynamic calculations: BuildModel and AnalyseComplex."""

import shutil
import subprocess
from pathlib import Path
from typing import Optional, Tuple

from .config import PipelineConfig


def run_buildmodel(config: PipelineConfig, structure_path: Path,
                   chain_id: str, ref_aa: str, position: int,
                   alt_aa: str, work_dir: Path) -> Optional[float]:
    """Run FoldX BuildModel. Returns average DDG or None on failure."""
    if config.foldx_binary is None or not config.foldx_binary.exists():
        return None

    foldx = config.foldx_binary
    n_runs = config.foldx_settings['n_runs']
    timeout = config.foldx_settings['timeout']
    rotabase = foldx.parent / 'rotabase.txt'

    work_dir.mkdir(parents=True, exist_ok=True)
    struct_name = structure_path.name
    shutil.copy2(structure_path, work_dir / struct_name)
    if rotabase.exists() and not (work_dir / 'rotabase.txt').exists():
        shutil.copy2(rotabase, work_dir / 'rotabase.txt')

    mut_str = f"{ref_aa}{chain_id}{position}{alt_aa};"
    (work_dir / 'individual_list.txt').write_text(mut_str + '\n')

    cmd = [
        str(foldx), '--command=BuildModel',
        f'--pdb={struct_name}',
        '--mutant-file=individual_list.txt',
        f'--numberOfRuns={n_runs}',
        f'--output-dir={work_dir}',
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=timeout, cwd=str(work_dir))
        if result.returncode != 0:
            return None

        # Parse Dif_*.fxout
        dif_file = None
        for f in work_dir.glob('Dif_*.fxout'):
            dif_file = f
            break
        if dif_file is None:
            return None

        ddg_values = []
        with open(dif_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('Pdb') or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    try:
                        ddg_values.append(float(parts[1]))
                    except ValueError:
                        continue
        return round(sum(ddg_values) / len(ddg_values), 4) if ddg_values else None

    except (subprocess.TimeoutExpired, Exception):
        return None


def run_analysecomplex(config: PipelineConfig, structure_path: Path,
                       chains: str, work_dir: Path) -> Optional[float]:
    """Run FoldX AnalyseComplex. Returns interaction energy or None."""
    if config.foldx_binary is None:
        return None

    foldx = config.foldx_binary
    timeout = config.foldx_settings['timeout']
    rotabase = foldx.parent / 'rotabase.txt'

    work_dir.mkdir(parents=True, exist_ok=True)
    struct_name = structure_path.name
    if not (work_dir / struct_name).exists():
        shutil.copy2(structure_path, work_dir / struct_name)
    if rotabase.exists() and not (work_dir / 'rotabase.txt').exists():
        shutil.copy2(rotabase, work_dir / 'rotabase.txt')

    cmd = [
        str(foldx), '--command=AnalyseComplex',
        f'--pdb={struct_name}',
        f'--analyseComplexChains={chains}',
        f'--output-dir={work_dir}',
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=timeout, cwd=str(work_dir))
        if result.returncode != 0:
            return None
        for f in work_dir.glob('Interaction_*_AC.fxout'):
            with open(f) as fh:
                for line in fh:
                    if line.startswith('Pdb') or line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        try:
                            return float(parts[5])
                        except ValueError:
                            continue
        return None
    except Exception:
        return None


def run_multimer_ddg(config: PipelineConfig, pdb_path: Path,
                     chain_gene: str, chain_partner: str,
                     ref_aa: str, position: int, alt_aa: str,
                     work_dir: Path) -> Tuple[Optional[float], Optional[float]]:
    """Compute multimer DDG: BuildModel + AnalyseComplex on WT and mutant.
    Returns (ddg_fold, ddg_binding)."""
    chains_str = f"{chain_gene},{chain_partner}"

    # WT interaction energy
    wt_dir = work_dir / 'wt'
    wt_dir.mkdir(exist_ok=True)
    ie_wt = run_analysecomplex(config, pdb_path, chains_str, wt_dir)

    # BuildModel (fold DDG)
    ddg_fold = run_buildmodel(config, pdb_path, chain_gene, ref_aa,
                              position, alt_aa, work_dir)

    # Find mutant PDB → AnalyseComplex
    mutant_pdb = None
    for f in work_dir.glob(f"{pdb_path.stem}_1_*.pdb"):
        mutant_pdb = f
        break
    if mutant_pdb is None:
        for f in work_dir.glob('*_1.pdb'):
            mutant_pdb = f
            break

    ie_mut = None
    if mutant_pdb is not None:
        mut_dir = work_dir / 'mut'
        mut_dir.mkdir(exist_ok=True)
        ie_mut = run_analysecomplex(config, mutant_pdb, chains_str, mut_dir)

    ddg_binding = None
    if ie_wt is not None and ie_mut is not None:
        ddg_binding = round(ie_mut - ie_wt, 4)

    return ddg_fold, ddg_binding

"""Structure loading, pLDDT extraction, contact counting, SASA, and DSSP."""

import subprocess
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch, ShrakeRupley

from .constants import THREE_TO_ONE, MAX_SASA

_pdb_parser = PDBParser(QUIET=True)
_cif_parser = MMCIFParser(QUIET=True)


def load_cif(path: Optional[Path]):
    if path and path.exists():
        try:
            return _cif_parser.get_structure('s', str(path))
        except Exception:
            pass
    return None


def load_pdb(path: Optional[Path]):
    if path and path.exists():
        try:
            return _pdb_parser.get_structure('s', str(path))
        except Exception:
            pass
    return None


def get_plddt(structure, chain_id: str = 'A') -> Dict[int, float]:
    """Per-residue pLDDT from B-factors. Returns only non-zero values."""
    plddt = {}
    if structure is None:
        return plddt
    model = structure[0]
    if chain_id not in model:
        for c in model:
            chain_id = c.id
            break
    if chain_id not in model:
        return plddt
    for res in model[chain_id].get_residues():
        if res.id[0] == ' ':
            p = None
            if 'CA' in res:
                p = res['CA'].bfactor
            else:
                for atom in res:
                    p = atom.bfactor
                    break
            if p is not None and p > 0:
                plddt[res.id[1]] = round(p, 2)
    return plddt


def get_residue_aa(structure, chain_id: str = 'A') -> Dict[int, str]:
    if structure is None:
        return {}
    model = structure[0]
    if chain_id not in model:
        for c in model:
            chain_id = c.id
            break
    if chain_id not in model:
        return {}
    return {
        r.id[1]: THREE_TO_ONE.get(r.resname, '?')
        for r in model[chain_id].get_residues() if r.id[0] == ' '
    }


def count_contacts(structure, chain_id: str = 'A',
                   distance: float = 5.0, seq_sep: int = 3) -> Dict[int, int]:
    """Intra-chain residue-residue contacts with sequence separation filter."""
    if structure is None:
        return {}
    model = structure[0]
    if chain_id not in model:
        for c in model:
            chain_id = c.id
            break
    if chain_id not in model:
        return {}

    residues = [r for r in model[chain_id].get_residues() if r.id[0] == ' ']
    contacts = {}
    for i, res in enumerate(residues):
        neighbor_set = set()
        for j, other in enumerate(residues):
            if abs(i - j) < seq_sep:
                continue
            for atom_i in res.get_atoms():
                found = False
                for atom_j in other.get_atoms():
                    if atom_i - atom_j < distance:
                        neighbor_set.add(other.id[1])
                        found = True
                        break
                if found:
                    break
        contacts[res.id[1]] = len(neighbor_set)
    return contacts


def count_interface(structure, my_chain: str, partner_chain: str,
                    distance: float = 5.0) -> Tuple[Dict[int, int], Set[int]]:
    """Inter-chain contacts. Returns (contact_map, interface_residues)."""
    if structure is None:
        return {}, set()
    model = structure[0]
    if my_chain not in model or partner_chain not in model:
        return {}, set()

    partner_atoms = list(model[partner_chain].get_atoms())
    if not partner_atoms:
        return {}, set()
    ns = NeighborSearch(partner_atoms)

    inter, iface = {}, set()
    for res in model[my_chain].get_residues():
        if res.id[0] != ' ':
            continue
        partner_residues = set()
        for atom in res.get_atoms():
            for nb in ns.search(atom.coord, distance, 'R'):
                if nb.id[0] == ' ':
                    partner_residues.add(nb.id[1])
        count = len(partner_residues)
        if count > 0:
            inter[res.id[1]] = count
            iface.add(res.id[1])
    return inter, iface


def get_accessibility(structure, chain_id: str = 'A') -> Dict[int, float]:
    """Relative solvent accessibility via Shrake-Rupley."""
    acc = {}
    if structure is None:
        return acc
    model = structure[0]
    if chain_id not in model:
        for c in model:
            chain_id = c.id
            break
    if chain_id not in model:
        return acc
    try:
        sr = ShrakeRupley()
        sr.compute(model, level='R')
        for res in model[chain_id].get_residues():
            if res.id[0] != ' ':
                continue
            aa = THREE_TO_ONE.get(res.resname, 'X')
            max_s = MAX_SASA.get(aa, 200)
            if max_s > 0 and hasattr(res, 'sasa'):
                acc[res.id[1]] = round(min(1.0, res.sasa / max_s), 9)
    except Exception:
        pass
    return acc


def get_secondary_structure(pdb_path: Optional[Path], chain_id: str = 'A',
                            dssp_binary: str = 'mkdssp') -> Dict[int, str]:
    """Secondary structure via DSSP."""
    ss_map = {}
    if pdb_path is None or not pdb_path.exists():
        return ss_map
    try:
        for cmd in [
            [dssp_binary, '--output-format', 'dssp', str(pdb_path)],
            [dssp_binary, '-i', str(pdb_path)],
            [dssp_binary, str(pdb_path)],
        ]:
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
                if result.returncode == 0 and len(result.stdout) > 100:
                    break
            except FileNotFoundError:
                return ss_map
        else:
            return ss_map
        if result.returncode != 0:
            return ss_map

        in_data = False
        for line in result.stdout.split('\n'):
            if '  #  RESIDUE' in line:
                in_data = True
                continue
            if not in_data or len(line) < 17:
                continue
            if line[13] == '!':
                continue
            try:
                chain = line[11]
                if chain != chain_id:
                    continue
                resnum = int(line[5:10].strip())
                sec = line[16] if len(line) > 16 and line[16] != ' ' else '-'
                ss_map[resnum] = sec
            except (ValueError, IndexError):
                continue
    except Exception:
        pass
    return ss_map


def load_monomer(config, gene: str):
    """Load monomer structure with pLDDT. Returns (plddt_map, structure, source, path)."""
    files = config.monomer_files(gene.lower())
    cif_path = files['cif']
    pdb_path = files['pdb']
    prefer_cif = files['plddt_source'] == 'cif'

    # Try preferred source first
    if prefer_cif and cif_path:
        struct = load_cif(cif_path)
        if struct:
            plddt = get_plddt(struct, 'A')
            if plddt:
                return plddt, struct, 'cif', cif_path

    if pdb_path:
        struct = load_pdb(pdb_path)
        if struct:
            plddt = get_plddt(struct, 'A')
            if plddt:
                return plddt, struct, 'pdb', pdb_path
            return {}, struct, 'pdb_no_plddt', pdb_path

    # Fallback: try CIF if we haven't yet
    if not prefer_cif and cif_path:
        struct = load_cif(cif_path)
        if struct:
            plddt = get_plddt(struct, 'A')
            if plddt:
                return plddt, struct, 'cif', cif_path

    return {}, None, None, None


def get_multimer_plddt(pdb_path, cif_path, chain_id):
    """Get pLDDT for multimer: PDB first for consistent chain ordering."""
    struct = load_pdb(pdb_path)
    if struct:
        plddt = get_plddt(struct, chain_id)
        if plddt:
            return plddt, 'pdb'
    struct = load_cif(cif_path)
    if struct:
        plddt = get_plddt(struct, chain_id)
        if plddt:
            return plddt, 'cif'
    return {}, None

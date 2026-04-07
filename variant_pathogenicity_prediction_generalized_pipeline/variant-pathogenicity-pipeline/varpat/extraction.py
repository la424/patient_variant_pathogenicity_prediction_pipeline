"""Monomer and multimer structural feature extraction."""

import pandas as pd
import numpy as np
from typing import Dict

from .config import PipelineConfig
from .constants import (
    get_grantham, classify_grantham, grantham_severity, get_property_changes,
    classify_plddt, classify_contacts, classify_burial,
)
from .structures import (
    load_monomer, load_pdb, load_cif, count_contacts, count_interface,
    get_accessibility, get_secondary_structure, get_residue_aa,
    get_plddt, get_multimer_plddt,
)


def extract_neighborhood(contact_map, plddt_map, variant_pos, config):
    """±N neighborhood contacts with distance-weighted sum."""
    t = config.thresholds
    radius = t.get('neighborhood_radius', 3)
    weights = t.get('neighborhood_weights', {0: 1.0, 1: 0.75, 2: 0.50, 3: 0.25})
    plddt_gate = t.get('plddt_confident', 50)

    result = {
        'contacts_weighted': None,
        'contacts_raw': None,
        'evaluable': False,
        'n_eval_positions': 0,
    }
    if not contact_map or not plddt_map:
        return result
    vp = plddt_map.get(variant_pos)
    if vp is None or vp < plddt_gate:
        return result

    result['evaluable'] = True
    weighted, raw, n_eval = 0.0, 0.0, 0
    for offset in range(-radius, radius + 1):
        pos = variant_pos + offset
        w = weights.get(abs(offset), 0)
        pp = plddt_map.get(pos)
        if pp is not None and pp >= plddt_gate:
            c = contact_map.get(pos, 0)
            weighted += w * c
            raw += c
            n_eval += 1
    result['contacts_weighted'] = round(weighted, 3)
    result['contacts_raw'] = raw
    result['n_eval_positions'] = n_eval
    return result


def extract_neighborhood_multimer(contact_map, inter_map, plddt_map,
                                   variant_pos, config):
    """±N neighborhood from multimer with inter-chain contacts."""
    t = config.thresholds
    radius = t.get('neighborhood_radius', 3)
    weights = t.get('neighborhood_weights', {0: 1.0, 1: 0.75, 2: 0.50, 3: 0.25})
    plddt_gate = t.get('plddt_confident', 50)

    result = {
        'contacts_weighted': None,
        'inter_weighted': None,
        'evaluable': False,
        'has_interface': False,
    }
    if not contact_map or not plddt_map:
        return result
    vp = plddt_map.get(variant_pos)
    if vp is None or vp < plddt_gate:
        return result

    result['evaluable'] = True
    wc, wi = 0.0, 0.0
    for offset in range(-radius, radius + 1):
        pos = variant_pos + offset
        w = weights.get(abs(offset), 0)
        pp = plddt_map.get(pos)
        if pp is not None and pp >= plddt_gate:
            wc += w * contact_map.get(pos, 0)
            inter = inter_map.get(pos, 0) if inter_map else 0
            wi += w * inter
            if inter > 0:
                result['has_interface'] = True
    result['contacts_weighted'] = round(wc, 3)
    result['inter_weighted'] = round(wi, 3)
    return result


def extract_monomer_features(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Extract monomer structural metrics for all variants."""
    print("Extracting monomer structural metrics...")
    t = config.thresholds
    dist = t.get('contact_distance', 5.0)
    seq_sep = t.get('sequence_separation', 3)

    rows = []
    for gene in df['gene'].unique():
        g = str(gene).lower()
        gene_vars = df[df['gene'].str.lower() == g]

        plddt_map, struct_plddt, plddt_src, plddt_path = load_monomer(config, g)

        # PDB for contacts/accessibility
        files = config.monomer_files(g)
        pdb_path = files['pdb']
        struct_pdb = load_pdb(pdb_path)
        struct_contacts = struct_pdb or struct_plddt

        contact_map = count_contacts(struct_contacts, 'A', dist, seq_sep) if struct_contacts else {}
        acc_map = get_accessibility(struct_pdb or struct_plddt, 'A')
        ss_map = get_secondary_structure(pdb_path, 'A', config.dssp_binary) if pdb_path else {}
        aa_map = get_residue_aa(struct_plddt or struct_pdb, 'A')

        has_struct = struct_plddt is not None or struct_pdb is not None
        n_plddt = sum(1 for _, r in gene_vars.iterrows() if plddt_map.get(int(r['position'])))
        print(f"  {gene}: struct={'YES' if has_struct else 'NO'} src={plddt_src} "
              f"pLDDT={n_plddt}/{len(gene_vars)}")

        for _, row in gene_vars.iterrows():
            pos = int(row['position'])
            gd = get_grantham(row['ref_aa'], row['alt_aa'])
            p = plddt_map.get(pos)
            c = contact_map.get(pos, 0) if contact_map else None
            a = acc_map.get(pos)
            sec = ss_map.get(pos, '-') if ss_map else None

            nbhd = extract_neighborhood(contact_map, plddt_map, pos, config)

            rows.append({
                'gene': gene, 'position': pos,
                'ref_aa': row['ref_aa'], 'alt_aa': row['alt_aa'],
                'grantham_distance': gd,
                'grantham_class': classify_grantham(gd),
                'substitution_severity': round(grantham_severity(gd), 2),
                'property_changes': get_property_changes(row['ref_aa'], row['alt_aa']),
                'monomer_plddt': p,
                'monomer_plddt_category': classify_plddt(p),
                'monomer_n_contacts': float(c) if c is not None else None,
                'monomer_contact_category': classify_contacts(c),
                'monomer_aa': aa_map.get(pos),
                'monomer_accessibility': a,
                'monomer_burial': classify_burial(a),
                'monomer_secondary_structure': sec,
                'monomer_contact_disruption': float(c) if c is not None else None,
                'nbhd_mono_contacts_weighted': nbhd['contacts_weighted'],
                'nbhd_mono_contacts_raw': nbhd['contacts_raw'],
                'nbhd_mono_evaluable': nbhd['evaluable'],
                'nbhd_mono_n_eval_positions': nbhd['n_eval_positions'],
            })

    result = pd.DataFrame(rows)
    print(f"✓ Monomer: {len(result)} variants, "
          f"pLDDT={result['monomer_plddt'].notna().sum()}/{len(result)}")
    return result


def extract_multimer_features(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Extract multimer structural metrics for all variants."""
    print("Extracting multimer metrics...")
    t = config.thresholds
    dist = t.get('contact_distance', 5.0)
    seq_sep = t.get('sequence_separation', 3)

    multi_data: Dict[tuple, dict] = {}
    variant_genes = set(df['gene'].str.lower().unique())
    all_partner_labels = set()

    for mdef in config.multimer_complexes:
        gene1 = mdef['gene'].lower()
        plabel = mdef['partner'].lower()
        chain1 = mdef['chain_gene']
        chain2 = mdef['chain_partner']
        cif_file = mdef.get('cif')
        pdb_file = mdef['pdb']

        pdb_path = config.find_file(pdb_file)
        cif_path = config.find_file(cif_file) if cif_file else None

        if pdb_path is None:
            print(f"  ⚠ {pdb_file} NOT FOUND — skipping {gene1}-{plabel}")
            continue

        struct_pdb = load_pdb(pdb_path)
        if struct_pdb is None:
            print(f"  ⚠ Failed to load {pdb_file}")
            continue

        # ── Forward: gene1 variants ──
        if gene1 in variant_genes:
            all_partner_labels.add(plabel)
            positions = set(df[df['gene'].str.lower() == gene1]['position'].values)

            plddt_a, _ = get_multimer_plddt(pdb_path, cif_path, chain1)
            contacts_a = count_contacts(struct_pdb, chain1, dist, seq_sep)
            inter_a, iface_a = count_interface(struct_pdb, chain1, chain2, dist)
            acc_a = get_accessibility(struct_pdb, chain1)
            ss_a = get_secondary_structure(pdb_path, chain1, config.dssp_binary)

            for pos in positions:
                key = (gene1, pos)
                if key not in multi_data:
                    multi_data[key] = {}
                pfx = f"multi_{plabel}"

                corrected, valid = config.apply_offset(gene1, plabel, pos)
                if not valid:
                    for sfx in ['_plddt', '_n_contacts', '_inter_contacts', '_is_interface',
                                '_accessibility', '_burial', '_sec_struct', '_disruption',
                                '_nbhd_contacts_weighted', '_nbhd_inter_weighted',
                                '_nbhd_evaluable', '_nbhd_has_interface']:
                        multi_data[key][f"{pfx}{sfx}"] = None
                    multi_data[key][f"{pfx}_is_interface"] = False
                    multi_data[key][f"{pfx}_nbhd_evaluable"] = False
                    multi_data[key][f"{pfx}_nbhd_has_interface"] = False
                    continue

                p = plddt_a.get(corrected)
                c = contacts_a.get(corrected, 0)
                ic = inter_a.get(corrected, 0)

                multi_data[key][f"{pfx}_plddt"] = p
                multi_data[key][f"{pfx}_n_contacts"] = float(c)
                multi_data[key][f"{pfx}_inter_contacts"] = float(ic)
                multi_data[key][f"{pfx}_is_interface"] = corrected in iface_a
                multi_data[key][f"{pfx}_accessibility"] = acc_a.get(corrected)
                multi_data[key][f"{pfx}_burial"] = classify_burial(acc_a.get(corrected))
                multi_data[key][f"{pfx}_sec_struct"] = ss_a.get(corrected, '-')
                multi_data[key][f"{pfx}_disruption"] = float(c)

                nbhd = extract_neighborhood_multimer(
                    contacts_a, inter_a, plddt_a, corrected, config)
                multi_data[key][f"{pfx}_nbhd_contacts_weighted"] = nbhd['contacts_weighted']
                multi_data[key][f"{pfx}_nbhd_inter_weighted"] = nbhd['inter_weighted']
                multi_data[key][f"{pfx}_nbhd_evaluable"] = nbhd['evaluable']
                multi_data[key][f"{pfx}_nbhd_has_interface"] = nbhd['has_interface']

            print(f"  FWD {gene1} → multi_{plabel}")

        # ── Reverse: partner gene variants ──
        pgm = config.partner_gene_map
        partner_gene = pgm.get(plabel, plabel)

        if partner_gene in variant_genes:
            all_partner_labels.add(plabel)
            positions = set(df[df['gene'].str.lower() == partner_gene]['position'].values)

            plddt_b, _ = get_multimer_plddt(pdb_path, cif_path, chain2)
            contacts_b = count_contacts(struct_pdb, chain2, dist, seq_sep)
            inter_b, iface_b = count_interface(struct_pdb, chain2, chain1, dist)
            acc_b = get_accessibility(struct_pdb, chain2)
            ss_b = get_secondary_structure(pdb_path, chain2, config.dssp_binary)

            for pos in positions:
                key = (partner_gene, pos)
                if key not in multi_data:
                    multi_data[key] = {}
                pfx = f"multi_{plabel}"

                if f"{pfx}_plddt" in multi_data[key] and multi_data[key][f"{pfx}_plddt"] is not None:
                    continue

                corrected, valid = config.apply_offset(partner_gene, plabel, pos)
                if not valid:
                    for sfx in ['_plddt', '_n_contacts', '_inter_contacts', '_is_interface',
                                '_accessibility', '_burial', '_sec_struct', '_disruption',
                                '_nbhd_contacts_weighted', '_nbhd_inter_weighted',
                                '_nbhd_evaluable', '_nbhd_has_interface']:
                        multi_data[key][f"{pfx}{sfx}"] = None
                    multi_data[key][f"{pfx}_is_interface"] = False
                    multi_data[key][f"{pfx}_nbhd_evaluable"] = False
                    continue

                p = plddt_b.get(corrected)
                c = contacts_b.get(corrected, 0)
                ic = inter_b.get(corrected, 0)

                multi_data[key][f"{pfx}_plddt"] = p
                multi_data[key][f"{pfx}_n_contacts"] = float(c)
                multi_data[key][f"{pfx}_inter_contacts"] = float(ic)
                multi_data[key][f"{pfx}_is_interface"] = corrected in iface_b
                multi_data[key][f"{pfx}_accessibility"] = acc_b.get(corrected)
                multi_data[key][f"{pfx}_burial"] = classify_burial(acc_b.get(corrected))
                multi_data[key][f"{pfx}_sec_struct"] = ss_b.get(corrected, '-')
                multi_data[key][f"{pfx}_disruption"] = float(c)

                nbhd = extract_neighborhood_multimer(
                    contacts_b, inter_b, plddt_b, corrected, config)
                multi_data[key][f"{pfx}_nbhd_contacts_weighted"] = nbhd['contacts_weighted']
                multi_data[key][f"{pfx}_nbhd_inter_weighted"] = nbhd['inter_weighted']
                multi_data[key][f"{pfx}_nbhd_evaluable"] = nbhd['evaluable']
                multi_data[key][f"{pfx}_nbhd_has_interface"] = nbhd['has_interface']

            print(f"  REV {partner_gene} ← multi_{plabel}")

    # Merge into df
    partner_labels = sorted(all_partner_labels)
    multi_df = pd.DataFrame.from_dict(multi_data, orient='index')
    multi_df.index = pd.MultiIndex.from_tuples(multi_df.index, names=['gene_lower', 'position'])
    multi_df = multi_df.reset_index()

    df['gene_lower'] = df['gene'].str.lower()
    df = df.merge(multi_df, on=['gene_lower', 'position'], how='left', suffixes=('', '_multi'))
    df = df.drop(columns=['gene_lower'], errors='ignore')

    # Compute summary columns
    df = _compute_multimer_summary(df, partner_labels)

    print(f"✓ Multimer: {len(partner_labels)} partners, "
          f"interface variants={df['is_interface_any'].sum()}/{len(df)}")
    return df


def _compute_multimer_summary(df: pd.DataFrame, partner_labels: list) -> pd.DataFrame:
    """Compute aggregate multimer summary columns."""
    for _, row in df.iterrows():
        idx = row.name
        n_complexes = 0
        partners_with_data = []
        iface_partners = []
        plddts = []
        contacts_list = []
        disruptions = []

        for pl in partner_labels:
            p = row.get(f"multi_{pl}_plddt")
            if pd.notna(p):
                n_complexes += 1
                partners_with_data.append(pl)
                plddts.append(float(p))
            c = row.get(f"multi_{pl}_n_contacts")
            if pd.notna(c):
                contacts_list.append(float(c))
            d = row.get(f"multi_{pl}_disruption")
            if pd.notna(d):
                disruptions.append(float(d))
            iface = row.get(f"multi_{pl}_is_interface")
            if iface is True:
                iface_partners.append(pl)

        df.at[idx, 'n_multimer_complexes'] = n_complexes
        df.at[idx, 'multimer_partners'] = ';'.join(partners_with_data) if partners_with_data else ''
        df.at[idx, 'is_interface_any'] = len(iface_partners) > 0
        df.at[idx, 'interface_partners'] = ';'.join(iface_partners) if iface_partners else ''
        df.at[idx, 'n_interface_partners'] = len(iface_partners)
        df.at[idx, 'multimer_plddt_max'] = max(plddts) if plddts else None
        df.at[idx, 'multimer_plddt_avg'] = np.mean(plddts) if plddts else None
        df.at[idx, 'multimer_contacts_max'] = max(contacts_list) if contacts_list else None
        df.at[idx, 'multimer_contacts_avg'] = np.mean(contacts_list) if contacts_list else None
        df.at[idx, 'multimer_disruption_max'] = max(disruptions) if disruptions else None
        df.at[idx, 'multimer_disruption_avg'] = np.mean(disruptions) if disruptions else None

    # best_plddt
    df['best_plddt'] = df.apply(
        lambda r: max(
            float(r['monomer_plddt']) if pd.notna(r.get('monomer_plddt')) else 0,
            float(r['multimer_plddt_max']) if pd.notna(r.get('multimer_plddt_max')) else 0,
        ) or None, axis=1
    )
    return df
